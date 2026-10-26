#!/usr/bin/env python3
"""Chronological QC coverage and exported PG.Quantity summaries.

PG.Quantity is the exported, software-processed protein-group signal, not
raw instrument intensity. Optional RT measurements and instrument events
must come from actual acquisition records; missing records are reported.
"""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import qc_analysis as qa
from metadata_io import read_sample_metadata, match_metadata


def write_event_records(out_dir: Path, events: pd.DataFrame | None) -> None:
    """Keep the exported event table consistent with this run's inputs."""
    path = out_dir / 'provided_acquisition_events.csv'
    if events is None:
        path.unlink(missing_ok=True)
    else:
        events.to_csv(path, index=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-dir', type=Path, required=True)
    parser.add_argument('--metadata', type=Path, required=True)
    parser.add_argument('--out-dir', type=Path, required=True)
    parser.add_argument('--events', type=Path, help='Actual event CSV: injection_order,event_type,description')
    parser.add_argument('--rt-metrics', type=Path, help='Measured QC RT CSV: sample_id,rt_metric')
    parser.add_argument('--rt-label', default='Provided retention-time metric')
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    reports = qa.load_all_reports(args.data_dir)
    pd.concat([report.sample_metadata for report in reports], ignore_index=True)[
        ['source_quantity_column', 'source_raw_file_identifier', 'sample_id']
    ].to_csv(args.out_dir / 'input_label_mapping.csv', index=False)
    matrix, _, _ = qa.build_combined_matrix(reports)
    key = match_metadata(read_sample_metadata(args.metadata), matrix.columns)
    metrics = qa.sample_qc_summary_from_reports(reports).set_index('sample_id')
    qc = key.loc[key.sample_type.eq('QC')].join(
        metrics[['detected_protein_groups', 'missing_percent', 'total_intensity', 'median_intensity']])
    qc = qc.sort_values('injection_order')
    qc['log2_total_exported_quantity'] = np.log2(qc.total_intensity)
    qc['log2_median_exported_quantity'] = np.log2(qc.median_intensity)
    # Mouse protein groups with >=20 observed pooled-QC values; no completion.
    qc_report = next(report for report in reports if report.is_qc_report)
    selection = qa.qc_regression_selection(qc_report)
    selection.to_csv(args.out_dir / 'qc_regression_selection_audit.csv', index=False)
    observed = qc_report.quantity.loc[selection.RegressionIncluded.to_numpy()].reindex(columns=qc.index)
    x = qc.injection_order.to_numpy(float)
    span = float(np.ptp(x))
    slopes = []
    for row in np.log2(observed.to_numpy(float)):
        valid = np.isfinite(row)
        slopes.append(float(np.polyfit(x[valid], row[valid], 1)[0]))
    ratios = np.exp2(np.abs(slopes) * span)
    pd.DataFrame({'accession': observed.index, 'slope': slopes,
                  'fitted_high_low_ratio': ratios}).to_csv(args.out_dir / 'qc_protein_injection_order_regression.csv', index=False)
    regression_summary = {'n': len(ratios), 'span': span,
                          'median_ratio': float(np.median(ratios)),
                          'p95_ratio': float(np.percentile(ratios, 95)),
                          'ge1p5': int((ratios >= 1.5).sum()),
                          'pct_ge1p5': float(100 * np.mean(ratios >= 1.5)),
                          'limitation': 'Processed quantities; no nonlinear-drift or causal batch inference.'}
    regression_summary['cv_eligible_entries'] = int(selection.CVEligible.sum())
    regression_summary['excluded_from_regression'] = int((selection.CVEligible & ~selection.RegressionIncluded).sum())
    regression_summary['organism_filter'] = "PG.Organisms contains Mus musculus"
    (args.out_dir / 'qc_protein_injection_order_regression.json').write_text(json.dumps(regression_summary, indent=2), encoding='utf-8')
    events = None
    if args.events:
        events = pd.read_csv(args.events)
        if not {'injection_order', 'event_type', 'description'}.issubset(events):
            raise ValueError('Events require injection_order,event_type,description')
        order = pd.to_numeric(events.injection_order, errors='coerce')
        if order.isna().any() or not np.isfinite(order).all() or not order.between(key.injection_order.min(), key.injection_order.max()).all():
            raise ValueError('Event order must be finite and within the acquired sequence')
        events['injection_order'] = order
    rt = None
    if args.rt_metrics:
        rt = pd.read_csv(args.rt_metrics).set_index('sample_id')
        if 'rt_metric' not in rt or rt.index.duplicated().any():
            raise ValueError('RT metrics require unique sample_id and rt_metric')
        if not rt.index.isin(qc.index).all():
            raise ValueError('RT metrics contain non-QC or unmatched sample IDs')
        rt['rt_metric'] = pd.to_numeric(rt.rt_metric, errors='coerce')
        if not np.isfinite(rt.rt_metric).all():
            raise ValueError('RT metrics must be finite measured values')
        qc = qc.join(rt[['rt_metric']])
    qc.to_csv(args.out_dir / 'qc_time_order_metrics.csv')
    fig, axes = plt.subplots(5 if rt is not None else 4, 1, figsize=(8, 9), sharex=True, layout='constrained')
    order = qc.injection_order
    axes[0].plot(order, qc.log2_total_exported_quantity, 'o-', label='log2 sum(PG.Quantity)')
    axes[0].set_ylabel('log2 sum(PG.Quantity)')
    axes[0].legend(fontsize=8)
    axes[1].plot(order, qc.log2_median_exported_quantity, 's-', color='darkorange')
    axes[1].set_ylabel('log2 median(PG.Quantity)')
    axes[2].plot(order, qc.detected_protein_groups, 'o-')
    axes[2].set_ylabel('Detected protein groups')
    axes[3].plot(order, qc.missing_percent, 'o-')
    axes[3].set_ylabel('Within-report missingness (%)')
    if rt is not None:
        axes[4].plot(order, qc.rt_metric, 'o-')
        axes[4].set_ylabel(args.rt_label)
    for ax in axes:
        ax.grid(alpha=0.25)
        if events is not None:
            for event in events.itertuples():
                ax.axvline(event.injection_order, color='darkred', linestyle='--', alpha=0.5)
    if events is not None:
        for event in events.itertuples():
            axes[0].text(event.injection_order, 1, event.event_type, transform=axes[0].get_xaxis_transform(),
                         rotation=90, fontsize=7, va='top')
    write_event_records(args.out_dir, events)
    axes[-1].set_xlabel('LC-MS/MS injection order')
    fig.suptitle('Chronological pooled-QC summary\nExported PG.Quantity; RT/events shown only when provided')
    fig.savefig(args.out_dir / 'figure_QC_time_order_summary.pdf', dpi=600)
    fig.savefig(args.out_dir / 'figure_QC_time_order_summary.png', dpi=300)
    plt.close(fig)
    availability = {'qc_injections': len(qc), 'exported_quantity_and_coverage_available': True,
                    'raw_instrument_intensity_verified': False,
                    'measured_rt_provided': rt is not None,
                    'qc_with_measured_rt': int(qc.rt_metric.notna().sum()) if rt is not None else 0,
                    'acquisition_events_provided': events is not None,
                    'acquisition_date_or_elapsed_time_present': any(c in key for c in ('acquisition_date', 'elapsed_time')),
                    'limitation': 'No absolute intensity-drift bound or causal batch correction is inferred.'}
    (args.out_dir / 'qc_evidence_availability.json').write_text(json.dumps(availability, indent=2), encoding='utf-8')
    print(json.dumps(availability, indent=2))


if __name__ == '__main__':
    main()
