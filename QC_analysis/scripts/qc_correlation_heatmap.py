#!/usr/bin/env python3
"""QC pairwise correlations ordered by injection sequence (descriptive)."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from figure_2F_injection_order import build_figure_2f_correlation, attach_injection_order


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--data-dir', type=Path, required=True)
    p.add_argument('--metadata', type=Path, required=True)
    p.add_argument('--out-dir', type=Path, required=True)
    args = p.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    corr, _ = build_figure_2f_correlation(args.data_dir)
    corr, key = attach_injection_order(corr, args.metadata)
    key = key.loc[key.sample_type.eq('QC')]
    ids = list(key.index)
    if len(ids) < 3:
        raise ValueError('At least three QC injections are required')
    qc = corr.loc[ids, ids]
    qc.to_csv(args.out_dir / 'qc_pairwise_correlation_matrix.csv')
    pairs = pd.DataFrame([{'sample_a': ids[i], 'sample_b': ids[j],
                           'order_distance': int(abs(key.loc[ids[i], 'injection_order'] - key.loc[ids[j], 'injection_order'])),
                           'correlation': float(qc.iloc[i, j])}
                          for i in range(len(ids)) for j in range(i+1, len(ids))])
    pairs.to_csv(args.out_dir / 'qc_pairwise_correlation_pairs.csv', index=False)
    fit = stats.linregress(pairs.order_distance, pairs.correlation)
    pd.DataFrame([{'QCCount': len(ids), 'PairCount': len(pairs), 'MedianCorrelation': pairs.correlation.median(),
                   'MinimumCorrelation': pairs.correlation.min(), 'MaximumCorrelation': pairs.correlation.max(),
                   'SlopePerOrderUnit_Descriptive': fit.slope, 'R2_Descriptive': fit.rvalue**2,
                   'RegressionPValue': np.nan,
                   'InferenceNote': 'Pairs share QC injections; ordinary independent-pair regression P is not reported'}]
                 ).to_csv(args.out_dir / 'qc_pairwise_correlation_metrics.csv', index=False)
    labels = [f"QC {i+1} | #{key.loc[s, 'injection_order']}" for i, s in enumerate(ids)]
    fig, ax = plt.subplots(figsize=(10, 9), layout='constrained')
    image = ax.imshow(qc, vmin=float(pairs.correlation.min())-0.001, vmax=1, cmap='Blues')
    ax.set_xticks(range(len(ids)), labels, rotation=90, fontsize=6)
    ax.set_yticks(range(len(ids)), labels, fontsize=6)
    for i in range(len(ids)):
        for j in range(len(ids)):
            ax.text(j, i, f'{qc.iloc[i,j]:.4f}', ha='center', va='center', fontsize=4.3,
                    color='white' if qc.iloc[i,j] > 0.994 else 'black')
    ax.set_title('Pooled-QC Pearson correlations by injection order')
    fig.colorbar(image, ax=ax, label='Pearson r')
    fig.savefig(args.out_dir / 'figure_QC_pairwise_correlation.pdf', dpi=600)
    fig.savefig(args.out_dir / 'figure_QC_pairwise_correlation.png', dpi=300)
    plt.close(fig)


if __name__ == '__main__':
    main()
