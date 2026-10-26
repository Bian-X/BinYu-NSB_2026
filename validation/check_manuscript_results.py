"""Check calculated panel data against the selected manuscript's numeric claims."""
import argparse,json,sys
from pathlib import Path
import numpy as np
import pandas as pd

def check_outputs(root,scope,strict=True,selection_sensitivity=False):
    checks=[]
    def test(name,actual,expected,tolerance=0):
        ok=bool(np.isclose(actual,expected,atol=tolerance,rtol=0)) if isinstance(expected,(int,float)) else actual==expected
        checks.append({'claim':name,'actual':actual,'expected':expected,'passed':ok})
    if scope in ('all','figure2'):
        q=root/'QC_analysis/outputs'
        samples=pd.read_csv(q/'sample_qc_observed.csv');cv=pd.read_csv(q/'qc_cv_observed.csv')
        eligible=cv.loc[cv.qc_detected_n.ge(20)]
        test('all acquisitions',len(samples),289)
        test('minimum detected protein groups',int(samples.detected_protein_groups.min()),6143)
        test('maximum detected protein groups',int(samples.detected_protein_groups.max()),9294)
        test('median detected protein groups',float(samples.detected_protein_groups.median()),8940)
        test('QC CV eligible exported entries',len(eligible),8908)
        test('median QC CV percent',float(eligible.qc_cv_percent.median()),9.89,.005)
        test('QC CV <=20 percent',int(eligible.qc_cv_percent.le(20).sum()),6957)
        test('QC CV <=30 percent',int(eligible.qc_cv_percent.le(30).sum()),8050)
        test('QC CV >100 percent',int(eligible.qc_cv_percent.gt(100).sum()),4)
        regression=json.loads((q/'qc_protein_injection_order_regression.json').read_text())
        # Reference metadata follows the retained RAW acquisition chronology (QC_13 at 141).
        for key,value,tol in [('n',8880,0),('span',288,0),('median_ratio',1.066,.0005),('p95_ratio',1.369,.0005),('ge1p5',249,0),('pct_ge1p5',2.80,.005),('excluded_from_regression',28,0)]:
            test('QC regression '+key,regression[key],value,tol)
        audit=pd.read_csv(q/'qc_regression_selection_audit.csv')
        test('regression exclusion count is CV minus regression',int(audit.CVEligible.sum()-audit.RegressionIncluded.sum()),28)
        correlations=pd.read_csv(q/'qc_pairwise_correlation_matrix.csv',index_col=0)
        upper=correlations.to_numpy()[np.triu_indices(len(correlations),1)]
        test('QC median pairwise Pearson r',float(np.median(upper)),.9912,.00005)
    if scope in ('all','figure3'):
        r=root/'sample_analysis/results'
        pl=pd.read_csv(r/'PL_G5_vs_G1_differential_results_BH_FC1p5.csv')
        test('PL tested protein groups',len(pl),8669)
        test('PL higher in G5',int(pl.Status.eq('High in G5').sum()),153)
        test('PL lower in G5',int(pl.Status.eq('Low in G5').sum()),98)
        go=pd.read_csv(r/'BLA_G4_vs_G1_unique/BLA_G4_vs_G1_GO_enrichment_results_BH_FC1p5.csv')
        test('BLA GO eligible terms',len(go),9097)
        test('BLA GO BH significant terms',int(go.AdjPValue_BH.lt(.05).sum()),37)
        test('BLA GO foreground genes',int(go.ForegroundSize.iloc[0]),400)
        test('BLA GO background genes',int(go.BackgroundSize.iloc[0]),8368)
        expanded=pd.read_csv(r/'expanded_mapping_sensitivity/PL_G5_vs_G1_GO_enrichment_results_BH_FC1p5.csv')
        test('PL expanded-mapping GO BH terms',int(expanded.AdjPValue_BH.lt(.05).sum()),17)
        for mode,n,pc1,pc2 in [('exclude_outlier',263,37.32,15.59),('include_all',264,37.59,15.51)]:
            scores=pd.read_csv(r/f'pca/brain_region_pca_{mode}_scores.csv')
            variance=pd.read_csv(r/f'pca/brain_region_pca_{mode}_variance.csv')
            test('PCA '+mode+' sample count',len(scores),n)
            test('PCA '+mode+' PC1 percent',float(variance.ExplainedVariancePercent.iloc[0]),pc1,.005)
            test('PCA '+mode+' PC2 percent',float(variance.ExplainedVariancePercent.iloc[1]),pc2,.005)
        z=pd.read_csv(r/'multi_region_G2_vs_G4_heatmap_zscore_matrix.csv',index_col=0)
        test('exploratory heatmap rows',len(z),36)
        corr=pd.read_csv(r/'brain_region_log2fc_correlation_matrix_G5_G1.csv',index_col=0)
        test('regional effect correlation dimensions',list(corr.shape),[11,11])
        screen=r/'workflow_example_comparison'
        if selection_sensitivity:
            if not (screen/'screen_metadata.csv').is_file():
                raise FileNotFoundError(f'Expected optional screen output: {screen / "screen_metadata.csv"}')
            screen_meta=pd.read_csv(screen/'screen_metadata.csv')
            test('exploratory screen comparisons',int(screen_meta.Comparisons.iloc[0]),165)
            test('exploratory screen GO tests',int(screen_meta.GOTests.iloc[0]),1510215)
            summary=pd.read_csv(screen/'comparison_summary.csv')
            bla=summary.loc[summary.Region.eq('BLA') & summary.GroupA.eq('G1') & summary.GroupB.eq('G4')]
            test('BLA GO terms after across-screen BH',int(bla.GO_AcrossSearch_BH.iloc[0]),23)
            complete=screen/'shortlist_no_imputation/BLA_G4_vs_G1_GO.csv'
            if not complete.is_file():
                raise FileNotFoundError(f'Expected optional complete-case output: {complete}')
            nofill=pd.read_csv(complete)
            selected=set(nofill.loc[nofill.AdjPValue_BH.lt(.05),'GOID'])
            main=set(go.loc[go.AdjPValue_BH.lt(.05),'GOID'])
            test('BLA complete-case GO terms',len(selected),41)
            test('BLA complete-case retains main GO terms',len(selected & main),37)
    manifest={'scope':scope,'strict_manuscript_match':strict,'selection_sensitivity':selection_sensitivity,'checks':checks,'all_passed':all(c['passed'] for c in checks)}
    folder=root/'validation';folder.mkdir(parents=True,exist_ok=True)
    (folder/'manuscript_result_checks.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    return manifest

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output-root',type=Path,required=True)
    p.add_argument('--figures',choices=['all','figure2','figure3'],default='all')
    p.add_argument('--environment',choices=['manuscript','compatible'],default='manuscript')
    p.add_argument('--selection-sensitivity',action='store_true');a=p.parse_args()
    result=check_outputs(a.output_root,a.figures,a.environment=='manuscript',a.selection_sensitivity)
    print(json.dumps(result,indent=2))
    if not result['all_passed'] and result['strict_manuscript_match']:sys.exit(1)
if __name__=='__main__':main()
