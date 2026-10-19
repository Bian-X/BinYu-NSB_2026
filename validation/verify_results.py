"""Independent numeric checks on the current study outputs; requires input CSVs."""
import argparse,json,re,sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import hypergeom
from scipy.spatial.distance import pdist
from sklearn.decomposition import PCA

REPO=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(REPO/'QC_analysis/scripts'))
sys.path.insert(0,str(REPO/'sample_analysis/scripts/python'))
import qc_analysis as qa
import brain_region_pca as production

def bh(values):
    p=np.asarray(values,dtype=float);out=np.full(len(p),np.nan);ix=np.flatnonzero(np.isfinite(p))
    order=ix[np.argsort(p[ix])];n=len(order)
    out[order]=np.minimum(1,np.minimum.accumulate((p[order]*n/np.arange(1,n+1))[::-1])[::-1])
    return out

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-dir',type=Path,required=True)
    parser.add_argument('--output-root',type=Path,default=REPO)
    args=parser.parse_args();output=args.output_root;result=output/'sample_analysis/results';report={}
    paths=production.collect_region_files(args.data_dir)
    matrices=[];meta=[];differences={}
    exclusion=pd.read_csv(REPO/'sample_analysis/sample_exclusions.csv')
    excluded=set(exclusion.sample_id)
    for region,path in paths.items():
        # Independently reconstruct the calculation from raw observations.
        raw=pd.read_csv(path,low_memory=False)
        keep=raw['PG.Organisms'].fillna('').str.contains('Mus musculus',regex=False)
        keep &= ~raw['PG.ProteinAccessions'].fillna('').map(lambda s:'iRT-Kit_WR_fusion' in [x.strip() for x in s.split(';')])
        q=raw.loc[keep].set_index('PG.ProteinAccessions').filter(regex=r'\.raw\.PG\.Quantity$')
        q.columns=[re.sub(r'\.raw\.PG\.Quantity$','',re.sub(r'^\[\d+\]\s+','',c)) for c in q.columns]
        q=q.drop(columns=list(set(q.columns)&excluded)).apply(pd.to_numeric,errors='coerce')
        q=q.mask((q<=0)|~np.isfinite(q));q=q.loc[q.notna().sum(axis=1)>=20]
        x=np.log2(q);offset=x.median(axis=0);x=x-offset+offset.median()
        missing=x.isna();med=x.median(axis=1);x=x.T.fillna(med).T
        assert np.isfinite(x.to_numpy()).all()
        for col in x:
            assert np.allclose(x.loc[missing[col],col],med.loc[missing[col]],rtol=0,atol=1e-12)
        py,m,_=production.read_and_preprocess_region(path,excluded)
        pd.testing.assert_frame_equal(x,py,check_names=False,atol=1e-12,rtol=0)
        r=pd.read_csv(output/'validation/R_matrices'/f'{region}.csv',index_col=0)
        r.columns=[production.clean_sample_name(c) for c in r.columns]
        r=r.loc[x.index,x.columns]
        maximum=float(np.max(np.abs(r.to_numpy()-x.to_numpy())))
        assert maximum<1e-10,(region,maximum)
        differences[region]={'proteins':len(x),'samples':x.shape[1],'R_Python_max_abs':maximum,'imputed_cells':int(missing.sum().sum())}
        matrices.append(x);meta.append(m)
    report['regional_matrices']=differences
    joint=pd.concat(matrices,axis=1,join='inner')
    fitted=PCA(n_components=10,svd_solver='full').fit(joint.T)
    saved=pd.read_csv(result/'pca/brain_region_pca_exclude_outlier_variance.csv')
    assert np.allclose(fitted.explained_variance_ratio_,saved.ExplainedVarianceRatio,rtol=1e-10,atol=1e-12)
    centered=joint.T.to_numpy();centered=centered-centered.mean(axis=0)
    eigenvalues=np.linalg.eigvalsh(centered@centered.T)[::-1]
    assert np.allclose(eigenvalues[:10]/eigenvalues.sum(),saved.ExplainedVarianceRatio,atol=1e-12,rtol=1e-10)
    report['PCA']={'shape':list(joint.shape),'variance_percent':(fitted.explained_variance_ratio_[:2]*100).tolist()}
    differential=pd.read_csv(result/'all_regions_G5_vs_G1_limma_results_BH_FC1p5.csv')
    counts={}
    for region,table in differential.groupby('Region'):
        assert np.allclose(bh(table.PValue),table.AdjPValue,equal_nan=True,atol=1e-12,rtol=1e-10)
        selected=(table.AdjPValue<.05)&(table.Log2FC.abs()>np.log2(1.5))
        assert np.array_equal(selected,table.Status.isin(['High in G5','High in G1']))
        counts[region]=int(selected.sum())
    report['BH_candidates']=counts
    heatmap_results=pd.read_csv(result/'multi_region_G2_vs_G4_differential_results.csv')
    adjusted=[c for c in heatmap_results if c.startswith('AdjPValue')]
    assert len(adjusted)==1,adjusted
    for _,g in heatmap_results.groupby('Region'):
        assert np.allclose(bh(g.PValue),g[adjusted[0]],equal_nan=True,atol=1e-12)
    before=pd.read_csv(result/'multi_region_G2_vs_G4_heatmap_log2_group_means.csv',index_col=0)
    after=pd.read_csv(result/'multi_region_G2_vs_G4_heatmap_log2_group_means_imputed.csv',index_col=0)
    expected=before.T.fillna(before.median(axis=1)).T
    pd.testing.assert_frame_equal(expected,after,atol=1e-12,rtol=0)
    report['heatmap_median_completion']=True
    report['GO']={}
    for sub,label,prefix in [(result,'unique','PL_G5_vs_G1'),
                             (result/'expanded_mapping_sensitivity','expanded','PL_G5_vs_G1'),
                             (result/'BLA_G4_vs_G1_unique','BLA_G4_vs_G1','BLA_G4_vs_G1')]:
        go=pd.read_csv(sub/f'{prefix}_GO_enrichment_results_BH_FC1p5.csv')
        p=hypergeom.sf(go.Count-1,go.BackgroundSize,go.TermSize,go.ForegroundSize)
        assert np.allclose(p,go.PValue,atol=1e-13,rtol=1e-9)
        assert np.allclose(bh(p),go.AdjPValue_BH,atol=1e-12,rtol=1e-9)
        assert (go.loc[go.Count.eq(0),'PValue']==1).all()
        report['GO'][label]={'terms':len(go),'zero_overlap':int(go.Count.eq(0).sum()),'BH_significant':int(go.AdjPValue_BH.lt(.05).sum())}
    reports=qa.load_all_reports(args.data_dir);matrix,metadata,annotation=qa.build_combined_matrix(reports)
    q=matrix.mask((matrix<=0)|~np.isfinite(matrix));q=q.loc[q.notna().mean(axis=1)>=.8]
    x=np.log2(q);offset=x.median(axis=0);x=x-offset+offset.median();x=x.T.fillna(x.median(axis=1)).T
    pd.testing.assert_frame_equal(x,qa.prepare_complete_matrix(matrix),atol=1e-12,rtol=0)
    ordered=pd.read_csv(output/'QC_analysis/outputs/injection_order_correlation_matrix.csv',index_col=0)
    expected=x.corr().loc[ordered.index,ordered.columns]
    delta=float(np.max(np.abs(expected.to_numpy()-ordered.to_numpy())))
    assert delta<1e-12,delta
    qc=pd.read_csv(output/'QC_analysis/outputs/qc_pairwise_correlation_matrix.csv',index_col=0)
    assert np.allclose(qc,x[qc.columns].corr().loc[qc.index,qc.columns],atol=1e-12)
    report['QC']={'matrix_shape':list(x.shape),'ordered_correlation_max_abs':delta,'pairwise_consistent':True}
    target=output/'validation/numeric_verification.json';target.write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(report,ensure_ascii=False,indent=2))

if __name__=='__main__':main()
