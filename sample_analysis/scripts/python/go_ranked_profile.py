"""Describe ranked GO associations while retaining the full BH conclusion.

This display does not select significant terms or change the tested universe.
"""
from pathlib import Path
import argparse,json,textwrap
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

PROJECT=Path(__file__).resolve().parents[2]

def ranked_terms(table:pd.DataFrame,top:int=10)->pd.DataFrame:
    if top<1:raise ValueError('top must be positive')
    valid=table.loc[table.Count.gt(0)&np.isfinite(table.PValue)&table.PValue.gt(0)].copy()
    valid=valid.sort_values(['PValue','Count','GOID'],ascending=[True,False,True]).head(top)
    valid.insert(0,'NominalRank',np.arange(1,len(valid)+1))
    valid['ExpectedCount']=valid.ForegroundSize*valid.TermSize/valid.BackgroundSize
    valid['FoldEnrichment']=valid.Count/valid.ExpectedCount
    valid['DisplayInterpretation']='Descriptive ranking; not an FDR-significant set'
    return valid

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--results-dir',type=Path,default=PROJECT/'results')
    p.add_argument('--figure-dir',type=Path,default=PROJECT/'figures')
    p.add_argument('--top',type=int,default=10)
    args=p.parse_args();args.figure_dir.mkdir(parents=True,exist_ok=True)
    table=pd.read_csv(args.results_dir/'PL_G5_vs_G1_GO_enrichment_results_BH_FC1p5.csv')
    if set(table.MappingPolicy)!={'unique'}:raise ValueError('This profile requires the unique-mapping primary analysis')
    ranked=ranked_terms(table,args.top)
    ranked.to_csv(args.results_dir/'PL_G5_vs_G1_GO_ranked_descriptive.csv',index=False)
    significant=int(table.AdjPValue_BH.lt(.05).sum())
    plt.rcParams.update({'font.family':'Arial','font.size':8,'pdf.fonttype':42})
    fig,ax=plt.subplots(figsize=(8.2,5.6))
    fig.subplots_adjust(left=.49,right=.82,bottom=.18,top=.80)
    g=ranked.iloc[::-1];y=np.arange(len(g))
    sc=ax.scatter(g.FoldEnrichment,y,c=-np.log10(g.PValue),s=18+8*g.Count,
                  cmap='viridis',edgecolors='#333333',linewidths=.5)
    ax.set_yticks(y,['\n'.join(textwrap.wrap(f'{r.TERM} [{r.Ontology}]',38)) for r in g.itertuples()],fontsize=7.2)
    ax.set_xlabel('Fold enrichment (observed / expected)')
    ax.set_xlim(0,max(1,float(g.FoldEnrichment.max()))*1.13)
    ax.axvline(1,color='#aaaaaa',ls='--',lw=.7,zorder=0)
    ax.grid(axis='x',color='#eeeeee');ax.set_axisbelow(True)
    ax.spines[['top','right']].set_visible(False)
    cax=fig.add_axes([.86,.39,.018,.30]);fig.colorbar(sc,cax=cax,label='-log10 raw P')
    sizes=sorted(set(int(v) for v in [g.Count.min(),g.Count.median(),g.Count.max()]))
    handles=[ax.scatter([],[],s=18+8*n,c='#c4c9ce',edgecolors='#333333',linewidths=.5,label=str(n)) for n in sizes]
    ax.legend(handles=handles,title='Gene count',loc='lower right',frameon=False,fontsize=7,title_fontsize=7)
    fig.suptitle('PL G5-G1: exploratory GO functional profile',fontsize=12,y=.95)
    fig.text(.5,.875,f'{significant} / {len(table):,} GO terms passed BH < 0.05',ha='center',fontsize=10,fontweight='bold')
    bh_range=f'BH-adjusted P of shown terms: {g.AdjPValue_BH.min():.2f}'
    if g.AdjPValue_BH.max()!=g.AdjPValue_BH.min():bh_range+=f' to {g.AdjPValue_BH.max():.2f}'
    fig.text(.5,.835,bh_range,ha='center',fontsize=8)
    fig.text(.5,.055,f'Top {len(ranked)} terms ranked by raw P, displayed for description only.\nUnique-gene mapping; the full GO test family and BH values are unchanged.',ha='center',fontsize=8,color='#555555')
    name=args.figure_dir/'PL_GO_ranked_descriptive_G5_vs_G1'
    fig.savefig(name.with_suffix('.pdf'));fig.savefig(name.with_suffix('.png'),dpi=400);plt.close(fig)
    metadata={'mapping':'unique','display':'top terms ranked by raw P; descriptive only','terms_tested':len(table),
              'BH_significant':significant,'displayed':len(ranked),'minimum_raw_P':float(table.PValue.min()),
              'nominal_P_below_005':int(table.PValue.lt(.05).sum()),'BH_values_shown':sorted(ranked.AdjPValue_BH.unique().tolist()),
              'fold_enrichment_definition':'Count / (ForegroundSize * TermSize / BackgroundSize)'}
    (args.results_dir/'PL_G5_vs_G1_GO_descriptive_display.json').write_text(json.dumps(metadata,indent=2),encoding='utf-8')
    print(json.dumps(metadata,indent=2))

if __name__=='__main__':main()
