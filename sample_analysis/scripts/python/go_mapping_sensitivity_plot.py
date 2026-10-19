"""Display expanded-mapping GO results with their protein-group provenance."""
from pathlib import Path
import argparse,json,textwrap
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

PROJECT=Path(__file__).resolve().parents[2]

def sensitivity_terms(results_dir:Path):
    expanded=results_dir/'expanded_mapping_sensitivity'
    full=pd.read_csv(expanded/'PL_G5_vs_G1_GO_enrichment_results_BH_FC1p5.csv')
    if set(full.MappingPolicy)!={'expanded'}:raise ValueError('Expected expanded-mapping GO results')
    primary=pd.read_csv(results_dir/'PL_G5_vs_G1_GO_enrichment_results_BH_FC1p5.csv')
    if set(primary.MappingPolicy)!={'unique'}:raise ValueError('Expected unique-mapping primary GO results')
    audit=pd.read_csv(expanded/'PL_G5_vs_G1_GO_mapping_audit.csv')
    foreground=audit.loc[audit.SelectedForeground & audit.IncludedByMapping].copy()
    foreground_ids=[set(str(x).split(';')) for x in foreground.MappedEntrezIDs]
    selected=full.loc[full.AdjPValue_BH.lt(.05)&full.Count.gt(0)].sort_values(
        ['AdjPValue_BH','PValue','Count','GOID'],ascending=[True,True,False,True]).copy()
    supports=[];ambiguous=[];accessions=[]
    for row in selected.itertuples():
        ids=set(str(row.EntrezIDs).split(';'))
        hit=np.array([bool(ids & mapped) for mapped in foreground_ids])
        groups=foreground.loc[hit]
        supports.append(len(groups));ambiguous.append(int(groups.MappedGeneCount.gt(1).sum()))
        accessions.append(' | '.join(groups.ProteinAccessions))
    selected['ContributingProteinGroups']=supports
    selected['AmbiguousProteinGroups']=ambiguous
    selected['ContributingAccessions']=accessions
    metadata={'MappingPolicy':'expanded','terms_tested':len(full),'BH_significant':len(selected),
              'primary_unique_BH_significant':int(primary.AdjPValue_BH.lt(.05).sum()),
              'foreground_genes':int(full.ForegroundSize.iloc[0]),'background_genes':int(full.BackgroundSize.iloc[0]),
              'interpretation':'Mapping sensitivity analysis; expanded gene IDs are not separate quantitative protein-group observations'}
    return selected,metadata

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--results-dir',type=Path,default=PROJECT/'results')
    p.add_argument('--figure-dir',type=Path,default=PROJECT/'figures')
    p.add_argument('--top',type=int,default=10)
    args=p.parse_args()
    if args.top<1:raise ValueError('top must be positive')
    selected,meta=sensitivity_terms(args.results_dir)
    expanded=args.results_dir/'expanded_mapping_sensitivity'
    selected.to_csv(expanded/'PL_GO_expanded_gene_to_protein_group_support.csv',index=False)
    (expanded/'PL_GO_expanded_display.json').write_text(json.dumps(meta,indent=2),encoding='utf-8')
    args.figure_dir.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.family':'Arial','font.size':8,'pdf.fonttype':42})
    fig,ax=plt.subplots(figsize=(8.2,5.7));fig.subplots_adjust(left=.48,right=.81,bottom=.20,top=.79)
    g=selected.head(args.top).iloc[::-1];y=np.arange(len(g))
    if g.empty:
        ax.axis('off');ax.text(.5,.5,'No expanded-mapping terms passed BH < 0.05',ha='center',transform=ax.transAxes)
    else:
        sc=ax.scatter(g.RichFactor,y,c=-np.log10(np.maximum(g.AdjPValue_BH,1e-300)),s=12+3*g.Count,
                      cmap='viridis',edgecolors='#333333',linewidths=.5)
        ax.set_yticks(y,['\n'.join(textwrap.wrap(r.TERM,35))+f" [{r.ContributingProteinGroups} {'PG' if r.ContributingProteinGroups==1 else 'PGs'}]" for r in g.itertuples()],fontsize=7)
        ax.set_xlabel('Rich factor (expanded IDs / term IDs)')
        ax.set_xlim(0,min(1.,max(.1,float(g.RichFactor.max())*1.15)))
        ax.spines[['top','right']].set_visible(False);ax.grid(axis='x',color='#eeeeee');ax.set_axisbelow(True)
        cax=fig.add_axes([.85,.40,.018,.29]);fig.colorbar(sc,cax=cax,label='−log10 BH-adjusted P')
        counts=sorted({int(g.Count.min()),int(g.Count.median()),int(g.Count.max())})
        handles=[ax.scatter([],[],s=12+3*n,c='#c4c9ce',edgecolors='#333333',linewidths=.5,label=str(n)) for n in counts]
        ax.legend(handles=handles,title='Expanded IDs',loc='lower right',frameon=False,fontsize=7,title_fontsize=7)
    fig.suptitle('PL G5−G1: expanded-mapping sensitivity analysis',fontsize=11,y=.965)
    fig.text(.5,.90,f"{meta['BH_significant']} / {meta['terms_tested']:,} terms passed BH < 0.05 under expanded mapping",ha='center',fontsize=9)
    fig.text(.5,.855,f"Primary unique mapping: {meta['primary_unique_BH_significant']} terms passed BH < 0.05",ha='center',fontsize=9,fontweight='bold')
    fig.text(.5,.07,f'Top {len(g)} terms ranked by BH-adjusted P; PGs = contributing protein-group entries.\nExpanded IDs can share one quantity and are not separate protein measurements.',ha='center',fontsize=8,color='#555555')
    target=args.figure_dir/'PL_GO_expanded_mapping_sensitivity'
    fig.savefig(target.with_suffix('.pdf'));fig.savefig(target.with_suffix('.png'),dpi=400);plt.close(fig)
    print(json.dumps(meta,indent=2))

if __name__=='__main__':main()
