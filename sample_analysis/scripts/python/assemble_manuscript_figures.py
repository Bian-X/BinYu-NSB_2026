"""Export separate scientific manuscript panels from regenerated result tables.

Publication scope is Figure 2D-F and Figure 3. Behavioral panels and image
assets are outside this package. All panels are calculated from current inputs.
"""
from pathlib import Path
import argparse
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import gridspec
from PIL import Image

RELEASE=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(RELEASE/'QC_analysis/scripts'))
import qc_analysis as qa
from brain_region_pca import REGION_PALETTE,add_confidence_ellipse,square_pca_axes
from go_ranked_profile import ranked_terms
from go_mapping_sensitivity_plot import sensitivity_terms

def panel(ax,label,y=1.08):
    ax.text(-.12*min(1.,.45/ax.get_position().width),y,label,transform=ax.transAxes,fontsize=12,fontweight='bold',va='top')

def save_formats(fig, stem, bbox=None, dpi=400):
    """Vector labels/lines; configurable PNG/TIFF resolution (default 400 dpi)."""
    for extension in ('pdf','svg','png'):
        fig.savefig(stem.with_suffix('.'+extension),dpi=dpi,bbox_inches=bbox)
    with Image.open(stem.with_suffix('.png')) as im:
        im.save(stem.with_suffix('.tif'),compression='tiff_lzw',dpi=(dpi,dpi))


def save_panels(fig, groups, output_dir):
    """Export original plot artists and vector text, one panel per file.

    Each group includes its axes, colorbars and inset axes. Bounding boxes
    include tick labels, legends and panel letters; no bitmap is sliced.
    """
    from matplotlib.transforms import Bbox
    for name, axes in groups.items():
        for ax in fig.axes:ax.set_visible(ax in axes)
        fig.canvas.draw()
        renderer=fig.canvas.get_renderer()
        boxes=[ax.get_tightbbox(renderer) for ax in axes]
        bbox=Bbox.union([box for box in boxes if box is not None]).transformed(fig.dpi_scale_trans.inverted())
        padding=.08
        bbox=Bbox.from_extents(bbox.x0-padding,bbox.y0-padding,bbox.x1+padding,bbox.y1+padding)
        if name=='Figure3C_PCA':
            side=max(bbox.width,bbox.height)
            xmid=(bbox.x0+bbox.x1)/2;ymid=(bbox.y0+bbox.y1)/2
            bbox=Bbox.from_extents(xmid-side/2,ymid-side/2,xmid+side/2,ymid+side/2)
        save_formats(fig, output_dir/name, bbox=bbox)
    for ax in fig.axes:ax.set_visible(True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--data-dir',type=Path,required=True)
    p.add_argument('--results-dir',type=Path,default=RELEASE/'sample_analysis/results')
    p.add_argument('--figure',choices=['all','figure2','figure3'],default='all')
    p.add_argument('--output-dir',type=Path,default=RELEASE/'manuscript_figures')
    p.add_argument('--go-display',choices=['bla-g4-g1','expanded-sensitivity','unique'],default='bla-g4-g1')
    p.add_argument('--font-family',default='Arial')
    args=p.parse_args();args.output_dir.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.family':args.font_family,'font.size':7,'axes.labelsize':7,'xtick.labelsize':6,'ytick.labelsize':6,'pdf.fonttype':42})
    if args.figure in ('all','figure3'):
        assemble_figure3(args, args.results_dir)
    if args.figure in ('all','figure2'):
        assemble_figure2(args)


def assemble_figure3(args, result):
    pl=pd.read_csv(result/'PL_G5_vs_G1_differential_results_BH_FC1p5.csv')
    go=pd.read_csv(result/'PL_G5_vs_G1_GO_selected_adjusted.csv')
    scores=pd.read_csv(result/'pca/brain_region_pca_exclude_outlier_scores.csv')
    variance=pd.read_csv(result/'pca/brain_region_pca_exclude_outlier_variance.csv')
    z=pd.read_csv(result/'multi_region_G2_vs_G4_heatmap_zscore_matrix.csv',index_col=0)
    info=pd.read_csv(result/'multi_region_G2_vs_G4_heatmap_selected_proteins.csv').set_index('ProteinAccessions')
    corr=pd.read_csv(result/'brain_region_log2fc_correlation_matrix_G5_G1.csv',index_col=0)
    from scipy.cluster.hierarchy import linkage,leaves_list
    from scipy.spatial.distance import squareform
    clustered=leaves_list(linkage(squareform(1-corr.to_numpy(),checks=False),method='complete'))
    corr=corr.iloc[clustered,clustered]
    fig=plt.figure(figsize=(7.2,10.0));gs=gridspec.GridSpec(3,2,figure=fig,height_ratios=[1,1.55,.8],hspace=.52,wspace=.55)
    groups={};start=len(fig.axes)
    ax=fig.add_subplot(gs[0,0]);panel(ax,'A')
    colors=np.where(pl.Status.eq('High in G5'),'#D94B3D',np.where(pl.Status.eq('Low in G5'),'#2878B5','#C3C7CB'))
    ax.scatter(pl.Log2FC,-np.log10(np.maximum(pl.AdjPValue,1e-300)),c=colors,s=3,alpha=.7,
               edgecolors='none',linewidths=0,rasterized=True)
    ax.axhline(-np.log10(.05),ls='--',lw=.6,color='grey')
    for x in [-np.log2(1.5),np.log2(1.5)]:ax.axvline(x,ls='--',lw=.6,color='grey')
    ax.set(xlabel='log2 fold change (G5/G1)',ylabel='−log10 BH-adjusted P')
    ax.text(.02,.98,f"PL: {pl.Status.eq('High in G5').sum()} higher / {pl.Status.eq('Low in G5').sum()} lower",transform=ax.transAxes,va='top',fontsize=6)
    groups['Figure3A_PL_volcano']=fig.axes[start:];start=len(fig.axes)
    ax=fig.add_subplot(gs[0,1]);panel(ax,'B')
    if args.go_display=='bla-g4-g1':
        selected=pd.read_csv(result/'BLA_G4_vs_G1_unique/BLA_G4_vs_G1_GO_selected_adjusted.csv')
        if not selected.MappingPolicy.eq('unique').all():raise ValueError('BLA example requires unique mapping')
        g=selected.sort_values(['AdjPValue_BH','PValue','GOID']).head(6).iloc[::-1];y=np.arange(len(g))
        import textwrap
        if g.empty:
            ax.axis('off');ax.text(.5,.5,'No BLA G4-G1 terms passed BH < 0.05',ha='center',transform=ax.transAxes,fontsize=6)
        else:
            sc=ax.scatter(g.RichFactor,y,c=-np.log10(np.maximum(g.AdjPValue_BH,1e-300)),s=10+g.Count,
                          cmap='viridis',edgecolors='#333333',linewidths=.3)
            ax.set_yticks(y,['\n'.join(textwrap.wrap(t,25)) for t in g.TERM],fontsize=4.7)
            ax.set_xlabel('Rich factor (gene count / term size)',fontsize=5.5)
            ax.set_xlim(0,float(g.RichFactor.max())*1.20)
            counts=sorted({int(g.Count.min()),int(g.Count.max())})
            handles=[ax.scatter([],[],s=10+n,c='#c4c9ce',edgecolors='#333333',linewidths=.3,label=str(n)) for n in counts]
            ax.legend(handles=handles,title='Unique IDs',loc='lower right',frameon=False,fontsize=4.5,title_fontsize=4.5,handletextpad=.4)
            fig.colorbar(sc,ax=ax,fraction=.04,pad=.03,label='−log10 BH P')
        ax.set_title(f'BLA G4−G1 GO: unique mapping\n{len(selected)} terms passed BH < 0.05',fontsize=5.4,pad=3)
    elif args.go_display=='expanded-sensitivity':
        expanded,meta=sensitivity_terms(result)
        g=expanded.head(6).iloc[::-1];y=np.arange(len(g))
        import textwrap
        if g.empty:
            ax.axis('off');ax.text(.5,.5,'No expanded-mapping terms passed BH < 0.05',ha='center',transform=ax.transAxes,fontsize=6)
        else:
            sc=ax.scatter(g.RichFactor,y,c=-np.log10(np.maximum(g.AdjPValue_BH,1e-300)),s=10+1.5*g.Count,
                          cmap='viridis',edgecolors='#333333',linewidths=.3)
            ax.set_yticks(y,['\n'.join(textwrap.wrap(r.TERM,25))+f" [{r.ContributingProteinGroups} {'PG' if r.ContributingProteinGroups==1 else 'PGs'}]" for r in g.itertuples()],fontsize=4.7)
            ax.set_xlabel('Rich factor (expanded IDs / term IDs)',fontsize=5.5)
            ax.set_xlim(0,min(1.,max(.1,float(g.RichFactor.max())*1.15)))
            counts=sorted({int(g.Count.min()),int(g.Count.max())})
            handles=[ax.scatter([],[],s=10+1.5*n,c='#c4c9ce',edgecolors='#333333',linewidths=.3,label=str(n)) for n in counts]
            ax.legend(handles=handles,title='Expanded IDs',loc='lower left',frameon=False,fontsize=4.5,title_fontsize=4.5,handletextpad=.4)
            fig.colorbar(sc,ax=ax,fraction=.04,pad=.03,label='−log10 BH P')
        ax.set_title(f"Expanded-mapping sensitivity\n{meta['BH_significant']} BH terms; primary unique: {meta['primary_unique_BH_significant']}",fontsize=5.4,pad=3)
    elif go.empty:
        full_go=pd.read_csv(result/'PL_G5_vs_G1_GO_enrichment_results_BH_FC1p5.csv')
        g=ranked_terms(full_go,6).iloc[::-1];y=np.arange(len(g))
        sc=ax.scatter(g.FoldEnrichment,y,c=-np.log10(g.PValue),s=10+4*g.Count,
                      cmap='viridis',edgecolors='#333333',linewidths=.3)
        import textwrap
        ax.set_yticks(y,['\n'.join(textwrap.wrap(t,25)) for t in g.TERM],fontsize=4.7)
        ax.set_xlabel('Fold enrichment',fontsize=6)
        ax.set_xlim(0,max(1,float(g.FoldEnrichment.max()))*1.16)
        ax.axvline(1,color='#aaaaaa',ls='--',lw=.5,zorder=0)
        shown_bh=f'{g.AdjPValue_BH.min():.2f}'
        if g.AdjPValue_BH.max()!=g.AdjPValue_BH.min():shown_bh+=f'–{g.AdjPValue_BH.max():.2f}'
        ax.set_title(f'Exploratory GO profile\nNo terms passed BH < 0.05; shown BH P = {shown_bh}',fontsize=5.4,pad=3)
        counts=sorted({int(g.Count.min()),int(g.Count.max())})
        handles=[ax.scatter([],[],s=10+4*n,c='#c4c9ce',edgecolors='#333333',linewidths=.3,label=str(n)) for n in counts]
        ax.legend(handles=handles,title='Genes',loc='lower right',frameon=False,fontsize=4.5,title_fontsize=4.5,handletextpad=.4)
        fig.colorbar(sc,ax=ax,fraction=.04,pad=.03,label='-log10 raw P')
    else:
        g=go.head(10).sort_values('RichFactor');y=np.arange(len(g))
        sc=ax.scatter(g.RichFactor,y,c=-np.log10(g.AdjPValue_BH),s=10+3*g.Count,cmap='RdBu_r',edgecolors='#333333',linewidths=.3)
        import textwrap
        ax.set_yticks(y,['\n'.join(textwrap.wrap(t,25)) for t in g.TERM],fontsize=5.5)
        ax.set_xlabel('Rich factor');fig.colorbar(sc,ax=ax,fraction=.04,pad=.02,label='−log10 BH P')
    groups['Figure3B_BLA_GO']=fig.axes[start:];start=len(fig.axes)
    ax=fig.add_subplot(gs[1,0]);panel(ax,'C')
    for region,g in scores.groupby('Region',sort=False):
        add_confidence_ellipse(ax,g.PC1.to_numpy(),g.PC2.to_numpy(),REGION_PALETTE[region])
        ax.scatter(g.PC1,g.PC2,s=3.5,marker='o',color=REGION_PALETTE[region],edgecolors='white',linewidths=.2,label=region,alpha=.88,zorder=3)
    square_pca_axes(ax)
    ax.set_title('95% data ellipses',fontsize=6,pad=3)
    ax.set_xlabel(f'PC1 ({variance.ExplainedVariancePercent.iloc[0]:.2f}%)');ax.set_ylabel(f'PC2 ({variance.ExplainedVariancePercent.iloc[1]:.2f}%)')
    ax.legend(title='Brain region',fontsize=5.5,title_fontsize=5.8,ncol=1,
              loc='upper left',bbox_to_anchor=(1.04,1),frameon=False,
              borderaxespad=0,handletextpad=.3,labelspacing=.4)
    groups['Figure3C_PCA']=fig.axes[start:];start=len(fig.axes)
    ax=fig.add_subplot(gs[1,1]);panel(ax,'D')
    im=ax.imshow(z.to_numpy(),aspect='auto',vmin=-2.5,vmax=2.5,cmap='RdBu_r',interpolation='nearest')
    ax.set_xticks(np.arange(z.shape[1]),z.columns,rotation=90,fontsize=5.5)
    labels=info.reindex(z.index).Genes.fillna(pd.Series(z.index,index=z.index)).str.split(';').str[0]
    ax.set_yticks(np.arange(len(z)),labels,fontsize=5.4);ax.yaxis.tick_right();ax.tick_params(length=0)
    cax=ax.inset_axes([-.16,.15,.025,.7]);fig.colorbar(im,cax=cax,label='Row Z-score');cax.yaxis.set_ticks_position('left');cax.yaxis.set_label_position('left')
    groups['Figure3D_heatmap']=fig.axes[start:];start=len(fig.axes)
    ax=fig.add_subplot(gs[2,:]);panel(ax,'E')
    im=ax.imshow(corr.to_numpy(),cmap='RdBu_r',vmin=-1,vmax=1,aspect='auto')
    ax.set_xticks(range(len(corr)),corr.columns,rotation=90,fontsize=6)
    ax.set_yticks(range(len(corr)),corr.index,fontsize=6)
    for i in range(len(corr)):
        for j in range(len(corr)):ax.text(j,i,f'{corr.iloc[i,j]:.2f}',ha='center',va='center',fontsize=4.5,color='white' if abs(corr.iloc[i,j])>.65 else 'black')
    fig.colorbar(im,ax=ax,fraction=.018,pad=.02,label='Pearson r')
    groups['Figure3E_region_correlation']=fig.axes[start:]
    fig.subplots_adjust(left=.10,right=.86,bottom=.06,top=.95)
    save_panels(fig, groups, args.output_dir)
    plt.close(fig)



def assemble_figure2(args):
    """Reference-style independent QC panels; unchanged numerical calculations."""
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    from matplotlib.ticker import FuncFormatter
    from matplotlib.colors import LinearSegmentedColormap,Normalize,to_rgba
    from matplotlib.cm import ScalarMappable
    reports=qa.load_all_reports(args.data_dir)
    matrix,meta,annotation=qa.build_combined_matrix(reports)
    sq=qa.sample_qc_summary_from_reports(reports)
    cv=qa.qc_cv_table(matrix,meta,annotation);cv=cv.loc[cv.qc_detected_n.ge(20)].copy()
    complete=qa.prepare_complete_matrix(matrix);correlation=complete.corr()
    order=meta.assign(rank=meta.region.map({v:i for i,v in enumerate(qa.HEATMAP_REGION_ORDER)})).sort_values(['rank','group','replicate']).sample_id
    correlation=correlation.loc[order,order]
    legend_order=['QC','PL','IL','SSC','MSC','BLA','CEA','DHPC','VHPC','PVH','PMD','LPB']
    handles=[Patch(color=qa.REGION_COLORS[r],label=r) for r in legend_order]

    fig=plt.figure(figsize=(9.6,4.5))
    detected=fig.add_axes([.08,.54,.82,.39])
    missing=fig.add_axes([.08,.10,.82,.39],sharex=detected)
    x=np.arange(len(sq));colors=[qa.REGION_COLORS[r] for r in sq.region]
    maximum=10000
    missing_max=40
    detected.bar(x,sq.detected_protein_groups,color=colors,width=.82,linewidth=0,antialiased=False,rasterized=True,zorder=3)
    detected.set(xlim=(-.5,len(sq)-.5),ylim=(0,maximum),ylabel='Detected')
    detected.set_yticks(np.arange(0,maximum+1,2500))
    detected.yaxis.set_major_formatter(FuncFormatter(lambda value,position:f'{value:.0f}'))
    missing.scatter(x,sq.missing_percent,c=colors,s=4,linewidths=0,zorder=3)
    missing.set(ylim=(0,missing_max),ylabel='Missing')
    missing_ticks=np.arange(0,missing_max+1,10)
    missing.set_yticks(missing_ticks)
    missing.yaxis.set_major_formatter(FuncFormatter(lambda value,position:f'{value:.0f}%'))
    # Retrospective primary exclusion threshold: more than 15% missingness.
    missing.axhline(15,color='#b2182b',ls=(0,(5,5)),lw=1,zorder=2)
    missing.text(1.006,15,'15%',color='#b2182b',transform=missing.get_yaxis_transform(),va='center',fontsize=8)
    boundaries=np.flatnonzero(sq.region.to_numpy()[1:]!=sq.region.to_numpy()[:-1])+.5
    for ax in (detected,missing):
        ax.set_xticks([]);ax.spines[['top','right']].set_visible(False)
        ax.spines['bottom'].set_color('#dddddd');ax.spines['bottom'].set_linewidth(.8)
        ax.spines['left'].set_linewidth(.8);ax.tick_params(axis='y',labelsize=8,length=3,width=.7)
        ax.set_ylabel(ax.get_ylabel(),fontsize=9)
        for boundary in boundaries:ax.axvline(boundary,color='#999999',lw=.8,ls=(0,(4,5)),zorder=1)
    missing.set_axisbelow(True);missing.yaxis.grid(color='#dddddd',lw=.8)
    for boundary in boundaries:
        fraction=.08+.82*(boundary+.5)/len(sq)
        fig.add_artist(Line2D([fraction,fraction],[.49,.54],transform=fig.transFigure,color='#999999',lw=.8,ls=(0,(4,5))))
    outlier=int(sq.missing_percent.argmax())
    detected.annotate('',xy=(outlier,maximum*1.012),xytext=(outlier,maximum*1.10),annotation_clip=False,
                      arrowprops={'arrowstyle':'->','color':'#e84b35','lw':1.2})
    detected.legend(handles=handles,loc='upper left',bbox_to_anchor=(1.022,.91),frameon=False,
                    fontsize=8,handlelength=.7,handletextpad=.4,borderaxespad=0,labelspacing=.25)
    fig.text(.008,.995,'D',fontsize=15,va='top')
    save_formats(fig,args.output_dir/'Figure2D_detection_missing');plt.close(fig)

    fig=plt.figure(figsize=(4.8,4.3));ax=fig.add_axes([.14,.14,.73,.76])
    low=cv.qc_cv_percent.le(100)
    ax.scatter(cv.loc[low,'qc_mean_log2_intensity'],cv.loc[low,'qc_cv_percent'],s=.30,c='#2c7fb8',
               alpha=.30,linewidths=0,rasterized=True)
    ax.scatter(cv.loc[~low,'qc_mean_log2_intensity'],np.full(int((~low).sum()),100),s=.30,c='#2c7fb8',
               alpha=.12,linewidths=0,rasterized=True)
    xmin=float(cv.qc_mean_log2_intensity.min());xmax=float(cv.qc_mean_log2_intensity.max())
    ax.set(xlim=(xmin,xmax),ylim=(0,100),xlabel='Mean log2 intensity in QC samples',ylabel='QC CV (%)')
    ax.set_xticks(np.linspace(xmin,xmax,6));ax.xaxis.set_major_formatter(FuncFormatter(lambda value,position:f'{value:.1f}'))
    ax.set_yticks(np.arange(0,101,20));ax.yaxis.set_major_formatter(FuncFormatter(lambda value,position:f'{value:.0f}%'))
    ax.set_axisbelow(True);ax.yaxis.grid(color='#dddddd',lw=.8)
    ax.spines[['top','right']].set_visible(False)
    for threshold in (20,30):
        ax.axhline(threshold,color='#b2182b',ls=(0,(5,5)),lw=1)
        ax.text(1.01,threshold,f'CV {threshold}%',transform=ax.get_yaxis_transform(),color='#b2182b',fontsize=7.5,va='center')
    ax.tick_params(labelsize=7.5,length=3);ax.set_xlabel(ax.get_xlabel(),fontsize=9);ax.set_ylabel('QC CV (%)',fontsize=9)
    fig.text(.012,.99,'E',fontsize=15,va='top')
    save_formats(fig,args.output_dir/'Figure2E_QC_CV');plt.close(fig)

    # Keep each correlation cell constant (nearest-neighbor display).
    # Bake 0.035 pt white separators into the 1200 dpi data layer so PDF
    # viewers cannot force each hairline to a full screen pixel.
    # Exact original color formula: r=.85 -> #deebf7; r=1 -> #084594.
    fig=plt.figure(figsize=(5,4.4));ax=fig.add_axes([.055,.04,.77,.875])
    values=correlation.to_numpy();fraction=np.clip((values-.85)/.15,0,1)
    light=np.array([222,235,247],dtype=float);dark=np.array([8,69,148],dtype=float)
    rgb=(light[None,None,:]*(1-fraction[:,:,None])+dark[None,None,:]*fraction[:,:,None]).astype(np.uint8)
    ax.imshow(rgb,origin='upper',interpolation='nearest',resample=False,aspect='equal')
    from matplotlib.collections import LineCollection
    n=len(correlation);boundaries=np.arange(.5,n-.5,1)
    segments=[[(v,-.5),(v,n-.5)] for v in boundaries]+[[(-.5,v),(n-.5,v)] for v in boundaries]
    ax.add_collection(LineCollection(segments,colors='white',linewidths=.035,antialiaseds=True,capstyle='butt',zorder=4))
    ax.set_rasterization_zorder(5)
    ax.set_xticks([]);ax.set_yticks([])
    for spine in ax.spines.values():spine.set_linewidth(.45);spine.set_color('#333333')
    region_by_id=meta.set_index('sample_id').region
    band=np.array([to_rgba(qa.REGION_COLORS[region_by_id[sample]]) for sample in correlation.columns])
    top=fig.add_axes([.055,.919,.77,.0105]);top.imshow(band[np.newaxis,:,:],aspect='auto',interpolation='nearest',resample=False);top.axis('off')
    left=fig.add_axes([.0405,.04,.0092,.875]);left.imshow(band[:,np.newaxis,:],aspect='auto',interpolation='nearest',resample=False);left.axis('off')
    top.add_collection(LineCollection([[(v,-.5),(v,.5)] for v in boundaries],colors='white',linewidths=.035,antialiaseds=True,capstyle='butt',zorder=4))
    left.add_collection(LineCollection([[(-.5,v),(.5,v)] for v in boundaries],colors='white',linewidths=.035,antialiaseds=True,capstyle='butt',zorder=4))
    top.set_rasterization_zorder(5);left.set_rasterization_zorder(5)
    cmap=LinearSegmentedColormap.from_list('reference_blues',['#deebf7','#084594'])
    cax=fig.add_axes([.858,.647,.022,.245]);bar=fig.colorbar(ScalarMappable(norm=Normalize(.85,1),cmap=cmap),cax=cax,ticks=[.85,1])
    bar.ax.set_yticklabels(['0.85','1.00']);bar.ax.tick_params(length=0,labelsize=8,pad=3);bar.outline.set_linewidth(.4)
    cax.set_title('r',loc='left',fontsize=9,pad=3)
    fig.legend(handles=handles,loc='upper left',bbox_to_anchor=(.850,.605),frameon=False,
               fontsize=8,handlelength=.7,handletextpad=.4,borderaxespad=0,labelspacing=.25)
    fig.text(.007,.99,'F',fontsize=15,va='top')
    save_formats(fig,args.output_dir/'Figure2F_sample_correlation',dpi=1200);plt.close(fig)


if __name__=='__main__':main()
