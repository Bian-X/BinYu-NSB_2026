#!/usr/bin/env python3
"""Generate and validate Figure 2D-F and Figure 3 from proteomics reports."""
from __future__ import annotations
import argparse,csv,datetime,hashlib,importlib.metadata,json,os,platform,subprocess,sys
from pathlib import Path
from validation.check_inputs import report_path
ROOT=Path(__file__).resolve().parent
PYTHON_PACKAGES={'numpy':'2.5.4','pandas':'3.0.1','scipy':'1.18.1','scikit-learn':'1.9.1','matplotlib':'3.11.2','openpyxl':'3.1.5','reportlab':'4.4.9','Pillow':'12.3.0'}
R_PACKAGES={'limma':'3.68.5','AnnotationDbi':'1.74.0','org.Mm.eg.db':'3.23.0','GO.db':'3.23.1','statmod':'1.5.2'}
def sha256(path):
    with Path(path).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()


def snapshot_outputs(folders):
    """Record existing files so a reused directory cannot publish stale results."""
    existing = {}
    for folder in folders:
        for path in folder.rglob('*'):
            if path.is_file():
                info = path.stat()
                existing[path] = (info.st_mtime_ns, info.st_size)
    return existing


def collect_current_outputs(output_root, folders, before):
    """Hash only files created or rewritten by the active analysis branches."""
    after = snapshot_outputs(folders)
    return [{'path': path.relative_to(output_root).as_posix(),
             'bytes': signature[1], 'sha256': sha256(path)}
            for path, signature in sorted(after.items())
            if before.get(path) != signature]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--data-dir',type=Path,required=True)
    p.add_argument('--metadata',type=Path,default=ROOT/'metadata/Sample_metadata.xlsx')
    p.add_argument('--rscript',default='Rscript');p.add_argument('--r-library',type=Path)
    p.add_argument('--output-root',type=Path,default=ROOT)
    p.add_argument('--figures',choices=['all','figure2','figure3'],default='all')
    p.add_argument('--environment',choices=['manuscript','compatible'],default='manuscript')
    p.add_argument('--font-family',default='Arial');p.add_argument('--selection-sensitivity',action='store_true')
    p.add_argument('--events',type=Path);p.add_argument('--rt-metrics',type=Path);p.add_argument('--rt-label',default='Provided retention-time metric')
    a=p.parse_args();a.data_dir=a.data_dir.resolve(strict=True);a.output_root=a.output_root.resolve();a.output_root.mkdir(parents=True,exist_ok=True)
    validation=a.output_root/'validation';validation.mkdir(exist_ok=True)
    env=os.environ.copy();env.update(PYTHONDONTWRITEBYTECODE='1',PYTHONUTF8='1')
    if platform.system()=='Windows':env['LC_ALL']='English_United States.utf8'
    if a.r_library:env['R_LIBS_USER']=str(a.r_library.resolve(strict=True))
    run={'package_version':'2026.10.26','scope':a.figures,'selection_sensitivity':a.selection_sensitivity,'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'python':platform.python_version(),'platform':platform.platform(),'environment_policy':a.environment,'commands':[],'status':'running','behavioral_analysis_included':False}
    def record():(validation/'run_manifest.json').write_text(json.dumps(run,indent=2),encoding='utf-8')
    def execute(executable,*arguments):
        command=[str(executable),*[str(x) for x in arguments]];print('+ '+subprocess.list2cmdline(command),flush=True)
        result=subprocess.run(command,cwd=ROOT,env=env);run['commands'].append({'argv':command,'returncode':result.returncode});record()
        if result.returncode:raise RuntimeError(f'Command failed with exit code {result.returncode}: {command[1]}')
    def py(path,*args):execute(sys.executable,ROOT/path,*args)
    def r(path,*args):execute(a.rscript,ROOT/path,*args)
    record()
    try:
        py('validation/check_inputs.py','--data-dir',a.data_dir)
        with (ROOT/'metadata/input_manifest.csv').open(encoding='utf-8',newline='') as handle:
            expected_reports=list(csv.DictReader(handle))
        report_files=[report_path(a.data_dir,row['filename']) for row in expected_reports]
        run['input_reports']=[{'filename':f.name,'sha256':sha256(f),'bytes':f.stat().st_size} for f in report_files]
        run['python_packages']={k:importlib.metadata.version(k) for k in PYTHON_PACKAGES}
        differences=[f'{k}: expected {v}, found {run["python_packages"][k]}' for k,v in PYTHON_PACKAGES.items() if run['python_packages'][k]!=v]
        if sys.version_info[:2]!=(3,12):differences.append(f'Python: expected 3.12.x (validated reference 3.12.14), found {platform.python_version()}')
        run['validated_reference_python']='3.12.14'
        run['python_patch_differs_from_reference']=platform.python_version()!='3.12.14'
        run['font_family']=a.font_family
        if a.figures in ('all','figure3'):
            probe=subprocess.run([a.rscript,str(ROOT/'validation/capture_R_environment.R')],cwd=ROOT,env=env,capture_output=True,text=True,encoding='utf-8',errors='replace',check=True)
            info=dict(line.split('\t',1) for line in probe.stdout.splitlines() if '\t' in line);run['R_environment']=info
            if info.get('R')!='4.6.1':differences.append(f'R: expected 4.6.1, found {info.get("R")}')
            differences += [f'{k}: expected {v}, found {info.get(k)}' for k,v in R_PACKAGES.items() if info.get(k)!=v]
            if a.selection_sensitivity and info.get('Matrix')!='1.7.5':differences.append(f'Matrix: expected 1.7.5 for --selection-sensitivity, found {info.get("Matrix")}')
        run['environment_differences']=differences;record()
        if differences:
            message='Environment differs from manuscript:\n'+'\n'.join(differences)
            if a.environment=='manuscript':raise RuntimeError(message+'\nInstall the recorded versions or explicitly select --environment compatible.')
            print('WARNING: '+message,flush=True)
        py('validation/test_regressions.py')
        if a.figures in ('all','figure3'):r('validation/test_regressions.R')
        result=a.output_root/'sample_analysis/results';figures=a.output_root/'sample_analysis/figures';qc=a.output_root/'QC_analysis/outputs';assembled=a.output_root/'manuscript_figures'
        for folder in (result,figures,qc,assembled):folder.mkdir(parents=True,exist_ok=True)
        active_folders=[assembled]
        if a.figures in ('all','figure3'):active_folders.extend((result,figures))
        if a.figures in ('all','figure2'):active_folders.append(qc)
        output_before=snapshot_outputs(active_folders)
        if a.figures in ('all','figure3'):
            common=['--data-dir',a.data_dir,'--output-dir',result,'--font-family',a.font_family]
            r('sample_analysis/scripts/R/volcano_PL_G5_vs_G1.R',*common,'--output-pdf',figures/'PL_volcano_G5_vs_G1.pdf')
            r('sample_analysis/scripts/R/go_enrichment_PL_G5_vs_G1.R',*common,'--output-pdf',figures/'PL_GO_enrichment_G5_vs_G1.pdf')
            expanded=result/'expanded_mapping_sensitivity';expanded.mkdir(exist_ok=True)
            r('sample_analysis/scripts/R/go_enrichment_PL_G5_vs_G1.R','--data-dir',a.data_dir,'--output-dir',expanded,'--font-family',a.font_family,'--go-mapping','expanded','--go-pvalue','adjusted','--output-pdf',expanded/'GO_PL_expanded_adjusted.pdf')
            for script in ('go_mapping_sensitivity_plot.py','go_ranked_profile.py'):py('sample_analysis/scripts/python/'+script,'--results-dir',result,'--figure-dir',figures)
            bla=result/'BLA_G4_vs_G1_unique';bla.mkdir(exist_ok=True)
            r('sample_analysis/scripts/R/go_enrichment_BLA_G4_vs_G1.R','--data-dir',a.data_dir,'--output-dir',bla,'--font-family',a.font_family,'--go-mapping','unique','--go-pvalue','adjusted','--output-pdf',figures/'BLA_GO_enrichment_G4_vs_G1.pdf')
            py('sample_analysis/scripts/python/brain_region_pca.py','--data-dir',a.data_dir,'--output-dir',result/'pca','--figure-dir',figures,'--font-family',a.font_family)
            r('sample_analysis/scripts/R/region_heatmap_G4_vs_G2.R',*common,'--output-pdf',figures/'multi_region_G2_vs_G4_heatmap.pdf')
            r('sample_analysis/scripts/R/region_log2fc_correlation_G5_G1.R',*common,'--output-pdf',figures/'brain_region_log2FC_correlation_G5_G1.pdf','--min-pairwise-proteins','30')
            r('sample_analysis/environment/capture_session_info.R',result/'provenance')
            if a.selection_sensitivity:
                screen=result/'workflow_example_comparison'
                r('sample_analysis/scripts/R/compare_workflow_examples.R','--data-dir',a.data_dir,'--output-dir',screen,'--reference-dir',result)
                r('sample_analysis/scripts/R/validate_workflow_shortlist.R','--data-dir',a.data_dir,'--screen-dir',screen,'--output-dir',result,'--targets','BLA:G1:G4')
            r('validation/export_R_matrices.R',a.data_dir,a.output_root)
        if a.figures in ('all','figure2'):
            a.metadata=a.metadata.resolve(strict=True);run['metadata_sha256']=sha256(a.metadata)
            for script in ('run_qc_pdf.py','figure_2F_injection_order.py','qc_drift_curve.py','qc_correlation_heatmap.py','qc_time_order_summary.py'):
                extra=[] if script=='run_qc_pdf.py' else ['--metadata',a.metadata]
                if script=='qc_time_order_summary.py':
                    if a.events:extra+=['--events',a.events.resolve(strict=True)]
                    if a.rt_metrics:extra+=['--rt-metrics',a.rt_metrics.resolve(strict=True),'--rt-label',a.rt_label]
                py('QC_analysis/scripts/'+script,'--data-dir',a.data_dir,'--out-dir',qc,*extra)
        py('sample_analysis/scripts/python/assemble_manuscript_figures.py','--data-dir',a.data_dir,'--results-dir',result,'--output-dir',assembled,'--figure',a.figures,'--font-family',a.font_family)
        if a.figures=='all':py('validation/verify_results.py','--data-dir',a.data_dir,'--output-root',a.output_root)
        check_args=['--output-root',a.output_root,'--figures',a.figures,'--environment',a.environment]
        if a.selection_sensitivity:check_args.append('--selection-sensitivity')
        py('validation/check_manuscript_results.py',*check_args)
        run['status']='passed';run['finished_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat()
        run['output_files']=collect_current_outputs(a.output_root,active_folders,output_before)
        record();print('Completed '+a.figures+'. Run provenance: '+str(validation),flush=True)
    except Exception as exc:run['status']='failed';run['error']=str(exc);record();raise
if __name__=='__main__':main()
