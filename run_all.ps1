param(
 [Parameter(Mandatory=$true)][string]$DataDir,
 [string]$Metadata='', [string]$Python='python', [string]$Rscript='Rscript',
 [string]$RLibrary='', [string]$OutputRoot='',
 [ValidateSet('all','figure2','figure3')][string]$Figures='all',
 [ValidateSet('manuscript','compatible')][string]$Environment='manuscript',
 [switch]$SelectionSensitivity, [string]$Events='', [string]$RtMetrics='',
 [string]$RtLabel='Provided retention-time metric', [string]$FontFamily='Arial'
)
$ErrorActionPreference='Stop'
$arguments=@((Join-Path $PSScriptRoot 'run_all.py'),'--data-dir',$DataDir,'--rscript',$Rscript,'--figures',$Figures,'--environment',$Environment,'--font-family',$FontFamily)
if($Metadata){$arguments+=@('--metadata',$Metadata)}
if($RLibrary){$arguments+=@('--r-library',$RLibrary)}
if($OutputRoot){$arguments+=@('--output-root',$OutputRoot)}
if($SelectionSensitivity){$arguments+='--selection-sensitivity'}
if($Events){$arguments+=@('--events',$Events)}
if($RtMetrics){$arguments+=@('--rt-metrics',$RtMetrics,'--rt-label',$RtLabel)}
& $Python @arguments
if($LASTEXITCODE -ne 0){throw "Proteomics workflow failed (exit $LASTEXITCODE)."}
