param(
    [string]$DataDir = "",
    [string]$Python = "python",
    [string]$Rscript = "Rscript",
    [string]$FontFamily = "Arial"
)

$ErrorActionPreference = "Stop"
$Root = (Resolve-Path $PSScriptRoot).Path
$FigureDir = Join-Path $Root "figures"
$ResultDir = Join-Path $Root "results"

function Test-RegionData {
    param([string]$Path)
    if (-not (Test-Path -LiteralPath $Path)) {
        return $false
    }
    $matches = Get-ChildItem -LiteralPath $Path -Filter "*.csv" |
        Where-Object { $_.Name -match "_(PMD|LPB|VHPC|PVH|BLA|DHPC|CEA|MSC|SSC|IL|PL)_DIA_LFQ.*\.csv$" }
    return ($matches.Count -ge 11)
}

if ([string]::IsNullOrWhiteSpace($DataDir)) {
    $candidates = @(
        (Join-Path $Root "data"),
        (Join-Path $Root "data\data1"),
        (Join-Path $Root "..\data\data1"),
        (Join-Path $Root "..\data")
    )
    foreach ($candidate in $candidates) {
        if (Test-RegionData $candidate) {
            $DataDir = (Resolve-Path $candidate).Path
            break
        }
    }
}
if ([string]::IsNullOrWhiteSpace($DataDir) -or -not (Test-RegionData $DataDir)) {
    throw "Could not find the 11 brain-region CSV files. Pass -DataDir <folder>."
}

$LocalRRoot = Join-Path $Root "..\.tools\r-base"
if (Test-Path -LiteralPath $LocalRRoot) {
    $env:PATH = "$LocalRRoot;$(Join-Path $LocalRRoot 'Library\bin');$(Join-Path $LocalRRoot 'Scripts');$env:PATH"
    $LocalRscript = Join-Path $LocalRRoot "Scripts\Rscript.exe"
    if ($Rscript -eq "Rscript" -and (Test-Path -LiteralPath $LocalRscript)) {
        $Rscript = $LocalRscript
    }
}
$env:PROTEOMICS_FIGURE_FONT = $FontFamily

New-Item -ItemType Directory -Force -Path $FigureDir, $ResultDir | Out-Null

function Invoke-Checked {
    param(
        [Parameter(Mandatory = $true)][string]$Executable,
        [Parameter(Mandatory = $true)][string[]]$Arguments
    )
    & $Executable @Arguments
    if ($LASTEXITCODE -ne 0) {
        throw "Command failed with exit code ${LASTEXITCODE}: $Executable $($Arguments -join ' ')"
    }
}

$RCommonArgs = @("--data-dir", $DataDir, "--output-dir", $ResultDir, "--font-family", $FontFamily)

Invoke-Checked $Rscript (@((Join-Path $Root "scripts\R\volcano_BLA_G4_vs_G2.R")) + $RCommonArgs + @("--output-pdf", (Join-Path $FigureDir "BLA_volcano_G4_vs_G2.pdf")))
Invoke-Checked $Rscript (@((Join-Path $Root "scripts\R\go_enrichment_BLA_G4_vs_G2.R")) + $RCommonArgs + @("--output-pdf", (Join-Path $FigureDir "BLA_GO_enrichment_G4_vs_G2.pdf")))
Invoke-Checked $Python @(
    (Join-Path $Root "scripts\python\brain_region_pca.py"),
    "--data-dir", $DataDir,
    "--output-dir", (Join-Path $ResultDir "pca"),
    "--output-pdf", (Join-Path $FigureDir "brain_region_pca.pdf"),
    "--font-family", $FontFamily
)
Invoke-Checked $Rscript (@((Join-Path $Root "scripts\R\region_paired_heatmap_G2_G4.R")) + $RCommonArgs + @("--output-pdf", (Join-Path $FigureDir "multi_region_G2_vs_G4_heatmap.pdf")))
Invoke-Checked $Rscript (@((Join-Path $Root "scripts\R\region_log2fc_correlation_G2_G4.R")) + $RCommonArgs + @("--output-pdf", (Join-Path $FigureDir "brain_region_log2FC_correlation.pdf"), "--min-pairwise-proteins", "30"))
Invoke-Checked $Rscript @((Join-Path $Root "environment\capture_session_info.R"))

$required = @(
    "BLA_volcano_G4_vs_G2.pdf",
    "BLA_GO_enrichment_G4_vs_G2.pdf",
    "brain_region_pca.pdf",
    "multi_region_G2_vs_G4_heatmap.pdf",
    "brain_region_log2FC_correlation.pdf"
)
foreach ($name in $required) {
    $path = Join-Path $FigureDir $name
    if (-not (Test-Path -LiteralPath $path)) {
        throw "Missing expected figure: $path"
    }
}

Write-Host "Figure generation completed."
Write-Host "Data directory: $DataDir"
Write-Host "Figures: $FigureDir"
Write-Host "Results: $ResultDir"
