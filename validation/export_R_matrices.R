# Run from any directory: Rscript export_R_matrices.R <input CSV directory>.
# Set R_LIBS_USER before starting R if dependencies are in a custom library.
args <- commandArgs(trailingOnly=TRUE)
script <- sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE)[1])
repo <- normalizePath(file.path(dirname(script),'..'),winslash='/',mustWork=TRUE)
data_dir <- if(length(args)) normalizePath(args[1],winslash='/',mustWork=TRUE) else file.path(repo,'..','data')
source(file.path(repo,'sample_analysis/scripts/R/common.R'))
output_root <- if(length(args)>=2L) normalizePath(args[2],winslash='/',mustWork=TRUE) else repo
parse_cli_options <- function(...) list(exclusions=file.path(repo,'sample_analysis/sample_exclusions.csv'),output_dir=file.path(output_root,'sample_analysis/results'))
out <- file.path(output_root,'validation/R_matrices')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
for(region in REGION_CODES) {
  p <- preprocess_region_matrix(region,data_dir)
  write.csv(data.frame(ProteinAccessions=rownames(p$normalized),p$normalized,check.names=FALSE),file.path(out,paste0(region,'.csv')),row.names=FALSE)
}
for(f in list.files(file.path(repo,'sample_analysis/scripts/R'),pattern='\\.R$',full.names=TRUE)) parse(f)
cat('All production R sources parse; eleven validation matrices exported.\n')
