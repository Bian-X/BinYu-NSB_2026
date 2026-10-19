# Tested release environment: R 4.6.1 / Bioconductor 3.23.
if(as.character(getRversion())!='4.6.1') stop('Use R 4.6.1 for the validated release environment')
if(.Platform$OS.type=='windows') {
  Sys.setenv(LC_ALL='English_United States.utf8',LANG='English_United States.utf8')
  Sys.setlocale('LC_CTYPE','English_United States.utf8')
}
target <- Sys.getenv('R_LIBS_USER',unset=file.path(getwd(),'.R_library'))
dir.create(target,recursive=TRUE,showWarnings=FALSE)
target <- normalizePath(target,winslash='/',mustWork=TRUE)
.libPaths(c(target,.libPaths()));Sys.setenv(R_LIBS_USER=target)
args <- commandArgs(trailingOnly=TRUE)
if(any(!args %in% '--check-only')) stop('Supported option: --check-only')
expected <- c(limma='3.68.5',AnnotationDbi='1.74.0',org.Mm.eg.db='3.23.0',GO.db='3.23.1',statmod='1.5.2',Matrix='1.7.5')
versions <- function() vapply(names(expected),function(p) {
  # Read installed metadata without loading DLLs that Windows cannot replace.
  # read.dcf also handles non-ASCII library paths before namespace loading.
  tryCatch(as.character(package_version(read.dcf(file.path(find.package(p),'DESCRIPTION'),fields='Version')[1,1])),
           error=function(e) 'MISSING')
},character(1))
actual <- versions()
if(!identical(actual,expected) && !('--check-only' %in% args)) {
  options(repos=c(CRAN='https://cloud.r-project.org'),timeout=600)
  if(!requireNamespace('BiocManager',quietly=TRUE)) install.packages('BiocManager',lib=target)
  BiocManager::install(version='3.23',lib=target,ask=FALSE,update=FALSE)
  # Matrix is used by the optional exhaustive comparison and complete-case scripts.
  needed <- names(expected)[actual!=expected]
  BiocManager::install(needed,lib=target,ask=FALSE,update=FALSE,force=TRUE)
  actual <- versions()
}
if(!identical(actual,expected)) {
  bad <- names(expected)[actual!=expected]
  stop(paste0('R dependencies differ from the validated release: ',
    paste(paste0(bad,': expected ',expected[bad],', found ',actual[bad]),collapse='; '),
    '. Matching archived packages may be needed if upstream versions have changed.'))
}
loadable <- vapply(names(expected),requireNamespace,logical(1),quietly=TRUE)
if(!all(loadable)) stop('R dependency namespaces failed to load: ',paste(names(expected)[!loadable],collapse=', '))
message('Validated release R dependency versions verified in the active library paths.')
