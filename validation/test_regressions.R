# Base-R regression tests: run without study data or Bioconductor packages.
script <- sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value=TRUE)[1])
repo <- normalizePath(file.path(dirname(script), '..'), winslash='/', mustWork=TRUE)
if (.Platform$OS.type == 'windows') {
  Sys.setlocale('LC_CTYPE', 'English_United States.utf8')
}
cases <- 0L
for (module in c('common.R', 'common_two_group.R')) {
  test_env <- new.env(parent=globalenv())
  source(file.path(repo, 'sample_analysis/scripts/R', module), local=test_env)
  test_env$capabilities <- function(...) FALSE
  for (requested in c('Arial', 'Courier')) {
    path <- tempfile(fileext='.pdf')
    messages <- character()
    withCallingHandlers(test_env$open_figure_pdf(path, width=2, height=2, family=requested),
      warning=function(w) {
        messages <<- c(messages, conditionMessage(w))
        invokeRestart('muffleWarning')
      })
    expected <- if (requested == 'Arial') 'Helvetica' else 'Courier'
    stopifnot(any(grepl(paste0('font family: ', expected), messages, fixed=TRUE)))
    plot.new()
    text(.5, .5, 'PDF fallback')
    dev.off()
    stopifnot(identical(rawToChar(readBin(path, 'raw', 4)), '%PDF'))
    unlink(path)
    cases <- cases + 1L
  }
}
cat('Base-R PDF font regression tests passed:', cases, 'cases\n')
