script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_file <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[1]) else NA_character_
script_dir <- if (!is.na(script_file) && file.exists(script_file)) dirname(normalizePath(script_file, winslash = "/", mustWork = TRUE)) else getwd()
root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)

packages <- c("limma", "AnnotationDbi", "org.Mm.eg.db", "GO.db", "statmod")
invisible(lapply(packages, function(pkg) {
  suppressPackageStartupMessages(require(pkg, character.only = TRUE))
}))

args <- commandArgs(trailingOnly=TRUE)
output_dir <- if(length(args)) args[1] else file.path(root, "results", "provenance")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
capture.output(sessionInfo(), file = file.path(output_dir, "R_sessionInfo.txt"))
message(file.path(output_dir, "R_sessionInfo.txt"))
