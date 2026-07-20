script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_file <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[1]) else NA_character_
script_dir <- if (!is.na(script_file) && file.exists(script_file)) dirname(normalizePath(script_file, winslash = "/", mustWork = TRUE)) else getwd()
root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)

packages <- c("limma", "AnnotationDbi", "org.Mm.eg.db", "GO.db")
invisible(lapply(packages, function(pkg) {
  suppressPackageStartupMessages(require(pkg, character.only = TRUE))
}))

dir.create(file.path(root, "environment"), recursive = TRUE, showWarnings = FALSE)
capture.output(sessionInfo(), file = file.path(root, "environment", "sessionInfo.txt"))
message(file.path(root, "environment", "sessionInfo.txt"))
