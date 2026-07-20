required_bioc <- c("limma", "AnnotationDbi", "org.Mm.eg.db", "GO.db")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

BiocManager::install(version = "3.20", ask = FALSE, update = FALSE)

missing <- required_bioc[!vapply(required_bioc, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  BiocManager::install(missing, ask = FALSE, update = FALSE)
}

message("R package check completed.")
