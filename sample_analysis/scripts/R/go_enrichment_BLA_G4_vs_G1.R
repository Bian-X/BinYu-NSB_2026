#!/usr/bin/env Rscript
# Selected GO workflow example: BLA G4-G1, six-group model and unique Entrez mapping.
# Kept separate from the PL G5-G1 volcano example.

options(stringsAsFactors = FALSE)
rm(list = ls())

suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(AnnotationDbi))
suppressPackageStartupMessages(library(org.Mm.eg.db))
has_go_db <- requireNamespace("GO.db", quietly = TRUE)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_file <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[1]) else NA_character_
script_dir <- if (!is.na(script_file) && file.exists(script_file)) dirname(normalizePath(script_file, winslash = "/", mustWork = TRUE)) else getwd()
source(file.path(script_dir, "common.R"))

paths <- analysis_paths("BLA_GO_enrichment_G4_vs_G1.pdf")
data_dir <- paths$data_dir
output_dir <- paths$output_dir
output_pdf <- paths$output_pdf
figure_font <- paths$font_family

min_go_size <- 5L
max_go_size <- 500L
# BH family: all eligible 5-500-background-gene terms pooled across BP/CC/MF.
# Zero-overlap terms are tested with P=1 before candidate selection.
options_cli <- parse_cli_options()
go_pvalue <- as.character(option_value(options_cli, "go_pvalue", "adjusted"))
go_cutoff <- as.numeric(option_value(options_cli, "go_cutoff", 0.05))
if (!go_pvalue %in% c("raw", "adjusted") || !is.finite(go_cutoff) || go_cutoff <= 0 || go_cutoff > 1) {
  stop("Use --go-pvalue raw|adjusted and --go-cutoff in (0,1]")
}
go_field <- if (go_pvalue == "adjusted") "AdjPValue_BH" else "PValue"
top_n <- 10L

analysis <- run_region_differential("BLA", "BLA", data_dir, group_a = "G1", group_b = "G4")
result <- analysis$result
result$Status[result$Status == "High in G1"] <- "Low in G4"
write.csv(
  result,
  file = file.path(output_dir, "BLA_G4_vs_G1_differential_results_BH_FC1p5.csv"),
  row.names = FALSE,
  na = ""
)

split_symbols <- function(x) {
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0L) {
    return(character(0))
  }
  y <- unlist(strsplit(x, split = ";", fixed = TRUE))
  y <- trimws(y)
  unique(y[!is.na(y) & y != ""])
}

map_symbol_to_entrez <- function(symbols) {
  symbols <- unique(symbols[!is.na(symbols) & symbols != ""])
  if (length(symbols) == 0L) {
    return(data.frame(SYMBOL = character(0), ENTREZID = character(0)))
  }
  map <- suppressMessages(AnnotationDbi::select(
    org.Mm.eg.db,
    keys = symbols,
    keytype = "SYMBOL",
    columns = c("SYMBOL", "ENTREZID")
  ))
  unique(map[!is.na(map$ENTREZID) & map$ENTREZID != "", c("SYMBOL", "ENTREZID")])
}

wrap_label <- function(x, width = 48) {
  vapply(x, function(z) paste(strwrap(z, width = width), collapse = "\n"), character(1))
}

background_rows <- is.finite(result$AdjPValue)
foreground_rows <- is.finite(result$AdjPValue) &
  result$AdjPValue < PVALUE_CUTOFF &
  abs(result$Log2FC) > LOG2FC_CUTOFF

source(file.path(script_dir, "go_mapping.R"))
mapping_policy <- as.character(option_value(options_cli, "go_mapping", "unique"))
mapped <- prepare_go_gene_sets(result, background_rows, foreground_rows, mapping_policy)
background_map <- mapped$background
foreground_map <- mapped$foreground
write.csv(mapped$audit, file.path(output_dir, "BLA_G4_vs_G1_GO_mapping_audit.csv"), row.names=FALSE)
background_entrez <- unique(background_map$ENTREZID)
foreground_entrez <- intersect(unique(foreground_map$ENTREZID), background_entrez)

if (length(background_entrez) == 0L) {
  stop("No background genes could be mapped to Entrez IDs.")
}
# An empty foreground is valid: all eligible terms have P=1 and the figure reports no candidates.

write.csv(
  foreground_map[foreground_map$ENTREZID %in% foreground_entrez, ],
  file = file.path(output_dir, "BLA_G4_vs_G1_GO_foreground_genes_BH_FC1p5.csv"),
  row.names = FALSE
)

go_map <- suppressMessages(AnnotationDbi::select(
  org.Mm.eg.db,
  keys = background_entrez,
  keytype = "ENTREZID",
  columns = c("ENTREZID", "GOALL", "ONTOLOGYALL")
))
go_map <- go_map[!is.na(go_map$GOALL) & go_map$GOALL != "", ]
go_map <- unique(go_map[, c("ENTREZID", "GOALL", "ONTOLOGYALL")])
colnames(go_map) <- c("ENTREZID", "GOID", "Ontology")

go_size <- aggregate(ENTREZID ~ GOID + Ontology, data = go_map, FUN = function(x) length(unique(x)))
colnames(go_size)[3] <- "TermSize"
go_size <- go_size[go_size$TermSize >= min_go_size & go_size$TermSize <= max_go_size, ]
go_map <- merge(go_map, go_size[, c("GOID", "Ontology")], by = c("GOID", "Ontology"))
if (nrow(go_map) == 0L) {
  stop("No GO annotations remained after term-size filtering.")
}

universe_n <- length(background_entrez)
foreground_n <- length(foreground_entrez)
go_list <- split(go_map$ENTREZID, go_map$GOID)

enrich_list <- lapply(names(go_list), function(go_id) {
  term_genes <- unique(go_list[[go_id]])
  overlap_genes <- intersect(foreground_entrez, term_genes)
  overlap_n <- length(overlap_genes)
  term_n <- length(term_genes)
  p_value <- phyper(q = overlap_n - 1, m = term_n, n = universe_n - term_n, k = foreground_n, lower.tail = FALSE)
  data.frame(
    GOID = go_id,
    Count = overlap_n,
    ForegroundSize = foreground_n,
    TermSize = term_n,
    BackgroundSize = universe_n,
    GeneRatio = if (foreground_n > 0L) overlap_n / foreground_n else 0,
    RichFactor = overlap_n / term_n,
    PValue = p_value,
    EntrezIDs = paste(overlap_genes, collapse = ";"),
    stringsAsFactors = FALSE
  )
})
enrich <- do.call(rbind, enrich_list)
if (is.null(enrich) || nrow(enrich) == 0L) {
  stop("No GO terms overlapped with the foreground gene set.")
}

if (has_go_db) {
  term_info <- suppressMessages(AnnotationDbi::select(
    GO.db::GO.db,
    keys = unique(enrich$GOID),
    keytype = "GOID",
    columns = c("GOID", "TERM", "ONTOLOGY")
  ))
  term_info <- unique(term_info)
} else {
  go_lookup <- unique(go_map[, c("GOID", "Ontology")])
  term_info <- data.frame(
    GOID = unique(enrich$GOID),
    TERM = unique(enrich$GOID),
    ONTOLOGY = go_lookup$Ontology[match(unique(enrich$GOID), go_lookup$GOID)],
    stringsAsFactors = FALSE
  )
}

enrich <- merge(enrich, term_info, by = "GOID", all.x = TRUE)
enrich$Ontology <- enrich$ONTOLOGY
enrich$TERM[is.na(enrich$TERM)] <- enrich$GOID[is.na(enrich$TERM)]

entrez_to_symbol <- unique(background_map[, c("ENTREZID", "SYMBOL")])
enrich$GeneSymbols <- vapply(
  strsplit(enrich$EntrezIDs, split = ";", fixed = TRUE),
  function(ids) paste(unique(entrez_to_symbol$SYMBOL[entrez_to_symbol$ENTREZID %in% ids]), collapse = ";"),
  character(1)
)

enrich$AdjPValue_BH <- p.adjust(enrich$PValue, method = "BH")
enrich$MappingPolicy <- mapping_policy
enrich$AdjustmentFamily <- "all eligible GO terms; BP/CC/MF pooled"
enrich <- enrich[, c(
  "GOID", "TERM", "Ontology", "Count", "ForegroundSize", "TermSize",
  "BackgroundSize", "GeneRatio", "RichFactor", "PValue", "AdjPValue_BH", "MappingPolicy", "AdjustmentFamily", "GeneSymbols", "EntrezIDs"
)]
enrich <- enrich[order(enrich$PValue, -enrich$Count), ]
write.csv(enrich, file = file.path(output_dir, "BLA_G4_vs_G1_GO_enrichment_results_BH_FC1p5.csv"), row.names = FALSE, na = "")

plot_data <- enrich[enrich[[go_field]] < go_cutoff & enrich$Count > 0L, ]
plot_data <- plot_data[order(plot_data[[go_field]], -plot_data$Count), ]
write.csv(plot_data, file.path(output_dir, paste0("BLA_G4_vs_G1_GO_selected_", go_pvalue, ".csv")), row.names = FALSE)
write.csv(data.frame(TermsTested=nrow(enrich), ForegroundGenes=foreground_n,
                    BackgroundGenes=universe_n, PValueField=go_field, Cutoff=go_cutoff,
                    EligibleCandidates=nrow(plot_data), MappingPolicy=mapping_policy,
                    ForegroundProteinGroups=mapped$foreground_groups, BackgroundProteinGroups=mapped$background_groups, BHFamily="BP/CC/MF pooled including zero-overlap terms"),
          file.path(output_dir, "BLA_G4_vs_G1_GO_analysis_summary.csv"), row.names = FALSE)
if (nrow(plot_data) == 0L) {
  open_figure_pdf(output_pdf, width=4.15, height=3.35, family=figure_font)
  plot.new()
  text(0.5, 0.55, "No GO terms passed the selected threshold", cex=0.8)
  text(0.5, 0.4, paste(go_field, "<", go_cutoff), cex=0.7)
  dev.off()
} else {
plot_data <- plot_data[1:min(top_n, nrow(plot_data)), ]
plot_data <- plot_data[order(plot_data$RichFactor), ]
plot_data$MinusLog10PValue <- -log10(pmax(plot_data[[go_field]], .Machine$double.xmin))
plot_data$Label <- wrap_label(paste0(plot_data$TERM, " [", plot_data$Ontology, "]"), width = 36)

color_value <- plot_data$MinusLog10PValue
color_palette <- colorRampPalette(c("#2B6CB0", "#88BBD6", "#F7F7F7", "#F4A261", "#B2182B"))(100)
color_index <- if (length(unique(color_value)) == 1L) rep(60, length(color_value)) else as.integer(cut(color_value, breaks = 100, include.lowest = TRUE))
point_color <- color_palette[color_index]
count_cex <- function(count) 0.75 + 1.35 * sqrt(count / max(plot_data$Count))
point_cex <- count_cex(plot_data$Count)

open_figure_pdf(output_pdf, width = 8, height = 5.8, family = figure_font)
par(mar = c(3.2, 15, 1.6, 3.3), mgp = c(1.95, 0.52, 0), las = 1, family = figure_font, cex.axis = 0.75, cex.lab = 0.85)
y_pos <- seq_len(nrow(plot_data))
plot(
  plot_data$RichFactor,
  y_pos,
  pch = 21,
  bg = point_color,
  col = "#27313A",
  lwd = 0.75,
  cex = point_cex,
  xlim = c(0, max(plot_data$RichFactor) * 1.24),
  ylim = c(0.5, nrow(plot_data) + 0.5),
  yaxt = "n",
  xlab = "Rich factor",
  ylab = "",
  main = "",
  sub = "",
  bty = "l"
)
axis(2, at = y_pos, labels = FALSE, tick = FALSE)
# axis() suppresses colliding multi-line labels; place every selected label explicitly.
label_usr <- par("usr")
text(label_usr[1] - diff(label_usr[1:2]) * 0.035, y_pos,
     labels = plot_data$Label, adj = c(1, 0.5), cex = 0.70, xpd = NA)
mtext("BLA G4-G1: unique-gene GO enrichment", side=3, line=0.35, cex=0.9)
grid(nx = NA, ny = NULL, col = "#EEF1F4", lty = 1, lwd = 0.50)

legend_counts <- sort(unique(plot_data$Count))
if (length(legend_counts) > 3L) {
  legend_counts <- unique(round(seq(min(legend_counts), max(legend_counts), length.out = 3)))
}
legend_cex <- count_cex(legend_counts)
legend(
  "bottomright",
  legend = paste0("Count = ", legend_counts),
  pt.cex = legend_cex,
  pt.bg = "#D8DCE1",
  col = "#343A40",
  pch = 21,
  bty = "n",
  title = "Gene count",
  cex = 0.52,
  y.intersp = 2.0,
  inset = c(0.01, 0.03)
)

usr <- par("usr")
x0 <- usr[2] + diff(usr[1:2]) * 0.03
x1 <- usr[2] + diff(usr[1:2]) * 0.055
y0 <- usr[3] + diff(usr[3:4]) * 0.20
y1 <- usr[3] + diff(usr[3:4]) * 0.54
legend_y <- seq(y0, y1, length.out = length(color_palette) + 1)
for (i in seq_along(color_palette)) {
  rect(x0, legend_y[i], x1, legend_y[i + 1], col = color_palette[i], border = NA, xpd = TRUE)
}
axis_ticks <- pretty(range(color_value), n = 3)
axis_ticks <- axis_ticks[axis_ticks >= min(color_value) & axis_ticks <= max(color_value)]
if (length(axis_ticks) > 0L) {
  axis_y <- if (length(unique(color_value)) == 1L) rep((y0 + y1) / 2, length(axis_ticks)) else y0 + (axis_ticks - min(color_value)) / diff(range(color_value)) * (y1 - y0)
  segments(x1, axis_y, x1 + diff(usr[1:2]) * 0.01, axis_y, xpd = TRUE, col = "#343A40", lwd = 0.55)
  text(x1 + diff(usr[1:2]) * 0.015, axis_y, labels = signif(axis_ticks, 3), adj = c(0, 0.5), cex = 0.50, xpd = TRUE)
}
text((x0 + x1) / 2, y1 + diff(usr[3:4]) * 0.035, labels = if (go_pvalue == "adjusted") "-log10(BH P)" else "-log10(raw P)", cex = 0.52, xpd = TRUE)
dev.off()

}
