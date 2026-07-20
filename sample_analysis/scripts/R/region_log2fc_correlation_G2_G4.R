#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
rm(list = ls())

suppressPackageStartupMessages(library(limma))

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_file <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[1]) else NA_character_
script_dir <- if (!is.na(script_file) && file.exists(script_file)) dirname(normalizePath(script_file, winslash = "/", mustWork = TRUE)) else getwd()
source(file.path(script_dir, "common.R"))

options_cli <- parse_cli_options()
paths <- analysis_paths("brain_region_log2FC_correlation.pdf", options_cli)
data_dir <- paths$data_dir
output_dir <- paths$output_dir
output_pdf <- paths$output_pdf
figure_font <- paths$font_family
cor_method <- "pearson"
min_pairwise_proteins <- as.integer(option_value(options_cli, "min_pairwise_proteins", 30L))

region_info <- data.frame(
  FileCode = c("PMD", "LPB", "VHPC", "PVH", "BLA", "DHPC", "CEA", "MSC", "SSC", "IL", "PL"),
  DisplayName = c("PMD", "LPB", "vHPC", "PVH", "BLA", "dHPC", "CEA", "MSC", "SSC", "IL", "PL"),
  stringsAsFactors = FALSE
)

all_result <- do.call(rbind, lapply(seq_len(nrow(region_info)), function(i) {
  analysis <- run_region_differential(region_info$FileCode[i], region_info$DisplayName[i], data_dir)
  result <- analysis$result
  result$Significant <- is.finite(result$PValue) &
    result$PValue < PVALUE_CUTOFF &
    abs(result$Log2FC_G4_vs_G2) > LOG2FC_CUTOFF
  result[, c(
    "Region", "RegionFileCode", "ProteinAccessions", "Genes",
    "Log2FC_G4_vs_G2", "PValue", "Status", "Significant"
  )]
}))

all_result$Label <- make_label(all_result$Genes, all_result$ProteinAccessions)
all_result <- all_result[order(all_result$Region, all_result$PValue, -abs(all_result$Log2FC_G4_vs_G2), na.last = TRUE), ]
write.csv(all_result, file = file.path(output_dir, "all_regions_G2_vs_G4_limma_results_pvalue_FC1p5.csv"), row.names = FALSE, na = "")

sig_result <- all_result[all_result$Significant, ]
if (nrow(sig_result) == 0L) {
  stop("No differential protein groups passed raw P value < 0.05 and FC > 1.5.")
}

sig_proteins <- sort(unique(sig_result$ProteinAccessions))
log2fc_matrix <- matrix(
  NA_real_,
  nrow = length(sig_proteins),
  ncol = nrow(region_info),
  dimnames = list(sig_proteins, region_info$DisplayName)
)
for (region in region_info$DisplayName) {
  sub <- all_result[all_result$Region == region, ]
  sub <- sub[match(sig_proteins, sub$ProteinAccessions), ]
  log2fc_matrix[, region] <- sub$Log2FC_G4_vs_G2
}
log2fc_matrix <- log2fc_matrix[rowSums(is.finite(log2fc_matrix)) >= 2L, , drop = FALSE]

label_map <- all_result[!duplicated(all_result$ProteinAccessions), c("ProteinAccessions", "Genes", "Label")]
label_map <- label_map[match(rownames(log2fc_matrix), label_map$ProteinAccessions), ]
write.csv(
  data.frame(ProteinAccessions = rownames(log2fc_matrix), Genes = label_map$Genes, Label = label_map$Label, log2fc_matrix, check.names = FALSE),
  file = file.path(output_dir, "union_differential_protein_log2fc_matrix.csv"),
  row.names = FALSE,
  na = ""
)

count_table <- as.data.frame.matrix(table(sig_result$Region, sig_result$Status))
count_table$Region <- rownames(count_table)
count_table <- count_table[, c("Region", setdiff(colnames(count_table), "Region")), drop = FALSE]
write.csv(count_table, file = file.path(output_dir, "significant_protein_counts_by_region.csv"), row.names = FALSE, na = "")

cor_matrix <- cor(log2fc_matrix, use = "pairwise.complete.obs", method = cor_method)
if (any(!is.finite(cor_matrix))) {
  stop("Correlation matrix contains non-finite values.")
}

pair_n <- outer(
  seq_len(ncol(log2fc_matrix)),
  seq_len(ncol(log2fc_matrix)),
  Vectorize(function(i, j) sum(is.finite(log2fc_matrix[, i]) & is.finite(log2fc_matrix[, j])))
)
dimnames(pair_n) <- list(colnames(log2fc_matrix), colnames(log2fc_matrix))
off_diagonal_pair_n <- pair_n[upper.tri(pair_n)]
if (any(off_diagonal_pair_n < min_pairwise_proteins)) {
  low_pairs <- which(pair_n < min_pairwise_proteins & upper.tri(pair_n), arr.ind = TRUE)
  low_pair_labels <- apply(low_pairs, 1, function(index) {
    paste0(rownames(pair_n)[index[1]], "-", colnames(pair_n)[index[2]], " (n=", pair_n[index[1], index[2]], ")")
  })
  stop("Insufficient pairwise protein coverage for correlation: ", paste(low_pair_labels, collapse = ", "), ".")
}

write.csv(data.frame(Region = rownames(cor_matrix), cor_matrix, check.names = FALSE), file = file.path(output_dir, "brain_region_log2fc_correlation_matrix.csv"), row.names = FALSE, na = "")
write.csv(data.frame(Region = rownames(pair_n), pair_n, check.names = FALSE), file = file.path(output_dir, "brain_region_log2fc_correlation_pairwise_n.csv"), row.names = FALSE, na = "")

plot_order <- hclust(as.dist(1 - cor_matrix), method = "complete")$order
cor_plot <- cor_matrix[plot_order, plot_order, drop = FALSE]

heat_colors <- colorRampPalette(COLOR_DIVERGING)(101)
breaks <- seq(-1, 1, length.out = length(heat_colors) + 1)

open_figure_pdf(output_pdf, width = 7.65, height = 2.55, family = figure_font)
par(mar = c(2.05, 2.45, 0.75, 3.80), mgp = c(1.25, 0.36, 0), family = figure_font, xpd = NA)
image(
  x = seq_len(ncol(cor_plot)),
  y = seq_len(nrow(cor_plot)),
  z = t(cor_plot[nrow(cor_plot):1, , drop = FALSE]),
  col = heat_colors,
  breaks = breaks,
  axes = FALSE,
  xlab = "",
  ylab = "",
  main = ""
)
axis(1, at = seq_len(ncol(cor_plot)), labels = colnames(cor_plot), las = 2, cex.axis = 0.50, tick = FALSE, line = 0.05)
axis(2, at = seq_len(nrow(cor_plot)), labels = rev(rownames(cor_plot)), las = 1, cex.axis = 0.50, tick = FALSE, line = 0.05)
abline(v = seq(0.5, ncol(cor_plot) + 0.5), col = "#FFFFFF", lwd = 0.8)
abline(h = seq(0.5, nrow(cor_plot) + 0.5), col = "#FFFFFF", lwd = 0.8)
box(lwd = 0.85, col = "#2F363D")

for (i in seq_len(nrow(cor_plot))) {
  for (j in seq_len(ncol(cor_plot))) {
    value <- cor_plot[i, j]
    text_color <- if (abs(value) >= 0.62) "#FFFFFF" else "#1F252B"
    text(j, nrow(cor_plot) - i + 1, labels = sprintf("%.2f", value), cex = 0.36, col = text_color)
  }
}

usr <- par("usr")
legend_x0 <- usr[2] + diff(usr[1:2]) * 0.035
legend_x1 <- usr[2] + diff(usr[1:2]) * 0.070
legend_y0 <- usr[3] + diff(usr[3:4]) * 0.12
legend_y1 <- usr[3] + diff(usr[3:4]) * 0.88
legend_y <- seq(legend_y0, legend_y1, length.out = length(heat_colors) + 1)
for (i in seq_along(heat_colors)) {
  rect(legend_x0, legend_y[i], legend_x1, legend_y[i + 1], col = heat_colors[i], border = NA)
}
rect(legend_x0, legend_y0, legend_x1, legend_y1, border = "#2F363D", lwd = 0.6)
legend_ticks <- seq(-1, 1, by = 0.5)
legend_tick_y <- legend_y0 + (legend_ticks + 1) / 2 * (legend_y1 - legend_y0)
segments(legend_x1, legend_tick_y, legend_x1 + diff(usr[1:2]) * 0.012, legend_tick_y, col = "#2F363D", lwd = 0.55)
text(legend_x1 + diff(usr[1:2]) * 0.015, legend_tick_y, labels = sprintf("%.1f", legend_ticks), adj = c(0, 0.5), cex = 0.42)
text((legend_x0 + legend_x1) / 2, legend_y1 + diff(usr[3:4]) * 0.060, labels = "Pearson r", cex = 0.46)
dev.off()
