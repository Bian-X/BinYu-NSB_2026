#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
rm(list = ls())

suppressPackageStartupMessages(library(limma))

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_file <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[1]) else NA_character_
script_dir <- if (!is.na(script_file) && file.exists(script_file)) dirname(normalizePath(script_file, winslash = "/", mustWork = TRUE)) else getwd()
source(file.path(script_dir, "common.R"))

paths <- analysis_paths("BLA_volcano_G4_vs_G2.pdf")
data_dir <- paths$data_dir
output_dir <- paths$output_dir
output_pdf <- paths$output_pdf
figure_font <- paths$font_family
max_label_each_direction <- 7L

analysis <- run_region_differential("BLA", "BLA", data_dir)
result <- analysis$result
result$Status[result$Status == "High in G2"] <- "Low in G4"

write.csv(
  result,
  file = file.path(output_dir, "BLA_G2_vs_G4_differential_results_pvalue_FC1p5.csv"),
  row.names = FALSE,
  na = ""
)

plot_data <- result[is.finite(result$Log2FC_G4_vs_G2) & is.finite(result$PValue), ]
if (nrow(plot_data) == 0L) {
  stop("No finite differential-analysis results are available for plotting.")
}

plot_data$MinusLog10PValue <- -log10(pmax(plot_data$PValue, .Machine$double.xmin))
plot_data$LabelScore <- abs(plot_data$Log2FC_G4_vs_G2) * plot_data$MinusLog10PValue
plot_data$Label <- first_candidate(plot_data$Genes)
empty_label <- is.na(plot_data$Label) | plot_data$Label == ""
plot_data$Label[empty_label] <- first_candidate(plot_data$ProteinAccessions[empty_label])

up_data <- plot_data[plot_data$Status == "High in G4", ]
down_data <- plot_data[plot_data$Status == "Low in G4", ]
if (nrow(up_data) > max_label_each_direction) {
  up_data <- up_data[order(-up_data$LabelScore, up_data$PValue, -abs(up_data$Log2FC_G4_vs_G2)), ][1:max_label_each_direction, ]
}
if (nrow(down_data) > max_label_each_direction) {
  down_data <- down_data[order(-down_data$LabelScore, down_data$PValue, -abs(down_data$Log2FC_G4_vs_G2)), ][1:max_label_each_direction, ]
}
label_data <- rbind(up_data, down_data)

color_map <- c(
  "Not significant" = COLOR_NOT_SIGNIFICANT,
  "Low in G4" = COLOR_HIGH_G2,
  "High in G4" = COLOR_HIGH_G4
)
point_color <- color_map[plot_data$Status]
x_lim <- max(abs(plot_data$Log2FC_G4_vs_G2), LOG2FC_CUTOFF, na.rm = TRUE) * 1.12
y_lim <- max(plot_data$MinusLog10PValue, -log10(PVALUE_CUTOFF), na.rm = TRUE) * 1.15

open_figure_pdf(output_pdf, width = 4.15, height = 3.35, family = figure_font)
layout(matrix(c(1, 2), nrow = 1), widths = c(4.9, 1.35))
par(
  mar = c(3.4, 3.6, 1.15, 0.25),
  mgp = c(2.0, 0.55, 0),
  las = 1,
  family = figure_font,
  cex.axis = 0.62,
  cex.lab = 0.72
)

plot(
  plot_data$Log2FC_G4_vs_G2,
  plot_data$MinusLog10PValue,
  pch = 16,
  cex = 0.42,
  col = adjustcolor(point_color, alpha.f = 0.78),
  xlim = c(-x_lim, x_lim),
  ylim = c(0, y_lim),
  xlab = "log2 fold change (G4/G2)",
  ylab = "-log10 P value",
  main = "",
  sub = "",
  bty = "l"
)
abline(v = c(-LOG2FC_CUTOFF, LOG2FC_CUTOFF), lty = 2, col = "#555B61", lwd = 1)
abline(h = -log10(PVALUE_CUTOFF), lty = 2, col = "#555B61", lwd = 1)

if (nrow(label_data) > 0L) {
  high_label <- label_data[label_data$Status == "High in G4", ]
  low_label <- label_data[label_data$Status == "Low in G4", ]

  if (nrow(high_label) > 0L) {
    high_label <- high_label[order(high_label$MinusLog10PValue, decreasing = TRUE), ]
    high_y_top <- min(y_lim * 0.93, max(high_label$MinusLog10PValue) + 0.12)
    high_y_bottom <- max(-log10(PVALUE_CUTOFF) + 0.06, min(high_label$MinusLog10PValue) - 0.04)
    high_label_y <- if (nrow(high_label) == 1L) high_label$MinusLog10PValue else seq(high_y_top, high_y_bottom, length.out = nrow(high_label))
    high_label_x <- max(LOG2FC_CUTOFF + 0.55, max(high_label$Log2FC_G4_vs_G2) + 0.28)
    high_label_x <- min(high_label_x, x_lim - 0.35)
    segments(
      x0 = high_label$Log2FC_G4_vs_G2,
      y0 = high_label$MinusLog10PValue,
      x1 = rep(high_label_x - 0.05, nrow(high_label)),
      y1 = high_label_y,
      col = adjustcolor(color_map["High in G4"], alpha.f = 0.70),
      lwd = 0.7
    )
    text(rep(high_label_x, nrow(high_label)), high_label_y, labels = high_label$Label, adj = c(0, 0.5), cex = 0.36, col = color_map["High in G4"], xpd = TRUE)
  }

  if (nrow(low_label) > 0L) {
    low_label <- low_label[order(low_label$MinusLog10PValue, decreasing = TRUE), ]
    low_y_top <- min(y_lim * 0.93, max(low_label$MinusLog10PValue) + 0.06)
    low_y_bottom <- max(-log10(PVALUE_CUTOFF) + 0.06, min(low_label$MinusLog10PValue) - 0.04)
    low_label_y <- if (nrow(low_label) == 1L) low_label$MinusLog10PValue else seq(low_y_top, low_y_bottom, length.out = nrow(low_label))
    low_label_x <- min(-LOG2FC_CUTOFF - 0.55, min(low_label$Log2FC_G4_vs_G2) - 0.28)
    low_label_x <- max(low_label_x, -x_lim + 0.35)
    segments(
      x0 = low_label$Log2FC_G4_vs_G2,
      y0 = low_label$MinusLog10PValue,
      x1 = rep(low_label_x + 0.05, nrow(low_label)),
      y1 = low_label_y,
      col = adjustcolor(color_map["Low in G4"], alpha.f = 0.70),
      lwd = 0.7
    )
    text(rep(low_label_x, nrow(low_label)), low_label_y, labels = low_label$Label, adj = c(1, 0.5), cex = 0.36, col = color_map["Low in G4"], xpd = TRUE)
  }
}

par(mar = c(3.4, 0.05, 1.15, 0.05), family = figure_font)
plot.new()
legend(
  "center",
  legend = c(
    paste0("High in G4 (n=", sum(result$Status == "High in G4"), ")"),
    paste0("Low in G4 (n=", sum(result$Status == "Low in G4"), ")"),
    paste0("Not significant (n=", sum(result$Status == "Not significant"), ")")
  ),
  col = c(color_map["High in G4"], color_map["Low in G4"], color_map["Not significant"]),
  pch = 16,
  pt.cex = 0.58,
  bty = "n",
  cex = 0.52
)
dev.off()
