#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
rm(list = ls())

suppressPackageStartupMessages(library(limma))

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_file <- if (length(script_arg) > 0) sub("^--file=", "", script_arg[1]) else NA_character_
script_dir <- if (!is.na(script_file) && file.exists(script_file)) dirname(normalizePath(script_file, winslash = "/", mustWork = TRUE)) else getwd()
source(file.path(script_dir, "common.R"))

paths <- analysis_paths("multi_region_G2_vs_G4_heatmap.pdf")
data_dir <- paths$data_dir
output_dir <- paths$output_dir
output_pdf <- paths$output_pdf
figure_font <- paths$font_family

# The selected source script displays up to 3 high-G2 and 3 high-G4 proteins
# for each of the six requested regions, for a maximum of 36 rows.
proteins_per_direction_per_region <- 3L
max_heatmap_proteins <- 36L
min_observed_heatmap_columns <- 6L

region_info <- data.frame(
  FileCode = c("BLA", "CEA", "DHPC", "VHPC", "PL", "IL"),
  DisplayName = c("BLA", "CEA", "dHPC", "vHPC", "PL", "IL"),
  stringsAsFactors = FALSE
)

region_results <- lapply(seq_len(nrow(region_info)), function(i) {
  run_region_differential(region_info$FileCode[i], region_info$DisplayName[i], data_dir)
})
names(region_results) <- region_info$DisplayName

diff_result <- do.call(rbind, lapply(region_results, function(x) x$result))
mean_result <- do.call(rbind, lapply(region_results, function(x) x$mean_table))
write.csv(diff_result, file = file.path(output_dir, "multi_region_G2_vs_G4_differential_results.csv"), row.names = FALSE, na = "")

sig_result <- diff_result[
  is.finite(diff_result$PValue) &
    diff_result$PValue < PVALUE_CUTOFF &
    abs(diff_result$Log2FC_G4_vs_G2) > LOG2FC_CUTOFF,
]
if (nrow(sig_result) == 0L) {
  stop("No significant protein groups passed raw P value < 0.05 and FC > 1.5 in the selected regions.")
}

best_p <- tapply(sig_result$PValue, sig_result$ProteinAccessions, min, na.rm = TRUE)
max_abs_fc <- tapply(abs(sig_result$Log2FC_G4_vs_G2), sig_result$ProteinAccessions, max, na.rm = TRUE)
primary_region <- vapply(names(best_p), function(protein) {
  sub <- sig_result[sig_result$ProteinAccessions == protein, ]
  sub$RegionOrder <- match(sub$Region, region_info$DisplayName)
  sub <- sub[order(sub$PValue, -abs(sub$Log2FC_G4_vs_G2), sub$RegionOrder), ]
  sub$Region[1]
}, character(1))

rank_table <- data.frame(
  ProteinAccessions = names(best_p),
  BestPValue = as.numeric(best_p),
  MaxAbsLog2FC = as.numeric(max_abs_fc),
  SignificantRegions = as.integer(table(sig_result$ProteinAccessions)[names(best_p)]),
  PrimaryRegion = primary_region,
  stringsAsFactors = FALSE
)
rank_table$PrimaryRegionOrder <- match(rank_table$PrimaryRegion, region_info$DisplayName)

column_info <- data.frame(
  Region = rep(region_info$DisplayName, each = 2),
  Treatment = rep(c("G2", "G4"), times = nrow(region_info)),
  stringsAsFactors = FALSE
)
column_info$ColumnName <- paste(column_info$Region, column_info$Treatment, sep = "_")

candidate_proteins <- rank_table$ProteinAccessions
candidate_mean_matrix <- matrix(
  NA_real_,
  nrow = length(candidate_proteins),
  ncol = nrow(column_info),
  dimnames = list(candidate_proteins, column_info$ColumnName)
)
for (i in seq_len(nrow(column_info))) {
  region <- column_info$Region[i]
  treatment <- column_info$Treatment[i]
  sub <- mean_result[mean_result$Region == region, ]
  sub <- sub[match(candidate_proteins, sub$ProteinAccessions), ]
  candidate_mean_matrix[, i] <- sub[[treatment]]
}

candidate_imputed_matrix <- candidate_mean_matrix
candidate_observed_columns <- rowSums(is.finite(candidate_mean_matrix))
for (i in seq_len(nrow(candidate_imputed_matrix))) {
  finite_value <- candidate_imputed_matrix[i, is.finite(candidate_imputed_matrix[i, ])]
  missing_index <- !is.finite(candidate_imputed_matrix[i, ])
  if (length(finite_value) > 0L && any(missing_index)) {
    candidate_imputed_matrix[i, missing_index] <- mean(finite_value)
  }
}

candidate_z_matrix <- row_zscore(candidate_imputed_matrix)
candidate_valid_row <- rowSums(is.finite(candidate_z_matrix)) >= 2L
candidate_z_matrix <- candidate_z_matrix[candidate_valid_row, , drop = FALSE]

zscore_region_rank <- do.call(rbind, lapply(region_info$DisplayName, function(region) {
  region_sig <- sig_result[
    sig_result$Region == region &
      sig_result$ProteinAccessions %in% rownames(candidate_z_matrix),
  ]
  if (nrow(region_sig) == 0L) {
    return(NULL)
  }
  g2_name <- paste(region, "G2", sep = "_")
  g4_name <- paste(region, "G4", sep = "_")
  region_sig$SelectedRegionZDiff <- abs(
    candidate_z_matrix[region_sig$ProteinAccessions, g4_name] -
      candidate_z_matrix[region_sig$ProteinAccessions, g2_name]
  )
  other_regions <- setdiff(region_info$DisplayName, region)
  other_pair_zdiff <- sapply(other_regions, function(other_region) {
    abs(
      candidate_z_matrix[region_sig$ProteinAccessions, paste(other_region, "G4", sep = "_")] -
        candidate_z_matrix[region_sig$ProteinAccessions, paste(other_region, "G2", sep = "_")]
    )
  })
  if (is.null(dim(other_pair_zdiff))) {
    other_pair_zdiff <- matrix(other_pair_zdiff, nrow = nrow(region_sig))
  }
  region_sig$OtherMaxZDiff <- apply(other_pair_zdiff, 1, max, na.rm = TRUE)
  region_sig$SelectedRegionSpecificity <- region_sig$SelectedRegionZDiff - region_sig$OtherMaxZDiff
  region_sig$ObservedHeatmapColumns <- candidate_observed_columns[region_sig$ProteinAccessions]
  region_sig$CoveragePass <- region_sig$ObservedHeatmapColumns >= min_observed_heatmap_columns
  region_sig$RegionDominant <- region_sig$SelectedRegionSpecificity >= 0
  region_sig$VisualScore <-
    region_sig$SelectedRegionZDiff *
      sqrt(pmin(region_sig$ObservedHeatmapColumns, ncol(candidate_mean_matrix)) / ncol(candidate_mean_matrix)) +
      pmax(region_sig$SelectedRegionSpecificity, 0) * 0.35
  region_sig <- region_sig[
    order(
      -as.integer(region_sig$CoveragePass),
      -as.integer(region_sig$RegionDominant),
      -region_sig$VisualScore,
      -region_sig$SelectedRegionZDiff,
      region_sig$PValue,
      -abs(region_sig$Log2FC_G4_vs_G2)
    ),
  ]
  region_sig[!duplicated(region_sig$ProteinAccessions), ]
}))
if (is.null(zscore_region_rank) || nrow(zscore_region_rank) == 0L) {
  stop("No significant protein groups could be ranked after candidate-level Z-score transformation.")
}

direction_info <- data.frame(
  Direction = c("Down_in_G4", "Up_in_G4"),
  DirectionOrder = c(1L, 2L),
  stringsAsFactors = FALSE
)
selected_records <- list()
selected_seen <- character(0)
for (region in region_info$DisplayName) {
  region_sig <- zscore_region_rank[zscore_region_rank$Region == region, ]
  if (nrow(region_sig) == 0L) {
    next
  }
  for (i in seq_len(nrow(direction_info))) {
    direction <- direction_info$Direction[i]
    direction_order <- direction_info$DirectionOrder[i]
    direction_sig <- if (direction == "Up_in_G4") {
      region_sig[region_sig$Log2FC_G4_vs_G2 > 0, ]
    } else {
      region_sig[region_sig$Log2FC_G4_vs_G2 < 0, ]
    }
    direction_sig <- direction_sig[!(direction_sig$ProteinAccessions %in% selected_seen), ]
    if (nrow(direction_sig) == 0L) {
      next
    }
    selected_sub <- direction_sig[seq_len(min(proteins_per_direction_per_region, nrow(direction_sig))), ]
    selected_sub$SelectedRegion <- region
    selected_sub$SelectedRegionOrder <- match(region, region_info$DisplayName)
    selected_sub$SelectedDirection <- direction
    selected_sub$SelectedDirectionOrder <- direction_order
    selected_sub$SelectedDirectionRank <- seq_len(nrow(selected_sub))
    selected_sub$SelectedRegionRank <- (selected_sub$SelectedDirectionOrder - 1L) * proteins_per_direction_per_region + selected_sub$SelectedDirectionRank
    selected_sub$SelectedRegionPValue <- selected_sub$PValue
    selected_sub$SelectedRegionLog2FC <- selected_sub$Log2FC_G4_vs_G2
    selected_records[[paste(region, direction, sep = "_")]] <- selected_sub[, c(
      "ProteinAccessions", "SelectedRegion", "SelectedRegionOrder",
      "SelectedDirection", "SelectedDirectionOrder", "SelectedDirectionRank",
      "SelectedRegionRank", "SelectedRegionPValue", "SelectedRegionLog2FC",
      "SelectedRegionZDiff", "SelectedRegionSpecificity", "OtherMaxZDiff",
      "ObservedHeatmapColumns", "CoveragePass", "RegionDominant", "VisualScore"
    )]
    selected_seen <- c(selected_seen, selected_sub$ProteinAccessions)
  }
}
if (length(selected_records) == 0L) {
  stop("No region-specific significant protein groups were selected for heatmap display.")
}

selected_order_table <- do.call(rbind, selected_records)
selected_order_table <- selected_order_table[
  order(
    selected_order_table$SelectedRegionOrder,
    selected_order_table$SelectedDirectionOrder,
    selected_order_table$SelectedDirectionRank,
    selected_order_table$SelectedRegionPValue,
    -abs(selected_order_table$SelectedRegionLog2FC)
  ),
]
selected_proteins <- selected_order_table$ProteinAccessions
if (length(selected_proteins) > max_heatmap_proteins) {
  selected_proteins <- selected_proteins[seq_len(max_heatmap_proteins)]
}

annotation_map <- diff_result[!duplicated(diff_result$ProteinAccessions), c("ProteinAccessions", "Genes", "ProteinDescriptions")]
selected_info <- merge(rank_table[rank_table$ProteinAccessions %in% selected_proteins, ], selected_order_table, by = "ProteinAccessions", all.x = TRUE, sort = FALSE)
selected_info <- merge(selected_info, annotation_map, by = "ProteinAccessions", all.x = TRUE, sort = FALSE)
selected_info$Label <- make_label(selected_info$Genes, selected_info$ProteinAccessions)
selected_info <- selected_info[match(selected_proteins, selected_info$ProteinAccessions), ]

write.csv(rank_table, file = file.path(output_dir, "multi_region_G2_vs_G4_significant_union.csv"), row.names = FALSE, na = "")
write.csv(selected_info, file = file.path(output_dir, "multi_region_G2_vs_G4_heatmap_selected_proteins.csv"), row.names = FALSE, na = "")

mean_matrix <- matrix(
  NA_real_,
  nrow = length(selected_proteins),
  ncol = nrow(column_info),
  dimnames = list(selected_proteins, column_info$ColumnName)
)
for (i in seq_len(nrow(column_info))) {
  region <- column_info$Region[i]
  treatment <- column_info$Treatment[i]
  sub <- mean_result[mean_result$Region == region, ]
  sub <- sub[match(selected_proteins, sub$ProteinAccessions), ]
  mean_matrix[, i] <- sub[[treatment]]
}
write.csv(data.frame(ProteinAccessions = rownames(mean_matrix), mean_matrix, check.names = FALSE), file = file.path(output_dir, "multi_region_G2_vs_G4_heatmap_log2_group_means.csv"), row.names = FALSE, na = "")

imputed_mean_matrix <- mean_matrix
imputed_count <- 0L
imputation_value <- rep(NA_real_, nrow(imputed_mean_matrix))
names(imputation_value) <- rownames(imputed_mean_matrix)
for (i in seq_len(nrow(imputed_mean_matrix))) {
  finite_value <- imputed_mean_matrix[i, is.finite(imputed_mean_matrix[i, ])]
  missing_index <- !is.finite(imputed_mean_matrix[i, ])
  if (length(finite_value) > 0L && any(missing_index)) {
    fill_value <- mean(finite_value)
    imputed_mean_matrix[i, missing_index] <- fill_value
    imputation_value[i] <- fill_value
    imputed_count <- imputed_count + sum(missing_index)
  }
}
write.csv(data.frame(ProteinAccessions = rownames(imputed_mean_matrix), imputed_mean_matrix, check.names = FALSE), file = file.path(output_dir, "multi_region_G2_vs_G4_heatmap_log2_group_means_imputed.csv"), row.names = FALSE, na = "")

z_matrix <- row_zscore(imputed_mean_matrix)
valid_row <- rowSums(is.finite(z_matrix)) >= 2L
z_matrix <- z_matrix[valid_row, , drop = FALSE]
mean_matrix <- mean_matrix[valid_row, , drop = FALSE]
imputed_mean_matrix <- imputed_mean_matrix[valid_row, , drop = FALSE]
selected_info <- selected_info[match(rownames(z_matrix), selected_info$ProteinAccessions), ]
imputation_value <- imputation_value[rownames(z_matrix)]
if (nrow(z_matrix) < 2L) {
  stop("Fewer than two protein groups remained for heatmap plotting after Z-score filtering.")
}
write.csv(data.frame(ProteinAccessions = rownames(z_matrix), z_matrix, check.names = FALSE), file = file.path(output_dir, "multi_region_G2_vs_G4_heatmap_zscore_matrix.csv"), row.names = FALSE, na = "")

write.csv(
  data.frame(
    ProteinAccessions = rownames(z_matrix),
    ImputationValue = imputation_value,
    stringsAsFactors = FALSE
  ),
  file = file.path(output_dir, "multi_region_G2_vs_G4_heatmap_imputation_summary.csv"),
  row.names = FALSE,
  na = ""
)

# Keep the selected-protein order: region-balanced, then direction-balanced.
z_plot <- z_matrix[seq_len(nrow(z_matrix)), , drop = FALSE]
plot_col_labels <- colnames(z_plot)
row_label_plot <- selected_info$Label[seq_len(nrow(z_plot))]
empty_row_labels <- is.na(row_label_plot) | trimws(row_label_plot) == ""
row_label_plot[empty_row_labels] <- rownames(z_plot)[empty_row_labels]
plot_row_labels <- rev(row_label_plot)

heat_colors <- colorRampPalette(c("#2C7BB6", "#F7F7F7", "#D7191C"))(101)
breaks <- seq(-2.5, 2.5, length.out = length(heat_colors) + 1)
z_plot_capped <- pmax(pmin(z_plot, max(breaks)), min(breaks))

open_figure_pdf(output_pdf, width = 4.90, height = 4.25, family = figure_font)
layout(
  matrix(c(0, 2, 1, 2), nrow = 2, byrow = TRUE),
  widths = c(0.45, 4.20),
  heights = c(3.35, 0.65)
)
par(mar = c(0.45, 0.05, 0.05, 0.05), family = figure_font, xpd = NA)
plot.new()
plot.window(xlim = c(0, 1), ylim = c(min(breaks) - 0.15, max(breaks) + 0.15))
legend_y <- seq(min(breaks), max(breaks), length.out = length(heat_colors) + 1)
for (i in seq_along(heat_colors)) {
  rect(0.46, legend_y[i], 0.80, legend_y[i + 1], col = heat_colors[i], border = NA)
}
rect(0.46, min(breaks), 0.80, max(breaks), border = "#1F252B", lwd = 0.6)
text(0.40, max(breaks), labels = sprintf("%.2f", max(breaks)), adj = c(1, 0.5), cex = 0.44)
text(0.40, min(breaks), labels = sprintf("%.2f", min(breaks)), adj = c(1, 0.5), cex = 0.44)

par(mar = c(0.65, 0.55, 2.85, 3.85), family = figure_font, xpd = NA)
image(
  x = seq_len(ncol(z_plot_capped)),
  y = seq_len(nrow(z_plot_capped)),
  z = t(z_plot_capped[nrow(z_plot_capped):1, , drop = FALSE]),
  col = heat_colors,
  breaks = breaks,
  axes = FALSE,
  xlab = "",
  ylab = "",
  main = ""
)
axis(side = 3, at = seq_len(ncol(z_plot_capped)), labels = plot_col_labels, las = 2, cex.axis = 0.46, col.axis = "#000000", line = 0.25, tick = FALSE)
axis(side = 4, at = seq_len(nrow(z_plot_capped)), labels = plot_row_labels, las = 1, cex.axis = 0.31, col.axis = "#000000", line = 0.25, tick = FALSE)
box(lwd = 0.65)
title(main = "", line = 1.8, cex.main = 0.50, font.main = 2)
dev.off()
