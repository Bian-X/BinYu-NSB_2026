# Regional preprocessing: >=20 detections; exclusion before filtering;
# log2; observed-value sample-median normalization; protein-wise median filling
# on the normalized log2 scale; all-six-group limma with robust/trend eBayes.
# Missingness/QC-CV summaries remain based on observed data without imputation.
options(stringsAsFactors = FALSE)

DEFAULT_FONT_FAMILY <- "Arial"
PVALUE_CUTOFF <- 0.05          # applied to the BH-adjusted P value (FDR)
FC_CUTOFF <- 1.5
LOG2FC_CUTOFF <- log2(FC_CUTOFF)
MIN_DETECTED_SAMPLES <- 20L    # fixed count: 20/24, or 20/23 for primary PVH
ADJUST_METHOD <- "BH"

cli_flag <- function(value) {
  if (is.logical(value)) return(isTRUE(value))
  value <- tolower(as.character(value))
  if (!value %in% c("true", "false", "yes", "no", "1", "0")) stop("Invalid boolean: ", value)
  value %in% c("true", "yes", "1")
}

sample_exclusions <- function() {
  path <- option_value(parse_cli_options(), "exclusions", file.path(project_root(), "sample_exclusions.csv"))
  policy <- read.csv(path, check.names = FALSE)
  if (!all(c("sample_id", "reason") %in% names(policy)) || anyNA(policy$sample_id) || anyDuplicated(policy$sample_id)) {
    stop("Invalid sample exclusion policy: ", path)
  }
  policy
}

clean_sample_name <- function(x) {
  sub("\\.raw\\.PG\\.Quantity$", "", sub("^\\[[0-9]+\\] +", "", x))
}

REGION_CODES <- c("PMD", "LPB", "VHPC", "PVH", "BLA", "DHPC", "CEA", "MSC", "SSC", "IL", "PL")
IRT_STANDARD_ACCESSIONS <- c("iRT-Kit_WR_fusion")

COLOR_HIGH_B <- "#D94B3D"
COLOR_HIGH_A <- "#2878B5"
COLOR_NOT_SIGNIFICANT <- "#B8BDC4"
COLOR_TEXT <- "#20262D"
COLOR_GRID <- "#E7E9EC"
COLOR_DIVERGING <- c("#2166AC", "#F7F7F7", "#B2182B")

script_directory <- function() {
  script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(script_arg) == 0) {
    return(normalizePath(getwd(), winslash = "/", mustWork = TRUE))
  }
  script_file <- sub("^--file=", "", script_arg[1])
  dirname(normalizePath(script_file, winslash = "/", mustWork = TRUE))
}

project_root <- function(script_dir = MODULE_DIR) {
  normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = TRUE)
}

common_source_file <- tryCatch(sys.frame(1)$ofile, error=function(e) NULL)
MODULE_DIR <- if (!is.null(common_source_file)) dirname(normalizePath(common_source_file, winslash="/", mustWork=TRUE)) else script_directory()
source(file.path(MODULE_DIR, "sample_validation.R"))

parse_cli_options <- function(args = commandArgs(trailingOnly = TRUE)) {
  values <- list()
  i <- 1L
  while (i <= length(args)) {
    token <- args[[i]]
    if (!startsWith(token, "--")) {
      stop("Unexpected positional argument: ", token)
    }
    token <- substring(token, 3L)
    if (grepl("=", token, fixed = TRUE)) {
      parts <- strsplit(token, "=", fixed = TRUE)[[1]]
      key <- parts[[1]]
      value <- paste(parts[-1], collapse = "=")
    } else {
      key <- token
      if (i == length(args) || startsWith(args[[i + 1L]], "--")) {
        value <- TRUE
      } else {
        i <- i + 1L
        value <- args[[i]]
      }
    }
    values[[gsub("-", "_", key, fixed = TRUE)]] <- value
    i <- i + 1L
  }
  values
}

option_value <- function(options, key, default) {
  value <- options[[key]]
  if (is.null(value)) default else value
}

has_region_csvs <- function(path) {
  if (!dir.exists(path)) {
    return(FALSE)
  }
  any(grepl("_(BLA|CEA|DHPC|IL|LPB|MSC|PL|PMD|PVH|SSC|VHPC)_DIA_LFQ.*\\.csv$", list.files(path)))
}

detect_data_dir <- function(root = project_root()) {
  candidates <- c(
    file.path(root, "data"),
    file.path(root, "data", "data1"),
    file.path(root, "..", "data", "data1"),
    file.path(root, "..", "data")
  )
  for (candidate in candidates) {
    if (has_region_csvs(candidate)) {
      return(normalizePath(candidate, winslash = "/", mustWork = TRUE))
    }
  }
  file.path(root, "data")
}

analysis_paths <- function(default_pdf_name, options = parse_cli_options()) {
  root <- project_root()
  data_dir <- normalizePath(
    option_value(options, "data_dir", Sys.getenv("PROTEOMICS_DATA_DIR", unset = detect_data_dir(root))),
    winslash = "/",
    mustWork = TRUE
  )
  output_dir <- normalizePath(
    option_value(options, "output_dir", file.path(root, "results")),
    winslash = "/",
    mustWork = FALSE
  )
  output_pdf <- normalizePath(
    option_value(options, "output_pdf", file.path(root, "figures", default_pdf_name)),
    winslash = "/",
    mustWork = FALSE
  )
  font_family <- as.character(option_value(
    options,
    "font_family",
    Sys.getenv("PROTEOMICS_FIGURE_FONT", unset = DEFAULT_FONT_FAMILY)
  ))

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(output_pdf), recursive = TRUE, showWarnings = FALSE)

  list(
    root = root,
    data_dir = data_dir,
    output_dir = output_dir,
    output_pdf = output_pdf,
    font_family = font_family
  )
}

open_figure_pdf <- function(filename, width, height, family = DEFAULT_FONT_FAMILY) {
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  if (capabilities("cairo")) {
    grDevices::cairo_pdf(filename, width = width, height = height, family = family, onefile = TRUE)
  } else {
    fallback_family <- if (family %in% names(grDevices::pdfFonts())) family else "Helvetica"
    warning("Cairo is unavailable; using the standard PDF device with font family: ", fallback_family)
    grDevices::pdf(file = filename, width = width, height = height, family = fallback_family, onefile = TRUE, useDingbats = FALSE)
  }
}

find_region_file <- function(region_code, data_dir) {
  pattern <- paste0("_", region_code, "_DIA_LFQ.*\\.csv$")
  files <- list.files(path = data_dir, pattern = pattern, full.names = TRUE)
  if (length(files) != 1L) {
    stop("Expected exactly one ", region_code, " CSV file in ", data_dir, "; found ", length(files), ".")
  }
  files[[1]]
}

mouse_protein_rows <- function(organisms) {
  !is.na(organisms) & grepl("Mus musculus", organisms, fixed = TRUE)
}

# Match the exact semicolon-delimited iRT accession, retaining mouse proteins
# such as Sirt1 and Teashirt.
is_irt_standard <- function(accessions) {
  vapply(
    strsplit(ifelse(is.na(accessions), "", accessions), ";", fixed = TRUE),
    function(ids) any(trimws(ids) %in% IRT_STANDARD_ACCESSIONS),
    logical(1)
  )
}

safe_row_mean <- function(x) {
  result <- rowMeans(x, na.rm = TRUE)
  result[!is.finite(result)] <- NA_real_
  result
}

first_candidate <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  vapply(strsplit(x, split = ";", fixed = TRUE), function(z) trimws(z[1]), character(1))
}

make_label <- function(gene, accession) {
  label <- as.character(gene)
  empty <- is.na(label) | trimws(label) == ""
  label[empty] <- accession[empty]
  gsub(";", "/", label, fixed = TRUE)
}

row_zscore <- function(mat) {
  row_mean <- rowMeans(mat, na.rm = TRUE)
  row_sd <- apply(mat, 1, sd, na.rm = TRUE)
  z <- sweep(mat, 1, row_mean, "-")
  z <- sweep(z, 1, row_sd, "/")
  z[!is.finite(z)] <- NA_real_
  z
}

read_region_table <- function(region_code, data_dir, required_columns) {
  data <- read.csv(
    file = find_region_file(region_code, data_dir),
    header = TRUE,
    check.names = FALSE,
    na.strings = c("NaN", "", "NA")
  )
  original_columns <- colnames(data)
  if (anyDuplicated(colnames(data))) stop("Duplicate input report columns")
  mapping_dir <- file.path(option_value(parse_cli_options(), "output_dir", file.path(project_root(), "results")), "input_label_mapping")
  dir.create(mapping_dir, recursive=TRUE, showWarnings=FALSE)
  write.csv(data.frame(source_column=original_columns, analysis_column=colnames(data)),
            file.path(mapping_dir, paste0(region_code, ".csv")), row.names=FALSE)
  missing_columns <- setdiff(required_columns, colnames(data))
  if (length(missing_columns) > 0L) {
    stop("Missing required columns in ", region_code, ": ", paste(missing_columns, collapse = ", "))
  }
  data
}

# ---------------------------------------------------------------------------
# Reference-workflow preprocessing for one region report.
# Returns the filtered annotation table plus the normalized log2 matrix
# (proteins x retained samples; exclusion precedes all preprocessing).
# ---------------------------------------------------------------------------
preprocess_region_matrix <- function(region_code, data_dir,
                                     min_detected = MIN_DETECTED_SAMPLES,
                                     include_all = cli_flag(option_value(parse_cli_options(), "include_all", FALSE)),
                                     imputation = option_value(parse_cli_options(), "imputation", "row-median")) {
  required_columns <- c(
    "PG.ProteinAccessions",
    "PG.Genes",
    "PG.Organisms",
    "PG.ProteinDescriptions",
    "PG.NrOfStrippedSequencesIdentified (Experiment-wide)"
  )
  data <- read_region_table(region_code, data_dir, required_columns)

  quantity_cols <- grep(paste0("_G[1-6]_[1-4]_", region_code, ".*\\.raw\\.PG\\.Quantity$"),
                        colnames(data), value = TRUE)
  if (length(quantity_cols) != 24L) {
    stop("Expected 24 biological quantity columns for ", region_code, "; found ", length(quantity_cols), ".")
  }
  validate_region_samples(quantity_cols, region_code)
  policy <- sample_exclusions()
  sample_ids <- clean_sample_name(quantity_cols)
  excluded <- sample_ids %in% policy$sample_id & !include_all
  relevant_policy <- policy$sample_id[endsWith(policy$sample_id, paste0("_", region_code))]
  if (any(!relevant_policy %in% sample_ids)) stop("Policy sample missing in ", region_code)
  manifest <- data.frame(Sample = sample_ids, Included = !excluded,
                         Reason = ifelse(sample_ids %in% policy$sample_id,
                           policy$reason[match(sample_ids, policy$sample_id)], ""))
  quantity_cols <- quantity_cols[!excluded]
  group_labels <- sub(paste0(".*_G([1-6])_[1-4]_", region_code, ".*"), "G\\1", quantity_cols)
  if (!identical(as.character(imputation), "row-median")) stop("Only --imputation row-median is supported")
  if (anyDuplicated(data$PG.ProteinAccessions) || anyNA(data$PG.ProteinAccessions)) stop("Duplicate or missing accessions in ", region_code)
  if (min_detected < 1 || min_detected > length(quantity_cols)) stop("Invalid detection threshold")

  # Step 1: mouse proteins only, iRT standard removed, non-positive -> NA.
  keep_rows <- mouse_protein_rows(data$PG.Organisms) & !is_irt_standard(data$PG.ProteinAccessions)
  data <- data[keep_rows, ]
  n_after_clean <- nrow(data)
  for (column in quantity_cols) {
    data[[column]] <- suppressWarnings(as.numeric(data[[column]]))
    data[[column]][!is.finite(data[[column]]) | data[[column]] <= 0] <- NA_real_
  }

  # Step 2: fixed >=20 detections among retained samples.
  detected <- rowSums(!is.na(as.matrix(data[, quantity_cols, drop = FALSE])))
  data <- data[detected >= min_detected, ]
  observed_counts <- as.matrix(data[, quantity_cols, drop = FALSE])
  n_after_filter <- nrow(data)
  if (n_after_filter == 0L) {
    stop("No proteins passed the ", min_detected, "/", length(quantity_cols), " detection filter for ", region_code, ".")
  }

  # Step 3: log2 and observed-value sample-median normalization.
  log2_matrix <- log2(as.matrix(data[, quantity_cols, drop = FALSE]))
  storage.mode(log2_matrix) <- "double"
  sample_medians <- apply(log2_matrix, 2, median, na.rm = TRUE)
  if (any(!is.finite(sample_medians))) stop("A retained sample has no finite observations")
  normalized <- sweep(log2_matrix, 2, sample_medians, FUN = "-")
  normalized <- sweep(normalized, 2, median(sample_medians), FUN = "+")
  # Step 4: one protein-wise median from observed normalized log2 quantities.
  imputed_count <- sum(is.na(normalized))
  protein_medians <- apply(normalized, 1, median, na.rm = TRUE)
  if (any(!is.finite(protein_medians))) stop("A retained protein group has no observations")
  missing <- which(is.na(normalized), arr.ind = TRUE)
  if (nrow(missing)) normalized[missing] <- protein_medians[missing[, 1]]
  rownames(normalized) <- as.character(data$PG.ProteinAccessions)

  list(
    data = data,
    normalized = normalized,
    quantity_cols = quantity_cols,
    group_labels = group_labels,
    n_input = n_after_clean,
    n_after_filter = n_after_filter,
    imputed_count = imputed_count,
    manifest = manifest,
    imputation = imputation,
    include_all = include_all,
    min_detected = min_detected,
    observed_counts = observed_counts
  )
}

# Generalized two-group differential analysis on the preprocessed matrix.
# The limma model is fitted on all six retained groups (shared variance
# model, as in the reference pipeline); the contrast is group_b - group_a.
run_region_differential <- function(region_code, display_name, data_dir,
                                    group_a = "G1", group_b = "G5",
                                    adjust_method = ADJUST_METHOD) {
  prep <- preprocess_region_matrix(region_code, data_dir)
  data <- prep$data
  normalized <- prep$normalized
  output_dir <- option_value(parse_cli_options(), "output_dir", file.path(project_root(), "results"))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  mode <- if (prep$include_all) "include_all" else "primary"
  write.csv(prep$manifest, file.path(output_dir, paste0(region_code, "_", mode, "_sample_manifest.csv")), row.names = FALSE)
  if (!all(c(group_a, group_b) %in% prep$group_labels) || group_a == group_b) stop("Invalid contrast groups")

  group <- factor(prep$group_labels, levels = paste0("G", 1:6))
  design <- model.matrix(~ 0 + group)
  colnames(design) <- levels(group)

  fit <- limma::lmFit(normalized, design)
  fit <- limma::contrasts.fit(fit, limma::makeContrasts(contrasts = paste(group_b, "-", group_a), levels = design))
  fit <- limma::eBayes(fit, trend = TRUE, robust = TRUE)

  a_columns <- prep$quantity_cols[prep$group_labels == group_a]
  b_columns <- prep$quantity_cols[prep$group_labels == group_b]
  mean_a <- rowMeans(normalized[, a_columns, drop = FALSE])
  mean_b <- rowMeans(normalized[, b_columns, drop = FALSE])

  log2fc <- as.numeric(fit$coefficients[, 1])
  pvalue <- as.numeric(fit$p.value[, 1])
  adj_pvalue <- p.adjust(pvalue, method = adjust_method)

  status <- rep("Not significant", nrow(data))
  status[adj_pvalue < PVALUE_CUTOFF & log2fc > LOG2FC_CUTOFF] <- paste("High in", group_b)
  status[adj_pvalue < PVALUE_CUTOFF & log2fc < -LOG2FC_CUTOFF] <- paste("High in", group_a)

  result <- data.frame(
    Region = display_name,
    RegionFileCode = region_code,
    ProteinAccessions = data$PG.ProteinAccessions,
    Genes = data$PG.Genes,
    Organisms = data$PG.Organisms,
    ProteinDescriptions = data$PG.ProteinDescriptions,
    StrippedSequences = data[["PG.NrOfStrippedSequencesIdentified (Experiment-wide)"]],
    MeanLog2_A = mean_a,
    MeanLog2_B = mean_b,
    Log2FC = log2fc,
    AverageExpression = as.numeric(fit$Amean),
    ModeratedT = as.numeric(fit$t[, 1]),
    PValue = pvalue,
    AdjPValue = adj_pvalue,
    LogOdds = as.numeric(fit$lods[, 1]),
    Status = status,
    Contrast = paste(group_b, "-", group_a),
    AnalysisSampleMode = mode,
    Samples_A = length(a_columns),
    Samples_B = length(b_columns),
    ModelSamples = ncol(normalized),
    DetectionThreshold = prep$min_detected,
    Imputation = "row-median; normalized log2 scale",
    Observed_A = rowSums(is.finite(prep$observed_counts[, a_columns, drop = FALSE])),
    Observed_B = rowSums(is.finite(prep$observed_counts[, b_columns, drop = FALSE])),
    stringsAsFactors = FALSE
  )
  result <- result[order(result$AdjPValue, -abs(result$Log2FC), na.last = TRUE), ]

  mean_table <- data.frame(
    ProteinAccessions = data$PG.ProteinAccessions,
    Genes = data$PG.Genes,
    ProteinDescriptions = data$PG.ProteinDescriptions,
    Region = display_name,
    stringsAsFactors = FALSE
  )
  mean_table[[group_a]] <- mean_a
  mean_table[[group_b]] <- mean_b

  list(
    result = result,
    mean_table = mean_table,
    group_a = group_a,
    group_b = group_b,
    qc = data.frame(
      Region = display_name,
      ProteinsAfterCleaning = prep$n_input,
      ProteinsAfterDetectionFilter = prep$n_after_filter,
      MedianImputedCells = prep$imputed_count,
      TestedProteins = nrow(result),
      SignificantProteins = sum(status != "Not significant"),
      stringsAsFactors = FALSE
    )
  )
}
