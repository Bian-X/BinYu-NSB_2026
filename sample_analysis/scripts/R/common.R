options(stringsAsFactors = FALSE)

DEFAULT_FONT_FAMILY <- "Arial"
PVALUE_CUTOFF <- 0.05
FC_CUTOFF <- 1.5
LOG2FC_CUTOFF <- log2(FC_CUTOFF)
MIN_VALID_NUMBER <- 3L

REGION_CODES <- c("PMD", "LPB", "VHPC", "PVH", "BLA", "DHPC", "CEA", "MSC", "SSC", "IL", "PL")

COLOR_HIGH_G4 <- "#D94B3D"
COLOR_HIGH_G2 <- "#2878B5"
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

project_root <- function(script_dir = script_directory()) {
  normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = TRUE)
}

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
    warning("Cairo is unavailable; using the standard PDF device.")
    grDevices::pdf(file = filename, width = width, height = height, family = family, onefile = TRUE, useDingbats = FALSE)
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
  missing_columns <- setdiff(required_columns, colnames(data))
  if (length(missing_columns) > 0L) {
    stop("Missing required columns in ", region_code, ": ", paste(missing_columns, collapse = ", "))
  }
  data
}

run_region_differential <- function(region_code, display_name, data_dir) {
  required_columns <- c(
    "PG.ProteinAccessions",
    "PG.Genes",
    "PG.Organisms",
    "PG.ProteinDescriptions",
    "PG.NrOfStrippedSequencesIdentified (Experiment-wide)"
  )
  data <- read_region_table(region_code, data_dir, required_columns)

  g2_columns <- grep(paste0("_G2_[1-4]_", region_code, ".*\\.raw\\.PG\\.Quantity$"), colnames(data), value = TRUE)
  g4_columns <- grep(paste0("_G4_[1-4]_", region_code, ".*\\.raw\\.PG\\.Quantity$"), colnames(data), value = TRUE)
  if (length(g2_columns) != 4L || length(g4_columns) != 4L) {
    stop(
      "Expected four G2 and four G4 samples for ", region_code,
      "; found G2=", length(g2_columns), ", G4=", length(g4_columns), "."
    )
  }

  sample_columns <- c(g2_columns, g4_columns)
  data <- data[mouse_protein_rows(data$PG.Organisms), ]

  for (column in sample_columns) {
    data[[column]] <- suppressWarnings(as.numeric(data[[column]]))
    data[[column]][data[[column]] <= 0] <- NA_real_
  }

  exp_log2 <- log2(as.matrix(data[, sample_columns]))
  storage.mode(exp_log2) <- "double"

  g2_exp <- exp_log2[, g2_columns, drop = FALSE]
  g4_exp <- exp_log2[, g4_columns, drop = FALSE]
  n_g2 <- rowSums(!is.na(g2_exp))
  n_g4 <- rowSums(!is.na(g4_exp))

  mean_g2 <- safe_row_mean(g2_exp)
  mean_g4 <- safe_row_mean(g4_exp)
  log2fc <- mean_g4 - mean_g2
  keep <- n_g2 >= MIN_VALID_NUMBER & n_g4 >= MIN_VALID_NUMBER

  group <- factor(c(rep("G2", length(g2_columns)), rep("G4", length(g4_columns))), levels = c("G2", "G4"))
  design <- model.matrix(~ 0 + group)
  colnames(design) <- levels(group)

  fit <- limma::lmFit(exp_log2[keep, , drop = FALSE], design)
  fit <- limma::contrasts.fit(fit, limma::makeContrasts(G4 - G2, levels = design))
  fit <- limma::eBayes(fit, trend = TRUE)

  pvalue <- rep(NA_real_, nrow(data))
  t_value <- rep(NA_real_, nrow(data))
  ave_expr <- rep(NA_real_, nrow(data))
  log_odds <- rep(NA_real_, nrow(data))

  log2fc[keep] <- fit$coefficients[, 1]
  pvalue[keep] <- fit$p.value[, 1]
  t_value[keep] <- fit$t[, 1]
  ave_expr[keep] <- fit$Amean
  log_odds[keep] <- fit$lods[, 1]

  status <- rep("Not tested", nrow(data))
  status[is.finite(pvalue)] <- "Not significant"
  status[is.finite(pvalue) & pvalue < PVALUE_CUTOFF & log2fc > LOG2FC_CUTOFF] <- "High in G4"
  status[is.finite(pvalue) & pvalue < PVALUE_CUTOFF & log2fc < -LOG2FC_CUTOFF] <- "High in G2"

  result <- data.frame(
    Region = display_name,
    RegionFileCode = region_code,
    ProteinAccessions = data$PG.ProteinAccessions,
    Genes = data$PG.Genes,
    Organisms = data$PG.Organisms,
    ProteinDescriptions = data$PG.ProteinDescriptions,
    StrippedSequences = data[["PG.NrOfStrippedSequencesIdentified (Experiment-wide)"]],
    N_G2 = n_g2,
    N_G4 = n_g4,
    MeanLog2_G2 = mean_g2,
    MeanLog2_G4 = mean_g4,
    Log2FC_G4_vs_G2 = log2fc,
    AverageExpression = ave_expr,
    ModeratedT = t_value,
    PValue = pvalue,
    LogOdds = log_odds,
    Status = status,
    stringsAsFactors = FALSE
  )
  result <- result[order(result$PValue, -abs(result$Log2FC_G4_vs_G2), na.last = TRUE), ]

  mean_table <- data.frame(
    ProteinAccessions = data$PG.ProteinAccessions,
    Genes = data$PG.Genes,
    ProteinDescriptions = data$PG.ProteinDescriptions,
    Region = display_name,
    G2 = mean_g2,
    G4 = mean_g4,
    stringsAsFactors = FALSE
  )

  list(result = result, mean_table = mean_table)
}
