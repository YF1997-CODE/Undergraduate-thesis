# =============================================================================
# GTEx ASE Pipeline: CYP3A Gene Family Allelic Imbalance Analysis
# =============================================================================
# Analyses allele-specific expression (ASE) data from GTEx v8 for three
# CYP3A cytochrome P450 genes (CYP3A4 / CYP3A5 / CYP3A7) across liver and
# small intestine, identifies cis-regulatory variants via chi-square testing,
# and produces publication-quality figures.
#
# Dependencies: tidyverse, ggplot2, patchwork, rvest, xml2, dplyr
# Author : <your name>
# License: MIT
# =============================================================================

library(tidyverse)
library(ggplot2)
library(patchwork)
library(rvest)
library(xml2)

# -----------------------------------------------------------------------------
# 0.  Configuration
# -----------------------------------------------------------------------------
DATA_DIR  <- "D:/data"
ASE_DIR   <- file.path(DATA_DIR, "GTEx_Analysis_v8_ASE_WASP_counts_by_subject")
OUT_DIR   <- file.path(DATA_DIR, "results")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

GENES <- c(
  CYP3A4 = "ENSG00000160868",
  CYP3A5 = "ENSG00000106258",
  CYP3A7 = "ENSG00000160870"
)

TISSUES <- c("LIVER", "SNTTRM")   # GTEx tissue codes

VARIANT_KEEP <- c(
  "stop_lost", "synonymous_variant",
  "3_prime_UTR_variant", "5_prime_UTR_variant",
  "missense_variant"
)

# Genomic coordinates for gene-body annotations (hg38)
GENE_COORDS <- list(
  CYP3A4 = list(start = 99756960L, end = 99784248L, mid = 99770604L),
  CYP3A5 = list(start = 99648194L, end = 99679998L, mid = 99664096L),
  CYP3A7 = list(start = 99684957L, end = 99735196L, mid = 99710077L)
)

# Significance threshold
ALPHA <- 0.05

# -----------------------------------------------------------------------------
# 1.  Load & merge GTEx ASE files
# -----------------------------------------------------------------------------
load_ase_data <- function(ase_dir, genes, tissues) {
  files <- list.files(ase_dir, full.names = TRUE)
  message(sprintf("Found %d ASE subject files.", length(files)))

  bind_rows(lapply(files, function(f) {
    df <- read_csv(f, show_col_types = FALSE)
    df |>
      filter(gene_id %in% genes, tissue %in% tissues)
  }))
}

ase_raw <- load_ase_data(ASE_DIR, GENES, TISSUES)

# Harmonise column names to snake_case
ase <- ase_raw |>
  rename_with(tolower) |>
  rename(
    gene_id    = gene_id,
    ref_count  = ref_count,
    alt_count  = alt_count,
    binom_p    = binom_p,
    variant_annotation = variant_annotation
  ) |>
  mutate(
    gene_name  = names(GENES)[match(gene_id, GENES)],
    sig_label  = if_else(binom_p < ALPHA, "<0.05", "≥0.05"),
    sig_label  = factor(sig_label, levels = c("<0.05", "≥0.05"))
  )

message(sprintf("Retained %d ASE records across %d tissues.", nrow(ase), n_distinct(ase$tissue)))

# -----------------------------------------------------------------------------
# 2.  Helper: subset by tissue × gene
# -----------------------------------------------------------------------------
subset_ase <- function(data, tissue_code, gene_ens) {
  data |> filter(tissue == tissue_code, gene_id == gene_ens)
}

# -----------------------------------------------------------------------------
# 3.  Filter functional variants (potential cis-markers)
# -----------------------------------------------------------------------------
ase_functional <- ase |>
  filter(variant_annotation %in% VARIANT_KEEP)

# Split by gene for downstream plotting
func_by_gene <- split(ase_functional, ase_functional$gene_name)

# -----------------------------------------------------------------------------
# 4.  Visualisation helpers
# -----------------------------------------------------------------------------
PALETTE_SIG <- c("<0.05" = "#E53935", "≥0.05" = "#1E88E5")

#' Scatter plot of REF vs ALT counts, coloured by Binomial P significance.
ase_scatter <- function(df, gene_label, tissue_label, add_lm = TRUE) {
  p <- ggplot(df, aes(x = ref_count, y = alt_count, colour = sig_label)) +
    geom_point(alpha = 0.65, size = 1.8) +
    scale_colour_manual(values = PALETTE_SIG, name = "Binom. P") +
    labs(
      title   = sprintf("%s — %s", gene_label, tissue_label),
      x       = "REF allele count",
      y       = "ALT allele count"
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "right")

  if (add_lm && nrow(df) > 2) {
    p <- p + geom_smooth(
      aes(group = 1), method = "lm", se = TRUE,
      colour = "grey40", linewidth = 0.7, linetype = "dashed"
    )
  }
  p
}

#' Manhattan-style plot: –log10(p) vs genomic position with gene annotation.
manhattan_plot <- function(df, gene_name, coords, colour_var) {
  g <- coords[[gene_name]]

  df |>
    filter(!is.na(pval)) |>
    mutate(logp = -log10(pval)) |>
    ggplot(aes(x = pos, y = logp, colour = .data[[colour_var]])) +
    geom_point(alpha = 0.7, size = 1.6) +
    # gene body arrow
    geom_segment(
      aes(x = g$start, xend = g$end, y = -0.08, yend = -0.08),
      colour = "black", linewidth = 1.1, inherit.aes = FALSE,
      arrow = arrow(type = "closed", length = unit(0.12, "cm"))
    ) +
    annotate(
      "text", x = g$mid, y = -log10(4) + 0.5,
      label = gene_name, fontface = "bold", size = 2.8, colour = "black"
    ) +
    scale_colour_manual(values = PALETTE_SIG) +
    labs(x = "Genomic position (hg38)", y = expression(-log[10](p))) +
    theme_bw(base_size = 11) +
    theme(
      axis.line   = element_line(arrow = arrow(length = unit(0.4, "cm"))),
      legend.position = "right"
    )
}

# -----------------------------------------------------------------------------
# 5.  Section A — All ASE scatter plots (liver & small intestine × 3 genes)
# -----------------------------------------------------------------------------
scatter_panels <- list()
for (tissue in TISSUES) {
  tissue_label <- if (tissue == "LIVER") "Liver" else "Small intestine"
  for (gname in names(GENES)) {
    key <- paste(tissue, gname, sep = "_")
    df  <- subset_ase(ase, tissue, GENES[[gname]])
    scatter_panels[[key]] <- ase_scatter(df, gname, tissue_label)
  }
}

# 2-row × 3-col panel (liver top, intestine bottom)
panel_scatter <- (
  scatter_panels[["LIVER_CYP3A4"]]  |
  scatter_panels[["LIVER_CYP3A5"]]  |
  scatter_panels[["LIVER_CYP3A7"]]
) / (
  scatter_panels[["SNTTRM_CYP3A4"]] |
  scatter_panels[["SNTTRM_CYP3A5"]] |
  scatter_panels[["SNTTRM_CYP3A7"]]
) + plot_annotation(
  title    = "CYP3A allele-specific expression: REF vs ALT counts",
  subtitle = "Red = Binomial P < 0.05  |  Blue = P ≥ 0.05",
  theme    = theme(plot.title = element_text(face = "bold", size = 13))
)

ggsave(file.path(OUT_DIR, "Fig1_ASE_scatter_all.pdf"),
       panel_scatter, width = 14, height = 9)

# -----------------------------------------------------------------------------
# 6.  Section B — Functional-variant scatter plots
# -----------------------------------------------------------------------------
func_panels <- lapply(names(GENES), function(gname) {
  df <- func_by_gene[[gname]]
  if (is.null(df) || nrow(df) == 0) return(NULL)
  ase_scatter(df, gname, "functional variants", add_lm = FALSE)
})
names(func_panels) <- names(GENES)

panel_func <- (func_panels[["CYP3A4"]] | func_panels[["CYP3A5"]] | func_panels[["CYP3A7"]]) +
  plot_annotation(
    title = "Functional variants (stop_lost, missense, UTR, synonymous)",
    theme = theme(plot.title = element_text(face = "bold", size = 13))
  )

ggsave(file.path(OUT_DIR, "Fig2_functional_variants.pdf"),
       panel_func, width = 14, height = 5)

# -----------------------------------------------------------------------------
# 7.  Section C — Linear model: ALT ~ REF (per gene × tissue)
# -----------------------------------------------------------------------------
run_lm <- function(data, tissue_code, gene_ens) {
  df <- subset_ase(data, tissue_code, gene_ens)
  if (nrow(df) < 5) return(NULL)
  fit <- lm(alt_count ~ ref_count, data = df)
  broom::tidy(fit) |> mutate(gene_id = gene_ens, tissue = tissue_code)
}

lm_results <- bind_rows(lapply(TISSUES, function(t) {
  lapply(GENES, function(g) run_lm(ase, t, g))
}))

write_csv(lm_results, file.path(OUT_DIR, "lm_results.csv"))
message("Linear model results saved.")

# -----------------------------------------------------------------------------
# 8.  Section D — Chi-square test: AEI status vs variant genotype
# -----------------------------------------------------------------------------
#' For each row (variant) in a genotype matrix, test whether the proportion
#' of AEI samples differs between HET and HOM carriers.
#'
#' @param geno_matrix  Matrix: rows = variants, cols = samples.
#'                     Last row must be AEI status ("AEI" / "NON").
#'                     Cells are genotype calls ("HET" / "HOM").
#' @param pos_vector   Numeric vector of genomic positions (length = nrow - 1).
#' @return             Data frame with columns: pos, chi_sq_p.

chisq_ase_scan <- function(geno_matrix, pos_vector) {
  n_variants <- nrow(geno_matrix) - 1L
  aei_row    <- as.character(geno_matrix[nrow(geno_matrix), ])
  aei_factor <- factor(aei_row, levels = c("AEI", "NON"))

  p_values <- vapply(seq_len(n_variants), function(i) {
    geno_factor <- factor(as.character(geno_matrix[i, ]),
                          levels = c("HET", "HOM"))
    tbl <- table(data.frame(aei = aei_factor, geno = geno_factor))
    if (any(dim(tbl) < 2)) return(NA_real_)
    suppressWarnings(chisq.test(tbl)$p.value)
  }, numeric(1))

  data.frame(pos = pos_vector, chi_sq_p = p_values)
}

# Example usage (assumes ud3a4, ud3a5, ud3a7 are pre-loaded genotype matrices
# and ud3a4$POS, ud3a5$POS, ud3a7$POS are the position vectors):
#
#   result_3a4 <- chisq_ase_scan(as.matrix(ud3a4[, -1]), ud3a4$POS)
#   result_3a5 <- chisq_ase_scan(as.matrix(ud3a5[, -1]), ud3a5$POS)
#   result_3a7 <- chisq_ase_scan(as.matrix(ud3a7[, -1]), ud3a7$POS)

# -----------------------------------------------------------------------------
# 9.  Section E — Manhattan plots (cis-regulatory loci)
# -----------------------------------------------------------------------------
# Read pre-processed chi-square results (exported from Section 8)
read_chisq_result <- function(path, tissue_label) {
  read_csv(path, show_col_types = FALSE) |>
    filter(!is.na(pval)) |>
    mutate(
      logp    = -log10(pval),
      tissue  = tissue_label,
      sig_label = if_else(pval < ALPHA, "<0.05", "≥0.05"),
      sig_label = factor(sig_label, levels = c("<0.05", "≥0.05"))
    )
}

graph4 <- read_chisq_result(file.path(DATA_DIR, "4.csv"), "Liver CYP3A4")
graph5 <- read_chisq_result(file.path(DATA_DIR, "5.csv"), "Liver CYP3A5")
graph7 <- read_chisq_result(file.path(DATA_DIR, "7.csv"), "Intestine CYP3A7")

p_manhattan_3a4 <- manhattan_plot(graph4, "CYP3A4", GENE_COORDS, "sig_label")
p_manhattan_3a5 <- manhattan_plot(graph5, "CYP3A5", GENE_COORDS, "sig_label")
p_manhattan_3a7 <- manhattan_plot(graph7, "CYP3A7", GENE_COORDS, "sig_label")

panel_manhattan <- (p_manhattan_3a4 / p_manhattan_3a5 / p_manhattan_3a7) +
  plot_annotation(
    title    = "Cis-regulatory variant scan — CYP3A gene family",
    subtitle = "Chi-square test: AEI status vs HET/HOM genotype",
    theme    = theme(plot.title = element_text(face = "bold", size = 13))
  )

ggsave(file.path(OUT_DIR, "Fig3_manhattan_cis.pdf"),
       panel_manhattan, width = 10, height = 11)

# -----------------------------------------------------------------------------
# 10. Section F — Web scraping: SNP annotation from NCBI dbSNP
# -----------------------------------------------------------------------------
#' Fetch functional annotation for a single rsID from dbSNP.
#' Returns a one-row tibble; NA columns on failure.
fetch_snp_annotation <- function(rsid, base_url = "https://www.ncbi.nlm.nih.gov/snp/") {
  url <- paste0(base_url, rsid)
  Sys.sleep(0.5)   # be polite to NCBI

  tryCatch({
    page  <- read_html(url, encoding = "utf-8")
    nodes <- page |>
      html_elements("#main_content main div.summary-box.usa-grid-full dl:nth-child(2)")

    clinical_sig <- nodes |> html_element("dd:nth-child(4) div") |> html_text(trim = TRUE)
    gene_info    <- nodes |> html_element("dd:nth-child(6)")      |> html_text(trim = TRUE)

    tibble(rsid = rsid, clinical_sig = clinical_sig, gene_info = gene_info)
  }, error = function(e) {
    warning(sprintf("Failed to scrape %s: %s", rsid, conditionMessage(e)))
    tibble(rsid = rsid, clinical_sig = NA_character_, gene_info = NA_character_)
  })
}

#' Batch-scrape a vector of rsIDs, return combined tibble.
fetch_snp_batch <- function(rsids) {
  bind_rows(lapply(rsids, fetch_snp_annotation))
}

# Usage:
#   snp_annotations <- fetch_snp_batch(data5n$snp_id)
#   # Prioritise SNPs with citation count >= 2
#   high_priority <- snp_annotations |> filter(cit >= 2) |> arrange(desc(cit))

# -----------------------------------------------------------------------------
# 11. Session info (reproducibility)
# -----------------------------------------------------------------------------
writeLines(capture.output(sessionInfo()),
           file.path(OUT_DIR, "session_info.txt"))
message("\n✓ Pipeline complete. Outputs written to: ", OUT_DIR)
