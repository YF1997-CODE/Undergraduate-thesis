# CYP3A ASE Pipeline

> Allele-specific expression analysis and *cis*-eQTL scanning for *CYP3A4*, *CYP3A5*, and *CYP3A7*  
> across liver and small intestine using GTEx v8 WASP-filtered ASE data.

![R](https://img.shields.io/badge/R-4.x-276DC3?style=flat-square&logo=r)
![License](https://img.shields.io/badge/license-MIT-green?style=flat-square)
![Genome](https://img.shields.io/badge/genome-hg38-orange?style=flat-square)
![GTEx](https://img.shields.io/badge/data-GTEx%20v8-blueviolet?style=flat-square)

---

## Background

The *CYP3A* subfamily accounts for ~25% of total hepatic CYP450 protein and metabolises approximately 60% of clinically used drugs. Inter-individual expression variability exceeds 30-fold, driven in part by *cis*-acting regulatory SNPs. This pipeline identifies candidate *cis*-eQTLs by detecting allelic expression imbalance (AEI) in heterozygous GTEx subjects, then running a chi-square association scan across ±250 kb of each gene's transcription start site.

---

## Target genes

| Gene | Ensembl ID | TSS (hg38) | Markers | Significant positions |
|------|-----------|-----------|---------|----------------------|
| *CYP3A4* | ENSG00000160868 | chr7:99,756,960 | 17 | 16 / 6,747 |
| *CYP3A5* | ENSG00000106258 | chr7:99,648,194 | 109 | 467 / 6,747 |
| *CYP3A7* | ENSG00000160870 | chr7:99,684,957 | 116 | 42 / 6,726 |

---

## Pipeline overview

```
GTEx ASE files (×838)
        │
        ▼
01  load_ase_data()        Filter to CYP3A genes × LIVER / SNTTRM
        │
        ▼
02  AEI classification     adjusted binomial p < 0.05 → AEI
                           CDS / 3′-UTR / 5′-UTR marker SNPs
        │
        ▼
03  run_lm()               ALT ~ REF linear model per gene × tissue
                           Larger p → greater allelic scatter → more cis-eQTL activity
        │
        ▼
04  chisq_ase_scan()       χ² test: HET/HOM × AEI/NON
                           All variants within ±250 kb of TSS
        │
        ▼
05  fetch_snp_batch()      Rate-limited dbSNP scrape → rsID annotation
        │
        ▼
    Figures + TSV/CSV outputs
```

---

## Installation

```r
install.packages(c(
  "tidyverse",
  "ggplot2",
  "patchwork",
  "rvest",
  "xml2",
  "broom"
))
```

R ≥ 4.0 is required. No Bioconductor packages are needed.

---

## Quick start

```r
# 1. Edit the paths at the top of the script
DATA_DIR <- "path/to/your/data"
ASE_DIR  <- file.path(DATA_DIR, "GTEx_Analysis_v8_ASE_WASP_counts_by_subject")
OUT_DIR  <- file.path(DATA_DIR, "results")

# 2. Run
source("cyp3a_ase_pipeline.R")
```

All outputs are written to `OUT_DIR/` automatically.

---

## Configuration

All tuneable parameters live at the top of `cyp3a_ase_pipeline.R`:

```r
GENES <- c(
  CYP3A4 = "ENSG00000160868",
  CYP3A5 = "ENSG00000106258",
  CYP3A7 = "ENSG00000160870"
)

TISSUES <- c("LIVER", "SNTTRM")

VARIANT_KEEP <- c(
  "stop_lost", "synonymous_variant",
  "3_prime_UTR_variant", "5_prime_UTR_variant", "missense_variant"
)

ALPHA <- 0.05   # AEI significance threshold (adjusted binomial p)
```

---

## Data requirements

```
data/
├── GTEx_Analysis_v8_ASE_WASP_counts_by_subject/
│   ├── GTEx-XXXX-LIVER.csv        # × 838 subject files
│   └── ...
├── 4.csv    # pre-computed chi-square results for CYP3A4
├── 5.csv    # pre-computed chi-square results for CYP3A5
└── 7.csv    # pre-computed chi-square results for CYP3A7
```

GTEx protected-access data requires dbGaP approval (`phs000424`).  
Genotype VCF files are accessed separately via the GTEx portal.

---

## Key functions

### `load_ase_data(ase_dir, genes, tissues)`

Reads all 838 subject-level ASE CSV files via `lapply` + `bind_rows`,
filtering to target gene IDs and tissue codes. Returns a single tidy
data frame in snake_case column names.

### `ase_scatter(df, gene_label, tissue_label, add_lm = TRUE)`

Scatter plot of REF vs ALT allele counts. Points are coloured by
binomial p significance (red = AEI, blue = non-AEI). Optionally adds
a dashed linear smoother. Returns a `ggplot` object for use with
`patchwork`.

### `manhattan_plot(df, gene_name, coords, colour_var)`

Plots −log₁₀(p) against hg38 genomic position. Annotates the gene body
with a directional arrow and bold label at the midpoint coordinate
defined in `GENE_COORDS`.

### `chisq_ase_scan(geno_matrix, pos_vector)`

Vectorised chi-square scan using `vapply`. Each row of `geno_matrix`
is a variant; the final row encodes AEI/NON status for the same
samples. Returns a data frame of `(pos, chi_sq_p)` pairs.

```r
# Example
result_3a5 <- chisq_ase_scan(
  geno_matrix = as.matrix(ud3a5[, -1]),
  pos_vector  = ud3a5$POS
)
```

### `run_lm(data, tissue_code, gene_ens)`

Fits `ALT_COUNT ~ REF_COUNT` via `lm()` and returns tidy coefficients
from `broom::tidy()`. A larger regression p-value indicates greater
scatter from cis-eQTL-driven allelic imbalance.

### `fetch_snp_batch(rsids)`

Wraps `fetch_snp_annotation()` in `lapply`. Each call to the inner
function reads the NCBI dbSNP page for one rsID using `rvest`,
extracts clinical significance and gene annotation fields, and sleeps
0.5 s between requests. Failures are caught by `tryCatch` and returned
as `NA` rows rather than stopping the batch.

```r
snp_info <- fetch_snp_batch(c("rs780822", "rs35987562", "rs10953293"))
```

---

## Outputs

| File | Description |
|------|-------------|
| `hits.tsv` | Every marker position passing the AEI cutoff |
| `lm_results.csv` | Linear model coefficients per gene × tissue |
| `Fig1_ASE_scatter_all.pdf` | 2 × 3 REF vs ALT panel, both tissues |
| `Fig2_functional_variants.pdf` | Missense / UTR / stop-lost marker scatter |
| `Fig3_manhattan_cis.pdf` | −log₁₀(p) vs hg38 position, stacked panel |
| `session_info.txt` | Full R session info for reproducibility |

---

## Top candidate *cis*-eQTLs

### *CYP3A4* — top signal p = 0.0178

| rsID | Position (hg38) | Locus |
|------|----------------|-------|
| rs35987562 | 99,855,044 | *CYP3A43* |
| rs517284 | 99,857,132 | intergenic |
| rs480596 | 99,858,788 | intergenic |
| rs680055 | 99,859,982 | intergenic |
| rs10278040 | 99,543,750 | upstream |

### *CYP3A5* — top signal p = 2.04 × 10⁻⁸

| rsID | Position (hg38) | Locus |
|------|----------------|-------|
| rs780822 | 99,622,460 | *ZSCAN5* |
| rs13362853 | 99,625,524 | *ZSCAN5* |
| rs6859590 | 99,629,549 | *ZSCAN5* |
| rs10229552 | 99,637,276 | *ZSCAN5* |
| rs7780328 | 99,659,221 | *CYP3A5* / *ZSCAN5* |

### *CYP3A7* — top signal p = 0.0026

| rsID | Position (hg38) | Locus |
|------|----------------|-------|
| rs10953293 | 99,981,409 | *AZGP1P1* |
| rs117268080 | 99,955,443 | intergenic |
| rs3843540 | 99,529,017 | *ZKSCAN5* |
| rs10264022 | 99,610,679 | *TMEM225B* |
| rs522415 | 99,853,750 | *CYP3A43* |

---

## Methods summary

**AEI classification** uses the GTEx adjusted binomial p-value
(`BINOM_P_ADJUSTED`). Samples with REF ratio ∈ [0.43, 0.56] are
consistently non-significant under this criterion, validating it over
the nominal p-value.

**Linear model comparison** interprets a larger regression p-value as
evidence of cis-regulatory AEI: a *cis*-eQTL shifts a subset of
heterozygotes onto an alternate REF–ALT slope, degrading the goodness
of fit of a single-line model.

**Chi-square scan** tests 2 × 2 contingency tables (HET/HOM ×
AEI/NON) at every biallelic variant within ±250 kb of the TSS.
Positions with χ² p < 0.05 are reported as candidate *cis*-eQTL sites.

**Web scraping** uses `rvest` / `xml2` to query NCBI dbSNP. A
`Sys.sleep(0.5)` call between requests respects NCBI rate limits.
SNPs with ≥ 2 citations are flagged for priority follow-up.

---

## Citation

If you use this pipeline in your research, please cite:

```
[Author names]. Identification of cis-eQTLs regulating allelic
expression imbalance of human CYP3A4, CYP3A5, and CYP3A7 using
GTEx allele-specific expression data. [Journal], [Year].
```

---

## License

MIT © [Author]

---

## Acknowledgements

GTEx data: dbGaP `phs000424`.  
*CYP3A* coordinate reference: Ensembl release 109, hg38.
