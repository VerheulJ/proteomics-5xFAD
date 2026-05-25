# =============================================================================
# 02_differential_expression.R
# Proteomics Analysis — 5xFAD Mouse Model
#
# Key difference vs 3xTg pipeline:
#   - p-value threshold: 0.01 (not 0.05)
#   - FC threshold: 50% change  →  log2(1.5) ≈ 0.58
#
# Input:
#   - data/processed/ajustado.xlsx
#   - data/processed/entera_5xfad.xlsx
#   - data/raw/sample-metadata.xlsx
#
# Output:
#   - results/tables/vulcanot_5xfad.xlsx
#   - results/tables/df_desreguladas_5xfad.xlsx
# =============================================================================

# --- Libraries ----------------------------------------------------------------
library(readxl)
library(openxlsx)
library(dplyr)
library(tibble)
library(effsize)

options(scipen = 999)

# --- Load data ----------------------------------------------------------------
entera_experimento_5xfad <- read_excel("data/processed/entera_5xfad.xlsx")
sample_metadata           <- read_xlsx("data/raw/sample-metadata.xlsx")
ajustado_raw              <- read_xlsx("data/processed/ajustado.xlsx")

accesion_col    <- dplyr::select(entera_experimento_5xfad, Accession)
gene_symbol_col <- dplyr::select(entera_experimento_5xfad, `Gene Symbol`)
indice_proteina <- cbind(accesion_col, gene_symbol_col)

ajustado <- as.data.frame(ajustado_raw)
rownames(ajustado) <- ajustado_raw$...1
ajustado$...1 <- NULL

# --- Separate samples by genotype ---------------------------------------------
prot_control <- sample_metadata %>% filter(Genotype == "WT")    %>% pull(`Sample-id`)
prot_5xfad   <- sample_metadata %>% filter(Genotype == "5xFAD") %>% pull(`Sample-id`)

tabla_control <- ajustado %>% select(all_of(prot_control))
tabla_5xfad   <- ajustado %>% select(all_of(prot_5xfad))

tabla_control_datos       <- tabla_control
tabla_control_datos$media <- rowSums(tabla_control, na.rm = TRUE) / ncol(tabla_control)

tabla_5xfad_datos       <- tabla_5xfad
tabla_5xfad_datos$media <- rowSums(tabla_5xfad, na.rm = TRUE) / ncol(tabla_5xfad)

nprot <- nrow(tabla_control)

# --- Statistical testing (protein-by-protein) ---------------------------------
p_5xfad         <- numeric(nprot)
normalidad_5xfad <- numeric(nprot)
normalidad_ctrl  <- numeric(nprot)

for (n in seq_len(nprot)) {
  c_5xfad  <- as.numeric(tabla_5xfad[n, ])
  c_control <- as.numeric(tabla_control[n, ])

  p_norm_5xfad <- if (length(unique(c_5xfad))   == 1) 1 else shapiro.test(c_5xfad)$p.value
  p_norm_ctrl  <- if (length(unique(c_control)) == 1) 1 else shapiro.test(c_control)$p.value

  normalidad_5xfad[n] <- p_norm_5xfad
  normalidad_ctrl[n]  <- p_norm_ctrl

  if (all(c_5xfad == c_control)) {
    p_5xfad[n] <- 1
  } else if (p_norm_5xfad > 0.05 && p_norm_ctrl > 0.05) {
    p_5xfad[n] <- t.test(c_5xfad, c_control, var.equal = FALSE)$p.value
  } else {
    p_5xfad[n] <- wilcox.test(c_5xfad, c_control)$p.value
  }
}

# --- Assemble results table ---------------------------------------------------
gene_symbols <- indice_proteina$`Gene Symbol`[
  match(rownames(tabla_5xfad), indice_proteina$Accession)
]

result_5xfad <- data.frame(
  Protein            = rownames(tabla_5xfad),
  Gene               = gene_symbols,
  P_Value            = p_5xfad,
  mean_control       = tabla_control_datos$media,
  mean_5xfad         = tabla_5xfad_datos$media,
  normalidad_5xfad   = normalidad_5xfad,
  normalidad_control = normalidad_ctrl
)

# --- Fold change and -log10 p-value -------------------------------------------
log2FC         <- log2(tabla_5xfad_datos$media / tabla_control_datos$media)
neg_log10_pval <- -log10(result_5xfad$P_Value)

descriptions_selected <- entera_experimento_5xfad$Description[
  match(result_5xfad$Protein, entera_experimento_5xfad$Accession)
]

vulcanot_5xfad <- data.frame(
  Accession          = result_5xfad$Protein,
  Gene.Symbol        = gene_symbols,
  log2cociente_5xfad = log2FC,
  log10pvalor_5xfad  = neg_log10_pval,
  P_Value            = result_5xfad$P_Value,
  Normality_5xfad    = normalidad_5xfad,
  Normality_control  = normalidad_ctrl,
  mean_5xfad         = tabla_5xfad_datos$media,
  mean_control       = tabla_control_datos$media,
  Description        = descriptions_selected
)

# --- Classification thresholds ------------------------------------------------
# NOTE: stricter than 3xTg pipeline — p < 0.01 and FC > 50%
p_umbral <- -log10(0.01)   # = 2
FC        <- 0.5
LFCP      <- log2(1 + FC)  # ≈  0.585
LFCN      <- log2(1 - FC)  # ≈ -1.000

up_regulated   <- filter(vulcanot_5xfad, log2cociente_5xfad >  LFCP & log10pvalor_5xfad > p_umbral)
down_regulated <- filter(vulcanot_5xfad, log2cociente_5xfad <  LFCN & log10pvalor_5xfad > p_umbral)

up_regulated$Status   <- "up-regulated"
down_regulated$Status <- "down-regulated"

df_desreguladas_5xfad <- rbind(down_regulated, up_regulated)

# --- Save outputs -------------------------------------------------------------
dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)
write.xlsx(vulcanot_5xfad,        "results/tables/vulcanot_5xfad.xlsx",        sheetName = "todas")
write.xlsx(df_desreguladas_5xfad, "results/tables/df_desreguladas_5xfad.xlsx", sheetName = "Desreguladas")

message(sprintf(
  "Differential expression complete.\n  Up-regulated:   %d\n  Down-regulated: %d",
  nrow(up_regulated), nrow(down_regulated)
))
