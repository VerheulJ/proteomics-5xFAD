# =============================================================================
# 01_preprocessing.R
# Proteomics Analysis — 5xFAD Mouse Model
#
# Description:
#   Loads raw LFQ abundance data (including isoforms), cleans gene symbols,
#   removes outliers, filters low-coverage proteins, imputes missing values,
#   and normalizes abundances by a sample correction factor.
#
# Key difference vs 3xTg pipeline:
#   - Coverage filter: >=4 valid values for 5xFAD males, >=2 for 5xFAD females
#   - Extended gene symbol dictionary (21 entries)
#
# Input:
#   - data/raw/230223-LFQ_Incl-Isoforms_FILTRADAS.xlsx
#   - data/raw/sample-metadata.xlsx
#
# Output:
#   - data/processed/ajustado.xlsx
#   - data/processed/entera_5xfad.xlsx
# =============================================================================

# --- Libraries ----------------------------------------------------------------
library(readxl)
library(openxlsx)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(factoextra)
library(FactoMineR)
library(missMDA)
library(VIM)

options(scipen = 999)

# --- Load data ----------------------------------------------------------------
entera_5xfad   <- read_excel("data/raw/230223-LFQ_Incl-Isoforms_FILTRADAS.xlsx")
sample_metadata <- read_xlsx("data/raw/sample-metadata.xlsx")

# --- Fix missing Gene Symbols -------------------------------------------------
gene_symbols <- c(
  "P47964" = "Rpl36",
  "P97461" = "Rps5",
  "Q8C3W1" = "P_Uncharacterized_protein_C1orf198_homolog",
  "Q91V76" = "P_Ester_hydrolase_C11orf54_homolog",
  "Q9CPR4" = "Rpl17",
  "P03987" = "Ig_gamma_3_chain_C_region",
  "Q64331" = "Myo6",
  "Q64436" = "Atp4a",
  "P22315" = "Fech",
  "Q6P4S6" = "Sik3",
  "Q80TK0" = "Btbd8",
  "Q80TT2" = "Baiap3",
  "O88456" = "Capns1",
  "O09167" = "Rpl21",
  "Q8BW41" = "Pomgnt2",
  "P97393" = "Arhgap5",
  "Q03157" = "Aplp1",
  "Q6DFV3" = "Arhgap21",
  "P97798" = "Neo1",
  "Q61285" = "Abcd2",
  "Q8C3W1" = "Uncharacterized_protein_C1orf198_homolog"
)

entera_5xfad <- entera_5xfad %>%
  mutate(
    `Gene Symbol` = ifelse(
      is.na(`Gene Symbol`) & Accession %in% names(gene_symbols),
      gene_symbols[Accession],
      `Gene Symbol`
    )
  )

# Clean species tag from descriptions
entera_5xfad$Description <- gsub("\\[OS=Mus musculus\\]", "", entera_5xfad$Description)
entera_5xfad <- entera_5xfad %>% arrange(Accession)

# --- Extract abundance columns ------------------------------------------------
df_experimento_5xfad <- entera_5xfad

accesion_col    <- dplyr::select(df_experimento_5xfad, Accession)
gene_symbol_col <- dplyr::select(df_experimento_5xfad, `Gene Symbol`)
indice_proteina <- cbind(accesion_col, gene_symbol_col)

colnam  <- colnames(entera_5xfad)
indices <- grep("Abundance: ", colnam)
abundancias_experimento_5xfad <- df_experimento_5xfad[, indices]

# Standardize column names to "mouseXX" format
colnam <- gsub(".*mouse(\\d+).*", "\\1", colnames(abundancias_experimento_5xfad))
colnam <- paste("mouse", gsub("\\D+", "", colnam), sep = "")
colnames(abundancias_experimento_5xfad) <- colnam

id <- accesion_col
df_experimento_5xfad <- cbind(accesion_col, gene_symbol_col, abundancias_experimento_5xfad)
rownames(df_experimento_5xfad) <- id$Accession

df_experimento_5xfad_rowname <- df_experimento_5xfad[, -(1:2)]

# --- QC: Missing values per sample --------------------------------------------
missing_by_sample <- colSums(is.na(df_experimento_5xfad_rowname))
df_missing_sample <- data.frame(
  Muestra = names(missing_by_sample),
  NAs     = missing_by_sample
)

ggplot(df_missing_sample, aes(x = Muestra, y = NAs)) +
  geom_bar(stat = "identity", fill = "red", alpha = 0.7) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Missing values per sample", x = "Sample", y = "Number of NAs")

# --- QC: Preliminary PCA (complete cases) -------------------------------------
tabla_sin_na   <- df_experimento_5xfad_rowname[complete.cases(df_experimento_5xfad_rowname), ]
pca_preliminar <- prcomp(tabla_sin_na, scale. = TRUE, center = TRUE)

pca_df_pre <- as.data.frame(pca_preliminar$x)
ggplot(pca_df_pre, aes(x = PC1, y = PC2)) +
  geom_point() +
  ggtitle("Preliminary PCA (no imputation, complete cases)")

# PCA by sample with metadata
pca_s    <- prcomp(t(tabla_sin_na), scale. = TRUE)
pca_df_s <- as.data.frame(pca_s$x)
pca_df_s$ID        <- sample_metadata$`Sample-id`
pca_df_s$Treatment <- sample_metadata$Genotype
pca_df_s$Sex       <- sample_metadata$sex

ggplot(pca_df_s, aes(x = PC1, y = PC2, color = Treatment, shape = Sex, label = ID)) +
  geom_point(size = 5) +
  geom_text_repel(size = 4, max.overlaps = 15) +
  scale_shape_manual(values = c("male" = 3, "female" = 1)) +
  theme_minimal() +
  labs(
    title = "PCA of protein expression",
    x = paste0("PC1 (", round(summary(pca_s)$importance[2, 1] * 100, 2), "%)"),
    y = paste0("PC2 (", round(summary(pca_s)$importance[2, 2] * 100, 2), "%)")
  )

# --- Separate samples by group ------------------------------------------------
prot_control_hembras <- sample_metadata %>%
  filter(Genotype == "WT",    sex == "female") %>% pull(`Sample-id`)
prot_control_machos  <- sample_metadata %>%
  filter(Genotype == "WT",    sex == "male")   %>% pull(`Sample-id`)
prot_5xfad_machos    <- sample_metadata %>%
  filter(Genotype == "5xFAD", sex == "male")   %>% pull(`Sample-id`)
prot_5xfad_hembras   <- sample_metadata %>%
  filter(Genotype == "5xFAD", sex == "female") %>% pull(`Sample-id`)

tabla_control_hembras <- df_experimento_5xfad_rowname %>% select(all_of(prot_control_hembras))
tabla_control_machos  <- df_experimento_5xfad_rowname %>% select(all_of(prot_control_machos))
tabla_5xfad_machos    <- df_experimento_5xfad_rowname %>% select(all_of(prot_5xfad_machos))
tabla_5xfad_hembras   <- df_experimento_5xfad_rowname %>% select(all_of(prot_5xfad_hembras))

# --- Outlier removal (3×IQR per protein per group) ----------------------------
handle_outliers <- function(x, factor = 3) {
  Q1 <- quantile(x, 0.25, na.rm = TRUE)
  Q3 <- quantile(x, 0.75, na.rm = TRUE)
  IQR_val <- Q3 - Q1
  x[x < (Q1 - factor * IQR_val) | x > (Q3 + factor * IQR_val)] <- NA
  return(x)
}

tabla_control_hembras <- as.data.frame(t(apply(tabla_control_hembras, 1, handle_outliers)))
tabla_control_machos  <- as.data.frame(t(apply(tabla_control_machos,  1, handle_outliers)))
tabla_5xfad_machos    <- as.data.frame(t(apply(tabla_5xfad_machos,    1, handle_outliers)))
tabla_5xfad_hembras   <- as.data.frame(t(apply(tabla_5xfad_hembras,   1, handle_outliers)))

# --- Coverage filter ----------------------------------------------------------
# Note: asymmetric thresholds reflect unequal group sizes in the 5xFAD cohort
tabla_control_hembras <- tabla_control_hembras[rowSums(!is.na(tabla_control_hembras)) >= 3, ]
tabla_control_machos  <- tabla_control_machos[ rowSums(!is.na(tabla_control_machos))  >= 3, ]
tabla_5xfad_machos    <- tabla_5xfad_machos[   rowSums(!is.na(tabla_5xfad_machos))    >= 4, ]
tabla_5xfad_hembras   <- tabla_5xfad_hembras[  rowSums(!is.na(tabla_5xfad_hembras))   >= 2, ]

# --- Keep only proteins detected in all four groups ---------------------------
proteinas_comunes <- Reduce(
  intersect,
  list(
    rownames(tabla_control_hembras),
    rownames(tabla_control_machos),
    rownames(tabla_5xfad_machos),
    rownames(tabla_5xfad_hembras)
  )
)

tabla_control_hembras <- tabla_control_hembras[proteinas_comunes, ]
tabla_control_machos  <- tabla_control_machos[ proteinas_comunes, ]
tabla_5xfad_machos    <- tabla_5xfad_machos[   proteinas_comunes, ]
tabla_5xfad_hembras   <- tabla_5xfad_hembras[  proteinas_comunes, ]

# --- Imputation (within-group median) -----------------------------------------
impute_median <- function(x) {
  x[is.na(x)] <- median(x, na.rm = TRUE)
  return(x)
}

tabla_control_hembras <- as.data.frame(t(apply(tabla_control_hembras, 1, impute_median)))
tabla_control_machos  <- as.data.frame(t(apply(tabla_control_machos,  1, impute_median)))
tabla_5xfad_machos    <- as.data.frame(t(apply(tabla_5xfad_machos,    1, impute_median)))

colnames(tabla_control_hembras) <- prot_control_hembras
colnames(tabla_control_machos)  <- prot_control_machos
colnames(tabla_5xfad_machos)    <- prot_5xfad_machos

imputados <- cbind(tabla_control_hembras, tabla_control_machos,
                   tabla_5xfad_machos, tabla_5xfad_hembras)

# --- Normalization by correction factor ---------------------------------------
media                <- colSums(imputados, na.rm = TRUE) / nrow(imputados)
valor_medio_minimo   <- min(media)
factor_correlaccion  <- media / valor_medio_minimo
factor_correlaccion_rep <- as.data.frame(
  do.call(cbind, lapply(factor_correlaccion, rep, times = nrow(imputados)))
)

ajustado <- imputados / factor_correlaccion_rep
rownames(ajustado) <- rownames(imputados)

# --- Save outputs -------------------------------------------------------------
dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)
write.xlsx(ajustado, "data/processed/ajustado.xlsx", rowNames = TRUE)

entera_5xfad_filtrada <- entera_5xfad %>% filter(Accession %in% rownames(ajustado))
write.xlsx(entera_5xfad_filtrada, "data/processed/entera_5xfad.xlsx", rowNames = TRUE)

message("Preprocessing complete. Files saved to data/processed/")
