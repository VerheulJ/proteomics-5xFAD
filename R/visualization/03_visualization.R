# =============================================================================
# 03_visualization.R
# Proteomics Analysis — 5xFAD Mouse Model
#
# Generates: PCA, Volcano plot, Heatmap, KEGG enrichment,
#            GO enrichment, Protein-protein interaction network
#
# Input:
#   - data/processed/ajustado.xlsx
#   - data/processed/entera_5xfad.xlsx
#   - data/raw/sample-metadata.xlsx
#   - results/tables/vulcanot_5xfad.xlsx
#   - results/tables/df_desreguladas_5xfad.xlsx
#   - data/processed/Enrichment_Results_5xfad_005.xlsx
#   - data/processed/RED_5xfad_005.xlsx
#
# Output: results/figures/
# =============================================================================

# --- Libraries ----------------------------------------------------------------
library(readxl)
library(openxlsx)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(igraph)
library(ggraph)
library(factoextra)

options(scipen = 999)

dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)

# Helper: save in TIFF, SVG, PDF and PNG
save_plot <- function(plot_obj, basename, width = 7.5, height = 7.5) {
  base <- file.path("results/figures", basename)
  tiff(paste0(base, ".tif"), width = width, height = height, units = "in", res = 600)
  print(plot_obj); dev.off()
  svg( paste0(base, ".svg"), width = width, height = height)
  print(plot_obj); dev.off()
  pdf( paste0(base, ".pdf"), width = width, height = height)
  print(plot_obj); dev.off()
  ggsave(paste0(base, ".png"), plot = plot_obj, width = width, height = height)
  message("Saved: ", basename)
}

# --- Load data ----------------------------------------------------------------
entera_experimento_5xfad <- read_excel("data/processed/entera_5xfad.xlsx")
sample_metadata           <- read_xlsx("data/raw/sample-metadata.xlsx")
ajustado_raw              <- read_xlsx("data/processed/ajustado.xlsx")
vulcanot_5xfad            <- read_xlsx("results/tables/vulcanot_5xfad.xlsx")
df_desreguladas_5xfad     <- read_xlsx("results/tables/df_desreguladas_5xfad.xlsx")

ajustado <- as.data.frame(ajustado_raw)
rownames(ajustado) <- ajustado_raw$...1
ajustado$...1 <- NULL

accesion_col    <- dplyr::select(entera_experimento_5xfad, Accession)
gene_symbol_col <- dplyr::select(entera_experimento_5xfad, `Gene Symbol`)
indice_proteina <- cbind(accesion_col, gene_symbol_col)

sample_metadata <- as.data.frame(sample_metadata)
rownames(sample_metadata) <- sample_metadata[, 1]

prot_control <- sample_metadata %>% filter(Genotype == "WT")    %>% pull(`Sample-id`)
prot_5xfad   <- sample_metadata %>% filter(Genotype == "5xFAD") %>% pull(`Sample-id`)

# ==============================================================================
# 1. PCA
# ==============================================================================
pca    <- prcomp(t(ajustado), scale. = TRUE)
pca_df <- as.data.frame(pca$x)
pca_df$ID        <- sample_metadata$`Sample-id`
pca_df$Treatment <- sample_metadata$Genotype
pca_df$Sex       <- sample_metadata$sex

pca_var <- summary(pca)$importance

g_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Treatment, shape = Sex, label = ID)) +
  geom_point(size = 5) +
  geom_text_repel(size = 4, max.overlaps = 15) +
  scale_shape_manual(values = c("male" = 3, "female" = 1)) +
  theme_minimal() +
  labs(
    title = "PCA of protein expression — 5xFAD vs WT",
    x = paste0("PC1 (", round(pca_var[2, 1] * 100, 2), "%)"),
    y = paste0("PC2 (", round(pca_var[2, 2] * 100, 2), "%)")
  )

save_plot(g_pca, "PCA_5xfad")

# ==============================================================================
# 2. Volcano plot
# ==============================================================================
# Stricter thresholds than 3xTg: p < 0.01, FC > 50%
p_umbral <- -log10(0.01)   # = 2
FC        <- 0.5
LFCP      <- log2(1 + FC)
LFCN      <- log2(1 - FC)

up_regulated   <- filter(vulcanot_5xfad, log2cociente_5xfad >  LFCP & log10pvalor_5xfad > p_umbral)
down_regulated <- filter(vulcanot_5xfad, log2cociente_5xfad <  LFCN & log10pvalor_5xfad > p_umbral)
resto          <- anti_join(vulcanot_5xfad, up_regulated) %>% anti_join(down_regulated)

up_regulated$Status   <- "up-regulated"
down_regulated$Status <- "down-regulated"

# Label selection
vulcanot_sig <- filter(vulcanot_5xfad, log10pvalor_5xfad > p_umbral)
top_sig      <- vulcanot_sig %>% arrange(desc(log10pvalor_5xfad)) %>% head(3)
top_fc_up    <- filter(vulcanot_sig, log2cociente_5xfad > LFCP) %>% arrange(desc(log2cociente_5xfad)) %>% head(3)
top_fc_down  <- filter(vulcanot_sig, log2cociente_5xfad < LFCN) %>% arrange(log2cociente_5xfad) %>% head(3)
subset_etiquetas <- bind_rows(top_sig, top_fc_up, top_fc_down) %>% distinct(Accession, .keep_all = TRUE)

a    <- length(prot_control) + length(prot_5xfad)
ref  <- 4.5; punt <- -10

g_volcano <- ggplot() +
  geom_point(data = up_regulated,   aes(x = log2cociente_5xfad, y = log10pvalor_5xfad, color = Status), size = 0.75) +
  geom_point(data = down_regulated, aes(x = log2cociente_5xfad, y = log10pvalor_5xfad, color = Status), size = 0.75) +
  geom_point(data = resto,          aes(x = log2cociente_5xfad, y = log10pvalor_5xfad), color = "grey", size = 0.5) +
  geom_vline(xintercept = c(LFCP, LFCN), linetype = "dashed", color = "black") +
  geom_hline(yintercept = p_umbral, linetype = "dashed", color = "black") +
  geom_label_repel(
    data = subset_etiquetas,
    aes(x = log2cociente_5xfad, y = log10pvalor_5xfad, label = Gene.Symbol),
    color = "black", box.padding = 0.2, point.padding = 0.1, force = 10
  ) +
  annotate("text", x = punt, y = ref + 2.095, label = "p-value(0.01) = 2",                        size = 4, hjust = 0) +
  annotate("text", x = punt, y = ref + 1.945, label = paste("N (proteins) =", nrow(vulcanot_5xfad)), size = 4, hjust = 0) +
  annotate("text", x = punt, y = ref + 1.795, label = paste("N (mice) =",     a),                  size = 4, hjust = 0) +
  annotate("text", x = punt, y = ref + 1.645, label = paste("down-regulated =", nrow(down_regulated)), size = 4, hjust = 0) +
  annotate("text", x = punt, y = ref + 1.495, label = paste("up-regulated =",   nrow(up_regulated)),   size = 4, hjust = 0) +
  scale_color_manual(values = c("up-regulated" = "red", "down-regulated" = "blue")) +
  coord_cartesian(xlim = c(-10, 10), ylim = c(0, 10)) +
  labs(
    x     = bquote(log[2](Fold~Change)),
    y     = bquote(-log[10](p~value)),
    color = ""
  ) +
  theme(
    panel.background = element_rect(fill = "white"),
    axis.line        = element_line(color = "black"),
    legend.text      = element_text(size = 12)
  )

save_plot(g_volcano, "efecto_5xfad")

# ==============================================================================
# 3. Heatmap
# ==============================================================================
desreguladas    <- df_desreguladas_5xfad$Accession
todas_muestras  <- c(prot_5xfad, prot_control)

matriz_elegidas <- ajustado[desreguladas, todas_muestras]
log_matriz      <- t(apply(matriz_elegidas, 1, log10))

desreguladas_nombres <- indice_proteina[match(rownames(log_matriz), indice_proteina$Accession), ]
desreguladas_nombres$Gene_Unique <- ifelse(
  grepl("-", desreguladas_nombres$Accession),
  paste0(desreguladas_nombres$`Gene Symbol`, "_iso",
         sub(".*-", "", desreguladas_nombres$Accession)),
  desreguladas_nombres$`Gene Symbol`
)
rownames(log_matriz) <- desreguladas_nombres$Gene_Unique

custom_colors <- colorRampPalette(c("blue", "white", "red"))(100)
breaks        <- seq(range(log_matriz, na.rm = TRUE)[1],
                     range(log_matriz, na.rm = TRUE)[2],
                     length.out = 101)

annotation_data <- sample_metadata[, "Genotype", drop = FALSE]

p_heatmap <- pheatmap(
  log_matriz,
  annotation_col       = annotation_data,
  color                = custom_colors,
  breaks               = breaks,
  border_color         = NA,
  show_rownames        = TRUE,
  show_colnames        = TRUE,
  fontsize_row         = 6,
  scale                = "none",
  angle_col            = 90,
  fontsize_col         = 10,
  clustering_distance_rows = "euclidean",
  clustering_method    = "complete",
  cluster_cols         = TRUE,
  silent               = TRUE
)

base_hm <- file.path("results/figures", "heatmap_resultado_5xfad")
tiff(paste0(base_hm, ".tif"), width = 11, height = 16, units = "in", res = 600); print(p_heatmap); dev.off()
svg( paste0(base_hm, ".svg"), width = 11, height = 16);                          print(p_heatmap); dev.off()
pdf( paste0(base_hm, ".pdf"), width = 11, height = 16);                          print(p_heatmap); dev.off()
png( paste0(base_hm, ".png"), width = 11, height = 16, units = "in", res = 150); print(p_heatmap); dev.off()
message("Saved: heatmap_resultado_5xfad")

# ==============================================================================
# 4. KEGG Enrichment
# ==============================================================================
Enrichment_Results <- read_excel("data/processed/Enrichment_Results_5xfad_005.xlsx")

termKEGG <- Enrichment_Results %>%
  filter(category == "KEGG", description != "Metabolic pathways") %>%
  mutate(count = sapply(strsplit(as.character(preferred_names), ","), length))

top_kegg <- termKEGG %>% group_by(category) %>% top_n(20, -fdr) %>% ungroup()

g_kegg <- ggplot(top_kegg, aes(x = reorder(description, -fdr), y = count, fill = fdr)) +
  geom_bar(stat = "identity", width = 0.6) +
  coord_flip() +
  scale_fill_gradient(low = "blue", high = "red") +
  facet_wrap(~ category, scales = "free_y", ncol = 1, strip.position = "top") +
  labs(title = "KEGG terms — 5xFAD", x = "Pathway", y = "Count", fill = "FDR") +
  theme_minimal() +
  theme(
    axis.text.y   = element_text(size = 10),
    strip.text    = element_text(size = 14, face = "bold"),
    panel.spacing = unit(3, "lines"),
    plot.title    = element_text(size = 16, face = "bold", hjust = 0.5)
  )

save_plot(g_kegg, "enrichkegg_5xfad", width = 10, height = 10)

# ==============================================================================
# 5. GO Enrichment
# ==============================================================================
termGO <- Enrichment_Results %>%
  filter(category %in% c("Process", "Component", "Function")) %>%
  mutate(
    num_genes = str_count(preferred_names, ",") + 1,
    count     = sapply(strsplit(as.character(preferred_names), ","), length)
  )

top_go <- termGO %>%
  group_by(category) %>% top_n(20, num_genes) %>%
  arrange(category, fdr) %>% ungroup()

g_go <- ggplot(top_go, aes(x = reorder(description, -fdr), y = count, fill = fdr)) +
  geom_bar(stat = "identity", width = 0.6) +
  coord_flip() +
  scale_fill_gradient(low = "blue", high = "red") +
  facet_wrap(~ category, scales = "free_y", ncol = 1, strip.position = "top") +
  labs(title = "GO terms — 5xFAD", x = "GO term", y = "Count", fill = "FDR") +
  theme_minimal() +
  theme(
    axis.text.y   = element_text(size = 7),
    strip.text    = element_text(size = 14, face = "bold"),
    panel.spacing = unit(5, "lines"),
    plot.title    = element_text(size = 16, face = "bold", hjust = 0.5)
  )

save_plot(g_go, "enrichGO_5xfad", width = 10, height = 10)

# ==============================================================================
# 6. Protein–Protein Interaction Network
# ==============================================================================
brEtOH          <- read_excel("data/processed/RED_5xfad_005.xlsx")
brinteracciones <- brEtOH[, 3:4]

replacements <- c(
  "Rps10l1" = "Rps10",       "LOC102555453" = "Rpl12",    "Qk" = "Qki",
  "LRRTM1"  = "Lrrtm1",     "LOC100361558" = "H3-3b",    "LOC100360449" = "Rpl9",
  "LOC100361756" = "Rps26",  "LOC103690091" = "Pcsk1n",   "LOC100363469" = "Rps24",
  "Rpl36al" = "Rpl36a",     "Slc9a3r2" = "Nherf2",      "Serpina3c" = "Serpina3k",
  "Car1" = "Ca1",            "Hmgb1-2" = "Hmgb1",        "COX2" = "Mtco2",
  "Srp54a" = "Srp54",       "Eif5b-2" = "Eif5b",        "H1f5" = "H1-5",
  "Cox6c" = "Cox6c2",       "Calm3" = "Calm2",           "LOC103693015" = "Vkorc1l1",
  "Tuba3b" = "Tuba3a; Tuba3b", "ATP8" = "Mt-atp8",      "Rack1-2" = "Rack1",
  "LOC100912380" = "Capns1", "LOC679739" = "Ndufs6",     "LOC100911575" = "Rplp2",
  "Hist1h2ah" = "P_Histone_H2a_Type1e", "H2ac1" = "P_Histone_H2a_Type4",
  "Hist2h2aa2" = "H2ac18",  "Hbb-b1" = "P_Hemoglobin_Subunit_Beta-2",
  "Atad3a" = "Atad3",       "Slc9a3r1" = "Nherf1",      "Timm8a1" = "Timm8a",
  "ND5" = "Mtnd5",          "Atp5md" = "Atp5mk",        "LOC100911110" = "Eif3h",
  "Rps21-2" = "Rps21",      "LOC100911402" = "Cend1",    "LOC103693564" = "Tagln3",
  "LOC100911769" = "Epb41l1","Rpl7a-2" = "Rpl7a",       "Acy1" = "Acy1a",
  "Hist1h2bc" = "P_Histone_H2b_Type1", "Hist2h4" = "H4c2-Hist1h4m-H4c16",
  "Rps18l1" = "Rps18",      "LOC100359951" = "Rps20",   "Lnpep-2" = "Lnpep",
  "LOC100360117" = "Rpl8",  "Car2" = "Ca2",             "F1M7L9_RAT" = "Pkp2",
  "Phb" = "Phb1",           "Ppial4d" = "Ppia",         "LOC100911372" = "Rps6",
  "LOC103689983" = "Hprt1", "LOC100909548" = "M6pr",    "Txn1" = "Txn",
  "Gm10053" = "Cycs",       "Atp5d" = "Atp5f1d",        "Atp5h" = "Atp5pd",
  "Ddx39" = "Ddx39a",       "H2aw" = "H2ac25",          "Mpp6" = "Pals2"
)

for (old_name in names(replacements)) {
  brinteracciones <- replace(brinteracciones, brinteracciones == old_name, replacements[[old_name]])
}

proteinas_con_interacciones <- unique(c(brinteracciones[[1]], brinteracciones[[2]]))

df_net <- df_desreguladas_5xfad %>%
  filter(Gene.Symbol %in% proteinas_con_interacciones) %>%
  group_by(Gene.Symbol) %>%
  filter(!(any(!grepl("-[0-9]+$", Accession)) & grepl("-[0-9]+$", Accession))) %>%
  ungroup() %>%
  distinct()

proteinas_all <- tibble(
  Gene.Symbol = df_net$Gene.Symbol,
  P_Value     = df_net$log10pvalor_5xfad,
  Cociente    = df_net$log2cociente_5xfad
)

g_net <- graph_from_data_frame(d = brinteracciones, vertices = proteinas_all, directed = FALSE)
V(g_net)$cociente  <- proteinas_all$Cociente
V(g_net)$p_valor   <- proteinas_all$P_Value
V(g_net)$categoria <- ifelse(V(g_net)$cociente > 0, "Up-regulated", "Down-regulated")

mi_paleta <- colorRampPalette(c("blue", "white", "red"))(100)

g_network <- ggraph(g_net, layout = "fr", niter = 10000) +
  geom_edge_link() +
  geom_node_point(aes(size = p_valor, fill = cociente, shape = categoria), colour = "black") +
  geom_node_label(aes(label = name), family = "sans", size = 3, fill = "white", color = "black") +
  theme_graph(base_family = "sans") +
  scale_shape_manual(values = c(21, 22)) +
  scale_fill_gradientn(colors = mi_paleta, limits = c(-5, 5)) +
  scale_size(range = c(5, 12)) +
  labs(
    fill  = "Log2 Fold Change",
    shape = "",
    size  = "-log10(p-value)",
    title = "Interaction network of deregulated proteins in 5xFAD"
  ) +
  theme(plot.title = element_text(family = "sans", size = 16, face = "bold"))

save_plot(g_network, "RED_5xfad", width = 10, height = 10)

message("All figures saved to results/figures/")
