## 03_atlas_celltypes_granja.R
## UMAP of the Granja et al. 2019 atlas (MPAL + healthy hematopoiesis),
## split into its two independent annotation layers so they are not
## conflated in one color scale:
##   A. MPAL disease classification (author ProjectClassification) --
##      6 categories, assigned only to the malignant MPAL blasts.
##   B. Healthy reference fine cell types (author BioClassification,
##      pulled from the separate scRNA-Healthy-Hematopoiesis-191120.rds
##      file and merged in by Group+Barcode) -- 26 categories, assigned
##      only to the healthy comparator cells.
## The two layers are complementary: every cell has an assignment in
## exactly one of them (celltype_fine in data/granja_seurat.rds already
## combines both, added during dataset preparation).
##
## Input:  data/granja_seurat.rds (celltype, celltype_fine, umap)
## Output: results/atlas_celltypes_granja.png / .pdf

suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(dplyr); library(patchwork)
})

proj_root <- getwd()
dir_data  <- file.path(proj_root, "data")
dir_out   <- file.path(proj_root, "results")

seu <- readRDS(file.path(dir_data, "granja_seurat.rds"))
emb <- Embeddings(seu, "umap")[, 1:2]

mpal_labels <- c("Progenitor_Like","Lymphoid_Like","Myeloid_Like","Healthy_Like","Erythroid_Like","TNK_Like")

df <- data.frame(UMAP_1 = emb[, 1], UMAP_2 = emb[, 2],
                  mpal_class = as.character(seu$celltype),
                  fine_class = as.character(seu$celltype_fine))
df <- df %>% filter(!is.na(mpal_class))

## ---- Panel A: MPAL disease classification layer ----
dfA <- df %>% mutate(group = ifelse(mpal_class %in% mpal_labels, mpal_class, NA_character_))
labA <- dfA %>% filter(!is.na(group)) %>% group_by(group) %>%
  summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2), .groups = "drop")
palA <- scales::hue_pal()(length(mpal_labels)); names(palA) <- mpal_labels

pA <- ggplot() +
  geom_point(data = dfA %>% filter(is.na(group)), aes(UMAP_1, UMAP_2), color = "grey88", size = 0.25, alpha = 0.4) +
  geom_point(data = dfA %>% filter(!is.na(group)), aes(UMAP_1, UMAP_2, color = group), size = 0.4, alpha = 0.75) +
  ggrepel::geom_text_repel(data = labA, aes(UMAP_1, UMAP_2, label = group), color = "black", size = 3,
                            fontface = "bold", bg.color = "white", bg.r = 0.12, max.overlaps = 20) +
  scale_color_manual(values = palA) + guides(color = "none") +
  labs(title = "A. MPAL disease classification",
       subtitle = "colored = MPAL blasts (author ProjectClassification)\ngrey = healthy reference cells (unclassified in this layer)") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 8, color = "grey30"),
        panel.grid.minor = element_blank())

## ---- Panel B: healthy reference fine classification layer ----
dfB <- df %>% mutate(group = ifelse(!mpal_class %in% mpal_labels, fine_class, NA_character_))
labB <- dfB %>% filter(!is.na(group)) %>% group_by(group) %>%
  summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2), .groups = "drop")
fine_types <- sort(unique(dfB$group[!is.na(dfB$group)]))
palB <- colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(length(fine_types)); names(palB) <- fine_types

pB <- ggplot() +
  geom_point(data = dfB %>% filter(is.na(group)), aes(UMAP_1, UMAP_2), color = "grey88", size = 0.25, alpha = 0.4) +
  geom_point(data = dfB %>% filter(!is.na(group)), aes(UMAP_1, UMAP_2, color = group), size = 0.4, alpha = 0.75) +
  ggrepel::geom_text_repel(data = labB, aes(UMAP_1, UMAP_2, label = group), color = "black", size = 2.4,
                            bg.color = "white", bg.r = 0.1, max.overlaps = 30, segment.size = 0.2) +
  scale_color_manual(values = palB) + guides(color = "none") +
  labs(title = "B. Healthy reference fine cell types",
       subtitle = "colored = healthy comparator cells (BioClassification, 26 types)\ngrey = MPAL blasts (unclassified in this layer)") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 8, color = "grey30"),
        panel.grid.minor = element_blank())

p <- (pA | pB) +
  plot_annotation(title = "Granja et al. 2019 -- MPAL + healthy hematopoiesis: two independent annotation layers",
                   theme = theme(plot.title = element_text(face = "bold", size = 13)))

ggsave(file.path(dir_out, "atlas_celltypes_granja.png"), p, width = 15, height = 7.5, dpi = 300)
ggsave(file.path(dir_out, "atlas_celltypes_granja.pdf"), p, width = 15, height = 7.5)
cat("Saved results/atlas_celltypes_granja.png / .pdf\n")
