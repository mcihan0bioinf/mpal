## 06_granja_signature_violin_topN.R
## Score LT-PLC top50 / top100 / top500 upregulated human signatures (maximal
## ortholog mapping) per cell in Granja et al. 2019, then show the score
## distribution per cell type across the same two annotation layers as
## 04_atlas_celltype_plots.R / 05_granja_signature_violin.R:
##   A. MPAL disease classification (author ProjectClassification)
##   B. Healthy reference fine cell types (author BioClassification)
##
## Cell-type order on the x-axis is fixed to the up_all-signature median
## order from 05_granja_signature_violin.R, so all panels are comparable.
##
## Input:  data/granja_seurat.rds
##         data/human_signature_sets.rds (maximal mapping; top50/top100/top500)
## Output: results/granja_signature_violin_topN_layerA.png / .pdf
##         results/granja_signature_violin_topN_layerB.png / .pdf

set.seed(1234)
suppressPackageStartupMessages({
  library(Seurat); library(UCell); library(ggplot2); library(dplyr); library(patchwork)
})

proj_root <- getwd()
dir_data  <- file.path(proj_root, "data")
dir_out   <- file.path(proj_root, "results")

seu <- readRDS(file.path(dir_data, "granja_seurat.rds"))
human_sets <- readRDS(file.path(dir_data, "human_signature_sets.rds"))
DefaultAssay(seu) <- "RNA"
all_genes <- rownames(seu)

variants <- c("up_top50", "up_top100", "up_top500")
sig <- lapply(variants, function(v) intersect(human_sets[[v]], all_genes))
names(sig) <- variants
for (nm in names(sig)) cat(sprintf("  %-12s n = %d\n", nm, length(sig[[nm]])))

## ---- score in memory-bounded batches ----
batch_size <- 10000
cells <- colnames(seu)
batches <- split(cells, ceiling(seq_along(cells) / batch_size))
scores_list <- vector("list", length(batches))
for (i in seq_along(batches)) {
  sub <- subset(seu, cells = batches[[i]])
  sc <- ScoreSignatures_UCell(GetAssayData(sub, layer = "counts"), features = sig, maxRank = 3000)
  scores_list[[i]] <- sc
  rm(sub); gc(verbose = FALSE)
  cat(sprintf("  batch %d/%d scored\n", i, length(batches)))
}
scores <- do.call(rbind, scores_list)[cells, , drop = FALSE]
colnames(scores) <- sub("_UCell$", "", colnames(scores))

mpal_labels <- c("Progenitor_Like","Lymphoid_Like","Myeloid_Like","Healthy_Like","Erythroid_Like","TNK_Like")

meta <- data.frame(cell = cells,
                    mpal_class = as.character(seu$celltype),
                    fine_class = as.character(seu$celltype_fine))

df <- cbind(meta, as.data.frame(scores)) %>% filter(!is.na(mpal_class))

size_labels <- c(up_top50 = "Top 50 genes", up_top100 = "Top 100 genes", up_top500 = "Top 500 genes")
df_long <- df %>%
  tidyr::pivot_longer(cols = all_of(colnames(scores)), names_to = "signature", values_to = "score") %>%
  mutate(size = factor(size_labels[signature], levels = size_labels))

## fixed cell-type order = up_all median order already established
orderA <- c("Healthy_Like","TNK_Like","Lymphoid_Like","Erythroid_Like","Myeloid_Like","Progenitor_Like")

dfA <- df_long %>% filter(mpal_class %in% mpal_labels) %>%
  mutate(mpal_class = factor(mpal_class, levels = orderA))

pA <- ggplot(dfA, aes(x = mpal_class, y = score, fill = mpal_class)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.85) +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", alpha = 0.6) +
  facet_wrap(~size, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = scales::hue_pal()(length(orderA)), guide = "none") +
  labs(title = "A. MPAL disease classification",
       x = NULL, y = "UCell score") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1),
        plot.title = element_text(face = "bold", size = 12),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"))

ggsave(file.path(dir_out, "granja_signature_violin_topN_layerA.png"), pA, width = 12, height = 6, dpi = 300)
ggsave(file.path(dir_out, "granja_signature_violin_topN_layerA.pdf"), pA, width = 12, height = 6)

## for the fine healthy layer, order by the top500 median (closest proxy to up_all)
orderB <- df_long %>% filter(!mpal_class %in% mpal_labels, size == "Top 500 genes") %>%
  group_by(fine_class) %>% summarise(m = median(score), .groups = "drop") %>% arrange(m) %>% pull(fine_class)

dfB <- df_long %>% filter(!mpal_class %in% mpal_labels) %>%
  mutate(fine_class = factor(fine_class, levels = orderB))

pB <- ggplot(dfB, aes(x = fine_class, y = score, fill = fine_class)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.85) +
  geom_boxplot(width = 0.25, outlier.shape = NA, fill = "white", alpha = 0.6) +
  facet_wrap(~size, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(length(orderB)), guide = "none") +
  labs(title = "B. Healthy reference fine cell types",
       x = NULL, y = "UCell score") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 7),
        plot.title = element_text(face = "bold", size = 12),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"))

ggsave(file.path(dir_out, "granja_signature_violin_topN_layerB.png"), pB, width = 15, height = 6.5, dpi = 300)
ggsave(file.path(dir_out, "granja_signature_violin_topN_layerB.pdf"), pB, width = 15, height = 6.5)

cat("\nSaved results/granja_signature_violin_topN_layerA.png/.pdf and layerB.png/.pdf\n")

cat("\nMedian UCell score, layer A:\n")
print(dfA %>% group_by(size, mpal_class) %>% summarise(median_score = median(score), .groups = "drop") %>%
        arrange(size, desc(median_score)))
cat("\nMedian UCell score, layer B:\n")
print(dfB %>% group_by(size, fine_class) %>% summarise(median_score = median(score), .groups = "drop") %>%
        arrange(size, desc(median_score)))
