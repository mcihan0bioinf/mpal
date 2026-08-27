## 05_granja_signature_violin.R
## Score the LT-PLC "up_all" human signature (all significantly upregulated
## genes, padj<0.05, maximal ortholog mapping) per cell in the Granja et al.
## 2019 atlas via UCell, then show the score distribution per cell type,
## split into the same two annotation layers as 04_atlas_celltype_plots.R:
##   A. MPAL disease classification (author ProjectClassification)
##   B. Healthy reference fine cell types (author BioClassification)
##
## Input:  data/granja_seurat.rds        (celltype, celltype_fine columns)
##         data/human_signature_sets.rds (maximal mapping; up_all set)
## Output: results/granja_signature_violin.png / .pdf

suppressPackageStartupMessages({
  library(Seurat); library(UCell); library(ggplot2); library(dplyr); library(patchwork)
})

proj_root <- getwd()
dir_data  <- file.path(proj_root, "data")
dir_out   <- file.path(proj_root, "results")

seu <- readRDS(file.path(dir_data, "granja_seurat.rds"))
human_sets <- readRDS(file.path(dir_data, "human_signature_sets.rds"))
sig <- list(LT_PLC_up = human_sets$up_all)
cat(sprintf("Signature 'up_all': %d human genes; %d present in atlas\n",
            length(sig$LT_PLC_up), sum(sig$LT_PLC_up %in% rownames(seu))))

## ---- score in memory-bounded batches (UCell on full 53k x 20k can spike RAM) ----
DefaultAssay(seu) <- "RNA"
batch_size <- 10000
cells <- colnames(seu)
batches <- split(cells, ceiling(seq_along(cells) / batch_size))
scores_list <- vector("list", length(batches))
for (i in seq_along(batches)) {
  sub <- subset(seu, cells = batches[[i]])
  sc <- ScoreSignatures_UCell(GetAssayData(sub, layer = "counts"), features = sig,
                               maxRank = 3000)
  scores_list[[i]] <- sc
  rm(sub); gc(verbose = FALSE)
  cat(sprintf("  batch %d/%d scored\n", i, length(batches)))
}
scores <- do.call(rbind, scores_list)[cells, , drop = FALSE]
seu$LT_PLC_up_UCell <- scores[, "LT_PLC_up_UCell"]

mpal_labels <- c("Progenitor_Like","Lymphoid_Like","Myeloid_Like","Healthy_Like","Erythroid_Like","TNK_Like")

df <- data.frame(score = seu$LT_PLC_up_UCell,
                  mpal_class = as.character(seu$celltype),
                  fine_class = as.character(seu$celltype_fine)) %>%
  filter(!is.na(mpal_class))

## ---- Panel A: MPAL disease classification ----
dfA <- df %>% filter(mpal_class %in% mpal_labels) %>%
  mutate(mpal_class = reorder(mpal_class, score, median))
ordA <- levels(dfA$mpal_class)

pA <- ggplot(dfA, aes(x = mpal_class, y = score, fill = mpal_class)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.85) +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", alpha = 0.6) +
  scale_fill_manual(values = scales::hue_pal()(length(ordA)), guide = "none") +
  labs(title = "A. MPAL disease classification",
       x = NULL, y = "UCell score (LT-PLC up signature)") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1),
        plot.title = element_text(face = "bold", size = 12),
        panel.grid.minor = element_blank())

## ---- Panel B: healthy reference fine cell types ----
dfB <- df %>% filter(!mpal_class %in% mpal_labels) %>%
  mutate(fine_class = reorder(fine_class, score, median))
ordB <- levels(dfB$fine_class)

pB <- ggplot(dfB, aes(x = fine_class, y = score, fill = fine_class)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.85) +
  geom_boxplot(width = 0.25, outlier.shape = NA, fill = "white", alpha = 0.6) +
  scale_fill_manual(values = colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(length(ordB)), guide = "none") +
  labs(title = "B. Healthy reference fine cell types",
       x = NULL, y = "UCell score (LT-PLC up signature)") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 8),
        plot.title = element_text(face = "bold", size = 12),
        panel.grid.minor = element_blank())

p <- (pA | pB) +
  plot_annotation(title = "Granja et al. 2019: LT-PLC upregulated signature score by cell type",
                   subtitle = sprintf("UCell, %d genes (%d matched in atlas); ordered by median score within each panel",
                                       length(sig$LT_PLC_up), sum(sig$LT_PLC_up %in% rownames(seu))),
                   theme = theme(plot.title = element_text(face = "bold", size = 13),
                                 plot.subtitle = element_text(size = 9.5, color = "grey30")))

ggsave(file.path(dir_out, "granja_signature_violin.png"), p, width = 15, height = 7, dpi = 300)
ggsave(file.path(dir_out, "granja_signature_violin.pdf"), p, width = 15, height = 7)
cat("Saved results/granja_signature_violin.png / .pdf\n")

## print median score per cell type for a quick numeric read
cat("\nMedian UCell score, MPAL classification layer:\n")
print(dfA %>% group_by(mpal_class) %>% summarise(median_score = median(score), n = n()) %>% arrange(desc(median_score)))
cat("\nMedian UCell score, healthy fine cell types:\n")
print(dfB %>% group_by(fine_class) %>% summarise(median_score = median(score), n = n()) %>% arrange(desc(median_score)))

saveRDS(seu, file.path(dir_data, "granja_seurat.rds"))
