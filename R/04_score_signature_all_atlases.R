## 04_score_signature_all_atlases.R
## Score the LT-PLC "up_all" human signature (UCell, maximal ortholog
## mapping) per cell in all four atlases, saving the score back into each
## Seurat object (data/<atlas>_seurat.rds$LT_PLC_up_UCell) so downstream
## scripts can reuse it without rescoring. Also produces one violin figure
## per atlas showing the score distribution by cell type / disease group.
##
## Granja gets a two-panel figure since it has two independent annotation
## layers (see 03_atlas_celltypes_granja.R); the other three atlases get a
## single cell-type panel plus a disease-group panel.
##
## Input:  data/granja_seurat.rds, data/caron_seurat.rds,
##         data/vangalen_seurat.rds, data/bonemarrowmap_seurat.rds
##         data/human_signature_sets.rds (maximal mapping; up_all set)
## Output: results/granja_signature_violin.png / .pdf
##         results/caron_signature_violin_celltype.png / .pdf
##         results/caron_signature_violin_group.png / .pdf
##         results/vangalen_signature_violin_celltype.png / .pdf
##         results/vangalen_signature_violin_group.png / .pdf
## (BoneMarrowMap is scored here too but has no dedicated figure -- it is
##  only used as the healthy baseline in 06_cross_atlas_comparison.R)

suppressPackageStartupMessages({
  library(Seurat); library(UCell); library(ggplot2); library(dplyr); library(patchwork)
})

proj_root <- getwd()
dir_data  <- file.path(proj_root, "data")
dir_out   <- file.path(proj_root, "results")

human_sets <- readRDS(file.path(dir_data, "human_signature_sets.rds"))
sig <- list(LT_PLC_up = human_sets$up_all)

## ---- shared helper: score a Seurat object in memory-bounded batches ----
score_atlas <- function(seu, batch_size = 10000) {
  DefaultAssay(seu) <- "RNA"
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
  seu$LT_PLC_up_UCell <- scores[, "LT_PLC_up_UCell"]
  seu
}

make_violin <- function(seu, celltype_col, title, out_png, out_pdf, angle = 45) {
  df <- data.frame(score = seu$LT_PLC_up_UCell, celltype = as.character(seu[[celltype_col, drop = TRUE]])) %>%
    filter(!is.na(celltype), celltype != "") %>%
    mutate(celltype = reorder(celltype, score, median))
  n_types <- length(unique(df$celltype))
  pal <- if (n_types <= 12) scales::hue_pal()(n_types) else
    colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(n_types)

  p <- ggplot(df, aes(x = celltype, y = score, fill = celltype)) +
    geom_violin(scale = "width", trim = TRUE, alpha = 0.85) +
    geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", alpha = 0.6) +
    scale_fill_manual(values = pal, guide = "none") +
    labs(title = title, x = NULL, y = "UCell score (LT-PLC up_all signature)") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = angle, hjust = 1),
          plot.title = element_text(face = "bold", size = 12.5),
          panel.grid.minor = element_blank())

  ggsave(out_png, p, width = 10, height = 6.5, dpi = 300)
  ggsave(out_pdf, p, width = 10, height = 6.5)
  cat(sprintf("Median UCell score by %s:\n", celltype_col))
  print(df %>% group_by(celltype) %>% summarise(median_score = median(score), n = n()) %>% arrange(desc(median_score)))
  df
}

## ================= Granja 2019 (MPAL) =================
cat("=== Granja 2019 ===\n")
granja <- readRDS(file.path(dir_data, "granja_seurat.rds"))
granja <- score_atlas(granja)
saveRDS(granja, file.path(dir_data, "granja_seurat.rds"))

mpal_labels <- c("Progenitor_Like","Lymphoid_Like","Myeloid_Like","Healthy_Like","Erythroid_Like","TNK_Like")
df <- data.frame(score = granja$LT_PLC_up_UCell,
                  mpal_class = as.character(granja$celltype),
                  fine_class = as.character(granja$celltype_fine)) %>%
  filter(!is.na(mpal_class))

dfA <- df %>% filter(mpal_class %in% mpal_labels) %>% mutate(mpal_class = reorder(mpal_class, score, median))
pA <- ggplot(dfA, aes(x = mpal_class, y = score, fill = mpal_class)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.85) +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", alpha = 0.6) +
  scale_fill_manual(values = scales::hue_pal()(length(levels(dfA$mpal_class))), guide = "none") +
  labs(title = "A. MPAL disease classification", x = NULL, y = "UCell score (LT-PLC up signature)") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1),
        plot.title = element_text(face = "bold", size = 12), panel.grid.minor = element_blank())

dfB <- df %>% filter(!mpal_class %in% mpal_labels) %>% mutate(fine_class = reorder(fine_class, score, median))
pB <- ggplot(dfB, aes(x = fine_class, y = score, fill = fine_class)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.85) +
  geom_boxplot(width = 0.25, outlier.shape = NA, fill = "white", alpha = 0.6) +
  scale_fill_manual(values = colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(length(levels(dfB$fine_class))), guide = "none") +
  labs(title = "B. Healthy reference fine cell types", x = NULL, y = "UCell score (LT-PLC up signature)") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 8),
        plot.title = element_text(face = "bold", size = 12), panel.grid.minor = element_blank())

p_granja <- (pA | pB) +
  plot_annotation(title = "Granja et al. 2019: LT-PLC upregulated signature score by cell type",
                   subtitle = sprintf("UCell, %d genes (%d matched in atlas); ordered by median score within each panel",
                                       length(sig$LT_PLC_up), sum(sig$LT_PLC_up %in% rownames(granja))),
                   theme = theme(plot.title = element_text(face = "bold", size = 13),
                                 plot.subtitle = element_text(size = 9.5, color = "grey30")))
ggsave(file.path(dir_out, "granja_signature_violin.png"), p_granja, width = 15, height = 7, dpi = 300)
ggsave(file.path(dir_out, "granja_signature_violin.pdf"), p_granja, width = 15, height = 7)
cat("Saved results/granja_signature_violin.png / .pdf\n")
rm(granja, df, dfA, dfB, p_granja); gc()

## ================= Caron 2020 (B-ALL / T-ALL) =================
cat("\n=== Caron 2020 ===\n")
caron <- readRDS(file.path(dir_data, "caron_seurat.rds"))
caron <- score_atlas(caron)
saveRDS(caron, file.path(dir_data, "caron_seurat.rds"))
make_violin(caron, "celltype", "Caron et al. 2020: LT-PLC score by cell/sample type",
            file.path(dir_out, "caron_signature_violin_celltype.png"),
            file.path(dir_out, "caron_signature_violin_celltype.pdf"), angle = 40)
make_violin(caron, "group", "Caron et al. 2020: LT-PLC score by disease group",
            file.path(dir_out, "caron_signature_violin_group.png"),
            file.path(dir_out, "caron_signature_violin_group.pdf"), angle = 0)
rm(caron); gc()

## ================= van Galen 2019 (AML) =================
cat("\n=== van Galen 2019 ===\n")
vg <- readRDS(file.path(dir_data, "vangalen_seurat.rds"))
vg <- score_atlas(vg)
saveRDS(vg, file.path(dir_data, "vangalen_seurat.rds"))
make_violin(vg, "celltype", "van Galen et al. 2019: LT-PLC score by cell state",
            file.path(dir_out, "vangalen_signature_violin_celltype.png"),
            file.path(dir_out, "vangalen_signature_violin_celltype.pdf"), angle = 45)
make_violin(vg, "group", "van Galen et al. 2019: LT-PLC score by disease group",
            file.path(dir_out, "vangalen_signature_violin_group.png"),
            file.path(dir_out, "vangalen_signature_violin_group.pdf"), angle = 0)
rm(vg); gc()

## ================= BoneMarrowMap (healthy reference) =================
## scored but not separately plotted -- used as the large independent
## healthy baseline in 06_cross_atlas_comparison.R
cat("\n=== BoneMarrowMap ===\n")
bmm <- readRDS(file.path(dir_data, "bonemarrowmap_seurat.rds"))
bmm <- score_atlas(bmm)
saveRDS(bmm, file.path(dir_data, "bonemarrowmap_seurat.rds"))
cat(sprintf("BoneMarrowMap scored: %d cells\n", ncol(bmm)))
rm(bmm); gc()

cat("\nDone: all four atlases scored with the LT-PLC up_all signature.\n")
