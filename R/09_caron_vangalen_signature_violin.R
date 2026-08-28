## 09_caron_vangalen_signature_violin.R
## Score the LT-PLC "up_all" human signature (UCell) per cell in Caron et al.
## 2020 (B-ALL/T-ALL) and van Galen et al. 2019 (AML), mirroring
## 05_granja_signature_violin.R, then show the score distribution per cell
## type/state in each atlas.
##
## Input:  data/caron_seurat.rds, data/vangalen_seurat.rds
##         data/human_signature_sets.rds (maximal mapping; up_all set)
## Output: results/caron_signature_violin.png / .pdf
##         results/vangalen_signature_violin.png / .pdf

suppressPackageStartupMessages({
  library(Seurat); library(UCell); library(ggplot2); library(dplyr)
})

proj_root <- getwd()
dir_data  <- file.path(proj_root, "data")
dir_out   <- file.path(proj_root, "results")

human_sets <- readRDS(file.path(dir_data, "human_signature_sets.rds"))
sig <- list(LT_PLC_up = human_sets$up_all)

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

## ================= Caron 2020 (B-ALL / T-ALL) =================
cat("=== Caron 2020 ===\n")
caron <- readRDS(file.path(dir_data, "caron_seurat.rds"))
caron <- score_atlas(caron)
saveRDS(caron, file.path(dir_data, "caron_seurat.rds"))

df_caron_ct <- make_violin(caron, "celltype",
  "Caron et al. 2020: LT-PLC score by cell/sample type",
  file.path(dir_out, "caron_signature_violin_celltype.png"),
  file.path(dir_out, "caron_signature_violin_celltype.pdf"), angle = 40)

df_caron_grp <- make_violin(caron, "group",
  "Caron et al. 2020: LT-PLC score by disease group",
  file.path(dir_out, "caron_signature_violin_group.png"),
  file.path(dir_out, "caron_signature_violin_group.pdf"), angle = 0)

rm(caron); gc()

## ================= van Galen 2019 (AML) =================
cat("\n=== van Galen 2019 ===\n")
vg <- readRDS(file.path(dir_data, "vangalen_seurat.rds"))
vg <- score_atlas(vg)
saveRDS(vg, file.path(dir_data, "vangalen_seurat.rds"))

df_vg_ct <- make_violin(vg, "celltype",
  "van Galen et al. 2019: LT-PLC score by cell state",
  file.path(dir_out, "vangalen_signature_violin_celltype.png"),
  file.path(dir_out, "vangalen_signature_violin_celltype.pdf"), angle = 45)

df_vg_grp <- make_violin(vg, "group",
  "van Galen et al. 2019: LT-PLC score by disease group",
  file.path(dir_out, "vangalen_signature_violin_group.png"),
  file.path(dir_out, "vangalen_signature_violin_group.pdf"), angle = 0)

rm(vg); gc()

cat("\nDone.\n")
