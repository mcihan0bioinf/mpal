## 06_cross_atlas_comparison.R
## Question C: across leukemia entities, where is the LT-PLC program
## strongest -- AML, B-ALL, T-ALL, or MPAL? This is the direct, headline
## answer to the reviewer's request (add MPAL, split ALL into B-ALL/T-ALL).
##
## Pools per-cell LT-PLC UCell scores (up_all signature, already computed
## by 04_score_signature_all_atlases.R) across the four single-cell atlases,
## using each atlas's own healthy reference as the local baseline:
##   Granja 2019    -> MPAL, Healthy (PBMC/BMMC/CD34 reference)
##   Caron 2020     -> B-ALL, T-ALL, Healthy (PBMMC)
##   van Galen 2019 -> AML, Healthy (BM), CellLine
##   BoneMarrowMap  -> Healthy (independent large healthy reference)
##
## Both a per-cell view (all cells pooled) and a per-patient/sample view
## (median score per sample -- the correct unit for statistics) are shown,
## since patient is the true biological replicate.
##
## Input:  data/granja_seurat.rds, data/caron_seurat.rds,
##         data/vangalen_seurat.rds, data/bonemarrowmap_seurat.rds
##         (all must already carry $LT_PLC_up_UCell -- run
##          04_score_signature_all_atlases.R first)
## Output: results/cross_atlas_comparison_percell.png / .pdf
##         results/cross_atlas_comparison_perpatient.png / .pdf
##         results/cross_atlas_patient_scores.csv

suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(dplyr)
})

proj_root <- getwd()
dir_data  <- file.path(proj_root, "data")
dir_out   <- file.path(proj_root, "results")

extract_scores <- function(path, group_col = "group", sample_col = "sample", atlas_name) {
  seu <- readRDS(path)
  if (!"LT_PLC_up_UCell" %in% colnames(seu@meta.data))
    stop(paste(path, "not scored yet -- run 04_score_signature_all_atlases.R first"))
  df <- data.frame(
    atlas = atlas_name,
    score = seu$LT_PLC_up_UCell,
    group = as.character(seu[[group_col, drop = TRUE]]),
    sample = as.character(seu[[sample_col, drop = TRUE]])
  )
  rm(seu); gc(verbose = FALSE)
  df
}

df_granja   <- extract_scores(file.path(dir_data, "granja_seurat.rds"), atlas_name = "Granja")
df_caron    <- extract_scores(file.path(dir_data, "caron_seurat.rds"), atlas_name = "Caron")
df_vangalen <- extract_scores(file.path(dir_data, "vangalen_seurat.rds"), atlas_name = "van Galen")
df_bmm      <- extract_scores(file.path(dir_data, "bonemarrowmap_seurat.rds"), atlas_name = "BoneMarrowMap")

df_all <- bind_rows(df_granja, df_caron, df_vangalen, df_bmm) %>%
  filter(!is.na(group))

## harmonize group labels: keep each atlas's own healthy reference separate
## (different platforms/pipelines -- pooling them would hide batch effects)
df_all <- df_all %>% mutate(
  group_label = case_when(
    group == "Healthy" ~ paste0("Healthy (", atlas, ")"),
    TRUE ~ group
  )
)

write.csv(df_all %>% group_by(atlas, group, sample) %>%
            summarise(median_score = median(score), n_cells = n(), .groups = "drop"),
          file.path(dir_out, "cross_atlas_patient_scores.csv"), row.names = FALSE)

## ---- disease-group order for the main comparison ----
disease_order <- c("Healthy (Caron)", "Healthy (van Galen)", "Healthy (Granja)", "Healthy (BoneMarrowMap)",
                    "CellLine", "AML", "B-ALL", "T-ALL", "MPAL")
df_plot <- df_all %>% filter(group_label %in% disease_order) %>%
  mutate(group_label = factor(group_label, levels = disease_order))

pal <- c("Healthy (Caron)" = "grey75", "Healthy (van Galen)" = "grey65", "Healthy (Granja)" = "grey55",
         "Healthy (BoneMarrowMap)" = "grey45", "CellLine" = "#756bb1",
         "AML" = "#e6550d", "B-ALL" = "#3182bd", "T-ALL" = "#31a354", "MPAL" = "#c51b8a")

## ---- per-cell view ----
p_cell <- ggplot(df_plot, aes(x = group_label, y = score, fill = group_label)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.85) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", alpha = 0.6) +
  scale_fill_manual(values = pal, guide = "none") +
  labs(title = "Cross-atlas comparison: LT-PLC signature score by disease entity (per cell)",
       subtitle = "UCell, up_all signature; each atlas's own healthy reference shown alongside its disease group(s)",
       x = NULL, y = "UCell score (LT-PLC up_all)") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 9, color = "grey30"),
        panel.grid.minor = element_blank())

ggsave(file.path(dir_out, "cross_atlas_comparison_percell.png"), p_cell, width = 11, height = 6.5, dpi = 300)
ggsave(file.path(dir_out, "cross_atlas_comparison_percell.pdf"), p_cell, width = 11, height = 6.5)

## ---- per-patient/sample view (correct statistical unit) ----
df_patient <- df_plot %>% group_by(atlas, group_label, sample) %>%
  summarise(median_score = median(score), n_cells = n(), .groups = "drop") %>%
  filter(n_cells >= 20)

p_patient <- ggplot(df_patient, aes(x = group_label, y = median_score, fill = group_label)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.85) +
  geom_jitter(width = 0.12, size = 1.8, alpha = 0.7) +
  scale_fill_manual(values = pal, guide = "none") +
  labs(title = "Cross-atlas comparison: LT-PLC signature score by disease entity (per patient/sample)",
       subtitle = "one point = median UCell score across cells for one sample (min 20 cells); patient is the biological replicate",
       x = NULL, y = "Median UCell score per sample") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 9, color = "grey30"),
        panel.grid.minor = element_blank())

ggsave(file.path(dir_out, "cross_atlas_comparison_perpatient.png"), p_patient, width = 11, height = 6.5, dpi = 300)
ggsave(file.path(dir_out, "cross_atlas_comparison_perpatient.pdf"), p_patient, width = 11, height = 6.5)

cat("\nSaved cross_atlas_comparison_percell.png/.pdf, cross_atlas_comparison_perpatient.png/.pdf, cross_atlas_patient_scores.csv\n")

cat("\nMedian per-cell score by group:\n")
print(df_plot %>% group_by(group_label) %>% summarise(median = median(score), n_cells = n()) %>% arrange(desc(median)))

cat("\nMedian per-patient score by group (n = number of samples, min 20 cells each):\n")
print(df_patient %>% group_by(group_label) %>% summarise(median = median(median_score), n_samples = n()) %>% arrange(desc(median)))

## disease vs its own matched healthy reference, Wilcoxon on per-patient medians
cat("\nDisease vs local healthy reference, per-patient Wilcoxon:\n")
compare_pairs <- list(
  AML   = c("AML", "Healthy (van Galen)"),
  BALL  = c("B-ALL", "Healthy (Caron)"),
  TALL  = c("T-ALL", "Healthy (Caron)"),
  MPAL  = c("MPAL", "Healthy (Granja)")
)
for (nm in names(compare_pairs)) {
  g <- compare_pairs[[nm]]
  a <- df_patient$median_score[df_patient$group_label == g[1]]
  b <- df_patient$median_score[df_patient$group_label == g[2]]
  if (length(a) >= 2 && length(b) >= 2) {
    wt <- wilcox.test(a, b)
    cat(sprintf("%-6s (n=%d) vs %-25s (n=%d): median %.4f vs %.4f, Wilcoxon p=%.3f\n",
                g[1], length(a), g[2], length(b), median(a), median(b), wt$p.value))
  }
}
