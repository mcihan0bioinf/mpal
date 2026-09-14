## 10_entity_pseudobulk_comparison.R
## Direct extension of the manuscript's Fig. 6A-B logic (mouse signature vs.
## human leukemia entities, correlated/scored per patient, compared with a
## rank test) rather than a new method: same idea (patient is the unit,
## entity is the grouping variable), same UCell instrument already used
## throughout this pipeline, just pointed at four disease entities instead
## of cell types within one atlas.
##
## Three changes from the earlier per-cell/per-cell-type analyses:
##   1. Score per-patient PSEUDOBULK (sum of counts across a patient's
##      cells), not per cell -- removes most of the library-depth confound.
##   2. Group by disease ENTITY (AML, B-ALL, T-ALL, MPAL, Normal), not by
##      cell type/state.
##   3. Add a size-matched random-gene-set null per patient, so the
##      absolute UCell value has a reference point ("higher than a random
##      gene set of the same size would give"), reported as a z-score.
##
## Headline signature: up_lfc1 (padj<0.05 & LFC>=1) -- a more defensible
## core set than up_all (2380 genes is large enough to partly reflect
## general transcriptional output). All seven variants from
## 09_granja_multisig_robustness.R are still scored side by side for
## robustness.
##
## Atlases / entities (four entities the manuscript needs, Huang 2026 and
## Zeng 2025 still not obtained -- see rebuttal doc TBD notes):
##   Granja 2019    -> MPAL (all 6 disease-classified compartments pooled
##                      per patient), Normal (4 healthy donors)
##   Caron 2020     -> B-ALL, T-ALL, Normal (PBMMC)
##   van Galen 2019 -> AML, Normal (BM); CellLine excluded (not a patient)
##   BoneMarrowMap  -> Normal (45 independent healthy donors)
##
## Input:  data/granja_seurat.rds, data/caron_seurat.rds,
##         data/vangalen_seurat.rds, data/bonemarrowmap_seurat.rds
##         data/human_signature_sets.rds (maximal mapping)
## Output: results/entity_pseudobulk_scores.csv
##         results/entity_pseudobulk_headline.png / .pdf   (up_lfc1 + null z)
##         results/entity_pseudobulk_multisig.png / .pdf   (all 7 variants)

suppressPackageStartupMessages({
  library(Seurat); library(UCell); library(ggplot2); library(dplyr); library(tidyr); library(patchwork)
})

set.seed(1234)
proj_root <- getwd()
dir_data  <- file.path(proj_root, "data")
dir_out   <- file.path(proj_root, "results")

human_sets <- readRDS(file.path(dir_data, "human_signature_sets.rds"))

variant_order <- c("up_all", "up_lfc1", "up_lfc2", "up_top50", "up_top100", "up_top500", "up_and_down")
variant_labels <- c(
  up_all      = "Up: all (padj<0.05)",
  up_lfc1     = "Up: LFC>=1 (headline)",
  up_lfc2     = "Up: LFC>=2",
  up_top50    = "Up: top 50",
  up_top100   = "Up: top 100",
  up_top500   = "Up: top 500",
  up_and_down = "Up (all) + Down (LFC<=-1)"
)
variants <- list(
  up_all      = human_sets$up_all,
  up_lfc1     = human_sets$up_lfc1,
  up_lfc2     = human_sets$up_lfc2,
  up_top50    = human_sets$up_top50,
  up_top100   = human_sets$up_top100,
  up_top500   = human_sets$up_top500,
  up_and_down = c(human_sets$up_all, paste0(human_sets$down_lfc1, "-"))
)

entity_order <- c("Normal", "AML", "B-ALL", "T-ALL", "MPAL")
entity_pal <- c(Normal = "grey60", AML = "#e6550d", `B-ALL` = "#3182bd",
                 `T-ALL` = "#31a354", MPAL = "#c51b8a")

## ================= pseudobulk helper =================
## sums raw counts across each patient's cells; drops patients below
## min_cells (pseudobulk from a handful of cells is not trustworthy)
pseudobulk_atlas <- function(seu, entity_col, patient, entity_map, min_cells = 50) {
  entity <- entity_map[as.character(seu@meta.data[[entity_col]])]
  keep <- !is.na(entity)
  seu <- seu[, keep]; entity <- entity[keep]; patient <- patient[keep]

  counts <- GetAssayData(seu, assay = "RNA", layer = "counts")
  patient_f <- factor(patient)
  ## sparse group-sum via indicator matrix: counts (genes x cells) %*% (cells x patients)
  ind <- Matrix::sparse.model.matrix(~ 0 + patient_f)
  colnames(ind) <- levels(patient_f)
  pb <- as.matrix(counts %*% ind)  # genes x patients

  n_cells <- table(patient_f)
  meta <- data.frame(patient = levels(patient_f),
                      entity = entity[match(levels(patient_f), patient)],
                      n_cells = as.integer(n_cells[levels(patient_f)]))
  meta <- meta %>% filter(n_cells >= min_cells)
  pb <- pb[, meta$patient, drop = FALSE]
  list(matrix = pb, meta = meta)
}

## size-matched random-gene-set null -> per-patient z-score for one signature
null_z <- function(pb_matrix, real_scores, sig_genes, K = 100, seed = 1) {
  background <- rownames(pb_matrix)[Matrix::rowSums(pb_matrix) > 0]
  n <- length(intersect(sig_genes, rownames(pb_matrix)))
  if (n < 5 || n > length(background)) return(rep(NA_real_, ncol(pb_matrix)))
  set.seed(seed)
  rand_sets <- lapply(seq_len(K), function(i) sample(background, n))
  names(rand_sets) <- paste0("rand", seq_len(K))
  rand_scores <- ScoreSignatures_UCell(pb_matrix, features = rand_sets, maxRank = 4500)
  mu <- rowMeans(rand_scores); sdv <- apply(rand_scores, 1, sd)
  (real_scores - mu) / sdv
}

## ================= build per-atlas pseudobulk matrices =================
cat("=== Granja 2019 ===\n")
granja <- readRDS(file.path(dir_data, "granja_seurat.rds"))
g_patient <- sub("_?T[0-9]+$", "", granja$sample)
pb_granja <- pseudobulk_atlas(granja, "group", g_patient,
                               c(MPAL = "MPAL", Healthy = "Normal"))
cat(sprintf("  %d patients kept: %s\n", nrow(pb_granja$meta), paste(pb_granja$meta$patient, collapse=", ")))
rm(granja); gc(verbose = FALSE)

cat("\n=== Caron 2020 ===\n")
caron <- readRDS(file.path(dir_data, "caron_seurat.rds"))
pb_caron <- pseudobulk_atlas(caron, "group", caron$sample,
                              c(`B-ALL` = "B-ALL", `T-ALL` = "T-ALL", Healthy = "Normal"))
cat(sprintf("  %d patients kept: %s\n", nrow(pb_caron$meta), paste(pb_caron$meta$patient, collapse=", ")))
rm(caron); gc(verbose = FALSE)

cat("\n=== van Galen 2019 ===\n")
vg <- readRDS(file.path(dir_data, "vangalen_seurat.rds"))
vg_patient <- sub("-D[0-9]+$", "", vg$orig.ident)  # AML707B-D113 -> AML707B; BM1 unchanged
pb_vg <- pseudobulk_atlas(vg, "group", vg_patient,
                           c(AML = "AML", Healthy = "Normal"))  # CellLine excluded (not in map -> NA -> dropped)
cat(sprintf("  %d patients kept: %s\n", nrow(pb_vg$meta), paste(pb_vg$meta$patient, collapse=", ")))
rm(vg); gc(verbose = FALSE)

cat("\n=== BoneMarrowMap ===\n")
bmm <- readRDS(file.path(dir_data, "bonemarrowmap_seurat.rds"))
pb_bmm <- pseudobulk_atlas(bmm, "group", bmm$sample, c(Healthy = "Normal"))
cat(sprintf("  %d patients kept\n", nrow(pb_bmm$meta)))
rm(bmm); gc(verbose = FALSE)

## ================= score each atlas's pseudobulk matrix =================
score_pb <- function(pb, atlas_name) {
  sc <- ScoreSignatures_UCell(pb$matrix, features = variants, maxRank = 4500)
  colnames(sc) <- sub("_UCell$", "", colnames(sc))
  z <- null_z(pb$matrix, sc[, "up_lfc1"], variants$up_lfc1, K = 100, seed = 42)
  cbind(atlas = atlas_name, pb$meta, as.data.frame(sc), z_up_lfc1 = z)
}

df_all <- bind_rows(
  score_pb(pb_granja, "Granja"),
  score_pb(pb_caron,  "Caron"),
  score_pb(pb_vg,     "van Galen"),
  score_pb(pb_bmm,    "BoneMarrowMap")
) %>% mutate(entity = factor(entity, levels = entity_order))

write.csv(df_all, file.path(dir_out, "entity_pseudobulk_scores.csv"), row.names = FALSE)
cat(sprintf("\nSaved results/entity_pseudobulk_scores.csv (%d patients total)\n", nrow(df_all)))
cat("\nPatients per entity:\n"); print(table(df_all$entity))

## ================= headline figure: up_lfc1 score + size-matched null z =================
p_score <- ggplot(df_all, aes(x = entity, y = up_lfc1, fill = entity)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.85, width = 0.55) +
  geom_jitter(width = 0.12, size = 1.6, alpha = 0.7) +
  scale_fill_manual(values = entity_pal, guide = "none") +
  labs(title = "A. LT-PLC score (up_lfc1) by entity", x = NULL, y = "Pseudobulk UCell score") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 12), panel.grid.minor = element_blank())

## IMPORTANT: the random background (and so the null mean/sd) is drawn
## separately per atlas -- each atlas has its own detected-gene universe and
## its own pseudobulk depth/composition (e.g. BoneMarrowMap's donor
## pseudobulks are far deeper and more heterogeneous than Granja's sorted
## MPAL blasts). A z-score is therefore only interpretable *within* the
## atlas it was computed in -- pooling z across atlases into one scale
## (as an earlier version of this figure did) is invalid and produced a
## misleading reversal (e.g. Granja MPAL looked "more below random" than
## its own PBMC/BMMC reference, driven entirely by that atlas's healthy
## samples being compositionally heterogeneous mixtures, not by a lower
## LT-PLC signal). Facet by atlas so only atlas-matched entities/local
## Normal are ever compared on the same axis.
p_z <- ggplot(df_all %>% filter(!is.na(z_up_lfc1)), aes(x = entity, y = z_up_lfc1, fill = entity)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_boxplot(outlier.shape = NA, alpha = 0.85, width = 0.55) +
  geom_jitter(width = 0.12, size = 1.6, alpha = 0.7) +
  facet_wrap(~atlas, nrow = 1, scales = "free_x") +
  scale_fill_manual(values = entity_pal, guide = "none") +
  labs(title = "B. Size-matched null: z-score vs. 100 random gene sets of the same size",
       subtitle = "z is only comparable WITHIN a panel (own atlas background/null) -- not across atlases",
       x = NULL, y = "z = (real - mean(random)) / sd(random)") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 8.5, color = "grey30"),
        strip.text = element_text(face = "bold", size = 9),
        panel.grid.minor = element_blank())

p_headline <- p_score / p_z + patchwork::plot_annotation(
  title = "Pseudobulk (per-patient) LT-PLC score across disease entities",
  subtitle = "one point = one patient; up_lfc1 signature (padj<0.05 & LFC>=1); Huang 2026 / Zeng 2025 arms not yet obtained",
  theme = theme(plot.title = element_text(face = "bold", size = 13),
                plot.subtitle = element_text(size = 9.5, color = "grey30")))

ggsave(file.path(dir_out, "entity_pseudobulk_headline.png"), p_headline, width = 8, height = 9, dpi = 300)
ggsave(file.path(dir_out, "entity_pseudobulk_headline.pdf"), p_headline, width = 8, height = 9)
cat("\nSaved results/entity_pseudobulk_headline.png / .pdf\n")

## ================= multi-variant grid (robustness across all 7 signatures) =================
df_multisig <- df_all %>%
  select(entity, all_of(variant_order)) %>%
  pivot_longer(cols = all_of(variant_order), names_to = "variant", values_to = "score") %>%
  mutate(variant = factor(variant, levels = variant_order, labels = variant_labels[variant_order]))

p_multisig <- ggplot(df_multisig, aes(x = entity, y = score, fill = entity)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.85, width = 0.6) +
  geom_jitter(width = 0.1, size = 1, alpha = 0.6) +
  facet_wrap(~variant, ncol = 4, scales = "free_y") +
  scale_fill_manual(values = entity_pal, guide = "none") +
  labs(title = "Pseudobulk LT-PLC score by entity, across all signature variants",
       subtitle = "one point = one patient; y-axis free per panel (raw scale differs by signature size)",
       x = NULL, y = "Pseudobulk UCell score") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        strip.text = element_text(face = "bold", size = 9),
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 9, color = "grey30"),
        panel.grid.minor = element_blank())

ggsave(file.path(dir_out, "entity_pseudobulk_multisig.png"), p_multisig, width = 15, height = 8, dpi = 300)
ggsave(file.path(dir_out, "entity_pseudobulk_multisig.pdf"), p_multisig, width = 15, height = 8)
cat("Saved results/entity_pseudobulk_multisig.png / .pdf\n")

## ================= statistics: patient-level, entity vs Normal + Kruskal-Wallis =================
cat("\n=== Kruskal-Wallis across all 5 entities (up_lfc1, patient-level) ===\n")
kw <- kruskal.test(up_lfc1 ~ entity, data = df_all)
print(kw)

cat("\n=== Each entity vs Normal, patient-level Wilcoxon (up_lfc1) ===\n")
for (e in setdiff(entity_order, "Normal")) {
  a <- df_all$up_lfc1[df_all$entity == e]
  b <- df_all$up_lfc1[df_all$entity == "Normal"]
  if (length(a) >= 2 && length(b) >= 2) {
    wt <- wilcox.test(a, b)
    cat(sprintf("%-6s (n=%d) vs Normal (n=%d): median %.4f vs %.4f, Wilcoxon p=%.4f\n",
                e, length(a), length(b), median(a), median(b), wt$p.value))
  }
}

cat("\n=== MPAL vs B-ALL / T-ALL / AML, patient-level Wilcoxon (up_lfc1) ===\n")
mpal <- df_all$up_lfc1[df_all$entity == "MPAL"]
for (e in c("AML", "B-ALL", "T-ALL")) {
  b <- df_all$up_lfc1[df_all$entity == e]
  if (length(mpal) >= 2 && length(b) >= 2) {
    wt <- wilcox.test(mpal, b)
    cat(sprintf("MPAL (n=%d) vs %-6s (n=%d): median %.4f vs %.4f, Wilcoxon p=%.4f\n",
                length(mpal), e, length(b), median(mpal), median(b), wt$p.value))
  }
}
