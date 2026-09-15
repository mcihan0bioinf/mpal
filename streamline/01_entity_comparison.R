## 01_entity_comparison.R
## Comment 1 (Reviewer 1, key analysis): "not clear why the comparison does not
## include MPAL... no distinction is made between B- vs T-ALL". Re-run the
## AML / B-ALL / T-ALL / MPAL / Normal comparison with the LT-PLC "up,
## padj<0.05 & LFC>=1" signature (Table S5), mapped to human orthologs with
## the STRICT one-to-one mapping (not the maximal mapping used elsewhere in
## this repo's R/ pipeline).
##
## Every atlas is scored the same way -- UCell per cell, then median per
## patient/sample -- so AML/B-ALL/T-ALL/MPAL/Normal are the same kind of
## quantity everywhere (no mixing with pseudobulk). This also makes the
## per-patient scores directly comparable to 02_granja_mpal_patient_level.R's
## Section A, which scores Granja the identical way.
##
## Entity and dataset are partly confounded (AML only from van Galen, B/T-ALL
## only from Caron, MPAL only from Granja -- different platforms/pipelines),
## so raw cross-atlas scores are shown but not over-interpreted. The primary
## comparison is a dataset-normalized enrichment: each patient's median score
## minus that same dataset's own Normal median. B-ALL vs T-ALL (same dataset,
## Caron) and MPAL vs its own Normal (same dataset, Granja) are the
## least-confounded comparisons and are called out separately from
## cross-dataset ones (e.g. AML vs MPAL), which are reported with a caveat.
##
## MPAL5R is the relapse sample of patient MPAL5, not a 6th independent
## patient -- it is excluded from the primary n=5 MPAL set and reported
## separately as a sensitivity point.
##
## Input:  data/granja_seurat.rds, data/caron_seurat.rds,
##         data/vangalen_seurat.rds, data/bonemarrowmap_seurat.rds
##         data/human_signature_sets_strict.rds (1:1 orthologs; up_lfc1)
## Output: streamline/results/entity_comparison_raw.png / .pdf
##         streamline/results/entity_comparison_normalized.png / .pdf
##         streamline/results/entity_patient_scores.csv

suppressPackageStartupMessages({
  library(Seurat); library(UCell); library(ggplot2); library(dplyr); library(patchwork)
})

proj_root <- getwd()
dir_data  <- file.path(proj_root, "data")
dir_out   <- file.path(proj_root, "streamline", "results")
if (!dir.exists(dir_out)) dir.create(dir_out, recursive = TRUE)

human_strict <- readRDS(file.path(dir_data, "human_signature_sets_strict.rds"))
sig <- list(LT_PLC_up_lfc1 = human_strict$up_lfc1)
cat(sprintf("Signature: up_lfc1, strict 1:1 human orthologs, %d genes\n\n", length(sig$LT_PLC_up_lfc1)))

## ---- shared helper: UCell-score a raw counts matrix per cell, in memory-bounded
## batches. Takes the counts matrix directly (not a Seurat object) and slices it
## by column -- Seurat's subset() on a Seurat object recreates a new object with
## its full metadata/graph overhead on every call, which is far more memory-
## hungry than slicing the sparse matrix directly, and was enough to OOM-kill
## this 15GB machine on BoneMarrowMap's 263k cells. batch_size is capped lower
## for very large atlases to keep peak per-batch memory down. ----
score_matrix <- function(counts, batch_size = 4000) {
  n_present <- sum(sig$LT_PLC_up_lfc1 %in% rownames(counts))
  cat(sprintf("  gene coverage: %d / %d (%.0f%%) signature genes detected in this atlas\n",
              n_present, length(sig$LT_PLC_up_lfc1), 100 * n_present / length(sig$LT_PLC_up_lfc1)))
  cells <- colnames(counts)
  batches <- split(cells, ceiling(seq_along(cells) / batch_size))
  scores_list <- vector("list", length(batches))
  for (i in seq_along(batches)) {
    sc <- ScoreSignatures_UCell(counts[, batches[[i]], drop = FALSE], features = sig, maxRank = 3000)
    scores_list[[i]] <- sc
    gc(verbose = FALSE, full = TRUE)
    if (i %% 5 == 0 || i == length(batches)) cat(sprintf("    batch %d/%d scored\n", i, length(batches)))
  }
  scores <- do.call(rbind, scores_list)[cells, , drop = FALSE]
  as.numeric(scores[, "LT_PLC_up_lfc1_UCell"])
}

## ---- per-patient median helper ----
patient_medians <- function(score, patient, entity, min_cells = 20) {
  data.frame(score = score, patient = patient, entity = entity) %>%
    filter(!is.na(entity)) %>%
    group_by(patient, entity) %>%
    summarise(median_score = median(score), n_cells = n(), .groups = "drop") %>%
    filter(n_cells >= min_cells)
}

## Each block below: read the Seurat object, immediately pull out the raw counts
## matrix and the small bits of metadata needed (sample/patient, group), then rm()
## the Seurat object itself before scoring -- so only the counts matrix (not the
## whole object plus its graphs/reductions/etc.) is resident during UCell scoring.

## ================= Granja 2019 (MPAL, Normal) =================
cat("=== Granja 2019 ===\n")
granja <- readRDS(file.path(dir_data, "granja_seurat.rds"))
DefaultAssay(granja) <- "RNA"
granja_counts <- GetAssayData(granja, layer = "counts")
granja_patient <- sub("_?T[0-9]+$", "", granja$sample)
granja_entity <- ifelse(granja$group == "Healthy", "Normal", ifelse(granja$group == "MPAL", "MPAL", NA))
rm(granja); gc(verbose = FALSE, full = TRUE)

granja_score <- score_matrix(granja_counts)
rm(granja_counts); gc(verbose = FALSE, full = TRUE)

df_granja <- patient_medians(granja_score, granja_patient, granja_entity) %>% mutate(dataset = "Granja")
## split MPAL5 (primary) from MPAL5R (relapse -- not an independent patient)
df_granja <- df_granja %>% mutate(
  is_relapse = patient == "MPAL5R",
  patient_primary = ifelse(is_relapse, "MPAL5", patient)
)

## ================= Caron 2020 (B-ALL, T-ALL, Normal) =================
cat("\n=== Caron 2020 ===\n")
caron <- readRDS(file.path(dir_data, "caron_seurat.rds"))
DefaultAssay(caron) <- "RNA"
caron_counts <- GetAssayData(caron, layer = "counts")
caron_sample <- caron$sample
caron_entity <- ifelse(caron$group == "Healthy", "Normal", ifelse(caron$group %in% c("B-ALL", "T-ALL"), as.character(caron$group), NA))
rm(caron); gc(verbose = FALSE, full = TRUE)

caron_score <- score_matrix(caron_counts)
rm(caron_counts); gc(verbose = FALSE, full = TRUE)

df_caron <- patient_medians(caron_score, caron_sample, caron_entity) %>%
  mutate(dataset = "Caron", is_relapse = FALSE, patient_primary = patient)

## ================= van Galen 2019 (AML, Normal) =================
cat("\n=== van Galen 2019 ===\n")
vg <- readRDS(file.path(dir_data, "vangalen_seurat.rds"))
DefaultAssay(vg) <- "RNA"
vg_counts <- GetAssayData(vg, layer = "counts")
vg_patient <- sub("-D[0-9]+$", "", vg$orig.ident)
vg_entity <- ifelse(vg$group == "Healthy", "Normal", ifelse(vg$group == "AML", "AML", NA))  # CellLine excluded
rm(vg); gc(verbose = FALSE, full = TRUE)

vg_score <- score_matrix(vg_counts)
rm(vg_counts); gc(verbose = FALSE, full = TRUE)

df_vg <- patient_medians(vg_score, vg_patient, vg_entity) %>%
  mutate(dataset = "van Galen", is_relapse = FALSE, patient_primary = patient)

## ================= BoneMarrowMap (independent Normal reference) =================
cat("\n=== BoneMarrowMap ===\n")
bmm <- readRDS(file.path(dir_data, "bonemarrowmap_seurat.rds"))
DefaultAssay(bmm) <- "RNA"
bmm_counts <- GetAssayData(bmm, layer = "counts")
bmm_sample <- bmm$sample
bmm_entity <- ifelse(bmm$group == "Healthy", "Normal", NA)
rm(bmm); gc(verbose = FALSE, full = TRUE)

bmm_score <- score_matrix(bmm_counts, batch_size = 2500)  # largest atlas -- smaller batches
rm(bmm_counts); gc(verbose = FALSE, full = TRUE)

df_bmm <- patient_medians(bmm_score, bmm_sample, bmm_entity) %>%
  mutate(dataset = "BoneMarrowMap", is_relapse = FALSE, patient_primary = patient)

## ================= combine =================
df_all <- bind_rows(df_granja, df_caron, df_vg, df_bmm)

## primary set: drop relapse samples from the independent-patient pool
df_primary <- df_all %>% filter(!is_relapse)

entity_order <- c("Normal", "AML", "B-ALL", "T-ALL", "MPAL")
entity_pal <- c(Normal = "grey60", AML = "#e6550d", `B-ALL` = "#3182bd", `T-ALL` = "#31a354", MPAL = "#c51b8a")

cat("\nPatients per dataset x entity (primary set, relapse excluded):\n")
print(df_primary %>% count(dataset, entity) %>% arrange(dataset, entity))

## ---- dataset-normalized enrichment: patient median - that dataset's own Normal median ----
normal_ref <- df_primary %>% filter(entity == "Normal") %>%
  group_by(dataset) %>% summarise(normal_median = median(median_score), normal_sd = sd(median_score), .groups = "drop")

df_norm <- df_primary %>% left_join(normal_ref, by = "dataset") %>%
  mutate(enrichment = median_score - normal_median,
         z_vs_normal = (median_score - normal_median) / normal_sd)

write.csv(df_all %>% left_join(normal_ref, by = "dataset") %>%
            mutate(enrichment = median_score - normal_median, z_vs_normal = (median_score - normal_median) / normal_sd),
          file.path(dir_out, "entity_patient_scores.csv"), row.names = FALSE)
cat("\nSaved streamline/results/entity_patient_scores.csv (includes relapse row, flagged is_relapse=TRUE)\n")

## ================= Figure A: raw per-patient scores, colored by dataset =================
df_plot <- df_primary %>% mutate(entity = factor(entity, levels = entity_order))

p_raw <- ggplot(df_plot, aes(x = entity, y = median_score, color = dataset)) +
  geom_boxplot(aes(group = entity), outlier.shape = NA, alpha = 0, color = "grey50", width = 0.5) +
  geom_jitter(width = 0.12, size = 2, alpha = 0.85) +
  labs(title = "A. Raw per-patient LT-PLC score (up_lfc1, strict 1:1 orthologs) by entity",
       subtitle = "one point = one patient (median UCell score, ≥20 cells); color = source dataset -- entity and\ndataset are confounded (AML only from van Galen, B/T-ALL only from Caron, MPAL only from Granja),\nso only within-dataset comparisons (e.g. B-ALL vs T-ALL) are directly interpretable here",
       x = NULL, y = "Median UCell score per patient", color = "Dataset") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12.5),
        plot.subtitle = element_text(size = 8.5, color = "grey30"),
        panel.grid.minor = element_blank())

## ================= Figure B: dataset-normalized enrichment =================
p_norm <- ggplot(df_norm %>% mutate(entity = factor(entity, levels = entity_order)),
                  aes(x = entity, y = enrichment, fill = entity)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_boxplot(outlier.shape = NA, alpha = 0.85, width = 0.55) +
  geom_jitter(width = 0.12, size = 2, alpha = 0.8) +
  scale_fill_manual(values = entity_pal, guide = "none") +
  labs(title = "B. Dataset-normalized enrichment (patient median − that dataset's own Normal median)",
       subtitle = "removes cross-platform offset; each entity is shown relative to its own dataset's healthy baseline",
       x = NULL, y = "Enrichment over dataset-matched Normal") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12.5),
        plot.subtitle = element_text(size = 8.5, color = "grey30"),
        panel.grid.minor = element_blank())

p_combined <- p_raw / p_norm + plot_annotation(
  title = "LT-PLC signature score across human leukemia entities (AML / B-ALL / T-ALL / MPAL / Normal)",
  subtitle = "Table S5 up-signature (padj<0.05, LFC≥1), strict 1:1 human orthologs, UCell per cell, median per patient; MPAL5R (relapse) excluded from n",
  theme = theme(plot.title = element_text(face = "bold", size = 13),
                plot.subtitle = element_text(size = 9, color = "grey30")))

ggsave(file.path(dir_out, "entity_comparison_raw.png"), p_raw, width = 9, height = 5.5, dpi = 300)
ggsave(file.path(dir_out, "entity_comparison_raw.pdf"), p_raw, width = 9, height = 5.5)
ggsave(file.path(dir_out, "entity_comparison_normalized.png"), p_norm, width = 9, height = 5.5, dpi = 300)
ggsave(file.path(dir_out, "entity_comparison_normalized.pdf"), p_norm, width = 9, height = 5.5)
ggsave(file.path(dir_out, "entity_comparison_combined.png"), p_combined, width = 9, height = 10, dpi = 300)
ggsave(file.path(dir_out, "entity_comparison_combined.pdf"), p_combined, width = 9, height = 10)
cat("\nSaved streamline/results/entity_comparison_raw / _normalized / _combined .png/.pdf\n")

## ================= statistics =================
cat("\nMedian raw score and normalized enrichment by entity (across all datasets pooled -- note confound):\n")
print(df_norm %>% group_by(entity) %>%
        summarise(n_patients = n(), median_raw = median(median_score), median_enrichment = median(enrichment)) %>%
        arrange(desc(median_enrichment)))

get_scores <- function(df, ent) df$median_score[df$entity == ent]

cat("\n--- Least-confounded, within-dataset comparisons ---\n")
cat("\nB-ALL vs T-ALL (both from Caron), Wilcoxon on raw per-patient median:\n")
a <- get_scores(df_primary %>% filter(dataset == "Caron"), "B-ALL")
b <- get_scores(df_primary %>% filter(dataset == "Caron"), "T-ALL")
print(wilcox.test(a, b)); cat(sprintf("median B-ALL=%.4f (n=%d) vs T-ALL=%.4f (n=%d)\n", median(a), length(a), median(b), length(b)))

cat("\nMPAL vs its own Normal (both from Granja), Wilcoxon on raw per-patient median:\n")
a <- get_scores(df_primary %>% filter(dataset == "Granja"), "MPAL")
b <- get_scores(df_primary %>% filter(dataset == "Granja"), "Normal")
print(wilcox.test(a, b)); cat(sprintf("median MPAL=%.4f (n=%d) vs Granja-Normal=%.4f (n=%d)\n", median(a), length(a), median(b), length(b)))

cat("\nAML vs its own Normal (both from van Galen), Wilcoxon on raw per-patient median:\n")
a <- get_scores(df_primary %>% filter(dataset == "van Galen"), "AML")
b <- get_scores(df_primary %>% filter(dataset == "van Galen"), "Normal")
print(wilcox.test(a, b)); cat(sprintf("median AML=%.4f (n=%d) vs vanGalen-Normal=%.4f (n=%d)\n", median(a), length(a), median(b), length(b)))

cat("\n--- Cross-dataset comparisons (entity and dataset confounded -- interpret cautiously) ---\n")
cat("Kruskal-Wallis across AML/B-ALL/T-ALL/MPAL, dataset-normalized enrichment:\n")
df_kw <- df_norm %>% filter(entity != "Normal")
print(kruskal.test(enrichment ~ entity, data = df_kw))

cat("\nPairwise Wilcoxon, dataset-normalized enrichment (MPAL vs each ALL/AML entity):\n")
mpal_e <- df_norm$enrichment[df_norm$entity == "MPAL"]
for (e in c("AML", "B-ALL", "T-ALL")) {
  be <- df_norm$enrichment[df_norm$entity == e]
  if (length(mpal_e) >= 2 && length(be) >= 2) {
    wt <- wilcox.test(mpal_e, be)
    cat(sprintf("MPAL (n=%d) vs %-6s (n=%d): median enrichment %.4f vs %.4f, Wilcoxon p=%.4f [cross-dataset]\n",
                length(mpal_e), e, length(be), median(mpal_e), median(be), wt$p.value))
  }
}

## ---- MPAL5R sensitivity check ----
cat("\n--- MPAL5R (relapse) sensitivity check, not part of primary n ---\n")
print(df_all %>% filter(dataset == "Granja") %>% select(patient, entity, median_score, n_cells, is_relapse))

cat("\nDone.\n")
