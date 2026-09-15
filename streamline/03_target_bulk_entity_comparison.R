## 03_target_bulk_entity_comparison.R
## Comment 1 (Reviewer 1, Fig. 6 fix) -- the definitive version. Replaces the
## small, underpowered single-cell entity comparison in 01_entity_comparison.R
## (T-ALL n=2) with one coherent bulk RNA-seq cohort: NCI TARGET, compiled
## exactly as in Mumme et al. 2025 (Nat Commun, PedSCAtlas paper) from
## TARGET-AML + TARGET-ALL-P2 + TARGET-ALL-P3 (GDC), giving AML / B-ALL /
## T-ALL / B-My-MPAL / T-My-MPAL / Normal all from ONE processing pipeline
## (same GDC STAR-Counts workflow, same TPM normalization) -- no cross-study
## platform confound like the scRNA-seq comparison had.
##
## Entity labels come from each project's own open-access GDC "Clinical
## Supplement" Excel file (verified during development, not guessed):
##   TARGET-AML            -> entity = AML (primary tumor); Normal (BM/blood normal)
##   TARGET-ALL-P2 Validation clinical file, "Cell of Origin" column
##                          -> B Cell ALL / T Cell ALL
##   TARGET-ALL-P3 clinical file (the ALAL/MPAL-focused phase), "WHO Dx" column
##                          -> B/M = B-My-MPAL, T/M = T-My-MPAL
##
## Scoring: NOT UCell (bulk, not single-cell) -- same dependency-light method
## already used for Mulet-Lazaro (07_mulet_lazaro_alal_bulk.R): mean of
## per-gene z-scores (log2(TPM+1), z-scored across all patients in this
## cohort) over the strict 1:1-ortholog up_lfc1 signature genes present in
## the data. This is the "same basic type of patient-level signature
## enrichment framework as the existing Fig. 6", per the agreed plan.
##
## Downloads are capped per entity (AML/Normal/B-ALL/T-ALL: n=150 random
## subsample, seeded) to keep this tractable on a resource-constrained
## machine; both MPAL subtypes are small in TARGET so ALL available
## diagnosis samples with RNA-seq are used. Everything is cached to
## data/target_bulk/ so reruns after the first do not re-download.
##
## Input:  GDC API (public, open-access STAR-Counts gene expression + open
##         clinical supplement Excel files -- verified open access for all
##         three projects during development)
##         data/human_signature_sets_strict.rds (1:1 orthologs; up_lfc1)
## Output: data/target_bulk/target_bulk_expr.rds, target_bulk_meta.rds (cache)
##         streamline/results/target_entity_scores.csv
##         streamline/results/target_entity_comparison.png / .pdf

suppressPackageStartupMessages({
  library(httr); library(jsonlite); library(readxl); library(dplyr); library(ggplot2)
})

set.seed(1234)
proj_root <- getwd()
dir_data  <- file.path(proj_root, "data", "target_bulk")
dir_out   <- file.path(proj_root, "streamline", "results")
for (d in c(dir_data, dir_out)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

expr_cache <- file.path(dir_data, "target_bulk_expr.rds")
meta_cache <- file.path(dir_data, "target_bulk_meta.rds")

GDC_FILES <- "https://api.gdc.cancer.gov/files?format=JSON"
GDC_DATA  <- "https://api.gdc.cancer.gov/data"

gdc_query <- function(filters, fields, size = 5000) {
  r <- POST(GDC_FILES,
            body = toJSON(list(filters = filters, fields = fields, size = size), auto_unbox = TRUE),
            content_type_json(), accept_json())
  stop_for_status(r)
  d <- fromJSON(content(r, as = "text", encoding = "UTF-8"), simplifyVector = FALSE)
  d$data$hits
}

filt_project_type_sampletype <- function(project, sample_types) {
  list(op = "and", content = list(
    list(op = "in", content = list(field = "cases.project.project_id", value = list(project))),
    list(op = "in", content = list(field = "data_type", value = list("Gene Expression Quantification"))),
    list(op = "in", content = list(field = "analysis.workflow_type", value = list("STAR - Counts"))),
    list(op = "in", content = list(field = "cases.samples.sample_type", value = as.list(sample_types)))
  ))
}

## file_id + case submitter_id (TARGET USI) + sample_type for one project
list_files <- function(project, sample_types) {
  hits <- gdc_query(filt_project_type_sampletype(project, sample_types),
                     fields = "file_id,cases.submitter_id,cases.samples.sample_type")
  do.call(rbind, lapply(hits, function(h) {
    data.frame(file_id = h$file_id,
               case_id = h$cases[[1]]$submitter_id,
               sample_type = h$cases[[1]]$samples[[1]]$sample_type,
               stringsAsFactors = FALSE)
  }))
}

## ---- bulk download a set of file_ids, parse STAR-Counts TSVs, return a
## gene x sample TPM matrix (gene_name rows, file_id columns) ----
download_and_parse <- function(file_ids, batch_size = 60) {
  if (length(file_ids) == 0) return(NULL)
  batches <- split(file_ids, ceiling(seq_along(file_ids) / batch_size))
  mats <- vector("list", length(batches))
  for (i in seq_along(batches)) {
    ids <- batches[[i]]
    tmp_tar <- tempfile(fileext = ".tar.gz")
    tmp_dir <- tempfile()
    dir.create(tmp_dir)
    r <- POST(GDC_DATA, body = toJSON(list(ids = as.list(ids)), auto_unbox = FALSE),
              content_type_json(), write_disk(tmp_tar, overwrite = TRUE))
    stop_for_status(r)
    untar(tmp_tar, exdir = tmp_dir)
    tsvs <- list.files(tmp_dir, pattern = "\\.tsv$", recursive = TRUE, full.names = TRUE)
    ## each tsv path is <file_id>/<uuid>.rna_seq...tsv -- recover file_id from parent dir
    batch_mat <- lapply(tsvs, function(f) {
      fid <- basename(dirname(f))
      df <- read.delim(f, skip = 1, header = TRUE, stringsAsFactors = FALSE)
      df <- df[!df$gene_id %in% c("N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous") &
                 !is.na(df$gene_name) & df$gene_name != "" & !duplicated(df$gene_name), ]
      setNames(list(setNames(df$tpm_unstranded, df$gene_name)), fid)
    })
    batch_mat <- do.call(c, batch_mat)
    all_genes <- Reduce(union, lapply(batch_mat, names))
    mats[[i]] <- sapply(batch_mat, function(x) x[all_genes])
    rownames(mats[[i]]) <- all_genes
    unlink(tmp_tar); unlink(tmp_dir, recursive = TRUE)
    cat(sprintf("  downloaded+parsed batch %d/%d (%d files)\n", i, length(batches), length(ids)))
  }
  all_genes <- Reduce(union, lapply(mats, rownames))
  combined <- sapply(mats, function(m) m[all_genes, , drop = FALSE])
  do.call(cbind, combined)
}

## ================= build sample manifest with entity labels =================
if (file.exists(expr_cache) && file.exists(meta_cache)) {

  cat("Cached TARGET bulk expression + metadata found -- skipping download.\n")
  expr <- readRDS(expr_cache)
  meta <- readRDS(meta_cache)

} else {

  cat("=== Querying GDC for candidate files ===\n")
  aml_tumor  <- list_files("TARGET-AML", c("Primary Blood Derived Cancer - Bone Marrow",
                                            "Primary Blood Derived Cancer - Peripheral Blood"))
  aml_normal <- list_files("TARGET-AML", c("Bone Marrow Normal", "Blood Derived Normal"))
  all2_tumor <- list_files("TARGET-ALL-P2", c("Primary Blood Derived Cancer - Bone Marrow",
                                               "Primary Blood Derived Cancer - Peripheral Blood"))
  all3_tumor <- list_files("TARGET-ALL-P3", c("Primary Blood Derived Cancer - Bone Marrow",
                                               "Primary Blood Derived Cancer - Peripheral Blood"))
  cat(sprintf("AML tumor: %d, AML normal: %d, ALL-P2 tumor: %d, ALL-P3 tumor: %d\n",
              nrow(aml_tumor), nrow(aml_normal), nrow(all2_tumor), nrow(all3_tumor)))

  cat("\n=== Downloading clinical label files ===\n")
  clin_dir <- file.path(dir_data, "clinical")
  if (!dir.exists(clin_dir)) dir.create(clin_dir, recursive = TRUE)
  fetch_clin <- function(file_id, name) {
    dest <- file.path(clin_dir, name)
    if (!file.exists(dest)) {
      r <- GET(paste0(GDC_DATA, "/", file_id), write_disk(dest, overwrite = TRUE))
      stop_for_status(r)
    }
    dest
  }
  ## file UUIDs identified and verified open-access during development (see script header)
  all2_clin_path <- fetch_clin("c6a63e40-f4e5-4002-9b87-3e6f3b3dfab4", "TARGET_ALL_P2_Validation.xlsx")
  all3_clin_path <- fetch_clin("d475576c-eac8-4f24-bfc7-67797b02b109", "TARGET_ALL_ClinicalData_Phase_III.xlsx")

  all2_clin <- read_excel(all2_clin_path) %>%
    transmute(case_id = `TARGET USI`, origin = `Cell of Origin`) %>%
    filter(origin %in% c("B Cell ALL", "T Cell ALL")) %>%
    mutate(entity = ifelse(origin == "B Cell ALL", "B-ALL", "T-ALL"))

  all3_clin <- read_excel(all3_clin_path) %>%
    transmute(case_id = `TARGET USI`, who_dx = trimws(`WHO Dx`)) %>%
    filter(who_dx %in% c("B/M", "T/M")) %>%
    mutate(entity = ifelse(who_dx == "B/M", "B-My-MPAL", "T-My-MPAL"))

  ## ---- join entity labels onto file manifests, cap large groups for tractability ----
  cap_n <- 150
  cap_sample <- function(df) if (nrow(df) > cap_n) slice_sample(df, n = cap_n) else df
  aml_lab  <- aml_tumor %>% mutate(entity = "AML") %>% distinct(case_id, .keep_all = TRUE) %>% cap_sample()
  norm_lab <- aml_normal %>% mutate(entity = "Normal") %>% distinct(case_id, .keep_all = TRUE) %>% cap_sample()
  ball_lab <- all2_tumor %>% inner_join(all2_clin %>% filter(entity == "B-ALL"), by = "case_id") %>%
    distinct(case_id, .keep_all = TRUE) %>% cap_sample()
  tall_lab <- all2_tumor %>% inner_join(all2_clin %>% filter(entity == "T-ALL"), by = "case_id") %>%
    distinct(case_id, .keep_all = TRUE) %>% cap_sample()
  bmpal_lab <- all3_tumor %>% inner_join(all3_clin %>% filter(entity == "B-My-MPAL"), by = "case_id") %>%
    distinct(case_id, .keep_all = TRUE)  ## small population -- keep all
  tmpal_lab <- all3_tumor %>% inner_join(all3_clin %>% filter(entity == "T-My-MPAL"), by = "case_id") %>%
    distinct(case_id, .keep_all = TRUE)

  manifest <- bind_rows(aml_lab, norm_lab, ball_lab, tall_lab, bmpal_lab, tmpal_lab) %>%
    select(file_id, case_id, entity)
  cat("\nManifest sample counts by entity:\n")
  print(table(manifest$entity))

  cat("\n=== Downloading + parsing gene expression files (this takes a while) ===\n")
  expr <- download_and_parse(manifest$file_id)
  meta <- manifest %>% filter(file_id %in% colnames(expr))

  saveRDS(expr, expr_cache)
  saveRDS(meta, meta_cache)
  cat(sprintf("\nCached expression matrix (%d genes x %d samples) and metadata to data/target_bulk/\n",
              nrow(expr), ncol(expr)))
}

## ================= score with strict up_lfc1 signature =================
cat("\n=== Scoring ===\n")
human_strict <- readRDS(file.path(proj_root, "data", "human_signature_sets_strict.rds"))
sig_genes <- human_strict$up_lfc1
n_present <- sum(sig_genes %in% rownames(expr))
cat(sprintf("Signature: up_lfc1, strict 1:1 human orthologs, %d genes (%d present in TARGET expression matrix, %.0f%%)\n",
            length(sig_genes), n_present, 100 * n_present / length(sig_genes)))

logexpr <- log2(expr + 1)
gene_var <- apply(logexpr, 1, var, na.rm = TRUE)
z <- (logexpr[gene_var > 0, ] - rowMeans(logexpr[gene_var > 0, ], na.rm = TRUE)) /
  apply(logexpr[gene_var > 0, ], 1, sd, na.rm = TRUE)

score_genes <- intersect(sig_genes, rownames(z))
patient_score <- colMeans(z[score_genes, , drop = FALSE], na.rm = TRUE)

df <- meta %>% mutate(score = patient_score[file_id]) %>% filter(!is.na(score))

## ---- sensitivity check: z-score computed WITHOUT Normal in the reference
## population -- tests whether the "every entity below Normal" pattern in the
## main (pooled) z-score is a genuine biological pattern or an artifact of
## z-scoring against a compositionally distinct Normal reference pooled in
## with heterogeneous-purity tumor samples (see RESULTS.md discussion). Same
## method as the main score, just restricted to leukemic samples only, matching
## 07_mulet_lazaro_alal_bulk.R's convention (leukemia-only z-score reference).
leuk_files <- meta$file_id[meta$entity != "Normal"]
logexpr_leuk <- logexpr[, leuk_files, drop = FALSE]
gene_var_leuk <- apply(logexpr_leuk, 1, var, na.rm = TRUE)
z_leuk <- (logexpr_leuk[gene_var_leuk > 0, ] - rowMeans(logexpr_leuk[gene_var_leuk > 0, ], na.rm = TRUE)) /
  apply(logexpr_leuk[gene_var_leuk > 0, ], 1, sd, na.rm = TRUE)
score_genes_leuk <- intersect(sig_genes, rownames(z_leuk))
patient_score_leuk <- colMeans(z_leuk[score_genes_leuk, , drop = FALSE], na.rm = TRUE)
df <- df %>% mutate(score_leuk_only = ifelse(file_id %in% names(patient_score_leuk),
                                              patient_score_leuk[file_id], NA_real_))

write.csv(df, file.path(dir_out, "target_entity_scores.csv"), row.names = FALSE)
cat(sprintf("\nSaved streamline/results/target_entity_scores.csv (%d patients)\n", nrow(df)))

cat("\nPatients per entity:\n")
print(table(df$entity))

## ================= figure =================
entity_order <- c("Normal", "AML", "B-ALL", "T-ALL", "B-My-MPAL", "T-My-MPAL")
entity_pal <- c(Normal = "grey60", AML = "#e6550d", `B-ALL` = "#3182bd", `T-ALL` = "#31a354",
                 `B-My-MPAL` = "#c51b8a", `T-My-MPAL` = "#756bb1")
df_plot <- df %>% filter(entity %in% entity_order) %>% mutate(entity = factor(entity, levels = entity_order))

p <- ggplot(df_plot, aes(x = entity, y = score, fill = entity)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.85, width = 0.6) +
  geom_jitter(width = 0.15, size = 1, alpha = 0.5) +
  scale_fill_manual(values = entity_pal, guide = "none") +
  labs(title = "TARGET bulk RNA-seq: LT-PLC signature score across leukemia entities",
       subtitle = sprintf("mean gene-wise z-score (log2 TPM+1), up_lfc1 signature (strict 1:1 orthologs, %d genes, %d present); one point = one patient",
                           length(sig_genes), n_present),
       x = NULL, y = "Mean z-score (LT-PLC up_lfc1)") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1),
        plot.title = element_text(face = "bold", size = 12.5),
        plot.subtitle = element_text(size = 8.5, color = "grey30"),
        panel.grid.minor = element_blank())

ggsave(file.path(dir_out, "target_entity_comparison.png"), p, width = 9, height = 6, dpi = 300)
ggsave(file.path(dir_out, "target_entity_comparison.pdf"), p, width = 9, height = 6)
cat("Saved streamline/results/target_entity_comparison.png / .pdf\n")

## ================= statistics =================
cat("\nMedian score by entity:\n")
print(df_plot %>% group_by(entity) %>% summarise(n = n(), median_score = median(score)) %>% arrange(desc(median_score)))

cat("\nKruskal-Wallis across all entities:\n")
print(kruskal.test(score ~ entity, data = df_plot))

cat("\nEach entity vs Normal, Wilcoxon:\n")
norm_scores <- df_plot$score[df_plot$entity == "Normal"]
for (e in setdiff(entity_order, "Normal")) {
  a <- df_plot$score[df_plot$entity == e]
  if (length(a) >= 2 && length(norm_scores) >= 2) {
    wt <- wilcox.test(a, norm_scores)
    cat(sprintf("%-10s (n=%d) vs Normal (n=%d): median %.4f vs %.4f, Wilcoxon p=%.4g\n",
                e, length(a), length(norm_scores), median(a), median(norm_scores), wt$p.value))
  }
}

cat("\nB-ALL vs T-ALL, Wilcoxon:\n")
print(wilcox.test(df_plot$score[df_plot$entity == "B-ALL"], df_plot$score[df_plot$entity == "T-ALL"]))

cat("\nB-My-MPAL vs T-My-MPAL, Wilcoxon:\n")
print(wilcox.test(df_plot$score[df_plot$entity == "B-My-MPAL"], df_plot$score[df_plot$entity == "T-My-MPAL"]))

cat("\nT-My-MPAL vs T-ALL, Wilcoxon (both T-lineage):\n")
print(wilcox.test(df_plot$score[df_plot$entity == "T-My-MPAL"], df_plot$score[df_plot$entity == "T-ALL"]))

cat("\nT-My-MPAL vs AML, Wilcoxon:\n")
print(wilcox.test(df_plot$score[df_plot$entity == "T-My-MPAL"], df_plot$score[df_plot$entity == "AML"]))

## ================= sensitivity: leukemia-only z-score reference =================
cat("\n\n=== SENSITIVITY CHECK: z-score computed WITHOUT Normal in the reference ===\n")
leuk_entities <- c("AML", "B-ALL", "T-ALL", "B-My-MPAL", "T-My-MPAL")
df_leuk <- df %>% filter(entity %in% leuk_entities, !is.na(score_leuk_only)) %>%
  mutate(entity = factor(entity, levels = leuk_entities))

cat("\nMedian leukemia-only-reference score by entity:\n")
print(df_leuk %>% group_by(entity) %>% summarise(n = n(), median_score = median(score_leuk_only)) %>% arrange(desc(median_score)))

p_leuk <- ggplot(df_leuk, aes(x = entity, y = score_leuk_only, fill = entity)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.85, width = 0.6) +
  geom_jitter(width = 0.15, size = 1, alpha = 0.5) +
  scale_fill_manual(values = entity_pal[leuk_entities], guide = "none") +
  labs(title = "Sensitivity check: LT-PLC score with z-score reference restricted to leukemia entities only",
       subtitle = "same up_lfc1 signature and method as the main figure, but Normal excluded from the z-score reference population\n(matches 07_mulet_lazaro_alal_bulk.R's convention) -- tests whether entity ranking depends on including Normal in scoring",
       x = NULL, y = "Mean z-score (LT-PLC up_lfc1), leukemia-only reference") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1),
        plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 8, color = "grey30"),
        panel.grid.minor = element_blank())

ggsave(file.path(dir_out, "target_entity_comparison_leukonly.png"), p_leuk, width = 9, height = 6, dpi = 300)
ggsave(file.path(dir_out, "target_entity_comparison_leukonly.pdf"), p_leuk, width = 9, height = 6)
cat("Saved streamline/results/target_entity_comparison_leukonly.png / .pdf\n")

cat("\nKruskal-Wallis across leukemia entities (leukemia-only reference):\n")
print(kruskal.test(score_leuk_only ~ entity, data = df_leuk))

cat("\nB-ALL vs T-ALL (leukemia-only reference):\n")
print(wilcox.test(df_leuk$score_leuk_only[df_leuk$entity == "B-ALL"], df_leuk$score_leuk_only[df_leuk$entity == "T-ALL"]))

cat("\nT-My-MPAL vs T-ALL (leukemia-only reference):\n")
print(wilcox.test(df_leuk$score_leuk_only[df_leuk$entity == "T-My-MPAL"], df_leuk$score_leuk_only[df_leuk$entity == "T-ALL"]))

cat("\nT-My-MPAL vs AML (leukemia-only reference):\n")
print(wilcox.test(df_leuk$score_leuk_only[df_leuk$entity == "T-My-MPAL"], df_leuk$score_leuk_only[df_leuk$entity == "AML"]))

cat("\nT-My-MPAL vs B-My-MPAL (leukemia-only reference):\n")
print(wilcox.test(df_leuk$score_leuk_only[df_leuk$entity == "T-My-MPAL"], df_leuk$score_leuk_only[df_leuk$entity == "B-My-MPAL"]))

cat("\nRank correlation between main (pooled) score and leukemia-only-reference score (leukemic patients only):\n")
print(cor.test(df_leuk$score[match(df_leuk$file_id, df$file_id)], df_leuk$score_leuk_only, method = "spearman"))

cat("\nDone.\n")
