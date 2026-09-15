## 06_target_bulk_full.R
## Same analysis as 03_target_bulk_entity_comparison.R (Fig. 6 fix), but using
## the ENTIRE available TARGET cohort per entity instead of the n=150 cap used
## there for tractability -- all AML tumor + Normal samples, all TARGET-ALL-P2
## samples matching a B-ALL/T-ALL clinical label, plus the already-uncapped
## MPAL subtypes (small populations in TARGET regardless). This is the
## maximum-power version of the same comparison; 03's capped version stays in
## place as the faster/reproducible baseline.
##
## Cached separately from 03 (target_bulk_expr_full.rds / _meta_full.rds,
## target_entity_*_full.csv/.png) so neither run overwrites the other.
##
## Input:  GDC API (same open-access STAR-Counts + clinical Excel files as 03)
##         data/human_signature_sets_strict.rds (1:1 orthologs; up_lfc1)
## Output: data/target_bulk/target_bulk_expr_full.rds, target_bulk_meta_full.rds (cache)
##         streamline/results/target_entity_scores_full.csv
##         streamline/results/target_entity_comparison_full.png / .pdf
##         streamline/results/target_entity_comparison_full_leukonly.png / .pdf

suppressPackageStartupMessages({
  library(httr); library(jsonlite); library(readxl); library(dplyr); library(ggplot2)
})

set.seed(1234)
proj_root <- getwd()
dir_data  <- file.path(proj_root, "data", "target_bulk")
dir_out   <- file.path(proj_root, "streamline", "results")
for (d in c(dir_data, dir_out)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

expr_cache <- file.path(dir_data, "target_bulk_expr_full.rds")
meta_cache <- file.path(dir_data, "target_bulk_meta_full.rds")

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

## batch_size=100 (vs 60 in 03) -- fewer, larger batches to cut per-batch HTTP/untar
## overhead given the much larger total file count here.
download_and_parse <- function(file_ids, batch_size = 100) {
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
    gc(verbose = FALSE, full = TRUE)
  }
  all_genes <- Reduce(union, lapply(mats, rownames))
  combined <- sapply(mats, function(m) m[all_genes, , drop = FALSE])
  do.call(cbind, combined)
}

## ================= build sample manifest with entity labels (NO CAPS) =================
if (file.exists(expr_cache) && file.exists(meta_cache)) {

  cat("Cached full-TARGET bulk expression + metadata found -- skipping download.\n")
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

  ## ---- join entity labels onto file manifests -- NO cap this time, use everything ----
  aml_lab  <- aml_tumor %>% mutate(entity = "AML") %>% distinct(case_id, .keep_all = TRUE)
  norm_lab <- aml_normal %>% mutate(entity = "Normal") %>% distinct(case_id, .keep_all = TRUE)
  ball_lab <- all2_tumor %>% inner_join(all2_clin %>% filter(entity == "B-ALL"), by = "case_id") %>%
    distinct(case_id, .keep_all = TRUE)
  tall_lab <- all2_tumor %>% inner_join(all2_clin %>% filter(entity == "T-ALL"), by = "case_id") %>%
    distinct(case_id, .keep_all = TRUE)
  bmpal_lab <- all3_tumor %>% inner_join(all3_clin %>% filter(entity == "B-My-MPAL"), by = "case_id") %>%
    distinct(case_id, .keep_all = TRUE)
  tmpal_lab <- all3_tumor %>% inner_join(all3_clin %>% filter(entity == "T-My-MPAL"), by = "case_id") %>%
    distinct(case_id, .keep_all = TRUE)

  manifest <- bind_rows(aml_lab, norm_lab, ball_lab, tall_lab, bmpal_lab, tmpal_lab) %>%
    select(file_id, case_id, entity)
  cat("\nManifest sample counts by entity (full TARGET set, no caps):\n")
  print(table(manifest$entity))

  cat(sprintf("\n=== Downloading + parsing %d gene expression files (this will take a while) ===\n", nrow(manifest)))
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

## leukemia-only z-score reference (Normal excluded) -- the version 03's
## sensitivity check determined is the one to trust/report as primary.
leuk_files <- meta$file_id[meta$entity != "Normal"]
logexpr_leuk <- logexpr[, leuk_files, drop = FALSE]
gene_var_leuk <- apply(logexpr_leuk, 1, var, na.rm = TRUE)
z_leuk <- (logexpr_leuk[gene_var_leuk > 0, ] - rowMeans(logexpr_leuk[gene_var_leuk > 0, ], na.rm = TRUE)) /
  apply(logexpr_leuk[gene_var_leuk > 0, ], 1, sd, na.rm = TRUE)
score_genes_leuk <- intersect(sig_genes, rownames(z_leuk))
patient_score_leuk <- colMeans(z_leuk[score_genes_leuk, , drop = FALSE], na.rm = TRUE)

df <- meta %>% mutate(score = patient_score_leuk[file_id]) %>% filter(!is.na(score), entity != "Normal")
write.csv(df, file.path(dir_out, "target_entity_scores_full.csv"), row.names = FALSE)
cat(sprintf("\nSaved streamline/results/target_entity_scores_full.csv (%d leukemia patients)\n", nrow(df)))

cat("\nPatients per entity (leukemia-only, full TARGET set):\n")
print(table(df$entity))

## ================= figure =================
entity_order <- c("AML", "B-ALL", "T-My-MPAL", "B-My-MPAL", "T-ALL")
entity_pal <- c(AML = "#e6550d", `B-ALL` = "#3182bd", `T-My-MPAL` = "#756bb1",
                 `B-My-MPAL` = "#c51b8a", `T-ALL` = "#31a354")
df_plot <- df %>% mutate(entity = factor(entity, levels = entity_order))

cat("\nMedian score by entity (full TARGET set, leukemia-only reference):\n")
print(df_plot %>% group_by(entity) %>% summarise(n = n(), median_score = median(score)) %>% arrange(desc(median_score)))

p <- ggplot(df_plot, aes(x = entity, y = score, fill = entity)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.85, width = 0.6) +
  geom_jitter(width = 0.15, size = 0.6, alpha = 0.35) +
  scale_fill_manual(values = entity_pal, guide = "none") +
  labs(title = "TARGET bulk RNA-seq (full cohort): LT-PLC score across leukemia entities",
       subtitle = sprintf("mean gene-wise z-score (log2 TPM+1), up_lfc1 signature (strict 1:1 orthologs, %d genes, %d present);\nleukemia-only z-score reference; entire available TARGET-AML/ALL-P2/ALL-P3 cohort, no per-entity cap; one point = one patient",
                           length(sig_genes), n_present),
       x = NULL, y = "Mean z-score (LT-PLC up_lfc1), leukemia-only reference") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1),
        plot.title = element_text(face = "bold", size = 12.5),
        plot.subtitle = element_text(size = 8.5, color = "grey30"),
        panel.grid.minor = element_blank())

ggsave(file.path(dir_out, "target_entity_comparison_full.png"), p, width = 9, height = 6, dpi = 300)
ggsave(file.path(dir_out, "target_entity_comparison_full.pdf"), p, width = 9, height = 6)
cat("Saved streamline/results/target_entity_comparison_full.png / .pdf\n")

## ================= statistics =================
cat("\nKruskal-Wallis across all entities:\n")
print(kruskal.test(score ~ entity, data = df_plot))

cat("\nT-My-MPAL vs T-ALL, Wilcoxon (both T-lineage):\n")
print(wilcox.test(df_plot$score[df_plot$entity == "T-My-MPAL"], df_plot$score[df_plot$entity == "T-ALL"]))

cat("\nT-My-MPAL vs AML, Wilcoxon:\n")
print(wilcox.test(df_plot$score[df_plot$entity == "T-My-MPAL"], df_plot$score[df_plot$entity == "AML"]))

cat("\nT-My-MPAL vs B-My-MPAL, Wilcoxon:\n")
print(wilcox.test(df_plot$score[df_plot$entity == "T-My-MPAL"], df_plot$score[df_plot$entity == "B-My-MPAL"]))

cat("\nB-ALL vs T-ALL, Wilcoxon:\n")
print(wilcox.test(df_plot$score[df_plot$entity == "B-ALL"], df_plot$score[df_plot$entity == "T-ALL"]))

cat("\nAML vs T-ALL, Wilcoxon:\n")
print(wilcox.test(df_plot$score[df_plot$entity == "AML"], df_plot$score[df_plot$entity == "T-ALL"]))

## ---- compare with the capped (n<=150) result from 03, where entities overlap ----
capped_path <- file.path(dir_out, "target_entity_scores.csv")
if (file.exists(capped_path)) {
  capped <- read.csv(capped_path) %>% filter(entity != "Normal") %>%
    group_by(entity) %>% summarise(n_capped = n(), median_capped = median(score_leuk_only), .groups = "drop")
  full_summary <- df_plot %>% group_by(entity) %>% summarise(n_full = n(), median_full = median(score), .groups = "drop")
  cat("\nCapped (n<=150) vs full-cohort median score, by entity:\n")
  print(full_summary %>% left_join(capped, by = "entity"))
}

cat("\nDone.\n")
