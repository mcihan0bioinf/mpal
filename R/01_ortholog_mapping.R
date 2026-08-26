## 01_ortholog_mapping.R
## Build the LT-PLC mouse signature gene set variants from the DESeq2
## diff-expression table (LT-PLC vs non-LT-PLC, LSK primary leukemia, n=5),
## then map each set to human orthologs in two modes:
##   - strict:  1:1 orthologs only (conservative, unambiguous)
##   - maximal: every valid ortholog hit kept, incl. one-to-many/paralogs
##              (maximizes gene coverage; downstream scoring uses this
##              unless a gene set specifically needs the strict version)
##
## Input:  Diff_genes_annotated.csv (project root; not tracked in git)
##         columns: mgi_symbol, log2FoldChange, padj (+ others)
## Output: data/mouse_signature_sets.rds
##         data/human_signature_sets.rds         (maximal mapping; use this downstream)
##         data/human_signature_sets_strict.rds  (strict 1:1 mapping)
##         results/ortholog_mapping_summary.csv

set.seed(1234)
suppressPackageStartupMessages({
  library(dplyr)
  library(gprofiler2)
})

proj_root <- getwd()
stopifnot("Diff_genes_annotated.csv" %in% list.files(proj_root))
dir_data <- file.path(proj_root, "data")
dir_out  <- file.path(proj_root, "results")
for (d in c(dir_data, dir_out)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

## ---- 1. Load DEG table and build mouse gene set variants ----

deg <- read.csv2(file.path(proj_root, "Diff_genes_annotated.csv"),
                  row.names = 1, stringsAsFactors = FALSE)
# read.csv2: ';' sep, ',' decimal (matches the file's EU-locale export)

stopifnot(all(c("log2FoldChange", "padj", "mgi_symbol") %in% names(deg)))

deg <- deg %>%
  filter(!is.na(mgi_symbol), mgi_symbol != "", !is.na(padj), !is.na(log2FoldChange)) %>%
  distinct(mgi_symbol, .keep_all = TRUE)

cat(sprintf("Loaded Diff_genes_annotated.csv: %d genes with usable mgi_symbol/padj/log2FC\n", nrow(deg)))

sig_sets <- list(
  up_all    = deg %>% filter(log2FoldChange > 0, padj < 0.05) %>% pull(mgi_symbol),
  up_lfc1   = deg %>% filter(log2FoldChange >= 1, padj < 0.05) %>% pull(mgi_symbol),
  up_lfc2   = deg %>% filter(log2FoldChange >= 2, padj < 0.05) %>% pull(mgi_symbol),
  up_top50  = deg %>% filter(padj < 0.05) %>% arrange(desc(log2FoldChange)) %>% slice_head(n = 50)  %>% pull(mgi_symbol),
  up_top100 = deg %>% filter(padj < 0.05) %>% arrange(desc(log2FoldChange)) %>% slice_head(n = 100) %>% pull(mgi_symbol),
  down_lfc1 = deg %>% filter(log2FoldChange <= -1, padj < 0.05) %>% pull(mgi_symbol)
)
for (nm in names(sig_sets)) cat(sprintf("  %-10s n = %d\n", nm, length(sig_sets[[nm]])))

saveRDS(sig_sets, file.path(dir_data, "mouse_signature_sets.rds"))
saveRDS(deg,       file.path(dir_data, "deg_table_clean.rds"))

## ---- 2. Map each set to human orthologs via gprofiler2::gorth ----
## Two modes computed from the same gorth() call:
##   strict  = keep only mouse genes with exactly one human ortholog, and
##             only human genes claimed by exactly one mouse gene
##   maximal = keep every valid ortholog hit (all paralogs / one-to-many
##             mappings included); maximizes the number of genes retained

clean_ortholog_names <- function(x) {
  x <- unique(x)
  x[x != "N/A" & !is.na(x) & x != ""]
}

map_orthologs <- function(genes, set_name) {
  if (length(genes) == 0) return(list(strict = character(0), maximal = character(0)))
  orth <- gorth(query = genes, source_organism = "mmusculus", target_organism = "hsapiens",
                mthreshold = Inf, filter_na = TRUE)

  maximal <- clean_ortholog_names(orth$ortholog_name)

  n_per_mouse <- table(orth$input)
  one_to_one_mouse <- names(n_per_mouse[n_per_mouse == 1])
  orth1 <- orth %>% filter(input %in% one_to_one_mouse)
  n_per_human <- table(orth1$ortholog_name)
  ambiguous_human <- names(n_per_human[n_per_human > 1])
  orth1 <- orth1 %>% filter(!ortholog_name %in% ambiguous_human)
  strict <- clean_ortholog_names(orth1$ortholog_name)

  cat(sprintf("  %-10s: %d mouse genes -> %d human orthologs (strict 1:1) / %d human orthologs (maximal)\n",
              set_name, length(genes), length(strict), length(maximal)))
  list(strict = strict, maximal = maximal)
}

mapped <- lapply(names(sig_sets), function(nm) map_orthologs(sig_sets[[nm]], nm))
names(mapped) <- names(sig_sets)

human_sets_strict  <- lapply(mapped, `[[`, "strict")
human_sets_maximal <- lapply(mapped, `[[`, "maximal")

saveRDS(human_sets_maximal, file.path(dir_data, "human_signature_sets.rds"))
saveRDS(human_sets_strict,  file.path(dir_data, "human_signature_sets_strict.rds"))

summary_df <- data.frame(set = names(sig_sets),
                          n_mouse = sapply(sig_sets, length),
                          n_human_strict  = sapply(human_sets_strict, length),
                          n_human_maximal = sapply(human_sets_maximal, length))
write.csv(summary_df, file.path(dir_out, "ortholog_mapping_summary.csv"), row.names = FALSE)
print(summary_df)
