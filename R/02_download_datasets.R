## 02_download_datasets.R
## Download and prepare the four labeled single-cell leukemia atlases used
## in this project. Idempotent: skips downloading a file that already
## exists. Each section downloads the raw published files (author-provided
## cell type / disease labels, no manual relabeling) and builds a
## harmonized Seurat object with common metadata columns:
##   celltype  -- author cell type / classification label
##   group     -- harmonized disease group (e.g. "MPAL", "AML", "B-ALL",
##                "T-ALL", "Healthy")
##   sample    -- sample / donor identifier
## Granja additionally gets celltype_fine (see section 1), combining the
## MPAL disease classification with a finer healthy-reference cell typing.
## saved to data/<atlas>_seurat.rds (not tracked in git).
##
## Atlases:
##   granja        -- Granja et al. 2019 Nat Biotechnol (MPAL + healthy ref)
##   bonemarrowmap -- Zeng et al. 2025 Blood Cancer Discov (healthy BM ref)
##   caron         -- Caron et al. 2020 Sci Rep (B-ALL, T-ALL, healthy)
##   vangalen      -- van Galen et al. 2019 Cell (AML, healthy, cell lines)

set.seed(1234)
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(data.table)
  library(SummarizedExperiment)
})

proj_root <- getwd()
dir_data  <- file.path(proj_root, "data")
if (!dir.exists(dir_data)) dir.create(dir_data, recursive = TRUE)

download_if_missing <- function(url, dest) {
  if (!file.exists(dest)) {
    cat("Downloading", basename(dest), "...\n")
    dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
    download.file(url, dest, mode = "wb", quiet = FALSE)
  } else {
    cat("Already present:", basename(dest), "\n")
  }
}

## ================= 1. Granja et al. 2019 (MPAL) =================
cat("\n=== Granja 2019 ===\n")
granja_dir <- file.path(dir_data, "granja_2019")
granja_file <- file.path(granja_dir, "scRNA-All-Hematopoiesis-MPAL-191120.rds")
download_if_missing(
  "https://jeffgranja.s3.amazonaws.com/MPAL-10x/Supplementary_Data/Healthy-Disease-Data/scRNA-All-Hematopoiesis-MPAL-191120.rds",
  granja_file)

se <- readRDS(granja_file)
counts <- assay(se, "counts")
rownames(counts) <- make.unique(rowData(se)$gene_name)
cd <- as.data.frame(colData(se))
seu <- CreateSeuratObject(counts = counts, meta.data = cd, project = "granja_2019")
seu <- NormalizeData(seu, verbose = FALSE)
umap_mat <- as.matrix(cd[, c("ProjectUMAP1", "ProjectUMAP2")])
colnames(umap_mat) <- c("UMAP_1", "UMAP_2"); rownames(umap_mat) <- colnames(seu)
seu[["umap"]] <- CreateDimReducObject(embeddings = umap_mat, key = "UMAP_", assay = DefaultAssay(seu))
seu$celltype <- seu$ProjectClassification
seu$group <- ifelse(grepl("^MPAL", seu$Group), "MPAL", "Healthy")
seu$sample <- seu$Group

## seu$celltype above only labels the malignant MPAL blasts (6 categories);
## healthy comparator cells are just "Reference" in that column. Pull the
## author's finer BioClassification (26 normal cell states) for those cells
## from the companion healthy-only object and merge in by Group+Barcode, so
## every cell ends up with a label in exactly one of two complementary
## annotation layers: celltype (MPAL disease classification) or
## celltype_fine (MPAL classification for blasts, BioClassification for
## healthy cells) -- see 03_atlas_celltypes_granja.R.
healthy_file <- file.path(granja_dir, "scRNA-Healthy-Hematopoiesis-191120.rds")
download_if_missing(
  "https://jeffgranja.s3.amazonaws.com/MPAL-10x/Supplementary_Data/Healthy-Data/scRNA-Healthy-Hematopoiesis-191120.rds",
  healthy_file)
se_h <- readRDS(healthy_file)
cd_h <- as.data.frame(colData(se_h))
cd_h$key <- paste(cd_h$Group, cd_h$Barcode, sep = "_")
key_m <- paste(seu$Group, seu$Barcode, sep = "_")
bioclass <- cd_h$BioClassification[match(key_m, cd_h$key)]
seu$celltype_fine <- ifelse(seu$celltype == "Reference" & !is.na(bioclass),
                             sub("^[0-9]+_", "", bioclass),  # "01_HSC" -> "HSC"
                             as.character(seu$celltype))
cat(sprintf("Granja: merged BioClassification for %d/%d healthy reference cells\n",
            sum(seu$celltype == "Reference" & !is.na(bioclass)), sum(seu$celltype == "Reference")))
rm(se_h, cd_h, bioclass, key_m)

saveRDS(seu, file.path(dir_data, "granja_seurat.rds"))
cat(sprintf("Granja: %d cells x %d genes -> data/granja_seurat.rds\n", ncol(seu), nrow(seu)))
rm(se, counts, cd, seu); gc()

## ================= 2. BoneMarrowMap (healthy reference) =================
cat("\n=== BoneMarrowMap ===\n")
bmm_dir <- file.path(dir_data, "bonemarrowmap")
bmm_file <- file.path(bmm_dir, "BoneMarrowMap_Annotated_Dataset_expandedFeatures.rds")
download_if_missing(
  "https://bonemarrowmap.s3.us-east-2.amazonaws.com/BoneMarrowMap_Annotated_Dataset_expandedFeatures.rds",
  bmm_file)

seu <- readRDS(bmm_file)
DefaultAssay(seu) <- "RNA"
seu[["SCENIC_TFs"]] <- NULL
seu$celltype <- seu$CellType
seu$group    <- "Healthy"
seu$sample   <- seu$Donor
saveRDS(seu, file.path(dir_data, "bonemarrowmap_seurat.rds"))
cat(sprintf("BoneMarrowMap: %d cells x %d genes -> data/bonemarrowmap_seurat.rds\n", ncol(seu), nrow(seu)))
rm(seu); gc()

## ================= 3. Caron et al. 2020 (B-ALL, T-ALL) =================
cat("\n=== Caron 2020 ===\n")
caron_dir <- file.path(dir_data, "caron_2020")
download_if_missing("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE132nnn/GSE132509/suppl/GSE132509_RAW.tar",
                     file.path(caron_dir, "GSE132509_RAW.tar"))
download_if_missing("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE132nnn/GSE132509/suppl/GSE132509_cell_annotations.tsv.gz",
                     file.path(caron_dir, "GSE132509_cell_annotations.tsv.gz"))
raw_dir <- file.path(caron_dir, "raw")
if (!dir.exists(raw_dir) || length(list.files(raw_dir)) == 0) {
  dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
  untar(file.path(caron_dir, "GSE132509_RAW.tar"), exdir = raw_dir)
}

samples <- unique(sub("\\.matrix\\.mtx\\.gz$", "", list.files(raw_dir, pattern = "\\.matrix\\.mtx\\.gz$")))
read_one <- function(prefix) {
  mtx <- readMM(file.path(raw_dir, paste0(prefix, ".matrix.mtx.gz")))
  genes <- read.delim(file.path(raw_dir, paste0(prefix, ".genes.tsv.gz")), header = FALSE, col.names = c("ensembl", "symbol"))
  barcodes <- readLines(file.path(raw_dir, paste0(prefix, ".barcodes.tsv.gz")))
  rownames(mtx) <- make.unique(genes$symbol); colnames(mtx) <- barcodes
  mtx
}
mats <- lapply(samples, read_one); names(mats) <- samples

ann <- read.delim(file.path(caron_dir, "GSE132509_cell_annotations.tsv.gz"), stringsAsFactors = FALSE)
ann <- ann %>% filter(cell_id != "", orig.ident != "orig.ident")
sample_tags <- sort(unique(ann$orig.ident))

## dem/anno sample tags share the same disease-prefix + index ordering as
## the GSM submission order (verified: identical counts per prefix)
gsm_prefix <- gsub("^GSM[0-9]+_|_[0-9]+$", "", samples)
tag_prefix <- gsub("\\.[0-9]+$", "", sample_tags)
gsm_prefix_norm <- gsub("-", ".", gsm_prefix)
tag_prefix_norm <- gsub("-", ".", tag_prefix)
t1 <- table(gsm_prefix_norm); names(dimnames(t1)) <- NULL
t2 <- table(tag_prefix_norm); names(dimnames(t2)) <- NULL
stopifnot(identical(sort(t1), sort(t2)))

map_df <- data.frame(gsm = samples, gsm_prefix = gsm_prefix_norm, stringsAsFactors = FALSE) %>%
  group_by(gsm_prefix) %>% mutate(idx = row_number()) %>% ungroup()
tag_df <- data.frame(tag = sample_tags, tag_prefix = tag_prefix_norm, stringsAsFactors = FALSE) %>%
  group_by(tag_prefix) %>% mutate(idx = row_number()) %>% ungroup()
map_df <- map_df %>% left_join(tag_df, by = c("gsm_prefix" = "tag_prefix", "idx"))

all_genes <- Reduce(union, lapply(mats, rownames))
combined_list <- list()
for (i in seq_along(mats)) {
  m <- mats[[i]]
  tag <- map_df$tag[map_df$gsm == names(mats)[i]]
  colnames(m) <- paste0(tag, "_", sub("-1$", "", colnames(m)))
  missing <- setdiff(all_genes, rownames(m))
  if (length(missing) > 0) {
    pad <- Matrix(0, nrow = length(missing), ncol = ncol(m), sparse = TRUE, dimnames = list(missing, colnames(m)))
    m <- rbind(m, pad)
  }
  combined_list[[i]] <- m[all_genes, , drop = FALSE]
}
counts <- do.call(cbind, combined_list)
rm(mats, combined_list); gc()

rownames(ann) <- ann$cell_id
common_cells <- intersect(colnames(counts), rownames(ann))
counts <- counts[, common_cells]; ann <- ann[common_cells, ]

seu <- CreateSeuratObject(counts = counts, meta.data = ann, project = "caron_2020")
seu <- NormalizeData(seu, verbose = FALSE)
umap_mat <- as.matrix(ann[, c("UMAP1", "UMAP2")]); colnames(umap_mat) <- c("UMAP_1", "UMAP_2"); rownames(umap_mat) <- colnames(seu)
seu[["umap"]] <- CreateDimReducObject(embeddings = umap_mat, key = "UMAP_", assay = DefaultAssay(seu))
# ann$celltype (author-provided) already carried through via meta.data as
# seu$celltype: literal cell types for healthy PBMMC samples, sample ID for
# leukemic samples (blasts not further sub-typed in the original annotation)
seu$group <- case_when(
  grepl("^ETV6|^HHD", seu$orig.ident) ~ "B-ALL",
  grepl("^PRE.T", seu$orig.ident) ~ "T-ALL",
  grepl("^PBMMC", seu$orig.ident) ~ "Healthy",
  TRUE ~ NA_character_
)
seu$sample <- seu$orig.ident
saveRDS(seu, file.path(dir_data, "caron_seurat.rds"))
cat(sprintf("Caron: %d cells x %d genes -> data/caron_seurat.rds\n", ncol(seu), nrow(seu)))
rm(counts, ann, seu); gc()

## ================= 4. van Galen et al. 2019 (AML) =================
cat("\n=== van Galen 2019 ===\n")
vg_dir <- file.path(dir_data, "vangalen_2019")
download_if_missing("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE116nnn/GSE116256/suppl/GSE116256_RAW.tar",
                     file.path(vg_dir, "GSE116256_RAW.tar"))
vg_raw <- file.path(vg_dir, "raw")
if (!dir.exists(vg_raw) || length(list.files(vg_raw)) == 0) {
  dir.create(vg_raw, recursive = TRUE, showWarnings = FALSE)
  untar(file.path(vg_dir, "GSE116256_RAW.tar"), exdir = vg_raw)
}

## NOTE: do not use the pre-packaged Figshare "Seurat_AML.rds" object for
## this atlas -- it stores counts as a DENSE matrix and needs ~15GB just to
## deserialize. The raw per-sample GEO files below are small and sparse.
dem_files <- list.files(vg_raw, pattern = "\\.dem\\.txt\\.gz$", full.names = TRUE)
anno_files <- list.files(vg_raw, pattern = "\\.anno\\.txt\\.gz$", full.names = TRUE)
anno_sample_id <- sub("^GSM[0-9]+_", "", sub("\\.anno\\.txt\\.gz$", "", basename(anno_files)))
names(anno_files) <- anno_sample_id

read_sample <- function(dem_path) {
  sample_id <- sub("^GSM[0-9]+_", "", sub("\\.dem\\.txt\\.gz$", "", basename(dem_path)))
  dt <- fread(dem_path, sep = "\t", header = TRUE, data.table = TRUE)
  genes <- dt[[1]]; mat <- as.matrix(dt[, -1, with = FALSE]); rownames(mat) <- genes
  mat_sparse <- Matrix(mat, sparse = TRUE)
  rm(dt, mat); gc(verbose = FALSE)
  list(sample_id = sample_id, mat = mat_sparse)
}
read_anno <- function(p) { a <- fread(p, sep = "\t", header = TRUE, data.table = FALSE); rownames(a) <- a$Cell; a }

all_mats <- list(); all_anno <- list()
for (f in dem_files) {
  s <- read_sample(f)
  all_mats[[s$sample_id]] <- s$mat
  ap <- anno_files[[s$sample_id]]
  if (!is.null(ap) && file.exists(ap)) all_anno[[s$sample_id]] <- read_anno(ap)
}
all_genes <- Reduce(union, lapply(all_mats, rownames))
padded <- lapply(all_mats, function(m) {
  missing <- setdiff(all_genes, rownames(m))
  if (length(missing) > 0) {
    pad <- Matrix(0, nrow = length(missing), ncol = ncol(m), sparse = TRUE, dimnames = list(missing, colnames(m)))
    m <- rbind(m, pad)
  }
  m[all_genes, , drop = FALSE]
})
counts <- do.call(cbind, padded)
rm(all_mats, padded); gc()

ann_all <- bind_rows(all_anno); rownames(ann_all) <- ann_all$Cell
common_cells <- intersect(colnames(counts), rownames(ann_all))
counts <- counts[, common_cells]; ann_all <- ann_all[common_cells, ]

seu <- CreateSeuratObject(counts = counts, meta.data = ann_all, project = "vangalen_2019")
rm(counts, ann_all); gc()
seu$sample <- sub("_[ACGT]+$", "", colnames(seu))
seu$group <- case_when(
  grepl("^AML", seu$sample) ~ "AML",
  grepl("^BM", seu$sample) ~ "Healthy",
  grepl("^MUTZ3|^OCI", seu$sample) ~ "CellLine",
  TRUE ~ NA_character_
)
seu$celltype <- seu$CellType

## no existing embedding in the raw GEO files -- compute one (HVG-restricted
## scale/PCA to bound memory; scaling all genes OOM-killed this on 15GB RAM)
seu <- NormalizeData(seu, verbose = FALSE)
seu <- FindVariableFeatures(seu, verbose = FALSE)
seu <- ScaleData(seu, features = VariableFeatures(seu), verbose = FALSE)
seu <- RunPCA(seu, features = VariableFeatures(seu), npcs = 30, verbose = FALSE)
seu <- RunUMAP(seu, dims = 1:30, verbose = FALSE, seed.use = 1234)

saveRDS(seu, file.path(dir_data, "vangalen_seurat.rds"))
cat(sprintf("van Galen: %d cells x %d genes -> data/vangalen_seurat.rds\n", ncol(seu), nrow(seu)))
rm(seu); gc()

cat("\nAll atlases downloaded and prepared.\n")
