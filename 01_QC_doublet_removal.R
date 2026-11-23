# ==============================================================================
# 01_QC_doublet_removal.R
#
# Quality control and per-sample doublet removal (DoubletFinder) for four public
# human lung scRNA-seq datasets: GSE168191, GSE178360, GSE196341, GSE196638.
#
# Input : data/<GSE>_noQC.rds
#         Unfiltered Seurat objects, one per GEO series, with metadata columns
#         `sample_id` (GSM accession) and `disease`. Each object was created by
#         reading the authors' original count matrices deposited in GEO into
#         Seurat (CreateSeuratObject) and merging all samples of that series,
#         without any filtering or QC.
# Output: results/01_QC/singlets_COPD_Healthy.rds
#
# Requires Seurat v5 and DoubletFinder >= 2.0.4
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(DoubletFinder)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
})

set.seed(1234)

data_dir <- "data"
out_dir  <- file.path("results", "01_QC")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Parameters
# ------------------------------------------------------------------------------
datasets   <- c("GSE168191", "GSE178360", "GSE196341", "GSE196638")
file_paths <- file.path(data_dir, paste0(datasets, "_noQC.rds"))

# Expected doublet rate per 1,000 recovered cells (10x Genomics, ~0.8%)
doublet_rate_per_1000 <- c(GSE168191 = 0.008,
                           GSE178360 = 0.008,
                           GSE196341 = 0.008,
                           GSE196638 = 0.008)

df_PCs        <- 1:15   # PCs used for per-sample clustering and DoubletFinder
df_res        <- 0.3    # clustering resolution for homotypic doublet modelling
min_cells     <- 300    # samples with fewer cells are not run through DoubletFinder
remove_levels <- c("Doublet_high", "Doublet_low")

# ------------------------------------------------------------------------------
# Load and merge datasets
# ------------------------------------------------------------------------------
objs <- lapply(file_paths, readRDS)
for (i in seq_along(objs)) {
  DefaultAssay(objs[[i]]) <- "RNA"
  objs[[i]]$dataset <- datasets[i]
}

combined <- merge(x = objs[[1]], y = objs[-1], add.cell.ids = datasets)
rm(objs); gc()

combined$orig.ident <- combined$sample_id

# Recovered cells per sample before QC (used for expected doublet rate)
n_recovered <- table(combined$orig.ident)

# ------------------------------------------------------------------------------
# QC metrics and filtering
# ------------------------------------------------------------------------------
combined <- PercentageFeatureSet(combined, "^MT-",          col.name = "percent.mt")
combined <- PercentageFeatureSet(combined, "^RP[SL]",       col.name = "percent_ribo")
combined <- PercentageFeatureSet(combined, "^HB[^(P|E|S)]", col.name = "percent_hb")
combined <- PercentageFeatureSet(combined, "^(PPBP|PF4)$",  col.name = "percent_plat")

qc_feats <- c("nFeature_RNA", "nCount_RNA", "percent.mt",
              "percent_ribo", "percent_hb", "percent_plat")

p <- VlnPlot(combined, group.by = "orig.ident", split.by = "disease",
             features = qc_feats, pt.size = 0.1, ncol = 3)
ggsave(file.path(out_dir, "QC_violin_before_filtering.png"),
       p, width = 12, height = 8, dpi = 300)

combined <- subset(
  combined,
  subset = nFeature_RNA >= 200 & nFeature_RNA <= 7000 &
           nCount_RNA   >= 500 & nCount_RNA   <= 60000 &
           percent.mt   <= 20
)

# ------------------------------------------------------------------------------
# Per-sample doublet detection (DoubletFinder)
# ------------------------------------------------------------------------------
combined    <- JoinLayers(combined)
seurat_list <- SplitObject(combined, split.by = "orig.ident")
rm(combined); gc()

seurat_list <- lapply(seurat_list, function(object) {
  object <- NormalizeData(object, normalization.method = "LogNormalize", scale.factor = 10000)
  object <- FindVariableFeatures(object, selection.method = "vst", nfeatures = 2000)
  object <- ScaleData(object)
  object <- RunPCA(object, npcs = 30, features = VariableFeatures(object), verbose = FALSE)
  object <- FindNeighbors(object, reduction = "pca", dims = df_PCs)
  object <- FindClusters(object, resolution = df_res)
  object <- RunUMAP(object, dims = df_PCs)
  object
})

run_doubletfinder <- function(seurat_obj, n_rec, rate_per_1000,
                              PCs = df_PCs, res = df_res, seed = 1234) {
  set.seed(seed)

  # pK selection
  sweep.res   <- paramSweep(seurat_obj, PCs = PCs, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
  bcmvn       <- find.pK(sweep.stats)
  pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))

  # Expected number of doublets
  n_cells  <- ncol(seurat_obj)
  rate     <- rate_per_1000 * n_rec / 1000
  nExp_poi <- round(rate * n_cells)

  homotypic.prop <- modelHomotypic(seurat_obj@meta.data[[paste0("RNA_snn_res.", res)]])
  nExp_adj       <- round(nExp_poi * (1 - homotypic.prop))

  # Run 1: unadjusted nExp
  old_cols   <- colnames(seurat_obj@meta.data)
  seurat_obj <- doubletFinder(seurat_obj, PCs = PCs, pN = 0.25, pK = pK,
                              nExp = nExp_poi, sct = FALSE)
  new1     <- setdiff(colnames(seurat_obj@meta.data), old_cols)
  pann_col <- grep("^pANN", new1, value = TRUE)
  df1      <- grep("^DF.classifications", new1, value = TRUE)
  seurat_obj$DF_low <- seurat_obj@meta.data[[df1]]
  seurat_obj$pANN   <- seurat_obj@meta.data[[pann_col]]

  # Run 2: homotypic-adjusted nExp, reusing pANN from run 1
  if (nExp_adj > 0 && nExp_adj < nExp_poi) {
    old_cols   <- colnames(seurat_obj@meta.data)
    seurat_obj <- doubletFinder(seurat_obj, PCs = PCs, pN = 0.25, pK = pK,
                                nExp = nExp_adj, reuse.pANN = pann_col, sct = FALSE)
    new2 <- setdiff(colnames(seurat_obj@meta.data), old_cols)
    df2  <- grep("^DF.classifications", new2, value = TRUE)
    seurat_obj$DF_high <- seurat_obj@meta.data[[df2]]
  } else {
    seurat_obj$DF_high <- seurat_obj$DF_low
  }

  seurat_obj$Doublet <- ifelse(seurat_obj$DF_high == "Doublet", "Doublet_high",
                        ifelse(seurat_obj$DF_low  == "Doublet", "Doublet_low", "Singlet"))

  seurat_obj$DF_pK      <- pK
  seurat_obj$DF_nExp    <- nExp_poi
  seurat_obj$DF_nExpAdj <- nExp_adj
  drop <- grep("^(pANN_|DF.classifications_)", colnames(seurat_obj@meta.data), value = TRUE)
  seurat_obj@meta.data[drop] <- NULL

  message(sprintf("%s: %d cells, pK = %.3f, nExp = %d, nExp_adj = %d",
                  unique(seurat_obj$orig.ident), n_cells, pK, nExp_poi, nExp_adj))
  seurat_obj
}

results <- lapply(names(seurat_list), function(s) {
  obj <- seurat_list[[s]]
  if (ncol(obj) < min_cells) {
    message(s, ": only ", ncol(obj), " cells, DoubletFinder skipped")
    obj$Doublet <- "Not_tested"
    return(obj)
  }
  run_doubletfinder(obj,
                    n_rec         = as.numeric(n_recovered[s]),
                    rate_per_1000 = doublet_rate_per_1000[[unique(obj$dataset)]])
})
names(results) <- names(seurat_list)
rm(seurat_list); gc()

# Summary table and per-sample UMAPs of doublet calls
df_summary <- do.call(rbind, lapply(results, function(obj) {
  tab <- table(factor(obj$Doublet,
                      levels = c("Singlet", "Doublet_low", "Doublet_high", "Not_tested")))
  data.frame(sample      = unique(obj$orig.ident),
             dataset     = unique(obj$dataset),
             n_cells     = ncol(obj),
             t(as.matrix(tab)),
             pct_removed = round(100 * sum(obj$Doublet %in% remove_levels) / ncol(obj), 2))
}))
write.csv(df_summary, file.path(out_dir, "DoubletFinder_summary.csv"), row.names = FALSE)

umap_plots <- lapply(results, function(obj) {
  DimPlot(obj, group.by = "Doublet", pt.size = 0.3,
          cols = c(Singlet = "grey80", Doublet_low = "orange",
                   Doublet_high = "red", Not_tested = "grey40")) +
    ggtitle(unique(obj$orig.ident)) +
    theme(legend.position = "none")
})
ggsave(file.path(out_dir, "DoubletFinder_UMAP_per_sample.pdf"),
       wrap_plots(umap_plots, ncol = 4),
       width = 16, height = 4 * ceiling(length(umap_plots) / 4), limitsize = FALSE)

# Keep singlets and merge
singlet_list <- lapply(results, function(obj) {
  obj <- obj[, !(obj$Doublet %in% remove_levels)]
  DietSeurat(obj, layers = "counts", assays = "RNA")
})
combined <- merge(singlet_list[[1]], y = singlet_list[-1])
rm(results, singlet_list); gc()

p <- VlnPlot(combined, group.by = "orig.ident", features = qc_feats,
             pt.size = 0.1, ncol = 1)
ggsave(file.path(out_dir, "QC_violin_after_filtering.png"),
       p, width = 20, height = 24, dpi = 300, limitsize = FALSE)

# ------------------------------------------------------------------------------
# Clinical group per sample (from GEO sample annotations)
# ------------------------------------------------------------------------------
group_map <- c(
  # GSE168191
  GSM5133588 = "COPD",    GSM5133589 = "COPD",    GSM5133590 = "COPD",
  GSM5133591 = "COPD",    GSM5133592 = "COPD",    GSM5133593 = "COPD",
  GSM5133594 = "COPD",    GSM5133595 = "Healthy", GSM5133596 = "Healthy",
  GSM5133597 = "Smoker",  GSM5133598 = "Smoker",  GSM5133599 = "Healthy",
  GSM5133600 = "Healthy", GSM5133603 = "Smoker",  GSM5133604 = "Healthy",
  # GSE196341
  GSM5873135 = "COPD", GSM5873136 = "COPD", GSM5873137 = "COPD",
  GSM5873138 = "COPD", GSM5873139 = "COPD", GSM5873140 = "COPD",
  GSM5873141 = "COPD", GSM5873142 = "COPD", GSM5873143 = "COPD",
  GSM5873144 = "COPD", GSM5873145 = "COPD", GSM5873146 = "COPD",
  # GSE178360
  GSM5388411 = "Smoker", GSM5388412 = "Healthy", GSM5388413 = "Smoker",
  # GSE196638
  GSM5897133 = "Healthy", GSM5897134 = "Healthy", GSM5897135 = "Healthy",
  GSM5897136 = "COPD",    GSM5897137 = "COPD",    GSM5897138 = "COPD"
)

combined$clinical_group <- unname(group_map[as.character(combined$sample_id)])
stopifnot(!any(is.na(combined$clinical_group)))

# Keep COPD and healthy non-smoker samples only
combined <- subset(combined, subset = clinical_group %in% c("COPD", "Healthy"))

# ------------------------------------------------------------------------------
# Sample and cell counts per dataset
# ------------------------------------------------------------------------------
meta <- combined@meta.data

sample_counts <- meta %>%
  dplyr::distinct(dataset, disease, orig.ident) %>%
  dplyr::count(dataset, disease, name = "n") %>%
  tidyr::pivot_wider(names_from = disease, values_from = n,
                     values_fill = 0, names_prefix = "n_samples_")

cell_counts <- meta %>%
  dplyr::count(dataset, disease, name = "n") %>%
  tidyr::pivot_wider(names_from = disease, values_from = n,
                     values_fill = 0, names_prefix = "n_cells_")

sample_summary <- dplyr::left_join(sample_counts, cell_counts, by = "dataset") %>%
  dplyr::arrange(dataset)
sample_summary <- dplyr::bind_rows(
  sample_summary,
  sample_summary %>%
    dplyr::summarise(dataset = "Total",
                     dplyr::across(where(is.numeric), ~ sum(.x, na.rm = TRUE)))
)
write.csv(sample_summary, file.path(out_dir, "sample_summary.csv"), row.names = FALSE)

saveRDS(combined, file.path(out_dir, "singlets_COPD_Healthy.rds"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
