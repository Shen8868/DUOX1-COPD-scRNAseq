# ==============================================================================
# 02_integration_annotation.R
#
# Normalization, Harmony integration, clustering and cell type annotation.
#
# Input : results/01_QC/singlets_COPD_Healthy.rds
# Output: results/02_integration/combined_annotated.rds
#         results/02_integration/epithelial_annotated.rds
#
# Note: cluster numbers used for annotation below depend on the exact package
# versions (see sessionInfo.txt). With other versions, re-check marker genes.
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(ggplot2)
})

in_file <- file.path("results", "01_QC", "singlets_COPD_Healthy.rds")
out_dir <- file.path("results", "02_integration")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Map cluster IDs to labels; clusters not listed keep their ID
label_clusters <- function(clusters, label_list) {
  lut <- unlist(lapply(names(label_list), function(lab) {
    setNames(rep(lab, length(label_list[[lab]])), as.character(label_list[[lab]]))
  }))
  clusters <- as.character(clusters)
  out <- unname(lut[clusters])
  out[is.na(out)] <- clusters[is.na(out)]
  out
}

# Dot plot of one gene across cell types, split by disease
dotplot_by_disease <- function(object, gene, group.by = "hand_celltype") {
  df <- DotPlot(object, features = gene, group.by = group.by, split.by = "disease")$data
  cell_order    <- unique(as.character(object@meta.data[[group.by]]))
  disease_order <- unique(as.character(object$disease))

  df$celltype <- NA_character_
  df$disease  <- NA_character_
  for (d in disease_order) {
    idx <- endsWith(as.character(df$id), paste0("_", d))
    df$disease[idx]  <- d
    df$celltype[idx] <- substr(as.character(df$id[idx]), 1,
                               nchar(as.character(df$id[idx])) - nchar(d) - 1)
  }
  df$celltype <- factor(df$celltype, levels = cell_order)
  df$disease  <- factor(df$disease, levels = rev(disease_order))

  ggplot(df, aes(x = celltype, y = disease)) +
    geom_point(aes(size = pct.exp, color = avg.exp)) +
    scale_size(range = c(0, 9), name = "Percent expressed") +
    scale_color_gradient(low = "lightgrey", high = "red", name = "Average expression") +
    labs(x = NULL, y = NULL, title = gene) +
    theme_classic(base_size = 13) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.line = element_blank(), axis.ticks = element_blank())
}

# ------------------------------------------------------------------------------
# Normalization and Harmony integration
# ------------------------------------------------------------------------------
combined <- readRDS(in_file)
DefaultAssay(combined) <- "RNA"
combined[["RNA"]] <- split(combined[["RNA"]], f = combined$sample_id)

set.seed(1234)
combined <- NormalizeData(combined)
combined <- CellCycleScoring(combined,
                             g2m.features = cc.genes$g2m.genes,
                             s.features   = cc.genes$s.genes)
combined$CC.Difference <- combined$S.Score - combined$G2M.Score

combined <- FindVariableFeatures(combined, nfeatures = 2000)
combined <- ScaleData(combined,
                      vars.to.regress = c("nCount_RNA", "percent.mt", "CC.Difference"),
                      verbose = FALSE)
combined <- RunPCA(combined)

combined <- RunHarmony(combined, group.by.vars = "orig.ident", reduction = "pca",
                       dims.use = 1:50, reduction.save = "harmony", verbose = TRUE)
combined <- JoinLayers(combined)

ggsave(file.path(out_dir, "ElbowPlot.png"), ElbowPlot(combined, ndims = 50),
       width = 6, height = 4, dpi = 300)

set.seed(1234)
combined <- FindNeighbors(combined, reduction = "harmony", dims = 1:30)
set.seed(1234)
combined <- FindClusters(combined, resolution = 1, cluster.name = "harmony_clusters")
set.seed(1234)
combined <- RunUMAP(combined, reduction = "harmony", dims = 1:30,
                    reduction.name = "umap.harmony")

saveRDS(combined, file.path(out_dir, "combined_integrated.rds"))

# ------------------------------------------------------------------------------
# Major cell type annotation
# ------------------------------------------------------------------------------
p <- FeaturePlot(combined,
                 features = c("FXYD3", "EPCAM", "CLDN5", "PECAM1", "CD53",
                              "PTPRC", "COL1A2", "DCN", "MFAP4"),
                 raster = TRUE, reduction = "umap.harmony")
ggsave(file.path(out_dir, "FeaturePlot_major_lineage_markers.png"),
       p, width = 12, height = 10, dpi = 300)

main_labels <- list(
  Epithelial  = c(1, 8, 9, 14, 16, 23, 25),
  Endothelial = c(4, 10, 13, 17, 20),
  Immune      = c(0, 2, 3, 6, 7, 12, 18, 21, 24, 26, 27),
  Stromal     = c(5, 11, 19, 28),
  doublet     = c(15, 22, 29:51)   # clusters co-expressing markers of multiple lineages
)
combined$main_celltype <- label_clusters(combined$harmony_clusters, main_labels)
combined <- subset(combined, subset = main_celltype != "doublet")

# ------------------------------------------------------------------------------
# Epithelial sub-clustering and annotation
# ------------------------------------------------------------------------------
epi_combined <- subset(combined, subset = main_celltype == "Epithelial")

set.seed(1234)
epi_combined <- FindClusters(epi_combined, resolution = 0.3,
                             cluster.name = "harmony_clusters")

epi_labels <- list(
  AT2          = c(0, 6),
  AT1          = 2,
  AT2_AT1      = 5,
  Basal        = 7,
  `Sec/Goblet` = 4,
  RAS          = 3,
  Ciliated     = 1,
  doublet      = 8
)
epi_combined$hand_celltype <- label_clusters(epi_combined$harmony_clusters, epi_labels)

# Re-cluster the AT2_AT1 cluster to separate AT2 and transitional AT2/AT1 cells
at2_at1 <- subset(epi_combined, subset = hand_celltype == "AT2_AT1")
set.seed(1234)
at2_at1 <- FindClusters(at2_at1, resolution = 0.1, cluster.name = "harmony_clusters")
at2_at1$hand_celltype <- label_clusters(at2_at1$harmony_clusters,
                                        list(`AT2/AT1` = 0, AT2 = 1))

epi_combined@meta.data[colnames(at2_at1), "hand_celltype"] <-
  at2_at1@meta.data[colnames(at2_at1), "hand_celltype"]
rm(at2_at1)

# Transfer labels to the full object (non-epithelial cells keep their major type)
combined$hand_celltype <- combined$main_celltype
combined@meta.data[colnames(epi_combined), "hand_celltype"] <-
  epi_combined@meta.data[colnames(epi_combined), "hand_celltype"]

saveRDS(combined,     file.path(out_dir, "combined_annotated.rds"))
saveRDS(epi_combined, file.path(out_dir, "epithelial_annotated.rds"))

# ------------------------------------------------------------------------------
# Figures
# ------------------------------------------------------------------------------
ggsave(file.path(out_dir, "UMAP_main_celltype.png"),
       DimPlot(combined, reduction = "umap.harmony", group.by = "main_celltype", label = TRUE),
       width = 10, height = 8, dpi = 300)

ggsave(file.path(out_dir, "FeaturePlot_lineage_markers.png"),
       FeaturePlot(combined, features = c("EPCAM", "PECAM1", "PTPRC", "PDGFRA", "PDGFRB"),
                   reduction = "umap.harmony", ncol = 5),
       width = 25, height = 5, dpi = 300, limitsize = FALSE)

for (g in c("DUOX1", "DUOXA1")) {
  ggsave(file.path(out_dir, paste0("UMAP_", g, ".png")),
         FeaturePlot(combined, features = g, reduction = "umap.harmony"),
         width = 8, height = 7, dpi = 300)
  ggsave(file.path(out_dir, paste0("UMAP_epi_", g, "_by_disease.png")),
         FeaturePlot(epi_combined, features = g, split.by = "disease",
                     reduction = "umap.harmony"),
         width = 12, height = 6, dpi = 300)
  ggsave(file.path(out_dir, paste0("DotPlot_epi_", g, ".png")),
         dotplot_by_disease(epi_combined, g),
         width = 12, height = 3.5, dpi = 300)
}

ggsave(file.path(out_dir, "UMAP_epi_hand_celltype.png"),
       DimPlot(epi_combined, reduction = "umap.harmony", group.by = "hand_celltype", label = TRUE),
       width = 10, height = 8, dpi = 300)

epi_markers <- c(
  "KRT15", "KRT17", "KRT5", "KRT19",                               # Basal
  "FOXJ1", "RSPH1", "C20orf85", "C9orf24",                         # Ciliated
  "BPIFB1", "LCN2", "TSPAN8", "MUC5B",                             # Secretory/Goblet
  "AGER", "CLIC5", "SPOCK2", "PDPN", "CAV1", "HOPX",               # AT1
  "MFSD2A", "C11orf96", "SFTPB", "SFTPA2", "SFTA2", "NAPSA",
  "SFTPD", "SFTPC",                                                # AT2
  "KLK11", "SOX4", "MGP", "SCGB3A2", "SCGB1A1", "SCGB3A1"          # RAS
)
p <- DoHeatmap(epi_combined, features = epi_markers, cells = 1:500,
               size = 4, angle = 90) + NoLegend()
ggsave(file.path(out_dir, "Heatmap_epithelial_markers.svg"), p, width = 8, height = 8)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
