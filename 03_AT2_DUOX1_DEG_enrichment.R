# ==============================================================================
# 03_AT2_DUOX1_DEG_enrichment.R
#
# AT2 cells: DUOX1_functional (DUOX1+ DUOXA1+) vs DUOX1_dysfunctional (others)
#   DEG A: AT2, COPD vs Control
#   DEG B: COPD AT2, DUOX1_dysfunctional vs DUOX1_functional
#   Four-quadrant plot of A vs B, then GO / KEGG / WikiPathways enrichment
#   of genes down-regulated in both comparisons.
#
# Input : results/02_integration/epithelial_annotated.rds
# Output: results/03_AT2_DUOX1/
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
})

set.seed(1234)

in_file <- file.path("results", "02_integration", "epithelial_annotated.rds")
out_dir <- file.path("results", "03_AT2_DUOX1")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Parameters
# ------------------------------------------------------------------------------
lfc_cut   <- 0.25    # |avg_log2FC| threshold
logp_cut  <- 5       # -log10(P) threshold
use_adj   <- TRUE    # TRUE: Bonferroni-adjusted P (Seurat p_val_adj); FALSE: raw P
min_pct   <- 0.05    # minimum fraction of cells expressing a gene in either group
def_genes <- c("DUOX1", "DUOXA1")   # grouping genes, excluded downstream

# ------------------------------------------------------------------------------
# 1. AT2 subset and DUOX1 status
# ------------------------------------------------------------------------------
epi_combined <- readRDS(in_file)

at2 <- subset(epi_combined, subset = hand_celltype == "AT2")
DefaultAssay(at2) <- "RNA"
stopifnot(all(def_genes %in% rownames(at2)))

cnt <- LayerData(at2, assay = "RNA", layer = "counts")[def_genes, , drop = FALSE]
at2$DUOX1_status <- ifelse(cnt["DUOX1", ] > 0 & cnt["DUOXA1", ] > 0,
                           "DUOX1_functional", "DUOX1_dysfunctional")
at2$DUOX1_status <- factor(at2$DUOX1_status,
                           levels = c("DUOX1_functional", "DUOX1_dysfunctional"))

frac_tab <- at2@meta.data %>%
  dplyr::group_by(orig.ident, disease) %>%
  dplyr::summarise(n_cells        = dplyr::n(),
                   n_functional   = sum(DUOX1_status == "DUOX1_functional"),
                   pct_functional = round(100 * n_functional / n_cells, 2),
                   .groups = "drop")
write.csv(frac_tab, file.path(out_dir, "DUOX1_functional_fraction_per_sample.csv"),
          row.names = FALSE)

p_frac <- ggplot(frac_tab, aes(disease, pct_functional, fill = disease)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, alpha = 0.6) +
  geom_jitter(width = 0.15, size = 2) +
  labs(x = NULL, y = "% DUOX1_functional in AT2 (per sample)") +
  theme_classic() +
  theme(legend.position = "none")
ggsave(file.path(out_dir, "DUOX1_functional_fraction_per_sample.pdf"),
       p_frac, width = 3.5, height = 4)

# ------------------------------------------------------------------------------
# 2. Differential expression
# ------------------------------------------------------------------------------
# DEG A: AT2, COPD vs Control
Idents(at2) <- "disease"
deg_A <- FindMarkers(at2, ident.1 = "COPD", ident.2 = "Control",
                     test.use = "wilcox", logfc.threshold = 0, min.pct = min_pct)

# DEG B: COPD AT2, DUOX1_dysfunctional vs DUOX1_functional
at2_copd <- subset(at2, subset = disease == "COPD")
Idents(at2_copd) <- "DUOX1_status"
deg_B <- FindMarkers(at2_copd, ident.1 = "DUOX1_dysfunctional", ident.2 = "DUOX1_functional",
                     test.use = "wilcox", logfc.threshold = 0, min.pct = min_pct)

tidy_deg <- function(df) {
  p <- if (use_adj) df$p_val_adj else df$p_val
  p[p == 0] <- .Machine$double.xmin   # avoid -log10(0) = Inf
  df$logP <- -log10(p)
  df$gene <- rownames(df)
  df
}
deg_A <- tidy_deg(deg_A)
deg_B <- tidy_deg(deg_B)

write.csv(deg_A, file.path(out_dir, "DEG_A_AT2_COPD_vs_Control.csv"), row.names = FALSE)
write.csv(deg_B, file.path(out_dir, "DEG_B_COPD_AT2_dysfunctional_vs_functional.csv"),
          row.names = FALSE)

# ------------------------------------------------------------------------------
# 3. Quadrant classification
# ------------------------------------------------------------------------------
merged <- dplyr::inner_join(
  deg_A %>% dplyr::select(gene, avg_log2FC, logP, pct.1, pct.2),
  deg_B %>% dplyr::select(gene, avg_log2FC, logP, pct.1, pct.2),
  by = "gene", suffix = c("_A", "_B")
) %>%
  dplyr::filter(!gene %in% def_genes)

sig <- merged$logP_A > logp_cut & merged$logP_B > logp_cut
fcA <- merged$avg_log2FC_A
fcB <- merged$avg_log2FC_B

merged$category <- "NS"
merged$category[sig & fcA >  lfc_cut & fcB >  lfc_cut] <- "Both_up"
merged$category[sig & fcA < -lfc_cut & fcB < -lfc_cut] <- "Both_down"
merged$category[sig & fcA >  lfc_cut & fcB < -lfc_cut] <- "A_up_B_down"
merged$category[sig & fcA < -lfc_cut & fcB >  lfc_cut] <- "A_down_B_up"
merged$category <- factor(merged$category,
                          levels = c("Both_up", "Both_down", "A_up_B_down", "A_down_B_up", "NS"))

write.csv(merged, file.path(out_dir, "DEG_A_B_merged_quadrant.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------
# 4. Four-quadrant plot
# ------------------------------------------------------------------------------
cat_cols <- c(Both_up = "#D7301F", Both_down = "#2C7FB8",
              A_up_B_down = "#F39C12", A_down_B_up = "#7B3294", NS = "grey80")

label_df <- merged %>%
  dplyr::filter(category != "NS") %>%
  dplyr::mutate(score = abs(avg_log2FC_A) + abs(avg_log2FC_B)) %>%
  dplyr::group_by(category) %>%
  dplyr::slice_max(score, n = 10, with_ties = FALSE) %>%
  dplyr::ungroup()

n_cat <- table(merged$category)
x_rng <- range(merged$avg_log2FC_A)
y_rng <- range(merged$avg_log2FC_B)

p_quad <- ggplot(merged, aes(avg_log2FC_A, avg_log2FC_B)) +
  geom_point(data = subset(merged, category == "NS"), color = cat_cols["NS"], size = 0.6) +
  geom_point(data = subset(merged, category != "NS"), aes(color = category), size = 1.3) +
  geom_hline(yintercept = 0, color = "grey40") +
  geom_vline(xintercept = 0, color = "grey40") +
  geom_hline(yintercept = c(-lfc_cut, lfc_cut), linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed", color = "grey60") +
  geom_text_repel(data = label_df, aes(label = gene, color = category),
                  size = 3, max.overlaps = 50, show.legend = FALSE,
                  segment.color = "grey50", box.padding = 0.3) +
  annotate("text", x = x_rng[2], y = y_rng[2], hjust = 1, vjust = 1,
           label = paste0("Both up: ", n_cat["Both_up"]), color = cat_cols["Both_up"]) +
  annotate("text", x = x_rng[1], y = y_rng[1], hjust = 0, vjust = 0,
           label = paste0("Both down: ", n_cat["Both_down"]), color = cat_cols["Both_down"]) +
  annotate("text", x = x_rng[2], y = y_rng[1], hjust = 1, vjust = 0,
           label = paste0("A up / B down: ", n_cat["A_up_B_down"]), color = cat_cols["A_up_B_down"]) +
  annotate("text", x = x_rng[1], y = y_rng[2], hjust = 0, vjust = 1,
           label = paste0("A down / B up: ", n_cat["A_down_B_up"]), color = cat_cols["A_down_B_up"]) +
  scale_color_manual(values = cat_cols, drop = FALSE) +
  labs(x = "log2FC  AT2: COPD vs Control (A)",
       y = "log2FC  COPD AT2: DUOX1_dysfunctional vs functional (B)",
       color = NULL,
       caption = sprintf("|log2FC| > %.2f and -log10(%s) > %g in both comparisons",
                         lfc_cut, ifelse(use_adj, "adj.P", "P"), logp_cut)) +
  theme_classic(base_size = 12)

ggsave(file.path(out_dir, "Quadrant_A_vs_B.pdf"), p_quad, width = 8, height = 6.5)
ggsave(file.path(out_dir, "Quadrant_A_vs_B.png"), p_quad, width = 8, height = 6.5, dpi = 300)

# ------------------------------------------------------------------------------
# 5. Enrichment of genes down-regulated in both comparisons
# ------------------------------------------------------------------------------
down_genes <- merged$gene[merged$category == "Both_down"]
writeLines(down_genes, file.path(out_dir, "Both_down_genes.txt"))

# Background: all genes tested in both comparisons
universe_symbol <- merged$gene
id_map <- bitr(universe_symbol, fromType = "SYMBOL", toType = "ENTREZID",
               OrgDb = org.Hs.eg.db)
down_entrez     <- unique(id_map$ENTREZID[id_map$SYMBOL %in% down_genes])
universe_entrez <- unique(id_map$ENTREZID)

save_enrich <- function(res, name, n_show = 20, split_ont = FALSE) {
  if (is.null(res) || nrow(as.data.frame(res)) == 0) {
    message(name, ": no significant terms")
    return(invisible(NULL))
  }
  write.csv(as.data.frame(res), file.path(out_dir, paste0(name, ".csv")), row.names = FALSE)
  p <- if (split_ont) {
    dotplot(res, showCategory = 10, split = "ONTOLOGY") +
      facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free_y")
  } else {
    dotplot(res, showCategory = n_show)
  }
  ggsave(file.path(out_dir, paste0(name, "_dotplot.pdf")), p + ggtitle(name),
         width = 8, height = ifelse(split_ont, 10, 7))
  invisible(res)
}

ego <- enrichGO(gene = down_genes, universe = universe_symbol,
                OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "ALL",
                pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
save_enrich(ego, "GO_Both_down", split_ont = TRUE)

ekegg <- enrichKEGG(gene = down_entrez, universe = universe_entrez, organism = "hsa",
                    pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
if (!is.null(ekegg)) ekegg <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
save_enrich(ekegg, "KEGG_Both_down")

ewp <- enrichWP(gene = down_entrez, universe = universe_entrez, organism = "Homo sapiens",
                pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
if (!is.null(ewp)) ewp <- setReadable(ewp, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
save_enrich(ewp, "WikiPathways_Both_down")

saveRDS(list(deg_A = deg_A, deg_B = deg_B, merged = merged,
             GO = ego, KEGG = ekegg, WP = ewp),
        file.path(out_dir, "AT2_DUOX1_results.rds"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
