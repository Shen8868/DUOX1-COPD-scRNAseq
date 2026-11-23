# DUOX1 in COPD alveolar type 2 cells: scRNA-seq analysis code

Analysis code for the single-cell RNA-seq part of:

> *Defective epithelial cell clearance drives emphysema and distinct pulmonary vascular phenotype in COPD.* Edma Loku1,*, Cheng-Yu Wu1,*, Vinita Sharma1, Chudamani Nuguri1, Anis Cilic1, Vidya Lakshmi1, Dana Cojocaru1, Henrike Tries1, Zixiong Shen1, Simone Kraut1, Baktybek Kojonazarov1, Marek Bartkuhn1, Andreas Guenther1, Claudia F Garcia Castro1, Werner Seeger1, Marija Gredic1, Hossein A Ghofrani1, Natascha Sommer1, Norbert Weissmann1 and Stefan Hadzic1,#. iScience (2026).

The scripts integrate four public human lung scRNA-seq datasets, annotate
epithelial cell types, and compare DUOX1/DUOXA1 double-positive
(`DUOX1_functional`) and other (`DUOX1_dysfunctional`) alveolar type 2 (AT2)
cells in COPD.

## Datasets

| GEO accession | Samples used |
|---|---|
| [GSE168191](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE168191) | COPD, healthy |
| [GSE178360](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE178360) | healthy |
| [GSE196341](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE196341) | COPD |
| [GSE196638](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE196638) | COPD, healthy |

Samples from smokers without COPD were excluded. The sample-to-group assignment
(by GSM accession) is listed in `01_QC_doublet_removal.R`.

## Repository structure

```
.
├── 01_QC_doublet_removal.R          # QC, per-sample doublet removal, sample selection
├── 02_integration_annotation.R      # normalization, Harmony integration, clustering, annotation
├── 03_AT2_DUOX1_DEG_enrichment.R    # DUOX1 grouping in AT2, DEG, quadrant plot, enrichment
├── data/                            # input Seurat objects (not included, see below)
└── results/                         # created by the scripts
```

## Input data

The count matrices deposited by the original authors were downloaded from GEO
and read into Seurat (`CreateSeuratObject`) without any filtering. All samples
of one GEO series were merged into one object and saved as:

```
data/GSE168191_noQC.rds
data/GSE178360_noQC.rds
data/GSE196341_noQC.rds
data/GSE196638_noQC.rds
```

Each object needs two metadata columns:

- `sample_id`: GSM accession of the sample
- `disease`: `COPD` or `Control`

## How to run

Run the scripts in order from the repository root:

```r
source("01_QC_doublet_removal.R")
source("02_integration_annotation.R")
source("03_AT2_DUOX1_DEG_enrichment.R")
```

Each script reads the output of the previous one, so they can also be run in
separate R sessions. Output files are written to `results/01_QC/`,
`results/02_integration/` and `results/03_AT2_DUOX1/`.

| Step | Main outputs |
|---|---|
| 01 | QC violin plots, DoubletFinder summary, sample summary, filtered object |
| 02 | Integrated and annotated objects, UMAPs, marker plots |
| 03 | DEG tables, four-quadrant plot, GO / KEGG / WikiPathways results |

## Requirements

R (version) with the following packages:

| Package | Version |
|---|---|
| Seurat | 5.x.x |
| SeuratObject | 5.x.x |
| DoubletFinder | 2.0.x (GitHub: chris-mcginnis-ucsf/DoubletFinder) |
| harmony | x.x.x |
| clusterProfiler | x.x.x |
| org.Hs.eg.db | x.x.x |
| enrichplot | x.x.x |
| ggplot2, ggrepel, patchwork, dplyr, tidyr | |

The full `sessionInfo()` used for the analysis is in `sessionInfo.txt`.
Installing [presto](https://github.com/immunogenomics/presto) is recommended
to speed up `FindMarkers`.

KEGG and WikiPathways enrichment download annotation data online, so results
may differ slightly with later database releases.

## Notes on reproducibility

- Random seeds are fixed (`set.seed(1234)`) before all stochastic steps.
- Cell type annotation in `02_integration_annotation.R` assigns labels to
  cluster numbers. Cluster numbering can change with different package
  versions. If you use other versions, check the marker gene plots before
  applying the labels.

## License

This code is released under the MIT License (see `LICENSE`).

## Contact

For questions, please open an issue in this repository.
