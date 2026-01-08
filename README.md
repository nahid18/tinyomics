# tinyomics 🧬

> *Tiny omics datasets, big possibilities.*

A lightweight R package providing instant access to curated omics datasets for teaching, prototyping, and benchmarking. All data loads as **tidyomics-compatible** container classes.

## Installation

```r
devtools::install_github("nahid18/tinyomics")
```

## Quick Start

```r
library(tinyomics)

# The simplest API - just two arguments
pbmc <- tiny("scrna", "pbmc")

# Or use the path-style syntax
pbmc <- tiny("scrna/pbmc")

# Modality-specific shortcuts
pbmc <- tiny_scrna("pbmc")
airway <- tiny_bulk("airway")
dlpfc <- tiny_spatial("dlpfc")
cytof <- tiny_cytof("levine32")
```

## Available Datasets

```r
tiny_list()  # See everything
```

### Single-cell RNA-seq (`scrna`)
| Dataset | Description |
|---------|-------------|
| `zeisel` | Mouse brain (3,005 cells, Zeisel et al.) |
| `darmanis` | Human brain (466 cells, Darmanis et al.) |
| `pbmc` | Human PBMC CITE-seq (Mair et al.) |
| `baron` | Human pancreas (Baron et al.) |
| `muraro` | Human pancreas (Muraro et al.) |
| `campbell` | Mouse hypothalamus (Campbell et al.) |

### Bulk RNA-seq (`bulk`)
| Dataset | Description |
|---------|-------------|
| `airway` | Airway smooth muscle (8 samples) |
| `tcga_brca` | TCGA Breast Cancer |
| `tcga_luad` | TCGA Lung Adenocarcinoma |
| `tcga_paad` | TCGA Pancreatic Cancer |

### Spatial Transcriptomics (`spatial`)
| Dataset | Description |
|---------|-------------|
| `dlpfc` | Human DLPFC Visium (Maynard et al.) |
| `mouse_embryo` | Mouse embryo seqFISH |
| `brain_libd` | LIBD Human DLPFC (47k spots) |
| `ovarian` | Mouse brain coronal Visium |

### CyTOF / Mass Cytometry (`cytof`)
| Dataset | Description |
|---------|-------------|
| `levine32` | Bone marrow (32 markers, 265k cells) |
| `levine13` | Bone marrow (13 markers) |
| `samusik` | Mouse bone marrow (39 markers) |
| `bodenmiller` | BCR-XL stimulation (16 samples) |
| `krieg` | Melanoma anti-PD-1 response |

### ChIP-seq (`chipseq`)
| Dataset | Description |
|---------|-------------|
| `tamoxifen` | ER ChIP tamoxifen resistance |

### Genotype (`genotype`)
| Dataset | Description |
|---------|-------------|
| `chr22` | 1000 Genomes chr22 variants |
| `structural` | Structural variants example |

### Multi-omics (`multiomics`)
| Dataset | Description |
|---------|-------------|
| `tcga_brca_multi` | BRCA RNA + Methylation + RPPA |
| `tcga_acc_multi` | ACC RNA + RPPA |

## Design Philosophy

### 🪶 Lightweight
Dependencies are loaded **dynamically**. Installing `tinyomics` doesn't pull in every Bioconductor package - only what you actually use.

### 🎯 Consistent API
Every dataset loads with the same pattern:
```r
data <- tiny("modality", "name")
```

### 📦 Tidyomics-Ready
All data returns as proper Bioconductor containers compatible with tidyomics:
- `SingleCellExperiment` for single-cell
- `SummarizedExperiment` for bulk/CyTOF
- `SpatialExperiment` for spatial
- `MultiAssayExperiment` for multi-omics
- `VCF` for genotype data

### 🔍 Discoverable
```r
tiny_list()              # What's available?
tiny_list("scrna")       # Just scRNA-seq
tiny_info("scrna/pbmc")  # Tell me more
```

## Example Workflows

### Quick scRNA-seq exploration
```r
library(tinyomics)
library(tidySingleCellExperiment)
library(dplyr)

tiny_scrna("pbmc") |>
  filter(label.fine != "Unknown") |>
  group_by(label.fine) |>
  summarise(n = n(), mean_counts = mean(nCount_RNA))
```

### Spatial visualization
```r
library(tinyomics)
library(ggspavis)

spe <- tiny_spatial("dlpfc")
plotSpots(spe, annotate = "ground_truth")
```

### CyTOF analysis
```r
library(tinyomics)
library(CATALYST)

cytof <- tiny_cytof("bodenmiller")
# Ready for clustering, differential analysis, etc.
```

## Contributing

Found a dataset that should be included? Open an issue or PR!

## License

MIT
