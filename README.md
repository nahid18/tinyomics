# tinyomics 🧬

> *Tiny omics datasets, big possibilities.*

A lightweight R package providing instant access to curated omics datasets for teaching, prototyping, and benchmarking. All data loads as **tidyomics-compatible** container classes.

## Installation

```r
devtools::install_github("nahid18/tinyomics")
```

## Usage

```r
library(tinyomics)

sce <- scrna$pbmc_2020()

se <- bulk$airway_2014()

spe <- spatial$dlpfc_2021()

cytof_data <- cytof$levine32_2015()
```

## Contributing

Found a dataset that should be included? Open an issue or PR!

## License

MIT
