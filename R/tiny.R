#' @keywords internal
"_PACKAGE"

#' Load any omics dataset with a simple interface
#'
#' @description The main entry point for loading datasets. Use `tiny("modality", "dataset")`
#' or `tiny("modality/dataset")` syntax for quick access to any dataset.
#'
#' @param modality The data modality (e.g., "scrna", "bulk", "spatial", "cytof").
#'   Can also be a path-style string like "scrna/zeisel".
#' @param name The dataset name within the modality. Optional if using path syntax.
#' @param verbose Print information about the dataset being loaded.
#'
#' @return A Bioconductor container object (SummarizedExperiment, SingleCellExperiment, 
#'   SpatialExperiment, MultiAssayExperiment, or VCF).
#'
#' @examples
#' \dontrun{
#' # These are all equivalent:
#' pbmc <- tiny("scrna", "pbmc")
#' pbmc <- tiny("scrna/pbmc")
#' 
#' # Other examples:
#' airway <- tiny("bulk", "airway")
#' dlpfc <- tiny("spatial", "dlpfc")
#' cytof <- tiny("cytof", "levine32")
#' 
#' # List available datasets:
#' tiny_list()
#' tiny_list("scrna")
#' }
#'
#' @export
tiny <- function(modality, name = NULL, verbose = TRUE) {
  # Handle path-style input: "scrna/zeisel"
  if (is.null(name) && grepl("/", modality)) {
    parts <- strsplit(modality, "/")[[1]]
    modality <- parts[1]
    name <- parts[2]
  }
  
  if (is.null(name)) {
    stop("Please specify a dataset name. Use tiny_list('", modality, "') to see options.")
  }
  
  # Get the dataset entry
  entry <- .get_entry(modality, name)
  
  if (verbose) {
    message(sprintf("Loading %s/%s from %s...", modality, name, entry$package))
  }
  
  # Load and return the data
  .require_pkg(entry$package)
  data <- entry$loader()
  
  if (verbose) {
    message(sprintf("Loaded: %s", entry$description))
    message(sprintf("Class: %s", class(data)[1]))
  }
  
  return(data)
}


#' Load single-cell RNA-seq datasets
#'
#' @description Convenience function for scRNA-seq data.
#'
#' @param name Dataset name: "zeisel", "darmanis", "pbmc", "baron", "muraro", "campbell"
#' @param verbose Print loading info
#'
#' @return SingleCellExperiment object
#'
#' @examples
#' \dontrun{
#' pbmc <- tiny_scrna("pbmc")
#' brain <- tiny_scrna("zeisel")
#' }
#'
#' @export
tiny_scrna <- function(name = "pbmc", verbose = TRUE) {
  tiny("scrna", name, verbose)
}


#' Load bulk RNA-seq datasets
#'
#' @description Convenience function for bulk RNA-seq data.
#'
#' @param name Dataset name: "airway", "tcga_brca", "tcga_luad", "tcga_paad"
#' @param verbose Print loading info
#'
#' @return RangedSummarizedExperiment or MultiAssayExperiment object
#'
#' @examples
#' \dontrun{
#' airway <- tiny_bulk("airway")
#' brca <- tiny_bulk("tcga_brca")
#' }
#'
#' @export
tiny_bulk <- function(name = "airway", verbose = TRUE) {
  tiny("bulk", name, verbose)
}


#' Load spatial transcriptomics datasets
#'
#' @description Convenience function for spatial data.
#'
#' @param name Dataset name: "dlpfc", "mouse_embryo", "brain_libd", "ovarian"
#' @param verbose Print loading info
#'
#' @return SpatialExperiment object
#'
#' @examples
#' \dontrun{
#' brain <- tiny_spatial("dlpfc")
#' embryo <- tiny_spatial("mouse_embryo")
#' }
#'
#' @export
tiny_spatial <- function(name = "dlpfc", verbose = TRUE) {
  tiny("spatial", name, verbose)
}


#' Load ATAC-seq datasets
#'
#' @description Convenience function for ATAC-seq data.
#'
#' @param name Dataset name
#' @param verbose Print loading info
#'
#' @return SummarizedExperiment or similar
#'
#' @export
tiny_atac <- function(name = "tcga_brca", verbose = TRUE) {
  tiny("atac", name, verbose)
}


#' Load CyTOF / mass cytometry datasets
#'
#' @description Convenience function for CyTOF data.
#'
#' @param name Dataset name: "levine32", "levine13", "samusik", "bodenmiller", "krieg"
#' @param verbose Print loading info
#'
#' @return SummarizedExperiment object
#'
#' @examples
#' \dontrun{
#' cytof <- tiny_cytof("levine32")
#' bcr <- tiny_cytof("bodenmiller")
#' }
#'
#' @export
tiny_cytof <- function(name = "levine32", verbose = TRUE) {
  
  tiny("cytof", name, verbose)
}


#' Load ChIP-seq datasets
#'
#' @description Convenience function for ChIP-seq data.
#'
#' @param name Dataset name: "tamoxifen", "tamoxifen_peaks"
#' @param verbose Print loading info
#'
#' @return SummarizedExperiment or DBA object
#'
#' @examples
#' \dontrun{
#' chip <- tiny_chipseq("tamoxifen")
#' }
#'
#' @export
tiny_chipseq <- function(name = "tamoxifen", verbose = TRUE) {
  tiny("chipseq", name, verbose)
}


#' Load perturbation / CRISPR screen datasets
#'
#' @description Convenience function for perturbation data.
#'
#' @param name Dataset name
#' @param verbose Print loading info
#'
#' @return SingleCellExperiment or similar
#'
#' @export
tiny_perturb <- function(name = "k562_crispr", verbose = TRUE) {
  tiny("perturb", name, verbose)
}


#' Load genotype / variant datasets
#'
#' @description Convenience function for VCF/genotype data.
#'
#' @param name Dataset name: "chr22", "structural"
#' @param verbose Print loading info
#'
#' @return VCF object
#'
#' @examples
#' \dontrun{
#' vcf <- tiny_genotype("chr22")
#' }
#'
#' @export
tiny_genotype <- function(name = "chr22", verbose = TRUE) {
  tiny("genotype", name, verbose)
}


#' Load multi-omics datasets
#'
#' @description Convenience function for multi-omics data.
#'
#' @param name Dataset name: "tcga_brca_multi", "tcga_acc_multi"
#' @param verbose Print loading info
#'
#' @return MultiAssayExperiment object
#'
#' @examples
#' \dontrun{
#' mae <- tiny_multiomics("tcga_brca_multi")
#' }
#'
#' @export
tiny_multiomics <- function(name = "tcga_brca_multi", verbose = TRUE) {
  tiny("multiomics", name, verbose)
}


#' List available datasets
#'
#' @description Show available modalities and datasets.
#'
#' @param modality Optional: show datasets for a specific modality only.
#'
#' @return Invisibly returns a data.frame of datasets. Also prints to console.
#'
#' @examples
#' \dontrun{
#' tiny_list()              # All datasets
#' tiny_list("scrna")       # Just scRNA-seq
#' tiny_list("spatial")     # Just spatial
#' }
#'
#' @export
tiny_list <- function(modality = NULL) {
  if (is.null(modality)) {
    modalities <- .get_modalities()
  } else {
    modalities <- modality
  }
  
  rows <- list()
  for (mod in modalities) {
    datasets <- .get_datasets(mod)
    for (nm in names(datasets)) {
      entry <- datasets[[nm]]
      rows[[length(rows) + 1]] <- data.frame(
        modality = mod,
        name = nm,
        class = entry$class,
        package = entry$package,
        description = entry$description,
        stringsAsFactors = FALSE
      )
    }
  }
  
  df <- do.call(rbind, rows)
  
  # Print nicely
  cat("\n")
  cat(sprintf("tinyomics: %d datasets across %d modalities\n", 
              nrow(df), length(unique(df$modality))))
  cat(strrep("-", 60), "\n\n")
  
  for (mod in unique(df$modality)) {
    cat(sprintf("== %s ==\n", toupper(mod)))
    sub <- df[df$modality == mod, ]
    for (i in seq_len(nrow(sub))) {
      cat(sprintf("  %-15s %s\n", sub$name[i], sub$description[i]))
    }
    cat("\n")
  }
  
  cat(strrep("-", 60), "\n")
  cat("Use: tiny('modality', 'name') or tiny('modality/name')\n")
  cat("Example: tiny('scrna', 'pbmc') or tiny('scrna/pbmc')\n\n")
  
  invisible(df)
}


#' Get detailed info about a dataset
#'
#' @description Print detailed information about a specific dataset.
#'
#' @param modality Modality name or path-style string
#' @param name Dataset name (optional if using path syntax)
#'
#' @return Invisibly returns the info list
#'
#' @examples
#' \dontrun{
#' tiny_info("scrna", "pbmc")
#' tiny_info("scrna/pbmc")
#' }
#'
#' @export
tiny_info <- function(modality, name = NULL) {
  # Handle path-style input
  if (is.null(name) && grepl("/", modality)) {
    parts <- strsplit(modality, "/")[[1]]
    modality <- parts[1]
    name <- parts[2]
  }
  
  entry <- .get_entry(modality, name)
  
  cat("\n")
  cat(sprintf("Dataset: %s/%s\n", modality, name))
  cat(strrep("-", 40), "\n")
  cat(sprintf("Description: %s\n", entry$description))
  cat(sprintf("Package:     %s\n", entry$package))
  cat(sprintf("Class:       %s\n", entry$class))
  cat("\n")
  cat(sprintf("Load with:   tiny('%s', '%s')\n", modality, name))
  cat(sprintf("         or: tiny('%s/%s')\n", modality, name))
  cat(sprintf("         or: tiny_%s('%s')\n", modality, name))
  cat("\n")
  
  invisible(entry)
}