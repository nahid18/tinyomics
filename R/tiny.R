#' @name tinyomics
#' @title tinyomics
#' @aliases scrna bulk spatial cytof chipseq atac perturb tiny_list tiny_search
#' @description Unified access to omics datasets
NULL

#' @export
scrna <- .create_modality("scrna")

#' @export
bulk <- .create_modality("bulk")

#' @export
spatial <- .create_modality("spatial")

#' @export
cytof <- .create_modality("cytof")

