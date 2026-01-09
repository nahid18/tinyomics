.require_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) utils::install.packages("BiocManager")
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

.require_github <- function(repo) {
  pkg <- basename(repo)
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (!requireNamespace("remotes", quietly = TRUE)) utils::install.packages("remotes")
    remotes::install_github(repo)
  }
}

.get_modalities <- function() names(.registry)

.get_datasets <- function(modality) {
  if (!modality %in% names(.registry)) stop("Unknown modality: ", modality)
  .registry[[modality]]
}

.get_entry <- function(modality, name) {
  datasets <- .get_datasets(modality)
  if (!name %in% names(datasets)) stop("Unknown dataset: ", name)
  datasets[[name]]
}

.create_modality <- function(modality) {
  structure(list(modality = modality), class = "tinyomics_modality")
}

#' @export
`$.tinyomics_modality` <- function(x, name) {
  modality <- x[["modality"]]
  entry <- .get_entry(modality, name)
  function(verbose = TRUE) {
    if (verbose) message("Loading ", modality, "/", name, "...")
    data <- entry()
    if (verbose) message("Done. Class: ", class(data)[1])
    data
  }
}

#' @export
names.tinyomics_modality <- function(x) names(.get_datasets(x[["modality"]]))

#' @export
print.tinyomics_modality <- function(x, ...) {
  modality <- x[["modality"]]
  datasets <- names(.get_datasets(modality))
  cat("tinyomics::", modality, " [", length(datasets), " datasets]\n", sep = "")
  cat(paste0("
## ", modality, "$", datasets, "()"), sep = "\n")
  invisible(x)
}