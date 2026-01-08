#' @title Dynamic Package Loading Utilities
#' @description Internal functions for lazy loading of Bioconductor packages.
#' @name utils
#' @keywords internal
NULL

#' Require a package, installing if necessary
#' 
#' @param pkg Package name
#' @param bioc Whether this is a Bioconductor package
#' @return TRUE if successful
#' @keywords internal
.require_pkg <- function(pkg, bioc = TRUE) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(sprintf("Installing %s...", pkg))
    if (bioc) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        utils::install.packages("BiocManager", repos = "https://cloud.r-project.org")
      }
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    } else {
      utils::install.packages(pkg, repos = "https://cloud.r-project.org")
    }
  }
  invisible(TRUE)
}

#' Load data from a package
#' 
#' @param pkg Package name
#' @param dataset Dataset name within the package
#' @param loader Optional custom loader function
#' @return The loaded data object
#' @keywords internal
.load_data <- function(pkg, dataset = NULL, loader = NULL) {
  .require_pkg(pkg)
  
  if (!is.null(loader)) {
    return(loader())
  }
  
  if (!is.null(dataset)) {
    e <- new.env()
    data(list = dataset, package = pkg, envir = e)
    return(get(dataset, envir = e))
  }
  
  stop("Either dataset or loader must be provided")
}

#' Format dataset info for display
#' 
#' @param info List containing dataset info
#' @return Formatted string
#' @keywords internal
.format_info <- function(info) {
  paste0(
    "\n", info$name, "\n",
    "  Description: ", info$description, "\n",
    "  Package: ", info$package, "\n",
    "  Class: ", info$class, "\n"
  )
}
