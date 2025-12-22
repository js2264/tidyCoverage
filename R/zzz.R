#' tidyCoverage: Tidyomics-based analysis of coverage tracks over genomic features
#'
#' @description
#' The tidyCoverage package provides a bridge between Bioconductor
#' SummarizedExperiment objects and the tidyomics project. 
#' 
#' @references
#' Serizay J, Koszul R., Epigenomics coverage data extraction and aggregation 
#' in R with tidyCoverage. Bioinformatics 40, btae487 (2024).
#' \doi{10.1093/bioinformatics/btae487}
#'
#' @seealso
#' Useful links:
#' \itemize{
#'   \item \url{https://github.com/js2264/tidyCoverage}
#'   \item Report bugs at \url{https://github.com/js2264/tidyCoverage/issues}
#' }
#'
#' @author Jacques Serizay 
#'
#' @docType package
#' @name tidyCoverage-package
#' @aliases tidyCoverage
#' @keywords internal
"_PACKAGE"

#' @importFrom utils packageDescription
.onAttach <- function(libname, pkgname) {
    attached <- tidyverse_attach()
    
    # Print loading message about printing
    cli::cli_alert_info(paste0(
        "tidyCoverage says: The tidy printing is now handled externally.\n",
        "If you want to visualize the data in a tidy way, do `library(tidyprint)`.\n",
        "See `https://github.com/tidyomics/tidyprint` for more information."
    ))
}
