#' aggregate 
#' 
#' Bin assays of an `CoverageExperiment` to a coarser resolution
#' 
#' @name aggregate
#' @rdname aggregate
#' 
#' @param x a `CoverageExperiment` or CoverageExperiment` object
#' @param ... ignored
#' @param bin an integer to bin each assay by. The `width` of the 
#' `CoverageExperiment` object should be a multiple of `bin`. 
#' @return an `CoverageExperiment` object
#' 
#' @importFrom S4Vectors aggregate
#' @examples 
#' data(ac)
#' aggregate(ac, bin = 10)
NULL 

#' @name aggregate
#' @export

setMethod("aggregate", signature(x = "CoverageExperiment"), 
    function(x, bin = 1, ...) {
        m0 <- matrix(
            list(), 
            nrow = length(features), ncol = length(tracks)
        )
        colnames(m0) <- colData(x)$track
        rownames(m0) <- rowData(x)$features
        assays <- c("coverage", "mean", "median", "min", "max", "sd", "se", "ci_low", "ci_high")
        l_assays <- vector("list", length = length(assays))
        names(l_assays) <- assays
        l_assays[['coverage']] <- assay(x, 'coverage')
        for(assay in assays[-1]) l_assays[[assay]] <- m0
        for (t in colData(x)$track) {
            for (f in rowData(x)$features) {
                m <- assay(x, 'coverage')[f, t][[1]]
                df <- .summarize_cov(m)
                for (assay in assays[-1]) {
                    vec <- .coarsen(df[[assay]], bin = bin, FUN = mean, na.rm = TRUE)
                    l_assays[[assay]][f, t][[1]] <- vec
                }
            }
        }
        assays(x) <- l_assays
        return(x)
    }
)
