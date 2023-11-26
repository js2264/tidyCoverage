#' aggregate 
#' 
#' Bin coverage contained in a `CoverageExperiment` into an
#' `AggregatedCoverage` object.
#' 
#' @name AggregatedCoverage
#' @aliases aggregate,CoverageExperiment-method
#' @rdname AggregatedCoverage
#' 
#' @param x a `CoverageExperiment` object
#' @param ... ignored
#' @param bin an integer to bin each assay by. The `width` of the 
#' `AggregatedCoverage` object should be a multiple of `bin`. 
#' @return an `AggregatedCoverage` object
#' 
#' @importFrom S4Vectors aggregate
#' @export
#' @examples 
#' data(ce)
#' aggregate(ce, bin = 10)
setMethod("aggregate", signature(x = "CoverageExperiment"), 
    function(x, bin = 1, ...) {
        m0 <- matrix(
            list(), 
            nrow = length(rowData(x)$features), ncol = length(colData(x)$track)
        )
        colnames(m0) <- colData(x)$track
        rownames(m0) <- rowData(x)$features
        assays <- c("mean", "median", "min", "max", "sd", "se", "ci_low", "ci_high")
        l_assays <- vector("list", length = length(assays))
        names(l_assays) <- assays
        for(assay in assays) l_assays[[assay]] <- m0
        for (t in colData(x)$track) {
            for (f in rowData(x)$features) {
                m <- assay(x, 'coverage')[f, t][[1]]
                df <- .summarize_cov(m)
                for (assay in assays) {
                    vec <- .coarsen(df[[assay]], bin = bin, FUN = mean, na.rm = TRUE)
                    l_assays[[assay]][f, t][[1]] <- vec
                }
            }
        }
        assays(x) <- l_assays
        AC <- methods::new(
            "AggregatedCoverage",
            SummarizedExperiment::SummarizedExperiment(
                rowRanges = rowRanges(x),
                colData = colData(x),
                assays = l_assays
            )
        )
        return(AC)
    }
)
