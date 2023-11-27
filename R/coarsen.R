#' @name coarsen
#' @aliases coarsen,CoverageExperiment-method
#' @rdname CoverageExperiment
#' 
#' @param x a `CoverageExperiment` object
#' @param window an integer to coarsen coverage by. 
#' @param BPPARAM Passed to BiocParallel.
#' 
#' @export
#' @examples 
#' data(ce)
#' coarsen(ce, window = 10)
setMethod("coarsen", signature(x = "CoverageExperiment"), 
    function(x, window = 1, BPPARAM = BiocParallel::bpparam()) {
        combs <- expand.grid(colData(x)$track, rowData(x)$features) |> 
            setNames(c("tracks", "features"))
        l <- BiocParallel::bplapply(seq_len(nrow(combs)), function(K) {
            t <- combs[K, "tracks"]
            f <- combs[K, "features"]
            mat <- .coarsen_mat(
                assay(x, 'coverage')[f, t][[1]], 
                bin = window, FUN = mean, na.rm = TRUE
            ) |> t()
        }, BPPARAM = BPPARAM)
        names(l) <- paste(combs$tracks, combs$features, sep = '^')
        for (t in colData(x)$track) {
            for (f in rowData(x)$features) {
                name <- paste(t, f, sep = '^')
                assays(x)[['coverage']][f, t][[1]] <- l[[name]]
            }
        }
        return(x)
    }
)
