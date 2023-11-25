#' as_tibble
#' 
#' Coerce an `CoverageExperiment` object into a `tibble`
#' 
#' @rdname as_tibble
#' 
#' @param x an `CoverageExperiment` object
#' @param ... ignored
#' @return `tibble`
#' 
#' @importFrom dplyr as_tibble
#' @export
#' @examples 
#' data(ac)
#' as_tibble(ac)

as_tibble.CoverageExperiment <- function(x, ...) {
    tracks <- colData(x)$track
    features <- rowData(x)$features
    assays <- names(assays(x))
    w <- width(rowRanges(x)[[1]])[[1]]
    bin <- w / length(assay(x, 1)[1, 1][[1]])
    lapply(tracks, function(t) {
        lapply(features, function(f) {
            lapply(assays(x), `[[`, f, t) |> 
                stats::setNames(assays) |> 
                dplyr::bind_cols() |> 
                dplyr::mutate(
                    coord = seq(-w/2, w/2-1, bin), 
                    track = t, features = f
                ) |> 
                dplyr::relocate(coord) |> 
                dplyr::relocate(features) |> 
                dplyr::relocate(track)
        }) |> dplyr::bind_rows()
    }) |> dplyr::bind_rows()
}

#' show
#' 
#' show method for `CoverageExperiment` objects
#' 
#' @rdname show
#' 
#' @param object a CoverageExperiment object
#' @return `tibble`
#' 
#' @export

setMethod("show", signature("CoverageExperiment"), function(object) {
    w <- width(rowRanges(object)[[1]][1])
    b <- w / ncol(assay(object, "coverage")[1, 1][[1]])
    if (
        isTRUE(x = getOption(x = "restore_SummarizedExperiment_show", default = FALSE)) | 
        isTRUE(x = getOption(x = "restore_CoverageExperiment_show", default = FALSE)) | 
        isFALSE("tidySummarizedExperiment" %in% .packages())
    ) {
        f <- getMethod(f = "show", signature = "SummarizedExperiment", 
            where = asNamespace(ns = "SummarizedExperiment"))
        f(object = object)
        cat(paste0("width: ", w, '\n'))
        cat(paste0("binning: ", b, '\n'))
    }
    else {
        print(object)
    }
})
