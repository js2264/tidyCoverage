#' as_tibble
#' 
#' Coerce an `CoverageExperiment` object into a `tibble`
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

as_tibble.AggregatedCoverage <- function(x, ...) {
    tracks <- colData(x)$track
    features <- rowData(x)$features
    assays <- names(assays(x))
    w <- width(rowRanges(x)[[1]])[[1]]
    bin <- w / length(assay(x, "mean")[1, 1][[1]])
    lapply(tracks, function(t) {
        lapply(features, function(f) {
            lapply(assays(x)[assays], `[[`, f, t) |> 
                stats::setNames(assays) |> 
                dplyr::bind_cols() |> 
                dplyr::mutate(
                    coord = seq(-w/2, w/2-1, by = bin), 
                    track = t, features = f
                ) |> 
                dplyr::relocate(coord) |> 
                dplyr::relocate(features) |> 
                dplyr::relocate(track)
        }) |> dplyr::bind_rows()
    }) |> 
        dplyr::bind_rows() |> 
        dplyr::left_join(colData(x) |> as.data.frame(), by = 'track')
}
