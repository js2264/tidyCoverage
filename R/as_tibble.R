#' as_tibble
#' 
#' Coerce an `CoverageExperiment` object into a `tibble`
#' 
#' @name as_tibble-methods
#' @rdname as_tibble-methods
#' @param x an `CoverageExperiment` object
#' @param ... ignored
#' @return `tibble`
#' 
#' @importFrom dplyr as_tibble
#' @examples 
#' data(ac)
#' as_tibble(ac)
NULL

#' @rdname as_tibble-methods
#' @export

as_tibble.CoverageExperiment <- function(x, ...) {
    tracks <- colData(x)$track
    features <- rowData(x)$features
    w <- width(rowRanges(x)[[1]])[[1]]
    bin <- w / ncol(assay(x, "coverage")[1, 1][[1]])
    lapply(features, function(f) {
        rr <- rowRanges(x)[[f]]
        rrdf <- as.data.frame(rr)
        coord <- lapply(seq_len(nrow(rrdf)), function(K) seq(rrdf[K, 'start'], rrdf[K, 'end'], by = bin)) |> unlist()
        lapply(tracks, function(t) {
            m <- assay(x, "coverage")[f, t][[1]] |> 
                as.data.frame() |>
                dplyr::mutate(
                    track = t, features = f, 
                    chr = as.vector(seqnames(rr)), 
                    ranges = as.character(rr), 
                    strand = as.vector(strand(rr))
                ) |> 
                dplyr::relocate(all_of(c("track", "features", "chr", "ranges", "strand")))
            d <- tidyr::pivot_longer(
                m, 
                !any_of(c("track", "features", "chr", "ranges", "strand")), 
                names_to = "coord", values_to = "coverage"
            )
            d$coord <- coord
            d
        }) |> dplyr::bind_rows()
    }) |> 
        dplyr::bind_rows() |> 
        dplyr::left_join(colData(x) |> as.data.frame(), by = 'track') |> 
        dplyr::group_by(track, features, ranges)
}

#' @rdname as_tibble-methods
#' @export

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