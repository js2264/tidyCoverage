#' Expand a CoverageExperiment object
#' 
#' @description
#' A `CoverageExperiment` object can be coerced into a `tibble` using the 
#' `tidySummarizedExperiment` package, but this will not turn 
#' each coverage matrix into a "long" format. The `expand` function 
#' provided here allows one to coerce a `CoverageExperiment`
#' object into a long data frame, and adds the `ranges` 
#' and `seqnames` to the resulting `tibble`. 
#' 
#' @name expand,CoverageExperiment
#' @aliases expand,CoverageExperiment-method
#' @rdname expand
#' 
#' @param x a `CoverageExperiment` object
#' @param ...,.name_repair ignored
#' @return a `tibble` object
#' 
#' @importFrom tidyr expand
#' @export
#' @examples 
#' data(ce)
#' ce
#' 
#' expand(ce)

expand.CoverageExperiment <- function(x, ..., .name_repair = NULL) {
    tracks <- colData(x)$track
    features <- rowData(x)$features
    w <- width(rowRanges(x)[[1]])[[1]]
    bin <- w / ncol(assay(x, "coverage")[1, 1][[1]])
    df <- lapply(features, function(f) {
        rr <- rowRanges(x)[[f]]
        rrdf <- as.data.frame(rr)
        coord <- lapply(seq_len(nrow(rrdf)), function(K) seq(rrdf[K, 'start'], rrdf[K, 'end'], by = bin)) |> unlist()
        coord.scaled <- lapply(seq_len(nrow(rrdf)), function(K) seq(-w/2, w/2-1, by = bin)) |> unlist()
        lapply(tracks, function(t) {
            m <- assay(x, "coverage")[f, t][[1]] |> 
                as.data.frame() |>
                dplyr::mutate(
                    track = t, features = f, 
                    chr = as.vector(seqnames(rr)), 
                    ranges = as.character(rr), 
                    strand = as.vector(strand(rr))
                ) |> 
                dplyr::relocate(tidyr::all_of(c("track", "features", "chr", "ranges", "strand")))
            d <- tidyr::pivot_longer(
                m, 
                !tidyr::any_of(c("track", "features", "chr", "ranges", "strand", "coord.scaled")), 
                names_to = "coord", values_to = "coverage"
            )
            d$coord <- coord
            d$coord.scaled <- coord.scaled
            d
        }) |> dplyr::bind_rows()
    }) |> 
        dplyr::bind_rows() |> 
        dplyr::left_join(colData(x) |> as.data.frame(), by = 'track') |> 
        dplyr::group_by(track, features, ranges)
    return(df)
}
