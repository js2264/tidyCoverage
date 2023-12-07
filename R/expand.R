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
#' @name expand
#' @aliases expand,CoverageExperiment-method
#' @rdname expand
#' 
#' @param x a `CoverageExperiment` object
#' @param ... ignored
#' @return a `tibble` object
#' 
#' @importFrom S4Vectors expand
#' @export
#' @examples 
#' data(ce)
#' ce
#' 
#' expand(ce)

setMethod("expand", signature(x = "CoverageExperiment"), function(x, ...) {
    tracks <- colData(x)$track
    features <- rowData(x)$features
    w <- width(rowRanges(x)[[1]])[[1]]
    bin <- w / ncol(assay(x, "coverage")[1, 1][[1]])
    df <- lapply(features, function(f) {
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
                dplyr::relocate(tidyr::all_of(c("track", "features", "chr", "ranges", "strand")))
            d <- tidyr::pivot_longer(
                m, 
                !tidyr::any_of(c("track", "features", "chr", "ranges", "strand")), 
                names_to = "coord", values_to = "coverage"
            )
            d$coord <- coord
            d
        }) |> dplyr::bind_rows()
    }) |> 
        dplyr::bind_rows() |> 
        dplyr::left_join(colData(x) |> as.data.frame(), by = 'track') |> 
        dplyr::group_by(track, features, ranges)
    return(df)
})
