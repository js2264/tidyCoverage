#' CoverageExperiment
#'
#' This function initiates a `CoverageExperiment` object. 
#'
#' @name CoverageExperiment
#' @rdname CoverageExperiment
#' 
#' @param tracks A genomic track imported as a `RleList` or a *named* list of 
#' genomic tracks.
#' @param features A set of features imported as `GRanges` or a *named* 
#' `GRangesList`. 
#' @param width Width to resize each set of genomic features
#' @param ignore.strand Logical, whether to not take the features strand
#' information
#' @param scale,center Logical, whether to scale and/or center tracks prior to 
#' summarization
#' @param BPPARAM Passed to BiocParallel
#' @return An `CoverageExperiment` object
#'
#' @import GenomicRanges
#' @import GenomeInfoDb
#' @import SummarizedExperiment
#' @import methods
#' @importFrom dplyr filter
#' @importFrom dplyr bind_rows
#' @importFrom dplyr bind_cols
#' @importFrom dplyr mutate
#' @importFrom dplyr relocate
#' @importFrom BiocParallel bpparam
#' @importFrom BiocParallel bplapply
#' @include AllGenerics.R
#'
#' @examples
#' library(rtracklayer)
#' library(purrr)
#' library(plyranges)
#' features <- system.file("extdata", "TSSs.bed", package = "CoverageExperiment") |> 
#'     import() |> filter(strand == '+') |> GRangesList(x = _)
#' tracks <- list(
#'     RNA_fwd = system.file("extdata", "SRR2045244.fwd.bw", package = "CoverageExperiment"),
#'     RNA_rev = system.file("extdata", "SRR2045244.rev.bw", package = "CoverageExperiment")
#' ) |> map(import, as = 'Rle')
#' CE <- CoverageExperiment(tracks, features, width = 5000)
#' AC
NULL

#' @rdname CoverageExperiment
#' @export

setMethod(
    "CoverageExperiment", 
    signature(tracks = "list", features = "GRangesList"), 
    function(
        tracks, features, width, 
        center = FALSE, scale = FALSE, 
        ignore.strand = TRUE, 
        BPPARAM = BiocParallel::bpparam()
    ) {
        ## Check that input args are valid
        stopifnot(length(unique(names(tracks))) == length(tracks)) # named tracks list
        stopifnot(length(unique(names(features))) == length(features)) # named features GRangesList

        ## Extend and filter features
        tracks <- .set_seqinfo(tracks)
        si <- seqinfo(tracks[[1]])[[1]]
        features <- lapply(features, .resize_granges, width = width, seqinfo = si) 

        ## Prepare cData and rData
        cData <- data.frame(
            track = names(tracks)
        )
        rData <- GRangesList(features)
        mcols(rData) <- data.frame(
            features = names(features)
        )

        ## Extract coverage scores
        combs <- expand.grid(names(tracks), names(features)) |> 
            setNames(c("tracks", "features"))
        l <- BiocParallel::bplapply(seq_len(nrow(combs)), function(K) {
            t <- combs[K, "tracks"]
            f <- combs[K, "features"]
            m <- .compute_cov(
                tracks[[t]], features[[f]], 
                scale = scale, center = center, ignore.strand = ignore.strand
            )
            return(m)
        }, BPPARAM = BPPARAM)
        names(l) <- paste(combs$tracks, combs$features, sep = '^')

        ## Fill out different assay matrices
        m0 <- matrix(
            list(), 
            nrow = length(features), ncol = length(tracks)
        )
        colnames(m0) <- names(tracks)
        rownames(m0) <- names(features)
        assays <- "coverage"
        l_assays <- vector("list", length = length(assays))
        names(l_assays) <- assays
        for (assay in assays) {
            m <- m0
            for (t in names(tracks)) {
                for (f in names(features)) {
                    name <- paste(t, f, sep = '^')
                    m[f, t][[1]] <- l[[name]]
                }
            }
            l_assays[[assay]] <- m
        }

        ## Instantiate and return the CoverageExperiment final object
        CE <- methods::new(
            "CoverageExperiment",
            SummarizedExperiment::SummarizedExperiment(
                rowRanges = rData,
                colData = cData,
                assays = l_assays
            )
        )
        return(CE)
    }
)


#' @rdname CoverageExperiment
#' @export

setMethod(
    "CoverageExperiment", 
    signature(tracks = "list", features = "list"), 
    function(
        tracks, features, width, 
        center = FALSE, scale = FALSE, 
        ignore.strand = TRUE, 
        BPPARAM = BiocParallel::bpparam()
    ) {
        features <- as(features, 'GRangesList')
        CoverageExperiment(
            tracks, features, width, 
            center, scale, 
            ignore.strand, 
            BPPARAM
        )
    }
)

#' @rdname CoverageExperiment
#' @export

setMethod(
    "CoverageExperiment", 
    signature(tracks = "RleList", features = "GRangesList"), 
    function(
        tracks, features, width, 
        center = FALSE, scale = FALSE, 
        ignore.strand = TRUE, 
        BPPARAM = BiocParallel::bpparam()
    ) {
        tracks <- list(track = tracks)
        CoverageExperiment(
            tracks, features, width, 
            center, scale, 
            ignore.strand, 
            BPPARAM
        )
    }
)

#' @rdname CoverageExperiment
#' @export

setMethod(
    "CoverageExperiment", 
    signature(tracks = "list", features = "GRanges"), 
    function(
        tracks, features, width, 
        center = FALSE, scale = FALSE, 
        ignore.strand = TRUE, 
        BPPARAM = BiocParallel::bpparam()
    ) {
        features <- GRangesList(features = features)
        CoverageExperiment(
            tracks, features, width, 
            center, scale, 
            ignore.strand, 
            BPPARAM
        )
    }
)

#' @rdname CoverageExperiment
#' @export

setMethod(
    "CoverageExperiment", 
    signature(tracks = "RleList", features = "GRanges"), 
    function(
        tracks, features, width, 
        center = FALSE, scale = FALSE, 
        ignore.strand = TRUE, 
        BPPARAM = BiocParallel::bpparam()
    ) {
        tracks <- list(track = tracks)
        features <- GRangesList(features = features)
        CoverageExperiment(
            tracks, features, width, 
            center, scale, 
            ignore.strand, 
            BPPARAM
        )
    }
)