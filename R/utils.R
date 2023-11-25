.set_seqinfo <- function(tracks) {
    n <- names(tracks)
    sls <- lapply(tracks, lengths)
    has_seqlengths <- any(!is.na(lapply(sls, sum)))
    if (!has_seqlengths) {
        sis <- lapply(tracks, function(rle) {
            Seqinfo(
                names(lengths(rle)), 
                lengths(rle)
            )
        }) |> unique()
        if (length(unique(sis)) > 1) 
            stop("More than 1 seqinfo inferred from the tracks.")
        si <- sis[[1]]
    } 
    else {
        sl <- sls[[which(!is.na(lapply(sls, sum)))[[1]]]]
        si <- list(Seqinfo(
            names(sl), 
            sl
        ))
    }
    tracks <- lapply(tracks, function(rle) {
        GenomeInfoDb::seqinfo(rle) <- si
        return(rle)
    })
    names(tracks) <- n
    return(tracks)
}

.set_seqinfo_bwfl <- function(tracks) {
    n <- names(tracks)
    sis <- lapply(tracks, seqinfo) |> unique()
    if (length(sis) > 1) 
        stop("More than 1 seqinfo inferred from the tracks.")
    si <- sis[[1]]
    return(tracks)
} 

.resize_granges <- function(gr, width, seqinfo) {
    GenomeInfoDb::seqlevels(gr, pruning.mode = "coarse") <- GenomeInfoDb::seqlevels(seqinfo)
    GenomeInfoDb::seqinfo(gr) <- seqinfo
    gr <- suppressWarnings(resize(gr, fix = 'center', width = width))
    w <- width
    gr <- trim(gr)
    gr <- gr[width(gr) == w]
}

#' @importFrom IRanges NumericList
#' @importFrom stats qt

.compute_cov <- function(rle, gr, center, scale, ignore.strand = TRUE) {
    scores <- rle[gr]
    scores <- IRanges::NumericList(scores)
    scores <- as.matrix(scores)
    scores <- t(scale(t(scores), center = center, scale = scale))
    scores
}

.compute_cov_bw <- function(bwf, gr, center, scale, ignore.strand = TRUE) {
    scores <- rtracklayer::import(
        bwf, 
        selection = rtracklayer::BigWigSelection(ranges = gr), 
        as = "NumericList", 
        format = "bigWig"
    )
    scores <- as.matrix(scores)
    scores <- t(scale(t(scores), center = center, scale = scale))
    scores
}

.summarize_cov <- function(scores) {
    mean <- colMeans(scores, na.rm = TRUE)
    median <- apply(scores, 2, median, na.rm = TRUE)
    min <- apply(scores, 2, min, na.rm = TRUE)
    max <- apply(scores, 2, max, na.rm = TRUE)
    sd <- apply(scores, 2, sd, na.rm = TRUE)
    se <- sd/sqrt(nrow(scores))
    ci_low <- mean - stats::qt(1 - (0.05 / 2), nrow(scores) - 1) * se
    ci_high <- mean + stats::qt(1 - (0.05 / 2), nrow(scores) - 1) * se
    data.frame(
        coord = seq(-ncol(scores)/2, ncol(scores)/2-1, 1), 
        mean = mean, median = median, min = min, max = max, 
        sd = sd, se = se, ci_low = ci_low, ci_high = ci_high
    )
}

.coarsen <- function(x, bin, FUN, ...) {
    stats::aggregate(
        x, 
        by = list(rep(seq(1, length(x)/bin), each = bin)), 
        FUN = FUN, ...
    )$x
}
