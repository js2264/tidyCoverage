#' show
#' 
#' show method for `CoverageExperiment` and `AggregatedCoverage` objects
#' 
#' @name show
#' @aliases show,CoverageExperiment-method
#' @aliases show,AggregatedCoverage-method
#' 
#' @param object a `CoverageExperiment` or `AggregatedCoverage` object
#' @param setup a setup object returned from [pillar::tbl_format_setup()].
#' @importFrom rlang names2
#' @importFrom pillar align
#' @importFrom pillar get_extent
#' @importFrom pillar style_subtle
#' @importFrom pillar tbl_format_header
#' @importFrom fansi strwrap_ctl
#' @importFrom purrr map_chr
#' @importFrom cli console_width
#' @importFrom cli symbol
#' @importFrom vctrs new_data_frame
#' @importFrom SummarizedExperiment assayNames
#' @importFrom stats setNames
#' @export
#' 
#' @inherit tibble::formatting
#' @return `Prints a message to the console describing
#' the contents of the `CoverageExperiment` or `AggregatedCoverage` objects.
#' @examples
#' data(ce)
#' print(ce)
#' data(ac)
#' print(ac)
NULL

#' @name show
#' @export

setMethod("show", signature("CoverageExperiment"), function(object) {
    w <- width(rowRanges(object)[[1]][1])
    # if (
    #     isTRUE(x = getOption(x = "restore_SummarizedExperiment_show", default = FALSE)) | 
    #     isTRUE(x = getOption(x = "restore_CoverageExperiment_show", default = FALSE)) | 
    #     isFALSE("tidySummarizedExperiment" %in% .packages())
    # ) {
        f <- getMethod(f = "show", signature = "SummarizedExperiment", 
            where = asNamespace(ns = "SummarizedExperiment"))
        f(object = object)
        cat(paste0("width: ", w, '\n'))
    # }
    # else {
    #     print(object)
    # }
})

#' @name show
#' @export

setMethod("show", signature("AggregatedCoverage"), function(object) {
    w <- width(rowRanges(object)[[1]][1])
    b <- w / length(assay(object, "mean")[1, 1][[1]])
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

#' @name show
#' @export

print.AggregatedCoverage <- function (x, ..., n = NULL) {
    .x <- x
    x <- as_tibble(.x) 
    x <- vctrs::new_data_frame(x, class=c("tidyAggregatedCoverage", "tbl"))
    attr(x, "width") <- width(rowRanges(.x)[[1]][1])
    attr(x, "binning") <- attr(x, "width") / length(assay(.x, 1)[1, 1][[1]])
    attr(x, "number_of_features") <- nrow(.x)
    attr(x, "number_of_tracks") <- ncol(.x)
    attr(x, "assay_names") <- names(assays(.x))
    attr(x, "named_header") <- sprintf(
        "%s %s %s", 
        nrow(x),
        cli::symbol$times,
        ncol(x)
    ) |>
    setNames("An AggregatedCoverage-tibble abstraction")
    print(x)
    invisible(x) 
}

#' @name show
#' @export

tbl_format_header.tidyAggregatedCoverage <- function(x, setup, ...) {
    width <- x |> attr("width")
    binning <- x |> attr("binning")
    number_of_features <- x |> attr("number_of_features")
    number_of_tracks <- x |> attr("number_of_tracks")
    named_header <- x |> attr("named_header")
    assay_names <- x |> attr("assay_names")
    if (all(rlang::names2(named_header) == "")) {
        header <- named_header
    } else {
        header <-
            paste0(
                align(paste0(rlang::names2(named_header), ":"), space="\U00A0"),
                " ",
                named_header
            ) |>
            # Add further info
            append(sprintf(
                "\033[90m features=%s | tracks=%s | assays=%s\033[39m",
                number_of_features,
                number_of_tracks,
                assay_names |> paste(collapse=", ")
            ), after = 1) |> 
            # Add further info re: width/binning
            append(sprintf(
                "\033[90m width=%s | binning=%s\033[39m",
                width,
                binning
            ))
    }
    pillar::style_subtle(.pillar___format_comment(header, width=setup$width))
}

.pillar___format_comment <- function (x, width) {
    if (length(x) == 0L) {
        return(character())
    }
    purrr::map_chr(x, .pillar___wrap, prefix="# ",
        width=min(width, cli::console_width()))
}

.pillar___strwrap2 <- function (x, width, indent) {
    fansi::strwrap_ctl(x, width=max(width, 0), indent=indent,
        exdent=indent + 2)
}

.pillar___wrap <- function (..., indent=0, prefix="", width) {
    x <- paste0(..., collapse="")
    wrapped <- .pillar___strwrap2(x, width - pillar::get_extent(prefix), indent)
    wrapped <- paste0(prefix, wrapped)
    wrapped <- gsub("\U00A0", " ", wrapped)
    paste0(wrapped, collapse="\n")
}
