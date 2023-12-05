#' Example `CoverageExperiment` and `AggregatedCoverage` objects
#'
#' Two example objects are provided in the `tidyCoverage` package: 
#' - `ce`: a `CoverageExperiment` dataset containing stranded RNA-seq coverage 
#' (forward and reverse) over Scc1 peaks (± 1kb).
#' - `ac`: an `AggregatedCoverage` object obtained with `aggregate(ce)`. 
#' 
#' Data was generated in yeast (S288c) and aligned to reference R64-1-1. 
#'
#' @name data
#' @rdname data
#' @format `CoverageExperiment` object containing 1 features set and 2 tracks.
#'
#' @usage data(ce)
"ce"

#' @rdname data
#' @format `AggregatedCoverage` object containing 1 features set and 2 tracks.
#' @usage data(ac)
"ac"
