methods::setClass(
    "CoverageExperiment", 
    contains = c("RangedSummarizedExperiment")
)
methods::setClass(
    "AggregatedCoverage", 
    contains = c("RangedSummarizedExperiment")
)
