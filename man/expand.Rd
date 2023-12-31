% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/expand.R
\name{expand,CoverageExperiment}
\alias{expand,CoverageExperiment}
\alias{expand.CoverageExperiment}
\alias{expand,CoverageExperiment-method}
\title{Expand a CoverageExperiment object}
\usage{
\method{expand}{CoverageExperiment}(data, ..., .name_repair = NULL)
}
\arguments{
\item{data}{a \code{CoverageExperiment} object}

\item{..., .name_repair}{ignored}
}
\value{
a \code{tibble} object
}
\description{
A \code{CoverageExperiment} object can be coerced into a \code{tibble} using the
\code{tidySummarizedExperiment} package, but this will not turn
each coverage matrix into a "long" format. The \code{expand} function
provided here allows one to coerce a \code{CoverageExperiment}
object into a long data frame, and adds the \code{ranges}
and \code{seqnames} to the resulting \code{tibble}.
}
\examples{
data(ce)
ce

expand(ce)
}
