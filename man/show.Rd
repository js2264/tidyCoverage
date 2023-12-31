% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/show.R
\name{show}
\alias{show}
\alias{show,CoverageExperiment-method}
\alias{show,AggregatedCoverage-method}
\alias{print.CoverageExperiment}
\alias{print.AggregatedCoverage}
\alias{tbl_format_header.tidyCoverageExperiment}
\alias{tbl_format_header.tidyAggregatedCoverage}
\title{show}
\usage{
\S4method{show}{CoverageExperiment}(object)

\S4method{show}{AggregatedCoverage}(object)

\method{print}{CoverageExperiment}(x, ..., n = NULL)

\method{print}{AggregatedCoverage}(x, ..., n = NULL)

\method{tbl_format_header}{tidyCoverageExperiment}(x, setup, ...)

\method{tbl_format_header}{tidyAggregatedCoverage}(x, setup, ...)
}
\arguments{
\item{object}{a \code{CoverageExperiment} or \code{AggregatedCoverage} object}

\item{x}{Object to format or print.}

\item{...}{Passed on to \code{\link[=tbl_format_setup]{tbl_format_setup()}}.}

\item{n}{Number of rows to show. If \code{NULL}, the default, will print all rows
if less than the \code{print_max} \link[pillar:pillar_options]{option}.
Otherwise, will print as many rows as specified by the
\code{print_min} \link[pillar:pillar_options]{option}.}

\item{setup}{a setup object returned from \code{\link[pillar:tbl_format_setup]{pillar::tbl_format_setup()}}.}
}
\value{
\verb{Prints a message to the console describing the contents of the }CoverageExperiment\code{or}AggregatedCoverage` objects.
}
\description{
show method for \code{CoverageExperiment} and \code{AggregatedCoverage} objects
}
\examples{
data(ce)
print(ce)
data(ac)
print(ac)
}
