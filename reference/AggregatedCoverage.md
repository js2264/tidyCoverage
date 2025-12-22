# aggregate

Bin coverage contained in a `CoverageExperiment` into an
`AggregatedCoverage` object.

## Usage

``` r
# S4 method for class 'CoverageExperiment'
aggregate(x, bin = 1, ...)
```

## Arguments

- x:

  a `CoverageExperiment` object

- bin:

  an integer to bin each assay by. The `width` of the
  `AggregatedCoverage` object should be a multiple of `bin`.

- ...:

  ignored

## Value

an `AggregatedCoverage` object

## Examples

``` r
data(ce)
aggregate(ce, bin = 10)
#> class: AggregatedCoverage 
#> dim: 1 2 
#> metadata(0):
#> assays(8): mean median ... ci_low ci_high
#> rownames(1): Scc1
#> rowData names(2): features n
#> colnames(2): RNA_fwd RNA_rev
#> colData names(1): track
#> width: 3000
#> binning: 100
```
