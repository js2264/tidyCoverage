# show method for `CoverageExperiment` and `AggregatedCoverage` objects

show method for `CoverageExperiment` and `AggregatedCoverage` objects

## Arguments

- object:

  a `CoverageExperiment` or `AggregatedCoverage` object

- setup:

  a setup object returned from
  [`pillar::tbl_format_setup()`](https://pillar.r-lib.org/reference/tbl_format_setup.html).

## Value

`Prints a message to the console describing the contents of the `CoverageExperiment`or`AggregatedCoverage\`
objects.

## Examples

``` r
data(ce)
print(ce)
#> class: CoverageExperiment 
#> dim: 1 2 
#> metadata(0):
#> assays(1): coverage
#> rownames(1): Scc1
#> rowData names(2): features n
#> colnames(2): RNA_fwd RNA_rev
#> colData names(1): track
#> width: 3000
data(ac)
print(ac)
#> class: AggregatedCoverage 
#> dim: 1 2 
#> metadata(0):
#> assays(8): mean median ... ci_low ci_high
#> rownames(1): Scc1
#> rowData names(1): features
#> colnames(2): RNA_fwd RNA_rev
#> colData names(1): track
#> width: 3000
#> binning: 1
```
