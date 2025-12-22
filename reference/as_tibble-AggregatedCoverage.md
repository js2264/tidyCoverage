# as_tibble

Coerce an `CoverageExperiment` or `AggregatedCoverage` object into a
`tibble`

## Usage

``` r
# S3 method for class 'AggregatedCoverage'
as_tibble(x, ...)
```

## Arguments

- x:

  A data frame, list, matrix, or other object that could reasonably be
  coerced to a tibble.

- ...:

  Unused, for extensibility.

## Value

`tibble`

## Row names

The default behavior is to silently remove row names.

New code should explicitly convert row names to a new column using the
`rownames` argument.

For existing code that relies on the retention of row names, call
`pkgconfig::set_config("tibble::rownames" = NA)` in your script or in
your package's [`.onLoad()`](https://rdrr.io/r/base/ns-hooks.html)
function.

## Life cycle

Using
[`as_tibble()`](https://tibble.tidyverse.org/reference/as_tibble.html)
for vectors is superseded as of version 3.0.0, prefer the more
expressive `as_tibble_row()` and `as_tibble_col()` variants for new
code.

## See also

[`tibble()`](https://tibble.tidyverse.org/reference/tibble.html)
constructs a tibble from individual columns.
[`enframe()`](https://tibble.tidyverse.org/reference/enframe.html)
converts a named vector to a tibble with a column of names and column of
values. Name repair is implemented using
[`vctrs::vec_as_names()`](https://vctrs.r-lib.org/reference/vec_as_names.html).

## Examples

``` r
data(ac)
as_tibble(ac)
#> # A tibble: 6,000 × 13
#>    .sample .feature track features coord    mean median   min   max    sd     se
#>    <chr>   <fct>    <chr> <fct>    <dbl>   <dbl>  <dbl> <dbl> <dbl> <dbl>  <dbl>
#>  1 RNA_fwd Scc1     RNA_… Scc1     -1500 -0.0921 -0.611 -3.25  10.5  1.42 0.0572
#>  2 RNA_fwd Scc1     RNA_… Scc1     -1499 -0.0915 -0.609 -3.25  10.5  1.42 0.0572
#>  3 RNA_fwd Scc1     RNA_… Scc1     -1498 -0.0898 -0.609 -3.25  10.5  1.42 0.0573
#>  4 RNA_fwd Scc1     RNA_… Scc1     -1497 -0.0914 -0.609 -3.25  10.5  1.42 0.0573
#>  5 RNA_fwd Scc1     RNA_… Scc1     -1496 -0.0915 -0.609 -3.25  10.5  1.42 0.0573
#>  6 RNA_fwd Scc1     RNA_… Scc1     -1495 -0.0912 -0.609 -3.25  10.5  1.42 0.0572
#>  7 RNA_fwd Scc1     RNA_… Scc1     -1494 -0.0912 -0.609 -3.25  10.5  1.42 0.0572
#>  8 RNA_fwd Scc1     RNA_… Scc1     -1493 -0.0915 -0.609 -3.25  10.5  1.42 0.0572
#>  9 RNA_fwd Scc1     RNA_… Scc1     -1492 -0.0907 -0.609 -3.25  10.5  1.42 0.0573
#> 10 RNA_fwd Scc1     RNA_… Scc1     -1491 -0.0903 -0.609 -3.25  10.5  1.42 0.0572
#> # ℹ 5,990 more rows
#> # ℹ 2 more variables: ci_low <dbl>, ci_high <dbl>
```
