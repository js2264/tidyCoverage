# Expand a CoverageExperiment object

A `CoverageExperiment` object can be coerced into a `tibble` using the
`tidySummarizedExperiment` package, but this will not turn each coverage
matrix into a "long" format. The `expand` function provided here allows
one to coerce a `CoverageExperiment` object into a long data frame, and
adds the `ranges` and `seqnames` to the resulting `tibble`.

## Usage

``` r
# S3 method for class 'CoverageExperiment'
expand(data, ..., .name_repair = NULL)
```

## Arguments

- data:

  A data frame.

- ...:

  \<[`data-masking`](https://tidyr.tidyverse.org/reference/tidyr_data_masking.html)\>
  Specification of columns to expand or complete. Columns can be atomic
  vectors or lists.

  - To find all unique combinations of `x`, `y` and `z`, including those
    not present in the data, supply each variable as a separate
    argument: `expand(df, x, y, z)` or `complete(df, x, y, z)`.

  - To find only the combinations that occur in the data, use `nesting`:
    `expand(df, nesting(x, y, z))`.

  - You can combine the two forms. For example,
    `expand(df, nesting(school_id, student_id), date)` would produce a
    row for each present school-student combination for all possible
    dates.

  When used with factors,
  [`expand()`](https://tidyr.tidyverse.org/reference/expand.html) and
  [`complete()`](https://tidyr.tidyverse.org/reference/complete.html)
  use the full set of levels, not just those that appear in the data. If
  you want to use only the values seen in the data, use
  `forcats::fct_drop()`.

  When used with continuous variables, you may need to fill in values
  that do not appear in the data: to do so use expressions like
  `year = 2010:2020` or `year = full_seq(year,1)`.

- .name_repair:

  One of `"check_unique"`, `"unique"`, `"universal"`, `"minimal"`,
  `"unique_quiet"`, or `"universal_quiet"`. See
  [`vec_as_names()`](https://vctrs.r-lib.org/reference/vec_as_names.html)
  for the meaning of these options.

## Value

a `tibble` object

## Grouped data frames

With grouped data frames created by
[`dplyr::group_by()`](https://dplyr.tidyverse.org/reference/group_by.html),
`expand()` operates *within* each group. Because of this, you cannot
expand on a grouping column.

## See also

[`complete()`](https://tidyr.tidyverse.org/reference/complete.html) to
expand list objects.
[`expand_grid()`](https://tidyr.tidyverse.org/reference/expand_grid.html)
to input vectors rather than a data frame.

## Examples

``` r
data(ce)
ce
#> class: CoverageExperiment 
#> dim: 1 2 
#> metadata(0):
#> assays(1): coverage
#> rownames(1): Scc1
#> rowData names(2): features n
#> colnames(2): RNA_fwd RNA_rev
#> colData names(1): track
#> width: 3000
expand(ce)
#> # A tibble: 368,400 × 8
#> # Groups:   track, features, ranges [1,228]
#>    track   features chr   ranges         strand coord coverage coord.scaled
#>    <chr>   <fct>    <chr> <chr>          <chr>  <dbl>    <dbl>        <dbl>
#>  1 RNA_fwd Scc1     II    II:4290-7289:+ +       4290   -0.257        -1500
#>  2 RNA_fwd Scc1     II    II:4290-7289:+ +       4300   -0.257        -1490
#>  3 RNA_fwd Scc1     II    II:4290-7289:+ +       4310   -0.257        -1480
#>  4 RNA_fwd Scc1     II    II:4290-7289:+ +       4320   -0.257        -1470
#>  5 RNA_fwd Scc1     II    II:4290-7289:+ +       4330   -0.257        -1460
#>  6 RNA_fwd Scc1     II    II:4290-7289:+ +       4340   -0.257        -1450
#>  7 RNA_fwd Scc1     II    II:4290-7289:+ +       4350   -0.257        -1440
#>  8 RNA_fwd Scc1     II    II:4290-7289:+ +       4360   -0.257        -1430
#>  9 RNA_fwd Scc1     II    II:4290-7289:+ +       4370   -0.257        -1420
#> 10 RNA_fwd Scc1     II    II:4290-7289:+ +       4380   -0.257        -1410
#> # ℹ 368,390 more rows
```
