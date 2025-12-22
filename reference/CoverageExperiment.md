# CoverageExperiment

`CoverageExperiment` objects store coverages for individual tracks over
different sets of features. The `coverage` assay contains a separate
matrix for each combination of track x features. `CoverageExperiment`
objects are instantiated using the `CoverageExperiment()` function, and
can be coarsened using the `coarsen()` function.

## Usage

``` r
CoverageExperiment(tracks, features, ...)

coarsen(x, window, ...)

# S4 method for class 'BigWigFileList,GRangesList'
CoverageExperiment(
  tracks,
  features,
  width = NULL,
  center = FALSE,
  scale = FALSE,
  ignore.strand = TRUE,
  window = 1,
  BPPARAM = BiocParallel::bpparam()
)

# S4 method for class 'BigWigFileList,GRanges'
CoverageExperiment(tracks, features, ...)

# S4 method for class 'BigWigFileList,list'
CoverageExperiment(tracks, features, ...)

# S4 method for class 'BigWigFile,GRangesList'
CoverageExperiment(tracks, features, ...)

# S4 method for class 'BigWigFile,GRanges'
CoverageExperiment(tracks, features, ...)

# S4 method for class 'BigWigFile,list'
CoverageExperiment(tracks, features, ...)

# S4 method for class 'list,GRangesList'
CoverageExperiment(
  tracks,
  features,
  width = NULL,
  center = FALSE,
  scale = FALSE,
  ignore.strand = TRUE,
  window = 1,
  BPPARAM = BiocParallel::bpparam()
)

# S4 method for class 'list,GRanges'
CoverageExperiment(tracks, features, ...)

# S4 method for class 'list,list'
CoverageExperiment(tracks, features, ...)

# S4 method for class 'RleList,GRangesList'
CoverageExperiment(tracks, features, ...)

# S4 method for class 'RleList,GRanges'
CoverageExperiment(tracks, features, ...)

# S4 method for class 'RleList,list'
CoverageExperiment(tracks, features, ...)

# S4 method for class 'CoverageExperiment'
coarsen(x, window = 1, BPPARAM = BiocParallel::bpparam())
```

## Arguments

- tracks:

  A genomic track imported as a `RleList` or a *named* list of genomic
  tracks.

- features:

  A set of features imported as `GRanges` or a *named* `GRangesList`.

- ...:

  Passed to the relevant method

- x:

  a `CoverageExperiment` object

- window:

  an integer to coarsen coverage by.

- width:

  Width to resize each set of genomic features

- scale, center:

  Logical, whether to scale and/or center tracks prior to summarization

- ignore.strand:

  Logical, whether to not take the features strand information

- BPPARAM:

  Passed to BiocParallel.

## Value

A `CoverageExperiment` object

## Examples

``` r
library(rtracklayer)
library(purrr)
#> 
#> Attaching package: ‘purrr’
#> The following object is masked from ‘package:GenomicRanges’:
#> 
#>     reduce
#> The following object is masked from ‘package:IRanges’:
#> 
#>     reduce
library(plyranges)
#> 
#> Attaching package: ‘plyranges’
#> The following objects are masked from ‘package:dplyr’:
#> 
#>     between, n, n_distinct
TSSs_bed <- system.file("extdata", "TSSs.bed", package = "tidyCoverage")
features <- import(TSSs_bed) |> filter(strand == '+')

#############################################################################
## 1. Creating a `CoverageExperiment` object from a single BigWigFile
#############################################################################

RNA_fwd <- system.file("extdata", "RNA.fwd.bw", package = "tidyCoverage")
tracks <- BigWigFile(RNA_fwd)
CoverageExperiment(tracks, features, width = 5000)
#> class: CoverageExperiment 
#> dim: 1 1 
#> metadata(0):
#> assays(1): coverage
#> rownames(1): features
#> rowData names(2): features n
#> colnames(1): track
#> colData names(1): track
#> width: 5000

#############################################################################
## 2. Creating a `CoverageExperiment` object from a BigWigFileList
#############################################################################

RNA_rev <- system.file("extdata", "RNA.rev.bw", package = "tidyCoverage")
tracks <- BigWigFileList(list(RNA_fwd = RNA_fwd, RNA_rev = RNA_rev))
CoverageExperiment(tracks, features, width = 5000)
#> class: CoverageExperiment 
#> dim: 1 2 
#> metadata(0):
#> assays(1): coverage
#> rownames(1): features
#> rowData names(2): features n
#> colnames(2): RNA_fwd RNA_rev
#> colData names(1): track
#> width: 5000

#############################################################################
## 3. Creating a `CoverageExperiment` object from imported bigwig files
#############################################################################

tracks <- list(
    RNA_fwd = system.file("extdata", "RNA.fwd.bw", package = "tidyCoverage"),
    RNA_rev = system.file("extdata", "RNA.rev.bw", package = "tidyCoverage")
) |> map(import, as = 'Rle')
CoverageExperiment(tracks, features, width = 5000)
#> class: CoverageExperiment 
#> dim: 1 2 
#> metadata(0):
#> assays(1): coverage
#> rownames(1): features
#> rowData names(2): features n
#> colnames(2): RNA_fwd RNA_rev
#> colData names(1): track
#> width: 5000

#############################################################################
## 4. Correct for strandness when recovering coverage
#############################################################################

TSSs_bed <- system.file("extdata", "TSSs.bed", package = "tidyCoverage")
features <- list(
    TSS_fwd = import(TSSs_bed) |> filter(strand == '+'), 
    TSS_rev = import(TSSs_bed) |> filter(strand == '-')
)
tracks <- list(
    RNA_fwd = system.file("extdata", "RNA.fwd.bw", package = "tidyCoverage"),
    RNA_rev = system.file("extdata", "RNA.rev.bw", package = "tidyCoverage")
) |> map(import, as = 'Rle')
CoverageExperiment(tracks, features, width = 5000, ignore.strand = FALSE)
#> class: CoverageExperiment 
#> dim: 2 2 
#> metadata(0):
#> assays(1): coverage
#> rownames(2): TSS_fwd TSS_rev
#> rowData names(2): features n
#> colnames(2): RNA_fwd RNA_rev
#> colData names(1): track
#> width: 5000

#############################################################################
## Aggregating a `CoverageExperiment` object
#############################################################################
data(ce)
coarsen(ce, window = 10)
#> class: CoverageExperiment 
#> dim: 1 2 
#> metadata(0):
#> assays(1): coverage
#> rownames(1): Scc1
#> rowData names(2): features n
#> colnames(2): RNA_fwd RNA_rev
#> colData names(1): track
#> width: 3000
```
