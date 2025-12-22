# Example `CoverageExperiment` and `AggregatedCoverage` objects

Two example objects are provided in the `tidyCoverage` package:

- `ce`: a `CoverageExperiment` dataset containing stranded RNA-seq
  coverage (forward and reverse) over Scc1 peaks (± 1kb).

- `ac`: an `AggregatedCoverage` object obtained with `aggregate(ce)`.

## Usage

``` r
data(ce)

data(ac)
```

## Format

`CoverageExperiment` object containing 1 features set and 2 tracks.

`AggregatedCoverage` object containing 1 features set and 2 tracks.

## Details

Data was generated in yeast (S288c) and aligned to reference R64-1-1.
