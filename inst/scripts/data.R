## ce
features <- list(
    Scc1 = system.file("extdata", "TSSs.bed", package = "tidyCoverage")
) |> map(import) |> map(filter, strand == '+') |> map(`[`, 1:1000)
tracks <- list(
    RNA_fwd = system.file("extdata", "RNA.fwd.bw", package = "tidyCoverage"),
    RNA_rev = system.file("extdata", "RNA.rev.bw", package = "tidyCoverage")
) |> map(import, as = 'Rle')
ce <- CoverageExperiment(tracks, features, width = 3000, scale = TRUE, center = TRUE)
usethis::use_data(ce, overwrite = TRUE)

## ac
ac <- aggregate(ce)
usethis::use_data(ac)
