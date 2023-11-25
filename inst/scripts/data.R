## ac
features <- list(
    Scc1 = system.file("extdata", "TSSs.bed", package = "AggregatedCoverage"),
    `Conv_transcription` = system.file("extdata", "conv_transcription_loci.bed", package = "AggregatedCoverage")
) |> map(import) |> map(filter, strand == '+') 
tracks <- list(
    RNA_fwd = system.file("extdata", "SRR2045244.fwd.CPM.bw", package = "AggregatedCoverage"),
    RNA_rev = system.file("extdata", "SRR2045244.rev.CPM.bw", package = "AggregatedCoverage"),
) |> map(import, as = 'Rle')
ac <- AggregatedCoverage(tracks, features, width = 3000, scale = TRUE, center = TRUE)
usethis::use_data(ac)
