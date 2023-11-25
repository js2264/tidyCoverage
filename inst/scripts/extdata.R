## conv_transcription_loci.bed
conv_transcription_loci <- rtracklayer::import('~/genomes/S288c/S288c.gtf') |> 
    filter(type == "transcript") |> 
    sort(ignore.strand = TRUE) |> 
    as_tibble() |> 
    group_by(seqnames) |>
    dplyr::group_modify(
        ~ dplyr::filter(.x, 
            strand == '+' & dplyr::lead(strand, n = 1) == '-' | 
            strand == '-' & dplyr::lag(strand, n = 1) == '+' 
        ) |> 
        dplyr::mutate(start2 = end, end = dplyr::lead(start), start = start2) |> 
        dplyr::filter(strand == '+' & dplyr::lead(strand, n = 1) == '-')
    ) |> 
    dplyr::filter(end > start) |> 
    makeGRangesFromDataFrame(keep.extra.columns = FALSE) |> 
    anchor_center() |> 
    mutate(width = 1)
si <- Seqinfo(genome = "R64")
seqlevelsStyle(conv_transcription_loci) <- 'Ensembl'
seqlevels(conv_transcription_loci) <- seqlevels(si)
seqinfo(conv_transcription_loci) <- si
export(conv_transcription_loci, 'inst/extdata/conv_transcription_loci.bed')

## TSSs.bed
library(plyranges)
txdb <- TxDb.Scerevisiae.UCSC.sacCer3.sgdGene::TxDb.Scerevisiae.UCSC.sacCer3.sgdGene
TSSs <- GenomicFeatures::genes(txdb) |> 
    # filter(strand == '+') |> 
    anchor_5p() |> 
    mutate(width = 1)
seqlevelsStyle(TSSs) <- "Ensembl"
export(TSSs, 'inst/extdata/TSSs.bed')

## Scc1-peaks.narrowPeak <--- /home/rsg/Projects/20220309_Christophe_GC-paper/data/WT/ChIP/peaks/CH224/CH224_vs-CH225_genome-S288c_YMT7BP_peaks.narrowPeak
file.copy(
    "/home/rsg/Projects/20220309_Christophe_GC-paper/data/WT/ChIP/peaks/CH224/CH224_vs-CH225_genome-S288c_YMT7BP_peaks.narrowPeak", 
    "inst/extdata/Scc1-peaks.narrowPeak"
)


## Scc1-vs-input.bw <-------- /home/rsg/Projects/20220309_Christophe_GC-paper/data/WT/ChIP/tracks/CH224/CH224^unmapped_CBS138^mapped_S288c^YMT7BP.vs-CH225.bw
## SRR2045244.fwd.CPM.bw <--- /home/rsg/Projects/20220309_Christophe_GC-paper/data/WT/RNA/tracks/SRR2045244/SRR2045244^mapped_S288c^KEYTQL.fwd.CPM.bw
## SRR2045244.rev.CPM.bw <--- /home/rsg/Projects/20220309_Christophe_GC-paper/data/WT/RNA/tracks/SRR2045244/SRR2045244^mapped_S288c^KEYTQL.rev.CPM.bw
## PolII.bw <---------------- /home/rsg/Projects/20220309_Christophe_GC-paper/data/WT/ChIP/tracks/CH244/CH244^mapped_S288c^CIJXLY.CPM.bw"
## MNase.bw <---------------- ~/Projects/20230517_Lea_MNase-timecourse/nuc_cov_Pneumo-time-course.bw
library(rtracklayer)
library(plyranges)
library(purrr)
t <- import("/home/rsg/Projects/20220309_Christophe_GC-paper/data/WT/ChIP/tracks/CH224/CH224^unmapped_CBS138^mapped_S288c^YMT7BP.vs-CH225.bw", as = 'Rle')
seqlengths(t) <- lengths(t)
t <- t[1:16]
si <- seqinfo(t)
seqlevels(si) <- seqlevels(si)[1:16]
bins <- plyranges::tile_ranges(si |> as("GRanges"), width = 50)
tracks <- list(
    `SRR2045244.fwd.bw` = "/home/rsg/Projects/20220309_Christophe_GC-paper/data/WT/RNA/tracks/SRR2045244/SRR2045244^mapped_S288c^KEYTQL.fwd.CPM.bw", 
    `SRR2045244.rev.bw` = "/home/rsg/Projects/20220309_Christophe_GC-paper/data/WT/RNA/tracks/SRR2045244/SRR2045244^mapped_S288c^KEYTQL.rev.CPM.bw", 
    `Scc1.bw` = "/home/rsg/Projects/20220309_Christophe_GC-paper/data/WT/ChIP/tracks/CH224/CH224^unmapped_CBS138^mapped_S288c^YMT7BP.vs-CH225.bw", 
    `PolII.bw` = "/home/rsg/Projects/20220309_Christophe_GC-paper/data/WT/ChIP/tracks/CH244/CH244^mapped_S288c^CIJXLY.CPM.bw", 
    `MNase.bw` = "~/Projects/20230517_Lea_MNase-timecourse/nuc_cov_Pneumo-time-course.bw" 
) |> map(import, as = 'Rle') 
tracks |> 
    map(function(x) {
        x <- binnedAverage(bins, x[1:16], 'mean', na.rm = FALSE) |> coverage(weigh = 'mean')
        x[c("II", "IV", "XVI")]
    }) |> 
    imap(~ export(.x, con = paste0("inst/extdata/", .y)))

