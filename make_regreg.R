#!/usr/bin/env Rscript

# make_regreg.R
#
# Define the default region for each gene. Merge the regions associated
# to each gene.
#
# Federico Marotta (federico.marotta@edu.unito.it)
# Jan,Apr,Jun 2020


suppressPackageStartupMessages({
    library(GenomicFeatures)
    library(data.table)
    library(docopt)
    library(parallel)
})

"make_regreg

Given the transcripts annotation, return a BED file of the regulatory regions
associated to each gene. By default, each transcript is assigned a region
(promoter) around its TSS. If a BED file with enhancers is provided, to each
gene it is also assigned its closest enhancer.

Usage:
  make_regreg [-u <bp>] [-d <bp>] [-e <BED>] [-c <chrN>] [-o <BED>] [-j <n>] <TRANSCRIPTS>
  make_regreg (-h | --help)
  make_regreg --version

Arguments:
  TRANSCRIPTS                 A GTF or GFF3 file with the transcripts annotation

Options:
  -u --upstream=<N>           Promoter start upstream of the TSS [default: 1500]
  -d --downstream=<N>         Promoter end downstream of the TSS [default: 500]
  -e --enhancers=<BED_FILE>   BED with the enhancers locations
  -c --chromosome=<chrN>      Consider only the specified chromosome [default: all]
  -o --output=<BED_FILE>      Path to the output file [default: STDOUT]
  -j --threads=<N>            Number of cores for parallel computations [default: 1]
  -h --help                   Show this message
  --version                   Show the version number

" -> doc
argv <- docopt(gsub(" \n\\s+", " ", x = doc, perl = T),
               version = 'make_regreg 0.0.1')

if (argv$output == "STDOUT") {
    argv$output <- ""
} else {
    if (!dir.exists(dirname(argv$output)))
        if (!dir.create(dirname(argv$output), recursive = T))
            stop("Failed to create the output directory. Please make sure ",
                 "that you have\nall the relevant permissions")
    if (!file.create(argv$output))
        stop("Failed to create the output file. Please make sure that the ",
             "directory\nexists and you have all the relevant permissions.")
}

# Transcripts
txdb <- makeTxDbFromGFF(argv$TRANSCRIPTS)
seqlevelsStyle(txdb) <- "UCSC"

# Enhancers
if (!is.null(argv$enhancers)) {
    enhancers <- as(fread(argv$enhancers,
                          select = 1:4,
                          col.names = c("CHR", "START", "END", "REGION")),
                    "GRanges")
    seqlevelsStyle(enhancers) <- "UCSC"
}

# Define the chromosomes
if (argv$chromosome != "all") {
    if ("UCSC" %in% seqlevelsStyle(argv$chromosome))
        chroms <- argv$chromosome
    else if ("NCBI" %in% seqlevelsStyle(argv$chromosome))
        chroms <- paste0("chr", argv$chromosome)
} else {
    chroms <- as.character(seqlevels(txdb))
}
if (!is.null(argv$enhancers)) {
    chroms <- intersect(chroms, as.character(seqlevels(enhancers)))
}

chroms <- chroms[order(nchar(chroms), chroms)]
for (chr in chroms) {
    message("Considering chromosome ", chr, "...")

    seqlevels(txdb) <- chr
    txgr <- transcriptsBy(txdb, by = "gene")

    regreg <- mclapply(seq_along(txgr), function(i) {
        # Assign the default regions to each transcript
        proms <- promoters(txgr[[i]],
                           upstream = as.integer(argv$upstream),
                           downstream = as.integer(argv$downstream))

        # This union of the regions with themselves is a trick to 'merge' them
        regreg <- union(proms, proms)

        # Assign the enhancers
        # To each gene, we assign the closest enhancers, as in (E) of Fig. 3 in
        # https://www.nature.com/articles/s41588-019-0538-0
        if (!is.null(argv$enhancers)) {
            enhs <- enhancers[nearest(genes(txdb)[i], enhancers)]
            regreg <- union(regreg, enhs)
        }

        data.table(
            CHR = as.character(seqnames(regreg)),
            START = format(start(regreg) - 1, scientific = FALSE, trim = TRUE),
            END = format(end(regreg), scientific = FALSE, trim = TRUE),
            NAME = paste0(names(txgr)[[i]], "_REG", 1:length(regreg)),
            # SCORE = ".",
            # STRAND = strand(gr), # doesn't make sense for regreg
            stringsAsFactors = FALSE
        )
    }, mc.cores = as.integer(argv$threads))

    # Exclude regions that gave an error (the "try-error" class is given automatically
    # by mclapply whenever a job fails; healthy regreg objects should be of class
    # "data.table" and "data.frame")
    regreg <- lapply(regreg, function(r) {
        if ("try-error" %in% class(r)) {
            r
        } else {
            message(r)
            NULL
        }
    })

    # Print the output
    df <- rbindlist(regreg)
    fwrite(df[order(nchar(CHR), CHR, as.integer(START))],
           argv$output,
           append = TRUE,
           col.names = FALSE,
           quote = FALSE,
           sep = "\t")

    # Reset the seqlevels of the txdb object
    seqlevels(txdb) <- seqlevels0(txdb)
    seqlevelsStyle(txdb) <- "UCSC"
}
