#!/usr/bin/env Rscript

# aggregate_tba.R
#
# Convert vcf_rider's output for the TBA in a matrix form. Moreover,
# exclude the regions whose indels are too large.
#
# Federico Marotta (federico.marotta@edu.unito.it)
# Oct 2019

suppressPackageStartupMessages({
    library(docopt)
    library(fplyr)
})

source("utils/utils.R")

"aggregate_tba

Aggregate the TBA into a matrix.

Usage:
  aggregate_tba [-l <FRACTION>] [-o <TBA_MATRIX>] [-j <N>] <SORTED_TBA> <SNPS_FILE>
  aggregate_tba (-h | --help)
  aggregate_tba --version

Arguments:
  SORTED_TBA                  Output of vcf_rider joined with the regreg and sorted.
  SNPS_FILE                   vcf_rider's SNPs output.

Options:
  -l --too-large=<FRACTION>   Fraction of indels considered too large [default: 0.1]
  -o --output=<TBA_MATRIX>    Path to the output file [default: STDOUT]
  -j --threads=<N>            Number of cores for parallel computations [default: 1]
  -h --help                   Show this message

" -> doc
argv <- docopt(gsub(" \n\\s+", " ", x = doc, perl = T))

if (argv$output == "STDOUT") {
    argv$output <- ""
} else {
    file.create(argv$output)
}

# Find the regions to be excluded due to length problems
invalid.regions <- exclude_from_indel_length(argv$SNPS_FILE, as.numeric(argv$`too-large`))

# Aggregate the TBA
ffply(argv$SORTED_TBA,
      argv$output,
      parallel = as.integer(argv$threads),
      header = FALSE,
      col.names = c("GENE", "REGION", "START", "END", "MODEL", "IID", "ALLELE", "TBA"),
      FUN = function(d, by) {
    # Exclude the problematic regions
    d <- d[!REGION %in% invalid.regions]
    if (nrow(d) == 0)
        return(NULL)

    # Sum the TBAs of the two alleles and of all regions
    d <- d[, .(TBA = log2(sum(as.numeric(TBA)))), by = c("MODEL", "IID")]
    dcast(d, IID ~ MODEL, value.var = "TBA")
})
