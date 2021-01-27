#!/usr/bin/env Rscript

# aggregate_delta_tba.R
#
# Compute the mean of all the REF and the mean of all the ALT alleles in the
# dataset. Store everything into a list of two matrices and save the
# corresponding .Rds. The difference between these two matrices is the
# delta TBA.
#
# Federico Marotta (federico.marotta@edu.unito.it)
# Aug-Sep 2020

suppressPackageStartupMessages({
    library(docopt)
    library(fplyr)
})

source("utils/utils.R")

"aggregate_delta_tba

Aggregate the delta TBA into a matrix.

Usage:
  aggregate_delta_tba [-l <FRACTION>] [-j <N>] -v <VCF> -t <DELTA_TBA> -s <SNPS_FILE> -o <DELTA_RDS>
  aggregate_delta_tba (-h | --help)

Options:
  -v --vcf=<VCF>              Reference LD dataset as a phased VCF.
  -t --tba=<DELTA_TBA>        TBA output of vcf_rider.
  -s --snps=<SNPS_FILE>       SNPs output of vcf_rider.
  -l --too-large=<FRACTION>   Fraction of indels considered too large [default: 0.1]
  -o --output=<DELTA_RDS>     Path to the output file
  -j --threads=<N>            Number of cores for parallel computations [default: 1]
  -h --help                   Show this message

" -> doc
argv <- docopt(gsub(" \n\\s+", " ", x = doc, perl = T))


# Read the genotypes (with impute = avg)
vcf <- fread(argv$vcf, drop = c(1:2, 4:9))
vcf <- melt(vcf, id.vars = "ID", variable.name = "IID", value.name = "GENOTYPE")
vcf$GENOTYPE <- tstrsplit(vcf$GENOTYPE, ":", fixed = T)[[1]]
vcf[, c("allele1", "allele2") := tstrsplit(vcf$GENOTYPE, "|", fixed = T)]
vcf <- melt(vcf[, !"GENOTYPE"], id.vars = c("ID", "IID"), variable.name = "ALLELE", value.name = "HAPLOTYPE")

# Find the regions to be excluded due to length problems
invalid_snps <- exclude_from_indel_length(argv$snps, as.numeric(argv$`too-large`))

# Read the tba
res <- flply(argv$tba,
             drop = c(2, 3),
             col.names = c("WINDOW", "MODEL", "IID", "ALLELE", "TBA"),
             parallel = as.integer(argv$threads),
             function(d) {
    d <- d[!WINDOW %in% invalid_snps, ]
    if (is.null(d) || nrow(d) == 0)
        return(NULL)

    d <- merge(d, vcf, by.x = c("WINDOW", "IID", "ALLELE"), by.y = c("ID", "IID", "ALLELE"))
    if (is.null(d) || nrow(d) == 0)
        return(NULL)

    ref <- d[HAPLOTYPE == 0, .(TBA = mean(TBA)), by = c("WINDOW", "MODEL")]
    ref <- dcast(ref, WINDOW ~ MODEL, value.var = "TBA")

    alt <- d[HAPLOTYPE == 1, .(TBA = mean(TBA)), by = c("WINDOW", "MODEL")]
    alt <- dcast(alt, WINDOW ~ MODEL, value.var = "TBA")

    list(ref, alt)
})

ref <- rbindlist(lapply(res, `[[`, 1))
nref <- ref$WINDOW
ref <- as.matrix(ref[, !"WINDOW"])
rownames(ref) <- nref

alt <- rbindlist(lapply(res, `[[`, 2))
nalt <- alt$WINDOW
alt <- as.matrix(alt[, !"WINDOW"])
rownames(alt) <- nalt

saveRDS(list(ref = ref, alt = alt), argv$output)
