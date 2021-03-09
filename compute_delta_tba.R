#!/usr/bin/env Rscript --vanilla

# compute_delta_tba.R
#
# For each SNP, compute the mean TBA for the individuals with genotype
# 0, 1, or 2. Store everything into a list of three matrices and save the
# corresponding .Rds. The difference between the 1 and the 0 matrices is the
# delta TBA.
#
# Federico Marotta (federico.marotta@edu.unito.it)
# Nov 2020

suppressPackageStartupMessages({
    library(GenomicRanges)
    library(docopt)
    library(fplyr)
})

source("utils/utils.R")


"compute_delta_tba

For each SNP, compute the mean TBA for the individuals with genotype 0,
1, or 2. Store everything into a list of three matrices and save the
corresponding .Rds. The difference between the 1 and the 0 matrices is
the delta TBA.

Usage:
  compute_delta_tba [-j <N>] [-l <FRACTION>] [-k <SAMPLES>] -s <SNPS_FILE> -v <VCF> -t <TBA> -b <BED> -o <DELTA_RDS>
  compute_delta_tba (-h | --help)

Options:
  -v --vcf=<VCF>              Reference LD dataset in VCF format
  -t --tba=<TBA>              TBA as computed by the compute_tba script
  -b --bed=<BED>              BED file with the regulatory regions
  -s --snps=<SNPS_FILE>       SNPs output of vcf_rider
  -k --keep=<SAMPLES_LIST>    Text file with a list of individuals to keep
  -l --too-large=<FRACTION>   Fraction of indels considered too large [default: 0.1]
  -o --output=<DELTA_RDS>     Path to the output file
  -j --threads=<N>            Number of cores for parallel computations [default: 1]
  -h --help                   Show this message

" -> doc
argv <- docopt(gsub(" \n\\s+", " ", x = doc, perl = T))

# Find the regions to be excluded due to length problems
invalid_regions <- exclude_from_indel_length(argv$snps, as.numeric(argv$too_large))
regreg <- fread(argv$bed,
                col.names = c("CHR", "START", "END", "ID"))
regreg <- regreg[!ID %in% invalid_regions]
regreg$END <- regreg$END - 1
regreg <- as(regreg, "GRanges")
seqlevelsStyle(regreg) <- "UCSC"

vcf <- fread(paste("zcat", argv$vcf), drop = 6:9)
bim <- GRanges(seqnames = vcf$`#CHROM`,
               ranges = IRanges(vcf$POS, vcf$POS),
               SNPID = vcf$ID,
               REF = vcf$REF,
               ALT = vcf$ALT)
seqlevelsStyle(bim) <- "UCSC"

ov <- findOverlaps(bim, regreg)
bim <- bim[unique(from(ov))]
vcf <- vcf[bim$SNPID, !c("#CHROM", "POS", "REF", "ALT"), on = "ID"]
vcf <- melt(vcf,
            id.vars = "ID",
            variable.name = "IID",
            value.name = "GENOTYPE",
            variable.factor = FALSE)
if (!is.null(argv$keep)) {
    samples <- readLines(argv$keep)
    vcf <- vcf[IID %in% samples]
}

vcf$GENOTYPE <- tstrsplit(vcf$GENOTYPE, ":", fixed = T)[[1]]
alleles <- tstrsplit(vcf$GENOTYPE, "|", fixed = T)
vcf$GENOTYPE <- as.numeric(alleles[[1]]) + as.numeric(alleles[[2]])

l <- flply(argv$tba, function(d) {
    gene <- d$GENE[1]
    rr <- regreg[grepl(gene, regreg$ID)]
    snps <- bim[unique(from(findOverlaps(bim, rr)))]$SNPID
    lapply(0:2, function(g) {
        m <- do.call(rbind, lapply(snps, function(id) {
            iids <- vcf[ID == id & GENOTYPE == g]$IID
            colMeans(d[IID %in% iids, !c("GENE", "IID")])
        }))
        if (length(m)) {
            rownames(m) <- snps
            return(m)
        }
        return(NULL)
    })
}, parallel = as.integer(argv$threads))

res <- list(bim = bim, tba_by_genotype = l)
saveRDS(res, argv$output)
