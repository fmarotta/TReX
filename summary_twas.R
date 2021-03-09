#!/usr/bin/env Rscript --vanilla

# summary_twas.R
#
# Perform the summary-level TWAS.
#
# Federico Marotta (federico.marotta@edu.unito.it)
# Feb,Jun 2020

# Part of the code is derived from
# `https://github.com/gusevlab/fusion_twas/blob/master/FUSION.assoc_test.R`
# in particular the allele.qc function and the high-level algorithm. See also
# `Gusev et al. “Integrative approaches for large-scale transcriptome-wide
# association studies” 2016 Nature Genetics`

# IMPORTANT: the ref LD dataset must have all the possible SNPs and all the SNPs
# must have a unique name. In this script, we try to match the sumstats to the
# ref LD, but the ref LD dominates. If the SNPs have the same position but different names,
# the names in the ref LD prevails over the name in the sumstat file.
#
# IMPORTANT: the regreg must end with _REG[[:digit:]]+.

suppressPackageStartupMessages({
    library(BEDMatrix)
    library(data.table)
    library(GenomicRanges)
    library(docopt)
    # library(plink2R)
})

source("../R/utils.R")


"summary_twas

Perform the association. The optional options are just required to compute the
delta-TBA from scratch; however, it is not recommended to compute the delta-TBA
through this script: it is better to use the dedicated 'compute_delta_tba.R'.

Usage:
  summary_twas [options]
  summary_twas (-h | --help)
  summary_twas --version

Options:
  -z --zscores=<FILE>             Summary statistics
  -w --weights=<FILE>             Regression coefficients (trex_cv_params.tsv)
  -l --ld=<FILE>                  Reference LD
  -r --regreg=<BED_FILE>          Regulatory regions
  -t --tba=<R_OBJECT>             Pre-computed TBAs

  -f --fasta=<FASTA_FILE>         Sequence of the reference genome (optional)
  -p --pwm=<TRANSFAC>             Positional Weight Matrices (optional)
  -b --background=<BG_FILE>       Background to compute the TBA (optional)

  -s --strip-ensg-version=<BOOL>  If the gene IDs are Ensembl IDs, remove the version.
                                  In the output file, the version in the 'z scores' dataset is
                                  added back. If the IDs are not in Ensemlb format, nothing is
                                  done even if true is specified. Possible values are 'true' and
                                  'false' (or their abbreviations). [default: true]
  -o --output=<FILE>              Output file
  -j --threads=<N>                Number of cores for parallel computations [default: 1]
  -v --verbose=<N>                Verbosity level (0, 1, or 2) [default: 1]
  -h --help                       Show this message

" -> doc
argv <- docopt(gsub(" \n\\s+", " ", x = doc, perl = T))


allele.qc = function(a1,a2,ref1,ref2) {
    # from https://github.com/gusevlab/fusion_twas
    a1 = toupper(a1)
    a2 = toupper(a2)
    ref1 = toupper(ref1)
    ref2 = toupper(ref2)

    ref = ref1
    flip = ref
    flip[ref == "A"] = "T"
    flip[ref == "T"] = "A"
    flip[ref == "G"] = "C"
    flip[ref == "C"] = "G"
    flip1 = flip

    ref = ref2
    flip = ref
    flip[ref == "A"] = "T"
    flip[ref == "T"] = "A"
    flip[ref == "G"] = "C"
    flip[ref == "C"] = "G"
    flip2 = flip;

    snp = list()
    snp[["keep"]] = !((a1 == "A" & a2 == "T") | (a1 == "T" & a2 == "A") | (a1 == "C" & a2 == "G") | (a1 == "G" & a2 == "C"))
    snp[["keep"]][ a1 != "A" & a1 != "T" & a1 != "G" & a1 != "C" ] = F
    snp[["keep"]][ a2 != "A" & a2 != "T" & a2 != "G" & a2 != "C" ] = F
    snp[["flip"]] = (a1 == ref2 & a2 == ref1) | (a1 == flip2 & a2 == flip1)

    snp
}

compute.tba <- function(fasta.file, background.file, pwm.file, regreg, snps, threads) {
    suppressPackageStartupMessages({
        library(BSgenome)
        library(MatrixRider)
        library(motifStack)
        library(stringi)
        library(universalmotif)
    })

    get_eta <- function(start, remaining, done, label = "") {
        elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
        eta <- remaining * elapsed / done
        eta.h <- floor(eta / 3600)
        eta.m <- floor(eta / 60 - eta.h * 60)
        eta.s <- floor(eta - 3600 * eta.h - 60 * eta.m)
        paste0(label,
               "  eta: ",
               stri_pad_left(eta.h, 2, "0"), ":",
               stri_pad_left(eta.m, 2, "0"), ":",
               stri_pad_left(eta.s, 2, "0"))
    }

    # Ref sequence
    fa <- readDNAStringSet(fasta.file)

    # PWM & background
    if (is.null(background.file)) {
        message("Using uniform background...")
        bg <- c(A = .25, C = .25, G = .25, T = .25)
    } else {
        bg <- as.numeric(strsplit(readLines(background.file), "\t")[[1]])
        names(bg) <- c("A", "C", "G", "T")
    }

    pfm <- motifStack::importMatrix(pwm.file, format = "transfac")
    pfm <- lapply(pfm, function(t) {
        t@background <- bg # Only this seems to work
        universalmotif::convert_motifs(pcm2pfm(t), "TFBSTools-PFMatrix")
    })

    tba <- list(
        ref = matrix(nrow = 0, ncol = length(pfm)),
        alt = matrix(nrow = 0, ncol = length(pfm))
    )
    for (chr in seqlevels(fa)) {
        rrchr <- chr
        seqlevelsStyle(rrchr) <- "NCBI"
        if (!(rrchr %in% seqlevels(regreg) && rrchr %in% seqlevels(snps)))
            next()
        rr <- keepSeqlevels(regreg, rrchr, pruning.mode = "coarse")
        rrsnps <- subsetByOverlaps(snps, rr)
        rrsnps <- rrsnps[max(length(rrsnps$A1), length(rrsnps$A2)) == 1] # Remove INDELS, because we use the "replaceletterat" function, not "insert/remove letter at".

        # Get ref and alt sequences
        # use "replaceletterat" twice: one to put effect alleles in the ref, one to put the other allele in the alt.
        # We don't care about ref and alt, we only care about effect and other.
        ref <- DNAStringSet(replaceLetterAt(DNAStringSet(fa[[chr]])[[1]],
                                            start(rrsnps),
                                            as.character(rrsnps$A2),
                                            if.not.extending = "replace",
                                            verbose = FALSE))
        alt <- DNAStringSet(replaceLetterAt(DNAStringSet(fa[[chr]])[[1]],
                                            start(rrsnps),
                                            as.character(rrsnps$A1),
                                            if.not.extending = "replace",
                                            verbose = FALSE))
        names(ref) <- names(alt) <- rrchr

        # compute the tba
        N <- length(pfm)
        for (dna in c("ref", "alt")) {
            message("\nConsidering ", chr, " ", dna, "...")
            start <- Sys.time()
            chr.tba <- mclapply(mc.cores = threads, seq_len(N), function(i) {
                tf <- pfm[[i]]
                l <- length(tf)
                neighb <- promoters(rrsnps, upstream = l - 1, downstream = l)

                # TODO: truncate the neighb at the boundaries of a regreg

                # TODO: combine the tbas of neighbourhoods that overlap

                seq <- getSeq(get(dna), neighb)
                names(seq) <- neighb$SNPID
                res <- sapply(seq, MatrixRider::getSeqOccupancy, tf, 0)

                eta <- get_eta(start, N - i, i, paste(dna, names(pfm)[i], sep = ", "))
                cat("\r", stri_pad_right(eta, 60, " "))

                res
            })
            chr.tba <- do.call(cbind, chr.tba)
            colnames(chr.tba) <- names(pfm)
            tba[[dna]] <- rbind(tba[[dna]], chr.tba)
        }
    }

    tba
}

read_plink_custom_fast <- function(root, impute = c('none', 'avg', 'random')) {
    # from https://github.com/gusevlab/fusion_twas/issues/13#issue-484106760
    if(impute == 'random') {
        stop("The 'impute' random option has not been implemented.", call. = FALSE)
    }
    if(!requireNamespace("data.table", quietly = TRUE)) {
        stop("Please install the 'data.table' package using install.packages('data.table').", call. = FALSE)
    }

    ## structure from https://github.com/gabraham/plink2R/blob/master/plink2R/R/plink2R.R
    proot <- path.expand(root)

    bedfile <- paste(proot, ".bed", sep="")
    famfile <- paste(proot, ".fam", sep="")
    bimfile <- paste(proot, ".bim", sep="")

    ## Could change this code to use data.table
    bim <- data.table::fread(
        bimfile,
        colClasses = list('character' = c(2, 5, 6), 'integer' = c(1, 3, 4)),
        col.names =  c("CHROM", "ID", "POS_CM", "POS", "A1", "A2"),
        showProgress = FALSE,
        # data.table = FALSE
    )
    fam <- data.table::fread(famfile,
                             colClasses = list('character' = 1:2, 'integer' = 3:6),
                             showProgress = FALSE,
                             # data.table = FALSE
    )
    ## Set the dimensions
    geno <- BEDMatrix::BEDMatrix(bedfile, n = nrow(fam), p = nrow(bim))

    ## Convert to a matrix
    geno <- as.matrix(geno)
    if(impute == 'avg') {
        ## Check if any are missing
        geno_na <- is.na(geno)
        if(any(geno_na)) {
            means <- colMeans(geno, na.rm = TRUE)
            geno[geno_na] <- rep(means, colSums(geno_na))
        }
    }
    ## Extract data using the data.table syntax
    ## in case fread(data.table = TRUE) was used (which is the default)
    colnames(geno) <- bim[[2]]
    rownames(geno) <- paste(fam[[1]], fam[[2]], sep=":")
    # colnames(geno) <- bim[,2]
    # rownames(geno) <- paste(fam[,1], fam[, 2], sep=":")

    ## If you used fread(data.table = TRUE) then you need to cast
    ## the objects using as.data.frame()
    list(bed=geno, fam=fam, bim=bim)
}


# Ref LD
genos <- read_plink_custom_fast(argv$ld, impute = "avg")

# Return all the SNPs of the ref LD dataset in GRanges format
snps <- GRanges(seqnames = genos$bim$CHROM,
                ranges = IRanges(genos$bim$POS, genos$bim$POS),
                SNPID = genos$bim$ID,
                A1 = genos$bim$A1,
                A2 = genos$bim$A2)
seqlevelsStyle(snps) <- "NCBI"

# Regulatory regions
regreg <- as(fread(argv$regreg,
                   select = 1:4,
                   col.names = c("CHR", "START", "END", "REGION")),
             "GRanges")
seqlevelsStyle(regreg) <- "NCBI"


# Compute the TBA of all the SNPS that fall in a regreg, even if we don't
# have their Zscore
if (!file.exists(argv$tba)) {
    warning("This method for computing the delta-TBA is deprecated. ",
            "You should use the dedicated script 'compute_delta_tba.sh' instead")
    message("Computing the 'delta TBA'...")
    tba <- compute.tba(argv$fasta, argv$background, argv$pwm, regreg, snps, argv$threads)
    dir.create(dirname(argv$tba), recursive = TRUE)
    saveRDS(tba, argv$tba)
} else {
    message("Using precomputed 'delta TBA'...")
    tba <- readRDS(argv$tba)
}

# Compute the delta TBA
if (argv$verbose > 0)
    message("The 'delta TBA' matrix has ", nrow(tba$ref),
            " SNPs and ", ncol(tba$ref), " TFs")
delta <- log(tba$alt) - log(tba$ref) # Assuming additive effects

# Sumstats
sumstat = fread(argv$zscores,
                header = TRUE,
                col.names = c("CHROM", "ID", "POS", "A1", "A2", "TYPE", "Z", "R2", "N", "P"))

# Find common SNPs and establish a common ordering
sumstat <- sumstat[!duplicated(sumstat$ID), ]
sumstat <- na.omit(sumstat)
sumstat$CHROM <- as.integer(sumstat$CHROM)
sumstat$POS <- as.integer(sumstat$POS)
genos$bim <- genos$bim[!duplicated(genos$bim$ID), ]

sumstat$A1 <- toupper(sumstat$A1)
sumstat$A2 <- toupper(sumstat$A2)
sumstat$REF <- toupper(sumstat$REF)
genos$bim$A1 <- toupper(genos$bim$A1)
genos$bim$A2 <- toupper(genos$bim$A2)

m <- merge(sumstat, genos$bim, by = c("CHROM", "POS"))
m <- m[(A1.x == A1.y & A2.x == A2.y) | (A1.x == A2.y & A2.x == A1.y), ]

sumstat <- sumstat[m$ID.x, c("ID", "A1", "A2", "Z"), on = "ID"]
genos$bed <- genos$bed[, m$ID.y]
genos$bim <- genos$bim[m$ID.y, , on = "ID"]

# QC / allele-flip the input and output
qc = allele.qc( sumstat$A1 , sumstat$A2 , genos$bim$A1 , genos$bim$A2 )

# Flip Z-scores for mismatching alleles
sumstat$Z[ qc$flip ] = -1 * sumstat$Z[ qc$flip ]
sumstat$A1[ qc$flip ] = genos$bim[qc$flip]$A1
sumstat$A2[ qc$flip ] = genos$bim[qc$flip]$A2

# Remove strand ambiguous SNPs (if any)
if ( sum(!qc$keep) > 0 ) {
    genos$bim = genos$bim[qc$keep,]
    genos$bed = genos$bed[,qc$keep]
    sumstat = sumstat[qc$keep,]
}

# Return the SNPs in GRanges format
snps <- GRanges(seqnames = genos$bim$CHROM,
                ranges = IRanges(genos$bim$POS, genos$bim$POS),
                SNPID = genos$bim$ID,
                A1 = sumstat$A1, # genos$bim$V5,
                A2 = sumstat$A2, # genos$bim$V6,
                Z = sumstat$Z)
seqlevelsStyle(snps) <- "NCBI"

# Weights
all_betas <- fread(argv$weights, select = c("data.name", colnames(delta)))
if (pmatch(tolower(argv$`strip-ensg-version`), "true", nomatch = FALSE)) {
    all_betas$data.name <- strip_ensg_version(all_betas$data.name)
}

# # Weights from fit
# fit <- readRDS(argv$weights)
# fit[sapply(fit, is.null)] <- NULL
# names(fit) <- strip_ensg_version(names(fit))

twas <- mclapply(mc.cores = argv$threads, all_betas$data.name, function(gene) {
# twas <- lapply(all_betas$data.name, function(gene) {
    if (argv$verbose > 0)
        message("Considering gene ", gene, "...")

    # We use all the SNPs that fall inside
    # one of the regulatory regions for that gene

    tryCatch({
        # Find the chromosome and the gene
        w <- regreg[grep(gene, regreg$REGION)]
        chr <- unique(as.character(seqnames(w)))
        if (is.null(chr) || length(chr) != 1) {
            warning("Couldn't find a regulatory region for ", gene)
            return(NULL)
        }
        twas_gene <- unique(sub("_REG[[:digit:]]+$", "", w$REGION))
        if (length(twas_gene) != 1) {
            warning("Couldn't find a matching gene for ", gene)
            return(NULL)
        }

        # Get the betas
        betas <- unlist(all_betas[data.name == gene, !"data.name"])

        # find the snps that affect that gene, i.e. those that fall in the regreg of that gene
        ov <- findOverlaps(regreg, snps)
        from <- which(grepl(gene, regreg$REGION))
        id <- snps[to(ov[from(ov) %in% from])]$SNPID

        if (argv$verbose > 1)
            message("This gene has ", length(id), " SNPs")
        if (!length(id)) {
            return(data.table(chr, twas_gene, 0, NA, NA))
        }

        # We use sapply just to ensure that weights and zetas have the same order
        # and are named vectors
        zetas <- sapply(id, function(snp) {
            snps[snps$SNPID == snp]$Z
        })

        # for each snp, the weight of that snp is the scalar product of beta and delta tba
        weights <- sapply(id, function(snp) {
            betas %*% delta[rownames(delta) == snp, names(betas)]
        })

        cur.genos = scale(genos$bed[, id])
        LD = t(cur.genos) %*% cur.genos / (nrow(cur.genos) - 1)

        n <- as.numeric(weights %*% zetas)
        d <- as.numeric(weights %*% LD %*% weights)

        twas.z <- n / sqrt(d)
        twas.p <- 2 * (pnorm(abs(twas.z), lower.tail = F))
        return(data.table(chr, twas_gene, length(id), twas.z, twas.p))
    },
    error = function(e) {
        message(e$message)
        return(NULL)
    })
})

twas <- rbindlist(twas)
names(twas) <- c("CHR", "GENE", "N_SNPS", "ZSCORE", "PVAL")
fwrite(twas, argv$output, sep = "\t", quote = FALSE)

