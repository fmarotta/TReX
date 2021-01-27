#!/usr/bin/env Rscript

# individual_twas.R
#
# Perform the individual-level TWAS.
#
# Federico Marotta (federico.marotta@edu.unito.it)
# Feb,Jun,Aug-Sep 2020

suppressPackageStartupMessages({
    library(docopt)
    library(fplyr)
})

"individual_twas

Associate the predicted expression to the phenotype. Use characters to force
factors. The first column must contain the individual ID, the second the
phenotype, and from the third onwards there must be the covariates.

Usage:
  individual_twas [options] <PHENO_FILE> <TREX_FILE>
  individual_twas (-h | --help)

Arguments:
  TREX_FILE                         The predicted TBA-regulated expression
  PHENO_FILE                        A file with phenotypes

Options:
  -t --twas-samples=<FILE>          [default: all of them]
  -c --covar=<FILE>                 Path to the covariates (optional)
  -q --pheno-quantile-normalize     Apply an inverse normal transformation to the
                                    quantitative phenotypes
  -o --outdir=<DIR>                 Path to the output file [default: ./]
  -j --threads=<N>                  Number of cores for parallel computations [default: 1]
  -v --verbose=<N>                  Level of verbosity (0, 1, or 2) [default: 1]
  -h --help                         Show this message

" -> doc
argv <- docopt(gsub(" \n\\s+", " ", x = doc, perl = T))


pheno <- fread(argv$PHENO_FILE,
               stringsAsFactors = TRUE,
               na.strings = c(getOption("datatable.na.strings", "NA"), ""))
w <- sapply(pheno, function(x) length(unique(na.omit(x))) > 1)
if (length(which(!w))) {
    warning("Removing ", paste(names(which(!w)), collapse = ", "),
            " from the phenotypes (all values are constant).")
    pheno <- pheno[, ..w]
}
names(pheno)[1] <- "IID"
pheno$IID <- as.character(pheno$IID)

if (argv$`twas-samples` != "all of them") {
    samples <- readLines(argv$`twas-samples`)
    pheno <- pheno[IID %in% samples]
}

d <- pheno

if (!is.null(argv$covar)) {
    covar <- fread(argv$covar,
                   stringsAsFactors = TRUE,
                   na.strings = c(getOption("datatable.na.strings", "NA"), ""))
    w <- sapply(covar, function(x) length(unique(na.omit(x))) > 1)
    if (length(which(!w))) {
        warning("removing ", paste(names(which(!w)), collapse = ", "),
                " from the covariates (all values are constant).")
        covar <- covar[, ..w]
    }
    names(covar)[1] <- "IID"
    covar$IID <- as.character(covar$IID)
    d <- merge(pheno, covar, by = "IID")
}

if (!dir.exists(argv$outdir))
    if (!dir.create(argv$outdir))
        stop("Could not create output directory.")
out <- paste0(argv$outdir, "/", names(pheno)[-1], ".trex_itwas.tsv")

fmply(argv$TREX_FILE, out, parallel = as.integer(argv$threads), function(expr) {
    d <- merge(expr[, c("IID", "TREX")], d, by = "IID")
    if (argv$verbose > 0)
        message("Considering gene ", expr$GENE[1], "...")
    r <- lapply(seq_len(length(pheno) - 1), function(ph) {
        ph <- ph + 1 # skip the first field
        if (!is.null(argv$covar)) {
            f <- as.formula(paste(names(pheno)[ph],
                                  paste("TREX",
                                        paste(names(covar)[-1], collapse = " + "),
                                        sep = " + "),
                                  sep = " ~ "))
        } else {
            f <- as.formula(paste(names(pheno)[ph], "TREX", sep = " ~ "))
        }
        if (argv$verbose > 1)
            message("Considering phenotype ", names(pheno)[ph], "...")
        if (is.factor(pheno[[ph]])) {
            tryCatch({
                l <- glm(f, family = binomial, d)
                allcoef <- coef(summary(l))
                excoef <- allcoef[rownames(allcoef) == "TREX"]
                if (length(excoef) == 0)
                    return(NULL)
                names(excoef) <- colnames(allcoef)
                t <- data.table(Gene = expr$GENE[1], t(excoef))
                t$Deviance <- l$deviance
                t$`DF Residual` <- l$df.residual
                t
            }, error = function(e) {
                message(e$message)
                return(NULL)
            }, warning = function(w) {
                message(w$message)
                return(NULL)
            })
        } else if (is.numeric(pheno[[ph]])) {
            tryCatch({
                if (argv$`pheno-quantile-normalize`) {
                    d[[names(pheno)[ph]]] <- qnorm(
                        (rank(d[[names(pheno)[ph]]], na.last = "keep") - 0.5) / sum(!is.na(d[[names(pheno)[ph]]]))
                    )
                }
                l <- lm(f, d)
                allcoef <- coef(summary(l))
                excoef <- allcoef[rownames(allcoef) == "TREX"]
                if (length(excoef) == 0)
                    return(NULL)
                names(excoef) <- colnames(allcoef)
                t <- data.table(Gene = expr$GENE[1], t(excoef))
                t$`R squared` <- summary(l)$r.squared
                t$`DF Residual` <- l$df.residual
                t
            }, error = function(e) {
                message(e$message)
                return(NULL)
            }, warning = function(w) {
                message(w$message)
                return(NULL)
            })
        } else {
            stop("Problem with PHENO_FILE: the type of the phenotype could not be determined.")
        }
    })

    r
})
