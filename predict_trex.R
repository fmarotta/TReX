#!/usr/bin/env Rscript

# predict_trex.R
#
# For each gene, predict the TBA-regulated expression.
#
# Federico Marotta (federico.marotta@edu.unito.it),
# Feb,Jun,Sep 2020

suppressPackageStartupMessages({
    library(docopt)
    library(fplyr)
})

source("../R/trex_models.R")


"predict_trex

Impute the expression. If the true expression values are provided, some
accuracy statistics are computed as well.

Usage:
fit_model [options] <FIT_OBJ> <TBA_FILE>
fit_model (-h | --help)
fit_model --version

Arguments:
  FIT_OBJ                         Rds with the model.
  TBA_FILE                        TBA matrix after aggregation.

Options:
  -t --test-samples=<FILE>        Samples for which the TReX will be computed. [default: all of them]
  -e --expression=<FILE>          True expression values (optional)
  -s --strip-ensg-version=<BOOL>  If the gene IDs are Ensembl IDs, remove the version.
                                  In the output file, the version in the 'prediction' dataset is
                                  added back. If the IDs are not in Ensemlb format, nothing is
                                  done even if true is specified. Possible values are 'true' and
                                  'false' (or their abbreviations). [default: true]
  -o --outdir=<DIR>               Path to the output directory [default: ./]
  -j --threads=<N>                Number of cores for parallel computations [default: 1]
  -v --verbose=<N>                Level of verbosity (0 or 1) [default: 1]
  -h --help                       Show this message

" -> doc
argv <- docopt(gsub(" \n\\s+", " ", x = doc, perl = T))


fit <- readRDS(argv$FIT_OBJ)
if (pmatch(tolower(argv$`strip-ensg-version`), "true", nomatch = FALSE)) {
    names(fit) <- strip_ensg_version(names(fit))
}


# We allow expression to be null. If so, don't compute statistics such as
# the cor between true and pred.
if (!is.null(argv$expression)) {
	# Read the expression file
    true_expr <- fread(argv$expression)
    names(true_expr)[1] <- "GENE"

    # Reshape it
    true_expr <- melt(true_expr,
                      id.vars = "GENE",
                      # measure.vars = samples, # samples may be NULL
                      value.name = "EXPRESSION",
                      variable.name = "IID",
                      variable.factor = FALSE)[order(GENE, IID), ]
}

# Output files
dir.create(argv$outdir)
out <- c(
	pred_expr = paste(argv$outdir, "trex_pred_expr.tsv", sep = "/"),
	pred_stat = if (exists("true_expr"))
	    paste(argv$outdir, "trex_pred_summary.tsv", sep = "/")
)
null <- lapply(out, file.remove, showWarnings = F)

# Prediction function
pred_function <- function(tba) {
    # Initialise the return value
    if (exists("true_expr"))
        r <- list(NULL, NULL)
    else
        r <- list(NULL)

    if (argv$`test-samples` != "all of them") {
        samples <- readLines(argv$`test-samples`)
        tba <- tba[IID %in% samples]
    }
    if (nrow(tba) == 0 || ncol(tba) == 0)
        return(r)

    if (pmatch(tolower(argv$`strip-ensg-version`), "true", nomatch = FALSE)) {
        gene <- strip_ensg_version(tba$GENE[1])
    } else {
        gene <- tba$GENE[1]
    }
    fit <- fit[[gene]]
    if (is.null(fit))
        return(r)

    if (argv$verbose > 0)
        message("Considering gene ", gene, "...")

    ypred <- predict(fit, as.matrix(tba[, !c("GENE", "IID")]))

    pred_expr <- data.table(GENE = tba$GENE, IID = tba$IID, TREX = ypred)
    names(pred_expr) <- c("GENE", "IID", "TREX")

    r[[1]] <- pred_expr

    if (exists("true_expr")) {
        d.test <- merge(pred_expr, true_expr)
        tryCatch({
            pearson <- compute_pearson(d.test$EXPRESSION, d.test$TREX)
            pearson <- cbind(
                data.name = gene,
                mse = cvmse(d.test$EXPRESSION, d.test$TREX),
                rsq = cvrsq(d.test$EXPRESSION, d.test$TREX),
                pearson[, c("parameter.df", "estimate.cor", "p.value")]
            )
            r[[2]] <- pearson
        }, error = function(e) {
            message(e$message)
            return(r)
        })
        return(r)
    }

    list(r)
}

fmply(argv$TBA_FILE, out, pred_function, parallel = as.integer(argv$threads))
