#!/usr/bin/env Rscript

# fit_model.R
#
# For each gene, perform a 5-fold cross-validation and save the R^2 in
# order to evaluate the performance of the model. Do a nested
# cross-validation to tune the parameters. Also fit the best model.
#
# Federico Marotta (federico.marotta@edu.unito.it),
# Feb,Jun,Sep 2020

suppressPackageStartupMessages({
    library(docopt)
    library(fplyr)
})

"fit_model

Train the model EXPR ~ TBA, using a nested cross-validation to evaluate
the performance and tune the parameters. After the cross-validation,
if the performance R^2 is greater than what specified by the `min-R2`
argument, the model with the optimised parameters is trained on the
whole data. Pass 'NA' to the --seed option in order to avoid setting the
seed.

Usage:
  fit_model [options] <EXPR_FILE> <TBA_FILE>
  fit_model (-h | --help)

Arguments:
  EXPR_FILE                  Gene expression.
  TBA_FILE                   TBA matrix after aggregation.

Options:
  -f --folds=<N>             Pass NULL or <= 1 to avoid doing the nested CV [default: 5]
  -n --nested-folds=<N>      Number of folds in the inner loop [default: 10]
  -m --min-R2=<U>            Minimum R^2 for which the model is fit [default: 0]
  -t --train-samples=<FILE>  List of training samples [default: 'all of them']
  -o --outdir=<DIR>          Path to the output directory [default: ./]
  -j --threads=<N>           Number of cores for parallel computations [default: 1]
  -s --seed=<N>              Seed to train the model [default: 2020]
  -v --verbose=<N>           Verbosity level from 0 to 4 [default: 1]
  -h --help                  Show this message

" -> doc
argv <- docopt(gsub(" \n\\s+", " ", x = doc, perl = T))

# Acquire the regression functions and parameters
source("utils/trex_models.R")


# Read the expression
expr <- fread(argv$EXPR_FILE)
names(expr)[1] <- "GENE"

# Read the samples
if (argv$`--train-samples` != "all of them") {
	samples <- readLines(argv$`--train-samples`)
	samples <- intersect(samples, names(expr))
} else {
    samples <- names(expr)[2:length(expr)]
}

# Reshape the expression file
expr <- melt(expr,
             id.vars = "GENE",
             measure.vars = samples,
             value.name = "EXPRESSION",
             variable.name = "IID",
             variable.factor = FALSE)[order(GENE, IID), ]

# Define output paths
# Create the 'perf' directory (with three files) if the outer CV is to be made;
# otherwise create only the 'cv' and 'fit' directories.
if (!is.null(argv$folds) & as.integer(argv$folds) > 1) {
    if (!dir.exists(paste(argv$outdir, "perf", sep = "/")))
        if (!dir.create(paste(argv$outdir, "perf", sep = "/"),
                        recursive = T, showWarnings = F))
            stop("Could not create output directory.")
    out <- c(
        perf_inner = paste(argv$outdir, "perf/trex_inner_cv.tsv", sep = "/"),
        perf_outer = paste(argv$outdir, "perf/trex_outer_cv_full.tsv", sep = "/"),
        perf_summ = paste(argv$outdir, "perf/trex_outer_cv_summary.tsv", sep = "/")
    )
} else {
    out <- c()
}

if (argv$`--min-R2` < 1) {
    if (!dir.exists(paste(argv$outdir, "cv", sep = "/")))
        if (!dir.create(paste(argv$outdir, "cv", sep = "/"),
                        recursive = T, showWarnings = F))
            stop("Could not create cv directory")
    if (!dir.exists(paste(argv$outdir, "fit", sep = "/")))
        if (!dir.create(paste(argv$outdir, "fit", sep = "/"),
                            recursive = T, showWarnings = F))
            stop ("Could not create fit directory")
    out <- c(
        out,
	    cv_params = paste(argv$outdir, "fit/trex_weights.tsv", sep = "/"),
	    cv_full = paste(argv$outdir, "cv/trex_cv_full.tsv", sep = "/"),
        cv_summ = paste(argv$outdir, "cv/trex_cv_summary.tsv", sep = "/")
    )
    out.cv_fit = paste(argv$outdir, "fit/trex_fit.Rds", sep = "/")
}

# Save the params
saveRDS(list(params.grid = params.grid, argv = argv),
        paste(argv$outdir, "trex_call.Rds", sep = "/"))

# Nested Cross-Validation
cv_function <- function(tba) {
    d.train <- merge(expr, tba)
    if (nrow(d.train) == 0)
        return(vector("list", length(out) + 1))

    if (argv$verbose > 0)
        message("Considering gene ", d.train$GENE[1], "...")

    if (argv$seed != "NA")
        set.seed(argv$seed)

    if (!is.null(argv$folds) && as.integer(argv$folds) > 1) {
        if (argv$verbose > 0)
            message("Nested Cross-Validation...")
        perf <- NestedCV(d.train,
                         params.grid = params.grid,
                         n.outer.folds = as.integer(argv$folds),
                         n.inner.folds = as.integer(argv$`--nested-folds`),
                         verbose = as.integer(argv$verbose) - 1)
        l <- list(perf$inner$full, perf$outer$full, perf$outer$summ)
    } else {
        l <- list()
    }

    # If the outer CV was not performed, or if the R2 of the outer CV is greater
    # than the user-specified threshold, fit the best model
    if ((is.null(argv$folds) || as.integer(argv$folds) <= 1) ||
        (!is.null(perf$outer$summ) && perf$outer$summ$pred_perf_R2 >= as.numeric(argv$`--min-R2`))) {
        if (argv$verbose > 0)
            message("Fitting best model...")
        cv <- CV(d.train,
                 NULL,
                 params.grid = params.grid,
                 n.folds = as.integer(argv$`--nested-folds`),
                 fit.best = TRUE,
                 verbose = as.integer(argv$verbose) - 1)
        if (!is.null(cv$fit)) {
		    betas <- cbind(data.name = d.train$GENE[1], Coef(cv$fit))
		    l <- append(l, list(betas, cv$full, cv$summ, cv$fit))
        } else {
            return(append(l, list(NULL, NULL, NULL, NULL)))
        }
    } else {
        return(append(l, list(NULL, NULL, NULL, NULL)))
    }

	l
}

fit <- fmply(argv$TBA_FILE,
             out,
             FUN = cv_function,
             parallel = as.integer(argv$threads),
             header = TRUE)

# Save the fit object
fit[sapply(fit, is.null)] <- NULL
if (length(fit))
    saveRDS(fit, out.cv_fit)
