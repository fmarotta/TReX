# Training models for T-REx. In principle, the user can customise the
# functions in this file, for instance to run PC regression instead of
# Ridge regression.
#
# Federico Marotta (federico.marotta@edu.unito.it)
# Feb 2020

# TODO: methods for getting the slots of autoparams objects.

suppressPackageStartupMessages({
    library(glmnet)
})

source("utils/utils.R")

# Define the parameters space
params.grid <- list(alpha = 0, lambda = 2^seq(4, -18, length = 100))

NestedCV <- function(d.train, params.grid,
                     outer.foldid = NULL, inner.foldid = NULL,
                     n.outer.folds = 5, n.inner.folds = 10,
                     verbose = 1) {

    # Define the folds
    if (is.null(outer.foldid)) {
        outer.foldid <- cut(1:nrow(d.train), breaks = n.outer.folds, labels = FALSE)
        outer.foldid <- sample(outer.foldid)
        # for (i in 1:n.outer.folds) {
        #     cat("Fold n.", i, ": ")
        #     cat(d.train$IID[foldid == i], "\n")
        # }
    }

    # Find the dimension of the parameter space
    n_params <- length(params.grid)

    nested.cv <- lapply(1:n.outer.folds, function(f) {
        if (verbose > 0)
            message("Outer fold", f, "--- Performing the CV on 4/5ths of the data")

        inner <- CV(d.train[outer.foldid != f], params.grid,
                    foldid = inner.foldid, n.folds = n.inner.folds,
                    fit.best = TRUE, verbose = verbose - 1)
        if (is.null(inner$fit))
            return(NULL)

        if (verbose > 0)
            message("Outer fold", f, "--- Testing on 1/5th of the data")

        outer <- Test(inner$fit, d.train[outer.foldid == f], verbose = verbose - 2)
        inner$fit <- NULL
        outer[, seq_len(n_params)] <- NULL

        list(inner = inner, outer = outer)
    })

    if (any(sapply(nested.cv, is.null))) {
        message(warnings())
        message(str(nested.cv))
        return(NULL)
    }

    # nested.cv is a list of five lists, each with two elements: inner and outer.
    # inner is itself, for each fold, a list of three elements: full, summ and fit.
    inner.cv <- lapply(nested.cv, "[[", "inner")
    outer.cv <- lapply(nested.cv, "[[", "outer")

    inner.cv.full <- lapply(seq_along(inner.cv), function(i) {
        inner.cv[[i]]$full$outer_fold <- i
        inner.cv[[i]]$full
    })
    inner.cv.full <- rbindlist(inner.cv.full)
    inner.cv.summ <- lapply(seq_along(inner.cv), function(i) {
        inner.cv[[i]]$summ$outer_fold <- i
        inner.cv[[i]]$summ
    })
    inner.cv.summ <- rbindlist(inner.cv.summ)

    outer.cv.full <- MergeFolds(outer.cv)
    outer.cv.summ <- AggregateFolds(outer.cv.full)

    # Return a list of two elements. each element is a list of two data.tables
    list(inner.cv = list(full = inner.cv.full, summ = inner.cv.summ),
         outer.cv = list(full = outer.cv.full, summ = outer.cv.summ))
}


CV <- function(d.train, params.grid,
               foldid = NULL, n.folds = 5,
               fit.best = FALSE, verbose = 1) {

    # Define the folds
    if (is.null(foldid)) {
        foldid <- cut(1:nrow(d.train), breaks = n.folds, labels = FALSE)
        foldid <- sample(foldid)
    }

    n_params <- length(params.grid)

    # For each fold, loop through all the parameters, then merge
    # everything at the end. This is more efficient than having the loop
    # over the parameters as the outer one. The loop over the parameters
    # must be done inside the user-specified function.
    kfold.cv <- lapply(1:n.folds, function(f) {
        if (verbose > 0)
            message("Fold", f, "--- Doing CV")

        Train(d.train[foldid != f],
              d.train[foldid == f],
              params.grid,
              verbose = verbose - 1)
    })

    if (any(sapply(kfold.cv, is.null))) {
        message(warnings())
        message(str(kfold.cv))
        return(NULL)
    }

    kfold.cv.full <- MergeFolds(kfold.cv)
    l <- list(full = kfold.cv.full, summ = AggregateFolds(kfold.cv.full))

    if (fit.best) {
        if (verbose > 0)
            message("CV --- Fitting the best model")

        # Train the model with the best parameters for this gene
        kfold.cv.params <- as.list(l$summ[which.min(pred_perf_pval), 1:n_params]) # TODO: use user-specified loss() to find best params
        kfold.cv.fit <- Train(d.train,
                              NULL,
                              kfold.cv.params,
                              return.fit = TRUE,
                              verbose = verbose - 1)
        l <- append(l, list(fit = kfold.cv.fit))
    }

    l
}


# TODO: check that d.test is null when return.fit

# TODO: check also that the number of params is one when return.fit

Train <- function(d.train, d.test = NULL, params.grid,
                  return.fit = FALSE, verbose = 0) {
    if (!length(params.grid$alpha))
        return(NULL)

    tryCatch({
        glm.fit <- glmnet(
            x = as.matrix(d.train[, !c("GENE", "IID", "EXPRESSION")]),
            y = d.train$EXPRESSION,
            alpha = as.numeric(params.grid$alpha),
            lambda = as.numeric(params.grid$lambda),
            standardize = FALSE
        )

        if (!is.null(d.test)) {
            return(Test(glm.fit, d.test, verbose = verbose))
        }

        if (return.fit)
            return(glm.fit)

    }, error = function(e) {
        message(e$message)
        warning(e$message)
        return(NULL)
    })
}


Test <- function(fit, d.test, verbose = 0) {
    glm.pred <- predict(fit, as.matrix(d.test[, !c("GENE", "IID", "EXPRESSION")]))

    pearsons <- lapply(seq_along(fit$lambda), function(p) {
        pearson <- compute_pearson(d.test$EXPRESSION, glm.pred[, p], by = d.test$GENE)
        cbind(
            alpha = params.grid$alpha,
            lambda = fit$lambda[p], # These need to be characters otherwise the merge doesn't work
            data.name = pearson$data.name,
            mse = cvmse(d.test$EXPRESSION, glm.pred[, p]),
            rsq = cvrsq(d.test$EXPRESSION, glm.pred[, p]),
            pearson[, c("parameter.df", "estimate.cor", "p.value")]
        )
    })
    rbindlist(pearsons)
}


Coef <- function(fit) {
    cbind(lambda = fit$lambda, as.data.table(t(as.matrix(coef(fit)))))
}
