#' Merge the folds to obtain a table with all the folds.
#' Operates on a list of "pearson" objects returned by Call_RF().
#'
#' @param cv.list A list of "pearson" data.tables.
#' @param n.params The number of parameters.
#'
#' @return A single data.table with all the folds cbinded.
MergeFolds <- function(cv.list) {
    n.params <- which(names(cv.list[[1]]) == "data.name") - 1
    i <- 1 # it may happend that the first fold is useless, so we go on until we find a viable one.
    while (i <= length(cv.list) && length(cv.list[[i]]) < n.params + 6) i <- i + 1
    if (i > length(cv.list))
        print(cv.list) # DEBUG
    cv.table <- cv.list[[i]]
    if (n.params > 0)
        cv.table[, (seq_len(n.params)) := lapply(.SD, as.character), .SDcols = seq_len(n.params)]
    names(cv.table)[(n.params + 2):(n.params + 6)] <- c("mse_1", "rsq_1", "parameter.df_1", "estimate.cor_1", "p.value_1")
    for (k in (i + 1):length(cv.list)) {
        cv.item <- cv.list[[k]]
        if (n.params > 0)
            cv.item[, (1:n.params) := lapply(.SD, as.character), .SDcols = 1:n.params]
        if (length(names(cv.item)) != n.params + 6)
            return(NULL)
        names(cv.item)[(n.params + 2):(n.params + 6)] <- c(paste0("mse_", k), paste0("rsq_", k), paste0("parameter.df_", k), paste0("estimate.cor_", k), paste0("p.value_", k))
        cv.table <- merge(cv.table, cv.item, all.x = TRUE, all.y = TRUE, sort = F)
    }
    cv.table
}

#' Compute the mean for all the folds.
#' Operates on the output of MergeFolds().
#'
#' @param cv.table The output of MergeFolds().
#' @param n.params The number of parameters.
#'
#' @return A single data.table with all the folds aggregated.
AggregateFolds <- function(cv.table) {
    n = length(cv.table)
    n.params <- which(names(cv.table) == "data.name") - 1
    cbind(
        cv.table[, 1:(n.params + 1)],
        cv.table[, .(mse_avg = rowMeans(.SD, na.rm = T),
                     mse_sd = rowSds(.SD)), .SDcols = seq(n.params + 2, n, 5)],
        cv.table[, .(rsq_avg = rowMeans(.SD, na.rm = T),
                     rsq_sd = rowSds(.SD)), .SDcols = seq(n.params + 3, n, 5)],
        cv.table[, .(rho_avg = rowMeans(.SD, na.rm = T),
                     rho_sd = rowSds(.SD),
                     pred_perf_R2 = rowMeans(.SD^2, na.rm = T)), .SDcols = seq(n.params + 5, n, 5)],
        cv.table[, .(rho_zscore = rowZscs(.SD)$zscore,
                     pred_perf_pval = rowZscs(.SD)$pvalue), .SDcols = seq(n.params + 6, n, 5)])
}

#' Compute the mean for all the genes.
#' Operates on the output of AggregateFolds().
#'
#' @param cv.aggr_table The output of AggregateFolds().
#' @param n.params The number of parameters.
#'
#' @return A single data.table with all the rows summarised.
SummariseGenes <- function(cv.aggr_table) {
    n.params <- which(names(cv.aggr_table) == "data.name") - 1
    cv.aggr_table[, .(
        prop_sign = mean(pred_perf_pval < 0.05, na.rm = TRUE),
        mean_rho_avg = mean(rho_avg, na.rm = TRUE),
        sd_rho_avg = stats::sd(rho_avg, na.rm = TRUE),
        mean_pred_perf_pval = mean(pred_perf_pval, na.rm = TRUE),
        sd_pred_perf_pval = stats::sd(pred_perf_pval, na.rm = TRUE),
        mean_rho_avg_of_sign = mean(rho_avg[pred_perf_pval < 0.05], na.rm = TRUE),
        sd_rho_avg_of_sign = stats::sd(rho_avg[pred_perf_pval < 0.05], na.rm = TRUE)
    ), by = eval(names(cv.aggr.table)[1:n.params])]
}


Sds <- function(...) {
    stats::sd(c(...), na.rm = T)
}
MapSds <- function(...) {
    mapply(Sds, ...)
}
#' Compute the standard deviations of the rows of a data frame.
#'
#' @param l A numeric data frame (or a list with equally long elements)
#' @return A list with the row-wise standard deviations.
rowSds <- function(l) {
    do.call("MapSds", l)
}

sumz <- function (p)
{
    na.fail(p)
    p <- as.numeric(p)
    keep <- (p > 0) & (p < 1)
    n <- sum(keep)
    if (n < 2)
        stop("Must have at least two valid p values")
    if (n != length(p))
        warning("Some studies omitted")
    zp <- as.numeric((qnorm(p[keep], lower.tail = FALSE) %*% rep(1, n)) / sqrt(n))
    res <- list(z = zp, p = pnorm(zp, lower.tail = FALSE), validp = p[keep])
    res
}
Zscs <- function(...) {
    pvals <- stats::na.omit(as.numeric(c(...)))
    rho_zscore <- NA
    pred_perf_pval <- NA
    if (length(pvals) > 2) {
        stouffer <- sumz(pvals)
        rho_zscore <- stouffer$z
        pred_perf_pval <- stouffer$p
    }

    c(zscore = rho_zscore,
      pvalue = pred_perf_pval)
}
MapZscs <- function(...) {
    t(mapply(Zscs, ...))
}
#' Compute row-wise Z-scores and p-valuse using Stouffer's method.
#'
#' @param l A numeric data frame of p-values (or a list with equally long
#'   elements)
#' @return A data frame with two columns: rho_zscore and pred_perf_pval
rowZscs <- function(l) {
    as.data.frame(do.call("MapZscs", l))
}

#' compute_pearson
#'
#' @param ytest Vector of true values.
#' @param ypred Vector of predicted values.
#' @param by Factor by which to split the data.
#' @param alternative See cor.test().
#'
#' @return A "pearson" data.table.
compute_pearson <- function(ytest, ypred, by = NULL, alternative = "g") {
    if (is.null(by)) {
        pearson <- t(unlist(cor.test(ytest, ypred, alternative = alternative)))
    } else {
        # Analyse the correlation by groups
        test.split <- split(ytest, as.factor(by))
        pred.split <- split(ypred, as.factor(by))
        pearson <- t(sapply(1:length(test.split), function(r) {
            l <- cor.test(test.split[[r]], pred.split[[r]], alternative = alternative)
            l$data.name <- names(test.split[r])
            unlist(l)
        }))
    }
    pearson <- as.data.table(pearson)
    pearson[, parameter.df := as.numeric(parameter.df)]
    pearson[, estimate.cor := as.numeric(estimate.cor)]
    pearson[, p.value := as.numeric(p.value)]

    pearson
}

#' MSE
#'
#' Mean squared error of the
#'
#' @param ytest A vector of the true response values.
#' @param ypred A vector of the predicted values.
#' @param p The number of predictors used by the model.
#'
#' @return The mean squared error of the model on the hold-out fold.
cvmse <- function(ytest, ypred) {
    mean((ytest - ypred)^2)
}

#' \ifelse{html}{\out{R<sup>2</sup>}}{\eqn{R^2}}
#'
#' @inheritParams cvmse
#'
#' @return The \ifelse{html}{\out{R<sup>2</sup>}}{\eqn{R^2}} of the model on the hold-out fold.
cvrsq <- function(ytest, ypred) {
    SS.total      <- sum((ytest - mean(ytest))^2)
    SS.residual   <- sum((ytest - ypred)^2)

    1 - SS.residual / SS.total
}

#' Fraction of Variability
#'
#' @inheritParams cvmse
#'
#' @return The fraction of variability explained by the model in the hold-out fold.
cvvarfrac <- function(ytest, ypred) {
    SS.total      <- sum((ytest - mean(ytest))^2)
    SS.regression <- sum((ypred - mean(ytest))^2)

    SS.regression / SS.total
}


#' Obtain the regions whose indels are too large
#'
#' @param snps_file The path to vcfrider's SNPs file
#' @param fraction Maximum length of indels relative to the length of the regreg
#'
#' @return A list of regions to be excluded
exclude_from_indel_length <- function(snps_file, fraction) {
    snps <- fread(snps_file,
                  header = FALSE,
                  sep = "\t",
                  col.names = c("REGION", "START", "END", "OV_MUT"))

    snps$CUTOFF <- (snps$END - snps$START) * fraction
    indels <- strsplit(snps$OV_MUT, ",")

    snps$MAX_INS <- sapply(indels, function(i) {
        i <- strsplit(i, "_")
        sum(sapply(i, function(j) {
            if (j[3] == "true" && j[2] > 0)
                as.numeric(j[2])
            else
                0
        }))
    })
    snps$MAX_DEL <- sapply(indels, function(i) {
        i <- strsplit(i, "_")
        sum(sapply(i, function(j) {
            if (j[3] == "true" && j[2] < 0)
                -as.numeric(j[2])
            else
                0
        }))
    })

    # snps[MAX_INS > CUTOFF | MAX_DEL > CUTOFF, c("REGION", "START", "END")]
    snps[MAX_INS > CUTOFF | MAX_DEL > CUTOFF]$REGION
}

#' Remove the version from Ensembl stable IDs
#'
#' @param ensg The gene IDs (possibly not in Ensembl format. Mixed formats are allowed.
#'
#' @return The IDs without the version
strip_ensg_version <- function(ensg) {
    # https://www.ensembl.org/info/genome/stable_ids/index.html
    ifelse(grepl("^ENS.*(E|FM|G|GT|P|R|T)[[:digit:]]{11}(\\.[[:digit:]]+)?$", ensg),
           tstrsplit(ensg, ".", fixed = TRUE)[[1]],
           ensg)
}
