# R/cor.ezpz2.R
#' Compute pairwise (weighted) correlations between all columns of a dataset with (bootstrapped) confidence intervals and return matrices, including significance stars.
#'
#' A helper function for psychologists to compute correlations,
#' p-values, significance stars, and confidence intervals (standard and bootstrap),
#' optionally weighted.
#'
#' @param data Data frame of numeric variables
#' @param method Correlation method ("pearson", "spearman", "kendall")
#' @param cis Logical, include confidence intervals?
#' @param weight Optional column name for weights
#' @param nboot Number of bootstrap samples
#' @param correction Correction to be applied c("none","bonferroni","fdr")
#' @param alpha Alpha level (default .05)
#' @return A list with correlation matrices, p-values, stars, and CIs
#' @import boot
#' @import weights
#' @export
cor.ezpz2 <- function(data,
                     method = "pearson",
                     cis = TRUE,
                     weight = NULL,
                     nboot = 1000,
                     correction = "none",
                     alpha = 0.05) {

  correction <- match.arg(correction)

  if (!is.null(weight)) {
    if (!(weight %in% colnames(data))) stop("Weight column not found")
    w <- data[[weight]]
    vars <- setdiff(colnames(data), weight)
  } else {
    vars <- colnames(data)
    w <- NULL
  }

  n <- length(vars)
  n_tests <- n * (n - 1) / 2  # number of unique correlations

  cor_matrix      <- matrix(NA, n, n, dimnames = list(vars, vars))
  p_matrix        <- matrix(NA, n, n, dimnames = list(vars, vars))
  ci_std_matrix   <- matrix("", n, n, dimnames = list(vars, vars))
  ci_boot_matrix  <- matrix("", n, n, dimnames = list(vars, vars))
  stars_matrix    <- matrix("", n, n, dimnames = list(vars, vars))

  assign_stars <- function(p) {
    if (is.na(p)) return("")
    if (p < .001) return("***")
    if (p < .01)  return("**")
    if (p < .05)  return("*")
    ""
  }

  # store unweighted boot objects for Bonferroni CIs
  unweighted_boot_list <- list()
  unweighted_idx_list <- c()

  idx_counter <- 1
  for (i in 1:n) {
    for (j in 1:i) {
      if (i == j) next

      x <- data[[vars[i]]]
      y <- data[[vars[j]]]

      if (is.null(weight)) {
        # unweighted correlation
        test <- stats::cor.test(x, y, method = method)
        r <- unname(test$estimate)
        pval <- test$p.value
        ci_std <- test$conf.int

        boot_fun <- function(d, idx) cor(d[idx,1], d[idx,2], method = method, use = "pairwise")
        boot_obj <- tryCatch(boot::boot(data.frame(x, y), boot_fun, R = nboot), error = function(e) NULL)

        ci_boot <- NA
        if (!is.null(boot_obj)) {
          ci_boot <- tryCatch({
            bci <- boot::boot.ci(boot_obj, type = "perc")
            if (!is.null(bci$percent)) bci$percent[4:5] else NA
          }, error = function(e) NA)
        }

        # store for Bonferroni adjustment
        if (!is.null(boot_obj)) {
          unweighted_boot_list[[length(unweighted_boot_list)+1]] <- boot_obj
          unweighted_idx_list <- c(unweighted_idx_list, idx_counter)
        }

      } else {
        # weighted correlation
        wc <- weights::wtd.cor(x, y, weight = w, mean1 = TRUE, bootse = FALSE)
        r <- unname(wc[1])
        pval <- unname(wc[4])
        ci_std <- NA

        boot_fun <- function(d, idx) {
          tryCatch(
            wtd.cor(d[idx,1], d[idx,2], weight = d[idx,3], mean1 = TRUE)[1],
            error = function(e) NA
          )
        }
        boot_obj <- tryCatch(boot::boot(cbind(x, y, w), boot_fun, R = nboot), error = function(e) NULL)

        ci_boot <- NA
        if (!is.null(boot_obj)) {
          ci_boot <- tryCatch({
            bci <- boot::boot.ci(boot_obj, type = "perc")
            if (!is.null(bci$percent)) bci$percent[4:5] else NA
          }, error = function(e) NA)
        }
      }

      # assign results
      cor_matrix[i,j] <- r
      p_matrix[i,j] <- pval
      stars_matrix[i,j] <- assign_stars(pval)

      if (!is.null(ci_std) && length(ci_std) == 2 && all(is.finite(ci_std))) {
        ci_std_matrix[i,j] <- paste0("[", round(ci_std[1],2), ", ", round(ci_std[2],2), "]")
      }

      if (!any(is.na(ci_boot))) {
        ci_boot_matrix[i,j] <- paste0("[", round(ci_boot[1],2), ", ", round(ci_boot[2],2), "]")
      }

      idx_counter <- idx_counter + 1
    }
  }

  # multiple testing correction (p-values only)
  if (correction != "none") {
    lower_idx <- which(lower.tri(p_matrix))
    pvals <- p_matrix[lower_idx]

    if (correction == "bonferroni") {
      p_adj <- p.adjust(pvals, method = "bonferroni")
      p_matrix[lower_idx] <- p_adj
      stars_matrix[lower_idx] <- sapply(p_adj, assign_stars)

      # Bonferroni-adjusted CIs for unweighted correlations only
      if (cis && length(unweighted_boot_list) > 0) {
        conf_level <- 1 - alpha / n_tests
        for (k in seq_along(unweighted_boot_list)) {
          idx <- unweighted_idx_list[k]
          boot_obj <- unweighted_boot_list[[k]]
          ci_b <- tryCatch({
            bci <- boot::boot.ci(boot_obj, type="perc", conf=conf_level)
            if (!is.null(bci$percent)) paste0("[", round(bci$percent[4],2), ", ", round(bci$percent[5],2), "]") else NA
          }, error = function(e) NA)
          if (!is.na(ci_b)) ci_boot_matrix[idx] <- ci_b
        }
      }

    } else if (correction == "fdr") {
      p_adj <- p.adjust(pvals, method = "fdr")
      p_matrix[lower_idx] <- p_adj
      stars_matrix[lower_idx] <- sapply(p_adj, assign_stars)
      # leave all CIs untouched
    }
  }

  cor_with_stars <- ifelse(is.na(cor_matrix), "", paste0(round(cor_matrix, 2), stars_matrix))

  return(list(
    cor = cor_matrix,
    p = p_matrix,
    stars = stars_matrix,
    cor_with_stars = cor_with_stars,
    ci_std = ci_std_matrix,
    ci_boot = ci_boot_matrix
  ))
}
