# R/cor.ezpz.R

#' Compute pairwise correlations with confidence intervals
#'
#' A helper function for psychologists to compute correlation matrices with
#' p-values, significance stars, and confidence intervals (standard and bootstrap),
#' optionally weighted.
#'
#' @param data Data frame of numeric variables
#' @param method Correlation method ("pearson", "spearman", "kendall")
#' @param cis Logical, include confidence intervals?
#' @param weight Optional column name for weights
#' @param nboot Number of bootstrap samples
#' @return A list with correlation matrices, p-values, stars, and CIs
#' @import boot
#' @import weights
#' @export
cor.ezpz <- function(data,
                     method = "pearson",
                     cis = TRUE,
                     weight = NULL,
                     nboot = 1000) {
  
  # Separate variables from weight column
  if (!is.null(weight)) {
    w <- data[[weight]]
    vars <- setdiff(colnames(data), weight)
  } else {
    vars <- colnames(data)
    w <- NULL
  }
  
  n <- length(vars)
  
  # matrices
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
  
  # loop over pairs
  for (i in 1:n) {
    for (j in 1:i) {
      if (i == j) next
      
      x <- data[[vars[i]]]
      y <- data[[vars[j]]]
      
      if (is.null(weight)) {
        test <- suppressWarnings(cor.test(x, y, method = method))
        r <- unname(test$estimate)
        pval <- test$p.value
        ci_std <- test$conf.int
        
        boot_fun <- function(d, idx) {
          cor(d[idx,1], d[idx,2], method = method, use = "pairwise")
        }
        boot_obj <- try(boot(data.frame(x, y), boot_fun, R = nboot), silent = TRUE)
        ci_boot <- tryCatch({
          bci <- boot.ci(boot_obj, type = "perc")
          if (!is.null(bci$percent)) bci$percent[4:5] else NA
        }, error = function(e) NA)
        
      } else {
        wc <- suppressWarnings(wtd.cor(x, y, weight = w, mean1 = TRUE, bootse = FALSE))
        r <- unname(wc[1])
        pval <- unname(wc[4])
        
        ci_std <- NA
        boot_fun <- function(d, idx) {
          wtd.cor(d[idx,1], d[idx,2], weight = d[idx,3], mean1 = TRUE)[1]
        }
        boot_obj <- try(boot(cbind(x, y, w), boot_fun, R = nboot), silent = TRUE)
        ci_boot <- tryCatch({
          bci <- boot.ci(boot_obj, type = "perc")
          if (!is.null(bci$percent)) bci$percent[4:5] else NA
        }, error = function(e) NA)
      }
      
      cor_matrix[i,j] <- r
      p_matrix[i,j]   <- pval
      stars_matrix[i,j] <- assign_stars(pval)
      
      if (!is.null(ci_std) && all(is.finite(ci_std))) {
        ci_std_matrix[i,j] <- paste0("[", round(ci_std[1],2), ", ", round(ci_std[2],2), "]")
      }
      if (!is.null(ci_boot) && all(is.finite(ci_boot))) {
        ci_boot_matrix[i,j] <- paste0("[", round(ci_boot[1],2), ", ", round(ci_boot[2],2), "]")
      }
    }
  }
  
  cor_with_stars <- ifelse(is.na(cor_matrix), "", paste0(round(cor_matrix, 2), stars_matrix))
  
  return(list(cor = cor_matrix,
              p = p_matrix,
              stars = stars_matrix,
              cor_with_stars = cor_with_stars,
              ci_std = ci_std_matrix,
              ci_boot = ci_boot_matrix))
}
