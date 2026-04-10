#' Bootstrap Confidence Intervals for Linear Model Coefficients
#'
#' Computes bootstrapped confidence intervals (percentile method) for
#' standardized regression coefficients using `boot::boot()`.
#' Optionally includes bootstrapped R-squared.
#'
#' @param data A data.frame
#' @param formula A model formula
#' @param boot_function Function used in bootstrapping (must return coefficients)
#' @param R Number of bootstrap resamples (default = 5000)
#' @param conf Confidence level (default = 0.95)
#' @param adjust Logical; whether to apply Bonferroni correction (default = FALSE)
#' @param include_intercept Logical; include intercept (default = TRUE)
#' @param include_r2 Logical; whether to include bootstrapped R² (default = FALSE)
#' @param ... Additional arguments passed to `boot_function`
#'
#' @return A data.frame with term, lower, and upper CI bounds
#'
#' @examples
#' \dontrun{
#' lm.boot.ci(data = df, formula = y ~ x1 + x2, boot_function = lm_boot)
#' }
#'
#' @export
lm.boot.ci <- function(data,
                       formula,
                       boot_function,
                       R = 5000,
                       conf = 0.95,
                       adjust = FALSE,
                       include_intercept = TRUE,
                       include_r2 = FALSE,
                       ...) {

  # --- Input checks ---
  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame.")
  }

  if (!inherits(formula, "formula")) {
    stop("`formula` must be a formula.")
  }

  # --- Run bootstrap for coefficients ---
  boot_out <- boot::boot(
    data = data,
    statistic = boot_function,
    R = R,
    formula = formula,
    include_intercept = include_intercept,
    ...
  )

  # --- Number of coefficients ---
  coef_names <- names(boot_out$t0)
  n_coefs <- length(coef_names)

  # --- Adjust confidence level if requested ---
  if (adjust) {
    conf_adj <- 1 - (1 - conf) / n_coefs
  } else {
    conf_adj <- conf
  }

  # --- Compute coefficient CIs ---
  ci_mat <- t(sapply(seq_len(n_coefs), function(i) {
    ci <- boot::boot.ci(
      boot_out,
      type = "perc",
      index = i,
      conf = conf_adj
    )

    if (is.null(ci$percent)) {
      return(c(NA, NA))
    }

    ci$percent[4:5]
  }))

  # --- Build output table ---
  out <- data.frame(
    term = coef_names,
    CI_lower = ci_mat[, 1],
    CI_upper = ci_mat[, 2],
    row.names = NULL
  )

  # --- Optionally compute R² ---
  if (include_r2) {

    boot_r2_fn <- function(data, indices, formula) {
      d <- data[indices, ]
      fit <- lm(formula, data = d)
      summary(fit)$r.squared
    }

    boot_r2 <- boot::boot(
      data = data,
      statistic = boot_r2_fn,
      R = R,
      formula = formula
    )

    probs <- c((1 - conf) / 2, 1 - (1 - conf) / 2)
    r2_ci <- quantile(boot_r2$t, probs = probs, na.rm = TRUE)

    r2_row <- data.frame(
      term = "R2",
      CI_lower = r2_ci[1],
      CI_upper = r2_ci[2]
    )

    out <- rbind(out, r2_row)
  }

  return(out)
}
