#' Bootstrap Confidence Intervals for Linear Model Coefficients
#'
#' Computes bootstrapped confidence intervals (percentile method) for
#' standardized regression coefficients using `boot::boot()`.
#'
#' @param data A data.frame
#' @param formula A model formula
#' @param R Number of bootstrap resamples (default = 5000)
#' @param conf Confidence level (default = 0.95)
#' @param adjust Logical; whether to apply Bonferroni correction (default = FALSE)
#' @param include_intercept Logical; include intercept (default = TRUE)
#' @param ... Additional arguments passed to `lm_boot`
#'
#' @return A data.frame with term, lower, and upper CI bounds
#'
#' @examples
#' \dontrun{
#' lm_boot_ci(data = df, formula = y ~ x1 + x2)
#' }
#'
#' @export
lm.boot.ci <- function(data,
                       formula,
                       R = 5000,
                       conf = 0.95,
                       adjust = FALSE,
                       include_intercept = TRUE,
                       ...) {

  # --- Input checks ---
  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame.")
  }

  if (!inherits(formula, "formula")) {
    stop("`formula` must be a formula.")
  }

  # --- Run bootstrap ---
  boot_out <- boot::boot(
    data = data,
    statistic = lm_boot,
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

  # --- Compute CIs ---
  ci_mat <- t(sapply(seq_len(n_coefs), function(i) {
    ci <- boot::boot.ci(
      boot_out,
      type = "perc",
      index = i,
      conf = conf_adj
    )

    # Handle potential NULL (rare but possible)
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

  return(out)
}
