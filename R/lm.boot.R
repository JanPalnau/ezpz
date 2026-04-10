#' Bootstrap Function for Standardized Linear Model Coefficients
#'
#' Internal helper function for use with `boot::boot()` to obtain
#' bootstrapped (standardized) regression coefficients.
#'
#' @param data A data.frame containing the variables
#' @param indices A vector of row indices provided by `boot`
#' @param formula A model formula (e.g., y ~ x1 + x2)
#' @param center Logical; whether to center variables (default = TRUE)
#' @param scale Logical; whether to scale variables (default = TRUE)
#' @param include_intercept Logical; whether to return intercept (default = TRUE)
#'
#' @return A numeric vector of regression coefficients
#'
#' @details
#' This function is designed to be passed to `boot::boot()` as the `statistic`
#' argument. Standardization is performed within each bootstrap sample.
#'
#' @examples
#' \dontrun{
#' library(boot)
#' boot_out <- boot(data = df, statistic = lm_boot, R = 1000,
#'                  formula = y ~ x1 + x2)
#' }
#'
#' @export
lm.boot <- function(data,
                    indices,
                    formula,
                    center = TRUE,
                    scale = TRUE,
                    include_intercept = TRUE) {

  # --- Input checks ---
  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame.")
  }

  if (!inherits(formula, "formula")) {
    stop("`formula` must be a valid formula.")
  }

  # --- Resample data ---
  d <- data[indices, , drop = FALSE]

  # --- Extract variables ---
  vars <- all.vars(formula)

  missing_vars <- setdiff(vars, names(d))
  if (length(missing_vars) > 0) {
    stop("Variables not found in data: ", paste(missing_vars, collapse = ", "))
  }

  # --- Standardize numeric variables only ---
  is_num <- sapply(d[vars], is.numeric)
  num_vars <- vars[is_num]

  if (length(num_vars) > 0) {
    d[num_vars] <- scale(d[num_vars], center = center, scale = scale)
  }

  # --- Fit model ---
  fit <- lm(formula, data = d)

  coefs <- coef(fit)

  # --- Optionally drop intercept ---
  if (!include_intercept) {
    coefs <- coefs[names(coefs) != "(Intercept)"]
  }

  return(coefs)
}
