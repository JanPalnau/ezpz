#' Fit Linear Model with Standardized Variables
#'
#' Fits a linear model (`lm`) using standardized (z-scored) numeric variables
#' from the supplied formula. Standardization is performed within the function.
#'
#' @param formula A model formula (e.g., y ~ x1 + x2)
#' @param data A data.frame containing the variables in the formula
#' @param center Logical; whether to center variables (default = TRUE)
#' @param scale Logical; whether to scale variables (default = TRUE)
#' @param return_data Logical; return standardized data as attribute (default = FALSE)
#'
#' @return An object of class `lm`. If `return_data = TRUE`, the standardized
#' data are attached as an attribute `"scaled_data"`.
#'
#' @details
#' Only numeric variables in the formula are standardized. Factors are left unchanged.
#'
#' @examples
#' df <- data.frame(
#'   y = rnorm(100),
#'   x1 = rnorm(100),
#'   x2 = rnorm(100),
#'   group = factor(sample(letters[1:2], 100, TRUE))
#' )
#'
#' model <- lm.beta(y ~ x1 + x2 + group, data = df)
#' summary(model)
#'
#' @export
lm.beta <- function(formula, data, center = TRUE, scale = TRUE, return_data = FALSE) {

  # --- Input checks ---
  if (!inherits(formula, "formula")) {
    stop("`formula` must be a valid formula.")
  }

  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame.")
  }

  # --- Extract variables used in formula ---
  vars <- all.vars(formula)

  # Check variables exist
  missing_vars <- setdiff(vars, names(data))
  if (length(missing_vars) > 0) {
    stop("Variables not found in data: ", paste(missing_vars, collapse = ", "))
  }

  # --- Copy data to avoid modifying original ---
  d <- data

  # --- Identify numeric variables only ---
  is_num <- sapply(d[vars], is.numeric)
  num_vars <- vars[is_num]

  # --- Standardize numeric variables ---
  if (length(num_vars) > 0) {
    d[num_vars] <- scale(d[num_vars], center = center, scale = scale)
  }

  # --- Fit model ---
  fit <- lm(formula, data = d)

  # --- Optionally attach standardized data ---
  if (return_data) {
    attr(fit, "scaled_data") <- d
  }

  return(fit)
}
