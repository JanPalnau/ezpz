#' Adjust p-values from a Linear Model
#'
#' Extracts p-values from an `lm` object and applies multiple-testing correction.
#'
#' @param model An object of class `lm`
#' @param method Adjustment method passed to `p.adjust` (e.g., "bonferroni", "holm", "BH")
#' @param include_intercept Logical; whether to include intercept in adjustment (default = FALSE)
#' @param terms Optional; character vector of coefficient names OR numeric indices to include
#' @param return_table Logical; if TRUE, returns a data.frame with estimates and adjusted p-values
#'
#' @return A named vector of adjusted p-values, or a data.frame if `return_table = TRUE`
#'
#' @examples
#' model <- lm(mpg ~ wt + hp, data = mtcars)
#'
#' # Adjust all predictors
#' lm_adjust(model, method = "holm")
#'
#' # Adjust only specific terms
#' lm_adjust(model, terms = c("wt", "hp"))
#'
#' @export
lm.p.adjust <- function(model,
                      method = "bonferroni",
                      include_intercept = FALSE,
                      terms = NULL,
                      return_table = FALSE) {

  # --- Input checks ---
  if (!inherits(model, "lm")) {
    stop("`model` must be an object of class 'lm'.")
  }

  # --- Extract coefficient table ---
  coefs <- summary(model)$coefficients

  # --- Handle intercept ---
  if (!include_intercept) {
    coefs <- coefs[rownames(coefs) != "(Intercept)", , drop = FALSE]
  }

  # --- Subset terms if specified ---
  if (!is.null(terms)) {

    if (is.character(terms)) {
      missing_terms <- setdiff(terms, rownames(coefs))
      if (length(missing_terms) > 0) {
        stop("Terms not found in model: ", paste(missing_terms, collapse = ", "))
      }
      coefs <- coefs[terms, , drop = FALSE]

    } else if (is.numeric(terms)) {
      if (any(terms > nrow(coefs))) {
        stop("Numeric indices exceed number of coefficients.")
      }
      coefs <- coefs[terms, , drop = FALSE]

    } else {
      stop("`terms` must be NULL, character vector, or numeric indices.")
    }
  }

  # --- Extract p-values ---
  p_vals <- coefs[, "Pr(>|t|)"]

  # --- Adjust p-values ---
  p_adj <- p.adjust(p_vals, method = method)

  # --- Return table if requested ---
  if (return_table) {
    out <- data.frame(
      term = rownames(coefs),
      estimate = coefs[, "Estimate"],
      std_error = coefs[, "Std. Error"],
      statistic = coefs[, "t value"],
      p_value = p_vals,
      p_adjusted = p_adj,
      row.names = NULL
    )
    return(out)
  }

  # --- Otherwise return named vector ---
  return(p_adj)
}
