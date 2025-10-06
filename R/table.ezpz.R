#' Export a table as a formatted Word document
#'
#' Quickly export a data frame or matrix to a Word document with consistent
#' formatting using the \pkg{flextable} and \pkg{officer} packages. Ideal for
#' generating tables for manuscripts, reports, or presentations.
#'
#' @param table A data frame or matrix to be turned into a formatted table.
#' @param name File name for the output Word document (default: "example_table.docx").
#' @param font_name Font family to use for all text (default: "Times New Roman").
#' @param font_size Font size for all table text (default: 8).
#' @param table_note Optional text note placed below the table.
#' @param align_all Default alignment for all table parts ("left", "center", or "right").
#' @param align_header Alignment for header text (default: "center").
#' @param column_names Optional named vector to rename columns (e.g., c("mean" = "Mean", "sd" = "SD")).
#'
#' @return A \code{flextable} object (and writes a .docx file to disk).
#' @import flextable
#' @import officer
#' @export
table.ezpz <- function(
    table = NULL,
    name = "example_table.docx",
    font_name = "Times New Roman",
    font_size = 8,
    table_note = "Example Note",
    align_all = "left",
    align_header = "center",
    column_names = NULL
) {
  if (is.null(table)) stop("You must provide a table or data frame.")
  
  # Add row names as a "variable" column for convenience
  table <- cbind(variable = rownames(table), table)
  
  ft <- flextable::flextable(as.data.frame(table))
  ft <- flextable::font(ft, fontname = font_name, part = "all")
  ft <- flextable::fontsize(ft, size = font_size, part = "all") 
  ft <- flextable::align(ft, align = align_all, part = "all")
  ft <- flextable::align(ft, align = align_header, part = "header")
  
  # Rename columns if requested
  if (!is.null(column_names)) {
    ft <- do.call(flextable::set_header_labels, c(list(ft), column_names))
  }
  
  # Add a note under the table
  formatted_note <- flextable::as_paragraph(flextable::as_i("Notes. "), table_note)
  ft <- flextable::add_footer_lines(ft, values = formatted_note)
  ft <- flextable::align(ft, align = "left", part = "footer")
  
  # Adjust width and structure
  ft <- flextable::set_table_properties(ft, width = 0.8)
  ft <- flextable::separate_header(ft)
  
  # Write to Word
  doc <- officer::read_docx()
  doc <- officer::body_add_flextable(doc, value = ft)
  print(doc, target = name)
  
  invisible(ft)
}
