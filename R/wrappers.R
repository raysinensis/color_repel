#' Wrapper to reorder ggplot colors to maximize color differences in space
#' @param g ggplot plot object
#' @param col colour or fill in ggplot
#' @return new ggplot object
#' @export
gg_color_repel <- function(g, col = "colour") {
  .f <- paste0("scale_", col, "_manual")
  newcols <<- color_repel(g, col = col)
  g + do.call(.f, c(values = list(newcols)))
}