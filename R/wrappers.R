#' Wrapper to reorder ggplot colors to maximize color differences in space
#' @param g ggplot plot object
#' @param col colour or fill in ggplot
#' @param verbose whether to print messages
#' @param downsample downsample when too many datapoints are present
#' @param seed sampling randomization seed
#' @return new ggplot object
#' @export
gg_color_repel <- function(g, col = "colour", verbose = FALSE, downsample = 10000, seed = 34) {
  .f <- paste0("scale_", col, "_manual")
  newcols <- color_repel(g, col = col, verbose = verbose, downsample = downsample, seed = seed)
  labs <- get_labs(g)
  if (all(is.na(labs))) {
    suppressMessages(g + do.call(.f, c(values = list(newcols))))
  } else {
    suppressMessages(g + do.call(.f, c(values = list(newcols), labels = list(labs))))
  }
}