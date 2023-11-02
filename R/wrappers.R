#' Wrapper to reorder ggplot colors to maximize color differences in space
#' @param g ggplot plot object
#' @param col colour or fill in ggplot
#' @param sim passing a colorbind simulation function if needed
#' @param severity severity of the color vision defect, between 0 and 1
#' @param verbose whether to print messages
#' @param downsample downsample when too many datapoints are present
#' @param nsamp how many random sampling color combinations to test, default 50000
#' @param seed sampling randomization seed
#' @param autoswitch try to switch between colour and fill automatically
#' @param out_worst output the worst combination instead of best
#' @return new ggplot object
#' @export
gg_color_repel <- function(g = ggplot2::last_plot(),
                           col = "colour",
                           sim = NULL,
                           severity = 0.5,
                           verbose = FALSE,
                           downsample = 5000,
                           nsamp = 50000,
                           seed = 34,
                           autoswitch = TRUE,
                           out_worst = FALSE) {
  newcols <- color_repel(g,
    col = col, verbose = verbose,
    downsample = downsample, nsamp = nsamp, seed = seed,
    sim = sim, severity = severity,
    autoswitch = autoswitch,
    out_worst = out_worst
  )

  if (autoswitch) {
    col <- check_colour_mapping(g, col = col, autoswitch = autoswitch)
  }
  .f <- paste0("scale_", col, "_manual")

  labs <- get_labs(g)

  if (all(is.na(labs))) {
    suppressMessages(g + do.call(.f, c(values = list(newcols))))
  } else {
    suppressMessages(g + do.call(.f, c(values = list(newcols), labels = list(labs))))
  }
}
