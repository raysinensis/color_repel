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
#' @param layer layer to detect color, defaults to first
#' @param out_orig output the original colors as named vector
#' @param out_worst output the worst combination instead of best
#' @param repel_label whether to add centroid labels with ggrepel
#' @param encircle whether to draw geom_encircle by cluster
#' @param encircle_alpha alpha argument passed to geom_encircle
#' @param encircle_expand expand argument passed to geom_encircle
#' @param encircle_shape shape/smoothing argument passed to geom_encircle
#' @param encircle_threshold threshold for removing outliers
#' @param encircle_nmin number of near neighbors for removing outliers
#' @param ... passed to repel_label
#' @examples
#' a <- ggplot2::ggplot(ggplot2::mpg, ggplot2::aes(displ, hwy)) +
#'   ggplot2::geom_point(ggplot2::aes(color = as.factor(cyl)))
#' b <- gg_color_repel(a, col = "colour")
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
                           layer = 1,
                           out_orig = FALSE,
                           out_worst = FALSE,
                           repel_label = FALSE,
                           encircle = FALSE,
                           encircle_alpha = 0.25,
                           encircle_expand = 0.02,
                           encircle_shape = 0.5,
                           encircle_threshold = 0.01,
                           encircle_nmin = 0.01,
                           ...) {
  newcols <- color_repel(g,
    col = col, verbose = verbose,
    downsample = downsample,
    nsamp = nsamp, seed = seed,
    sim = sim, severity = severity,
    autoswitch = autoswitch, layer = layer,
    out_orig = out_orig,
    out_worst = out_worst
  )

  if (autoswitch) {
    col <- check_colour_mapping(g, col = col, autoswitch = autoswitch, layer = layer)
  }
  .f <- paste0("ggplot2:::scale_", col, "_manual")

  labs <- get_labs(g)

  if (all(is.na(labs))) {
    g <- suppressMessages(g + do.call(eval(parse(text = .f)), c(values = list(newcols))))
  } else {
    g <- suppressMessages(g + do.call(eval(parse(text = .f)), c(values = list(newcols), labels = list(labs))))
  }

  if (encircle) {
    dat <- prep_encircle(g, threshold = encircle_threshold, nmin = encircle_nmin, downsample = downsample, seed = seed)
    g <- g + ggalt::geom_encircle(
      data = dat,
      ggplot2::aes(x = x, y = y, fill = group),
      expand = encircle_expand,
      s_shape = encircle_shape,
      alpha = encircle_alpha,
      show.legend = FALSE
    )
    g <- suppressMessages(g + do.call(eval(parse(text = "ggplot2:::scale_fill_manual")), c(values = list(newcols))))
  }

  if (repel_label) {
    g <- label_repel(g, ...)
  }

  g
}
