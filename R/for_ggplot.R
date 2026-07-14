#' Add a color-repel scale to a ggplot
#'
#' `scale_color_repel()` can be added to a ggplot with `+`. It calculates a
#' reordered color palette from the plot object, then adds a manual color scale
#' using those optimized colors.
#'
#' @param nsamp how many random sampling color combinations to test, default 50000
#' @param sim passing a colorbind simulation function if needed
#' @param severity severity of the color vision defect, between 0 and 1
#' @param verbose whether to print messages
#' @param downsample downsample when too many datapoints are present, or use chull
#' @param polychrome_recolor whether to replace the original colors with polychrome creation
#' @param seed sampling randomization seed
#' @param col colour or fill in ggplot
#' @param autoswitch try to switch between colour and fill automatically
#' @param layer layer to detect color, defaults to first
#' @param out_worst output the worst combination instead of best
#' @param ... passed to [ggplot2::scale_color_manual()] or
#'   [ggplot2::scale_fill_manual()]
#' @examples
#' a <- ggplot2::ggplot(ggplot2::mpg, ggplot2::aes(displ, hwy)) +
#'   ggplot2::geom_point(ggplot2::aes(color = as.factor(cyl)))
#' b <- a + scale_color_repel()
#' c <- ggplot2::ggplot(ggplot2::mpg, ggplot2::aes(class, fill = as.factor(cyl))) +
#'   ggplot2::geom_bar() +
#'   scale_fill_repel()
#' @return object added to a ggplot with `+`
#' @export
scale_color_repel <- function(nsamp = 50000,
                              sim = NULL,
                              severity = 0.5,
                              verbose = FALSE,
                              downsample = 5000,
                              polychrome_recolor = FALSE,
                              seed = 34,
                              col = "colour",
                              autoswitch = TRUE,
                              layer = 1,
                              out_worst = FALSE,
                              ...) {
  structure(
    list(
      repel_params = list(
        nsamp = nsamp,
        sim = sim,
        severity = severity,
        verbose = verbose,
        downsample = downsample,
        polychrome_recolor = polychrome_recolor,
        seed = seed,
        col = col,
        autoswitch = autoswitch,
        layer = layer,
        out_worst = out_worst
      ),
      scale_params = list(...)
    ),
    class = "scale_color_repel"
  )
}

#' @rdname scale_color_repel
#' @export
scale_colour_repel <- scale_color_repel

#' @rdname scale_color_repel
#' @export
scale_fill_repel <- function(nsamp = 50000,
                             sim = NULL,
                             severity = 0.5,
                             verbose = FALSE,
                             downsample = 5000,
                             polychrome_recolor = FALSE,
                             seed = 34,
                             col = "fill",
                             autoswitch = TRUE,
                             layer = 1,
                             out_worst = FALSE,
                             ...) {
  scale_color_repel(
    nsamp = nsamp,
    sim = sim,
    severity = severity,
    verbose = verbose,
    downsample = downsample,
    polychrome_recolor = polychrome_recolor,
    seed = seed,
    col = col,
    autoswitch = autoswitch,
    layer = layer,
    out_worst = out_worst,
    ...
  )
}

#' @importFrom ggplot2 ggplot_add
#' @method ggplot_add scale_color_repel
#' @export
ggplot_add.scale_color_repel <- function(object, plot, ...) {
  ggbuild <- ggplot2::ggplot_build(plot)
  col <- check_colour_mapping(
    plot,
    col = object$repel_params[["col"]],
    autoswitch = object$repel_params[["autoswitch"]],
    layer = object$repel_params[["layer"]],
    ggbuild = ggbuild
  )

  repel_params <- object$repel_params
  repel_params[["col"]] <- col
  repel_params[["autoswitch"]] <- FALSE

  newcols <- do.call(
    color_repel,
    c(
      list(
        g = plot,
        ggbuild = ggbuild
      ),
      repel_params
    )
  )

  scale_params <- object$scale_params
  scale_params[["values"]] <- newcols
  labs <- get_labs(plot, ggbuild = ggbuild)
  if (!all(is.na(labs)) && is.null(scale_params[["labels"]])) {
    scale_params[["labels"]] <- labs
  }

  scale_fun <- switch(col,
    fill = ggplot2::scale_fill_manual,
    colour = ggplot2::scale_color_manual,
    color = ggplot2::scale_color_manual,
    ggplot2::scale_color_manual
  )

  plot + do.call(scale_fun, scale_params)
}
