#' Prepare ggplot object to ggplotly-compatible layer and image layer
#' @param g ggplot plot object
#' @param repel_color whether to rearrange colors
#' @param repel_label whether to add centroid labels with ggrepel
#' @param encircle whether to draw geom_encircle by cluster
#' @param width plot width
#' @param height plot height
#' @param filename temp file location for saving image
#' @param draw_box if a colored background should be included
#' @param background if specified, use this ggplot object or file as background instead
#' @param background_alpha alpha value of background image
#' @param use_cairo whether to use cairo for saving plots, maybe needed for certain ggplot extensions
#' @param label_lim whether to limit labels to avoid edge fraction
#' @param ggbuild already built ggplot_built object if available
#' @param crop whether to call cropping of the background image to remove whitespace
#' @param size_nudge slight image size adjustment, default to none
#' @param ... arguments passed to gg_color_repel
#' @examples
#' a <- ggplot2::ggplot(ggplot2::mpg, ggplot2::aes(displ, hwy)) +
#'   ggplot2::geom_point(ggplot2::aes(color = as.factor(cyl)))
#' new_colors <- color_repel(a)
#' b <- ggplotly_background(a, filename = NULL)
#' @return plotly object with background image of layers unsupported by plotly
#' @export
ggplotly_background <- function(g, repel_color = TRUE, repel_label = TRUE, encircle = FALSE,
                                width = 5, height = 5,
                                filename = "temp.png", draw_box = NULL,
                                background = NULL, background_alpha = 1, use_cairo = FALSE, label_lim = 0.05,
                                ggbuild = NULL, crop = TRUE, size_nudge = 0,
                                ...) {
  a <- g
  if (is.null(ggbuild)) {
    c <- ggplot2::ggplot_build(a)
  } else {
    c <- ggbuild
  }

  xmin <- min(c$data[[1]]$x)
  xmax <- max(c$data[[1]]$x)
  ymin <- min(c$data[[1]]$y)
  ymax <- max(c$data[[1]]$y)
  b <- gg_color_repel(a,
    out_orig = !repel_color, repel_label = repel_label, encircle = encircle,
    force = 10, xlim = expand_lims(xmin, xmax, -label_lim), ylim = expand_lims(ymin, ymax, -label_lim),
    ...
  )

  if ((is.null(filename))) {
    return(plotly::ggplotly(a))
  }
  if (is.null((background))) {
    tempbg <- save_background(prep_background(remove_geom(b), 
                                              xmin, xmax, ymin, ymax, draw_box),
                              filename = filename, use_cairo = use_cairo)
    if (crop) {
      tempbg <- crop_background(tempbg)
    }
  } else {
    if (!("character" %in% class(background))) {
      tempbg <- save_background(prep_background(background, 
                                                expand_lims(xmin, xmax, label_lim)[1],
                                                expand_lims(xmin, xmax, label_lim)[2],
                                                expand_lims(ymin, ymax, label_lim)[1],
                                                expand_lims(ymin, ymax, label_lim)[2],
                                                draw_box),
                                filename = filename, use_cairo = use_cairo)
      if (crop) {
        tempbg <- crop_background(tempbg)
      }
    } else {
      tempbg <- background
    }
  }

  ggplotly_withbg(b, xmin, xmax, ymin, ymax, filename = tempbg, alpha = background_alpha, size_nudge = size_nudge)
}

remove_geom <- function(g, layer = 1) {
  g2 <- unserialize(serialize(g, NULL))

  if ("patchwork" %in% class(g2)) {
    g2[[1]]$layers[[layer]]$aes_params$alpha <- 0.01
  } else {
    g2$layers[[layer]]$aes_params$alpha <- 0.01
  }
  g2
}

prep_background <- function(g, xmin, xmax, ymin, ymax, draw_box = NULL) {
  g <- g + ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::coord_cartesian(expand = F, xlim = expand_lims(xmin, xmax, 0), ylim = expand_lims(ymin, ymax, 0), clip = "off") +
    ggplot2::theme(
      axis.line = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(), axis.text.y = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(), axis.title.y = ggplot2::element_blank(),
      axis.ticks.length = ggplot2::unit(0, "lines"), axis.ticks.margin = ggplot2::unit(0, "lines"),
      legend.position = "none",
      panel.background = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.margin = ggplot2::unit(0, "lines"),
      plot.title = ggplot2::element_blank(),
      plot.margin = ggplot2::unit(c(-1, -1, -1.5, -1.5), "lines")
    )

  if (is.null(draw_box)) {
    g + ggplot2::theme(plot.background = ggplot2::element_blank())
  } else {
    g + ggplot2::theme(plot.background = ggplot2::element_rect(fill = draw_box))
  }
}

save_background <- function(g, filename = "temp.png", width = 5, height = 5, use_cairo = FALSE) {
  if (!use_cairo) {
    ggplot2::ggsave(filename, g, width = width, height = height)
  } else {
    ggplot2::ggsave(filename, g, width = width, height = height, type = "cairo")
  }
  return(filename)
}

crop_background <- function(filename = "temp.png") {
  knitr::plot_crop(filename)
  return(filename)
}

ggplotly_withbg <- function(g, xmin, xmax, ymin, ymax, filename = "temp.png",
                            width = 5, height = 5, alpha = 1, legend = FALSE,
                            size_nudge = 0) {
  xmin <- expand_lims(xmin, xmax, size_nudge)[1]
  xmax <- expand_lims(xmin, xmax, size_nudge)[2]
  ymin <- expand_lims(ymin, ymax, size_nudge)[1]
  ymax <- expand_lims(ymin, ymax, size_nudge)[2]
  
  p <- plotly::ggplotly(g, width = width * 100, height = height * 100)
  p <- plotly::layout(p,
    autosize = F,
    margin = list(l = 0, r = 0, b = 0, t = 0, pad = 0, autoexpand = T),
    scene = list(aspectmode = "data"),
    xaxis = list(autorange = F, range = as.list(expand_lims(xmin, xmax))),
    yaxis = list(autorange = F, range = as.list(expand_lims(ymin, ymax)), scaleanchor = "x", scaleratio = (xmax - xmin) / (ymax - ymin)),
    images = list(
      source = plotly::raster2uri(grDevices::as.raster(png::readPNG(filename))),
      opacity = alpha,
      x = (xmax + xmin)/2, y = (ymax + ymin)/2,
      sizex = xmax - xmin, sizey = ymax - ymin,
      xref = "x1", yref = "y1",
      sizing = "stretch",
      xanchor = "center", yanchor = "middle", layer = "below"
    ),
    showlegend = legend
  )
  p <- plotly::config(p, displayModeBar = F)
  p
}
