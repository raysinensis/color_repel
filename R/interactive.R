#' Prepare ggplot object to ggplotly-compatible layer and image layer
#' @param g ggplot plot object
#' @param width plot width
#' @param height plot height
#' @param filename temp file location for saving image
#' @param draw_box if a colored background should be included
#' @param background_alpha alpha value of background image
#' @param ... arguments passed to gg_color_repel
#' @examples
#' a <- ggplot2::ggplot(ggplot2::mpg, ggplot2::aes(displ, hwy)) +
#'   ggplot2::geom_point(ggplot2::aes(color = as.factor(cyl)))
#' new_colors <- color_repel(a)
#' b <- ggplotly_background(a, filename = NULL)
#' @return plotly object with background image of layers unsupported by plotly
#' @export
ggplotly_background <- function(g, width = 5, height = 5, filename = "temp.png", draw_box = NULL, background_alpha = 1, ...) {
  a <- g
  b <- gg_color_repel(a, nudge_x = 2, nudge_y = 2, force = 10, ...)
  c <- ggplot2::ggplot_build(a)
  xmin <- min(c$data[[1]]$x)
  xmax <- max(c$data[[1]]$x)
  ymin <- min(c$data[[1]]$y)
  ymax <- max(c$data[[1]]$y)
  
  if ((is.null(filename))) {
    return(plotly::ggplotly(a))
  }
  tempbg <- crop_background(save_background(prep_background(remove_geom(b), xmin, xmax, ymin, ymax, draw_box),
                                            filename = filename))
  ggplotly_withbg(b, xmin, xmax, ymin, ymax, filename = tempbg, alpha = background_alpha)
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
  g <- g + ggplot2::labs(x=NULL, y=NULL) + 
    ggplot2::coord_cartesian(expand = F, xlim = c(xmin * 1.0, xmax * 1.0), ylim = c(ymin * 1.0, ymax * 1.0), clip = "off") +
    ggplot2::theme(axis.line = ggplot2::element_blank(), 
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
          plot.margin = ggplot2::unit(c(-1, -1, -1.5, -1.5), "lines"))
  if (is.null(draw_box)) {
    g + ggplot2::theme(plot.background = ggplot2::element_blank())
  } else {
    g + ggplot2::theme(plot.background = ggplot2::element_rect(fill=draw_box))
  }
}

save_background <- function(g, filename = "temp.png", width = 5, height = 5) {
  ggplot2::ggsave(filename, g, width = width, height = height)
  return(filename)
}

crop_background <- function(filename = "temp.png") {
  knitr::plot_crop(filename)
  return(filename)
}

ggplotly_withbg <- function(g, xmin, xmax, ymin, ymax, filename = "temp.png", width = 5, height = 5, alpha = 1) {
  p <- plotly::ggplotly(g, width = width * 100, height = height * 100)
  
  p <- plotly::layout(p, autosize = F, 
                      margin = list(l = 0, r = 0, b = 0, t = 0, pad = 0, autoexpand = T),
                      scene = list(aspectmode = "data"),
                      xaxis = list(autorange = F, range = list(xmin*1.0, xmax*1.0)),
                      yaxis = list(autorange = F, range = list(ymin*1.0, ymax*1.0), scaleanchor= 'x', scaleratio = (xmax - xmin)/(ymax - ymin)),
                      images = list(#source = base64enc::dataURI(file = filename), 
                                    source = raster2uri(as.raster(png::readPNG(filename))),
                                    opacity = alpha,
                                    x = xmin, y = ymin,
                                    sizex = xmax - xmin, sizey = ymax - ymin,
                                    xref = "x1", yref= "y1",
                                    sizing = "stretch",
                                    xanchor = "left", yanchor = "bottom", layer = "below"),
                      showlegend = FALSE)
  p <- plotly::config(p, displayModeBar = F)
  p
}