remove_geom <- function(g, layer = 1) {
  a <- check_patchwork(g)
  a <- unserialize(serialize(a, NULL))
  a$layers[[layer]]$aes_params$alpha <- 0.1
  a
}

prep_background <- function(g, xmin, xmax, ymin, ymax, draw_box = NULL) {
  g <- g + labs(x=NULL, y=NULL) + 
    coord_cartesian(expand = F, xlim = c(xmin * 1.0, xmax * 1.0), ylim = c(ymin * 1.0, ymax * 1.0), clip = "off") +
    theme(axis.line = element_blank(), 
          axis.text.x = element_blank(), axis.text.y = element_blank(),
          axis.ticks = element_blank(), 
          axis.title.x = element_blank(), axis.title.y = element_blank(), 
          axis.ticks.length = unit(0, "lines"), axis.ticks.margin = unit(0, "lines"), 
          legend.position = "none", 
          panel.background = element_blank(), 
          panel.border = element_blank(), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.margin = unit(0, "lines"), 
          plot.title = element_blank(), 
          plot.margin = unit(c(-1, -1, -1.5, -1.5), "lines"))
  if (is.null(draw_box)) {
    g + theme(plot.background = element_blank())
  } else {
    g + theme(plot.background = element_rect(fill=draw_box))
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

ggplotly_withbg <- function(g, filename = "temp.png", width = 5, height = 5) {
  ggplotly(g, width = width * 100, height = height * 100) %>%
    layout(autosize = F, 
           margin = list(l = 0, r = 0, b = 0, t = 0, pad = 0, autoexpand = T),
           scene = list(aspectmode = "data"),
           xaxis = list(autorange = F, range = list(xmin*1.0, xmax*1.0)),
           yaxis = list(autorange = F, range = list(ymin*1.0, ymax*1.0), scaleanchor= 'x', scaleratio = (xmax - xmin)/(ymax - ymin)),
           images = list(source = base64enc::dataURI(file = filename), 
                         x = 0, y = 0,
                         sizex = 1, sizey = 0.995,
                         xref = "xaxis", yref= "yaxis",
                         sizing = "stretch",
                         xanchor = "left", yanchor = "bottom", layer = "below"),
           showlegend = FALSE) %>% 
    config(displayModeBar = F)
}