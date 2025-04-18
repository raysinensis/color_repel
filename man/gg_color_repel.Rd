% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wrappers.R
\name{gg_color_repel}
\alias{gg_color_repel}
\title{Wrapper to reorder ggplot colors to maximize color differences in space}
\usage{
gg_color_repel(
  g = ggplot2::last_plot(),
  col = "colour",
  sim = NULL,
  severity = 0.5,
  verbose = FALSE,
  downsample = 5000,
  nsamp = 50000,
  polychrome_recolor = FALSE,
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
  mascarade = FALSE,
  ggbuild = NULL,
  ...
)
}
\arguments{
\item{g}{ggplot plot object}

\item{col}{colour or fill in ggplot}

\item{sim}{passing a colorbind simulation function if needed}

\item{severity}{severity of the color vision defect, between 0 and 1}

\item{verbose}{whether to print messages}

\item{downsample}{downsample when too many datapoints are present}

\item{nsamp}{how many random sampling color combinations to test, default 50000}

\item{polychrome_recolor}{whether to replace the original colors with polychrome creation}

\item{seed}{sampling randomization seed}

\item{autoswitch}{try to switch between colour and fill automatically}

\item{layer}{layer to detect color, defaults to first}

\item{out_orig}{output the original colors as named vector}

\item{out_worst}{output the worst combination instead of best}

\item{repel_label}{whether to add centroid labels with ggrepel}

\item{encircle}{whether to draw geom_encircle by cluster}

\item{encircle_alpha}{alpha argument passed to geom_encircle}

\item{encircle_expand}{expand argument passed to geom_encircle}

\item{encircle_shape}{shape/smoothing argument passed to geom_encircle}

\item{encircle_threshold}{threshold for removing outliers}

\item{encircle_nmin}{number of near neighbors for removing outliers}

\item{mascarade}{use mascarade package to outline clusters}

\item{ggbuild}{already built ggplot_built object if available}

\item{...}{passed to repel_label}
}
\value{
new ggplot object
}
\description{
Wrapper to reorder ggplot colors to maximize color differences in space
}
\examples{
a <- ggplot2::ggplot(ggplot2::mpg, ggplot2::aes(displ, hwy)) +
  ggplot2::geom_point(ggplot2::aes(color = as.factor(cyl)))
b <- gg_color_repel(a, col = "colour")
}
