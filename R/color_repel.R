#' Reorder ggplot colors to maximize color differences in space
#' @param g ggplot plot object
#' @param coord coordinates, default is inferred
#' @param groups groups corresponding to color/fill, default is inferred
#' @param nsamp how many random sampling color combinations to test, default 50000
#' @param sim passing a colorbind simulation function if needed
#' @param severity severity of the color vision defect, between 0 and 1
#' @param verbose whether to print messages
#' @param downsample downsample when too many datapoints are present
#' @param seed sampling randomization seed
#' @param col colour or fill in ggplot
#' @param autoswitch try to switch between colour and fill automatically
#' @param out_orig output the original colors as named vector
#' @param out_worst output the worst combination instead of best
#' @examples 
#' a <- ggplot2::ggplot(ggplot2::mpg, ggplot2::aes(displ, hwy)) +
#'  ggplot2::geom_point(ggplot2::aes(color = as.factor(cyl)))
#' new_colors <- color_repel(a)
#' b <- a + ggplot2::scale_color_manual(values = new_colors)
#' @return vector of reordered colors
#' @export
color_repel <- function(g,
                        coord = NULL,
                        groups = NULL,
                        nsamp = 50000,
                        sim = NULL,
                        severity = 0.5,
                        verbose = FALSE,
                        downsample = 5000,
                        seed = 34,
                        col = "colour",
                        autoswitch = TRUE,
                        out_orig = FALSE,
                        out_worst = FALSE) {
  g <- check_patchwork(g)
  if (verbose) {
    message("extract original colors...")
  }
  temp <- check_colour_mapping(g, col = col, return_col = TRUE, autoswitch = autoswitch)
  col <- temp[["col"]]
  cols <- temp[["cols"]]
  g2 <- ggplot2::ggplot_build(g)

  if (length(cols) <= 1) {
    warning("Did not detect multiple colors, did you specify the correct mapping? Trying to autoswitch...")
  }
  if (verbose) {
    message(cols)
  }
  orig_cols <- cols
  if (out_orig) {
    temp <- orig_cols
    names(temp) <- sort(unique(g$data[[ggplot2::as_label(g$mapping[[col]])]]))
    return(temp)
  }

  # deficiency simulation
  if (!is.null(sim)) {
    cols <- do.call(sim, c(list(cols), severity = list(severity)))
  }
  # rgb matrix
  colsm <- t(grDevices::col2rgb(cols))
  # convert to lab
  colslab <- grDevices::convertColor(colsm, from = "sRGB", to = "Lab")
  # euclidean distance
  coldist <- as.matrix(stats::dist(colslab))
  coldist[coldist == 0] <- Inf
  coldist <- (coldist - min(coldist[coldist != 0])) / 1000

  if (verbose) {
    message("extract plot distances...")
  }
  if (all(c("x", "y") %in% colnames(g2$data[[1]]))) {
    em <- dplyr::select(g2$data[[1]], x, y)
    # clustering info
    clust <- as.character(g2$data[[1]][[col]])
    clust <- as.character(as.numeric(factor(clust, levels = orig_cols)))
    if (nrow(em) > downsample) {
      frac <- downsample / nrow(em)
      res <- by_cluster_sampling(em, clust, frac, seed = seed)
      em <- res[[1]]
      clust <- res[[2]]
    }
    # min distance between clusters on plot
    cdist <- suppressMessages(calc_distance(em, clust))
    if (verbose) {
      message("extract plot distances (part 2)...")
    }
    rownames(cdist) <- as.character(1:nrow(cdist))
    cdist <- suppressMessages(average_clusters_rowwise(cdist, metadata = clust, if_log = FALSE, method = "min", output_log = F, trim = T))
    ord <- gtools::mixedorder(colnames(cdist))
    cdist <- cdist[ord, ord]
    cdist[cdist < max(cdist) / 100] <- max(cdist) / 100
    cdist[cdist > max(cdist) / 3] <- NA
  } else {
    cdist <- as.matrix(stats::dist(data.frame(x = unique(g2$data[[1]]$group))))
    cdist <- cdist^2
  }
  if (verbose) {
    message(dim(cdist))
  }
  if (is.null(nsamp)) {
    nsamp <- min(factorial(ncol(cdist)) * 5, 20000)
  }
  if (verbose) {
    message("iterate color combinations...")
  }
  res <- matrix2_score_n(1 / cdist, 1 / coldist, n = nsamp, verbose = verbose, seed = seed, out_worst = out_worst)
  temp <- orig_cols[res]
  names(temp) <- sort(unique(g$data[[ggplot2::as_label(g$mapping[[col]])]]))
  temp
}