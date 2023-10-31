#' Reorder ggplot colors to maximize color differences in space
#' @param g ggplot plot object
#' @param coord coordinates, default is inferred
#' @param groups groups corresponding to color/fill, default is inferred
#' @param nsamp how many random sampling color combinations to test, default 10000
#' @param sim passing a colorbind simulation function if needed
#' @param verbose whether to print messages
#' @param downsample downsample when too many datapoints are present
#' @param col colour or fill in ggplot
#' @return vector of reordered colors
#' @export
color_repel <- function(g, coord = NULL, groups = NULL, nsamp = NULL, sim = NULL, verbose = F, downsample = 10000, col = "colour") {
  if (verbose) {
    message("extract original colors...")
  }
  g2 <- ggplot2::ggplot_build(g)
  cols <- g2$data[[1]] %>%
    arrange(group) %>%
    pull(col) %>%
    unique()
  orig_cols <- cols

  # deficiency simulation
  if (!is.null(sim)) {
    cols <- do.call(sim, list(cols))
  }
  # cols <- colorspace::tritan(cols, severity = 0.5)
  # rgb matrix
  colsm <- grDevices::col2rgb(cols) %>% t()
  # convert to lab
  colslab <- grDevices::convertColor(colsm, from = "sRGB", to = "Lab")
  # euclidean distance
  coldist <- dist(colslab) %>% as.matrix()
  coldist[coldist == 0] <- Inf
  coldist <- (coldist - min(coldist[coldist != 0])) / 1000

  if (verbose) {
    message("extract plot distances...")
  }
  if (all(c("x", "y") %in% colnames(g2$data[[1]]))) {
    em <- g2$data[[1]] %>% select(x, y)
    # clustering info
    clust <- as.character(g2$data[[1]]$group)
    if (nrow(em) > downsample) {
      frac <- downsample / nrow(em)
      res <- by_cluster_sampling(em, clust, frac)
      em <- res[[1]]
      clust <- res[[2]]
    }
    # message(dim(em))
    # message(length(clust))
    # min distance between clusters on plot
    cdist <- clustifyr::calc_distance(em, clust)
    if (verbose) {
      message("extract plot distances (part 2)...")
    }
    cdist <- cdist %>% average_clusters_rowwise(metadata = clust, if_log = F, method = "min", output_log = F, trim = T)
    cdist[cdist < max(cdist) / 100] <- max(cdist) / 100
    cdist[cdist > max(cdist) / 3] <- NA
  } else {
    cdist <- data.frame(x = unique(g2$data[[1]]$group)) %>%
      dist() %>%
      as.matrix()
    cdist <- cdist^2
  }

  if (verbose) {
    message("iterate color combinations...")
  }
  if (is.null(nsamp)) {
    nsamp <- min(factorial(ncol(cdist)) * 10, 100000)
  }
  res <- matrix2_score_n(cdist, coldist, n = nsamp)
  orig_cols[res]
}

matrix2_score <- function(dist1, dist2) {
  1 / (dist1 * dist2) %>%
    rowSums(na.rm = T) %>%
    mean(na.rm = T)
}

matrix2_score_n <- function(dist1, dist2, n = min(factorial(ncol(dist2)) * 10, 100000)) {
  ord1 <- 1:ncol(dist2)
  score1 <- matrix2_score(dist1, dist2)
  s <- list()
  for (i in 1:n) {
    s[[i]] <- sample(1:ncol(dist2))
  }
  # message(length(s))
  s <- s %>% unique()
  for (i in (length(s) + 1):n) {
    s[[i]] <- sample(1:ncol(dist2))
  }
  s <- s %>% unique()
  # message(length(s))
  for (i in 1:length(s)) {
    ord_temp <- s[[i]]
    score_temp <- matrix2_score(dist1, dist2[, ord_temp])
    if (score_temp < score1) {
      ord1 <- ord_temp
      score1 <- score_temp
    }
  }
  ord_temp
}
