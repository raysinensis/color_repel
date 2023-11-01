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
                        autoswitch = TRUE) {
  if (verbose) {
    message("extract original colors...")
  }
  g2 <- ggplot2::ggplot_build(g)
  cols <- arrange(g2$data[[1]], group)
  cols <- unique(pull(cols, col))
  if (verbose) {
    message(cols)
  }
  if (autoswitch) {
    if (length(cols) <= 1) {
      if (col == "fill") {
        col <- "colour"
      } else {
        col <- "fill"
      }
      cols <- arrange(g2$data[[1]], group)
      cols <- unique(pull(cols, col))
    }
  }
  if (length(cols) <= 1) {
    warning("Did not detect multiple colors, did you specify the correct mapping? Trying to autoswitch...")
  }
  orig_cols <- cols

  # deficiency simulation
  if (!is.null(sim)) {
    cols <- do.call(sim, c(list(cols), severity = list(severity)))
  }
  # cols <- colorspace::tritan(cols, severity = 0.5)
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
    em <- select(g2$data[[1]], x, y)
    # clustering info
    clust <- as.character(as.numeric(as.factor(as.character(g2$data[[1]][[col]]))))
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
    cdist <- suppressMessages(average_clusters_rowwise(cdist,metadata = clust, if_log = F, method = "min", output_log = F, trim = T))
    cdist[cdist < max(cdist) / 100] <- max(cdist) / 100
    cdist[cdist > max(cdist) / 3] <- NA
  } else {
    cdist <- as.matrix(stats::dist(data.frame(x = unique(g2$data[[1]]$group))))
    cdist <- cdist^2
  }

  if (verbose) {
    message("iterate color combinations...")
  }
  if (is.null(nsamp)) {
    nsamp <- min(factorial(ncol(cdist)) * 5, 20000)
  }
  res <- matrix2_score_n(1/cdist, 1/coldist, n = nsamp, verbose = verbose, seed = seed)
  orig_cols[res]
}

# matrix2_score <- function(dist1, dist2) {
#   temp <- dist1 * dist2
#   # temp[temp == Inf] <- NA
#   # temp <- colSums(temp, na.rm = T)
#   sum(temp, na.rm = T)/(ncol(temp)^2-sum(is.na(temp)))
# }

matrix2_score <- function(dist1, dist2) {
  temp <- dist1 * dist2
  max(temp, na.rm = T)
}

matrix2_score_n <- function(dist1, 
                            dist2, 
                            n = min(factorial(ncol(dist2)) * 10, 20000),
                            verbose = F,
                            seed = 34) {
  dist1[dist1 == Inf] <- NA
  dist2[dist2 == Inf] <- NA
  len <- ncol(dist2)
  ord1 <- 1:ncol(dist2)
  score1 <- matrix2_score(dist1, dist2)
  score0 <- score1
  scoremax <- score1
  s <- vector("list", n)
  # set.seed(seed)
  dqrng::dqset.seed(seed)
  for (i in 1:n) {
    s[[i]] <- dqrng::dqsample.int(len, len)
  }
  s <- unique(s)
  if (verbose) {
    message("attempting ", length(s), " calcuations...")
    if (length(s) == min(factorial(len))) {
      message("all color combos covered")
    }
  }
  for (i in 1:length(s)) {
    ord_temp <- s[[i]]
    dist3 <- dist2[ord_temp, ord_temp]
    score_temp <- matrix2_score(dist1, dist3)
    if (score_temp > scoremax) {
      scoremax <- score_temp
    }
    if (score_temp < score1) {
      ord1 <- ord_temp
      score1 <- score_temp
    }
  }
  if (verbose) {
    scale1 <- 10 ^ floor(log10(score0))
    # message("scale: ", scale1)
    message(ord1)
    message("original score: ", score0 / scale1)
    message("worst score: ", scoremax / scale1)
    message("optimal score: ", score1 / scale1)
  }
  ord1
}

create_matrix_lookup <- function(mat1, mat2) {
  res <- list()
  for (i in 1:ncol(mat1)) {
    for (j in 1:ncol(mat2)) {
      res[[paste0(i, "-", j)]] <- 1/(mat1[, i, drop = T] * mat2[, j, drop = T])
    }
  }
  res
}

matrix_lookup <- function(mat1, mat2, s) {
   ml <- create_matrix_lookup(mat1, mat2)
   scores <- c()
   for (i in 1:length(s)) {
     l <- str_c(1:length(s[[i]]), "-", s[[i]])
     temp <- data.frame(ml[l])
     temp[temp == Inf] <- NA
     scores[i] <- mean(rowSums(temp, na.rm = T), na.rm = T)
   }
   scores
}
