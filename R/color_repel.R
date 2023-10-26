color_repel <- function(g, coord = NULL, groups = NULL, n = 100000, sim = NULL) {
  g2 <- ggplot2::ggplot_build(g)
  cols <- g2$data[[1]] %>% arrange(group) %>% pull(colour) %>% unique()
  orig_cols <- cols
  
  # deficiency simulation
  if (!is.null(sim)) {
    cols <- do.call(sim, list(cols))
  }
  # cols <- colorspace::tritan(cols, severity = 0.5)
  # rgb matrix
  colsm <- col2rgb(cols) %>% t()
  # convert to lab
  colslab <- grDevices::convertColor(colsm, from = "sRGB", to = "Lab")
  # euclidean distance
  coldist <- dist(colslab) %>% as.matrix()
  coldist[coldist == 0] <- Inf
  coldist <- (coldist - min(coldist[coldist != 0]))/1000
  
  # em <- Embeddings(so, "umap")
  em <- g2$data[[1]] %>% select(x,y)
  # clustering info
  # clust <- so$type
  clust <- as.character(g2$data[[1]]$group)
  # min distance between clusters on plot
  cdist <- clustifyr::calc_distance(em, clust) %>% t() %>% average_clusters_trim(metadata = clust, if_log = F, method = "min", output_log = F, trim = T)
  cdist[cdist < max(cdist)/100] <- max(cdist)/100
  cdist[cdist > max(cdist)/3] <- NA
  
  res <- matrix2_score_n(cdist, coldist, n = n)
  orig_cols[res]
}

matrix2_score <- function(dist1, dist2) {
  1/(dist1 * dist2) %>% rowSums(na.rm = T) %>% mean(na.rm = T)
}

matrix2_score_n <- function(dist1, dist2, n = min(factorial(ncol(dist2))*10, 100000)) {
  ord1 <- 1:ncol(dist2)
  score1 <- matrix2_score(dist1, dist2)
  s <- list()
  for (i in 1:n) {
    s[[i]] <- sample(1:ncol(dist2))
  } 
  s <- s %>% unique()
  for (i in (length(s) + 1):n) {
    s[[i]] <- sample(1:ncol(dist2))
  }
  s <- s %>% unique()
  for (i in 1:length(s)) {
    ord_temp <- s[[i]]
    score_temp <- matrix2_score(dist1, dist2[,ord_temp])
    if (score_temp < score1) {
      ord1 <- ord_temp
      score1 <- score_temp
    }
  }
  ord_temp
}

average_clusters_trim <- function (mat, metadata, cluster_col = "cluster", if_log = TRUE, 
                               cell_col = NULL, low_threshold = 0, method = "mean", output_log = TRUE, 
                               subclusterpower = 0, cut_n = NULL, trim = FALSE) 
{
  cluster_info <- metadata
  if (!(is.null(cell_col))) {
    if (!(all(colnames(mat) == cluster_info[[cell_col]]))) {
      mat <- mat[, cluster_info[[cell_col]]]
    }
  }
  if (is.null(colnames(mat))) {
    stop("The input matrix does not have colnames.\n", "Check colnames() of input object")
  }
  if (is.vector(cluster_info)) {
    if (ncol(mat) != length(cluster_info)) {
      stop("vector of cluster assignments does not match the number of columns in the matrix", 
           call. = FALSE)
    }
    cluster_ids <- split(colnames(mat), cluster_info)
  }
  else if (is.data.frame(cluster_info) & !is.null(cluster_col)) {
    if (!is.null(cluster_col) && !(cluster_col %in% colnames(metadata))) {
      stop("given `cluster_col` is not a column in `metadata`", 
           call. = FALSE)
    }
    cluster_info_temp <- cluster_info[[cluster_col]]
    if (is.factor(cluster_info_temp)) {
      cluster_info_temp <- droplevels(cluster_info_temp)
    }
    cluster_ids <- split(colnames(mat), cluster_info_temp)
  }
  else if (is.factor(cluster_info)) {
    cluster_info <- as.character(cluster_info)
    if (ncol(mat) != length(cluster_info)) {
      stop("vector of cluster assignments does not match the number of columns in the matrix", 
           call. = FALSE)
    }
    cluster_ids <- split(colnames(mat), cluster_info)
  }
  else {
    stop("metadata not formatted correctly,\n         supply either a vector or a dataframe", 
         call. = FALSE)
  }
  if (subclusterpower > 0) {
    cluster_ids <- overcluster(mat, cluster_ids, power = subclusterpower)
  }
  if (method == "mean") {
    out <- lapply(cluster_ids, function(cell_ids) {
      if (!all(cell_ids %in% colnames(mat))) {
        stop("cell ids not found in input matrix", call. = FALSE)
      }
      if (if_log) {
        mat_data <- expm1(mat[, cell_ids, drop = FALSE])
      }
      else {
        mat_data <- mat[, cell_ids, drop = FALSE]
      }
      res <- Matrix::rowMeans(mat_data, na.rm = TRUE)
      if (output_log) {
        res <- log1p(res)
      }
      res
    })
  }
  else if (method == "median") {
    out <- lapply(cluster_ids, function(cell_ids) {
      if (!all(cell_ids %in% colnames(mat))) {
        stop("cell ids not found in input matrix", call. = FALSE)
      }
      mat_data <- mat[, cell_ids, drop = FALSE]
      res <- matrixStats::rowMedians(as.matrix(mat_data), 
                                     na.rm = TRUE)
      res[is.na(res)] <- 0
      names(res) <- rownames(mat_data)
      res
    })
  }
  else if (method == "trimean") {
    out <- lapply(cluster_ids, function(cell_ids) {
      if (!all(cell_ids %in% colnames(mat))) {
        stop("cell ids not found in input matrix", call. = FALSE)
      }
      mat_data <- mat[, cell_ids, drop = FALSE]
      res1 <- matrixStats::rowQuantiles(as.matrix(mat_data), 
                                        probs = 0.25, na.rm = TRUE)
      res2 <- matrixStats::rowQuantiles(as.matrix(mat_data), 
                                        probs = 0.5, na.rm = TRUE)
      res3 <- matrixStats::rowQuantiles(as.matrix(mat_data), 
                                        probs = 0.75, na.rm = TRUE)
      res <- 0.5 * res2 + 0.25 * res1 + 0.25 * res3
      res[is.na(res)] <- 0
      names(res) <- rownames(mat_data)
      res
    })
  }
  else if (method == "truncate") {
    out <- lapply(cluster_ids, function(cell_ids) {
      if (!all(cell_ids %in% colnames(mat))) {
        stop("cell ids not found in input matrix", call. = FALSE)
      }
      mat_data <- mat[, cell_ids, drop = FALSE]
      res <- apply(mat_data, 1, function(x) mean(x, trim = 0.1, 
                                                 na.rm = TRUE))
      colnames(res) <- names(cell_ids)
      res
    })
  }
  else if (method == "min") {
    out <- lapply(cluster_ids, function(cell_ids) {
      if (!all(cell_ids %in% colnames(mat))) {
        stop("cell ids not found in input matrix", call. = FALSE)
      }
      mat_data <- mat[, cell_ids, drop = FALSE]
      if (trim) {
        res <- matrixStats::rowQuantiles(as.matrix(mat_data), 
                                         na.rm = TRUE, probs = 0.01)
      } else {
        res <- matrixStats::rowMins(as.matrix(mat_data), 
                                    na.rm = TRUE)
      }
      res[is.na(res)] <- 0
      names(res) <- rownames(mat_data)
      res
    })
  }
  else if (method == "max") {
    out <- lapply(cluster_ids, function(cell_ids) {
      if (!all(cell_ids %in% colnames(mat))) {
        stop("cell ids not found in input matrix", call. = FALSE)
      }
      mat_data <- mat[, cell_ids, drop = FALSE]
      res <- matrixStats::rowMaxs(as.matrix(mat_data), 
                                  na.rm = TRUE)
      res[is.na(res)] <- 0
      names(res) <- rownames(mat_data)
      res
    })
  }
  out <- do.call(cbind, out)
  if (low_threshold > 0) {
    fil <- vapply(cluster_ids, FUN = length, FUN.VALUE = numeric(1)) >= 
      low_threshold
    if (!all(as.vector(fil))) {
      message("The following clusters have less than ", 
              low_threshold, " cells for this analysis: ", 
              paste(colnames(out)[!as.vector(fil)], collapse = ", "), 
              ". They are excluded.")
    }
    out <- out[, as.vector(fil)]
  }
  else {
    fil <- vapply(cluster_ids, FUN = length, FUN.VALUE = numeric(1)) >= 
      10
    if (!all(as.vector(fil))) {
      message("The following clusters have less than ", 
              10, " cells for this analysis: ", paste(colnames(out)[!as.vector(fil)], 
                                                      collapse = ", "), ". Classification is likely inaccurate.")
    }
  }
  if (!(is.null(cut_n))) {
    expr_mat <- out
    expr_df <- as.matrix(expr_mat)
    df_temp <- t(matrixStats::colRanks(-expr_df, ties.method = "average"))
    rownames(df_temp) <- rownames(expr_mat)
    colnames(df_temp) <- colnames(expr_mat)
    expr_mat[df_temp > cut_n] <- 0
    out <- expr_mat
  }
  return(out)
}
