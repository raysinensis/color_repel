#' Score matrix distances
#' @param dist1 distanct matrix 1
#' @param dist2 distanct matrix 2
#' @return numeric score
matrix2_score <- function(dist1, dist2) {
  temp <- dist1 * dist2
  max(temp, na.rm = TRUE)
}

#' Score matrix distances in multiple combinations
#' @param dist1 distanct matrix 1
#' @param dist2 distanct matrix 2
#' @param n number of iterations
#' @param verbose whether to output more messages
#' @param seed random seed
#' @param out_worst instead of default output of best combination, output worst instead
#' @return reordered vector
matrix2_score_n <- function(dist1,
                            dist2,
                            n = min(factorial(ncol(dist2)) * 10, 20000),
                            verbose = FALSE,
                            seed = 34,
                            out_worst = FALSE) {
  dist1[dist1 == Inf] <- NA
  dist2[dist2 == Inf] <- NA
  len <- ncol(dist2)

  s <- vector("list", n)
  dqrng::dqset.seed(seed)
  for (i in 1:n) {
    s[[i]] <- dqrng::dqsample.int(len, len)
  }
  s <- unique(s)
  if (verbose) {
    message("attempting ", length(s), " calculations...")
    if (length(s) == min(factorial(len))) {
      message("all color combos covered")
    }
  }

  ord1 <- 1:ncol(dist2)
  score1 <- matrix2_score(dist1, dist2)
  score0 <- score1
  scoremax <- score1

  for (i in 1:length(s)) {
    ord_temp <- s[[i]]
    dist3 <- dist2[ord_temp, ord_temp]
    score_temp <- matrix2_score(dist1, dist3)
    if (score_temp > scoremax) {
      scoremax <- score_temp
      if (out_worst) {
        ord1 <- ord_temp
      }
    }
    if (score_temp < score1) {
      if (!out_worst) {
        ord1 <- ord_temp
      }
      score1 <- score_temp
    }
  }
  if (verbose) {
    scale1 <- 10^floor(log10(score0))
    message("original score (scaled): ", score0 / scale1)
    message("worst score: ", scoremax / scale1)
    message("optimal score: ", score1 / scale1)
  }
  ord1
}

#' Balanced downsampling of matrix/data.frame based on cluster assignment vector
#' @param df expression matrix or data.frame
#' @param vec vector of ids
#' @param frac fraction 0-1 to downsample to
#' @param seed sampling randomization seed
#' @return list with new downsampled matrix/data.frame and id vector
#' @examples
#' res <- by_cluster_sampling(data.frame(y = c(1, 2, 3, 4, 5, 6)),
#'   vec = c(1, 2, 1, 2, 1, 2), frac = 0.5
#' )
#' @export
by_cluster_sampling <- function(df, vec, frac, seed = 34) {
  dfs <- split(df, vec)
  vecout <- c()
  dflist <- list()
  set.seed(seed)
  for (x in names(dfs)) {
    df1 <- dfs[[x]]
    samp <- sample(1:nrow(df1), round((frac * nrow(df1))))
    em1 <- df1[samp, , drop = FALSE]
    vec1 <- rep(x, round((frac * nrow(df1))))
    vecout <- c(vecout, vec1)
    dflist[[x]] <- em1
  }
  dfout <- do.call(rbind, dflist)
  list(dfout, vecout)
}

#' Rowwise math from matrix/data.frame per cluster based on another vector/metadata,
#' similar to clustifyr::average_clusters but ids as rows
#' @param mat expression matrix
#' @param metadata data.frame or vector containing cluster assignments per cell.
#' Order must match column order in supplied matrix. If a data.frame
#' provide the cluster_col parameters.
#' @param if_log input data is natural log,
#' averaging will be done on unlogged data
#' @param cluster_col column in metadata with cluster number
#' @param cell_col if provided, will reorder matrix first
#' @param low_threshold option to remove clusters with too few cells
#' @param method whether to take mean (default), median, 10% truncated mean, or trimean,
#' max, min
#' @param output_log whether to report log results
#' @param cut_n set on a limit of genes as expressed, lower ranked genes
#' are set to 0, considered unexpressed
#' @param trim whether to remove 1 percentile when doing min caluculation
#' @return average expression matrix, with genes for row names, and clusters
#'  for column names
#' @examples
#' mat <- average_clusters_rowwise(data.frame(
#'   y = c(1, 2, 3, 4, 5, 6),
#'   x = c(1, 2, 3, 4, 5, 6)
#' ), metadata = c(1, 2, 1, 2, 1, 2), method = "min")
#' @export
average_clusters_rowwise <- function(mat, metadata, cluster_col = "cluster", if_log = FALSE,
                                     cell_col = NULL, low_threshold = 0, method = "mean", output_log = FALSE,
                                     cut_n = NULL, trim = FALSE) {
  cluster_info <- metadata
  if (!(is.null(cell_col))) {
    if (!(all(rownames(mat) == cluster_info[[cell_col]]))) {
      mat <- mat[cluster_info[[cell_col]], ]
    }
  }

  if (is.vector(cluster_info)) {
    if (nrow(mat) != length(cluster_info)) {
      stop("vector of cluster assignments does not match the number of rows in the matrix",
        call. = FALSE
      )
    }
    cluster_ids <- split(rownames(mat), cluster_info)
  } else if (is.data.frame(cluster_info) & !is.null(cluster_col)) {
    if (!is.null(cluster_col) && !(cluster_col %in% colnames(metadata))) {
      stop("given `cluster_col` is not a column in `metadata`",
        call. = FALSE
      )
    }
    cluster_info_temp <- cluster_info[[cluster_col]]
    if (is.factor(cluster_info_temp)) {
      cluster_info_temp <- droplevels(cluster_info_temp)
    }
    cluster_ids <- split(rownames(mat), cluster_info_temp)
  } else if (is.factor(cluster_info)) {
    cluster_info <- as.character(cluster_info)
    if (nrow(mat) != length(cluster_info)) {
      stop("vector of cluster assignments does not match the number of rows in the matrix",
        call. = FALSE
      )
    }
    cluster_ids <- split(rownames(mat), cluster_info)
  } else {
    stop("metadata not formatted correctly,\n         supply either a vector or a dataframe",
      call. = FALSE
    )
  }
  if (method == "mean") {
    out <- lapply(cluster_ids, function(cell_ids) {
      if (!all(cell_ids %in% colnames(mat))) {
        stop("cell ids not found in input matrix", call. = FALSE)
      }
      if (if_log) {
        mat_data <- expm1(mat[cell_ids, , drop = FALSE])
      } else {
        mat_data <- mat[cell_ids, , drop = FALSE]
      }
      res <- Matrix::colMeans(mat_data, na.rm = TRUE)
      if (output_log) {
        res <- log1p(res)
      }
      res
    })
  } else if (method == "median") {
    out <- lapply(cluster_ids, function(cell_ids) {
      if (!all(cell_ids %in% colnames(mat))) {
        stop("cell ids not found in input matrix", call. = FALSE)
      }
      mat_data <- mat[cell_ids, , drop = FALSE]
      res <- matrixStats::colMedians(as.matrix(mat_data),
        na.rm = TRUE
      )
      res[is.na(res)] <- 0
      names(res) <- colnames(mat_data)
      res
    })
  } else if (method == "trimean") {
    out <- lapply(cluster_ids, function(cell_ids) {
      if (!all(cell_ids %in% colnames(mat))) {
        stop("cell ids not found in input matrix", call. = FALSE)
      }
      mat_data <- mat[cell_ids, , drop = FALSE]
      res1 <- matrixStats::colQuantiles(as.matrix(mat_data),
        probs = 0.25, na.rm = TRUE
      )
      res2 <- matrixStats::colQuantiles(as.matrix(mat_data),
        probs = 0.5, na.rm = TRUE
      )
      res3 <- matrixStats::colQuantiles(as.matrix(mat_data),
        probs = 0.75, na.rm = TRUE
      )
      res <- 0.5 * res2 + 0.25 * res1 + 0.25 * res3
      res[is.na(res)] <- 0
      names(res) <- colnames(mat_data)
      res
    })
  } else if (method == "truncate") {
    out <- lapply(cluster_ids, function(cell_ids) {
      if (!all(cell_ids %in% rownames(mat))) {
        stop("cell ids not found in input matrix", call. = FALSE)
      }
      mat_data <- mat[cell_ids, , drop = FALSE]
      res <- apply(mat_data, 2, function(x) {
        mean(x,
          trim = 0.1,
          na.rm = TRUE
        )
      })
      rownames(res) <- names(cell_ids)
      res
    })
  } else if (method == "min") {
    out <- lapply(cluster_ids, function(cell_ids) {
      if (!all(cell_ids %in% rownames(mat))) {
        stop("cell ids not found in input matrix", call. = FALSE)
      }
      mat_data <- mat[cell_ids, , drop = FALSE]
      if (trim) {
        res <- matrixStats::colQuantiles(as.matrix(mat_data),
          na.rm = TRUE, probs = 0.01
        )
      } else {
        res <- matrixStats::colMins(as.matrix(mat_data),
          na.rm = TRUE
        )
      }
      res[is.na(res)] <- 0
      names(res) <- colnames(mat_data)
      res
    })
  } else if (method == "max") {
    out <- lapply(cluster_ids, function(cell_ids) {
      if (!all(cell_ids %in% rownames(mat))) {
        stop("cell ids not found in input matrix", call. = FALSE)
      }
      mat_data <- mat[cell_ids, , drop = FALSE]
      res <- matrixStats::colMaxs(as.matrix(mat_data),
        na.rm = TRUE
      )
      res[is.na(res)] <- 0
      names(res) <- colnames(mat_data)
      res
    })
  }
  out <- do.call(cbind, out)
  if (low_threshold > 0) {
    fil <- vapply(cluster_ids, FUN = length, FUN.VALUE = numeric(1)) >=
      low_threshold
    if (!all(as.vector(fil))) {
      message(
        "The following clusters have less than ",
        low_threshold, " cells for this analysis: ",
        paste(colnames(out)[!as.vector(fil)], collapse = ", "),
        ". They are excluded."
      )
    }
    out <- out[as.vector(fil), ]
  } else {
    fil <- vapply(cluster_ids, FUN = length, FUN.VALUE = numeric(1)) >=
      10
    if (!all(as.vector(fil))) {
      message(
        "The following clusters have less than ",
        10, " cells for this analysis: ", paste(rownames(out)[!as.vector(fil)],
          collapse = ", "
        ), ". Classification is likely inaccurate."
      )
    }
  }
  if (!(is.null(cut_n))) {
    expr_mat <- out
    expr_df <- as.matrix(expr_mat)
    df_temp <- t(matrixStats::rowRanks(-expr_df, ties.method = "average"))
    rownames(df_temp) <- rownames(expr_mat)
    colnames(df_temp) <- colnames(expr_mat)
    expr_mat[df_temp > cut_n] <- 0
    out <- expr_mat
  }
  return(out)
}

#' Extract custom labels from ggplot object
#' @param g ggplot object
#' @return named vector of labels
#' @examples
#' a <- ggplot2::ggplot(ggplot2::mpg, ggplot2::aes(displ, hwy)) +
#'   ggplot2::geom_point(ggplot2::aes(color = as.factor(cyl))) +
#'   ggplot2::geom_text(ggplot2::aes(label = model))
#' get_labs(a)
#' @export
get_labs <- function(g) {
  g2 <- ggplot2::ggplot_build(g)
  nlayer <- length(g2$plot$scales$scales)
  for (x in 1:nlayer) {
    ls <- g2$plot$scales$scales[[x]]$get_labels()
    if (length(ls) > 0) {
      return(ls)
    }
  }
}

check_colour_mapping <- function(g, col = "colour", return_col = FALSE, autoswitch = TRUE) {
  g2 <- ggplot2::ggplot_build(g)
  cols <- dplyr::arrange(g2$data[[1]], group)
  cols <- unique(cols[[col]])
  if (length(cols) <= 1) {
    if (!autoswitch) {
      stop("only 1 colour detected, please check mapping")
    }
    if (col == "fill") {
      col <- "colour"
    } else {
      col <- "fill"
    }
    cols <- dplyr::arrange(g2$data[[1]], group)
    cols <- unique(cols[[col]])
  }
  if (return_col) {
    list(col = col, cols = cols)
  } else {
    col
  }
}

#' Distance calculations for spatial coord
#' @param coord dataframe or matrix of spatial coordinates, cell barcode as rownames
#' @param metadata data.frame or vector containing cluster assignments per cell.
#' Order must match column order in supplied matrix. If a data.frame
#' provide the cluster_col parameters.
#' @param cluster_col column in metadata with cluster number
#' @param collapse_to_cluster instead of reporting min distance to cluster per cell, summarize to cluster level
#' @return min distance matrix
calc_distance <- function(
    coord,
    metadata,
    cluster_col = "cluster",
    collapse_to_cluster = FALSE) {
  distm <- distances::distances(coord)
  res <- average_clusters(distm,
    metadata,
    cluster_col,
    if_log = FALSE,
    output_log = FALSE,
    method = "min"
  )
  if (collapse_to_cluster) {
    res2 <- average_clusters(t(res),
      metadata,
      cluster_col,
      if_log = FALSE,
      output_log = FALSE,
      method = "min"
    )
    res2
  } else {
    res
  }
}

#' Average expression values per cluster
#'
#' @param mat expression matrix
#' @param metadata data.frame or vector containing cluster assignments per cell.
#' Order must match column order in supplied matrix. If a data.frame
#' provide the cluster_col parameters.
#' @param if_log input data is natural log,
#' averaging will be done on unlogged data
#' @param cluster_col column in metadata with cluster number
#' @param cell_col if provided, will reorder matrix first
#' @param low_threshold option to remove clusters with too few cells
#' @param method whether to take mean (default), median, 10% truncated mean, or trimean, max, min
#' @param output_log whether to report log results
#' @param cut_n set on a limit of genes as expressed, lower ranked genes
#' are set to 0, considered unexpressed
#' @return average or other desired calculation by group/cluster matrix
average_clusters <- function(mat,
                             metadata,
                             cluster_col = "cluster",
                             if_log = TRUE,
                             cell_col = NULL,
                             low_threshold = 0,
                             method = "mean",
                             output_log = TRUE,
                             cut_n = NULL) {
  cluster_info <- metadata
  if (!(is.null(cell_col))) {
    if (!(all(colnames(mat) == cluster_info[[cell_col]]))) {
      mat <- mat[, cluster_info[[cell_col]]]
    }
  }

  if (is.vector(cluster_info)) {
    if (ncol(mat) != length(cluster_info)) {
      stop("vector of cluster assignments does not match the number of columns in the matrix",
        call. = FALSE
      )
    }
    if (!is.null(colnames(mat))) {
      cluster_ids <- split(colnames(mat), cluster_info)
    } else {
      cluster_ids <- split(1:length(cluster_info), cluster_info)
    }
  } else if (is.data.frame(cluster_info) & !is.null(cluster_col)) {
    if (!is.null(cluster_col) &&
      !(cluster_col %in% colnames(metadata))) {
      stop("given `cluster_col` is not a column in `metadata`", call. = FALSE)
    }

    cluster_info_temp <- cluster_info[[cluster_col]]
    if (is.factor(cluster_info_temp)) {
      cluster_info_temp <- droplevels(cluster_info_temp)
    }
    cluster_ids <- split(colnames(mat), cluster_info_temp)
  } else if (is.factor(cluster_info)) {
    cluster_info <- as.character(cluster_info)
    if (ncol(mat) != length(cluster_info)) {
      stop("vector of cluster assignments does not match the number of columns in the matrix",
        call. = FALSE
      )
    }
    cluster_ids <- split(colnames(mat), cluster_info)
  } else {
    stop("metadata not formatted correctly,
         supply either a vector or a dataframe",
      call. = FALSE
    )
  }

  if (method == "mean") {
    out <- lapply(
      cluster_ids,
      function(cell_ids) {
        if (!all(cell_ids %in% colnames(mat))) {
          stop("cell ids not found in input matrix",
            call. = FALSE
          )
        }
        if (if_log) {
          mat_data <- expm1(mat[, cell_ids, drop = FALSE])
        } else {
          mat_data <- mat[, cell_ids, drop = FALSE]
        }
        res <- Matrix::rowMeans(mat_data, na.rm = TRUE)
        if (output_log) {
          res <- log1p(res)
        }
        res
      }
    )
  } else if (method == "median") {
    out <- lapply(
      cluster_ids,
      function(cell_ids) {
        if (!all(cell_ids %in% colnames(mat))) {
          stop("cell ids not found in input matrix",
            call. = FALSE
          )
        }
        mat_data <- mat[, cell_ids, drop = FALSE]
        # mat_data[mat_data == 0] <- NA
        res <- matrixStats::rowMedians(as.matrix(mat_data),
          na.rm = TRUE
        )
        res[is.na(res)] <- 0
        names(res) <- rownames(mat_data)
        res
      }
    )
  } else if (method == "trimean") {
    out <- lapply(
      cluster_ids,
      function(cell_ids) {
        if (!all(cell_ids %in% colnames(mat))) {
          stop("cell ids not found in input matrix",
            call. = FALSE
          )
        }
        mat_data <- mat[, cell_ids, drop = FALSE]
        # mat_data[mat_data == 0] <- NA
        res1 <- matrixStats::rowQuantiles(as.matrix(mat_data),
          probs = 0.25,
          na.rm = TRUE
        )
        res2 <- matrixStats::rowQuantiles(as.matrix(mat_data),
          probs = 0.5,
          na.rm = TRUE
        )
        res3 <- matrixStats::rowQuantiles(as.matrix(mat_data),
          probs = 0.75,
          na.rm = TRUE
        )
        res <- 0.5 * res2 + 0.25 * res1 + 0.25 * res3
        res[is.na(res)] <- 0
        names(res) <- rownames(mat_data)
        res
      }
    )
  } else if (method == "truncate") {
    out <- lapply(
      cluster_ids,
      function(cell_ids) {
        if (!all(cell_ids %in% colnames(mat))) {
          stop("cell ids not found in input matrix",
            call. = FALSE
          )
        }
        mat_data <- mat[, cell_ids, drop = FALSE]
        # mat_data[mat_data == 0] <- NA
        res <- apply(mat_data, 1, function(x) mean(x, trim = 0.1, na.rm = TRUE))
        colnames(res) <- names(cell_ids)
        res
      }
    )
  } else if (method == "min") {
    out <- purrr::map(
      cluster_ids,
      function(cell_ids) {
        mat_data <- mat[, cell_ids, drop = FALSE]
        res <- matrixStats::rowMins(mat_data,
          na.rm = TRUE
        )
        res[is.na(res)] <- 0
        names(res) <- rownames(mat_data)
        res
      }
    )
  } else if (method == "max") {
    out <- lapply(
      cluster_ids,
      function(cell_ids) {
        if (!all(cell_ids %in% colnames(mat))) {
          stop("cell ids not found in input matrix",
            call. = FALSE
          )
        }
        mat_data <- mat[, cell_ids, drop = FALSE]
        # mat_data[mat_data == 0] <- NA
        res <- matrixStats::rowMaxs(as.matrix(mat_data),
          na.rm = TRUE
        )
        res[is.na(res)] <- 0
        names(res) <- rownames(mat_data)
        res
      }
    )
  }

  out <- do.call(cbind, out)
  if (low_threshold > 0) {
    fil <- vapply(cluster_ids,
      FUN = length,
      FUN.VALUE = numeric(1)
    ) >= low_threshold
    if (!all(as.vector(fil))) {
      message(
        "The following clusters have less than ", low_threshold, " cells for this analysis: ",
        paste(colnames(out)[!as.vector(fil)], collapse = ", "),
        ". They are excluded."
      )
    }
    out <- out[, as.vector(fil)]
  } else {
    fil <- vapply(cluster_ids,
      FUN = length,
      FUN.VALUE = numeric(1)
    ) >= 10
    if (!all(as.vector(fil))) {
      message(
        "The following clusters have less than ", 10, " cells for this analysis: ",
        paste(colnames(out)[!as.vector(fil)], collapse = ", "),
        ". Classification is likely inaccurate."
      )
    }
  }
  if (!(is.null(cut_n))) {
    expr_mat <- out
    expr_df <- as.matrix(expr_mat)
    df_temp <- t(matrixStats::colRanks(-expr_df,
      ties.method = "average"
    ))
    rownames(df_temp) <- rownames(expr_mat)
    colnames(df_temp) <- colnames(expr_mat)
    expr_mat[df_temp > cut_n] <- 0
    out <- expr_mat
  }

  return(out)
}

#' ggrepel labeling of clusters
#' @param g ggplot object or data.frame
#' @param group_col column name in data.frame, default to "label" or "group" in ggplot data
#' @param x column name in data.frame for x
#' @param y column name in data.frame for y
#' @param txt_pt text size
#' @param remove_current whether to remove current text
#' @param layer text layer to remove, defaults to last
#' @param ... arguments passed to geom_text_repel
#' @return function, if data.frame input, or new ggplot object
#' @examples
#' g <- label_repel(ggplot2::ggplot(mtcars, ggplot2::aes(x = hp, y = wt, color = as.character(cyl))) +
#'   ggplot2::geom_point(), remove_current = FALSE)
#' @export
label_repel <- function(g, group_col = "auto", x = "x", y = "y",
                        txt_pt = 3, remove_current = "auto", layer = "auto", ...) {
  g_orig <- g
  if (is.data.frame(g)) {
    so_df <- g
  } else {
    g2 <- ggplot2::ggplot_build(g)
    if (layer == "auto") {
      layer <- length(g2$data)
    }
    so_df <- g2$data[[layer]]
  }
  if (group_col == "auto") {
    if ("label" %in% colnames(so_df)) {
      group_col <- "label"
    } else {
      group_col <- "group"
    }
  }
  if (is.numeric(so_df[[group_col]])) {
    temp_group <- get_labs(g)
    so_df[[group_col]] <- factor(so_df[[group_col]], labels = temp_group)
  }
  centers <- dplyr::group_by(so_df, !!dplyr::sym(group_col))
  centers <- dplyr::summarize(centers,
    t1 = stats::median(!!dplyr::sym(x)),
    t2 = stats::median(!!dplyr::sym(y)),
    a = 1
  )
  centers <- dplyr::ungroup(centers)

  labdata <- dplyr::select(
    so_df,
    !!dplyr::sym(group_col),
    !!dplyr::sym(x),
    !!dplyr::sym(y)
  )
  labdata[[1]] <- ""
  labdata$a <- 0
  colnames(labdata) <- colnames(centers)
  alldata <- rbind(labdata, centers)

  d <- ggrepel::geom_text_repel(
    data = alldata,
    color = "black",
    size = txt_pt,
    mapping = ggplot2::aes(
      x = !!dplyr::sym("t1"),
      y = !!dplyr::sym("t2"),
      # fill = NA,
      # alpha = !!dplyr::sym("a"),
      label = .data[[group_col]]
    ),
    point.padding = 0.5,
    box.padding = 0.5,
    max.iter = 50000,
    max.overlaps = 10000,
    ...
  )

  if (is.data.frame(g)) {
    d
  } else {
    if (remove_current == "auto") {
      remove_current <- check_labels(g_orig)
    }
    if (remove_current) {
      g <- remove_current_labels(g, layer = layer)
    }
    g + d
  }
}

check_patchwork <- function(g, layer = 1) {
  if ("patchwork" %in% class(g)) {
    g[[layer]]
  } else {
    g
  }
}

check_labels <- function(g, layer = "auto", text = "text|label") {
  g <- check_patchwork(g)
  if (layer == "auto") {
    layer <- length(g[["layers"]])
  }
  cs <- stringr::str_to_lower(class(g[["layers"]][[layer]][["geom"]]))
  any(stringr::str_detect(cs, text))
}

remove_current_labels <- function(g, layer = "auto") {
  g <- check_patchwork(g)
  if (layer == "auto") {
    layer <- length(g[["layers"]])
  }
  g[["layers"]][[layer]] <- NULL
  g
}

prep_encircle <- function(g, threshold = 0.01, nmin = 0.01, downsample = 5000, seed = 42) {
  g <- check_patchwork(g)
  g <- ggplot2::ggplot_build(g)

  em <- dplyr::select(g$data[[1]], c(x, y))
  clust <- g$data[[1]]$group
  if (nrow(em) > downsample) {
    frac <- downsample / nrow(em)
    res <- by_cluster_sampling(em, clust, frac, seed = seed)
    em <- res[[1]]
    clust <- res[[2]]
  }
  ems <- split(em, clust)
  dat <- purrr::map(1:length(ems), function(x) {
    em1 <- ems[[x]]
    distm1 <- distances::distances(em1)
    distm1 <- as.matrix(distm1)
    cut1 <- stats::quantile(unlist(distm1), probs = threshold)
    n1 <- colSums(distm1 <= cut1)
    sel1 <- n1 >= (ceiling(nrow(em1) * nmin))
    dat1 <- em1[sel1, ]
    if (nrow(dat1) <= 3) {
      message("too few points remain in group ", names(ems)[x])
    }
    dat1$group <- names(ems)[x]
    dat1
  })
  dplyr::bind_rows(dat)
}

expand_lims <- function(xmin, xmax, by = 0.1) {
  len <- xmax - xmin
  return(c(xmin - len*by, xmax + len*by))
}
