#' Balanced downsampling of matrix/data.frame based on cluster assignment vector
#' @param df expression matrix or data.frame
#' @param vec vector of ids
#' @param frac fraction 0-1 to downsample to
#' @param seed sampling randomization seed
#' @return list with new downsampled matrix/data.frame and id vector
#' @examples
#' res <- by_cluster_sampling(data.frame(y = c(1, 2, 3, 4, 5, 6)), vec = c(1, 2, 1, 2, 1, 2), frac = 0.5)
#' @export
by_cluster_sampling <- function(df, vec, frac, seed = 34) {
  dfs <- split(df, vec)
  vecout <- c()
  dflist <- list()
  for (x in names(dfs)) {
    df1 <- dfs[[x]]
    set.seed(seed)
    seed <- seed + 1
    samp <- sample(1:nrow(df1), round((frac * nrow(df1))))
    em1 <- df1[samp, , drop = F]
    vec1 <- rep(x, round((frac * nrow(df1))))
    vecout <- c(vecout, vec1)
    dflist[[x]] <- em1
  }
  dfout <- do.call(rbind, dflist)
  list(dfout, vecout)
}

#' Rowwise math from matrix/data.frame per cluster based on another vector/metadata, similar to clustifyr::average_clusters but ids as rows
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
#' @param subclusterpower whether to get multiple averages per original cluster
#' @param cut_n set on a limit of genes as expressed, lower ranked genes
#' are set to 0, considered unexpressed
#' @param trim whether to remove 1 percentile when doing min caluculation
#' @return average expression matrix, with genes for row names, and clusters
#'  for column names
#' @examples
#' mat <- average_clusters_rowwise(data.frame(y = c(1, 2, 3, 4, 5, 6), x = c(1, 2, 3, 4, 5, 6)), metadata = c(1, 2, 1, 2, 1, 2), method = "min")
#' @importFrom matrixStats rowMaxs rowMedians colRanks
#' @export
average_clusters_rowwise <- function(mat, metadata, cluster_col = "cluster", if_log = FALSE,
                                     cell_col = NULL, low_threshold = 0, method = "mean", output_log = FALSE,
                                     subclusterpower = 0, cut_n = NULL, trim = FALSE) {
  cluster_info <- metadata
  if (!(is.null(cell_col))) {
    if (!(all(rownames(mat) == cluster_info[[cell_col]]))) {
      mat <- mat[cluster_info[[cell_col]], ]
    }
  }
  if (is.null(rownames(mat))) {
    stop("The input matrix does not have rownames.\n", "Check colnames() of input object")
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
  if (subclusterpower > 0) {
    cluster_ids <- overcluster(mat, cluster_ids, power = subclusterpower)
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
#' @export
get_labs <- function(g) {
  g2 <- ggplot2::ggplot_build(g)
  g2$plot$scales$scales[[1]]$get_labels()
}

check_colour_fill <- function(g) {
  col <- "fill"
  g2 <- ggplot2::ggplot_build(g)
  cols <- arrange(g2$data[[1]], group)
  cols <- unique(pull(cols, col))
  if (length(cols) <= 1) {
    if (col == "fill") {
      col <- "colour"
    } else {
      col <- "fill"
    }
  }
  col
}