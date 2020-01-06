
#'
#' @title Scandal default UMAP configuration
#'
#' @description Returns a UMAP configuration object containing the default configuration which scandal uses to
#' compute the UMAP coordinates.
#'
#' @details The current configuration uses the Pearson's correlation as default distance metric.
#'
#' @return A \code{umap.defaults} object.
#'
#' @export
scandal_default_umap_config <- function() {

  config <- umap::umap.defaults
  config$metric <- "pearson"

  return (config)
}

#'
#' @title Setup analysis tools
#'
#' @description This function sets-up analysis-related tools such as computing a cell-correlation
#' matrix, t-SNE and UMAP to be used later in downstream analysis.
#'
#' @param object
#' @param computation_opts
#' @param cor_method
#' @param tsne_perplexity
#' @param tsne_pca
#' @param tsne_initial_dims
#' @param tsne_pca_center
#' @param tsne_pca_scale
#' @param tsne_normalize
#' @param umap_config
#'
#' @author Avishay Spitzer
#'
#' @export
scandal_setup_analysis <- function(object, computation_opts = c("correlation" = TRUE, "tsne" = TRUE, "umap" = TRUE),
                                   cor_method = "pearson",
                                   tsne_perplexity = 30, tsne_pca = TRUE, tsne_initial_dims = 50, tsne_pca_center = TRUE, tsne_pca_scale = FALSE, tsne_normalize = TRUE,
                                   umap_config = scandal_default_umap_config(),
                                   verbose = FALSE) {

  stopifnot(is_scandal_object(object))
  stopifnot(!is.null(computation_opts),
            is.vector(computation_opts),
            is.logical(computation_opts),
            base::all(names(computation_opts) %in% c("correlation", "tsne", "umap")) == TRUE)
  stopifnot(!is.null(cor_method), is.character(cor_method), length(cor_method) == 1, cor_method %in% c("pearson", "spearman", "kendall"))
  stopifnot(is.numeric(tsne_perplexity), is.logical(tsne_pca), is.numeric(tsne_initial_dims), is.logical(tsne_pca_center), is.logical(tsne_pca_scale), is.logical(tsne_normalize))
  stopifnot(!is.null(umap_config), class(umap_config) == class(scandal_default_umap_config()))

  if (.is_option_requested(computation_opts, "correlation")) {

    if (isTRUE(verbose))
      message(paste0("computing ", cor_method, " correlation matrix for ", nodeID(object)))

    reducedDim(object, "cor") <- .compute_cor(assay(object), cor_method = cor_method)
  }

  if (.is_option_requested(computation_opts, "tsne")) {

    if (isTRUE(verbose))
      message(paste0("computing t-SNE for ", nodeID(object)))

    reducedDim(object, "tsne") <- .compute_tsne(assay(object),
                                                perplexity = tsne_perplexity,
                                                pca = tsne_pca,
                                                initial_dims = tsne_initial_dims,
                                                check_duplicates = FALSE,
                                                pca_center = tsne_pca_center,
                                                pca_scale = tsne_pca_scale,
                                                normalize = tsne_normalize)
  }

  if (.is_option_requested(computation_opts, "umap")) {

    if (isTRUE(verbose))
      message(paste0("computing umap for ", nodeID(object)))

    reducedDim(object, "umap") <- .compute_umap(assay(object), umap_config = umap_config)
  }

  return (object)
}

#'
#' @importFrom FNN get.knn
#' @importFrom igraph graph_from_data_frame simplify cluster_louvain membership
#'
#' @export
compute_louvain_clusters <- function(data, k) {

  knn <- get.knn(as.matrix(data), k = k)

  knn <- data.frame(from = rep(1:nrow(knn$nn.index), k), to = as.vector(knn$nn.index), weight = 1/(1 + as.vector(knn$nn.dist)))

  nw <- graph_from_data_frame(knn, directed = FALSE)
  nw <- simplify(nw)

  lc <- cluster_louvain(nw)

  clusters <- membership(lc)
  clusters <- as.character(clusters)
  names(clusters) <- rownames(data)

  return (clusters)
}

#'
#' @title Compute cluster markers density
#'
#' @description This function computes the expression density of a set of marker genes
#' for each cluster. The expression density is defined as the mean expression of the
#' markers gene set accross all cells belonging to the same cluster.
#'
#' @param x gene expression matrix
#' @param clusters a character/factor vector of cluster assigment for each observation (column)
#' in \code{x}
#' @param markers a character vector of marker gene symbols corresponding to subset of genes (rows)
#' in \code{x}
#' @param return_sorted wether the result should be returned in descending order. Defaule is FALSE
#'
#' @return A \link{tibble} with two columns: Cluster and MCD (Markers Cluster Density).
#'
#' @author Avishay Spitzer
#'
#' @importFrom tibble tibble
#' @importFrom dplyr %>% group_by summarise arrange desc
#'
#' @export
scandal_cluster_markers_density <- function(x, clusters, mname, markers, return_sorted = FALSE) {

  stopifnot(is_valid_assay(x))

  markers <- markers[markers %in% rownames(x)]

  x <- center_matrix(as.matrix(x[markers, ]), by = "row", method = "mean", scale = FALSE)

  markers_mean <- apply(x, 2, mean)

  tbl <- tibble(Cluster = clusters, Mmean = markers_mean)

  tbl <- tbl %>%
            group_by(Cluster) %>%
            summarise (MCD = mean(Mmean))

  if (isTRUE(return_sorted))
    tbl <- tbl %>% arrange(desc(MCD))

  if (!is.null(mname))
    colnames(tbl)[colnames(tbl) == "MCD"] <- mname

  return (tbl)
}

#' @importFrom tibble is_tibble
#' @importFrom dplyr left_join
#'
#' @export
`%MCD%` <- function(t1, t2) {
  stopifnot(is_tibble(t1), is_tibble(t2))
  stopifnot("Cluster" %in% colnames(t1), "Cluster" %in% colnames(t2))
  stopifnot(nrow(t1) == nrow(t2))

  return (left_join(t1, t2, by = "Cluster"))
}

#' @importFrom tibble is_tibble
#' @importFrom dplyr left_join
#' @importFrom methods is
#'
#' @export
`%JCS%` <- function(o1, o2) {

  if (is_scandal_object(o1))
    o1 <- as_tibble(colData(o1), rownames = "CellID")
  else if (is_tibble(o1))
    stopifnot("CellID" %in% colnames(o1))

  if (is(o2, "ScandalMetaprograms"))
    o2 <- mpScores(o2, as_tibble = TRUE)
  else if (is_tibble(o2))
    stopifnot("CellID" %in% colnames(o2))

  return (left_join(o1, o2, by = "CellID"))
}

.is_option_requested <- function(opts_vec, opt) (opt %in% names(opts_vec)) && (isTRUE(opts_vec[opt]))

.compute_cor <- function(x, cor_method) {

  x <- center_matrix(x, by = "row", method = "mean", scale = FALSE)

  res <- stats::cor(as.matrix(x), use = "all.obs", method = cor_method)

  return (res)
}

.compute_tsne <- function(x, perplexity = 30, ...) {

  if (3 * perplexity >= ncol(x) - 1) {
    old_perplexity <- perplexity

    perplexity <- (ncol(x) - 1) / 3

    message(paste0("Perplexity ", old_perplexity, " is too large for ", ncol(x), " cells, setting to ", perplexity))
  }

  res <- Rtsne::Rtsne(t(as.matrix(x)), ...)

  return (res$Y)
}

.compute_umap <- function(x, umap_config) {

  res <- umap::umap(t(as.matrix(x)), config = umap_config)

  return (res$layout)
}
