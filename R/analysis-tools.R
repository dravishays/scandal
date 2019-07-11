
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
                                   tsne_perplexity = 30, tsne_pca = TRUE, tsne_initial_dims = 50, tsne_pca_center = TRUE, tsne_pca_scale = TRUE, tsne_normalize = TRUE,
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

.is_option_requested <- function(opts_vec, opt) (opt %in% names(opts_vec)) && (isTRUE(opts_vec[opt]))

.compute_cor <- function(x, cor_method) {

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
