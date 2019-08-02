
#'
#' @title Plot and save a figure
#'
#' @description Plots a figure to the graphics device and saves it to file if
#' requested.
#'
#' @param p a \code{ggplot} object.
#' @param show_plot whether the plot should be printed to the graphics device.
#' Defualt is TRUE.
#' @param save_to_file whether the plot should be saved to the disk.
#' Default is TRUE.
#' @param project_dir the root directory of the current project to which all
#' plots are saved. Default is NULL meaning that the plots will be saved to the
#' working directory.
#' @param plots_dir a subdirectory in \code{project_dir} to which all plots
#' are saved. Default is "plots".
#' @param filename a name for the file to be saved. If NULL and \code{save_to_file}
#' is set to TRUE then an error will be generated.
#' @param dirname a name for the subdirectory in \code{plots_dir} to save the file into.
#' If NULL and \code{save_to_file} is set to TRUE then an error will be generated.
#' @param device a device to use. See \link{ggsave} for details. Default is "png"
#'
#' @details This function is basically a wrapper for \link{ggsave} that takes care of
#' creating the directories to store plots in using the given arguments.
#'
#' @return Invisibly returns the plot \code{p}.
#'
#' @examples
#' # Generate some data
#' y <- rnorm(1000)
#'
#' # Create an ordered scatter plot of the data
#' p <- scandal_scatter_plot(x = NULL, y = y, plot_ordered = TRUE, order_by_axis = "y")
#'
#' # Will print p to the graphics device and save it to the file "scatter.png" in the current working directory
#' scandal_plot(p, show_plot = TRUE, save_to_file = TRUE, project_dir = ".", filename = "scatter.png", dirname = ".", plots_dir = ".")
#'
#' # In this case p will be saved into a "plots" directory which will be created in the current working directory
#' # without printing p to the graphics device
#' scandal_plot(p, show_plot = FALSE, save_to_file = TRUE, project_dir = ".", filename = "scatter.png", dirname = ".", plots_dir = "plots")
#'
#' @author Avishay Spitzer
#'
#' @import ggplot2
#'
#' @export
scandal_plot <- function(p, show_plot = TRUE, save_to_file = FALSE, project_dir = ".", plots_dir = "plots", filename = NULL, dirname = ".", device = "png") {

  stopifnot(is.ggplot(p), is.logical(show_plot), is.logical(save_to_file))

  if (isTRUE(show_plot))
    print(p)

  if(isTRUE(save_to_file)) {
    stopifnot(is.character(project_dir), is.character(plots_dir), is.character(filename), is.character(dirname), is.character(device))

    plotting_dir <- paste0(project_dir, "/", plots_dir, "/", dirname)

    if (!dir.exists(plotting_dir)) {
      .setup_plotting_dir(project_dir, plots_dir)
      message(paste0("Creating plotting directory - ", plotting_dir))
      base::dir.create(plotting_dir, showWarnings = TRUE)
    }

   ggsave(filename, plot = p, path = plotting_dir, device = device)
  }

  invisible(p)
}

#'
#' @title Create a scatter plot
#'
#' @description This function creates a scatter plot i.e. a two-dimensional plot
#' that uses dots to visualize the values of two different variables \code{x} and
#' \code{y}.
#'
#' @param x a numeric vector representing the variable to plot in the x axis. If
#' set to NULL then a vector \code{1:length(y)} will be used instead.
#' @param y a numeric vector representing the variable to plot in the y axis. If
#' set to NULL then a vector \code{1:length(x)} will be used instead.
#' @param labels an optional vector of character labels for each (x, y) point.
#' Used for color coding each point. Default is NULL.
#' @param color_legend_name the name of the legend, valid only if the labels vector
#' is provided. Default is NULL (can be set to NULL as well if labels vector
#' is provided).
#' @param title an optional string for the title of the plot. Default is NULL
#' (no title).
#' @param xlab a label for the x axis. Default is NULL (no label).
#' @param ylab a label for the y axis. Default is NULL (no label).
#' @param plot_ordered a logical indicating of the provided data should be
#' plotted ordered. Valid only if either \code{x} or \code{y} are provided
#' but not both of them. Default is FALSE.
#' @param order_by_axis the axis by which to order the data (either "x" or "y").
#' Default is "y".
#' @param title_text_size text size of the title. Default is 20.
#'
#' @return A \link{ggplot2} object representing the scatter plot.
#'
#' @examples
#'
#' # Example #1 - One variable of interest
#'
#' # Generate randomly 1000 data points
#' y <- rnorm(1000)
#'
#' # See the difference between the ordered and unordered scatter plots
#' scandal_scatter_plot(x = NULL, y = y, plot_ordered = TRUE, order_by_axis = "y", title = "Ordered plot")
#' scandal_scatter_plot(x = NULL, y = y, plot_ordered = FALSE, title = "Unordered plot")
#'
#' # Example #2 - Two variables with linear correlation
#'
#' library(MASS)
#' library(ggplot2)
#'
#' samples <- 1000
#' r <- 0.8
#'
#' # Generate 2 variables with a linear correlation between them (1000 data points with pearsnon's r 0.8)
#' data <- mvrnorm(n = samples, mu = c(0, 0), Sigma = matrix(c(1, r, r, 1), nrow = 2), empirical = TRUE)
#'
#' # Plot the two variables. As scandal_scatter_plot returns a ggplot object we can add to it a
#' # linear regression line showing the positive correlation
#' scandal_scatter_plot(x = data[, 1], y = data[, 2], plot_ordered = FALSE) +
#'     geom_smooth(method = "glm")
#'
#' @author Avishay Spitzer
#'
#' @export
scandal_scatter_plot <- function(x, y, labels = NULL, color_legend_name = NULL, title = NULL, xlab = NULL, ylab = NULL, plot_ordered = FALSE, order_by_axis = "y", title_text_size = 20) {

  if (is.null(x) & is.null(y))
    stop("Either x or y can be NULL but not both!")

  stopifnot(is.null(x) | (is.vector(x) & is.numeric(x)), is.null(y) | (is.vector(y) & is.numeric(y)))
  stopifnot(is.logical(plot_ordered), is.character(order_by_axis) & order_by_axis %in% c("x", "y"))

  if (!is.null(x) & !is.null(y) & isTRUE(plot_ordered))
    stop("Ordering is not allowed with two variables")

  if (is.null(x))
    x <- seq_len(length(y))
  else if (is.null(y))
    y <- seq_len(length(x))

  stopifnot(length(x) == length(y))
  stopifnot((is.null(labels) | (is.vector(labels) & is.character(labels))) | (!is.null(labels) & (length(labels) == length(x))))
  stopifnot(is.null(color_legend_name) | is.character(color_legend_name),
            is.null(title) | is.character(title),
            is.null(xlab) | is.character(xlab),
            is.null(ylab) | is.character(ylab),
            is.numeric(title_text_size) & title_text_size > 0)

  if (isTRUE(plot_ordered)) {
    if (order_by_axis == "x") {
      ord <- order(x)
      x <- x[ord]
    }
    else {
      ord <- order(y)
      y <- y[ord]
    }

    if (!is.null(labels))
      labels <- labels[ord]
  }

  p <- qplot(x = x, y = y, colour = labels, geom = "point", xlab = xlab, ylab = ylab, main = title) +
    theme_classic() +
    labs(colour = color_legend_name) +
    theme(plot.title = element_text(hjust = 0.5, size = title_text_size))

  return (p)
}

#'
#' @title Create a histogram plot
#'
#' @description This function creates a histogram plot i.e. uses bars to visualize
#' the frequency distribution of \code{data}.
#'
#' @param data a numeric vector of observations.
#' @param title an optional string for the title of the plot. Default is NULL
#' (no title).
#' @param xlab a label for the x axis. Default is NULL (no label).
#' @param ylab a label for the y axis. Default is "Frequency".
#' @param by controls the binwidth. Default is 0.2.
#' @param title_text_size text size of the title. Default is 20.
#' @param colour color of the border of the histogram bars. Default is "black".
#' @param fill fill color of the histogram bars. Default is "lightblue".
#' @param alpha transparency of the histogram bars fill color. Default is 0.2.
#'
#' @return A \link{ggplot2} object representing the histogram plot.
#'
#' @examples
#'
#' # Generate 10K normally distributed data points
#' y <- rnorm(10000, mean = 5, sd = 1)
#'
#' # Generate the histogram plot with bin width of 0.2
#' scandal_histogram_plot(y, by = .2)
#'
#' library(ggplot2)
#'
#' # Generate the histogram plot with bin width of 0.2. Ad the function returns a ggplot object
#' # we can attach to it a vertical line marking the mean of the distribution with regular
#' # ggplot arithmetics
#' scandal_histogram_plot(y, by = .2) +
#'     geom_vline(xintercept = mean(y), colour = "red", linetype = "dashed", size = 1)
#'
#' # Same example with bin width 0.5
#' scandal_histogram_plot(y, by = .5) +
#'     geom_vline(xintercept = mean(y), colour = "red", linetype = "dashed", size = 1)
#'
#' @author Avishay Spitzer
#'
#' @export
scandal_histogram_plot <- function(data, title = NULL, xlab = NULL, ylab = "Frequency", by = .2, title_text_size = 20, colour = "black", fill = "lightblue", alpha = .2) {

  stopifnot(!is.null(data), is.vector(data), is.numeric(data), is.null(dim(data)) | length(dim(data)) == 1)
  stopifnot(is.null(title) | is.character(title),
            is.null(xlab) | is.character(xlab),
            is.null(ylab) | is.character(ylab),
            is.numeric(title_text_size) & title_text_size > 0)


  df <- data.frame(x = seq_len(length(data)), y = data)

  p <- ggplot(df, aes(x = y)) +
    geom_histogram(breaks = seq(0, max(df$y), by = by), colour = colour, fill = fill, alpha = alpha) +
    theme_classic() +
    labs(x = xlab, y = ylab, title = title) +
    theme(plot.title = element_text(hjust = 0.5, size = title_text_size))

  return (p)
}

#'
#' @title Create a box-and-whiskers plot
#'
#' @description This function depicts the distribution of data within groups
#' using a box and whiskers diagram. Each box depicts the interquantile
#' distribution within a group while the whiskers depict the distribution in
#' the upper and lower quartiles.
#'
#' @param data a numeric vector of observations.
#' @param labels a vector of label per observation in \code{data} used to spilt
#' \code{data} into the different groups.
#' @param title an optional string for the title of the plot. Default is NULL
#' (no title).
#' @param xlab a label for the x axis. Default is NULL (no label).
#' @param ylab a label for the y axis. Default is NULL (no label).
#' @param palette an optional color palette used instead of \code{ggplot}'s
#' default palette. Default is NULL.
#' @param legend_name the name of the leged. Default is NULL (no label).
#' @param title_text_size text size of the title. Default is 20.
#'
#' @return A \link{ggplot2} object representing the histogram plot.
#'
#' @examples
#'
#' # Generate 1K normally distributed data points
#' y <- rnorm(1000, mean = 5, sd = 1)
#'
#' # Randomlly assign each data point to a group
#' labels <- paste0("Group", sample(x = 1:5, size = 1000, replace = TRUE))
#'
#' # Plot the distribution
#' scandal_whiskers_plot(y, labels)
#'
#' @author Avishay Spitzer
#'
#' @export
scandal_whiskers_plot <- function(data, labels, title = NULL, xlab = NULL, ylab = NULL, palette = NULL, legend_name = NULL, title_text_size = 20) {

  stopifnot(!is.null(data), is.vector(data), is.numeric(data), is.null(dim(data)) | length(dim(data)) == 1)
  stopifnot((is.null(labels) | (is.vector(labels) & is.character(labels))) | (!is.null(labels) & (length(labels) == length(data))))
  stopifnot(is.null(title) | is.character(title),
            is.null(xlab) | is.character(xlab),
            is.null(ylab) | is.character(ylab),
            is.numeric(title_text_size) & title_text_size > 0)

  df <- data.frame(x = labels, y = data, label = labels)

  p <- ggplot(df, aes(x = x, y = y, fill = label)) +
    geom_boxplot() +
    stat_summary(fun.y = mean, colour = "black", geom = "point", shape = 18, size = 3, show.legend = FALSE) +
    stat_summary(fun.y = mean, colour = "black", geom = "text", show.legend = FALSE, vjust = -0.7, aes(label = round(..y.., digits = 1))) +
    scale_fill_hue(name = legend_name) +
    labs(x = xlab, y = ylab, title = title) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, size = title_text_size), axis.text.x = element_blank(), axis.ticks.x = element_blank())

  if (!is.null(palette)) {
    if (!(palette %in% rownames(RColorBrewer::brewer.pal.info))) {
      warning(paste0("Unknown palette ", palette, ", setting to default palette (Paired)"))
      palette <- "Paired"
    }

    n <- RColorBrewer::brewer.pal.info[palette, "maxcolors"]

    colormap <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(n, palette))(length(table(labels)))
    names(colormap) <- names(table(labels))

    p <- p + scale_fill_brewer(name = legend_name, palette = palette)
  }

  return (p)
}

#'
#' @importFrom reshape2 melt
#' @importFrom scales squish
#'
#' @export
scandal_simple_heatmap_plot <- function(data, center = TRUE, cluster_rows = TRUE, cluster_columns = TRUE, is_corr_matrix = FALSE,
                                        legend_name = "Expression", low = "dodgerblue", high = "red", mid = "white", limits = c(NA, NA), midpoint = 0) {

  stopifnot(!is.null(data), is.matrix(data), is.numeric(data), length(dim(data)) == 2)

  if (isTRUE(center))
    data <- center_matrix(data, by = "row", method = "mean", scale = FALSE)

  if (isTRUE(cluster_columns)) {
    if (isFALSE(is_corr_matrix))
      ord <- scrabble::clusta(mat = t(data))
    else
      ord <- scrabble::clusta(CR = t(data))

    data <- data[ord$ORD, ]
  }

  if (isTRUE(cluster_rows)) {
    if (isFALSE(is_corr_matrix))
      ord <- scrabble::clusta(mat = data)
    else
      ord <- scrabble::clusta(CR = data)

    data <- data[, ord$ORD]
  }

  melted_data <- melt(data)

  p <- ggplot(data = melted_data, aes(Var2, Var1, fill = value)) +
    geom_raster() +
    scale_fill_gradient2(low = low, high = high, mid = mid, limits = limits, midpoint = midpoint, oob = squish, space = "Lab", name = legend_name) +
    theme_void()

  return (p)
}

#'
#' @title Create a t-SNE plot
#'
#' @description This function plots the t-SNE coordinates using a scatter plot.
#'
#' @param object a \links4class{ScandalDataSet} object.
#' @param tsne_labels  an optional vector of character labels for each (x, y) point.
#' Used for color coding each point. Default is NULL.
#' @param legend_name the name of the leged. Default is NULL (no label).
#' @param title an optional string for the title of the plot. Default is NULL
#' (no title).
#' @param title_text_size text size of the title. Default is 20.
#'
#' @return A \link{ggplot2} object representing the histogram plot.
#'
#' @seealso \link{Rtsne}
#'
#' @author Avishay Spitzer
#'
#' @export
scandal_tsne_plot <- function(object, tsne_labels = NULL, legend_name = NULL, title = DEFAULT_TITLE(object, "t-SNE plot"), title_text_size = 20) {

  stopifnot(is_scandal_object(object))
  stopifnot((is.null(tsne_labels) | (is.vector(tsne_labels) & is.character(tsne_labels))) | (!is.null(tsne_labels) & (length(tsne_labels) == ncol(object))))
  stopifnot(is.null(title) | is.character(title),
            is.null(legend_name) | is.character(legend_name),
            is.numeric(title_text_size) & title_text_size > 0)

  tsne_data <- reducedDim(object, "tsne")

  if (is.null(tsne_data))
    stop("t-SNE data not found")

  p <- scandal_scatter_plot(x = tsne_data[, 1],
                            y = tsne_data[, 2],
                            labels = tsne_labels,
                            color_legend_name = legend_name,
                            title = title,
                            xlab = "t-SNE dim1",
                            ylab = "t-SNE dim2",
                            plot_ordered = FALSE,
                            title_text_size = title_text_size)

  return (p)
}

#'
#' @title Create a UMAP plot
#'
#' @description This function plots the UMAP coordinates using a scatter plot.
#'
#' @param object a \links4class{ScandalDataSet} object.
#' @param umap_labels  an optional vector of character labels for each (x, y) point.
#' Used for color coding each point. Default is NULL.
#' @param legend_name the name of the leged. Default is NULL (no label).
#' @param title an optional string for the title of the plot. Default is NULL
#' (no title).
#' @param title_text_size text size of the title. Default is 20.
#'
#' @return A \link{ggplot2} object representing the histogram plot.
#'
#' @seealso \link{umap}
#'
#' @author Avishay Spitzer
#'
#' @export
scandal_umap_plot <- function(object, umap_labels = NULL, legend_name = NULL, title = DEFAULT_TITLE(object, "UMAP plot"), title_text_size = 20) {

  stopifnot(is_scandal_object(object))
  stopifnot((is.null(umap_labels) | (is.vector(umap_labels) & is.character(umap_labels))) | (!is.null(umap_labels) & (length(umap_labels) == ncol(object))))
  stopifnot(is.null(title) | is.character(title),
            is.null(legend_name) | is.character(legend_name),
            is.numeric(title_text_size) & title_text_size > 0)

  umap_data <- reducedDim(object, "umap")

  if (is.null(umap_data))
    stop("UMAP data not found")

  p <- scandal_scatter_plot(x = umap_data[, 1],
                            y = umap_data[, 2],
                            labels = umap_labels,
                            color_legend_name = legend_name,
                            title = title,
                            xlab = "UMAP dim1",
                            ylab = "UMAP dim2",
                            plot_ordered = FALSE,
                            title_text_size = title_text_size)

  return (p)
}

.setup_plotting_dir <- function(project_dir, plots_dir) {
  plotting_root_dir <- paste0(project_dir, "/", plots_dir)

  if (!dir.exists(project_dir)) {
    message(paste0("creating directory ", project_dir))
    dir.create(project_dir, showWarnings = FALSE)
  }

  if (!dir.exists(plotting_root_dir)) {
    message(paste0("creating directory ", plotting_root_dir))
    dir.create(plotting_root_dir, showWarnings = FALSE)
  }
}

DEFAULT_TITLE <- function(object, text) { paste0(nodeID(object), " - ", text) }
