
#'
#' @title
#'
#' @description
#'
#' @param p
#' @param show_plot
#' @param save_to_file
#' @param project_dir
#' @param filename
#' @param dirname
#' @param device
#'
#' @details
#'
#' @return
#'
#' @author Avishay Spitzer
#'
#' @export
scandal_plot <- function(p, show_plot = TRUE, save_to_file = FALSE, project_dir = NULL, plots_dir = "plots", filename = NULL, dirname = NULL, device = "png") {

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

    ggplot2::ggsave(filename, plot = p, path = plotting_dir, device = device)
  }

  invisible(p)
}

#' @export
generate_scatter_plot <- function(y_data, x_data = NULL, title = NULL, labels = NULL, name = "", xlab = "Cells", ylab = NULL, plot_ordered = TRUE, scale_y = FALSE, y_axis_breaks = 0.2) {

  if (is.null(y_data))
    stop("y_data is NULL")

  if (isTRUE(plot_ordered) && is.null(x_data)) {
    y_data <- y_data[order(y_data)]

    if (!is.null(labels))
      labels <- labels[names(y_data)]
  }

  if (is.null(x_data))
    x_data <- seq_len(length(y_data))

  if (is.null(labels))
    df <- data.frame(x = x_data, y = y_data)
  else
    df <- data.frame(x = x_data, y = y_data, label = labels)

  if (is.null(labels))
    p <- ggplot2::ggplot(df, aes(x = x, y = y))
  else
    p <- ggplot2::ggplot(df, aes(x = x, y = y, color = label)) +
      ggplot2::scale_color_hue(name = name)

  p <- p + ggplot2::geom_point()

  if (isTRUE(scale_y))
    p <- p +
    ggplot2::scale_y_continuous(breaks = seq(min(y_data), max(y_data), y_axis_breaks), minor_breaks = seq(min(y_data), max(y_data), y_axis_breaks / 10)) +
    ggplot2::theme_minimal()
  else
    p <- p + ggplot2::theme_classic()

  p <- p +
    ggplot2::ggtitle(title) +
    ggplot2::xlab(xlab) +
    ggplot2::ylab(ylab) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size=22))

  return (p)
}

#' @export
generate_histogram_plot <- function(data, title, xlab = "Cells", ylab = NULL, by = .2) {

  if (is.null(data))
    stop("data is NULL")

  df <- data.frame(x = seq_len(length(data)), y = data)

  p <- ggplot2::ggplot(df, aes(x = df$y)) +
    ggplot2::geom_histogram(breaks = seq(0, max(df$y), by = by), col = "black", fill = "lightblue", alpha = .2) +
    #geom_histogram(binwidth = nclass.Sturges, col = "black", fill = "lightblue", alpha = .2) +
    ggplot2::ggtitle(title) +
    ggplot2::theme_classic() +
    ggplot2::xlab(xlab) +
    ggplot2::ylab(ylab) +
    ggplot2::xlim(c(0, max(df$y))) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size=22))

  return (p)
}

#' @export
generate_whiskers_plot <- function(data, labels, title = NULL, xlab = NULL, ylab = NULL, palette = NULL, name = "") {

  df <- data.frame(x = labels, y = data, label = labels)

  p <- ggplot2::ggplot(df, aes(x = x, y = y, fill = label)) +
    ggplot2::geom_boxplot() +
    ggplot2::stat_summary(fun.y = mean, colour = "black", geom = "point", shape = 18, size = 3, show.legend = FALSE) +
    ggplot2::stat_summary(fun.y = mean, colour = "black", geom = "text", show.legend = FALSE, vjust = -0.7, aes(label = round(..y.., digits = 1))) +
    ggplot2::scale_fill_hue(name = name) +
    ggplot2::labs(x = xlab, y = ylab) +
    ggplot2::theme_classic()

  if (!is.null(palette)) {
    if (!(palette %in% rownames(RColorBrewer::brewer.pal.info))) {
      warning(paste0("Unknown palette ", palette, ", setting to default palette (Paired)"))
      palette = "Paired"
    }

    n <- RColorBrewer::brewer.pal.info[palette, "maxcolors"]

    colormap <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(n, palette))(length(table(labels)))
    names(colormap) <- names(table(labels))

    p <- p + ggplot2::scale_fill_brewer(name = name, palette = palette)
  }

  if (!is.null(title)) {
    p <- p +
      ggplot2::labs(title = title) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 22), axis.text.x = ggplot2::element_blank(), axis.ticks.x = element_blank())
  }

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
