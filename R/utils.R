
#'
#' @title Saves a data object as an RDS file
#'
#' @description Saves a data object as an RDS files, possibly into a specified sub-directory.
#' Creates the sub-directory if it is missing.
#'
#' @param object an object to be saved as a RDS file.
#' @param filename the name of the file **without .RDS suffix**.
#' @param dirname a sub-directory to which the object will be saved. Default is ".",
#' meaning no sub-directory will be used.
#' @param project_dir the root directory of the current project to which all
#' data files are saved. Default is NULL meaning that the files will be saved to the
#' working directory.
#' @param data_dir a subdirectory in \code{project_dir} to which all data files
#' are saved. Default is "data".
#' @param verbose suppresses all messages from this function. Default is FALSE.
#'
#' @details Saves a data object as an RDS files, possibly into a specified sub-directory.
#' Creates the sub-directory if it is missing.
#'
#' @return Invisibly returns \code{NULL}.
#'
#' @author Avishay Spitzer
#'
#' @export
scandal_save_data <- function(object, filename, dirname = ".", project_dir = ".", data_dir = "data", verbose = FALSE) {

  stopifnot(is.character(filename), is.character(dirname), is.character(project_dir), is.character(data_dir))

  target_dir <- paste0(project_dir, "/", data_dir, "/", dirname)

  if (!dir.exists(target_dir)) {
    .setup_data_dir(project_dir, data_dir)

    if (!dir.exists(target_dir)) {
      message(paste0("Creating data directory - ", target_dir))
      base::dir.create(target_dir, showWarnings = TRUE)
    }
  }

  if (isTRUE(verbose))
    message("Saving file ", paste0(filename, ".RDS"), " into ", target_dir)

  base::saveRDS(object = object, file = paste0(target_dir, "/", filename, ".RDS"))

  invisible(NULL)
}

#'
#' @title Loads an RDS file from the file system
#'
#' @description Loads an RDS file from the file system, possibly from a specified sub-directory.
#'
#' @param filename the name of the file **without .RDS suffix**.
#' @param dirname a sub-directory from which the object will be loaded. Default is ".",
#' meaning no sub-directory will be used.
#' @param project_dir the root directory of the current project to which all
#' data files are saved. Default is NULL meaning that the files will be loaded from the
#' working directory.
#' @param data_dir a subdirectory in \code{project_dir} to which all data files
#' are saved. Default is "data".
#' @param verbose suppresses all messages from this function. Default is FALSE.
#'
#' @details Loads an RDS file from the file system, possibly from a specified sub-directory.
#'
#' @return Returns the data object which was loaded from the RDS file.
#'
#' @author Avishay Spitzer
#'
#' @export
scandal_load_data <- function(filename, dirname = ".", project_dir = ".", data_dir = "data", verbose = FALSE) {

  stopifnot(is.character(filename), is.character(dirname), is.character(project_dir), is.character(data_dir))

  target_dir <- paste0(project_dir, "/", data_dir, "/", dirname)

  if (!dir.exists(target_dir))
    stop("Directory ", target_dir, " does not exist")

  file <- paste0(target_dir, "/", filename, ".RDS")

  if (!file.exists(file))
    stop("File ", filename, ".RDS does not exist in directory ", target_dir)

  res <- base::readRDS(file = file)

  return (res)
}

.setup_data_dir <- function(project_dir, data_dir) {
  data_root_dir <- paste0(project_dir, "/", data_dir)

  if (!dir.exists(project_dir)) {
    message(paste0("creating directory ", project_dir))
    dir.create(project_dir, showWarnings = FALSE)
  }

  if (!dir.exists(data_root_dir)) {
    message(paste0("creating directory ", data_root_dir))
    dir.create(data_root_dir, showWarnings = FALSE)
  }
}
