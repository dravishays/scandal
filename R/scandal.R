
#' @title scandal: a package for single-cell analysis
#'
#' @description The scandal package provides a framework for performing analysis of
#' single-cell sequencing data and is oriented for analyzing cancer-related data.
#'
#' @section Classes:
#' \describe{
#'   \item{\code{ScandalDataSet}}{The main data structure of the package.}
#'   \item{\code{PreprocConfig}}{Maintains preprocessing configuration parameters.}
#'   \item{\code{QCResults}}{Stores data related to the Quality Control procedure
#'   such as IDs of cells and genes that passed QC, some statistics regardign the
#'   QC process and the \code{PreprocConfig} object used in the process.}
#' }
#'
#' @section Functions:
#' This section lists and gives a short description of the main functions supplied
#' by the scandal package. This is in no way a complete list of the function.
#' \describe{
#'   \item{\code{load_dataset}}{Loads a dataset from file.}
#'   \item{\code{scandal_preprocess}}{Performs preprocessing of a dataset to enable
#'   downstream analysis.}
#'   \item{\code{scandal_plot_qc_metrics}}{Plots the various quality control metric.}
#' }
#'
#' @author Avishay Spitzer
#'
#' @docType package
#' @name scandal
NULL
