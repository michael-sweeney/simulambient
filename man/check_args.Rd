% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_args.R
\name{check_args}
\alias{check_args}
\title{Check User Arguments}
\usage{
check_args(
  dataset,
  contamination_levels,
  contamination = NULL,
  metadata = NULL,
  cellTypes = NULL,
  batch = NULL,
  bgBatch = NULL
)
}
\arguments{
\item{dataset}{An n x k matrix-like object where n = number of genes and k = number of nuclei to be decontaminated.}

\item{contamination_levels}{A numeric value or a vector of numeric values indicating the desired levels of contamination that the user would like to simulate. A value of 1 indicates the user wants to simulate the level of contamination detected in the input dataset; a value of 0.5 indicates that the user wants to simulate half of the contamination detected in the original dataset, etc.}

\item{contamination}{An n x j matrix-like object where n = number of genes and j = number of ambient mRNA droplets. The droplets in this matrix are believed to be by the user to be composed of ambient mRNA. Including droplets beneath a specific UMI cutoff (i.e. 100) is a basic idea to get ambient mRNA droplets.}

\item{metadata}{A matrix with k rows that corresponds to metadata of the k nuclei provided in \code{dataset}. Each column should represent a different aspect of metadata, i.e. sex, cell types, etc.}

\item{cellTypes}{A vector of length k should cell type information not be provided in \code{metadata}. If it is provided, this parameter should be a string of the column name that represents the cell type information in \code{metadata}. If cell type information is not desired to be included, default is NULL.}

\item{batch}{A vector of length k that provides the batch of the nucleus in \code{dataset}. If there is no batch effect, parameter can be left blank.}

\item{bgBatch}{A vector of length j corresponding to batch information of the ambient mRNA droplets in \code{contamination}.}
}
\value{
Validated metadata provided by the user.
}
\description{
check_args will take any cellTypes, batch, and sex arguments provided by the user and check them for validity.
}
