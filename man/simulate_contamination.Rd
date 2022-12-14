% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate_contamination.R
\name{simulate_contamination}
\alias{simulate_contamination}
\title{Simulate snRNA-seq Contamination with Multinomial Distribution}
\usage{
simulate_contamination(geneMeans, numCounts, batch = NULL, seed = 615)
}
\arguments{
\item{geneMeans}{Matrix-like object of ambient mRNA gene expression means for each metadata group. (Rows = genes, columns = metadata groups)}

\item{numCounts}{A vector of length k corresponding to the contamination counts needed for each nucleus.}

\item{batch}{A vector of length k corresponding to batch information of the ambient mRNA droplets in the nucleus \code{dataset}.}

\item{seed}{Seed for random number generation. Can be set by the user. Defaults to 615.}
}
\value{
A sparse matrix of size n x k corresponding to the contamination counts to be added to each simulated nucleus.
}
\description{
simulate_contamination simulates snRNA-seq ambient RNA contamination given parameters for the multinomial distribution.
}
