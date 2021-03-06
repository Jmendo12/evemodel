% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/twoThetaModel.R
\name{simTwoTheta}
\alias{simTwoTheta}
\title{Simulate expression data for the two-theta EVE model, i.e. two different thetas assigned to specific edges}
\usage{
simTwoTheta(n, tree, colSpecies, isTheta2edge, theta1, theta2, sigma2,
  alpha, beta)
}
\arguments{
\item{n}{Number of "genes" to simulate}

\item{tree}{Phylogeny}

\item{colSpecies}{A character vector with same length as columns in returned expression matrix, 
specifying the species, i.e. tip labels in the phylogeny, for the corresponding column.}

\item{isTheta2edge}{Logical vector with same length as number of edges in tree specifying whether
the corresponding edge theta parameter should be theta2 (TRUE) or theta1 (FALSE)}

\item{theta1}{Value of the theta1 parameter}

\item{theta2}{Value of the theta2 parameter}

\item{sigma2}{Value of the sigma2 parameter}

\item{alpha}{Value of the alpha parameter}

\item{beta}{Value of the beta parameter}
}
\value{
Matrix of simulated gene expression values with samples in columns and genes in rows
}
\description{
Simulate expression data for the two-theta EVE model, i.e. two different thetas assigned to specific edges
}
