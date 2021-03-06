% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/csim-main.R
\name{data.preprocess}
\alias{data.preprocess}
\title{A pre-processing function}
\usage{
data.preprocess(y = NULL, Tr, Xsc = NULL, X1D = NULL, optTr = NULL)
}
\arguments{
\item{y}{treatment outcomes, a n-by-1 vector}

\item{Tr}{treatment indicators, a n-by-1 vector; each element represents one of the K available treatment options}

\item{Xsc}{pre-treatment covarates, a n-by-p matrix}

\item{X1D}{a list of q number of 1 dimensional functional covariates, in which each of the functional covariates is evalauted on some grid points}

\item{optTr}{a n-by-1 vector; each element represents the optimal treatment option for the given covariate information; only available for simulated data.}
}
\value{
\item{y}{a n-by-1 vector of the re-ordered treatment outcomes}
\item{Tr}{a n-by-1 vector of the re-ordered treatment indicators}
\item{Xsc}{a n-by-p matrix of the re-ordered pretreatment covarate matrix}
\item{X1D}{a list of q number of the re-ordered 1 dimensional functional covariates, in which each of the functional covariates is evalauted on some grid points }
}
\description{
A pre-processing  function that can be  used to construct a data frame consist of a set of the observations  in the order of the treatment indicators; i.e., (y, Tr, X) with Tr=1 come first,  (y, Tr, X) with Tr=2 come second, ..., and (y, Tr, X) with Tr=K come last.
}
\seealso{
\code{fit.csim}
}
\author{
Park, Petkova, Tarpey, Ogden
}
