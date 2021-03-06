% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/csim-block-cd.R
\name{forw.cd}
\alias{forw.cd}
\title{A coordinate descent function}
\usage{
forw.cd(w.a, u, y, Tr, X.a, nbasis.t = c(6, 6, 8), rho.grid = c(0,
  0.25, 0.5), it.max = 60, eps = 10^-4, trace = F,
  linear.link = FALSE, ortho.constr = TRUE, i.fx = NULL)
}
\arguments{
\item{w.a}{an initial estimate of the single index coefficient vector.}

\item{u}{a vector of the current singe index variable, u = alpha'X.}

\item{y}{treatment outcomes, n-by-1 vector.}

\item{Tr}{treatment indicators, n-by-1 vector; each element represents one of the K available treatment options.}

\item{X.a}{a matrix of the pre-treatment covarates, associated with the single index coefficients.}

\item{nbasis.t}{a length K+1 vector; each element specifies the number of B-spline basis funtions used in approximating the treatment-specific link function; the last element is for the "main effect" link function.}

\item{rho.grid}{a grid vector of (ridge-type) smoothing parameters for approximating the link functions.}

\item{it.max}{an integer value specifying the maximum number of iterations for each coordinate.}

\item{eps}{a value specifying the converge criterion of algorithm.}

\item{trace}{if \code{TRUE}, show the trace of the fitting procedure; the default is \code{FALSE}.}

\item{linear.link}{if \code{TRUE}, restrict the link functions to be linear functions; the default is \code{FALSE}.}

\item{ortho.constr}{the constraint that separates the interaction effects from the main effect (without this, the interaction effect can be confounded by the main effect); the default is \code{TRUE}.}

\item{i.fx}{a value of the index to be fixed throughout the estimation for model identifiability; the default is \code{NULL}, hence it is estimated.}
}
\value{
a list of information of the optimized single index coefficients including
\item{w.a}{a vector of the estimated single index coefficients.} \item{u}{a n x 1 vector of the estimated single index variable.} \item{it}{the number of iterations.}
}
\description{
\code{forw.cd} implements an ordinary (not a block-) coordinate descent to opitmize the single index coefficients, given an initial single index variable u = alpha'X.
}
\details{
This function can be used to fit an un-regularized constrained single index model. This function is also used in obtaining the "post-selected" (re-fitted) single index coefficients, \code{cfs.pst}; see \code{fit.csim} and \code{block.cd}; see also Radchenko (2015) for the post-selected  single index coefficient estimator.
}
\seealso{
\code{pred.csim},  \code{fit.csim},  \code{fit.csim.cv}
}
