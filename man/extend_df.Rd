% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_svine_func.R
\name{extend_df}
\alias{extend_df}
\title{Helper function for fit_svine}
\usage{
extend_df(dat, var_order, k = 1)
}
\arguments{
\item{dat}{N x d matrix; each column corresponds to a
univariate vector}

\item{var_order}{a vector that contains the desired order of the variable
(If var_order = 1:d, it means the order of the variable stays the same
as the original order in dat)}

\item{k}{Markov order k >= 1, default is k = 1}
}
\value{
an extended (N−k)×d(k+1) matrix with variables in the order specified
        by var_order and Markov order k
}
\description{
Helper function for fit_svine
}
\examples{
test_d3 = matrix(data = runif(30), ncol = 3, nrow = 10)
extend_df(test_d3, var_order = 1:3)

}
