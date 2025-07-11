% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/svine-array.R
\name{corr_index}
\alias{corr_index}
\title{Helper function, return corresponding index}
\usage{
corr_index(diag_ele, index, is.next = FALSE)
}
\arguments{
\item{diag_ele}{Diagonal elements of a Markov order 1 vine array. A 2d vector.}

\item{index}{The index for which we want to get the corresponding index,
A vector}

\item{is.next}{Boolean. If TRUE, returns the next corresponding index.
Otherwise return the previous index. Default FALSE.}
}
\value{
Corresponding index, a vector with the same length as index.
}
\description{
Helper function, return corresponding index
}
\examples{
A <- c(3, 2, 1, 6, 5, 4)
B <- c(4, 5, 4)
corr_index(diag_ele=A, index=B)

}
