\name{calcprob}
\alias{calcprob}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ calculates item response probabilities }
\description{
  Calculates item response probabilities over a theta grid.
}
\usage{
calcprob(ipar, theta)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ipar}{ a data frame containing the following columns: a, cb1, cb2,..., cb(maxCat) }
  \item{theta}{ a grid of theta values, e.g., theta <- seq(-4,4,.1) }
}
\details{
  Calculates an array of item response probabilities according to the graded response model (GRM: Samejima, 1969)
  over a grid of theta values. Two required input objects are ipar and theta. ipar is a data frame containing
  item parameters in the following order: a, cb1, cb2,..., cb(maxCat). Items can have different numbers of
  categories. The variable maxCAT is determined by the function as the maximum number of
  category threshold parameters across all items plus 1. theta is a vector containing a grid of theta values.
}
\value{
  Returns an array of item response probabilities of dimension, c(nq, ni, maxCAT), where
  nq is the length of the theta grid, ni is the number of items in ipar, i.e., nrow(ipar), and maxCAT is the maximum
  number of response categories across all items.
}
\references{ Samejima, F. (1969). Estimation of latent ability using a response pattern of graded scores. Psychometrika Monograph, 17. }
\author{ Seung W. Choi <s-choi@northwestern.edu> }
\seealso{ \code{\link{probgrm}} }
\examples{
##item.par<-read.csv(fn,head=F,col.names=c("a","cb1","cb2","cb3"))
##theta <- seq(-4,4,.1)
\dontrun{calcprob(item.par,theta)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
