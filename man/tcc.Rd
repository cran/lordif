\name{tcc}
\alias{tcc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ computes a test characteristic curve (tcc) }
\description{
  computes a test characteristic curve (tcc) from input item parameters
}
\usage{
tcc(a, cb, theta)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{a}{ a vector of slope parameters }
  \item{cb}{ a matrix of category boundary/threshold parameters }
  \item{theta}{ a grid of theta values }
}
\details{
  The graded response model is assumed.
  The lowest score for each item is set to 0 not 1. 
  length(a) must be the same as nrow(cb)
}
\value{
  Retruns a vector representing a tcc
}
\references{ Samejima, F. (1969). Estimation of latent ability using a response pattern of graded scores. Psychometrika Monograph, 17. }
\author{ Seung W. Choi <s-choi@northwestern.edu> }
\seealso{ \code{\link{equate}}, \code{\link{probgrm}} }
\examples{
##
\dontrun{tcc(a,cb,theta)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
