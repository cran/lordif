\name{equate}
\alias{equate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ performs Stocking-Lord Equating }
\description{
  Computes linear transformation constants to equate a set of GRM/2PL item parameters to a target scale
  using a common-item test characteristic curve equating procedure (Stocking & Lord, 1983)
}
\usage{
equate(ipar.to, ipar.from, theta)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ipar.to}{ a data frame containing target item parameters in the following order: a, cb1, cb2,..., cb(ncat-1) }
  \item{ipar.from}{ a data frame containing to-be-equated item parameters in the following order: a, cb1, cb2,..., cb(ncat-1) }
  \item{theta}{ a theta grid (quadrature points) }
}
\details{
  Computes linear transformation constants that equate a set of item parameters (ipar.from) to the scale
  defined by a target item parameters (ipar.to) by minimizing the squared difference between the test characteristic curves
  as described in Stocking and Lord (1983). The minimization is performed by the nlminb function (in stats).
}
\value{
  returns a vector of two elements c(A, X) where A is a multiplicative constant and K is an additive constant
}
\references{ Stocking, M. L. & Lord, F. M. (1983). Developing a Common Metric in Item Response Theory. Applied Psychological Measurement, 7(2), 201-210. }
\author{ Seung W. Choi <s-choi@northwestern.edu> }
\note{ 
  The item parameters are assumed to be on the theta metric (0,1).
  The number of category threshold parameters can differ across items.
}
\seealso{ \code{\link{tcc}} }
\examples{
##ipar.to is a data frame containing "target" item parameters
##ipar.from is a data frame containing "to-be-equated" item parameters
\dontrun{AK <- equate(ipar.to,ipar.from)}
#AK[1] contains the multiplicative constant
#AK[2] contains the additive constant
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.