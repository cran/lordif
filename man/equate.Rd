\name{equate}
\alias{equate}
\title{ performs Stocking-Lord Equating }
\description{
  Computes linear transformation constants to equate a set of GRM/GPCM item parameters to a target scale
  using a test characteristic curve equating procedure (Stocking & Lord, 1983)
}
\usage{
  equate(ipar.to, ipar.from, theta, model = "GRM", start.AK = c(1, 0), 
  lower.AK = c(0.5, -2), upper.AK = c(2, 2))
}
\arguments{
  \item{ipar.to}{ a data frame containing target item parameters in the following order: a, cb1, cb2,..., cb(maxCat-1) }
  \item{ipar.from}{ a data frame containing to-be-equated item parameters in the following order: a, cb1, cb2,..., cb(maxCat-1) }
  \item{theta}{ a theta grid }
  \item{model}{ IRT model, either "GRM" or "GPCM" }
  \item{start.AK}{ a vector of starting values, c(A, K) where A is a multiplicative constant and K is an additive constant }
  \item{lower.AK}{ a vector of lower limits, c(A, K) where A is a multiplicative constant and K is an additive constant }
  \item{upper.AK}{ a vector of upper limits, c(A, K) where A is a multiplicative constant and K is an additive constant }
}
\details{
  Computes linear transformation constants (A and K) that equate a set of item parameters (ipar.from) to the scale
  defined by a target item parameters (ipar.to) by minimizing the squared difference between the test characteristic curves
  (Stocking & Lord, 1983). The minimization is performed by the nlminb function (in stats).
}
\value{
  returns a vector of two elements, c(A, K) where A is a multiplicative constant and K is an additive constant
}
\references{
  Stocking, M. L. & Lord, F. M. (1983). Developing a Common Metric in Item Response Theory. Applied Psychological Measurement, 7(2), 201-210.
}
\author{
  Seung W. Choi <choi.phd@gmail.com>
}
\note{
  The item parameters are assumed to be on the theta metric (0,1).
  The number of category threshold parameters may differ across items but not greater than (maxCat-1).
}
\seealso{
  \code{\link{tcc}}
}
\examples{
##ipar.to is a data frame containing "target" item parameters
##ipar.from is a data frame containing "to-be-equated" item parameters
\dontrun{AK <- equate(ipar.to,ipar.from)}
#AK[1] contains the multiplicative constant
#AK[2] contains the additive constant
}
