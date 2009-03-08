\name{calctheta}
\alias{calctheta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ calculates EAP theta and associated standard error }
\description{
  Calculates the Expected A Posteriori (EAP) theta and the
  associated posterior standard deviation (PSD) as standard error estimates.
}
\usage{
calctheta(ipar, resp.data, theta, prior.mean = 0, prior.sd = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ipar}{ a data frame containing the following columns: a, cb1, cb2,..., cb(maxCat) }
  \item{resp.data}{ a data frame containing item responses }
  \item{theta}{ a theta grid (quadrature points) }
  \item{prior.mean}{ prior mean }
  \item{prior.sd}{ prior standard deviation }
}
\details{
  Calculates EAP theta estimates and PSD standard error estimates based on the input item parameters (ipar)
  and the item response data (resp.data).
}
\value{
  A list object with the components
  \item{EAP}{Expected A Posteriori estimates of theta}
  \item{SE}{Standard Error estimates}
}
\references{ Bock, R. D. & Mislevy, R. J. (1982). Adaptive EAP Estimation of Ability in  
  a Microcomputer Environment. Applied Psychological Measurement, 6(4), 431-444. }
\author{ Seung W. Choi <s-choi@northwestern.edu> }
\note{
 Not all item responses may be present.
}
\seealso{ \code{\link{calcprob}},  \code{\link{probgrm}} }
\examples{
\dontrun{calctheta(ipar,resp.data)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
