\name{calctheta}
\alias{calctheta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ calculates EAP theta estimates and associated standard errors }
\description{
  Calculates Expected A Posteriori (EAP) theta estimates and
  associated standard error estimates (posterior standard deviations).
}
\usage{
  calctheta(ipar, resp.data, theta, prior.mean = 0, prior.sd = 1, model = "GRM")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ipar}{ a data frame containing the following columns: a, cb1, cb2,..., cb(maxCat) }
  \item{resp.data}{ a data frame containing item responses }
  \item{theta}{ a theta grid (quadrature points) }
  \item{prior.mean}{ prior mean }
  \item{prior.sd}{ prior standard deviation }
  \item{model}{ IRT model, either "GRM" or "GPCM") }
}
\details{
  Calculates EAP theta estimates and standard error estimates based on the input item parameters (ipar),
  the item response data (resp.data), and the IRT model specified ("GRM" or "GPCM").
}
\value{
  A list object with the following components
  \item{EAP}{Expected A Posteriori estimates of theta}
  \item{SE}{Standard Error estimates}
}
\references{ Bock, R. D. & Mislevy, R. J. (1982). Adaptive EAP Estimation of Ability in  
  a Microcomputer Environment. Applied Psychological Measurement, 6(4), 431-444. }
\author{ Seung W. Choi <choi.phd@gmail.com> }
\note{
  Some missing item responses (NA) are allowed.
}
\seealso{ \code{\link{calcprob}}, \code{\link{probgrm}}, \code{\link{probgpcm}} }
\examples{
  \dontrun{calctheta(ipar,resp.data,model="GPCM")}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
