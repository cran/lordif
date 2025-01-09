\name{tcc}
\alias{tcc}
\title{
  computes a test characteristic curve (tcc)
}
\description{
  computes a test characteristic curve (tcc) from input item parameters
}
\usage{
  tcc(a, cb, theta, model = "GRM")
}

\arguments{
  \item{a}{ a vector of slope parameters }
  \item{cb}{ a matrix of category boundary/threshold parameters }
  \item{theta}{ a grid of theta values }
  \item{model}{ IRT model, either "GRM" or "GPCM" }
}
\details{
  Either the graded response model (GRM) or the generalized partial credit model (GPCM) can be specified.
  The lowest response/score for each item is coded as 0 not 1.
}
\value{
  Retruns a vector of tcc values over a theta grid
}

\author{
  Seung W. Choi <choi.phd@gmail.com>
}

\seealso{
  \code{\link{probgrm}}, \code{\link{probgpcm}}
}
\examples{
\dontrun{tcc(a,cb,theta,model="GRM")}
\dontrun{tcc(a,cb,theta,model="GPCM")}
}
