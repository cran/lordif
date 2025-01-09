\name{probgrm}
\alias{probgrm}
\title{
  calculates item response probabilities according to GRM
}
\description{
  Calculates a matrix of item response probabilities over a grid of theta values for an item
}
\usage{
probgrm(theta, a, cb)
}
\arguments{
  \item{theta}{ a vector of theta values (e.g., quadrature points) }
  \item{a}{ a slope parameter value }
  \item{cb}{ a vector of category threshold values }
}
\details{
  The Graded Response Model (Samejima, 1969) is assumed.
}
\value{
  Returns a matrix of item response probabilities. The first dimension corresponds to the length of theta.
}
\references{
  Samejima, F. (1969). Estimation of latent ability using a response pattern of graded scores. Psychometrika Monograph, 17.
}
\author{
  Seung W. Choi <choi.phd@gmail.com>
}

\seealso{
  \code{\link{calcprob}}, \code{\link{probgpcm}}
}
\examples{
  probgrm(seq(-4,4,.1), 1.5, c(-1.2,0.5,1.5))
}
