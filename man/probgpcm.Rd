\name{probgpcm}
\alias{probgpcm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  calculates item response probabilities according to GPCM
}
\description{
  Calculates a matrix of item response probabilities over a grid of theta values for an item
}
\usage{
  probgpcm(theta, a, cb)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta}{ a vector of theta values (e.g., quadrature points) }
  \item{a}{ a slope parameter value }
  \item{cb}{ a vector of category threshold values }
}
\details{
  The Generalized Partial Credit Model (Muraki, 1992) is assumed.
}
\value{
  Returns a matrix of item response probabilities. The first dimension corresponds to the length of theta.
}
\references{
  Muraki, E. (1992). A generalized partial credit model: Application of an EM algorithm. Applied Psychological Measurement, 16, 159-176.
}
\author{
  Seung W. Choi <choi.phd@gmail.com>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{calcprob}}, \code{\link{probgrm}}
}
\examples{
probgrm(seq(-4,4,.1), 1.5, c(-1.2,0.5,1.5))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
