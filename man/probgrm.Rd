\name{probgrm}
\alias{probgrm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  calculates item response probabilities according to GRM
}
\description{
  Calculates a matrix of item response probabilities over a grid of theta values for an item
}
\usage{
probgrm(theta, a, cb)
}
%- maybe also 'usage' for other objects documented here.
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

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{calcprob}}, \code{\link{probgpcm}}
}
\examples{
  probgrm(seq(-4,4,.1), 1.5, c(-1.2,0.5,1.5))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
