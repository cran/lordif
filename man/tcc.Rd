\name{tcc}
\alias{tcc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	computes a test characteristic curve (tcc)
}
\description{
	computes a test characteristic curve (tcc) from input item parameters
}
\usage{
	tcc(a, cb, theta, model = "GRM")
}
%- maybe also 'usage' for other objects documented here.
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
\references{
  Samejima, F. (1969). Estimation of latent ability using a response pattern of graded scores. Psychometrika Monograph, 17.

  Muraki, E. (1992). A generalized partial credit model: Application of an EM algorithm. Applied Psychological Measurement, 16, 159-176.
}
\author{
  Seung W. Choi <choi.phd@gmail.com> 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{probgrm}}, \code{\link{probgpcm}}
}
\examples{
##
  \dontrun{tcc(a,cb,theta,model="GPCM")}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
