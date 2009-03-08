\name{getcutoff}
\alias{getcutoff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ determines a cutoff threshold }
\description{
  Determines cutoff thresholds for statistics generated from Monte Carlo simulations
}
\usage{
getcutoff(stat, alpha, reverse)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{stat}{ a vector containing statistics sampled from a Monte Carlo simulation }
  \item{alpha}{ a value determining the proportion to be cut off from the bottom (or top) of stat }
  \item{reverse}{ if TRUE it cuts off the right tail }
}
\details{
  if reversed is TRUE it cuts off the right tail (i.e., largest values)
}
\value{
  returns a scalar that cuts off the designated proportion (specified by alpha)
}
\author{ Seung W. Choi <s-choi@northwestern.edu> }
\seealso{ \code{\link{montecarlo}} }
\examples{
#top 1 percent
getcutoff(runif(1000),0.01,TRUE)
#bottom 1 percent
getcutoff(runif(1000),0.01,FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

