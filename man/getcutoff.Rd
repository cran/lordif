\name{getcutoff}
\alias{getcutoff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ determines a cutoff threshold }
\description{
  Determines an empirical cutoff value for statistics generated from a Monte Carlo simulation.
}
\usage{
  getcutoff(stat, alpha, reverse)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{stat}{ a vector containing statistics sampled from a Monte Carlo simulation }
  \item{alpha}{ a p-value specifying the quantile of the statistics to be determined, e.g., (alpha x 100)th percentile }
  \item{reverse}{ if TRUE, the quantile is determined by (1-alpha) }
}
\details{
  Calculates the quantile value of statistics sampled from a Monte Carlo simulation. 
  For example, when alpha = 0.1 the function determines the 1st percentile of the statistics.
  Conversely, if reverse is True, the function determines the 99th percentile, i.e., (1-alpha) x 100.
}
\value{
  returns a scalar corresponding to the quantile of the statistics determined by alpha or (1-alpha)
}
\author{ Seung W. Choi <choi.phd@gmail.com> }
\seealso{ 
  \code{\link{montecarlo}}, \code{\link{permute}}
}
\examples{
#top 1 percent
getcutoff(runif(1000),0.01,TRUE)
#bottom 1 percent
getcutoff(runif(1000),0.01,FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
