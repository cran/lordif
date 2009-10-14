\name{montecarlo}
\alias{montecarlo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ performs Monte Carlo simulations for empirical cutoff thresholds }
\description{
  Performs Monte Carlo simulations using multiple datasets without DIF.
}
\usage{
montecarlo(obj, alpha = 0.01, nr = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{ an object returned from lordif }
  \item{alpha}{ desired significance level (e.g., .01) }
  \item{nr}{ number of replications }
}
\details{
  The simulated datasets have the same dimensions as the empirical data.  
  Group differences (impact) in theta between groups are preserved in simulated datasets.  
  Returns empirical thresholds for various statistics and effect size measures.
}

\value{
  Returns a data frame with the following components:
  \item{chi12}{prob associated with the LR Chi-square test comparing Model 1 vs. 2}
  \item{chi13}{prob associated with the LR Chi-square test comparing Model 1 vs. 3}
  \item{chi23}{prob associated with the LR Chi-square test comparing Model 2 vs. 3}
  \item{pseudo12.CoxSnell}{Cox & Snell pseudo R-square change from Model 1 to 2}
  \item{pseudo13.CoxSnell}{Cox & Snell pseudo R-square change from Model 1 to 3}
  \item{pseudo23.CoxSnell}{Cox & Snell pseudo R-square change from Model 2 to 3}
  \item{pseudo12.Nagelkerke}{Nagelkerke pseudo R-square change from Model 1 to 2}
  \item{pseudo13.Nagelkerke}{Nagelkerke pseudo R-square change from Model 1 to 3}
  \item{pseudo23.Nagelkerke}{Nagelkerke pseudo R-square change from Model 2 to 3}
  \item{pseudo12.McFadden}{McFadden pseudo R-square change from Model 1 to 2}
  \item{pseudo13.McFadden}{McFadden pseudo R-square change from Model 1 to 3}
  \item{pseudo23.McFadden}{McFadden pseudo R-square change from Model 2 to 3}
  \item{beta12}{proportional beta change from Model 1 to 2}
  \item{alpha}{significance level}
  \item{nr}{number of replications}
  \item{cutoff}{thresholds for the statistics}
}
\references{ Choi, S. W., Gibbons, L. E., & Crane, P. K. (April, 2009). Development of freeware for an iterative hybrid ordinal logistic regression/IRT DIF. Paper Presented at the National Council of Measurement in Education. San Diego, CA. }
\author{ Seung W. Choi <s-choi@northwestern.edu> }
\note{ 
  nr must be a large number (e.g., 500) for smooth distributions.
}
\seealso{ \code{\link{lordif}} }
\examples{
##run lordif first
\dontrun{age.dif <- lordif(Anxiety[paste("R",1:29,sep="")],Anxiety$age)}
##the following takes a long time
\dontrun{mc1 <- montecarlo(age.dif,alpha=0.05,nr=500)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
