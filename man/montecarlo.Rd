\name{montecarlo}
\alias{montecarlo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  performs Monte Carlo simulations to generate empirical distributions
}
\description{
  performs Monte Carlo simulations under no-DIF conditions to generate empirical distributions of statistics
}
\usage{
montecarlo(obj, alpha = 0.01, nr = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{ an object returned from \code{\link{lordif}} }
  \item{alpha}{ desired significance level (e.g., .01) }
  \item{nr}{ number of replications }
}
\details{
  Simulated datasets are generated under no-DIF conditions and have the same dimensions as the empirical dataset.  
  Group n-sizes and differences (impact) in theta estimates are preserved in the simulated datasets.  
  Returns empirical distributions and thresholds for various statistics and effect size measures.
}
\value{
  Returns an object (list) of class "lordif.MC" with the following components:
  \item{call}{calling expression}
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
\references{
  Choi, S. W., Gibbons, L. E., Crane, P. K. (2011). lordif: An R Package for Detecting Differential Item Functioning Using Iterative Hybrid Ordinal Logistic Regression/Item Response Theory and Monte Carlo Simulations. Journal of Statistical Software, 39(8), 1-30. URL http://www.jstatsoft.org/v39/i08/.
}
\author{
  Seung W. Choi <choi.phd@gmail.com>
}
\note{
  nr must be a large integer (e.g., 500) to generate smooth distributions.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{permute}}, \code{\link{lordif}}
}
\examples{
##load PROMIS Anxiety sample data (n=766)
\dontrun{data(Anxiety)}
##age : 0=younger than 65 or 1=65 or older
##run age-related DIF on all 29 items (takes about a minute)
\dontrun{age.DIF <- lordif(Anxiety[paste("R",1:29,sep="")],Anxiety$age)}
##the following takes several minutes
\dontrun{age.DIF.MC <- montecarlo(age.DIF,alpha=0.01,nr=100)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
