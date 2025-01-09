\name{permute}
\alias{permute}
\title{
  performs permutation test for empirical cutoff thresholds
}
\description{
  performs permutation tests under no-DIF conditions to generate empirical distributions of DIF statistics
}
\usage{
permute(obj, alpha = 0.01, nr = 100)
}
\arguments{
  \item{obj}{ an object returned from \code{\link{lordif}} }
  \item{alpha}{ desired significance level (e.g., .01) }
  \item{nr}{ number of replications }
}
\details{
  The vector of group designations is randomly shuffled nr times to estimate the sampling distribution
  of the statistics when the null hypothesis is true.
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
  Seung W. Choi <choi.phd@gmail.com>}
\note{
  nr must be a large integer (e.g., 500) for smooth distributions.
}

\seealso{
  \code{\link{montecarlo}}, \code{\link{lordif}}
}
\examples{
##load PROMIS Anxiety sample data (n=766)
\dontrun{data(Anxiety)}
##age : 0=younger than 65 or 1=65 or older
##run age-related DIF on all 29 items (takes about a minute)
\dontrun{age.DIF <- lordif(Anxiety[paste("R",1:29,sep="")],Anxiety$age)}
##the following takes several minutes
\dontrun{age.DIF.MC <- permute(age.DIF,alpha=0.01,nr=100)}
}

