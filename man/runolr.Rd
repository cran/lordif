\name{runolr}
\alias{runolr}
\title{
  runs ordinal logistic regression models
}
\description{
  Runs ordinal logistic regression models and produces DIF statistics and effect size measures
}
\usage{
runolr(rv, ev, gr, wt)
}
\arguments{
  \item{rv}{ a response variable }
  \item{ev}{ an explanatory variable (e.g., conditioning variable) }
  \item{gr}{ a vector of group identifiers }
  \item{wt}{ a vector of optional sample weights }
}
\details{
  Model 1: ev
  
  Model 2: ev + gr
  
  Model 3: ev*gr or equivalently ev + gr + ev*gr
}
\value{
  Returns a list of the following components:
  \item{chi12}{prob for the LR Chi-square comparing Model 1 vs. Model 2}
  \item{chi13}{prob for the LR Chi-square comparing Model 1 vs. Model 3}
  \item{chi23}{prob for the LR Chi-square comparing Model 2 vs. Model 3}
  \item{beta12}{proportional change in the coefficient for ev}
  \item{pseudo1.CoxSnell}{Cox & Snell psudo R-square for Model 1}
  \item{pseudo2.CoxSnell}{Cox & Snell psudo R-square for Model 2}
  \item{pseudo3.CoxSnell}{Cox & Snell psudo R-square for Model 1}
  \item{pseudo1.Nagelkerke}{Nagelkerke psudo R-square for Model 1}
  \item{pseudo2.Nagelkerke}{Nagelkerke psudo R-square for Model 2}
  \item{pseudo3.Nagelkerke}{Nagelkerke psudo R-square for Model 3}
  \item{pseudo1.McFadden}{McFadden psudo R-square for Model 1}
  \item{pseudo2.McFadden}{McFadden psudo R-square for Model 2}
  \item{pseudo3.McFadden}{McFadden psudo R-square for Model 3}
  \item{pseudo12.CoxSnell}{Cox & Snell R-square change from Model 1 to Model 2}
  \item{pseudo13.CoxSnell}{Cox & Snell R-square change from Model 1 to Model 3}
  \item{pseudo23.CoxSnell}{Cox & Snell R-square change from Model 2 to Model 3}
  \item{pseudo12.Nagelkerke}{Nagelkerke R-square change from Model 1 to Model 2}
  \item{pseudo13.Nagelkerke}{Nagelkerke R-square change from Model 1 to Model 3}
  \item{pseudo23.Nagelkerke}{Nagelkerke R-square change from Model 2 to Model 3}
  \item{pseudo12.McFadden}{McFadden R-square change from Model 1 to Model 2}
  \item{pseudo13.McFadden}{McFadden R-square change from Model 1 to Model 3}
  \item{pseudo23.McFadden}{McFadden R-square change from Model 2 to Model 3}
  \item{df12}{df for the LR Chi-square comparing Model 1 and Model 2}
  \item{df13}{df for the LR Chi-square comparing Model 1 and Model 3}
  \item{df23}{df for the LR Chi-square comparing Model 2 and Model 3}
}
\references{
  Choi, S. W., Gibbons, L. E., Crane, P. K. (2011). lordif: An R Package for Detecting Differential Item Functioning Using Iterative Hybrid Ordinal Logistic Regression/Item Response Theory and Monte Carlo Simulations. Journal of Statistical Software, 39(8), 1-30. URL http://www.jstatsoft.org/v39/i08/.
  
  Crane, P. K., Gibbons, L. E., Jolley, L., & van Belle, G. (2006). Differential item functioning analysis with ordinal logistic regression techniques: DIF detect and difwithpar. Medical Care, 44(11 Suppl 3), S115-S123.
}
\author{
  Seung W. Choi <choi.phd@gmail.com>
}

\seealso{
  \code{\link{rundif}},  \code{\link{lordif}}, \code{\link[rms:lrm]{rms::lrm()}}
}
\examples{
\dontrun{runolr(rv, ev, gr)}
}
