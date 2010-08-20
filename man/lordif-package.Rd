\name{lordif-package}
\alias{lordif-package}
%\alias{lordif}
\docType{package}
\title{
LOgistic Regression Differential Item Functioning using IRT
}
\description{
Analysis of Differential Item Functioning (DIF) for dichotomous and polytomous items  
  using an iterative hybrid of (ordinal) logistic regression and item response theory (IRT).
}
\details{
\tabular{ll}{
Package: \tab lordif\cr
Type: \tab Package\cr
Version: \tab 0.1-6\cr
Date: \tab 2010-08-19\cr
License: \tab GPL (>=2)\cr
LazyLoad: \tab yes\cr
}
  Ordinal logistic regression (OLR) provides a flexible framework for detecting various  
  types of DIF. Previous efforts extended the framework by substituting the matching variable  
  based on sum scores with IRT based trait scores and by employing an iterative process of  
  purifying the matching variable with the use of group-specific item parameters (Crane et. al., 2006).  
  This package represents an effort to integrate both statistical and IRT procedures into  
  a single program. A Monte Carlo simulation approach was incorporated to derive empirical threshold 
  values for various DIF statistics and effect size measures.  
  The two most important functions are: \code{\link{lordif}} and \code{\link{montecarlo}}.
}
\author{
Seung W. Choi, with contributions from  
  Laura E. Gibbons and Paul K. Crane
Maintainer: Seung W. Choi <s-choi@northwestern.edu>
}
\references{
  Crane, P. K., Gibbons, L. E., Jolley, L., & van Belle, G. (2006). Differential item functioning analysis with ordinal logistic regression techniques: DIF detect and difwithpar. Medical Care, 44(11 Suppl 3), S115-S123.
}
\keyword{ package }
\seealso{
  \code{\link[ltm:ltm-package]{ltm}},  \code{\link[Design:Design]{Design}}
}
\examples{
##load PROMIS Anxiety sample data (n=766)
  data(Anxiety)
##age : 0=younger than 65 or 1=65 or older
##gender: 0=Male or 1=Female
##education: 0=some college or higher 1=high school of lower
##run age-related DIF on all 29 items (takes about a minute)
  age.dif <- lordif(Anxiety[paste("R",1:29,sep="")],Anxiety$age) 
##print output
##print(age.dif)
##print extended output
##summary(age.dif)
##generate plots for DIF items (reference group: <65)
##plot(age.dif,labels=c("Younger","Older")) 
##run Monte Carlo simulations for threshold values
##this may take several minutes
  age.dif.MC<-montecarlo(age.dif,alpha=0.05,nr=100)
##print output
##print(age.dif.MC)
##print extended output
##summary(age.dif.MC)
##generate plots for Monte Carlo threshold values
##plot(age.dif.MC)  
}
