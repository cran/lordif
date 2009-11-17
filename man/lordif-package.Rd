\name{lordif-package}
\alias{lordif-package}
%\alias{lordif}
\docType{package}
\title{
Logistic Regression Differential Item Functioning using IRT
}
\description{
Analysis of Differential Item Functioning (DIF) for dichotomous and polytomous items  
  using an iterative hybrid of (ordinal) logistic regression and item response theory (IRT).
}
\details{
\tabular{ll}{
Package: \tab lordif\cr
Type: \tab Package\cr
Version: \tab 0.1-4\cr
Date: \tab 2009-11-12\cr
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
  The three most important functions are: \code{\link{lordif}}, \code{\link{plotDIF}}, and \code{\link{montecarlo}}.
}
\author{
Seung W. Choi, with contributions from  
  Laura E. Gibbons and Paul K. Crane

Maintainer: Seung W. Choi <s-choi@northwestern.edu>
}
\references{
  Choi, S. W., Gibbons, L. E., & Crane, P. K. (under review). Development of an iterative hybrid ordinal logistic regression/IRT DIF: A Monte Carlo simulation approach for determining cutoff values. Behavior Research Methods. 
  
  Crane, P. K., Gibbons, L. E., Jolley, L., & van Belle, G. (2006). Differential item functioning analysis with ordinal logistic regression techniques: DIF detect and difwithpar. Medical Care, 44(11 Suppl 3), S115-S123.
}
\keyword{ package }
\seealso{
  \code{\link[ltm:ltm-package]{ltm}},  \code{\link[Design:Design]{Design}}
}
\examples{
##data(Anxiety)
##age : 0=younger than 65 and 1=65 or older
##gender: 0=Male and 1=Female
##education: 0=some college or higher 1=high school of lower
#runs age DIF on all 29 items (takes about a minute)
\dontrun{age.dif <- lordif(Anxiety[paste("R",1:29,sep="")],Anxiety$age)} 
#generates plots for DIF items, reference="Younger" & focal="Older"
\dontrun{plotDIF(sample.run,labels=c("Younger","Older"))} 
\dontrun{montecarlo(age.dif,alpha=0.05,nr=500)} 
##this may take a few minutes
}
