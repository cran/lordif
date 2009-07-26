\name{plotMC}
\alias{plotMC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ plots Monte Carlo output }
\description{
  generate plots for Monte Carlo output
}
\usage{
plotMC(obj, mfrow = c(3, 1), cexp = 1.0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{ output object from lordif }
  \item{mfrow}{ number of rows and columns per page for multi-fane plots }
  \item{cexp}{ character expansion factor, e.g., cexp=2 for two times a normal size }
}
\details{
  Generates the following graphs:
  1. thresholds for Chi-square probability for Model 1 vs. 2
  2. thresholds for Chi-square probability for Model 1 vs. 3
  3. thresholds for Chi-square probability for Model 2 vs. 3
  4. pseudo R-square change from Model 1 to 2
  5. pseudo R-square change from Model 1 to 3
  6. pseudo R-square change from Model 2 to 3
  7. proportional beta change from Model 1 to 2
}
\value{
  Returns no object.
}
\references{ 
  Choi, S. W., Gibbons, L. E., & Crane, P. K. (April, 2009). Development of freeware for an iterative hybrid ordinal logistic regression/IRT DIF. Paper Presented at the National Council of Measurement in Education. San Diego, CA.
}
\author{ Seung W. Choi <s-choi@northwestern.edu> }
\note{ 
  obj is an object output from montecarlo
}
\seealso{ \code{\link{lordif}} }
\examples{
##run lordif first
\dontrun{age.dif <- lordif(Anxiety[paste("R",1:29,sep="")],Anxiety$age)}
\dontrun{plotDIF(age.dif,labels=c("Younger","Older"))}
\dontrun{age.dif.MC <- montecarlo(age.dif,alpha=.05,nr=500)}
\dontrun{plotMC(age.dif.MC,mfrow=c(1,1))}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
