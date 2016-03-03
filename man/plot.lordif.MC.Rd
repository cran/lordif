\name{plot.lordif.MC}
\alias{plot.lordif.MC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Plot method for Monte Carlo simulation output
}
\description{
  Produces plots for Monte Carlo output
}
\usage{
  \method{plot}{lordif.MC}(x, mfrow = c(3, 1), width = 7, height = 7, \dots)
%plot.lordif.MC(x, mfrow = c(3, 1), width = 7, height = 7, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ an object of class lordif.MC returned from \code{\link{montecarlo}} }
  \item{mfrow}{ number of rows and columns per page for multi-fane plots }
  \item{width}{ width of the canvas of the plotting device in inches }
  \item{height}{ height of the canvas of the plotting device in inches }
  \item{...}{ extra graphical parameters }
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
  Choi, S. W., Gibbons, L. E., Crane, P. K. (2011). lordif: An R Package for Detecting Differential Item Functioning Using Iterative Hybrid Ordinal Logistic Regression/Item Response Theory and Monte Carlo Simulations. Journal of Statistical Software, 39(8), 1-30. URL http://www.jstatsoft.org/v39/i08/.

}
\author{
  Seung W. Choi <choi.phd@gmail.com>
}
\note{
  x is an object of class lordif.MC returned from \code{\link{montecarlo}}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{lordif}}, \code{\link{montecarlo}}, \code{\link{permute}}
}
\examples{
## s3 plot method for class 'lordif.MC':
\dontrun{age.dif <- lordif(Anxiety[paste("R",1:29,sep="")],Anxiety$age)}
\dontrun{age.dif.MC <- montecarlo(age.dif,alpha=.05,nr=500)}
\dontrun{plot(age.dif.MC,mfrow=c(1,1),width=8,height=7)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
