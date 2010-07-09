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
  \method{plot}{lordif.MC}(x, mfrow = c(3, 1), \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ an object of class lordif.MC returned from montecarlo }
  \item{mfrow}{ number of rows and columns per page for multi-fane plots }
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
\author{ Seung W. Choi <s-choi@northwestern.edu> }
\note{
  x is an object of class lordif.MC returned from montecarlo
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{lordif}} 
}
\examples{
  ## s3 plot method for class 'lordif.MC':
  ## age.dif <- lordif(Anxiety[paste("R",1:29,sep="")],Anxiety$age)
  ## age.dif.MC <- montecarlo(age.dif,alpha=.05,nr=500)
  ## plot(age.dif.MC,mfrow=c(1,1))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
