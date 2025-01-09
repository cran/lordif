\name{DFIT}
\alias{DFIT}
\title{ calculates DFIT statistics }
\description{
  Calculates DFIT statistics using an object of class "lordif"
}
\usage{
  DFIT(obj)
}
\arguments{
  \item{obj}{ an object of class "lordif" }
}
\details{
  Calculates DFIT statistics, including the compensatory differential item functioning (CDIF), 
  the non-compensatory differential item functioning (NCDIF), and the differential test functioning (DTF), 
  based on an object returned from lordif.
}
\value{
  \item{CDIF}{ a data frame of dimension ni by (ng-1), containing compensatory differential item functioning statistics for ni items and (ng-1) groups }
  \item{NCDIF}{ a data frame containing non-compensatory differential item functioning statistics }
  \item{DTF}{ the Differential Test Functioning (DTF) statistic (Raju, van der Linden, & Fleer, 1995) }
  \item{ipar}{ a list of item parameter estimates by group }
  \item{TCC}{ a list of test characteristic functions by group }
}
\references{
  Oshima, T., & Morris, S. (2008). Raju's differential functioning of items and tests (DFIT). Educational Measurement: Issues and Practice, 27, 43-50.
  
  Raju, N. S., van der Linden, W. J., & Fleer, P. F., (1995). An IRT-based internal measure of test bias with application of differential item functioning. Applied Psychological Measurement, 19, 353-368.  
}
\author{ Seung W. Choi <choi.phd@gmail.com> }
\seealso{ \code{\link{lordif}} }
\examples{
##load PROMIS Anxiety sample data (n=766)
\dontrun{data(Anxiety)}
##age : 0=younger than 65 or 1=65 or older
##run age-related DIF on all 29 items (takes about a minute)
\dontrun{age.DIF <- lordif(Anxiety[paste("R",1:29,sep="")],Anxiety$age)}
##run DFIT
\dontrun{age.DIF.DFIT <- DFIT(age.DIF)}
}
