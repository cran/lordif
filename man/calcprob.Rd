\name{calcprob}
\alias{calcprob}
\title{ calculates item response probabilities }
\description{
  Calculates item response probabilities over a theta grid according to either the GRM or the GPCM.
}
\usage{
  calcprob(ipar, theta, model = "GRM")
}
\arguments{
  \item{ipar}{ a data frame containing the following columns: a, cb1, cb2,..., cb(maxCat-1) }
  \item{theta}{ a grid of theta values, e.g., theta <- seq(-4,4,.1) }
  \item{model}{ IRT model, either "GRM" or "GPCM" }
}
\details{
  Calculates an array of item response probabilities according to either the Graded Response Model (GRM: Samejima, 1969)
  or the Generalized Partial Credit Model (GPCM: Muraki, 1992) over a grid of theta values. 
  The two required input objects are ipar and theta. ipar is a data frame containing
  item parameters in the following order: a, cb1, cb2,..., cb(maxCat-1). Items may have different numbers of
  categories. The variable maxCAT is the maximum number of response categories across all items. 
  theta is a vector containing a grid of theta values. The IRT model can be either "GRM" or "GPCM".
}
\value{
  Returns an array of item response probabilities of dimension, c(nq, ni, maxCAT-1), where
  nq is the length of the theta grid, ni is the number of items in ipar, i.e., nrow(ipar), and maxCAT is the maximum
  number of response categories across all items.
}
\references{ 
  Samejima, F. (1969). Estimation of latent ability using a response pattern of graded scores. Psychometrika Monograph, 17. 

  Muraki, E. (1992). A generalized partial credit model: Application of an EM algorithm. Applied Psychological Measurement, 16, 159-176.
}
\author{ Seung W. Choi <choi.phd@gmail.com> }
\seealso{
	\code{\link{probgrm}}, \code{\link{probgpcm}}
}
\examples{
##item.par<-read.csv(fn,head=F,col.names=c("a","cb1","cb2","cb3","cb4"))
##theta <- seq(-4,4,.1)
\dontrun{calcprob(item.par,theta,model="GPCM")}
}
