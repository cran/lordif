\name{lordif}
\alias{lordif}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ performs Logistic Regression Differential Item Functioning using IRT }
\description{
  Performs iterative hybrid ordinal logistic regression/IRT DIF.
}
\usage{lordif(resp.data, group, selection=NULL,  
  criterion="Chisqr", pseudo.R2="McFadden",  
  alpha=0.01, beta.change=0.1, R2.change=0.035,  
  maxIter=10, minCell=5, minTheta=-4.0, maxTheta=4.0, inc=0.1)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{resp.data}{ a data frame or matrix containing item responses }
  \item{group}{ a vector of group designations }
  \item{selection}{ a vector specifying a subset of items to be analyzed or NULL for all items}
  \item{criterion}{ criterion for flagging (i.e., "CHISQR", "R2", or "BETA") }
  \item{pseudo.R2}{ pseudo R-squared measure (i.e., "McFadden", "Nagelkerke", or "CoxSnell") }
  \item{alpha}{ significance level for Chi-squared criterion }
  \item{beta.change}{ proportional change for Beta criterion }
  \item{R2.change}{ R-squared change for pseudo R-squared criterion }
  \item{maxIter}{ maximum number of iterations for purification }
  \item{minCell}{ minimum cell frequncy }
  \item{minTheta}{ minimum for theta grid }
  \item{maxTheta}{ maximum for theta grid }
  \item{inc}{ increment for theta grid }
}
\details{
  Performs a ordinal (common odds-ratio) logistic regression DIF analysis using IRT theta estimates  
  as the conditioning variable. The graded response (GR) model is used for IRT trait estimation.  
  Flagged items are treated as separate items and group-specific item parameters are obtained.  
  The procedure runs iteratively until the same set of items is flagged over two iterations.
}
\value{
  \item{options}{Options used for the run}
  \item{selection}{all or a subset of items analyzed}
  \item{stats}{a matrix containing output statistics}
  \item{flag}{a logical vector of final flags indicating whether each item is displaying DIF or not}
  \item{recoded}{a data frame containing recoded item responses}
  \item{group}{a vector of group designation values}
  \item{ng}{a scalar for the number of groups}
  \item{ncat}{a vector of the number of response categories after collapsing/recoding}
  \item{calib}{a vector of theta estimates based on the overall (non-group-specific) item parameters}
  \item{calib.sparse}{a vector of theta estimates based on the group-specific item parameters (for DIF items)}
  \item{iteration}{a scalar for the number of iterations}
  \item{ipar}{a data frame of the overall (non-group-specific) item parameter estimates}
  \item{ipar.sparse}{a data frame of the group-specific item parameter estimates}
  \item{stats.raw}{a matrix containing output statistics (the same components as stats above but based on raw scores)}
  \item{meanraw}{a vector containing mean raw scores}
  \item{flag.raw}{a logical vector of DIF final flags based on raw scores}
}
\references{ 
  Choi, S. W., Gibbons, L. E., and Crane, P. K. (April, 2009). Development of freeware for an iterative hybrid ordinal logistic regression/IRT DIF. Paper Presented at the National Council of Measurement in Education. San Diego, CA. 
  
  Crane, P. K., Gibbons, L. E., Jolley, L., and van Belle, G. (2006). Differential item functioning analysis with ordinal logistic regression techniques: DIF detect and difwithpar. Medical Care, 44(11 Suppl 3), S115-S123.
}
\author{ Seung W. Choi <s-choi@northwestern.edu> }
\note{ 
  requires the \pkg{ltm} and \pkg{Design} packages
}
\seealso{ \code{\link{rundif}} }
\examples{
##data(Anxeity)
\dontrun{ resp.data <- Anxiety[paste("R",1:29,sep="")] }
\dontrun{ age <- Anxiety$age }
\dontrun{ age.DIF <- lordif(resp.data,age) }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
