\name{lordif}
\alias{lordif}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  performs LOgistic Regression Differential Item Functioning using IRT
}
\description{
  performs iterative hybrid ordinal logistic regression/IRT DIF
}
\usage{
  lordif(resp.data, group, selection = NULL,
  criterion = c("Chisqr", "R2", "Beta"),
  pseudo.R2 = c("McFadden", "Nagelkerke", "CoxSnell"), alpha = 0.01,
  beta.change = 0.1, R2.change = 0.02, maxIter = 10, minCell = 5,
  minTheta = -4, maxTheta = 4, inc = 0.1, NQ=41)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{resp.data}{ data frame or matrix containing item responses }
  \item{group}{ a vector of group designations }
  \item{selection}{ vector specifying a subset of items to be analyzed or NULL for all items }
  \item{criterion}{ criterion for flagging (i.e., "Chisqr", "R2", or "Beta") }
  \item{pseudo.R2}{ pseudo R-squared measure (i.e., "McFadden", "Nagelkerke", or "CoxSnell") }
  \item{alpha}{ significance level for Chi-squared criterion }
  \item{beta.change}{ proportionate change for Beta criterion }
  \item{R2.change}{ R-squared change for pseudo R-squared criterion }
  \item{maxIter}{ maximum number of iterations for purification }
  \item{minCell}{ minimum cell frequency to avoid collapsing }
  \item{minTheta}{ minimum for theta grid }
  \item{maxTheta}{ maximum for theta grid }
  \item{inc}{ increment for theta grid }
  \item{NQ}{ number of quadrature points for IRT parameter estimation (maximum of 61) }
}
\details{
  Performs a ordinal (common odds-ratio) logistic regression DIF analysis using IRT theta estimates  
  as the conditioning variable. The graded response (GR) model is used for IRT trait estimation.  
  Flagged items are treated as separate items and group-specific item parameters are obtained.  
  The procedure runs iteratively until the same set of items is flagged over two iterations.
}
\value{
  Returns a list of class "lordif" with the following components:
  \item{call}{calling expression}
  \item{options}{options used for the run}
  \item{selection}{all or a subset of items analyzed}
  \item{stats}{matrix containing output statistics}
  \item{flag}{logical vector of final flags indicating whether each item is displaying DIF or not}
  \item{recoded}{data frame containing recoded item responses}
  \item{group}{vector of group designation values}
  \item{ng}{scalar for the number of groups}
  \item{ncat}{vector of the number of response categories after collapsing/recoding}
  \item{calib}{vector of theta estimates based on the overall (non-group-specific) item parameters}
  \item{calib.sparse}{vector of theta estimates based on the group-specific item parameters (for DIF items)}
  \item{iteration}{scalar for the number of iterations}
  \item{ipar}{data frame of the overall (non-group-specific) item parameter estimates}
  \item{ipar.sparse}{data frame of the group-specific item parameter estimates}
  \item{stats.raw}{matrix containing output statistics (the same components as stats above but based on raw scores)}
  \item{meanraw}{vector containing mean raw scores}
  \item{flag.raw}{logical vector of DIF final flags based on raw scores}
}
\references{
  Choi, S. W., Gibbons, L. E., Crane, P. K. (2011). lordif: An R Package for Detecting Differential Item Functioning Using Iterative Hybrid Ordinal Logistic Regression/Item Response Theory and Monte Carlo Simulations. Journal of Statistical Software, 39(8), 1-30. URL http://www.jstatsoft.org/v39/i08/.
  
  Crane, P. K., Gibbons, L. E., Jolley, L., and van Belle, G. (2006). Differential item functioning analysis with ordinal logistic regression techniques: DIF detect and difwithpar. Medical Care, 44(11 Suppl 3), S115-S123.
}
\author{
  Seung W. Choi <s-choi@northwestern.edu>
}
\note{
  requires the \pkg{ltm} and \pkg{rms} packages
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{rundif}} }
\examples{
  \dontrun{data(Anxiety)}
  \dontrun{resp.data <- Anxiety[paste("R",1:29,sep="")]}
  \dontrun{age <- Anxiety$age}
  \dontrun{age.DIF <- lordif(resp.data,age)}
  \dontrun{print(age.DIF)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
