\name{rundif}
\alias{rundif}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ runs ordinal logistic regression DIF }
\description{
  Runs ordinal logistic regression DIF
}
\usage{
  rundif(item, resp, theta, gr, criterion, alpha, beta.change, 
         pseudo.R2, R2.change)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{item}{ a selection of items to be analyzed }
  \item{resp}{ a data frame containing item responses }
  \item{theta}{ a conditioning (matching) variable }
  \item{gr}{ a vector of group identifiers }
  \item{criterion}{ criterion for flagging (i.e., "CHISQR", "R2", or "BETA") }
  \item{alpha}{ significance level for Chi-squared criterion }
  \item{beta.change}{ proportional change for Beta criterion }
  \item{pseudo.R2}{ pseudo R-squared measure (i.e., "McFadden", "Nagelkerke", or "CoxSnell") }
  \item{R2.change}{ R-squared change for pseudo R-squared criterion }
}
\details{
  The argument item lists the column numbers of the data frame resp to be included in the analysis.
}
\value{
  Returns a list of the following components:
  \item{stats}{a data frame containing output statistics}
  \item{flag}{a logical vector of DIF flags}
}
\references{
  Crane, P. K., Gibbons, L. E., Jolley, L., and van Belle, G. (2006). Differential item functioning analysis with ordinal logistic regression techniques: DIF detect and difwithpar. Medical Care, 44(11 Suppl 3), S115-S123.
}
\author{ Seung W. Choi <s-choi@northwestern.edu> }
\seealso{ \code{\link{runolr}},  \code{\link{lordif}} }
\examples{
##
\dontrun{rundif(item,resp,theta,gr)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.