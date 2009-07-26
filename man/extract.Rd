\name{extract}
\alias{extract}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ extracts item parameters }
\description{
  extracts IRT item parameter estimates from an output returned from ltm:grm 
}
\usage{
extract(ipar)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ipar}{ output from the grm function in the ltm package }
}
\details{
  similar to the coef function in the ltm package
}
\value{
  a data frame containing item parameter estimates in the order of 
  a, cb1, cb2,..., cb(maxCat-1).
}
\author{ Seung W. Choi <s-choi@northwestern.edu> }
\seealso{ \code{\link{lordif}} }
\examples{
##calib.sparse <- grm(sparse.matrix,constrained=FALSE,IRT.param=TRUE) 
\dontrun{ipar.sparse <- extract(calib.sparse)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.