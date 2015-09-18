\name{extract}
\alias{extract}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ extracts IRT item parameters }
\description{
  Extracts IRT item parameter estimates from an output returned from \code{\link{mirt}}
}
\usage{
  extract(ipar)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ipar}{ output from the \code{\link{mirt}} function of the mirt package }
}
\details{
  similar to the coef function in the mirt package
}
\value{
  a data frame containing item parameter estimates in the order of 
  a, cb1, cb2,..., cb(maxCat-1).
}
\author{ Seung W. Choi <choi.phd@gmail.com> }
\seealso{ 
  \code{\link{lordif}}, \code{\link{mirt}}
}
\examples{
##calib.sparse <- mirt(sparse.matrix,1,itemtype="gpcm") 
\dontrun{ipar.sparse <- extract(calib.sparse)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
