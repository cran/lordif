\name{extract}
\alias{extract}
\title{ extracts IRT item parameters }
\description{
  Extracts IRT item parameter estimates from an output returned from \code{\link[mirt:mirt]{mirt::mirt()}}
}
\usage{
  extract(ipar)
}
\arguments{
  \item{ipar}{ output from the \code{\link[mirt:mirt]{mirt::mirt()}} function of the \pkg{mirt} package }
}
\details{
  similar to the coef function in the \pkg{mirt} package
}
\value{
  a data frame containing item parameter estimates in the order of 
  a, cb1, cb2,..., cb(maxCat-1).
}
\author{ Seung W. Choi <choi.phd@gmail.com> }
\seealso{ 
  \code{\link{lordif}}, \code{\link[mirt:mirt]{mirt::mirt()}}
}
\examples{
##calib.sparse <- mirt(sparse.matrix,1,itemtype="gpcm") 
\dontrun{ipar.sparse <- extract(calib.sparse)}
}
