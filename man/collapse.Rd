\name{collapse}
\alias{collapse}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ collapses response categories }
\description{
  Collapses response categories if cell frequencies by group are below a minimum
  threshold specified by minCell and returns a vector of recoded responses.
}
\usage{
  collapse(resp, group, minCell)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{resp}{ a vector of item responses }
  \item{group}{ a vector of group designations }
  \item{minCell}{ a minimum cell frequency }
}
\details{
  Collapses item response categories in resp if the two-way contingency table (resp x group)  
  has cell frequencies less than the minimum cell count threshold specified by minCell.
}
\value{
  a numeric vector of the same length as resp with collapsed/recoded values.
}
\author{ Seung W. Choi <choi.phd@gmail.com> }
\note{ 
  The lowest item response category/score should be coded as 1 not 0 (e.g., 1, 2, 3, 4, 5; and not 0, 1, 2, 3, 4).
  There must be at least two unique categories after collapsing/recoding.
}
\seealso{ \code{\link{recode}} }
\examples{
  r1 <- c(1,1,2,1,1,2,2,1,2,2,1,2,2,1,1,2,1,2,2,3,3,1,2,3)
  gr <- c(0,0,0,1,1,0,1,1,0,0,1,0,1,1,0,1,0,1,0,1,0,1,0,1) 
  collapse(r1,gr,2) #minCell=2
  ## returns  c(1,1,2,1,1,2,2,1,2,2,1,2,2,1,1,2,1,2,2,2,2,1,2,2)
  ## response categories 2 and 3 are collapsed
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
