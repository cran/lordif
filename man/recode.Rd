\name{recode}
\alias{recode}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ recodes item responses }
\description{
  Recodes item responses as specified by original and modified
}
\usage{
recode(vec, original, modified)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vec}{ a vector of item responses to be recoded }
  \item{original}{ original item response categories, e.g., c(0,1,2,3) }
  \item{modified}{ modified item response categories, e.g., c(1,2,3,4) }
}
\details{
  vec, original, and modified must be of the same mode.
  original and modified must be of the same length.
}
\value{
  Returns a vector of the same length and mode as vec with recoded values.
}
\author{ Seung W. Choi <s-choi@northwestern.edu> }
\examples{
x <- c(0,1,2,3,4,0,1,2,3,4)
y <- c(0,1,2,3,4)
z <- c(1,2,3,4,5)
recode(x,y,z)
##returns c(1,2,3,4,5,1,2,3,4,5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.