\name{recode}
\alias{recode}
\title{
  recodes item responses
}
\description{
  Recodes item responses as specified by original and modified
}
\usage{
  recode(vec, original, modified)
}
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
\author{
  Seung W. Choi <choi.phd@gmail.com>
}

\examples{
x <- c(0,1,2,3,4,0,1,2,3,4)
y <- c(0,1,2,3,4)
z <- c(1,2,3,4,5)
recode(x,y,z)
##returns c(1,2,3,4,5,1,2,3,4,5)
}
