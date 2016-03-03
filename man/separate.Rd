\name{separate}
\alias{separate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  splits item response vectors of DIF items by group
}
\description{
  Splits item response vectors of DIF items into separate group-specific vectors.
}
\usage{
  separate(resp, flag, gr)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{resp}{ a data frame (or matrix) of item responses }
  \item{flag}{ a vector of DIF flags }
  \item{gr}{ a vector of group identifiers }
}
\details{
  To obtain group specific item calibrations, response vectors of DIF items are split into
  multiple vectors by group and treated as separate items.
}
\value{
  Returns a data frame with item response vectors for non-DIF items followed by
  separated item response vectors for DIF items
}

\author{
  Seung W. Choi <choi.phd@gmail.com>
}

\examples{
\dontrun{separate(resp,flag,gr)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
