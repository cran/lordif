\name{separate}
\alias{separate}

\title{
  splits item response vectors of DIF items by group
}
\description{
  Splits item response vectors of DIF items into separate group-specific vectors.
}
\usage{
  separate(resp, flag, gr)
}

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
