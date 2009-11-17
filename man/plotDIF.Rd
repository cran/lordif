\name{plotDIF}
\alias{plotDIF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ plots DIF graphs for flagged items }
\description{
  plots diagnostic graphs for items identified as displaying DIF
}
\usage{
plotDIF(obj, labels = c("Reference", "Focal"), cexp = 0.8, lwd = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{ output from lordif }
  \item{labels}{ labels for group levels, e.g., c("Male","Female") }
  \item{cexp}{ character expansion factor, e.g., cexp=2 for two times a normal size }
  \item{lwd}{ line width, e.g., lwd=2 for thicker lines for focal groups }
}
\details{
  Generates the following graphs IF there is one or more DIF items:
  1. Trait Distributions - density graphs for groups
  2. Item True Score Functions - true score functions by theta for groups
  3. Differences in Item True Score Functions - unsigned differences
  4. Item Response Functions - item response function for groups
  5. Impact (Weighted by Density) - unsigned differences weighted by theta distributions for the focal group
  6. Comparison of TCCs by group based on all items and DIF items
  7. Impact on Theta Estimates - theta estimates by group before and after accounting for DIF
  
  Graphs 2-5 are generated for each DIF item
  Graphs 1, 6 and 7 are generated for the whole analysis
}
\value{
  Returns no object.
}
\references{ 
  Choi, S. W., Gibbons, L. E., & Crane, P. K. (under review). Development of an iterative hybrid ordinal logistic regression/IRT DIF: A Monte Carlo simulation approach for determining cutoff values. Behavior Research Methods.
}
\author{ Seung W. Choi <s-choi@northwestern.edu> }
\note{ 
  Runs only if obj contains DIF items.
}
\seealso{ \code{\link{lordif}} }
\examples{
##run lordif first
\dontrun{age.dif <- lordif(Anxiety[paste("R",1:29,sep="")],Anxiety$age)}
\dontrun{plotDIF(age.dif,labels=c("Younger","Older"))}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
