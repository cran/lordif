\name{plot.lordif}
\alias{plot.lordif}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plot method for lordif class }
\description{
  plots diagnostic graphs for items identified as displaying DIF
}
\usage{
  \method{plot}{lordif}(x, labels = c("Reference", "Focal"), \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ output from lordif }
  \item{labels}{ labels for group levels, e.g., c("Male","Female") }
  \item{...}{ extra graphical parameters }
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
\author{
  Seung W. Choi <s-choi@northwestern.edu>}
\note{
  Produces graphs only if x contains DIF items.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{lordif}} }
\examples{
##run lordif first
  \dontrun{age.dif <- lordif(Anxiety[paste("R",1:29,sep="")],Anxiety$age)}
  \dontrun{plot(age.dif,labels=c("Younger","Older"),cex=0.8,lwd=1)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.