\name{stresc}
\alias{stresc}
\title{Utility function to escape LaTeX special characters present in a string}
\description{
 This function returns a string in which LaTeX special characters are
 escaped, this is useful in conjunction with xtable.
}
\usage{
stresc(string)
}
\arguments{
\itemize{
  \item{string}{The string to deal with.}
  }
}
\value{
  Returns a string with escaped characters.
}
\references{
  \code{citation("seqinr")}
}
\author{J.R. Lobry}
\seealso{ \code{\link{s2c}} }
\examples{
  stresc("MISC_RNA")
}
\keyword{utilities}


