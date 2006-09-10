\name{stresc}
\alias{stresc}
\title{Utility function to escape LaTeX special characters present in a string}
\description{
 This function returns a string in which LaTeX special characters are
 escaped, this is useful in conjunction with xtable.
}
\usage{
stresc(string = "->\\<->{<->}<->$<->^<->_<->%<->#<->&<->~<-",
     fromchar = s2c("\\{}$^_%#&~"),
     tochar = c("$\\backslash$", "\\{", "\\}", "\\$", "\\^{}",
                                "\\_", "\\%", "\\#", "\\&", "\\~{}"))
}
\arguments{
\itemize{
  \item{string}{The string to deal with.}
  \item{fromchar}{The characters to be escaped.}
  \item{tochar}{The escaped values for the characters to be escaped.}
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
  stopifnot( stresc() == "->$\\backslash$<->\\{<->\\}<->\\$<->\\^{}<->\\_<->\\%<->\\#<->\\&<->\\~{}<-")
}
\keyword{utilities}


