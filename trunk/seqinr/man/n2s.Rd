\name{n2s}
\alias{n2s}
\title{ function to convert the numeric encoding of a DNA sequence
into a vector of characters }
\description{
By default, if no `levels' arguments is provided, this function
will just transform your vector of integer into a DNA sequence 
according to the lexical order: 1 -> "a", 2 -> "c", 3 -> "g",
4 -> "t", others -> NA
}
\usage{
n2s(nseq, levels = c("a", "c", "g", "t"))
}
\arguments{
  \item{nseq}{a vector of integers }
  \item{levels}{the translation vector }
}
\value{
a vector of characters
}
\author{ J.R. Lobry }
\note{EXPERIMENTAL}

\seealso{ \code{\link{s2n}} }

\examples{
nseq <- sample(x = 1:4, size = 100, replace = TRUE)
n2s(nseq) 
# Show what's happen with out-of-range and NA values :
nseq[1] <- NA
nseq[2] <- 777
n2s(nseq)[1:10]
}
\keyword{ utilities }
