\name{s2n}
\alias{s2n}
\title{ simple numerical encoding of a DNA sequence.
 }
\description{
By default, if no \code{levels} arguments is provided, this function will
just code your DNA sequence in integer values following the lexical
order (a > c > g > t), that is 1 for "a", 2 for "c", 3 for "g", 4 for
"t" and NA for ambiguous bases.
}
\usage{
s2n(seq, levels, ...)
}
\arguments{
  \item{seq}{ a vector of chars }
  \item{levels}{ allowed char values, by default a, c, g and t }
  \item{\dots}{ further arguments to factor }
}
\value{
  a vector of integers
}
\author{J.R. Lobry }
\note{ EXPERIMENTAL }


\seealso{ \code{\link{factor}}, \code{\link{codes}} }

\examples{
urn <- c("a","c","g","t")
seq <- sample( urn, 100, replace = TRUE )
seq
s2n(seq)
}
\keyword{ utilities }
