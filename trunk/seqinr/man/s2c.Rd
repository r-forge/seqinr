\name{s2c}
\alias{s2c}
\title{ conversion of a string into a vector of chars }
\description{
  This is a simple utility function to convert a string such
as "BigBang" into a vector of chars such as "B" "i" "g" "B" "a" "n" "g"
}
\usage{
s2c(string = "BigBang")
}

\arguments{
  \item{string}{ a string of chars }
}

\value{
  a vector of chars
}
\author{ J.R. Lobry }

\section{WARNING }{ Does not check its argument }

\seealso{ \code{\link{c2s}} }

\examples{
c2s( "BigBang" )
}
\keyword{ Programming }
\keyword{ utilities }


