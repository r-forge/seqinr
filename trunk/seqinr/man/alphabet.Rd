\name{alphabet}
\alias{alphabet}
\title{To display alphabet of specified length }
\description{
  This function displays all the words of the specified length that you can do with the alphabet (A,C,G,T).       
}
\usage{
alphabet(al)
}
\arguments{
  \item{al}{is the length of the words}
}
\value{
  \code{alphabet} return an array.  
}
\author{ D .Charif } 
\seealso{ \code{\link{count}} }
\examples{
alphabet(2)
alphabet(5)
as.vector(alphabet(3))
}
\keyword{utilities}
