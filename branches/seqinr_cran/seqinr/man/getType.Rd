\name{getType}
\alias{getType}
\title{ To Get Sequence Type }
\description{
 This function allows to get all the sequence's type present in the opened ACNUC databank.  
}
\usage{
getType(socket)
}
\arguments{
  \item{socket}{ An object of class \code{connection} returned by \code{choosebank}. }
}
\value{
  It returns a list containing for each sequence type, a short description.
}
\references{ ~put references to the literature/web site here ~ }
\author{ D. Charif }
\seealso{ \code{\link{choosebank}}, \code{\link{query}}, \code{\link{plot.SeqAcnucWeb}} }
\examples{
\dontrun{s = choosebank("genbank")}
\dontrun{geType(s$socket) }
}
\keyword{utilities}
