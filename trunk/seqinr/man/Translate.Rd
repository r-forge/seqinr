\name{Translate}
\alias{Translate}
\title{ Generic Function to translate a sequence }
\description{
   Methods for the class of sequences:  \code{SeqFastadna},  \code{SeqFastaAA},  \code{SeqAcnucWeb},  \code{SeqAcnucLocal},  \code{SeqFrag}. 
}
\usage{
Translate(seq, frame = 0, sens = "F", numcode = 1)
}
\arguments{
  \item{object}{ one of this object used to select a method }
}
\value{
 a vector of chars containing the sequence
}
\author{ D. Charif } 
\seealso{ \code{\link{translate}} }
\examples{
 methods(Translate)
}
\keyword{ utilities }