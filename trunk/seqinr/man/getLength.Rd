\name{getLength}
\alias{getLength}
\alias{getLength.SeqAcnucWeb}
\alias{getLength.SeqAcnucLocal}
\alias{getLength.SeqFastadna}
\alias{getLength.SeqFastaAA}
\alias{getLength.SeqFrag}
\alias{getLength.default}
\title{Generic Function for obtaining the length of a sequence}
\description{
  Methods for the class of sequences: \code{SeqFastadna},  \code{SeqFastaAA},  \code{SeqAcnucWeb},  \code{SeqAcnucLocal},  \code{SeqFrag}. 
}
\usage{
getLength(object)
}
\arguments{
   \item{object}{ one of this object used to select a method }
}
\value{
  a numeric vector giving the length of the sequence
}
\author{ D. Charif } 
\seealso{ \code{\link{getName}}, \code{\link{getFrag}} }
\examples{
 methods(getLength)
}
\keyword{ utilities }
