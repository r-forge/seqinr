\name{getLength}
\alias{getLength}
\alias{getLength.SeqAcnucWeb}
\alias{getLength.SeqAcnucLocal}
\alias{getLength.SeqFastadna}
\alias{getLength.SeqFastaAA}
\alias{getLength.SeqFrag}
\alias{getLength.default}
\title{Generic Function To obtain the length of a sequence}
\description{
  Methods for the classes of sequences: \code{SeqFastadna},  \code{SeqFastaAA},  \code{SeqAcnucWeb},  \code{SeqAcnucLocal},  \code{SeqFrag}. 
}
\usage{
getLength(object)
}
\arguments{
   \item{object}{ an object of class \code{SeqFastadna},  \code{SeqFastaAA},  \code{SeqAcnucWeb},  \code{SeqAcnucLocal} or \code{SeqFrag}. }
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
