\name{getSequence}
\alias{getSequence}
\alias{getSequence.SeqAcnucWeb}
\alias{getSequence.SeqAcnucLocal}
\alias{getSequence.SeqFastadna}
\alias{getSequence.SeqFastaAA}
\alias{getSequence.SeqFrag}
\alias{getSequence.default}
\alias{getNumber.socket}
\alias{getseq}
\alias{getSequenceSocket}
\alias{getseq2}

\title{ Generic Function for obtaining a sequence }
\description{
   Methods for the class of sequences:  \code{SeqFastadna},  \code{SeqFastaAA},  \code{SeqAcnucWeb},  \code{SeqAcnucLocal},  \code{SeqFrag}. 
}
\usage{
getSequence(object)
}
\arguments{
  \item{object}{ one of this object used to select a method }
}
\value{
 a vector of chars containing the sequence
}
\author{ D. Charif } 
\seealso{ \code{\link{getFrag}} }
\examples{
 methods(getSequence)
}
\keyword{ utilities }

