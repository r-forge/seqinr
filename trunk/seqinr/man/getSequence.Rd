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

\title{ Generic Function To Get a sequence }
\description{
   Method for the classes of sequences:  \code{SeqFastadna},  \code{SeqFastaAA},  \code{SeqAcnucWeb},  \code{SeqAcnucLocal} and \code{SeqFrag}. 
}
\usage{
getSequence(object)
}
\arguments{
  \item{object}{an object of class \code{SeqFastadna},  \code{SeqFastaAA},  \code{SeqAcnucWeb},  \code{SeqAcnucLocal} or \code{SeqFrag}.}
}
\value{
 A vector of chars containing the sequence
}
\author{ D. Charif } 
\seealso{ \code{\link{getFrag}}, \code{\link{c2s}} }
\examples{
 methods(getSequence)
}
\keyword{ utilities }
