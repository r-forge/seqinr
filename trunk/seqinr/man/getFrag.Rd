\name{getFrag}
\alias{getFrag}
\alias{getFrag.SeqAcnucWeb}
\alias{getFrag.SeqAcnucLocal}
\alias{getFrag.SeqFastadna}
\alias{getFrag.SeqFastaAA}
\alias{getFrag.SeqFrag}
\title{ Generic Function to obtain a fragment of a sequence }
\description{
  Methods for the class of sequences: \code{SeqFastadna},  \code{SeqFastaAA},  \code{SeqAcnucWeb},  \code{SeqAcnucLocal},  \code{SeqFrag}. 
}
\usage{
getFrag(object, begin, end)
}
\arguments{
  \item{object}{ one of this object used to select a method}
  \item{begin}{ First base }
  \item{end}{ Last base }
}
\value{
  An object of class SeqFrag, which is a vector of chars with many attributes:
  \item{seqMother}{ the name of the sequence from which the sequence comes}
  \item{begin}{ the position of the first base of the fragment on the mother sequence}
  \item{end}{ the position of the last base of the fragment on the mother sequence }
  \item{class}{ \code{SeqFrag} which is the newest class of the sequence } 
}
\author{ D. Charif } 
\seealso{ \code{\link{getSequence}} }
\examples{
  methods(getFrag)
}
\keyword{ utilities }

