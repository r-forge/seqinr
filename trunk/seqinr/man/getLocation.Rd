\name{getLocation}
\alias{getLocation}
\alias{getLocation.SeqAcnucWeb}
\alias{getLocation.SeqAcnucLocal}
\alias{getLocation.default}
\alias{getLocationSocket}
\alias{getExon}
\title{ Generic Function for obtaining the position of a sequence on the parent sequence }
\description{
  Methods for the class of sequences: \code{SeqAcnucWeb},  \code{SeqAcnucLocal}. If the sequence is a coding sequence, the function will return the position of each exon on the parent sequence.
}
\usage{
getLocation(object)
}
\arguments{
  \item{object}{ one of this object used to select a method }
}

\value{
  a list giving the position of the sequence on the parent sequence.
}
\author{ D. Charif } 
\seealso{ \code{\link{SeqAcnucWeb}}, \code{\link{SeqAcnucLocal}}  }
\examples{
  methods(getLocation)
}
\keyword{ utilities }
