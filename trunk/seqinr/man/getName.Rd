\name{getName}
\alias{getName}
\alias{getName.SeqAcnucWeb}
\alias{getName.SeqAcnucLocal}
\alias{getName.SeqFastadna}
\alias{getName.SeqFastaAA}
\alias{getName.SeqFrag}
\title{ Generic Function for obtaining the name of a sequence }
\description{
  Methods for the class of sequences: \code{SeqFastadna},  \code{SeqFastaAA},  \code{SeqAcnucWeb},  \code{SeqAcnucLocal},  \code{SeqFrag}. 
}
\usage{
getName(object)
}
\arguments{
  \item{object}{ one of this object used to select a method }
}

\value{
  a vector of chars containing the name of the sequence
}
\author{ D. Charif } 
\seealso{ \code{\link{getName}}, \code{\link{getFrag}} }
\examples{
  methods(getName)
}
\keyword{ utilities }
