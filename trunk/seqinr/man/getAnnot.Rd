\name{getAnnot}
\alias{getAnnot}
\alias{getAnnot.SeqAcnucWeb}
\alias{getAnnot.SeqAcnucLocal}
\alias{getAnnot.SeqFastadna}
\alias{getAnnot.SeqFastaAA}
\title{ Generic Function to obtain annotation for a sequence }
\description{
  Methods for the class of sequences: \code{SeqFastadna},  \code{SeqFastaAA},  \code{SeqAcnucWeb},  \code{SeqAcnucLocal}. 
}
\usage{
getAnnot(object, nl)
}
\arguments{
  \item{object}{ one of this object used to select a method}
  \item{nl}{ the number of line of annotation to read}
}
\value{
  a vector of sting of character containing the annotation. 
}
\author{ D. Charif } 
\seealso{ \code{\link{getProp}} }
\examples{
  methods(getAnnot)
}
\keyword{ utilities }
