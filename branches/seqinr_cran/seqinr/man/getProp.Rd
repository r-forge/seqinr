\name{getProp}
\alias{getProp}
\alias{getProp.SeqAcnucWeb}
\alias{getProp.SeqAcnucLocal}
\alias{getProp.SeqFastadna}
\alias{getProp.SeqFastaAA}
\alias{getProp.SeqFrag}
\alias{getProp.default}
\alias{AAprop}
\alias{getAttribut.socket}
\alias{getAttribut}
\title{ Generic Function to obtain some properties of a sequence }
\description{
  Method for the classes of sequence: \code{SeqFastadna},  \code{SeqFastaAA},  \code{SeqAcnucWeb},  \code{SeqAcnucLocal},  \code{SeqFrag}. 
}
\usage{
getProp(object)
}
\arguments{
  \item{object}{ an object of class code{SeqFastadna},  \code{SeqFastaAA},  \code{SeqAcnucWeb},  \code{SeqAcnucLocal} or \code{SeqFrag} }
}
\value{
  A list with two components:
  \item{frame }{the reading frame for the sequence in the case of a CDS}
  \item{code }{the genetic code to use in order to translate the sequence (number of genetic code according to ACNUC) }
}
\author{ D. Charif }
\section{Warning }{ Just available for coding sequences } 
\seealso{ \code{\link{getName}}, \code{\link{getFrag}}  }
\examples{
 methods(getProp)
}
\keyword{ utilities }
