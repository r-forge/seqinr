\name{getTrans}
\alias{getTrans}
\alias{getTrans.SeqAcnucWeb}
\alias{getTrans.SeqAcnucLocal}
\alias{getTrans.SeqFastadna}
\alias{getTrans.SeqFastaAA}
\alias{getTrans.SeqFrag}
\alias{getTrans.default}
\alias{translateCDS}
\title{ Generic Function to translate a sequence }
\description{
   Methods for the class of sequences:  \code{SeqFastadna},  \code{SeqFastaAA},  \code{SeqAcnucWeb},  \code{SeqAcnucLocal},  \code{SeqFrag}. 
}
\usage{
getTrans(seq, frame = 0, sens = "F", numcode = 1)
}
\arguments{
  \item{seq}{ an object of class \code{seq}.}
  \item{frame}{ Frame(s) (0,1,2) to translate. By default the frame \code{0} is used. }
  \item{sens}{ Sense to translate: \code{F} for forward sense and \code{R} for revers sense. }
  \item{numcode}{ The number of the code to use. By default the standard genetic code is used. }
}
\value{
 a vector of chars containing the sequence
}
\author{ D. Charif } 
\seealso{ \code{\link{translate}} }
\examples{
 methods(getTrans)
}
\keyword{ utilities } 
