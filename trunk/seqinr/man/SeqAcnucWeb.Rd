\name{SeqAcnucWeb}
\alias{SeqAcnucWeb}
\alias{as.SeqAcnucWeb}
\alias{is.SeqAcnucWeb}
\title{Sequence coming from ACNUC data bases located on the web}
\description{
    \code{as.SeqAcnucWeb} is called by many functions as \code{query} but not directly by the user. It creates an object of class \code{SeqAcnucWeb}. 
  \code{is.SeqAcnucWeb} returns TRUE if the object is of class \code{SeqAcnucWeb}. 
}
\usage{
as.SeqAcnucWeb(object,socket=F)
is.SeqAcnucWeb(object)
}
\arguments{
  \item{object}{ a string giving the name of a sequence presents in the data base}
  \item{socket}{ an object of class \code{socket}}
}
\value{
  \code{as.SeqAcnucWeb} returns an object sequence of class \code{SeqAcnucWeb}
}
\author{ D. Charif }
\examples{
  \dontrun{s = choosebank.socket("genbank")}
  \dontrun{query.socket(s$socket,"felis","sp=felis catus et t=cds et o=mitochondrion")}
  \dontrun{is.SeqAcnucWeb(felis$req[[1]])}
}
\keyword{utilities}
