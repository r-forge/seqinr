\name{SeqAcnucWeb}
\alias{SeqAcnucWeb}
\title{Sequence coming from ACNUC data bases }
\description{
    \code{as.SeqAcnucWeb} is called by many functions as \code{query} but not directly by the user. It creates an object of class \code{SeqAcnucWeb}. 
  \code{is.SeqAcnucWeb} returns TRUE if the object is of class \code{SeqAcnucWeb}. 
}
\usage{
as.SeqAcnucWeb(object)
is.SeqAcnucWeb(object)
}
\arguments{
  \item{object}{ a string giving the name of a sequence presents in the data base}
}
\value{
  \code{as.SeqAcnucWeb} returns an object sequence of class \code{SeqAcnucWeb}
}
\author{ D. Charif }
\examples{
  s = choosebank.socket("genbank")
  query.socket(s$socket,"felis","sp=felis catus et t=cds et o=mitochondrion")
  is.SeqAcnucWeb(felis$req[[1]])
}
\keyword{utilities}
