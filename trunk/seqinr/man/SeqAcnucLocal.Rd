\name{SeqAcnucLocal}
\alias{SeqAcnucLocal}
\alias{is.SeqAcnucLocal}
\alias{as.SeqAcnucLocal}
\alias{summary.SeqAcnucLocal}
\title{Sequence coming from local ACNUC data bases}
\description{
    \code{as.SeqAcnucLocal} is called by many functions as \code{query} but not directly by the user. It creates an object of class \code{SeqAcnucLocal}. 
  \code{is.SeqFastadna} returns TRUE if the object is of class \code{SeqAcnucLocal}. 
}
\usage{
as.SeqAcnucLocal(object)
is.SeqAcnucLocal(object)
summary.SeqAcnucLocal(object,...)
}
\arguments{
  \item{object}{ a string giving the name of a sequence presents in the data base}
}
\value{
  \code{as.SeqAcnucLocal} returns an object sequence of class \code{SeqAcnucLocal}
  \code{summary.SeqAcnucLocal} returns a list which the following components:
	\item{name}{the name of th sequence}
	\item{composition}{the base counting of the sequence} 
	\item{GC}{the percentage of G+C in the sequence}
}
\author{ D. Charif }
\examples{
  choosebank()
  query("coli","sp=esherichia coli")
  is.SeqAcnucLocal(coli$req[[1]])
}
\keyword{utilities}
