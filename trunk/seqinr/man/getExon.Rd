\name{getExon}
\alias{getExon}
\title{ To get the position of the exon }
\description{
   This function allows to get the position of the exon associated to a CDS from the class of sequence \code{SeqAcnucLocal}.
}
\usage{
getExon(object)
}
\arguments{
  \item{object}{ an object of class \code{SeqAcnucLocal}}
}
\value{
  a list giving the position of the first base and the position of the last base for each exon of the CDS.
}
\author{ D. Charif }
\note{  This function will be soon available for the class of sequence \code{SeqAcnucWeb}}
\examples{
 \dontrun{choosebank("embl")}
 \dontrun{query("s","sp=felis catus et t=cds")}
 \dontrun{getExon(s$req[[7]])}	
}
\keyword{manip}
