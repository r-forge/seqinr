\name{getExon}
\alias{getExon}
\alias{getExon.SeqAcnucWeb}
\alias{getExon.SeqAcnucLocal}
\alias{getExon.socket}
\alias{getExon.default}
\title{ Generic Function to get start and end positions of exon }
\description{
   This function allows to get start and end position of exon  whitin parent sequence. Methods for the class of sequences: \code{SeqAcnucLocal} or  \code{SeqAcnucWeb}.
}
\usage{
getExon(object)
}
\arguments{
  \item{object}{ an object of class \code{SeqAcnucLocal} or \code{SeqAcnucWeb} }
}
\value{
  a list of numeric vector giving the position of the first base and the position of the last base for each exon whitin the parent sequence.
}
\author{ D. Charif }
\examples{
 \dontrun{choosebank("embl")}
 \dontrun{query("s","sp=felis catus et t=cds")}
 \dontrun{getExon(s$req[[7]])}	
}
\keyword{ utilities }
