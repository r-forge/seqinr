\name{getKeyword}
\alias{getKeyword}
\alias{getKeyword.SeqAcnucWeb}
\alias{getKeyword.SeqAcnucLocal}
\alias{getKeyword.socket}
\alias{getKeyword.default}
\title{ Generic Function to get keyword associated to a sequence }
\description{
  This function allows to get keyword associated to a sequence. Methods for the class of sequences \code{SeqAcnucLocal} or \code{SeqAcnucWeb}. 
}
\usage{
getKeyword(object)
}
\arguments{
  \item{object}{an object of class \code{SeqAcnucLocal} or \code{SeqAcnucWeb}}
}
\value{
   A vector of strings containing the keyword(s).
}
\author{D. Charif }
\examples{
 \dontrun{choosebank("embl")}
 \dontrun{query("s","sp=felis catus et t=cds")}
 \dontrun{getKeyword(s$req[[7]])}
}
\keyword{ manip }
