\name{getKeyword}
\alias{getKeyword}
\title{ To get the keyword associated to a sequence }
\description{
  This function allows to get the keyword associated to a sequence of class SeqAcnucLocal. 
}
\usage{
getKeyword(object)
}
\arguments{
  \item{object}{an object of class \code{SeqAcnucLocal}
}
\value{
   A vector of strings containing the keyword(s).
}
\author{D. Charif }
\note{  This function will be soon available for the class of sequence \code{SeqAcnucWeb}}
\examples{
 \dontrun{choosebank("embl")}
 \dontrun{query("s","sp=felis catus et t=cds")}
 \dontrun{getKeyword(s$req[[7]])}
}
\keyword{ manip }
