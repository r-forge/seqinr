\name{getKeyword}
\alias{getKeyword}
\alias{getKeyword.SeqAcnucWeb}
\alias{getKeyword.SeqAcnucLocal}
\alias{getKeywordsocket}
\alias{getKeyword.default}
\title{ Generic Function to get keyword(s) associated to a sequence }
\description{
  This function allows to get keyword(s) associated to a sequence. Method for the class of sequences \code{SeqAcnucLocal} or \code{SeqAcnucWeb}. 
}
\usage{
getKeyword(object)
}
\arguments{
  \item{object}{an object of class \code{SeqAcnucLocal} or \code{SeqAcnucWeb}}
}
\value{
   A vector of strings containing the keyword(s) associated to a sequence.
}
\author{D. Charif }
\seealso{ \code{\link{query.socket}} }
\examples{
 \dontrun{s=choosebank.socket("embl")}
 \dontrun{query.socket(s$socket,"fc","sp=felis catus et t=cds")}
 \dontrun{getKeyword(fc$req[[7]])}
}
\keyword{ manip }
