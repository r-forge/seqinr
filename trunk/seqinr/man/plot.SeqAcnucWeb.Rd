\name{plot.SeqAcnucWeb}
\alias{plot.SeqAcnucWeb}
\alias{getType.socket}
\title{To plot subsequences on the parent sequence}
\description{
  This function allows to plot all the type of subsequences on a parent sequence. As an example, types could be
  CDS, TRNA, RRNA .... Subsequences are represented by colored rectangle on the parent sequence. 
}
\usage{
plot.SeqAcnucWeb(x, ...)
}
\arguments{
  \item{x}{A sequence of class \code{SeqAcnucWeb}}
  \item{\dots}{graphical parameters can also be specified as arguments}	
}
\details{
  User is asking for the number of the group of subsequences to plot. 	  
 }
\value{
 The plot
}
\references{ \url{http://pbil.univ-lyon1.fr/databases/acnuc.html}}
\author{D. Charif}
\seealso{\code{\link{getType.socket}} }
\examples{
\dontrun{ s=choosebank.socket("hovernucl")}
\dontrun{ query.socket(s$socket,"list","sp=homo sapiens et k=globin@")}	
\dontrun{ plot.SeqAcnucWeb(list$req[[22]]) }
}
\keyword{ hplot }
