\name{closebank.socket}
\alias{closebank.socket}
\title{ close the bank opening with choosebank.socket }
\description{
  It closes  the bank opening with \code{choosebank.socket}.
}
\usage{
closebank.socket(socket)
}
\arguments{
  \item{socket}{ An object of class socket coming from \code{choosebank.socket}} 
}
\value{
  THe following messsage is printed \code{OK ACNUC SOCKET STOPED}
}
\author{D .Charif }
\seealso{\code{\link{choosebank.socket}} }
\examples{
 \dontrun{s=choosebank.socket("genbank")}
 \dontrun{closebank.socket(s$socket)}
}
\keyword{utilities}
