\name{choosebank.socket}
\alias{choosebank.socket}
\title{To select a database structured under ACNUC via a socket protocol}
\description{
  This function allows to select one of the databases structured under ACNUC and located at the \bold{pbil}. Then you can use \code{query.socket} to compose your query and get some sequences.  
}
\usage{
choosebank.socket(bank)
}
\arguments{
  \item{bank}{ The name of the bank.}
}
\details{
  ACNUC databases are made of a series of flat text files (called divisions) containing annotations and sequences and a series of index files allowing efficient access to sequence data.
}
\value{
  \item{socket}{ an object of class socket }
  \item{bankname}{ the name of the bank }
}
\references{ 
\url{http://pbil.univ-lyon1.fr/databases/acnuc.html}\cr
Gouy, M., Milleret, F., Mugnier, C., Jacobzone, M., Gautier,C. (1984) ACNUC: a nucleic acid sequence data base and analysis system. 
\emph{Nucl. Acids Res.}, \bold{12}:121-127.\cr
Gouy, M., Gautier, C., Attimonelli, M., Lanave, C., Di Paola, G. (1985) 
ACNUC - a portable retrieval system for nucleic acid sequence databases:
logical and physical designs and usage.
\emph{Comput. Appl. Biosci.}, \bold{3}:167-172.\cr
Gouy, M., Gautier, C., Milleret, F. (1985) System analysis and nucleic acid sequence banks.
\emph{Biochimie}, \bold{67}:433-436.\cr
For further information about the socket, please see \code{connection} function. For more information about the socket communication protocol with ACNUC please get at \url{http://pbil.univ-lyon1.fr/databases/acnuc/remote_acnuc.html}.

}
\author{ D. Charif }
\note{ In order to know the current bankname, edit the object \code{bankname}.}
\seealso{ \code{\link{query.socket}}, \code{\link{mkbank}}, \code{\link{connection}} }
\examples{
  \dontrun{choosebank.socket("genbank")}
}
\keyword{ manip }
