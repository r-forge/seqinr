\name{choosebank.socket}
\alias{choosebank.socket}
\title{To select a database structured under ACNUC and located on the web}
\description{
  This function allows to select one of the databases structured under ACNUC and located on the web. Then, you can use \code{query.socket} to compose your query and get a list of sequence names.    
}
\usage{
choosebank.socket(bank = NA, host = "pbil.univ-lyon1.fr", port = 5558))
}
\arguments{
  \item{bank}{ The name of the bank. If NA, it returns, for all database known by the server, its name, availability (off means temporarily unavailable), and description.}
  \item{host}{ Host name for port}
  \item{port}{ The TCP port number}
}
\details{
  \code{choosebank.socket()} returns a list of all the databases located at the Pbil. See the reference section for more details. 
  ACNUC databases are made of a series of flat text files (called divisions) containing annotations and sequences and a series of index files allowing efficient access to sequence data. 
}
\value{
A list with five components
  \item{socket}{ an object of class \code{socket} }
  \item{bankname}{ the name of the bank }
  \item{totseqs}{ the total number of sequences present in the opened database }
  \item{totspecs}{ the total number of species present in the opened database }
  \item{totkeys}{ the total number of keywords present in the opened database }	
}
\references{ 
A description of all the databases located at the pbil is available at the followint url:
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
\seealso{ \code{\link{query.socket}}, \code{\link{connection}} }
\examples{
  \dontrun{choosebank.socket("genbank")}
}
\keyword{ manip }
