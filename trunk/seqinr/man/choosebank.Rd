\name{choosebank}
\alias{choosebank}
\title{To select a local database structured under ACNUC}
\description{
  This function allows to select one of the available local databases structured under ACNUC. Then you can use \code{query} to compose your query and get some sequences. If the database are not installed in local, you can create your own database with \code{mkbank} and query it. 
}
\usage{
choosebank(bank = "demo")
}
\arguments{
  \item{bank}{ The name of the bank. A local bank of demonstration is set by default if there is no bank name.}
}
\details{
  ACNUC databases are made of a series of flat text files (called divisions) containing annotations and sequences and a series of index files allowing efficient access to sequence data. Two environment variables \emph{acnuc} and \emph{gcgacnuc} are used by all ACNUC programs to define the name of the directories where index and flat files are located, respectively. This function set this two environment variables and open the database specified by \code{bank}.
}
\value{
  it returns a short description of the bank.
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
\emph{Biochimie}, \bold{67}:433-436.
}
\author{ D. Charif }
\note{ In order to know the current bankname, edit the object \code{bankname}.}

\seealso{ \code{\link{choosebank.socket}}, \code{\link{query}}, \code{\link{mkbank}} }
\examples{
  choosebank(bank="demo")
}
\keyword{ manip }
