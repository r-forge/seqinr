\name{zlibexseq}
\alias{zlibexseq}
\title{Low level function to extract sequences from socket in zlib compressed  mode}
\description{
Low level function to extract sequences from socket in zlib compressed  mode
}

\usage{

zlibexseq(soc,nlz,ver)
}

\arguments{
  \item{soc}{the socket}
  \item{nlz}{number of lines}
  \item{ver}{if \code{TRUE}, verbose mode is on}
}

\details{
 Under development.
}

}
\value{
  Sequence data.
}
\references{ 
\code{citation("seqinr")}
}
\author{S.Penel}
\seealso{
  \code{\link{extractseqs}}, 
}
}
\examples{
 \dontrun{
 # Need internet connection
 choosebank("swissprot")
 query("MyListName", "k=globin")
 MyList.fasta <- exseq("MyListName", verbose=T)
 #
 
 }
}
\keyword{utilities}
