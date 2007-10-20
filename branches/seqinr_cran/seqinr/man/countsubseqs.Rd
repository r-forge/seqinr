\name{countsubseqs}
\alias{countsubseqs}
\alias{css}
\title{Number of subsequences in an ACNUC list}
\description{
Returns the number of subsequences in the ACNUC list of rank \code{lrank}. 
}
\usage{
countsubseqs(lrank, socket = "auto")
css(lrank, socket = "auto")
}
\arguments{
  \item{lrank}{the rank of the ACNUC list to consider.}
  \item{socket}{an object of class \code{connection} returned by 
      \code{choosebank}. The default bank is used when socket == "auto".}
}

\value{
Numeric.
}
\references{ \url{http://pbil.univ-lyon1.fr/databases/acnuc.html}

\code{citation("seqinr")}
}
\author{J.R. Lobry}
\seealso{\code{\link{choosebank}}, \code{\link{query}}, \code{\link{glr}} to
get a list rank from its name.}
\examples{
\dontrun{ 
  ### Need internet connection
  choosebank("emblTP")
  query("mylist", "N=@", virtual=TRUE) # select all (seqs + subseqs)
  mylist$nelem   # 14138094 seqs + subseqs
  css(glr("mylist")) # 1604500 subsequences only
  }
}
\keyword{ utilities }
