\name{comp}
\alias{comp}
\title{ complement a sequence }
\description{
    Complement a sequence.
}
\usage{
comp(seq)
}
\arguments{
  \item{seq}{ an object of class \code{seq} }
}

\value{
  It returns a vector of character which is the complement of a sequence.
}
\author{ D. CHARIF }

\examples{
   a=s2c("acgtgctcgtqtqtqtqtgccc")
   comp(a)		
}

\keyword{ manip }

