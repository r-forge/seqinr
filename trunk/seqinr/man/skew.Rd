\name{skew}
\alias{skew}
\title{ Direct representation of AT and GC skews }
\description{
  The function \code{skew} returns an object of class skew from
a sequence.
}
\usage{
skew(seq, wsizekb = 10, wstepkb = 10)
\method{plot}{skew}(x, option = 1)
}

\arguments{
  \item{seq}{ a vector of DNA characters, if missing a demo
sequence will be used }
  \item{wsizekb}{ the window size in Kb }
  \item{wstepkb}{ the window step in Kb }
  \item{x}{ an object of class skew to be plotted }
  \item{option}{ an integer with value 1 to plot the GC skew and 2 to
plot the AT skew }
}

\value{
An objet of class skew
}
\references{ plot of skew objects are as in figure 1 from Lobry, J.R.
(1996) Asymmetric substitution patterns in the two DNA strands of bacteria.
\emph{Molecular Biology and Evolution},\bold{13}:660-665.\cr
\url{http://mbe.library.arizona.edu/data/1996/1305/4lobr.pdf}\cr}
\author{ J.R. Lobry }
\note{ This kind of direct representation is not very good, consider using
dnawalks instead. }

\section{WARNING }{Experimental !!!!}

\seealso{ \code{\link{oriloc}} }

\examples{
\dontrun{ skw <- skew(); plot(skw) }
}
\keyword{ utilities }
