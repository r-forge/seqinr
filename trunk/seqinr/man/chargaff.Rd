\name{chargaff}
\alias{chargaff}
\docType{data}
\title{Base composition in ssDNA for 7 bacterial DNA}
\description{
  ~~ A concise (1-5 lines) description of the dataset. ~~
  Sould be completed with accuracy.
}
\usage{data(chargaff)}
\format{
  A data frame with 7 observations on the following 4 variables.
  \describe{
    \item{[A]}{frequencies of A bases in percent}
    \item{[G]}{frequencies of G bases in percent}
    \item{[C]}{frequencies of C bases in percent}
    \item{[T]}{frequencies of T bases in percent}
  }
}
\details{
  ~~ If necessary, more details than the __description__ above ~~
}
\source{
Rudner, R., Karkas, J.D., Chargaff, E. (1969) Separation of microbial deoxyribonucleic acids into complementary strands. \emph{Proceedings of the National Academy of Sciences of the United States of America}, \bold{63}:152-159.
}
\references{
  ~~ possibly secondary sources and usages ~~ To be completed
  Figures not checked yet!

}

\examples{
data(chargaff)
op <- par(no.readonly = TRUE)
par(mfrow=c(4,4))
xlim <- c(0,100)
ylim <- xlim
par(mai=rep(0,4))
par(omd=rep(0.01,4))
par(xaxs="i")
par(yaxs="i")

for( i in 1:4 )
{
  for( j in 1:4 )
  {
    if( i == j )
    {
      plot(chargaff[,i], chargaff[,j],t="n", xlim=xlim, ylim=ylim,
      xlab="", ylab="", xaxt="n", yaxt="n")
      polygon(x=c(0,0,100,100),y=c(0,100,100,0), col="lightgrey")
      for( k in seq(0,100,by=10) )
      {
        lseg <- 3
        segments(k,0,k,lseg)
        segments(k,100-lseg,k,100)
        segments(0,k,lseg,k)
        segments(100-lseg,k,100,k)
      }
      string <- paste(names(chargaff)[i],"\n\n",xlim[1],"\% -",xlim[2],"\%")
      text(x=mean(xlim),y=mean(ylim), string, cex = 1.5)
    }
    else
    {
      plot(chargaff[,i], chargaff[,j], pch=20, xlim=xlim, ylim=ylim,
      xlab="",ylab="", xaxt="n", yaxt="n")
      iname <- names(chargaff)[i]
      jname <- names(chargaff)[j]
      direct <- function() segments(0,0,50,50, col="blue")
      invers <- function() segments(0,50,50,0, col="blue")
      PR2 <- function()
      {
        if( iname == "[A]" & jname == "[T]" ) { direct(); return() }
        if( iname == "[T]" & jname == "[A]" ) { direct(); return() }
        if( iname == "[C]" & jname == "[G]" ) { direct(); return() }
        if( iname == "[G]" & jname == "[C]" ) { direct(); return() }
        invers()
      }
      PR2()
      polygon(x=c(0,100,100), y=c(100,100,0), col="lightpink")
      polygon(x=c(0,0,100), y=c(0,100,0))
    }
  }
}
# Clean up
par(op)
rm("PR2")
rm("direct")
rm("i")
rm("iname")
rm ("invers")
rm("j")       
rm("jname")
rm("k")        
rm("lseg")
rm("op")
rm("string")  
rm("xlim")
rm("ylim")
}
\keyword{datasets}
