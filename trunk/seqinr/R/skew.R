#
# This should be removed after R_1.7
#

require(methods) 

#
# Class representation
#

setClass("skew", representation(
  bp = "numeric",
  Kb = "numeric", 
  Mb = "numeric", 
  ATskew  = "numeric",
  GCskew  = "numeric",
  wsizekb = "numeric",
  wstepkb = "numeric"   )
)

#
# Class constructor
#

skew <- function( seq, wsizekb = 10, wstepkb = 10)
{
  if( missing(seq) ) # demo file
  {
    filename = system.file("sequences/bb.fasta", package = "seqinr")
    seq <- tolower( read.fasta( filename )[[1]]) 
  }

  lseq <- length(seq)
  bp <- seq( from = 1, to = lseq, by = 1000*wstepkb )
  a <- ifelse( seq == "a", 1, 0)
  t <- ifelse( seq == "t", 1, 0)
  g <- ifelse( seq == "g", 1, 0)
  c <- ifelse( seq == "c", 1, 0)
  at <- function(pos)
  {
    subseq <- pos:(pos+1000*wsizekb)
    subseq <- subseq[ subseq <= lseq ]
    na <- sum(a[subseq])
    nt <- sum(t[subseq])
    return( 100*(na - nt)/(nt + na) )
  }
  gc <- function(pos)
  {
    subseq <- pos:(pos+1000*wsizekb)
    subseq <- subseq[ subseq <= lseq ]
    ng <- sum(g[subseq])
    nc <- sum(c[subseq])
    return( 100*(nc - ng)/(nc + ng) )
  }
  ATskew <- sapply(bp, at)
  GCskew <- sapply(bp, gc)
  new("skew", bp = bp, Kb = bp/1000, Mb = bp/1000000, 
    ATskew = ATskew, GCskew = GCskew, wsizekb = wsizekb,
    wstepkb = wstepkb )
}

#
# Plot DNA walk
#
setMethod("plot",
  signature(x = "skew", y = "missing"),
  definition = function(x, y, option = 1, ...)
  {
    fillpoly <- function(x, y)
    {
      for( i in 1:(length(x)-1 ) )
      {
        fill <- FALSE
        if( y[i] >= 0 & y[i+1] >= 0 )
        {
          xpol <- c(x[i],x[i],x[i+1],x[i+1])
          ypol <- c(0, y[i], y[i+1],0)
          fill <- TRUE
        }
        if( y[i] >= 0 & y[i+1] < 0 )
        {
          xinter <- x[i] + y[i]*(x[i+1]-x[i])/(y[i]-y[i+1])
          xpol <- c(x[i], x[i], xinter)
          ypol <- c(0, y[i], 0)
          fill <- TRUE
        }
        if( y[i] < 0 & y[i+1] >= 0 )
        {
          xinter <- x[i+1] - y[i+1]*(x[i+1]-x[i])/(y[i+1]-y[i])
          xpol <- c(xinter, x[i+1], x[i+1])
          ypol <- c(0, y[i+1], 0)
          fill <- TRUE
        }
        if( fill )
        {
          polygon( xpol, ypol, col="black" )
        }
      }
    }

    if( option == 1 ) # GC skew only
    {
      ylab <- paste("GC skew: 100*(C-G)/(C+G) in", x@wsizekb, "Kb window")
      xlab <- "Map position in Kb"
      plot(x@Kb, x@GCskew, ,ylab=ylab, xlab=xlab,type="l", ...)
      abline(h=0)
      fillpoly( x@Kb, x@GCskew )
    }
    if( option == 2 ) # AT skew only
    {
      ylab <- paste("AT skew: 100*(A-T)/(A+T) in", x@wsizekb, "Kb window")
      xlab <- "Map position in Kb"
      plot(x@Kb, x@ATskew, ,ylab=ylab, xlab=xlab,type="l", ...)
      abline(h=0)
      fillpoly( x@Kb, x@ATskew )
    }

  }
)
