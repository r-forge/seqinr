uco<-function( seq, frame = 0, freq = FALSE, as.data.frame = FALSE){
	sq<-splitseq(seq,frame)
	eff<-table(factor(sq,levels=.SEQINR.UTIL$CODON.AA$CODON))
	if(freq == TRUE) eff<-round(eff/(floor(length(seq)/3)),4)
	if(as.data.frame == TRUE){
	eff=as.data.frame(cbind(as.character(.SEQINR.UTIL$CODON.AA$AA),as.character(.SEQINR.UTIL$CODON.AA$CODON),as.vector(eff)))
	if(freq == TRUE) names(eff)=c("aa","codon","freq")
	else names(eff)=c("aa","codon","eff")
	}
	eff
} 

dotchart.uco <- function(x, numcode = 1, aa3 = TRUE, cex = 0.7, 
  alphabet = s2c("TCAG"), pch = 21, gpch = 20, bg = par("bg"), 
  color = par("fg"), gcolor = par("fg"), lcolor = "gray", 
  xlim, main = NULL, xlab = NULL, ylab = NULL, ...)
{
  if( is.null(names(x)) ) names(x) <- words( alphabet = alphabet )
#
# General sorting 
#
  x <- sort(x)
  labels <- names(x)
  stringlabel = paste(labels, sep="", collapse="")
  groups <- as.factor(translate(s2c(stringlabel), numcode =  numcode))
  gdata <- sapply(split(x, groups), sum)
#
# Now, sorting by aa order
#
  gordered <- rank(gdata)
  xidx <- numeric(64)

  for( i in 1:64 )
  {
    xidx[i] <- -0.01*i + gordered[groups[i]]
  }

  x <- x[order(xidx)]
  labels <- names(x)
  stringlabel = paste(labels, sep="", collapse="")
  aa <- translate(s2c(stringlabel), numcode =  numcode)
  groups <- factor(aa, levels = unique(aa))
  gdata <- sapply(split(x, groups), sum)

  if( missing(xlim) ) xlim <- c(0, max(gdata))
  if( aa3 )
  {
    levels(groups) <- aaa(levels(groups))
  }
  dotchart(x = x, labels = labels, groups = groups, gdata = gdata,
   cex = cex, pch = pch, gpch = gpch, bg = bg, color = color,
   gcolor = gcolor, lcolor = lcolor, xlim = xlim, main = main, 
   xlab = xlab, ylab = ylab, ...)
#
# Return invisibly for further plots
#
  result <- list(0)
  result$x <- x
  result$labels <- labels
  result$groups <- groups
  result$gdata <- gdata

  ypg <- numeric( length(levels(groups)) )
  i <- 1
  for( aa in levels(groups) )
  {
    ypg[i] <- length(which(groups == aa)) + 2
    i <- i + 1
  }
  ypg <- rev(cumsum(rev(ypg))) - 1
  names(ypg) <- levels(groups)
  result$ypg <- ypg

  ypi <- numeric( length(x) )
  for( i in 1:length(x) )
  {
    ypi[i] <- ypg[groups[i]]
  }
  antirank <- function(x) 
  {
    return( seq(length(x),1,by=-1 ))
  }
  ypi <- ypi - unlist(sapply(split(x, groups),antirank))
  names(ypi) <- labels
  result$ypi <- ypi

  return( invisible(result) ) 
}
