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
  xlim = c( 0, max(x)), main = NULL, xlab = NULL, ylab = NULL, ...)
{
  if( is.null(names(x)) ) names(x) <- words( alphabet = alphabet )
#
# General sorting 
#
  x <- sort(x)
  labels <- names(x)
  groups <- as.factor(translate(s2c(labels), numcode =  numcode))
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
  aa <- translate(s2c(labels), numcode =  numcode)
  groups <- factor(aa, levels = unique(aa))
  gdata <- sapply(split(x, groups), sum)

  xlim <- c(0, max(gdata))
  if( aa3 )
  {
    levels(groups) <- aaa(levels(groups))
  }

  dotchart(x = x, labels = labels, groups = groups, gdata = gdata,
   cex = cex, pch = pch, gpch = gpch, bg = bg, color = color,
   gcolor = gcolor, lcolor = lcolor, xlim = xlim, main = main, 
   xlab = xlab, ylab = ylab, ...) 
}
