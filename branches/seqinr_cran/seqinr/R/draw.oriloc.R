draw.oriloc <- function(ori){

  ta=ori$x
  cg=ori$y
  skew=ori$skew
  cdsskew=ori$CDS.excess

  ## wrapped cds

  start.kb=ori$start.kb
  end.kb=ori$end.kb
  

  wrapped=which(abs(ori$start.kb-ori$end.kb)>=50) 
  if(wrapped==length(cdsskew)){
    start.kb[wrapped]=max(start.kb[wrapped],end.kb[wrapped])
    end.kb[wrapped]=max(start.kb[wrapped],end.kb[wrapped])
  }
  if(wrapped==1){
    end.kb[wrapped]=min(start.kb[wrapped],end.kb[wrapped])
    start.kb[wrapped]=min(start.kb[wrapped],end.kb[wrapped])
  }
     
  meancoord=(start.kb+end.kb)/2

  ymin <- min(ta,cg,skew)
  ymax <- max(ta,cg,skew)
  xmin <- min(meancoord)
  xmax <- max(meancoord)
  
  ticks <- pretty(cdsskew)
  
  ticks.y <- (ymax-ymin)/(max(cdsskew)-min(cdsskew))*(ticks - min(cdsskew)) + ymin
  cds.y   <- (ymax-ymin)/(max(cdsskew)-min(cdsskew))*(cdsskew - min(cdsskew)) + ymin

 plot(meancoord, cg, type="l", xlab="Map position (gene index)",
      ylab = "Cumulated normalized skew",xlim=c(xmin,xmax),ylim=c(ymin,ymax),cex.lab=1.35,col="lightblue",main="Rearranged nucleotide skews",lwd=1.5)

  axis(side = 4, at = ticks.y, labels = ticks, col = "lightgreen", col.axis ="black")

  tmp <- pretty(meancoord)
  abline(v=tmp, col="grey", lty=3)
  tmp <- tmp[-length(tmp)] + diff(tmp)/2
  abline(v=tmp, col="grey", lty=3)
  
  lines(meancoord,ta, col="pink",lwd=1.5)

  lines(meancoord,skew, col="black",lwd=1.5)
  
  lines(meancoord,cds.y, col="lightgreen",lwd=1.5)
  
  mtext("Cumul. T-A skew", col="pink", adj=0)
  mtext("Cumul. C-G skew", col="lightblue")
  mtext("Cumul. CDS skew", col="lightgreen", adj=1)

 

}
