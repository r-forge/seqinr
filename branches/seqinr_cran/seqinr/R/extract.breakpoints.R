extract.breakpoints <- function(rearr.ori,type=c("atfw","atrev","gcfw","gcrev"),nbreaks,gridsize=100,it.max=500){

  if (!require(segmented)){
    stop("This functions requires the segmented package, but it couldn't be loaded.")
  }
  
  if(length(type)==0){
    stop("You must specify the type of skew: atfw, atrev, gcfw or gcrev. See ?extract.breakpoints for more details.")
  }
  if(sum(!type%in%c("atfw","atrev","gcfw","gcrev"))!=0){
    stop("The type of skew must be one of the following: atfw, atrev, gcfw or gcrev. See ?extract.breakpoints for more details.")
  }

  result=list()
  
  for(t in type){

    if(type=="gcfw"){
      x.breaks=rearr.ori$meancoord.rear[rearr.ori$strand.rear=="forward"]
      y.breaks=cumsum(rearr.ori$gcskew.rear[rearr.ori$strand.rear=="forward"])
     }
     if(type=="gcrev"){
      x.breaks=rearr.ori$meancoord.rear[rearr.ori$strand.rear=="reverse"]
      y.breaks=cumsum(rearr.ori$gcskew.rear[rearr.ori$strand.rear=="reverse"])
     }
    if(type=="atfw"){
      x.breaks=rearr.ori$meancoord.rear[rearr.ori$strand.rear=="forward"]
      y.breaks=cumsum(rearr.ori$atskew.rear[rearr.ori$strand.rear=="forward"])
     }
     if(type=="atrev"){
      x.breaks=rearr.ori$meancoord.rear[rearr.ori$strand.rear=="reverse"]
      y.breaks=cumsum(rearr.ori$atskew.rear[rearr.ori$strand.rear=="reverse"])
     }
    
    assign("x.breaks",x.breaks,envir=globalenv())
    assign("y.breaks",y.breaks,envir=globalenv())
    
    rss=numeric(0)
    starts=list()

    i=0
    while(i<gridsize){
      initpsi=runif(nbreaks[which(type==t)],min=min(x.breaks),max=max(x.breaks))
      if(exists("seg")){
        rm(seg)
      }
      try(seg <- segmented(lm(y~x),x,psi=initpsi,it.max=it.max),silent=T)
      
      if(exists("seg")){
        starts[[length(starts)+1]]=initpsi
        rss=c(rss, sum(seg$residuals^2))
        i=i+1
      }
      
      gc()
    }

    wmin=which.min(rss)
    
    starts=starts[[wmin]]
    
    seg=seg <- segmented(lm(y~x),x,psi=starts,it.max=it.max)
    
    breaks=round(seg$psi[,2])
    
    breaks=breaks[order(breaks)]

    sl=c(x[1],breaks,x[length(x)])-x[1]+1
    
    slopes=numeric(length(sl)-1)
    
    for(i in 1:(length(sl)-1)){
      u=(sl[i]:sl[i+1])+x[1]-1
      slopes[i]=summary(lm(y[sl[i]:sl[i+1]]~u))$coefficients[2,1]
    }
    
    slopes.left=slopes[-length(slopes)]
    slopes.right=slopes[-1]
    
    result[[t]]=list()

    result[[t]]$breaks=breaks
    result[[t]]$slopes.left=slopes.left
    result[[t]]$slopes.right=slopes.right

    real.coord=rearr.ori$meancoord.real[breaks]
    result[[t]]$real.coord=real.coord
    
    
  }

  return(result)

}
