splitseq<-function(sq,phase=0,nombre=3){
	if(nombre==1) stop("invalid number")
	l<-(floor((length(sq)-phase)/nombre)*nombre)+phase
	sq1=sq[seq(phase+1,l,nombre)]
	for(i in 2:nombre) sq1=paste(sq1,sq[seq(phase+i,l,nombre)],sep="")
	sq1
}