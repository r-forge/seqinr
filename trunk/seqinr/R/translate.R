translate<-function(sq,phase=0,sens="+",numcode=1)
{
	if(sens=="-") sq<-comp(invers(sq))
	else if(sens=="+") 
	sq<-replace(sq,sq=="A",2)
	sq<-replace(sq,sq=="C",1)
	sq<-replace(sq,sq=="G",3)
	sq<-replace(sq,sq=="T",0)
	sq<-as.numeric(sq) 
	l<-floor((length(sq)-phase)/3)*3
	a<-seq(phase+1,l+phase,3)
	b<-seq(phase+2,l+phase,3)
	c<-seq(phase+3,l+phase,3)
	tra<-16*sq[a]+4*sq[b]+sq[c]+1
	code<-unlist(strsplit(as.vector(CODES[numcode,]$V2),""))	
	code[tra]
}


