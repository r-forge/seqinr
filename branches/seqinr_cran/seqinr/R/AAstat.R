AAstat = function(seq, plot = TRUE){

	tutu=list(Tiny = which(object %in% SEQINR.UTIL$AA.PROPERTY$Tiny),
	Small = which(object %in% SEQINR.UTIL$AA.PROPERTY$Small),
	Aliphatic = which(object %in% SEQINR.UTIL$AA.PROPERTY$Aliphatic),
	Aromatic = which(object %in% SEQINR.UTIL$AA.PROPERTY$Aromatic),
	Non.polar = which(object %in% SEQINR.UTIL$AA.PROPERTY$Non.polar),
	Polar = which(object %in% SEQINR.UTIL$AA.PROPERTY$Polar),
	Charged = which(object %in% SEQINR.UTIL$AA.PROPERTY$Charged),
	Basic = which(object %in% SEQINR.UTIL$AA.PROPERTY$Basic),
	Acidic = which(object %in% SEQINR.UTIL$AA.PROPERTY$Acidic))
	
	if(plot == TRUE){
	coul = rainbow(length(tutu))
	plot(c(0,length(object)),c(0,length(tutu)+1),type="n",axes=FALSE,ann=FALSE,xlim=c(0,length(object)+1))
	title(xlab=paste("Position of the residues along the sequence",sep=""))
	axis(2, at = seq(1.5,10,1), lab = names(tutu), col.lab = "blue",las=1)
	axis(1,at = seq(0,length(object),15), lab = seq(0,length(object),15) , col.axis = "blue")
	m=lapply(1:length(tutu),function(x){
	segments(tutu[[x]],x,tutu[[x]],x+1,col=coul[x],lwd=2)
	rect(0,x,length(object),1,lwd=2)
	 }	
	)
	rect(0,length(tutu)+1,length(object),1,lwd=2)
 }
	res1 = lapply(tutu,function(x){ length(x)/length(object) })
	res2 = table(factor(object, levels = levels(SEQINR.UTIL$CODON.AA$L)))
	res3 = computePI(object)
	return(list(res2,res1,PI=res3))
}
