read.mase = function( File = system.file("sequences/test.mase", package = "seqinr") ){
	
	mase = .Call("read_mase",File)
	mase = lapply(mase,as.character)
	names(mase)=list("seq.number","site.number","sequence","comment","name")
	class(mase)=c("SeqMase")
	return( mase )

	}


dist.mase = function( mase , matrix = "similarity"){
	
	if (!inherits(mase, "SeqMase")) 
        stop("Object of class 'SeqMase' expected")	
	t = c("similarity","identity")
	m = grep(matrix,t)
	l=as.numeric(mase$seq.number)
	dist = .Call("distance",mase$sequence,l,m)
	mat=matrix(dist,l,l, byrow =TRUE)
        dimnames(mat)=list(mase$name,mase$name)
	return(as.dist(mat))

}
