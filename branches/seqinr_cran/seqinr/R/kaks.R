kaks = function(x){
	if(attr(x,"class") != "alignment") error("object x must be of class alignment")
	x$nb = as.numeric(x$nb)
	l = .Call("kaks",x$seq,x$nb,PACKAGE="seqinr")
	m = lapply(l,function(k){
	if(! is.null(x$nam)) tmp = matrix( k, x$nb, x$nb, byrow = TRUE, dimnames=list(x$nam,x$nam))
	else{
	 n = paste("seq",c(1:x$nb),sep="")
	 tmp = matrix( k, x$nb, x$nb, byrow = TRUE, dimnames=list(n,n))
	}
	 as.dist(t(tmp))
	})
	m = lapply(m,round,digits=4)
	names(m)=c("ks","ka","vks","vka")
	return(m)
}
