read.alignment = function( File, format)
{

	ali = switch( format,
	fasta = .Call("read_fasta_align",File), 
	mase = .Call("read_mase",File),
	phylip = .Call("read_phylip_align",File),
	msf = .Call("read_msf_align",File),
	clustal = .Call("read_clustal_align",File))
	ali = lapply(ali, as.character)
	if(format == "mase") ali = list(nb = ali[[1]], nam = ali[[2]], seq = ali[[3]], com =ali[[4]]) 
	else ali = list(nb = ali[[1]], nam = ali[[2]], seq = ali[[3]], com = NA)
	class(ali)="alignment"
	return(ali)

}


dist.alignment = function( x , matrix = "similarity" )
{
	if (!inherits(x, "alignment")) 
        stop("Object of class 'alignment' expected")
	t1 = c("similarity","identity")
	m1 = grep(matrix,t1)
	m2 = as.numeric(.Call("is_a_protein_seq",x$seq[1]) >= 0.8)
	l = as.numeric(x$nb)
	maxseq = length(x$seq)
	dist = .Call("distance",x$seq,l,m1,m2)
	mat = matrix(dist,l,l, byrow =TRUE)
        dimnames(mat) = list(x$nam,x$nam)
	return( as.dist(mat) )

}
