scanK = function(elemlist, motif){
	if(is.null(elemlist)) return(NULL)
	else if(any(i = grep(motif,elemlist,value=F,ignore.case=T,perl=F))) elemlist
	}

# prend en paramètre une liste de mnémonique et un motif et renvoie un vecteur de chaine de charactères contenant les  mnemos qui comportent le motif recherché dans leurs keywords.


scanKeyword = function(listmnemo, motif){
	tmp = lapply(listmnemo,getKeyword)
	names(tmp) = listmnemo
	tmp = lapply(tmp,scanK,motif)
	tmp = lapply(tmp,is.null)
	tmp = unlist(tmp)
	tmp = tmp[tmp == F]
	return(names(tmp))
}

#idem scanKeyword mais pred en parametre le nombre de ligne d'annotations sur lesquelles on désire efectuer la recherche de motif

scanAnnot = function(listmnemo, motif,nligne){
	tmp = lapply(listmnemo,getAnnots,nligne)
	names(tmp) = listmnemo
	tmp = lapply(tmp,scanK,motif)
	tmp = lapply(tmp,is.null)
	tmp = unlist(tmp)
	tmp = tmp[tmp == F]
	return(names(tmp))
}




