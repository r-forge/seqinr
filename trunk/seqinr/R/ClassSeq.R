	#####################################################################################
	# 		classes de séquences						    #
	#  toutes les classes doivent avoir exactement la même interface à savoir:	    #
	#  une fonction initNomClasse qui retourne une instance de la classe		    #
	#  des spécialisations des fonctions                                                #
	#	     getSequence(seq) retourne un vecteur de char	     		    #
	#            getFrag(seq,begin,end) retourne un vecteur de caractères               #
	#            getLength(seq) retourne un "entier"                                    #
	#            getName(seq) retourne une chaîne                                       #
	#            getProp(seq) retourne une liste nommée                                 # 
	#####################################################################################


getSequenceC=function(x){return(NULL)}
getFragC=function(x,y,z){return(NULL)}
getLengthC=function(x){return(0)}
getNameC=function(x){return(NULL)}
getPropC=function(x){return(list())}
getAnnotC=function(x,y){return(NULL)}

getFrag =  function(x,y,z) {
if(is.null(attr(x,"class"))) {getFragC(x,y,z)}
else UseMethod("getFrag")
}

getSequence = function(x){
if(is.null(attr(x,"class"))) {getSequenceC(x)}
else UseMethod("getSequence")
}


getLength =  function(x) {
if(is.null(attr(x,"class"))) {getLengthC(x,y,z)}
else UseMethod("getLength")
}

getName =  function(x) {
if(is.null(attr(x,"class"))) {getNameC(x,y,z)}
else UseMethod("getName")
}

getProp =  function(x) {
if(is.null(attr(x,"class"))) {getPropC(x,y,z)}
else UseMethod("getProp")
}


getAnnot = function(x,y) {
if(is.null(attr(x,"class"))) {getAnnotC(x,y,z)}
else UseMethod("getAnnot")
}






	########################################################################################################
	#		Classe de sequence SeqFastadna et ses méthodes:                                        #
	#	La classe de séquence SeqFasta pour les séquences résultants de la lecture d'un fichier au     # 
	#	format fasta.                                                                                  #
	########################################################################################################

	##################################################################################################
	# as.SeqFasta sera appelée au moment de la lecture d'un fichier au format fasta par read.fasta() #
	##################################################################################################

as.SeqFastadna = function(elemlist){
	class(elemlist)="SeqFastadna"	
        return(elemlist)
        }

is.SeqFastadna = function(x){
	
	inherits(x,"SeqFastadna")
}


getFrag.SeqFastadna = function( SeqFastadna, begin = 1, end = getLength(SeqFastadna)){
	if(end > getLength(SeqFastadna)) stop("invalid end")	
	return(SeqFastadna[begin:end])
	}

getLength.SeqFastadna = function( SeqFastadna){
	return(length(SeqFastadna))
	}

getName.SeqFastadna = function( SeqFastadna){
	return(attr(SeqFastadna,"name"))
}

getProp.SeqFastadna = function(Seqfastadna){
	return(list(seqtype="DNA"))
}

summary.SeqFastadna = function(SeqFastadna){
	compo=count(SeqFastadna,1)
	return(list(composition=compo,GC=GC(SeqFastadna)))
}


	###############################################################################
	#		Classe de sequences SeqFastaAA et ses méthodes:               #
	###############################################################################

as.SeqFastaAA = function(elemlist){
	class(elemlist)="SeqFastaAA"	
        return(elemlist)
        }

is.SeqFastaAA = function(x){
	
	inherits(x,"SeqFastaAA")
}


getFrag.SeqFastaAA = function( SeqFastaAA, begin = 1, end = getLength(SeqFastaAA)){
	if(end > getLength(SeqFastaAA)) stop("invalid end")	
	return(SeqFastaAA[begin:end])
	}


getLength.SeqFastaAA = function( SeqFastaAA){
	return(length(SeqFastaAA))
	}


getName.SeqFastaAA = function( SeqFastaAA){
	return(attr(SeqFastaAA,"name"))
}


getProp.SeqFastaAA = function(SeqfastaAA){
	return(list(seqtype="AA"))
}


summary.SeqFastaAA = function(SeqFastaAA){
	compo=table(factor(SeqFastaAA, levels = levels(SEQINR.UTIL$CODON.AA$L)))
	return(list(composition=compo/getLength(SeqFastaAA),AA.Property=AApropr(SeqFastaAA)))
}




	############################################################################################
	#	
	#		Classe de sequences SeqAcnucLocal et ses méthodes: 			   #
	# La classe de séquence SeqAcnucLocal pour les séquences provenant                         #
	# d' une requete dans les banques structurées sous ACNUC                                   #
	# Cette classe ne contiendra pas la sequence car le but est de ne pas la stocker.          #
	# On stockera le mnemo et le nom de la banque d'où elle provient.                          #
	#                                                                                          #
	############################################################################################


	###################################################################################################
	# l'initialisation se fera au moment de l'appel à la fonction getreq. On attribuera à chaque mnemo#
	# de la liste résultant de la requete la classe SeqAcnucLocal					  #	
	###################################################################################################



as.SeqAcnucLocal = function(name){
	class(name)="SeqAcnucLocal"
	return(name)
}


is.SeqAcnucLocal = function(x){

	inherits(x, "SeqAcnucLocal")
}


getFrag.SeqAcnucLocal = function(SeqAcnucLocal,born1=1,born2){
	b = getLength(SeqAcnucLocal)
	if((born2 > b) || (born1 > b)) stop("born out of limits")
	else{  
	s = .Call("getseq2",SeqReq,born1,born2)
	return(s2c(s))
	}
}


getSequence.SeqAcnucLocal = function(SeqAcnucLocal){
	return(getseq(SeqReq,as.string=F))
	}



getName.SeqAcnucLocal = function(SeqAcnucLocal){return(SeqAcnucLocal)}

getLength.SeqAcnucLocal = function(SeqAcnucLocal){
	return(getAttribut(SeqAcnucLocal)[[1]])
	}

getProp.SeqAcnucLocal = function(SeqAcnucLocal){
	return(getAttribut(SeqAcnucLocal)[2:3])
}

getAnnot.SeqAcnucLocal = function(SeqAcnucLocal,nbl){
	return(getAnnots(SeqAcnucLocal,nbl))
}


#summary.SeqAcnucLocal = function(SeqAcnucLocal){
# 	return(list(mnemo=getName(SeqAcnucLocal),GC.percent=GC(SeqAcnucLocal),base.count=count(SeqAcnucLocal,1)))
#	}



	######################
	# fonctions annexes  #
	######################



AApropr = function(SeqFastaAA){
	s=table(factor(SeqFastaAA, levels = levels(SEQINR.UTIL$CODON.AA$L)))
	t=sum(s)
	list(Tiny=sum(s[which(names(s) %in% SEQINR.UTIL$AA.PROPERTY$Tiny)])/t,Small=sum(s[which(names(s) %in% SEQINR.UTIL$AA.PROPERTY$Small)])/t,Aliphatic=sum(s[which(names(s) %in% SEQINR.UTIL$AA.PROPERTY$Aliphatic)])/t,Aromatic=sum(s[which(names(s) %in% SEQINR.UTIL$AA.PROPERTY$Aromatic)])/t,Non.polar=sum(s[which(names(s) %in% SEQINR.UTIL$AA.PROPERTY$Non.polar)])/t,Polar=sum(s[which(names(s) %in% SEQINR.UTIL$AA.PROPERTY$Polar)])/t,Charged=sum(s[which(names(s) %in% SEQINR.UTIL$AA.PROPERTY$Charged)])/t,Basic=sum(s[which(names(s) %in% SEQINR.UTIL$AA.PROPERTY$Basic)])/t,Acidic=sum(s[which(names(s) %in% SEQINR.UTIL$AA.PROPERTY$Acidic)])/t)
}






####################################################################################################
#												   #
#	Classe de Sequences SeqAcnucWeb                                                            #
#												   #			
####################################################################################################





as.SeqAcnucWeb = function( name, socket=F ){

	class(name)="SeqAcnucWeb"
	attributes(name)=list(class="SeqAcnucWeb",socket=socket)
	name
}


is.SeqAcnucWeb = function( x ){
	
	inherits(x,"SeqAcnucWeb")

}



getSequence.SeqAcnucWeb = function( SeqAcnucWeb){
	b=getLength( SeqAcnucWeb )
	getfrag( SeqAcnucWeb,1,b)
}



getFrag.SeqAcnucWeb = function( SeqAcnucWeb ,born1,born2 ){

	b = getLength(SeqAcnucWeb)
	if((born2 > b) || (born1 > b)) stop("born out of limits")  
	bb=born2-born1+1
	return( getsequence.socket(attr(SeqAcnucWeb,"socket"),SeqAcnucWeb,start=born1,length=bb))
}



getName.SeqAcnucWeb = function( SeqAcnucWeb ){	

	return( SeqAcnucWeb )

}

getLength.SeqAcnucWeb = function( SeqAcnucWeb ){

	return( getAttribut.socket(attr(SeqAcnucWeb,"socket"),SeqAcnucWeb)[[1]] )

}


getProp.SeqAcnucWeb = function( SeqAcnucWeb ){

	return( getAttribut.socket(attr(SeqAcnucWeb,"socket"),SeqAcnucWeb) ) 

}

getAnnot.SeqAcnucWeb = function( SeqAcnucWeb, nbl ){
		
	return( readAnnots.socket( socket= attr(SeqAcnucWeb,"socket"),name=SeqAcnucWeb, nl=nbl) ) 

}
