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


getSequence.default = function(x){
	if(length(x) == 1) x=s2c(x)	
 	xx = tolower(x)
 	if(length(grep("[acgtu]",xx)) != length(xx)) stop("Biological sequence is needed !")
 	if(y>length(xx) || z>length(xx) || y>z) stop("borns are not correct")	
 	else return(xx)		
}

getFrag.default = function(x,y,z){ 
	if(length(x) == 1) x=s2c(x)	
 	xx = tolower(x)
 	if(length(grep("[acgtu]",xx)) != length(xx)) stop("Biological sequence is needed !")
 	if(y>length(xx) || z>length(xx) || y>z) stop("borns are not correct")	
 	else return(xx[y:z])		
}

getLength.default = function(x){
	if(length(x) == 1) x=s2c(x)	
 	xx = tolower(x)
 	if(length(grep("[acgtu]",xx)) != length(xx)) stop("Biological sequence is needed !")
 	return(length(xx))
}

getName.default = function(x){
 	stop("no name")
}

getProp.default = function(x){
 	stop("no property")
}

getAnnot.default = function(x,y){ 
 	stop("no annotation for this sequence")
}


##################################################################

getFrag =  function(x,y,z) {
	if(! inherits(x,c("SeqFastadna","SeqFastaAA","SeqAcnucLocal","SeqAcnucWeb","SeqFrag"))) { getFrag.default(x,y,z) }
	else UseMethod("getFrag")
}

getSequence = function(x){
	if(! inherits(x,c("SeqFastadna","SeqFastaAA","SeqAcnucLocal","SeqAcnucWeb","SeqFrag"))) {getSequence.default(x)}
	else UseMethod("getSequence")
}


getLength =  function(x) {
	if(! inherits(x,c("SeqFastadna","SeqFastaAA","SeqAcnucLocal","SeqAcnucWeb","SeqFrag"))) {getLength.default(x)}
	else UseMethod("getLength")
}

getName =  function(x) {
	if(! inherits(x,c("SeqFastadna","SeqFastaAA","SeqAcnucLocal","SeqAcnucWeb","SeqFrag"))) {getName.default(x)}
	else UseMethod("getName")
}

getProp =  function(x) {
	if(! inherits(x,c("SeqFastadna","SeqFastaAA","SeqAcnucLocal","SeqAcnucWeb","SeqFrag"))) {getProp.default(x)}
	else UseMethod("getProp")
}


getAnnot = function(x,y) {
	if(! inherits(x,c("SeqFastadna","SeqFastaAA","SeqAcnucLocal","SeqAcnucWeb","SeqFrag"))) {getAnnot.default(x,y)}
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

getSequence.SeqFastadna = function(SeqFastadna){
	return(SeqFastadna)
	}

getFrag.SeqFastadna = function( SeqFastadna, begin, end){
	if(end > getLength(SeqFastadna)) stop("invalid end")	
	newSeq = SeqFastadna[begin:end]
	newSeq = as.SeqFrag(newSeq,begin,end,compl=T,name=getName(SeqFastadna))
	return(newSeq)
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

getAnnot.SeqFastadna = function(Seqfastadna){
	return(attr(SeqFastadna,"Annot"))
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

getSequence.SeqFastaAA = function(SeqFastAA){
	return(SeqFastaAA)
	}


getFrag.SeqFastaAA = function( SeqFastaAA, begin, end){
	if(end > getLength(SeqFastaAA)) stop("invalid end")	
	newSeq = SeqFastaAA[begin:end]
	newSeq = as.SeqFrag(newSeq,begin,end,compl=T,name=getName(SeqFastaAA))
	return(newSeq)
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

getAnnot.SeqFastadna = function(Seqfastadna){
	return(attr(SeqFastaAA,"Annot"))
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


getFrag.SeqAcnucLocal = function(SeqAcnucLocal,born1,born2){
	b = getLength(SeqAcnucLocal)
	if((born2 > b) || (born1 > b)) stop("born out of limits")
	else{  
	s = .Call("getseq2",SeqAcnucLocal,born1,(born2-born1+1))
	seq = s2c(s)
	return(as.SeqFrag(seq,born1,born2,compl=T,name=getName(SeqAcnucLocal)))
	}
}


getSequence.SeqAcnucLocal = function(SeqAcnucLocal){
	return(getseq(SeqAcnucLocal,as.string=F))
	}

getName.SeqAcnucLocal = function(SeqAcnucLocal){
	return(as.character(SeqAcnucLocal))
}

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
	getsequence.socket(attr(SeqAcnucWeb,"socket"),SeqAcnucWeb,start=1,length=b)
}



getFrag.SeqAcnucWeb = function( SeqAcnucWeb ,born1,born2 ){

	b = getLength(SeqAcnucWeb)
	if((born2 > b) || (born1 > b)) stop("born out of limits")  
	bb=born2-born1+1
	newSeq = getsequence.socket(attr(SeqAcnucWeb,"socket"),SeqAcnucWeb,start=born1,length=bb)
	newSeq = as.SeqFrag(newSeq,begin=born1,end=born2,compl=T,name=getName(SeqAcnucWeb))
	return(newSeq)
}



getName.SeqAcnucWeb = function( SeqAcnucWeb ){	

	return( as.character(SeqAcnucWeb) )

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



	############################################################################
	#		Classe de sequences SeqFrag et ses méthodes:               #
	############################################################################




as.SeqFrag = function(seq,begin,end,compl=F,name="frag"){
	if(compl){ attr(seq,"seqMother") = name }
	else attr(seq,"seqMother") = getName(seq)
        attr(seq,"begin") = begin
	attr(seq,"end") = end
	class(seq) = "SeqFrag"
        return(seq)
        }

getSequence.SeqFrag = function(seq){
	return(seq)
	}


getFrag.SeqFrag = function(seq,begin,end){
        if((end<begin) || (end>getLength(seq)))  stop("invalid end")
        newBegin = attr(seq,"begin")+begin-1
        newEnd = attr(seq,"begin")+end-1
	newSeq = seq[begin:end]
        newSeq = as.SeqFrag(seq=newSeq,begin=newBegin,end=newEnd,compl=T,name=getName(seq))
	return(newSeq)
        }

getLength.SeqFrag = function(seq){
	return(attr(seq,"end")-(attr(seq,"begin")+1))
}

getName.SeqFrag = function(seq){
	return(attr(seq,"seqMother"))
}


getprop.SeqFrag = function(seq){
	return(list())
}






	