######################################################################
#classes de séquences
# toutes les classes doivent aoir exactement la même interface à savoir:
# une fonction initNomClasse qui retourne une instance de la classe
# des spécialisations des fonctions
#                    getFrag(seq,begin=0,end=getLength(seq)) retourne un tableau de caractères
#                    getLength(seq) retourne un "entier"
#                    getName(seq) retourne une chaîne
#                    getProp(seq) retourne une liste nommée
######################################################################



getFragC=function(x,y,z){return(NULL)}
getLengthC=function(x){return(0)}
getNameC=function(x){return(NULL)}
getPropC=function(x){return(list())}

getFrag =  function(x,y,z) {
if(is.null(attr(x,"class"))) {getFragC(x,y,z)}
else UseMethod("getFrag")
}

getLength =  function(x,y,z) {
if(is.null(attr(x,"class"))) {getLengthC(x,y,z)}
else UseMethod("getLength")
}

getName =  function(x,y,z) {
if(is.null(attr(x,"class"))) {getNameC(x,y,z)}
else UseMethod("getName")
}

getProp =  function(x,y,z) {
if(is.null(attr(x,"class"))) {getPropC(x,y,z)}
else UseMethod("getProp")
}



######################################
#   Classe de sequence SeqSimple
####################################

initSeqSimple = function(s,name="seq") {
        r=list(seq=s,name=name)
        class(r)="SeqSimple"
        return(r)
        }

getFrag.SeqSimple=function(SeqSimple,begin=1,end=getLength(seq)){
        if(end<begin)return(NULL)
        return(s2c(substr(SeqSimple$seq,begin,end)))
        }
        
getLength.SeqSimple = function(SeqSimple){return(nchar(SeqSimple$seq))}

getName.SeqSimple = function(SeqSimple){return(SeqSimple$name)}

getprop.SeqSimple = function(SeqSimple){return(list())}


################################################################################
# Ce qui touche aux sous-sequences
#  		Classe SeqFrag : un fragment d'une séquence
#BEWARE seq must be a character string representing the name of the variable
################################################################################

initSeqFrag=function(seq,begin,end,name="Frag",compl=FALSE,num=0){
        if(compl){
                nom=paste(name,num,sep="")
                }
        else nom=name
        seq=as.symbol(seq)
        r=list(seqmother=seq,begin=begin,end=end,name=nom)
        class(r)="SeqFrag"
        return(r)
        }

getFrag.SeqFrag=function(seq,begin=1,end=getLength(seq)) {
        if(end<begin)return(NULL)
        s=eval(as.symbol(seq$seqmother),envir=globalenv())
        bm=seq$begin-1
        b=bm+begin
        em=bm+end
        e=min(em,s$end)
        return(getFrag(s,b,e))
        }

getLength.SeqFrag = function(seq){return(seq$end-seq$begin+1)}

getName.SeqFrag = function(seq){return(seq$name)}

getprop.SeqFrag = function(seq){return(list())}


##############################################
#
#	Class SeqDaugther : la concaténation de plusieurs séquences (par exemple des fragments)
#
##############################################



seqDC=function(seq,i){ # MUST NEVER BE CALLED
        if(class(seq)!="SeqDaughter") return(NULL)
        if(i<=seq$nbFrag)return(seq[1][[1]][i])
        else return(NULL)
        }

initSeqDaughter=function(name,...) {
        l=list(...)
        l=c(l,name=name)
        nbFrag=length(l[1][[1]])
        l=c(l,nbFrag=nbFrag)
        class(l)="SeqDaugther"
        length=0
        for(i in 1:l$nbFrag) length=length+getLength(SeqDC(l,i))
        l=c(l,length=length)
        return(l)
        }

getLength.SeqDaugther = function(seq){return(seq$length)}

getName.seqDaugther = function(seq){return(seq$name)}

getFrag.seqDaugther = function(seq,begin=1,end=getLength(seq)){
        if(end<begin)return(NULL)
        if(begin>getLength(seq))return(NULL)
        if(end>getLength(seq))end=getLength(seq)
        if(begin<1)begin=1
        k=0
# looking for first implied fragment
        i=1
        for(i in 1:seq$nbFrag) {
                k=k+getLength(seqDC(seq,i))
                if(begin<k)break
                }
        b=getLength(seqDC(seq,i))-k+begin #the begin in the current frag
        e=getLength(seqDC(seq,i))-k+end
        e1=min(e,getLength(seqDC(seq,i))) #the end in the current frag
        frag=getFrag(seqDC(seq,i),b,e1)
        e=e-e1
        while(e>0){
                i=i+1
                e1=min(e,getLength(seqDC(seq,i)))
                frag=c(frag,getFrag(seqDC(seq,i),1,e1))
                e=e-e1
                }
        return(frag)
        }


##############################################################################################
#		Classe de sequence SeqFasta et ses méthodes:
#La classe de séquence SeqFasta pour les séquences résultants de la lecture d'un fichier au
#format fasta. 
###########################################################################"

##################################
#initSeqFasta sera appelée au moment de la lecture d'un fichier au format fasta par read.fasta()
####################################

initSeqFasta = function(elemlist){
	class(elemlist)="SeqFasta"	
        return(elemlist)
        }

getFrag.SeqFasta = function( SeqFasta, begin = 1, end = getLength(SeqFasta)){
	if(end > getLength(SeqFasta)) stop("invalid end")	
	return(SeqFasta[[1]][begin:end])
	}

getLength.SeqFasta = function( SeqFasta){
	return(length(SeqFasta[[1]]))
	}

getName.SeqFasta = function( SeqFasta){
	return(names(SeqFasta))
}


#############################################################################################
#		Classe de sequence SeqReq et ses méthodes: 
#La classe de séquence SeqReq pour les séquences provenant 
#d' une requete dans les banques structurées sous ACNUC
#Cette classe ne contiendra pas la sequence car le but est de ne pas la stocker.
#On stockera le mnemo et le nom de la banque d'où elle provient.
#
#############################################################################################


######
# l'initialisation se fera au moment de l'appel à la fonction getreq. On attribuera à chaque mnemo
# de la liste résultant de la requete la classe SeqReq 
#EN COURS: INCORPORATION DU NOM DE LA BANQUE D'0RIGINE
######

initSeqReq = function(name)
	class(name)="SeqReq"
	return(name)
}

getFrag.SeqReq = function(SeqReq, born1 , born2 ){ 
	s = .Call("getseq2",SeqReq,born1,born2)
	return(toupper(s2c(s)))
}

getName.SeqReq = function(SeqReq){return(SeqReq)}

getLength.SeqReq = function(SeqReq){
	s = .Call("getseq",SeqReq)
	return(nchar(s))
	}




