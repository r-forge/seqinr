###################################################################################################
# Ces fonctions permettent d'interroger les banques structurées sous ACNUC au travers des sockets.#
###################################################################################################






choosebank.socket = function( bank ){
	
	# ouverture d'un "client socket" sur le serveur pbil et sur le port 5557 (non définitif :  plus tard 5558)
 
    	socket = socketConnection( host = "pbil.univ-lyon1.fr", port = 5557, server = F, blocking=T)
	rep1 = readLines(socket, n=1)

	# vérification du bon fonctionnement de la connection et ouverture de la banque
    
	request = paste("acnucopen&db=",bank,sep="")

	writeLines( request, socket, sep = "\n")
	rep2 = readLines(socket, n=1)
	
	if(parser.socket(rep2) == "1"){
		print("bank name incorrect")
		rm(socket)
	}
	else{ 
		assign("banknameSocket",bank,.GlobalEnv)
		return(list(socket=socket, bankname = bank))
	}
}




closebank.socket = function( socket ){
	
	# fermeture de la banque acnuc

	writeLines( "acnucclose", socket, sep="\n" )
	rep = readLines( socket ,n=1)
	if(parser.socket(rep)=="0") print("OK ACNUC SOCKET STOPED")
	close(socket)
	rm(socket)
}


parser.socket = function(p)
{
	p1=s2c(p)
	b=c(which(p1=="="))
  	a=c(which(p1=="&"),length(p1)+1)
	o=character(length(a))
	for(i in 1:length(a)) {o[i]=substr(p,b[i]+1,a[i]-1) }
	return(o)
}
	


getsequence.socket = function( socket, name, start, length){
	
	request2 = paste("gfrag&name=", name,"&start=", start, "&length=", length, sep= "")
	writeLines( request2, socket, sep="\n" )
	s = readLines(socket,n=1)

	if(length(s)==0){
		 print("invalid sequence name")
		}
	else{		
		s = s2c(s)
		sequence = s[(grep("&",s)+1):length(s)]
		return(sequence)
	}
}
	


getAttribut.socket = function( socket, name){
	
	# Récupération des attributs d'une séquence
	
	request = paste( "isenum&name=",name, sep="")
	writeLines( request, socket, sep="\n")
	res = readLines( socket, n=1 )
	p=parser.socket(res)
	l=list( length = as.numeric(p[2]),frame = as.numeric(p[3]), gencode =as.numeric(p[4]) )
	return(l)
}


readAnnots.socket = function( socket, name, nl){

	request= paste( "read_annots&name=", name, "&nl=", nl, sep= "")
	annots=character(nl)
	for(i in 1:nl){
	   writeLines( request , socket, sep="\n")
	   annots[i] =  readLines( socket , n=1 )
	}
	return(annots)
}




getNumber.socket = function( socket, name){
	
	request = paste( "isenum&name=",name, sep="")
	writeLines( request, socket, sep="\n")
	s = readLines(socket,n=1)
	return(parser.socket(s)[1])
}



query.socket = function( socket, listname , query ){
	
	writeLines( "prep_requete", socket , sep="\n")
	readLines( socket, n=1 )
	request = paste( "proc_requete&query=\"", query, "\"&name=\"", listname,"\"", sep="")
	writeLines( request, socket, sep= "\n")
	res = readLines( socket , n=1 )
	p=parser.socket(res)

	# initialisation des variables	

	lrank=p[2]
	first=1
	liste = character(as.numeric(p[3]))
	
	# récupération de tous les éléments de la liste
	
	for(i in 1:length(liste)){
		writeLines(paste("nexteltinlist&lrank=",lrank,"&first=",first,"&type=SQ",sep=""),socket, sep= "\n")
		res = readLines( socket , n=1 )
		r = parser.socket(res)
		first = r[1]
		liste[i] = r[2] 
	}

	liste=lapply(liste,as.SeqAcnucWeb,socket)
	result = list(call = match.call(), name= listname, req = as.list(liste), socket = socket)
	class(result) = c("qaw")
	assign(listname, result, env = .GlobalEnv)
	print(result)
}



print.qaw <- function(x, ...)
{

	cat("\n")
	cat("\n$socket: ")
	print(x$socket)
	cat("\n$banque: ")
	cat(get("banknameSocket",env=.GlobalEnv))
	cat("\n$call: ")
	print(x$call)
	cat("$name: ")
	print(x$name)
	cat("\n")
	sumry <- array("", c(1, 4), list(1, c("list", "length", "mode", "content")))
	sumry[1, ] <- c("$req",length(x$req),"character","sequences")
	class(sumry) <- "table"
	print(sumry)
	cat("\n")
}


	



getKeyword.socket = function( socket, name )
{

}

