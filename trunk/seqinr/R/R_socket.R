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
        res = parser.socket(rep2)

	if(res[1] != "0"){
		print("bank name incorrect")
		rm(socket)
	}
	else{ 
		assign("banknameSocket",bank,.GlobalEnv)
		return(list(socket=socket, bankname = bank, totseqs = res[3], totspecs=res[4], totkeys=res[5]))
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
	


getSequence.socket = function( socket, name, start, length){
	
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

	if(as.numeric(p[1]) != 0)  stop(paste("invalid request:",p[2],sep=""))

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
	liste=lapply(liste,noquote)
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



getKeyword.socket <- function( socket, name){

         writeLines(paste("isenum&name=",name,sep=""),socket,sep="\n")
         res = readLines( socket , n=1 )
         number = parser.socket(res)[1] 

         writeLines(paste("readsub&num=",number,sep=""),socket,sep="\n")
         res2 = readLines( socket , n=1 ) 
         rr = parser.socket(res2)
         
         writeLines(paste("readshrt&num=",rr[7],sep=""),socket,sep="\n")
         res3 = readLines( socket , n=1 ) 
         
	 p1=s2c(res3)
	 b=c(which(p1=="="))
  	 a=c(which(p1=="&"))
         d=c(which(p1==","),length(p1)+1)
	 o=character(length(a))
         o[1]=substr(res3,a[2]+1,d[1]-1)
	 s = seq(2,length(a)*2,by=2)
         for(i in 1:(length(s)-1)){o[i+1] = substr(res3,d[s[i]]+1,d[s[i]+1]-1)} 

	 lapply(o,function(x){
          	writeLines(paste("readkey&num=",x,sep=""),socket,sep="\n")	
	        res4 = readLines( socket , n=1 ) 
                parser.socket(res4)[2]
})

} 


getExon.socket <- function( socket, name){

	 writeLines(paste("isenum&name=",name,sep=""),socket,sep="\n")
         res = readLines( socket , n=1 )
         number = parser.socket(res)[1] 
	 
         writeLines(paste("readsub&num=",number,sep=""),socket,sep="\n")
         res2 = readLines( socket , n=1 ) 
         rr = parser.socket(res2)
  
	 # Test si subsequence           
     
	 l=list()	
         if(as.numeric(rr[5]) != 0) stop(" pas un CDS : sequence mere\n")
         else {
		i=1
 	 	writeLines(paste("readext&num=",rr[6],sep=""),socket,sep="\n")		
		res3 = readLines( socket , n=1 )
		r = parser.socket(res3)
		l[[i]] = c(r[3],r[4])
		n=r[5] 
	}
        while(as.numeric(n) != 0){
		i=i+1
		writeLines(paste("readext&num=",n,sep=""),socket,sep="\n")	 
		res4 = readLines( socket , n=1 )
		rrr = parser.socket(res4)
		l[[i]] = c(rrr[3],rrr[4])
		n=rrr[5]
  		}
	return(l)
}	



#
#getSpecies.socket <- function( socket, name){
#
#         writeLines(paste("isenum&name=",name,sep=""),socket,sep="\n")
#         res = readLines(socket , n=1)
#         number = parser.socket(res)[1] 
#
#         writeLines(paste("readsub&num=",number,sep=""),socket,sep="\n")
#         res2 = readLines( socket , n=1 ) 
#         rr = parser.socket(res2)
#	
#	 if(as.numeric(rr[5]) != 0){ num = rr[5] }
#		else{
#			writeLines(paste("readext&num=",rr[6],sep=""),socket,sep="\n")
#	 		p = readLines(socket , n=1 ) 
#			m = parser.socket(p)		
#         		writeLines(paste("readsub&num=",m[2],sep=""),socket,sep="\n")
#      	       	        res = readLines( socket , n=1 ) 
#      	                g = parser.socket(res)	
#	                num = g[5]
#		} 
#	
#		writeLines(paste("readloc&num=",num,sep=""),socket,sep="\n")
#	 	r = readLines(socket , n=1 ) 
#		k = parser.socket(r)
#		writeLines(paste("readspec&num=",k[4],sep=""),socket,sep="\n")
# 		o = readLines(socket , n=1 ) 
#	 	return(parser.socket(o)[2])
#				        
#} 
#
