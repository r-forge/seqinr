choosebank <- function( bankname )
{
  if( missing( bankname ) )
    stop("Mandatory argument is missing")

  ad <- Sys.getenv( bankname )
  if( ad == "" )
    stop(paste("Unable to get environment variable ", bankname))

  var <- unlist(strsplit(ad, " "))
  acnuc <- var[1]
  gcgacnuc <- var[2]
  .C("Racnucopen", acnuc, gcgacnuc)
}

choixbanque <- choosebank # Just an alias for compatibility


getreq <- function(nomliste,requete)
{	
	liste<-.Call("getreq",nomliste,requete)
	liste<-as.character(liste)
	liste<-strsplit(liste,"")
	p<-function(m) paste(m[m!=" "],collapse="")
	liste<-lapply(liste,p)
	toto=list(call=requete,name=nomliste,req=liste)
	class(toto)=c("requete")
	assign(nomliste,toto,env = .GlobalEnv)
	print.requete(toto)
}



getseq <- function(name, as.string = TRUE)
{
 	if(as.string == TRUE){
 	x=.Call("getseq",name)
	class(x)=sequence
	x
	}
	else{
  	tata<-.Call("getseq",name)
  	seq <- unlist(strsplit(tata,""))
	class(tata)=sequence
  	toupper(seq)
	}
}

getseq2 <- function(name,B1,B2,as.string = TRUE)
{
 	if(as.string == TRUE){
 	.Call("getseq2",name,B1,B2)
	}
	else{
 	tata<-.Call("getseq2",name,B1,B2)
 	seq <- unlist(strsplit(tata,""))
  	toupper(seq)
	}
}


changebanque <- function(nombanque)
{
	.C("Racnucclose")
	choixbanque.R(nombanque)
}


print.requete <- function(x)
{
	cat("class: ")
	cat(class(x))
	cat("\nbanque: ")
	cat(Sys.getenv("acnuc"))
	cat("   ")
	cat(Sys.getenv("gcgacnuc"))
	cat("\n")
	cat("\n$call: ")
	print(x$call)
	cat("\n$name: ")
	print(x$name)
	cat("\n")
	sumry <- array("", c(1, 4), list(1, c("list", "length", 
        "mode", "content")))
	sumry[1, ] <- c("$req",length(x$req),"character","mnemoniques")
	class(sumry) <- "table"
	print(sumry)
	cat("\n")
}

translateCDS <- function(name)
{
	.Call("translateCDS",name)
}

