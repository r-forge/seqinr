
mkbank = function(bankname,fileseq,bank="genbank",format="flatfiles"){

#
# Vérification du format pour le fichier de séquence
#
	if(any(grep("\\.",fileseq))) b=which(s2c(fileseq)==".") else(stop("you don't have any file extension !"))
	b=which(s2c(fileseq)==".")
	if(substr(fileseq,b+1,b+3)!="dat" & substr(fileseq,b+1,b+3)!="seq") stop("you don't have the good extension for the file of sequences !") 
#	
# Création des répertoires flatfiles et index qui vont accueillir la banque ACNUC
#
	path1=paste(getwd(),bankname,sep=.Platform$file.sep)
	path2=paste(path1,"flatfiles",sep=.Platform$file.sep)
	path3=paste(path1,"index",sep=.Platform$file.sep)
	path4=system.file("exec/initf",package="seqinr")
	path5=system.file("exec/gbemgener",package="seqinr")
	dir.create(path1)
	dir.create(path2)
	dir.create(path3)
	file.copy(fileseq,path2)
#
# Génération de la banque vide par le programme initf
#
	if(format=="gcg") form=format
	com=paste(paste(path4,bank,form,sep=" "),paste("div",substr(fileseq,1,b-1),sep="="),sep=" ")
	system(paste(paste("cd",path3,sep=" "),com,sep=" ;"))

#
# Définition des nouvelles variables d'environnements
#
	Sys.putenv("acnuc"=path3)
	Sys.putenv("gcgacnuc"=path2)
#
# Remplissage des fichiers indexs et génération de la banque structurée ACNUC par le programme gbemgener
#	
	system(paste(paste(path5,"d",sep=" "),substr(fileseq,1,b-1),sep=" "))

#
# Création d'une variable d'environnement portant le nom de la banque et indiquant l'emplacement des fichiers plats et indexs.
#
	x=paste(bankname, paste(path3,path2,sep=" "), sep = "=")
        invisible(.Internal(putenv(x)))
}

                                                                                                                                                                                                                                                                                                                                                                                                                                                         