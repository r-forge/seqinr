#include "requete_acnuc.h"
#include "dir_acnuc.h"
#include <Rinternals.h>
#include <R.h>
#include <Rdefines.h>





void Racnucopen(char **acnuc,char **gcgacnuc){

  chg_acnuc(*acnuc,*gcgacnuc);
 
  /*ouverture complete en lecture de la banque acnuc designées par les variables acnuc et gcgacnuc*/
  acnucopen(); 
  prep_acnuc_requete();
}

void Racnucclose(){

  dir_acnucclose();
}

void getdesc(char **name,char **seq,char **p){

  int tot,num, lseq, deb;
  
  /*obtenir le numero de la sequence dans acnuc à partir de son nom ou mnemo*/
  num=gsnuml(*name,&lseq,NULL,NULL);
  if(num==0) error("ce mnemo est invalide");
	
  /*acces à la séquence */
  tot=gfrag(num,1,lseq,*seq);
  *p=short_descr(num,*seq,80);
}



SEXP getreq(SEXP nom,SEXP req){

  char *nom_liste;
  char *requete;
  char message[100];
  int numliste,num,total,*debut_liste,erreur,i;
  SEXP chaine;
  
  nom_liste=CHAR(STRING_ELT(nom,0));
  requete=CHAR(STRING_ELT(req,0));
  
  erreur = proc_requete(requete, message, nom_liste, &numliste);
  
  if(erreur){
    printf("Message d'erreur:\n%s\n\n",message);
    error("renouvelez la requete");
  }

  total = defllen[numliste];
  //printf("\nNbre de seqs selectionnees = %d\n", total);

  PROTECT(chaine=allocVector(VECSXP,total));

  debut_liste = defbitlist + numliste*lenw;
  num = 1;   
  total=1; 
	while( (num = irbit(debut_liste, num, nseq)) != 0)  {
	 readsub(num);
	 SET_ELEMENT(chaine,total-1,mkChar(psub->name));
	 total++;
       }
	 UNPROTECT(1);
	 return(chaine);
	 }



SEXP getseq(SEXP nom){
  //char seq[10001];
  char *seq;
  char *name;
  int tot,num,deb,lseq;
  SEXP chaine;

  name=CHAR(STRING_ELT(nom,0));

  /*obtention du numero de la sequence dans acnuc, rendu avec sa longueur, sa frame (ici NULL) et son code genetique (ici NULL)*/
  num=gsnuml(name,&lseq,NULL,NULL);
  if(num==0) error("ce mnemo est invalide");

  seq=(char *)malloc(lseq*sizeof(char));
  PROTECT(chaine=NEW_CHARACTER(1));
  /* acces à la séquence de numéro num de la position 1 à lseq et stockage dans seq */
  tot=gfrag(num,1,lseq,seq);
  SET_ELEMENT(chaine,0,mkChar(seq));
  UNPROTECT(1);
  return(chaine);
}

SEXP getseq2(SEXP nom, SEXP b1, SEXP b2){
  //char seq[10001];
  char *seq;
  char *name;
  double borne1,borne2;
  int tot,num,deb,lseq;
  SEXP chaine;

  borne1=REAL(b1)[0];
  borne2=REAL(b2)[0];
  name=CHAR(STRING_ELT(nom,0));

  /*obtention du numero de la sequence dans acnuc, rendu avec sa longueur, sa frame (ici NULL) et son code genetique (ici NULL)*/
  num=gsnuml(name,&lseq,NULL,NULL);
  if(num==0) error("ce mnemo est invalide");

  seq=(char *)malloc(lseq*sizeof(char));
  PROTECT(chaine=NEW_CHARACTER(1));
  /* acces à la séquence de numéro num de la position 1 à lseq et stockage dans seq */
  tot=gfrag(num,borne1,borne2,seq);
  SET_ELEMENT(chaine,0,mkChar(seq));
  UNPROTECT(1);
  return(chaine);
}
