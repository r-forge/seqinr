/*################################################################################################# */
/* Programme: R_acnuc.c                                                                             */
/*                                                                                                  */
/* But: Fonctions pour l'interfaçage entre R et ACNUC                                               */
/*                                                                                                  */
/* Fonctions disponibles: Racnucopen, Racnucclose, getreq, getseq, getseq2, getAttribut, getExon,   */
/*                        translateCDS, getAnnots, s2c, getKey.                                     */
/*                                                                                                  */
/* Auteur: Delphine Charif                                                                          */
/*                                                                                                  */
/*################################################################################################# */




#include ".exec/requete_acnuc.h"
#include ".exec/dir_acnuc.h"
#include <Rinternals.h>
#include <R.h>
#include <Rdefines.h>
#define VRAI 1
#define FAUX 0
#define MAXMOT 100



/*################################################################################################# */
/*#  Ouverture du système ACNUC                                                                     */
/*################################################################################################# */


void Racnucopen(char **acnuc,char **gcgacnuc){

  chg_acnuc(*acnuc,*gcgacnuc);
 
  /*ouverture complete en lecture de la banque acnuc designées par les variables acnuc et gcgacnuc*/
  acnucopen(); 
  prep_acnuc_requete();
}



/*################################################################################################# */
/*#  Fermeture du système ACNUC                                                                   # */
/*################################################################################################# */



void Racnucclose(){

  dir_acnucclose();
}


/*################################################################################################# */
/*#   Effectuer une requete                                                                       # */
/*################################################################################################# */


SEXP getreq(SEXP nom,SEXP req){

  char *nom_liste;
  char *requete;
  char message[100];
  int numliste,num,total,*debut_liste,erreur,i;
  SEXP chaine;
  
  /*Passages des objets R (paramètres) dans des variables C */ 
  nom_liste=CHAR(STRING_ELT(nom,0));
  requete=CHAR(STRING_ELT(req,0));
  
  erreur = proc_requete(requete, message, nom_liste, &numliste);
  
  if(erreur){
    printf("Message d'erreur:\n%s\n\n",message);
    error("renouvelez la requete");
  }

  total = defllen[numliste];
  
  /*Création d'une liste de chaine de caractère, objet R destiné à contenir les mnémoniques*/ 
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

/*############################################################################################################################*/
/*#   Récupération des attributs d'une séquence à partir de son mnémonique : longueure, phase de lecture, code génétique     #*/
/*############################################################################################################################*/


SEXP getAttribut(SEXP nom){
  char *name;
  int lseq,num,frame,gencode;
  SEXP l;

  name=CHAR(STRING_ELT(nom,0));
  num=gsnuml(name,&lseq,&frame,&gencode);
  PROTECT(l = NEW_INTEGER(3));
  INTEGER(l)[0]=lseq;
  INTEGER(l)[1]=frame;
  INTEGER(l)[2]=get_ncbi_gc_number(gencode);  
  UNPROTECT(1);
  return(l);
}

/*############################################################*/
/*#   Récupération de la séquence à partir de son mnémonique #*/
/*############################################################*/


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

/*##########################################################################*/
/*#   Récupération d'un fragment de la séquence à partir de son mnémonique #*/
/*##########################################################################*/



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

  /*obtention du numero de la sequence dans acnuc, rendu avec sa longueure,
  sa frame (ici NULL) et son code genetique (ici NULL)*/
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


/*##############################*/
/*# Traduction de la séquence  #*/
/*##############################*/



SEXP translateCDS(SEXP nom){
  
  char *protein;
  char *name,*seq;
  int frame, code, num, lseq;
  SEXP chaine;

  name=CHAR(STRING_ELT(nom,0));
  num=gsnuml(name,&lseq,&frame,&code);
  seq=(char *)malloc(lseq*sizeof(char));
  PROTECT(chaine=NEW_CHARACTER(1));
  protein = translate_cds(num);
  SET_ELEMENT(chaine,0,mkChar(protein));
 
  UNPROTECT(1);
  return(chaine);

  }

/*####################################################*/
/*# Obtention des mots clefs associés un mnémonique  #*/
/*####################################################*/


SEXP getKey(SEXP nom){
  
  char *mnemo;
  int num, kw, point,total,i;
  
  SEXP chaine;
  total=0;
  i=0;

  mnemo = CHAR(STRING_ELT(nom,0));
 
  num = isenum(mnemo);
  
  readsub(num);
  point = psub->plkey;
  
  while(point !=0){
    total++;
    readshrt(point);
    kw = pshrt->val;
    point = pshrt->next;
  }

  PROTECT(chaine=allocVector(VECSXP,total));

  point = psub->plkey;

 while(point !=0 && i<=total){
    readshrt(point);
    kw = pshrt->val;
    readkey(kw);
    SET_ELEMENT(chaine,i,mkChar(pkey->name));
    point = pshrt->next;
    i++;
  }
  UNPROTECT(1);
  return(chaine);
}


/*#######################*/
/*# Position des Exons  #*/
/*#######################*/


SEXP getExon(SEXP nom){

 char *mnemo;
 int deb,fin,ss,nsub,i,tmp;
 int total=0;

 SEXP Exon;
 
 mnemo = CHAR(STRING_ELT(nom,0));

 nsub = isenum(mnemo);

 if(nsub == 0) error("Ce mnemo est invalide");
     
 /* lecture de l'enregistrement correspondant */
 readsub(nsub);

/* si le CDS est une mere, alors pas de decoupage en exon */

 if(psub->pext <= 0) {
   error(" pas un CDS : sequence mere\n");
 }

 /* sinon s'il s'agit d'une subseq, parcours du decoupage du CDS */
 else {
   
   ss = psub->pext;
   /*parcours de la liste chainée pour connaitre nombre exon*/
   while (ss != 0) {
     total++;
     readext(ss);
     ss =  pext->next;
   }

   PROTECT(Exon=NEW_INTEGER(2*total));

   /*réinitialisation de la variable ss*/
   ss = psub->pext;
   
   i=0;
   while (ss != 0 && i<=(2*total)) {
     readext(ss);
     deb = pext->deb;
     fin = pext->fin;
     if(deb > fin){
       tmp = fin;
       fin = deb;
       deb = tmp;
     }
     INTEGER(Exon)[i]=deb;
     i++;
     INTEGER(Exon)[i]=fin;
     i++;
     ss =  pext->next;
   }
 }
 UNPROTECT(1);
 return(Exon);
}

/*########################################################*/
/*# Obtention des annotations associées à un mnémonique  #*/
/*########################################################*/


SEXP getAnnots(SEXP name, SEXP nligne){
  
  int nsub,div,i,lignes;
  long annot;
  char *mnemo;
 
  SEXP chaine;
  
  lignes = INTEGER_VALUE(nligne);
  PROTECT(chaine=allocVector(VECSXP,lignes));
  mnemo = CHAR(STRING_ELT(name,0));
 
  nsub = isenum(mnemo);
  if(nsub == 0) error("Ce mnemo est invalide");
  
  seq_to_annots(nsub, &annot, &div);
  
  SET_ELEMENT(chaine,0,mkChar(read_annots(annot,div)));
	      
	      for(i=1;i<lignes;i++){
		SET_ELEMENT(chaine,i,mkChar(next_annots(NULL)));
	      }
	      UNPROTECT(1);
	      return(chaine);
}	
 







