#include <Rinternals.h>
#include <R.h>
#include <Rdefines.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>




/*##################################################*/
/*# Convertir une string en vecteur de caract�res  #*/
/*##################################################*/


SEXP s2c(SEXP seq){
  char *string;
  int lseq,i;
  char mot[2];
  
  
  SEXP chaine;

  string = CHAR(STRING_ELT(seq,0));
  
  lseq = strlen(string);
  
  PROTECT(chaine=NEW_CHARACTER(lseq));


  for(i=0;i<lseq;i++){  
    mot[0]=string[i];
    mot[1]='\0';
    SET_ELEMENT(chaine,i,mkChar(mot));
    }
  UNPROTECT(1);
  return(chaine);
}


/*##################################################*/
/*# Mettre de nouvelles variables d'environnement  #*/
/*##################################################*/


void Sysputenv(char **varname, char **valeur){
  
  int i,l;  
  char *vartot;
 
  l=strlen(*varname)+strlen(*valeur)+2;
  vartot = malloc(l);
  strcpy(vartot,*varname); 
  strcat(vartot,"=");
  strcat(vartot,*valeur);
  printf("%s\n",vartot);
  i=putenv(vartot);
  printf("%d\n",i);
}

 
/*#####################################################*/
/*# Tester si une s�quence est prot�ique ou nucl�ique #*/
/*#####################################################*/




SEXP is_a_protein_seq(SEXP sequence)
/* returns TRUE if seq looks like a protein sequence (less than 80% ACGTU) */
{

  SEXP res;
  char *seq;
  static char dna[]="ACGTU";
  int total=0, length=0;
  
  seq = CHAR(STRING_ELT(sequence,0)); 

  while(*seq != 0) {
    if(*seq != '-') {
      if( strchr(dna, toupper(*seq)) != NULL ) total++;
      length++; 
    }
    seq++;
  }
  
   PROTECT(res=NEW_NUMERIC(1));
   REAL(res)[0]=(float)(total) / length ;
   
   UNPROTECT(1);
   return ( res );
}


