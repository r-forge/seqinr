#include "mase.h"
#include <Rinternals.h>
#include <R.h>
#include <Rdefines.h>



/***********************************************************************************************************************/
/*  lit un fichier MASE, renvoie une liste (objet R) contenant les séquences, les commentaies et les noms des espèces. */
/***********************************************************************************************************************/


SEXP read_mase(SEXP nomfic)
{
  char *fic_name;
  FILE *fic;
  struct SEQMASE *aln;
  int nb_seq;
  int lg_max = 0, lg, lgs, lgc;
  char string[MAXSTRING + 1];
  char c1, c2;
  int i,ii, jj, kk, numline, maxcom = 0;
 
  SEXP listseq;
  SEXP essai;
  SEXP listcom;
  SEXP listmn;
  SEXP nombreseq;
  SEXP nombresite;

  /*Passages des objets R (paramètres) dans des variables C */ 
  fic_name=CHAR(STRING_ELT(nomfic,0)); 



  if((fic = fopen(fic_name, "r")) == NULL) {
    error("impossible to open file");
  }
  
  
  c1 = 0;
  nb_seq = 0;
  lg = lgc = 0;
  while(fgets(string, MAXSTRING, fic) != NULL) {
    string[MAXSTRING] = 0;
   
    lgs = strlen(string);
    
    if(lgs >= (MAXSTRING - 1)) {
      fprintf(stderr, "\n Too long line in alignment (> %d).\n", MAXSTRING);
      error("increase MAXSTRING and recompile.\n"); 
    }
    
    c2 = string[0];
    
    if(string[0] == ';' && string[1] != ';') {
      lgc += (lgs + 1);
    }
    
    
    if(c1 == ';' && c2 != c1) {
      nb_seq++;
     if(lg > lg_max) lg_max = lg;
     if(lgc > maxcom) maxcom = lgc;
     lg = lgc = 0;
    }
    
    else if(c2 != ';') lg += lgs;
    c1 = c2;
    
  }
  if(lg > lg_max) lg_max = lg;
  

  /******************************************/
  /* Création de 6 objets R qui seront      */
  /******************************************/

  PROTECT(listseq=allocVector(VECSXP,nb_seq));
  PROTECT(essai=allocVector(VECSXP,5));
  PROTECT(listcom=allocVector(VECSXP,nb_seq));
  PROTECT(listmn=allocVector(VECSXP,nb_seq));
  PROTECT(nombreseq=NEW_INTEGER(1));
  PROTECT(nombresite=NEW_INTEGER(1));
 
 aln = (struct SEQMASE *) check_alloc(nb_seq + 1, sizeof(struct SEQMASE));
 
 for(ii = 0; ii <= nb_seq; ii++) {
    aln[ii].seq = (char *) check_alloc(lg_max + 1, sizeof(char));
    aln[ii].com = (char *) check_alloc(maxcom + 1, sizeof(char));
    aln[ii].com[0] = 0;
 }
 

 rewind(fic);

 numline = 0;
 ii = -1;
 while(fgets(string, MAXSTRING, fic) != NULL) {
   numline++;
   string[MAXSTRING] = 0;
   
   
   c2 = string[0];
   if(string[0] == ';' && string[1] != ';') {
     strcat(aln[ii + 1].com, string);
        }
   
   if(c1 == ';' && c2 != c1) {
     ii++;
     kk = aln[ii].lg = 0;
	
     rem_blank(string);

     if((int) strlen(string) >= (MAXMNMASE - 1)) {
       fprintf(stderr, "\nEXIT: sequence name too long!\n");
       fprintf(stderr, "max  %d characters\n", MAXMNMASE);
       fprintf(stderr, "(line %d)\n", numline);
       exit(1);
     }
     strcpy(aln[ii].mn, string);
     
     lg = 0;
   }
   
   else if(c2 != ';') {
     for(jj = 0; jj < MAXSTRING; jj++) {
       if(string[jj] == 0) break;
       if(string[jj] == ' ') continue;
       if(string[jj] == '\n') continue;
       if(string[jj] == '\t') continue;
            aln[ii].seq[kk++] = string[jj];
            aln[ii].lg = kk;
     }
   }
   c1 = c2;
   
 }
 

 fclose(fic);
 
 lg_max = aln[0].lg;

 for(ii = 1; ii < nb_seq; ii++)
   if(aln[ii].lg > lg_max) lg_max = aln[ii].lg;
 
 fprintf(stderr,
	 "%d sequences, max length = %d sites\n", nb_seq, lg_max);
 fflush(stderr);
 

 INTEGER(nombreseq)[0]=(int)nb_seq;
 INTEGER(nombresite)[0]=(int)lg_max;

 for(i=0;i<nb_seq;i++){
   SET_ELEMENT(listseq,i,mkChar(aln[i].seq));
 }

 for(i=0;i<nb_seq;i++){
   SET_ELEMENT(listcom,i,mkChar(aln[i].com));
 }

for(i=0;i<nb_seq;i++){
   SET_ELEMENT(listmn,i,mkChar(aln[i].mn));
 }
 
 SET_ELEMENT(essai,0,nombreseq);
 SET_ELEMENT(essai,1,nombresite);
 SET_ELEMENT(essai,2,listseq);
 SET_ELEMENT(essai,3,listcom);
 SET_ELEMENT(essai,4,listmn);
  

 UNPROTECT(6);
 free_mase(aln,nb_seq);

 return(essai);

}


/********************** end read_mase ****************************/

  
char *check_alloc(int nbrelt, int sizelt)
{
  char *retval;
  if( (retval=calloc(nbrelt,sizelt)) != NULL ) return(retval);
  fprintf(stderr,"\nERROR: not enough memory (check_alloc).\n");
  exit(1);
}

/********************** end check_alloc ****************************/


void rem_blank(char *string)
{
  int ii;


  ii = strlen(string);

  for( ;ii >=0; ii--) {
    if(string[ii] == 0 || string[ii] == '\n' ||
       string[ii] == ' ' || string[ii] == '\t') string[ii] = 0;
    else break;
  }
  

}
/**************************** end rem_blank ************************/

void free_mase(struct SEQMASE * aln, int nbsq)
     
{
  int ii;


  for(ii = 0; ii <= nbsq; ii++) {
    free(aln[ii].seq);
    free(aln[ii].com);
  }
  
  free((char *) aln);
  
  
}

/******************************** end free_mase ************************/


/*************************************************************************/
/* Compute distance between two aligned sequences using different matrix */
/*************************************************************************/


SEXP distance(SEXP sequences,SEXP nbseq, SEXP matNumber){

  SEXP d;
  char **seq;
  double dist[MAXNSEQS][MAXNSEQS];
  int i, j, k, n,totseqs, seq_long, nbases;
  int ndiff[MAXNSEQS][MAXNSEQS];
  int mat_number;

  int mat_pos[]  = { 17, -1, 15, 0, 1, 12, 18,  4,  9, -1,  2, 10, 16, 5,
                       -1, 19,  6, 3, 7,  8, -1, 11, 13, -1, 14, -1 };
  const char Prot[] = "DEKRHNQSTILVFWYCMAGPX*-";
  int matp[20][20] = { {1, 1, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3},
		       {1, 1, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3},
		       {2, 2, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3},
		       {2, 2, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3},
		       {2, 2, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3},
		       {2, 2, 2, 2, 2, 1, 1, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3},
		       {2, 2, 2, 2, 2, 1, 1, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3},
		       {2, 2, 2, 2, 2, 2, 2, 1, 1, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3},
		       {2, 2, 2, 2, 2, 2, 2, 1, 1, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3},
		       {3, 3, 3, 3, 3, 3, 3, 3, 3, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3},
		       {3, 3, 3, 3, 3, 3, 3, 3, 3, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3},
		       {3, 3, 3, 3, 3, 3, 3, 3, 3, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3},
		       {3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 1, 1, 1, 2, 2, 3, 3, 3},
		       {3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 1, 1, 1, 2, 2, 3, 3, 3},
		       {3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 1, 1, 1, 2, 2, 3, 3, 3},
		       {3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 1, 2, 3, 3, 3},
		       {3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 1, 3, 3, 3},
		       {3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 1, 1, 3},
		       {3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 1, 1, 3},
		       {3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 1} };
  

 
 
  totseqs = INTEGER_VALUE(nbseq);
  mat_number= INTEGER_VALUE(matNumber);

  PROTECT(d=NEW_NUMERIC(totseqs*totseqs));

  seq = (char **)malloc(totseqs*sizeof(char *));
   
  for(i=0;i<=totseqs;i++){
    seq[i] = CHAR(STRING_ELT(sequences,i));
  }

  /**************************************/
  /* Initialising matrix                */
  /*************************************/


  if (mat_number == 1){
      for (i = 0; i < 20; i++){
	for (j = 0; j < 20; j++){
	  matp[i][j] = matp[i][j]/3.0;
	}
      }
    }
  

  /*********************************************************/
  /*     Computing distance between sequences i and j      */
  /*********************************************************/

    seq_long = (int)strlen(seq[0]);
    
    for (i = 0; i < totseqs; i++)
      {
        dist[i][i] = 0.0;
      }

    for (i = 0; i < totseqs; i++)
    {
      for (j = 0; j < i; j++)
        {
	  ndiff[j][i] = ndiff[i][j] = 0;
	  nbases = 0;
	  for (k = 0; k < seq_long; k++)
            {

	      if(mat_number == 1){

		/********************************************/
		/*Protein sequences with similarity matrix  */
		/********************************************/

		  if ((strchr(Prot, seq[i][k]) == NULL) || (strchr(Prot, seq[j][k]) == NULL))
		    seq[i][k] = seq[j][k] = 'X';
		  if (seq[i][k] != '-' && seq[i][k] != '*' && seq[i][k] != 'X' &&
		      seq[j][k] != '-' && seq[j][k] != '*' && seq[j][k] != 'X')
                    {
		      nbases++;
		      ndiff[i][j] = ndiff[i][j] + matp[mat_pos[seq[i][k] - 65]][mat_pos[seq[j][k] - 65]];
		      ndiff[j][i] = ndiff[i][j];
                
		    }
	      }
	      else if( mat_number == 2){

		/********************************************/
		/*  Protein sequences with identity matrix  */
		/********************************************/

		if ((strchr(Prot, seq[i][k]) == NULL) || (strchr(Prot, seq[j][k]) == NULL))
		  seq[i][k] = seq[j][k] = 'X';
		if (seq[i][k] != '-' && seq[i][k] != '*' && seq[i][k] != 'X' &&
		      seq[j][k] != '-' && seq[j][k] != '*' && seq[j][k] != 'X')
		  {
		      nbases++;
		      if (seq[i][k] != seq[j][k])
                        {
			  ndiff[j][i]++;
			  ndiff[i][j]++;
                        }
		  }	
		
	      }
	      
	      dist[i][j] = dist[j][i] = sqrt((double)ndiff[i][j]/nbases);
	    }
	}
    }
    
    
    /********************************************************************************/
    /* Remplissage de l'objet R (matrice de taille nb_seq * nb_seq  avec dist       */
    /********************************************************************************/

    n=0;
    
    for(i=0;i<totseqs;i++){
      for(j=0;j<totseqs;j++){
	REAL(d)[n+j]=dist[i][j];
      }
      n=n+totseqs;
    }

    UNPROTECT(1);
    
    return(d);
}
