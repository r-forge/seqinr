#include "dir_acnuc.h"

int tlist = 50;
#define MAXLISTES 100
char *deflnames[MAXLISTES]; /* noms listes pour def */
int defoccup[MAXLISTES];
int deflocus[MAXLISTES];
char defgenre[MAXLISTES];
int defllen[MAXLISTES];
int *defbitlist; /* listes de bits pour interf avec routine def */


/* local prototypes */
void prep_acnuc_requete(void);
void free_list(int num);


/* external */
void *mycalloc(int nbre, int size);


void prep_acnuc_requete(void)
/* fonction pour demarrer l'usage d'une banque acnuc */
{
int i;

defbitlist = (int *)mycalloc(tlist*lenw, sizeof(int));
for(i=0; i<MAXLISTES; i++) {
	defoccup[i]=FALSE;
	deflnames[i]=NULL;
	}
lngbit(3,defbitlist);  /* liste des valides */
for(i=0; i< lenw; i++) defbitlist[i]= ~defbitlist[i];
bit0(defbitlist,1);
quick_list_meres(defbitlist+lenw); /* liste des meres */
defoccup[0]=defoccup[1]=TRUE;
deflnames[0]=(char *)mycalloc(12,1);
strcpy(deflnames[0],"!VALIDSEQS!");
deflnames[1]=(char *)mycalloc(12,1);
strcpy(deflnames[1],"!MERES!");
}



void free_list(int num)
/* fonction pour liberer la liste de rang num */
{
defoccup[num]=FALSE;
free(deflnames[num]);
}


