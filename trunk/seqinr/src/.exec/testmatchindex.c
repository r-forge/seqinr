/*
testmatchindex
	Teste adequation fichiers indexs / fichiers plats
	Rend 0 si ok, != 0 sinon
	Accepte divisions d'updates dont des sequences peuvent etre
	indexees ailleurs
	Accepte sequences non indexees
*/
#include "dir_acnuc.h"

/* included functions */
char *last_mnemo(FILE *flat, long *flat_addr, char *marque);
int mal_indexee(int numseq, long addr, int div, char *marque, char *nom);

/* external functions and variables */
void gcgors(char *type_fic, int div, int stop_if_error);
extern FILE *divannot[];


int main(void)
{
int numdiv, numseq, div;
long index_addr, flat_addr;
FILE *flat;
char flat_name[100], *mnemo, marque[10];

dir_acnucopen("RO");
if(genbank) strcpy(marque, "LOCUS ");
else if(embl || swissprot) strcpy(marque, "ID ");
else strcpy(marque, "ENTRY ");
for(numdiv = 0; numdiv <= divisions; numdiv++) {
/* pour chaque division */
	gcgors("inf", numdiv, TRUE);
	flat = divannot[numdiv];
/* travailler jusqu'a la fin de la division */
	fseek(flat, 0, SEEK_END);
	flat_addr = ftell(flat);
	do	{
/* recherche de la derniere sequence avant pos flat_addr */
		mnemo = last_mnemo(flat, &flat_addr, marque);
		if(mnemo == NULL) break; /* pas de seq trouvee */
		numseq = isenum(mnemo);
		if(numseq == 0) { 
			/* si seq non indexee, passer a seq precedante */
			div = numdiv + 1;
			continue;
			}
		/* trouver adresse de cette seq dans les indexs */
		seq_to_annots(numseq, &index_addr, &div);
		/* arreter si seq trouvee est mal indexee */
		if( mal_indexee(numseq, index_addr, div, marque, mnemo) ) exit(ERREUR);
		}
	/* continuer en remontant dans le fichier en ignorant les seqs
	bien indexees dans d'autres divisions ou indexees plus tot dans la meme division
	*/
	while(div != numdiv || index_addr < flat_addr);
	if(mnemo == NULL || div != numdiv || index_addr != flat_addr) exit(ERREUR);
	}
return 0;
}

char *last_mnemo(FILE *flat, long *flat_addr, char *marque)
/* recherche du debut de derniere seq avant pos flat_addr */
{
#define MAXBUF 50000
static char mnemo[30], buffer[MAXBUF + 1];
long fpos;
int lu, trouve, alire;
char *debut, *p;

fpos = *flat_addr - 50;
if(fpos < 0) return NULL;
do	{
/* lecture par paquets de MAXBUF chevauchants de 50 en reculant */
	fpos = fpos - MAXBUF + 50; if(fpos < 0) fpos = 0;
	fseek(flat, fpos, SEEK_SET);
	fpos = ftell(flat);
	alire = MAXBUF;
	if(fpos + alire > *flat_addr) 
		alire = *flat_addr - fpos;
	lu = fread(buffer, 1, alire, flat);
	buffer[lu] = 0;
	debut = buffer - 1; 
/* recherche du dernier debut de seq dans le paquet */
	trouve = FALSE;
	while( ( p = strstr(debut + 1, marque) ) != NULL) {
		debut = p;
		trouve = TRUE;
		}
/* un vrai debut doit etre precede par \n ou \r ou etre 1er caractere du fichier */
	if(trouve) {
		if(fpos == 0 && debut == buffer) ;
		else if(debut == buffer || (*(debut - 1) != '\n' && *(debut - 1) != '\r') )
			trouve = FALSE;
		}
	}
while( !trouve && fpos > 0 );
if( !trouve ) return NULL;
/* calcul position de ce debut dans le fichier */
*flat_addr = fpos + (debut - buffer);
/* recherche du nom de la sequence */
p = debut;
while( *p != ' ') p++;
while( *p == ' ') p++;
debut = mnemo;
while( *p != ' ') *(debut++) = *(p++);
*debut = 0;
return mnemo;
#undef MAXBUF
}


int mal_indexee(int numseq, long addr, int div, char *marque, char *nom)
/* rend TRUE si seq numseq/nom est mal indexee par (addr , div) */
{
char *p;
int l;

p = read_annots(addr, div);
if(p == NULL || strncmp(p, marque, strlen(marque) ) != 0) return TRUE;
while( *p != ' ') p++;
while( *p == ' ') p++;
l = strlen(nom);
return ( p[l] != ' ' || strncmp(p, nom, l ) != 0 );
}


