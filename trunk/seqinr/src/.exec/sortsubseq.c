#include "dir_acnuc.h"
#ifdef unix
#include <unistd.h>
#endif

#define new_rsub(i) ( (i >=1 && i <= nseq) ? seqtable[i] : 0 )

/* globals */
int *seqtable;
int to_h_in_rsub; /* echange entre sortsubseq et sort_subseq_pointers */

/* included functions */
int sort_subseq_pointers(const void *p1, const void *p2);
int sort_record_pointers(const void *p1, const void *p2);
char **merge_table(char **table, int total, int fin_trie);
int sortsubseq(void);
void process_short_list(int point);
void *myalloc(size_t size);

/* external functions */
#ifdef unix
void write_quick_meres(void);
#endif

int main(int argc, char **argv)
{
int tloc, tlng, num, i, text, tacc, tbib;

dir_acnucopen("RW");


/* Sorting SUBSEQ in place */
printf("Sorting SUBSEQ\n"); fflush(stdout);
#ifdef unix
dir_set_mmap(kshrt);
#endif
seqtable = (int *) myalloc( (nseq+1) * sizeof(int) );
sortsubseq();
seqtable[0]=0;

/* processing LONGL */
printf("Processing LONGL\n"); fflush(stdout);
#ifdef unix
dir_set_mmap(klng);
#endif
tlng = read_first_rec(klng, NULL);
memset(plng, 0, lrlng); /* mettre liste effacees a vide */
plng->sub[0] = 1;
writelng(3);
for(num = 2; num <= tlng; num++) {
	readlng(num);
	for(i = 0; i < SUBINLNG; i++)
		plng->sub[i] = new_rsub(plng->sub[i]);
	writelng(num);
	}
dir_close(klng);

/* processing LOCUS */
printf("Processing LOCUS\n"); fflush(stdout);
#ifdef unix
dir_set_mmap(kloc);
#endif
tloc = read_first_rec(kloc,NULL);
for(num = 2; num <= tloc; num++) {
	readloc(num);
	ploc->sub = new_rsub(ploc->sub);
	writeloc(num);
	}
dir_close(kloc);

if( !(nbrf || swissprot) ) {
	/* processing EXTRACT */
	printf("Processing EXTRACT\n"); fflush(stdout);
	#ifdef unix
	dir_set_mmap(kext);
	#endif
	text = read_first_rec(kext,NULL);
	for(num = 2; num <= text; num++) {
		readext(num);
		pext->mere = new_rsub(pext->mere);
		writeext(num);
		}
	dir_close(kext);
	}

printf("Processing SHORTL\n"); fflush(stdout);

/* processing ACCESS */
tacc = read_first_rec(kacc,NULL);
for(num = 2; num <= tacc; num++) {
	readacc(num);
	process_short_list(pacc->plsub);
	}

/* processing BIBLIO */
tbib = read_first_rec(kbib,NULL);
for(num = 2; num <= tbib; num++) {
	readbib(num);
	process_short_list(pbib->plsub);
	}

/* processing sequence hashing */
for(num = 3; num <= (hsub + 1)/2 + 2; num++) {
	readshrt(num);
	pshrt->val = new_rsub(pshrt->val); pshrt->next = new_rsub(pshrt->next);
	writeshrt(num);
	}

dir_close(kshrt);
#ifdef unix
dir_close(kspec);
dir_close(kkey);
dir_close(ksmj);
dir_close(kacc);
dir_close(ktxt);
dir_close(kaut);
dir_close(kbib);
write_quick_meres();
#endif
printf("Normal end\n");
return 0;
}


int sort_subseq_pointers(const void *p1, const void *p2)
{
char *c1, *c2, *last, *deb1, *deb2;
int h1, h2;

c1 = *(char **)p1;
c2 = *(char **)p2;
deb1 = c1; 
deb2 = c2;
last = c1 + (L_MNEMO - 1);
while(*c1 == *c2) {
	if(*c1=='.') {
		/* rendre ordre des pointeurs vers location */
		h1 = *(int *)( deb1 + to_h_in_rsub);
		h2 = *(int *)( deb2 + to_h_in_rsub);
		if( (unsigned int) h1 > (unsigned int) h2 )
			return 1;
		else if ( (unsigned int) h1 < (unsigned int) h2 )
			return -1;
		else
			return 0;
		}
	if(c1 >= last) return 0;
	c1++; c2++;
	}
return *c1 - *c2;
}


int sort_record_pointers(const void *p1, const void *p2)
{
char *rec1, *rec2;
rec1 = *(char **)p1; rec2 = *(char **)p2;
return memcmp(rec1, rec2, L_MNEMO);
}


char **merge_table(char **table, int total, int fin_trie)
{
char **new_list;
int rank_1, rank_2, rank_new;

new_list = (char **)myalloc(total * sizeof(char *));
rank_1 = 0; rank_2 = fin_trie; rank_new = 0;
while(rank_new < total) {
	while( rank_1 >= fin_trie || memcmp(table[rank_1], 
			table[rank_2], L_MNEMO) > 0) {
		new_list[rank_new++] = table[rank_2++];
		if(rank_2 >= total) 
			break;
		}
	if(rank_new >= total) break;
	while( rank_2 >= total || memcmp(table[rank_1], 
			table[rank_2], L_MNEMO) < 0) {
		new_list[rank_new++] = table[rank_1++];
		if(rank_1 >= fin_trie) 
			break;
		}
	}
free(table);
return new_list;
}


int sortsubseq(void)
{
char *array, croix[L_MNEMO], *debut, *fin, **deb_record, *previous;
int insize, outsize, old_rank, new_rank, num, trie, fin_trie;
char *sub_buffer;
int *h_field, *htable;

/* offset vers champ h ds record rsub */
to_h_in_rsub = ( (char *)&(psub->h) - (char *)psub );
memset(croix, 'x', L_MNEMO);
insize = read_first_rec(ksub, NULL);
array = (char *)myalloc(insize * lrsub);
htable = (int *)myalloc((nseq + 1) * sizeof(int));
/* lecture en memoire de tout le fichier a trier */
dir_read(ksub, 2, insize - 1, array);
if(!( nbrf || swissprot) ) { 
/* cas du fichier SUBSEQ et banque nucleique:
on veut pour ordre des filles celui de leur citation dans la table des features
on remplace le champ h inutile de SUBSEQ par le pointeur vers feature
*/
/* offset vers 2 champs ds record rsub */
	int to_pext_in_rsub, to_plinf_in_rsub, *pext_field, *plinf_field;
	to_pext_in_rsub = ( (char *)&(psub->pext) - (char *)psub );
	to_plinf_in_rsub = ( (char *)&(psub->plinf) - (char *)psub );
	debut = array; fin = array + (insize - 1) * lrsub;
	old_rank = 2;
	while(debut < fin) {
		/* trouver champ h */
		h_field = (int *)(debut + to_h_in_rsub); 
		/* conserver ancienne valeur de h */
		htable[old_rank++] = *h_field;
		/* trouver champ pext */
		pext_field = (int *)(debut + to_pext_in_rsub); 
		if( *pext_field > 0 ) { /* pour une fille (pext > 0) */
			/* trouver champ plinf */
			plinf_field = (int *)(debut + to_plinf_in_rsub); 
			num = *plinf_field; /* lire champ plinf */
			/*lire ds SHORTL adresse de la location de cette fille*/
			readshrt(num); 
			/* stocker cette location ds champ h inutile */
			*h_field = pshrt->val; 
			}
		debut += lrsub;
		}
	}
deb_record = (char **)myalloc(sizeof(char *) * insize);
debut = array; fin = array + (insize - 1) * lrsub;
outsize = 0;
while(debut < fin) {
	if(strncmp(debut, croix, L_MNEMO) != 0) {
		deb_record[outsize++] = debut;
		}
	debut += lrsub;
	}
if( !( nbrf || swissprot ) ) {
	qsort(deb_record, outsize, sizeof(char *), sort_subseq_pointers);
	}
else	{
/* calculer longueur de la partie deja triee */
	previous = deb_record[0];
	for(fin_trie = 1; fin_trie < outsize; fin_trie++) {
		trie = (memcmp(deb_record[fin_trie], previous, L_MNEMO) > 0);
		if( !trie ) break;
		previous = deb_record[fin_trie];
		}
	if(fin_trie > outsize) fin_trie = outsize;
	if(fin_trie < 2) fin_trie = 0;
/* trier la partie non triee du tableau deb_record */
	if( outsize - fin_trie > 1) qsort(deb_record + fin_trie, 
		outsize - fin_trie, sizeof(char *), sort_record_pointers);
/* trier tout deb_record en mergeant ses 2 parties */
	if( outsize > fin_trie && fin_trie > 0 ) 
		deb_record = merge_table(deb_record, outsize, fin_trie);
	}
/* calculer seqtable */
memset(seqtable, 0, (insize + 1) * sizeof(int) );
seqtable[1] = 1;
new_rank = 1;
for(num = 0; num < outsize; num++) {
	old_rank = (deb_record[num] - array) / lrsub + 2;
	seqtable[old_rank] = ++new_rank;
	}

/* ecrire SUBSEQ sous forme triee */
num = (outsize + 1) * lrsub;
if(num > lrsub * 10000) num = lrsub * 10000;
sub_buffer = (char *)myalloc(num);
dir_resize_buff(ksub, sub_buffer, num);
write_first_rec(ksub, outsize + 1, outsize + 1);
new_rank = 1;
for(num = 0; num < outsize; num++) {
	old_rank = (deb_record[num] - array) / lrsub + 2;
	/* traiter champ h */
	h_field = (int *)(deb_record[num] + to_h_in_rsub); 
	if(! (swissprot || nbrf) ) {
		*h_field = new_rsub(htable[old_rank]);
		}
	else 	{
		*h_field = new_rsub(*h_field);
		}
	dir_write(ksub, ++new_rank, 1, deb_record[num]);
	}
#ifdef unix
ftruncate(ksub->fd, (outsize + 1) * ksub->record_length );
#endif
dir_close(ksub);
free(array); free(deb_record); free(htable);
return outsize + 1;
}


void process_short_list(int point)
{
while(point != 0) {
	readshrt(point);
	pshrt->val = new_rsub(pshrt->val);
	writeshrt(point);
	point = pshrt->next;
	}
}


void *myalloc(size_t size)
{
void *retval;

retval = malloc(size);
if(retval == NULL) {
	fprintf(stderr, "Not enough memory\nBank left unchanged\n");
	exit(ERREUR);
	}
return retval;
}


