#include "dir_acnuc.h"

#ifdef unix
void write_quick_meres(void);
#endif


DIR_FILE *klng2, *ksub2, *kan11, *ktxt2, *kshrt2, *kloc2, *kaut2;
int *acctable, *loctable, *seqtable, *bibtable, *autable, *smjtable, *htable,
	newsub, lastlng2, compsize, *lng_list_aux;
unsigned lastshrt2;
struct paire {
	int old_val, new_val;
	} *tabpairloc;
struct rspec ksrec;


void mem_message(char *text)
{
fprintf(stderr,"Error: %s\n",text);
exit(ERREUR);
}


void check_desc_pointers(DIR_FILE *kan, char *type)
{
int total, num;
char message[100];

total=read_first_rec(kan, NULL);
for(num = 2; num <= total; num++) {
	dir_read(kan, num, 1, &ksrec);
	if(ksrec.desc < 0) {
		sprintf(message, "isolated %s at rec # %d", type, num);
		mem_message(message);
		}
	}
}


void install_new_file(char *root, char *where)
{
char old_name[200], new_name[200];
int err;

strcpy(old_name, prepare_env_var(where));
strcat(old_name, root);
strcpy(new_name, old_name);
strcat(new_name, ".NEW");
err = remove(old_name);
if(err) {
	fprintf(stderr, "Error while deleting file %s\n", old_name);
	exit(ERREUR);
	}
err = rename(new_name, old_name);
if(err) {
	fprintf(stderr, "Error while renaming file %s to %s\n", 
		new_name, old_name);
	exit(ERREUR);
	}
return;
}


static void write_and_check(DIR_FILE *fich,int pos, int numrec, void *buffer)
{
static char text[200];
if( dir_write(fich,pos,numrec,buffer)!=0 ) {
	sprintf(text,"Error writing to file %s (probably disk space)",
		fich->filename);
	perror(text);
	exit(ERREUR);
	}
}


static void read_and_check(DIR_FILE *fich,int pos, int numrec, void *buffer)
{
static char text[200];
if( dir_read(fich,pos,numrec,buffer) == -1 ) {
	sprintf(text,"Error writing/reading to file %s (probably disk space)",
		fich->filename);
	perror(text);
	exit(ERREUR);
	}
}


static void close_and_check(DIR_FILE *fich)
{
static char text[200];
if( dir_close(fich)!=0 ) {
	sprintf(text,"Error closing file %s (probably disk space)",
		fich->filename);
	perror(text);
	exit(ERREUR);
	}
}


static void flush_and_check(DIR_FILE *fich)
{
static char text[200];
if( dir_flush(fich)!=0 ) {
	sprintf(text,"Error writing to file %s (probably disk space)",
		fich->filename);
	perror(text);
	exit(ERREUR);
	}
}


void make_lng(int point)
/* supprime trous, seq effacees et repetitions mais ne trie pas */
{
int i, num, count;
static struct rlng srlng;
memset(lng_list_aux,0,lenw*sizeof(int));
count= -1;
while(point) {
	readlng(point); point=plng->next;
	for(i=0; i<SUBINLNG; i++) {
		num=seqtable[plng->sub[i]];
		if(num && !testbit(lng_list_aux,num)) {
			bit1(lng_list_aux,num);
			++count;
			if(count>=SUBINLNG) {
				++lastlng2;
				srlng.next=lastlng2+1;
				write_and_check(klng2,lastlng2,1,&srlng);
				count=0;
				}
			srlng.sub[count]=num;
			}

		}
	}
while(++count<SUBINLNG) srlng.sub[count]=0;
srlng.next=0;
write_and_check(klng2,++lastlng2,1,&srlng);
}


static int int_sort(const void *i, const void *j)
{
return *(int *)i - *(int *)j;
}


void make_lng_sort(int point)
/* supprime trous, seq effacees et repetitions et trie */
{
int i, j, num, count, maxi, mini_0, maxi_0, mot, inlocal, totlocal;
static struct rlng srlng;
#define TRANCHE 500
static int local[TRANCHE + SUBINLNG];

totlocal = 0;
while(point) {
	readlng(point); point = plng->next;
	for(i = 0; i < SUBINLNG; i++) {
		if( (num = plng->sub[i]) != 0) {
			if( (num = seqtable[num]) != 0) local[totlocal++] = num;
			}
		}
	if(totlocal >= TRANCHE) break;
	}
inlocal = (point == 0);
if(inlocal) {
	if(totlocal > 1) qsort(local, totlocal, sizeof(int), int_sort);
	}
else	{
	mini_0 = maxi_0 = local[0];
	for(i = 1; i < totlocal; i++) {
		if(local[i] < mini_0)  mini_0 = local[i];
		if(local[i] > maxi_0)  maxi_0 = local[i];
		}
	mini_0 = (mini_0 - 1)/lmot;
	maxi_0 = (maxi_0 - 1)/lmot;
	i=mini_0; while(i <= maxi_0) lng_list_aux[i++] = 0;
	for(i=0; i<totlocal; i++) bit1(lng_list_aux, local[i] );
	while(point) {
		readlng(point); point = plng->next;
		for(i = 0; i < SUBINLNG; i++) {
			if( (num = plng->sub[i]) != 0) {
				if( (num = seqtable[num]) != 0) {
					mot = (num - 1)/lmot;
					if(mot < mini_0) {
						j = mot;
						while(j < mini_0)
							lng_list_aux[j++] = 0;
						mini_0 = mot;
						}
					else if(mot > maxi_0) {
						j = maxi_0 + 1;
						while(j <= mot)
							lng_list_aux[j++] = 0;
						maxi_0 = mot;
						}
					bit1(lng_list_aux, num );
					}
				}
			}
		}
	}

count= -1; /* ecriture de la nouvelle liste */
if(inlocal) {
	num = 0;
	for(i=0; i < totlocal; i++) {
		if(local[i] == num) continue; /* enlever les repetitions */
		num = local[i];
		++count;
		if(count >= SUBINLNG) {
			++lastlng2;
			srlng.next = lastlng2+1;
			write_and_check(klng2, lastlng2, 1, &srlng);
			count = 0;
			}
		srlng.sub[count] = num;
		}
	}
else	{
	num = mini_0 * lmot;
	maxi = (maxi_0 + 1) * lmot;
	while( (num = irbit(lng_list_aux, num, maxi) ) != 0) {
		++count;
		if(count >= SUBINLNG) {
			++lastlng2;
			srlng.next = lastlng2+1;
			write_and_check(klng2, lastlng2, 1, &srlng);
			count = 0;
			}
		srlng.sub[count] = num;
		}
	}
while( ++count < SUBINLNG) srlng.sub[count] = 0;
srlng.next = 0;
write_and_check(klng2, ++lastlng2, 1, &srlng);
#undef TRANCHE
}


void conv_shrt(int old)
{
while(old != 0) {
	dir_read(kshrt,old,1,pshrt);
	old=pshrt->next;
	++lastshrt2;
	if(old != 0) pshrt->next=lastshrt2+1;
	pshrt->val=seqtable[pshrt->val];
	write_and_check(kshrt2,lastshrt2,1,pshrt);
	}
}

void copy_shrt(int old)
{
while(old != 0) {
	dir_read(kshrt,old,1,pshrt);
	old=pshrt->next;
	++lastshrt2;
	if(old != 0) pshrt->next=lastshrt2+1;
	write_and_check(kshrt2,lastshrt2,1,pshrt);
	}
}

int comp_paires(const void *p1, const void *p2)
{
return ( ((struct paire *)p1)->new_val - ((struct paire *)p2)->new_val);
}

void *myalloc(size_t size)
{
void *retval;
retval=malloc(size);
if(retval==NULL)mem_message("Not enough memory");
return retval;
}


void make_desc_list(int next)
{
int ldes, pkey, suiv;

while(1) {
	read_and_check(kshrt2,next,1,pshrt);
	if(pshrt->next == 0) return;
	ldes=pshrt->next;
	pshrt->next=lastshrt2+1;
	write_and_check(kshrt2,next,1,pshrt);
	while(1) {
		dir_read(kshrt,ldes,1,pshrt);
		ldes=pshrt->next;
		dir_read(kshrt,pshrt->val,1,pshrt);
		pkey=pshrt->val;
		suiv=pshrt->next;
		dir_read(kan11,abs(pkey),1,&ksrec);
		next= ++lastshrt2;
		if(ksrec.desc <= 0) break;
		if(ldes == 0) {
			pshrt->val=ksrec.desc; pshrt->next=0;
			write_and_check(kshrt2,lastshrt2,1,pshrt);
			return;
			}
		pshrt->val=ksrec.desc; pshrt->next=lastshrt2+1;
		write_and_check(kshrt2,lastshrt2,1,pshrt);
		}
	pshrt->val=lastshrt2+1; pshrt->next=ldes;
	write_and_check(kshrt2,lastshrt2,1,pshrt);
	pshrt->val=pkey; pshrt->next=suiv;
	write_and_check(kshrt2,++lastshrt2,1,pshrt);
	ksrec.desc=lastshrt2;
	write_and_check(kan11,abs(pkey),1,&ksrec);
	make_desc_list(lastshrt2);
	}
}


int to_h_in_rsub; /* echange entre sortfile et sort_subseq_pointers */

int sort_subseq_pointers(const void *p1, const void *p2)
{
char *c1, *c2, *last, *deb1, *deb2;
int h1, h2;

c1 = *(char **)p1;
c2 = *(char **)p2;
deb1 = c1; 
deb2 = c2;
last = c1 + (compsize - 1);
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
return memcmp(rec1, rec2, compsize);
}


char **merge_table(char **table, int total, int fin_trie)
{
char **new_list;
int rank_1, rank_2, rank_new;

new_list = (char **)myalloc(total * sizeof(char *));
rank_1 = 0; rank_2 = fin_trie; rank_new = 0;
while(rank_new < total) {
	while( rank_1 >= fin_trie || memcmp(table[rank_1], 
			table[rank_2], compsize) > 0) {
		new_list[rank_new++] = table[rank_2++];
		if(rank_2 >= total) 
			break;
		}
	if(rank_new >= total) break;
	while( rank_2 >= total || memcmp(table[rank_1], 
			table[rank_2], compsize) < 0) {
		new_list[rank_new++] = table[rank_1++];
		if(rank_1 >= fin_trie) 
			break;
		}
	}
free(table);
return new_list;
}


int sortfile(DIR_FILE *kan, size_t rsize, int keysize, 
	char *outfname, char *env_var, DIR_FILE **outfkan, int *table)
{
char *array, croix[40], *debut, *fin, **deb_record, *previous;
int lcroix, insize, outsize, old_rank, new_rank, num, trie, fin_trie, num_lu;
size_t bigout;

printf("Sorting file %s\n",outfname); fflush(stdout);
lcroix = (keysize<40 ? keysize : 40);
num=0; while(num<lcroix) croix[num++]='x';
insize = read_first_rec(kan,NULL);
array = (char *)myalloc(insize * rsize);
/* lecture en memoire de tout le fichier a trier */
dir_set_normal(kan);
dir_resize_buff(kan, array, insize * rsize);
dir_read_buff(kan, 2, &num_lu);

if(kan == ksub && !( nbrf || swissprot) ) { 
/* cas du fichier SUBSEQ et banque nucleique:
on veut pour ordre des filles celui de leur citation dans la table des features
on remplace le champ h inutile de SUBSEQ par le pointeur vers feature
*/
/* offset vers 3 champs ds record rsub */
	int to_pext_in_rsub, to_plinf_in_rsub, *pext_field, *plinf_field, 
		*h_field;
	to_h_in_rsub = ( (char *)&(psub->h) - (char *)psub );
	to_pext_in_rsub = ( (char *)&(psub->pext) - (char *)psub );
	to_plinf_in_rsub = ( (char *)&(psub->plinf) - (char *)psub );
	debut = array; fin = array + (insize - 1) * rsize;
#ifdef unix
	dir_set_mmap(kshrt);
#endif
	while(debut < fin) {
		/* trouver champ pext */
		pext_field = (int *)(debut + to_pext_in_rsub); 
		if( *pext_field > 0 ) { /* pour une fille (pext > 0) */
			/* trouver champ plinf */
			plinf_field = (int *)(debut + to_plinf_in_rsub); 
			num = *plinf_field; /* lire champ plinf */
			/*lire ds SHORTL adresse de la location de cette fille*/
			readshrt(num); 
			/* trouver champ h */
			h_field = (int *)(debut + to_h_in_rsub); 
			/* stocker cette location ds champ h inutile */
			*h_field = pshrt->val; 
			}
		debut += rsize;
		}
#ifdef unix
	dir_set_normal(kshrt);
#endif
	}
deb_record = (char **)myalloc(sizeof(char *) * insize);
debut = array; fin = array + (insize - 1) * rsize;
outsize = 0;
while(debut < fin) {
	if(strncmp(debut, croix, lcroix) != 0) {
		deb_record[outsize++] = debut;
		}
	debut += rsize;
	}
compsize = keysize;
if( kan == ksub && !( nbrf || swissprot ) ) {
	qsort(deb_record, outsize, sizeof(char *), sort_subseq_pointers);
	}
else	{
/* calculer longueur de la partie deja triee */
	previous = deb_record[0];
	for(fin_trie = 1; fin_trie < outsize; fin_trie++) {
		trie = (memcmp(deb_record[fin_trie], previous, keysize) > 0);
		if( !trie ) break;
		previous = deb_record[fin_trie];
		}
	if(fin_trie > outsize) fin_trie = outsize;
	if(fin_trie < 2) fin_trie = 0;
/* trier la partie non triee du tableau deb_record */
	if( outsize - fin_trie > 0) qsort(deb_record + fin_trie, 
		outsize - fin_trie, sizeof(char *), sort_record_pointers);
/* trier tout deb_record en mergeant ses 2 parties */
	if( outsize > fin_trie && fin_trie > 0 ) 
		deb_record = merge_table(deb_record, outsize, fin_trie);
	}
bigout = (outsize + 1) * rsize;
if(bigout>1024000) bigout = 1024000;
*outfkan = dir_open(outfname, env_var, "w+", rsize, bigout);
if(*outfkan == NULL){
	char text[100];
	sprintf(text,"Cannot create new file %s",outfname);
	mem_message(text);
	}
write_first_rec(*outfkan, outsize + 1, outsize + 1);
memset(table, 0, (insize + 1) * sizeof(int) );
table[1]=1;
new_rank = 1;
for(num = 0; num < outsize; num++) {
	write_and_check(*outfkan, ++new_rank, 1, deb_record[num]);
	old_rank = (deb_record[num] - array) / rsize + 2;
	table[old_rank] = new_rank;
	}
flush_and_check(*outfkan);
free(deb_record);
	{char *tmp;
	tmp = (char *)myalloc(8*1024); 
	dir_resize_buff(kan, tmp, 8*1024);
	}
return outsize + 1;
}


int main(void)
{
int i,totspec, totkey, suiv, next, newns, tot12, is , isub, newloc, tloc,
	lastext2, old, newacc, point, newaut, newbib;
char buffer[100], where_new_files[40];

dir_acnucopen("RO");

tloc=read_first_rec(kloc,NULL);
#ifdef vms
strcpy(where_new_files, "BANK");
#else
strcpy(where_new_files, "acnuc");
#endif

/* Sorting SUBSEQ */
seqtable=(int *)myalloc((nseq+1)*sizeof(int));
newsub = sortfile(ksub, lrsub, L_MNEMO, "SUBSEQ.NEW", where_new_files, 
	&ksub2, seqtable);
seqtable[0]=0;
dir_close(ksub);

/* writing list of loci and unvalid seqs */
printf("Writing list of loci and unvalid seqs\n");
fflush(stdout);
#ifdef unix
dir_set_mmap(klng);
#endif
klng2 = dir_open("LONGL.NEW", where_new_files, "w+", lrlng, 1000*1024);
if(klng2==NULL)mem_message("Pb creating file LONGL.NEW");
lng_list_aux=(int *)myalloc(lenw*sizeof(int));
lastlng2=3;
dir_read(klng,2,1,plng);
if(plng->next==0) suiv=0;
else suiv=lastlng2+1;
i=0; while(i<SUBINLNG) {plng->sub[i]=seqtable[plng->sub[i]]; i++; }
next=plng->next;
plng->next=suiv;
write_and_check(klng2,2,1,plng);
i=1; while(i<SUBINLNG) plng->sub[i++]=0;
plng->next=0;
plng->sub[0]=1;
write_and_check(klng2,3,1,plng);
if(next != 0)make_lng(next);

/* sorting SMJYT.NEW */
old=read_first_rec(ksmj,NULL);
smjtable = (int *)myalloc((old+1)*sizeof(int));
smjtable[0]=0;
newns = sortfile(ksmj,lrsmj,20,"SMJYT.NEW", where_new_files, &kan11,smjtable);
ktxt2 = dir_open("TEXT.NEW", where_new_files, "w+", 60, 60*512);
if(ktxt2==NULL)mem_message("Pb creating file TEXT.NEW");
tot12=1;
for(is=2; is<=newns; is++) {
	dir_read(kan11,is,1,psmj);
	if(psmj->libel != 0) {
		dir_read(ktxt,psmj->libel,1,buffer);
		write_and_check(ktxt2,++tot12,1,buffer);
		psmj->libel=tot12;
		}
	if(psmj->plong != 0) {
		i=psmj->plong;
		psmj->plong=lastlng2+1;
		make_lng_sort(i);
		}
	write_and_check(kan11,is,1,psmj);
	}
close_and_check(kan11);
hashmn(psub->name); /* necessaire avant de fermer ksmj cause detection hashing algo */
dir_close(ksmj);

/* listes de sous-sequences, calcul hashing, extract.new */
if( !(nbrf || swissprot) ) 
	printf("Writing lists of sub-sequences and EXTRACT.NEW ");
printf("Computing sequence hashing\n");
fflush(stdout);
htable=(int *)myalloc((hsub+2)*sizeof(int));
loctable=(int *)myalloc((tloc+1)*sizeof(int));
i=0; while(i<=tloc) loctable[i++]=0;
i=1; while(i<=hsub+1) htable[i++]=0;
newloc=1;
if( !(nbrf || swissprot) ) {
	kan11 = dir_open("EXTRACT.NEW", where_new_files, "w+", lrext, 100*1024);
	if(kan11==NULL)mem_message("Pb creating file EXTRACT.NEW");
	lastext2=1;
#ifdef unix
	dir_set_mmap(kext);
#endif
	}
for(isub=2; isub<=newsub; isub++) {
	dir_read(ksub2,isub,1,psub);
	i=hashmn(psub->name); /* calcul hashing */
	psub->h=htable[i];
	htable[i]=isub;
	if(psub->pext<=0) { /* liste des filles */
		loctable[psub->plinf]= ++newloc;
		psub->plinf=newloc;
		if(psub->pext < 0) {
			i= -psub->pext;
			psub->pext= -(lastlng2+1);
			make_lng_sort(i);
			}
		}
	else	{ /* fabrication nouvel EXTRACT */
		old=psub->pext;
		psub->pext=lastext2+1;
		while(old != 0) {
			dir_read(kext,old,1,pext);
			old=pext->next;
			++lastext2;
			if(old != 0) pext->next= lastext2+1;
			else pext->next=0;
			pext->mere=seqtable[pext->mere];
			write_and_check(kan11,lastext2,1,pext);
			}
		}
	psub->type=smjtable[psub->type];
	write_and_check(ksub2,isub,1,psub);
	}
if( !( nbrf || swissprot) ) {
	write_first_rec(kan11,lastext2,0);
	close_and_check(kan11);
	dir_close(kext);
	}

/* SPECIES.NEW */
printf("Writing SPECIES.NEW\n");
fflush(stdout);
kan11 = dir_open("SPECIES.NEW", where_new_files, "w+", lrspec, 100*1024);
if(kan11==NULL)mem_message("Pb creating file SPECIES.NEW");
totspec=read_first_rec(kspec,NULL);
for(is=2; is<=totspec; is++ ) {
	dir_read(kspec,is,1,pspec);
	if(pspec->libel != 0) {
		dir_read(ktxt,pspec->libel,1,buffer);
		write_and_check(ktxt2,++tot12,1,buffer);
		pspec->libel=tot12;
		}
	i=pspec->plsub;
	if(i != 0) {
		pspec->plsub=lastlng2+1;
		make_lng_sort(i);
		}
	i=pspec->plhost;
	if(i != 0) {
		pspec->plhost=lastlng2+1;
		make_lng_sort(i);
		}
	pspec->desc= - pspec->desc;
	write_and_check(kan11,is,1,pspec);
	}
write_first_rec(kan11,totspec,0);
close_and_check(kan11);
dir_close(kspec);

/* KEYWORDS.NEW */
printf("Writing KEYWORDS.NEW\n");
fflush(stdout);
kan11 = dir_open("KEYWORDS.NEW", where_new_files, "w+", lrkey, 100*1024);
if(kan11==NULL)mem_message("Pb creating file KEYWORDS.NEW");
totkey=read_first_rec(kkey,NULL);
for(is=2; is<=totkey; is++ ) {
	dir_read(kkey,is,1,pkey);
	if(pkey->libel != 0) {
		dir_read(ktxt,pkey->libel,1,buffer);
		write_and_check(ktxt2,++tot12,1,buffer);
		pkey->libel=tot12;
		}
	i=pkey->plsub;
	if(i != 0) {
		pkey->plsub=lastlng2+1;
		make_lng_sort(i);
		}
	pkey->desc= - pkey->desc;
	write_and_check(kan11,is,1,pkey);
	}
write_first_rec(kan11,totkey,0);
close_and_check(kan11);
dir_close(kkey);
write_first_rec(klng2,lastlng2,0);
close_and_check(klng2);
free(lng_list_aux);
dir_close(klng);
write_first_rec(ktxt2,tot12,0);
close_and_check(ktxt2);
dir_close(ktxt);

/* sequence hashing data */
printf("Writing hashing data\n");
fflush(stdout);
kshrt2 = dir_open("SHORTL.NEW", where_new_files, "w+", lrshrt, 100*1024);
if(kshrt2==NULL)mem_message("Pb creating file SHORTL.NEW");
pshrt->val= -hsub; pshrt->next= -hkwsp;
write_and_check(kshrt2,2,1,pshrt);
old=2;
for(i=1; i<=hsub; i+=2) {
	pshrt->val=htable[i]; pshrt->next=htable[i+1];
	write_and_check(kshrt2,++old,1,pshrt);
	}
free((char *)htable);
#ifdef unix
dir_set_mmap(kshrt);
#endif

/* transfer of keyw and spec hashing data */
lastshrt2=2+(hsub+1)/2+hkwsp+1;
for(old=3+(hsub+1)/2; old<=lastshrt2; ++old) {
	dir_read(kshrt,old,1,pshrt);
	write_and_check(kshrt2,old,1,pshrt);
	}

/* short lists of keywords and info records */
printf("Short lists of keywords and info records\n");
fflush(stdout);
for(isub=2; isub<=newsub; isub++ ) {
	dir_read(ksub2,isub,1,psub);
	if(psub->plkey != 0) {
		i=psub->plkey;
		psub->plkey=lastshrt2+1;
		copy_shrt(i);
		}
	if(psub->pext > 0) {
		i=psub->plinf;
		psub->plinf=lastshrt2+1;
		copy_shrt(i);
		}
	write_and_check(ksub2,isub,1,psub);
	}
close_and_check(ksub2);

/* sorting ACCESS */
old=read_first_rec(kacc,NULL);
acctable=(int *)myalloc((old+1)*sizeof(int));
newacc = sortfile(kacc, kacc->record_length, ACC_LENGTH, "ACCESS.NEW",
		where_new_files, &kan11, acctable);
dir_close(kacc);
for(old=2; old<=newacc; old++) {
	dir_read(kan11,old,1,buffer);
	memcpy(&i, buffer + ACC_LENGTH, sizeof(int));
	if(i>0) {
		point=lastshrt2+1;
		memcpy(buffer + ACC_LENGTH, &point, sizeof(int));
		conv_shrt(i);
		}
	write_and_check(kan11,old,1,buffer);
	}
close_and_check(kan11);

/* sorting BIBLIO */
old=read_first_rec(kbib,NULL);
bibtable=(int *)myalloc((old+1)*sizeof(int));
newbib = sortfile(kbib, lrbib, 40, "BIBLIO.NEW", where_new_files, &kan11, bibtable);
dir_close(kbib);
printf("Writing LOCUS.NEW and lists of access#s and refers\n");
fflush(stdout);
kloc2 = dir_open("LOCUS.NEW", where_new_files, "w+", lrloc, 100*1024);
if(kloc2 == NULL)mem_message("Pb creating file LOCUS.NEW");
tabpairloc = (struct paire *)myalloc((newloc+1)*sizeof(struct paire));
suiv= 1;
i=2; while( i<=tloc ) {
	if(loctable[i] != 0) {
		suiv++;
		(tabpairloc+suiv)->old_val=i;
		(tabpairloc+suiv)->new_val=loctable[i];
		}
	i++;
	}
if(suiv != newloc)mem_message("erreur logique tabpairloc");
qsort(tabpairloc+2, newloc-1, sizeof(struct paire), comp_paires);
for(old=2; old<=newloc; old++) {
	dir_read(kloc,(tabpairloc+old)->old_val,1,ploc);
	ploc->bef=loctable[ploc->bef];
	ploc->next=loctable[ploc->next];
	ploc->molec=smjtable[ploc->molec];
	ploc->stat=smjtable[ploc->stat];
	ploc->org=smjtable[ploc->org];
	ploc->sub=seqtable[ploc->sub];
	if(ploc->placc != 0) {
		i=ploc->placc;
		ploc->placc=lastshrt2+1;
		while(i != 0) {
			dir_read(kshrt,i,1,pshrt);
			pshrt->val=acctable[pshrt->val];
			i=pshrt->next;
			++lastshrt2;
			if(i != 0) pshrt->next=lastshrt2+1;
			write_and_check(kshrt2,lastshrt2,1,pshrt);
			}
		}
	if(ploc->plref != 0) {
		i=ploc->plref;
		ploc->plref=lastshrt2+1;
		while(i != 0) {
			dir_read(kshrt,i,1,pshrt);
			pshrt->val=bibtable[pshrt->val];
			i=pshrt->next;
			++lastshrt2;
			if(i != 0) pshrt->next=lastshrt2+1;
			write_and_check(kshrt2,lastshrt2,1,pshrt);
			}
		}
	if( ( nbrf || swissprot ) && ploc->spec != 0 ) {
		i=ploc->spec;
		ploc->spec=lastshrt2+1;
		copy_shrt(i);
		}
	write_and_check(kloc2,old,1,ploc);
	}
write_first_rec(kloc2,newloc,0);
close_and_check(kloc2);
dir_close(kloc);
free((char *)tabpairloc);
free((char *)loctable);
free((char *)acctable);

/* sorting AUTHOR */
old=read_first_rec(kaut,NULL);
autable=(int *)myalloc((old+1)*sizeof(int));
newaut = sortfile(kaut, lraut, 20, "AUTHOR.NEW", where_new_files, 
	&kaut2, autable);
dir_close(kaut);
printf("Writing lists of seqs and authors for refers\n");
fflush(stdout);
for(old=2; old<=newbib; old++) {
	dir_read(kan11,old,1,pbib);
	if(pbib->plsub != 0) {
		i=pbib->plsub;
		pbib->plsub=lastshrt2+1;
		conv_shrt(i);
		}
	if(pbib->plaut != 0) {
		i=pbib->plaut;
		pbib->plaut=lastshrt2+1;
		while(i != 0) {
			dir_read(kshrt,i,1,pshrt);
			pshrt->val=autable[pshrt->val];
			i=pshrt->next;
			++lastshrt2;
			if(i != 0) pshrt->next=lastshrt2+1;
			write_and_check(kshrt2,lastshrt2,1,pshrt);
			}
		}
	pbib->j=smjtable[pbib->j];
	pbib->y=smjtable[pbib->y];
	write_and_check(kan11,old,1,pbib);
	}
close_and_check(kan11);
free((char *)smjtable);
free((char *)autable);
free((char *)seqtable);

/* writing AUTHOR.NEW */
printf("Writing lists of refers for authors\n");
fflush(stdout);
for(old=2; old<=newaut; old++) {
	dir_read(kaut2,old,1,paut);
	if(paut->plref != 0) {
		i=paut->plref;
		paut->plref=lastshrt2+1;
		while(i != 0) {
			dir_read(kshrt,i,1,pshrt);
			pshrt->val=bibtable[pshrt->val];
			i=pshrt->next;
			++lastshrt2;
			if(i != 0) pshrt->next=lastshrt2+1;
			write_and_check(kshrt2,lastshrt2,1,pshrt);
			}
		}
	write_and_check(kaut2,old,1,paut);
	}
close_and_check(kaut2);
free((char *)bibtable);

/* arbre des KEYWORDS */
printf("Writing tree structure of keywords\n");
fflush(stdout);
kan11 = dir_open("KEYWORDS.NEW", where_new_files, "r+", lrkey, totkey*lrkey);
if(kan11==NULL) mem_message("Pb re-opening KEYWORDS.NEW");
dir_read(kan11,2,1,&ksrec);
dir_read(kshrt,abs(ksrec.desc),1,pshrt);
write_and_check(kshrt2,++lastshrt2,1,pshrt);
ksrec.desc=lastshrt2;
write_and_check(kan11,2,1,&ksrec);
make_desc_list(lastshrt2);
check_desc_pointers(kan11, "keyword");
close_and_check(kan11);

/* arbre des SPECIES */
printf("Writing tree structure of species\n");
fflush(stdout);
kan11 = dir_open("SPECIES.NEW", where_new_files, "r+", lrspec, totspec*lrspec);
if(kan11==NULL) mem_message("Pb re-opening SPECIES.NEW");
dir_read(kan11,2,1,&ksrec);
dir_read(kshrt,abs(ksrec.desc),1,pshrt);
write_and_check(kshrt2,++lastshrt2,1,pshrt);
ksrec.desc=lastshrt2;
write_and_check(kan11,2,1,&ksrec);
make_desc_list(lastshrt2);
check_desc_pointers(kan11, "species");
close_and_check(kan11);

write_first_rec(kshrt2,lastshrt2,0);
close_and_check(kshrt2);
dir_close(kshrt);

/* renommage des fichiers */
printf("Replacing old index files by new ones\n");
install_new_file("ACCESS", where_new_files);
install_new_file("AUTHOR", where_new_files);
install_new_file("BIBLIO", where_new_files);
if( ! (nbrf || swissprot) ) install_new_file("EXTRACT", where_new_files);
install_new_file("KEYWORDS", where_new_files);
install_new_file("LOCUS", where_new_files);
install_new_file("LONGL", where_new_files);
install_new_file("SHORTL", where_new_files);
install_new_file("SMJYT", where_new_files);
install_new_file("SPECIES", where_new_files);
install_new_file("SUBSEQ", where_new_files);
install_new_file("TEXT", where_new_files);

#ifdef unix
write_quick_meres();
#endif

printf("Normal end\n");
return 0;
} /* end of main */
