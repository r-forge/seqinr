#include "dir_acnuc.h"
#include <ctype.h>

/* included functions */
void dir_acnucflush(void);
void dir_writeerr(DIR_FILE *fich, int recnum);
void writeacc(int recnum);
void write_first_rec(DIR_FILE *fp, int total, int endsort);
void suphsh(int numrec, DIR_FILE *kan);
void addhsh(int recnum, DIR_FILE *kan);
void delseq(int nsub);
void write_sorted_part(DIR_FILE *fp);
int crespecies(char *ascend, char *nom);
int crekeyword(char *ascend, char *nom);
void cre_new_division(char *name);
void update_div_size(int numdiv);
unsigned get_current_div_size(int numdiv);
int read_addr_loc_qualif(long faddr, int div, 
	char *location, int maxlocat, char *type, char *qualif, int maxqualif);
char *get_qualif_value(char *p);
int chgnam(int isub, char *qualif, int maxqual);
int loadtypes( int **type_ranks, char ***type_codes, char ***type_abbrevs);
#ifdef unix
/* Appeler cette fonction une fois au debut du programme
puis dans un endroit ou on peut aller vers un fin propre.
A cet endroit, envoyer vers fin propre si check_term rend TRUE
*/
int check_term(void);
void write_quick_meres(void);
#endif


/* prototypes of used functions */
void gcgors(char *type_fic, int div, int stop_if_error);

/* extern data used */
extern FILE *divannot[], *divseq[];
extern int use_div_sizes, max_divisions;
extern int annotopened[], seqopened[];
extern unsigned div_offset[];
extern char *gcgname[];


extern DIR_FILE *kinf, *knuc;

void dir_acnucflush(void)
{
dir_flush(kshrt);
dir_flush(ksub);
dir_flush(kloc);
if( !(nbrf || swissprot) )dir_flush(kext);
if(nbmrfa==lmot) {
	dir_flush(kinf);
	dir_flush(knuc);
	}
if(ksmj != NULL) {
	dir_flush(ksmj);
	dir_flush(kkey);
	dir_flush(kspec);
	dir_flush(kbib);
	dir_flush(kacc);
	dir_flush(ktxt);
	dir_flush(kaut);
	dir_flush(klng);
	}
return;
}


void dir_writeerr(DIR_FILE *fich, int recnum)
{
char text[150];
sprintf(text,"Error while writing in %s at record # %d",fich->filename,recnum);
perror(text);
exit(ERREUR);
}

extern int acc_length_on_disk; /* taille chaine ACCESSION sur le disque */

void writeacc(int recnum)
{
/* ici structure en memoire ne colle pas celle sur le disque
*/
static char point[ACC_LENGTH + 2*sizeof(int)]; /* trop grand expres */
memset(point, 0, sizeof(point) );
memcpy(point, pacc->name, acc_length_on_disk);
memcpy(point + acc_length_on_disk, &(pacc->plsub), sizeof(int));
if(dir_write(kacc,recnum,1,point)) dir_writeerr(kacc,recnum);
}



void write_first_rec(DIR_FILE *fp, int total, int endsort)
/* ecriture du 1er record
total: nbre total de records
endsort: fin de la partie triee: pas trie si 0
				SORTED si endsort==total
				1/2SOR si endsort<total
*/
{
static int buffer[64];
char *record = (char *)buffer;

memset(record,0,fp->record_length);
buffer[0]=total;
if(endsort == total)
	memcpy(record+sizeof(int),"SORTED",6);
else if(endsort != 0) {
	memcpy(record+sizeof(int),"1/2SOR",6);
	memcpy(record+sizeof(int)+6, (char *)&endsort,  sizeof(int) );
	}
if( dir_write(fp,1,1,record) )exit(ERREUR);
}



void suphsh(int numrec, DIR_FILE *kan)
{
int next, point, suiv, prec, h;
if(kan==ksub) {
	readsub(numrec);
	next=psub->h;
	h=hashmn((char *)psub);
	point=hoffst+(h+1)/2;
	}
else 	{
	dir_read(kan,numrec,1,pspec);
	next=pspec->h;
	h=hasnum((char *)pspec);
	if(kan==kkey) {
		point=hoffst+(hsub+1)/2;
		}
	else	{
		point=hoffst+(hsub+1)/2+(hkwsp+1)/2;
		}
	point += (h+1)/2;
	}
readshrt(point);
if(h%2) prec=pshrt->val;
else prec=pshrt->next;
if(prec==numrec) {
	if(h%2) pshrt->val=next;
	else pshrt->next=next;
	writeshrt(point);
	if(kan==ksub) {
		psub->h=0;
		writesub(numrec);
		}
	else 	{
		pspec->h=0;
		if(dir_write(kan,numrec,1,pspec))
			dir_writeerr(kan,numrec);
		}
	}
else if (prec!=0){
   	while(TRUE) {
		if(kan==ksub) {
			readsub(prec);
			suiv=psub->h;
			}
		else	{
			dir_read(kan,prec,1,pspec);
			suiv=pspec->h;
			}
		if(suiv==numrec) break;
		if(suiv==0)return;
		prec=suiv;
		}
	if(kan==ksub) {
		psub->h=next;
		writesub(prec);
		}
	else	{
		pspec->h=next;
		if(dir_write(kan,prec,1,pspec))
			dir_writeerr(kan,prec);
		}
	}
} /* end of suphsh */

void addhsh(int recnum, DIR_FILE *kan)
{
int point, h, next;
if(kan==ksub) {
	readsub(recnum);
	h=hashmn(psub->name);
	point=hoffst+(h+1)/2;
	}
else	{
	if(dir_read(kan,recnum,1,pspec) != 1) dir_readerr(kan,recnum);
	h=hasnum(pspec->name);
	point=hoffst+(hsub+1)/2+(h+1)/2;
	if(kan==kspec) point += (hkwsp+1)/2;
	}
readshrt(point);
if(h%2)	{
	next=pshrt->val;
	pshrt->val=recnum;
	}
else	{
	next=pshrt->next;
	pshrt->next=recnum;
	}
writeshrt(point);
	if(kan==ksub) {
	psub->h=next;
	writesub(recnum);
	}
else	{
	pspec->h=next;
	if(dir_write(kan,recnum,1,pspec))dir_writeerr(kan,recnum);
	}
} /* end of addhsh */


/* global for delseq+delsubseq */
static const char croix[]="xxxxxxxxxxxxxxxx";

static void delsubseq(int isub, int mere)
/* pour traiter la partie subseq de l'effacement d'une sequence
isub=num de la seq
mere=num de la mere qui est aussi effacee si on efface une fille a cause de 
     sa mere, =0 si on efface une fille ou une mere
*/
{
struct rsub rsub2, *psub2= &rsub2;
int valeur, point;
dir_read(ksub,isub,1,psub2);
addlng(3,isub);
suphsh(isub,ksub);
memset(psub,0,lrsub);
memcpy(psub->name,croix,16);
writesub(isub);
if(psub2->type) mdlng(ksmj,psub2->type,-1,isub,NULL);
point=psub2->plkey;
while(point) {
	readshrt(point);
	valeur = pshrt->val;
	point = pshrt->next;
	mdlng(kkey,valeur,-2,isub,NULL);
	}
point=psub2->pext;
if(point>0) {
	int autremere=1, premere;
	char name2[17], *pos, *fin=name2+15;
	memcpy(name2,psub2->name,16);
	name2[16]=0;
	pos=strchr(name2,'.');
	if(pos != NULL) {
		while(pos<=fin) *(pos++)=' ';
		}
	premere=0;
	while(point) {
		readext(point);
		point=pext->next;
		if(pext->mere==premere)continue;
		premere= pext->mere;
		if(premere!=mere) mdlng(ksub,premere,-3,isub,NULL);
		if(autremere) {
			readsub(premere);
			autremere= strncmp(psub->name,name2,16);
			}
		}
	if(autremere) {
		autremere=isenum(name2);
		if(autremere!=mere) mdlng(ksub,autremere,-3,isub,NULL);
		}
	}	

} /* end of delsubseq */


void delseq(int nsub)
{
struct rloc rloc2;
struct rlng rlng2;
int locus, list_fi, i, point;
readsub(nsub);
if(psub->pext > 0) goto lab1000;
locus=psub->plinf;
list_fi= - psub->pext;
readloc(locus);
if(ploc->molec) mdlng(ksmj,ploc->molec,-1,nsub,NULL);
if(ploc->stat) mdlng(ksmj,ploc->stat,-1,nsub,NULL);
if(ploc->org) mdlng(ksmj,ploc->org,-1,nsub,NULL);
suplng(2,nsub);
point=ploc->placc;
while(point) {
	readshrt(point);
	point=pshrt->next;
	mdshrt(kacc,pshrt->val,-1,nsub,NULL);
	}
if( ( nbrf || swissprot ) && ploc->spec != 0) {
	point=ploc->spec;
	while(point) {
		readshrt(point);
		point=pshrt->next;
		mdlng(kspec,pshrt->val,-2,nsub,NULL);
		}
	}
else if (ploc->spec!=0) {
	mdlng(kspec,ploc->spec,-2,nsub,NULL);
	}
if(ploc->host) {
	mdlng(kspec,ploc->host,-6,nsub,NULL);
	}
point=ploc->plref;
while(point) {
	int nbib;
	readshrt(point);
	point=pshrt->next;
	nbib=pshrt->val;
	mdshrt(kbib,nbib,-1,nsub,NULL);
	readbib(nbib);
	if(pbib->j) mdlng(ksmj,pbib->j,-1,nsub,NULL);
	if(pbib->y) mdlng(ksmj,pbib->y,-1,nsub,NULL);
	}
if(ploc->bef) {
	dir_read(kloc,ploc->bef,1,&rloc2);
	rloc2.next=0;
	dir_write(kloc,ploc->bef,1,&rloc2);
	}
if(ploc->next) {
	dir_read(kloc,ploc->next,1,&rloc2);
	rloc2.bef=0;
	dir_write(kloc,ploc->next,1,&rloc2);
	}
memset(&rloc2,0,lrloc);
memcpy(rloc2.date,croix,16);
dir_write(kloc,locus,1,&rloc2);

while(list_fi) {
	dir_read(klng,list_fi,1,&rlng2);
	list_fi=rlng2.next;
	for(i=0; i<SUBINLNG; i++ ) 
		if(rlng2.sub[i]) delsubseq(rlng2.sub[i],nsub);
	}
lab1000:
delsubseq(nsub,0);
point=read_first_rec(ksub,NULL);
write_first_rec(ksub,point,0);
} /* end of delseq */


/* a appeler en fin de programme qui allonge kaut, kbib, kacc ou ksmj
pour chaque fichier allonge' */

extern int fcode_nsortd[];
void write_sorted_part(DIR_FILE *fp)
{
int total, numf;
if (fp==kaut) { numf=0; }
else if (fp==kbib) { numf=1; }
else if (fp==kacc)  { numf=2; }
else if (fp==ksmj)  { numf=3; }
else return;
total=read_first_rec(fp,NULL);
write_first_rec(fp,total,fcode_nsortd[numf]);
}


static int cre_key_spec_node(DIR_FILE *kan, int pdes, char *nom)
{
unsigned totshrt, totspec, next;

totspec=read_first_rec(kan,NULL);
totspec++;
totshrt=read_first_rec(kshrt,NULL);
totshrt++;
write_first_rec(kshrt,totshrt+1,0);
memset(pspec,0,lrspec);
padtosize(pspec->name,nom,WIDTH_KS);
pspec->desc=totshrt;
dir_write(kan,totspec,1,pspec);
write_first_rec(kan,totspec,0);
addhsh(totspec,kan);
readshrt(pdes);
next=pshrt->next;
pshrt->next=totshrt+1;
writeshrt(pdes);
pshrt->val=totspec;
pshrt->next=0;
writeshrt(totshrt);
pshrt->val=totshrt;
pshrt->next=next;
writeshrt(totshrt+1);
return totspec;
}

int crespecies(char *ascend, char *nom)
/*
returns number of species nom if exists
if does not exist, creates it under species ascend 
(created as new root if does not exist)
ascend=NULL to create nom as new root
*/
{
int numascend, i, totspec, pdes;
char upnom[WIDTH_KS + 1];

padtosize(upnom,nom,WIDTH_KS);
for(i=0;i<WIDTH_KS;i++) upnom[i]=toupper(upnom[i]);
totspec=iknum(upnom,kspec);
while(totspec && pspec->syno >0) {
	totspec=pspec->syno;
	readspec(totspec);
	}
if(totspec) return totspec;
if(ascend==NULL) {
	numascend=2;
	}
else	{
	char upascend[WIDTH_KS + 1];
	padtosize(upascend,ascend,WIDTH_KS);
	for(i=0;i<WIDTH_KS;i++) upascend[i]=toupper(upascend[i]);
	numascend=iknum(upascend,kspec);
	if(numascend==0) 
		numascend=crespecies(NULL,ascend);
	else	{
		readspec(numascend);
		while(pspec->syno >0) {
			numascend=pspec->syno;
			readspec(numascend);
			}
		}
	}
readspec(numascend);
pdes=pspec->desc;
readshrt(pdes);
while(pshrt->next!=0) {
	pdes=pshrt->next;
	readshrt(pdes);
	}
return cre_key_spec_node(kspec,pdes,upnom);
}


static int *liste_racines, l_liste_racines;
static int end_of_root_chain;

static void init_liste_racines(void)
{
int totkey, next, val;
totkey=read_first_rec(kkey,NULL);
l_liste_racines = (totkey-1)/lmot + 1 + 50;
liste_racines = (int *)calloc(l_liste_racines, sizeof(int));
if(liste_racines == NULL) { /* brutal mais devrait etre suffisant */
	fprintf(stderr, "Error: no memory for init_liste_racines\n");
	exit(ERREUR);
	}
readkey(2);
next=pkey->desc;
readshrt(next);
while(pshrt->next) {
	next=pshrt->next;
	readshrt(next);
	val=pshrt->val;
	readshrt(val);
	bit1(liste_racines, abs(pshrt->val));
	readshrt(next);
	}
end_of_root_chain = next;
}

static void mem_new_root(int num)
{
int totkey, newmots;
totkey=read_first_rec(kkey,NULL);
newmots = (totkey-1)/lmot + 1;
if(newmots > l_liste_racines) {
	int *point;
	newmots += 50;
	point = (int *)calloc(newmots, sizeof(int));
/* si pas de memoire on ne memorise plus, les autres mots-cles a mettre
en descendant de quelque chose resteront racine */
	if(point == NULL) return;
	memcpy(point, liste_racines, l_liste_racines * sizeof(int));
	free(liste_racines);
	liste_racines = point;
	l_liste_racines = newmots;
	}
bit1(liste_racines, num);
}

static int is_keyw_root(int num)
{
if( num > l_liste_racines * lmot )
	return FALSE;
return testbit(liste_racines,num);
}


static void unsetroot(int numkey)
{
int val, next, pdes, pre;

readkey(numkey);
pdes=pkey->desc;
readkey(2);
pre=pkey->desc;
readshrt(pre);
next=pshrt->next;
while(next) {
	readshrt(next);
	if(pshrt->val==pdes) break;
	pre=next;
	next=pshrt->next;
	}
readshrt(pre);
val=pshrt->val;
next=pshrt->next;
readshrt(next);
pshrt->val=val;
writeshrt(pre);
end_of_root_chain=pre;
bit0(liste_racines,numkey);
}


static int fin_desc_list(int keynum)
/* rendre recnum ds SHORTL de la fin de la liste des descendants de keynum
version avec memorisation de plusieurs fins
*/
{
#define TOT_MEM_ASCEND 100
static int fin_list_from_ascend[TOT_MEM_ASCEND];
static int list_num_ascend[TOT_MEM_ASCEND];
static int tot_list = 0;
int num_list, point;

for(num_list = 0; num_list < tot_list; num_list++) 
	if(list_num_ascend[num_list] == keynum) break;
if(num_list >= tot_list) { /* ascendant non deja memorise */
	/* parcours complet des descendants jusqu'a la fin */
	readkey(keynum);
	point = pkey->desc;
	readshrt(point);
	while(pshrt->next) {
		point = pshrt->next;
		readshrt(point);
		}
	if(tot_list < TOT_MEM_ASCEND) { /* possible de memoriser plus */
		/* memorisation pour nouveau mot-cle */
		list_num_ascend[tot_list] = keynum;
		fin_list_from_ascend[tot_list] = point;
		tot_list++;
		}
	}
else	{
	/* utiliser et mettre a jour fin deja connue pour ce mot-cle */
	point = fin_list_from_ascend[num_list];	
	readshrt(point);
	while(pshrt->next) {
		point = pshrt->next;
		readshrt(point);
		}
	fin_list_from_ascend[num_list] = point;
	}
return point;
#undef TOT_MEM_ASCEND
}


static void setkeywbranch(int numascend, int num)
/* do not use for making a root 
mettre keyw num comme descendant de keyw numascend
*/
{
int point, pdes;
unsigned totshrt;

/* recherche de la fin de la liste des descendants de numascend */
point = fin_desc_list(numascend);
readkey(num);
pdes = pkey->desc;
totshrt = read_first_rec(kshrt,NULL);
totshrt++;
write_first_rec(kshrt,totshrt,0);
readshrt(point);
pshrt->next=totshrt;
writeshrt(point);
pshrt->val=pdes;
pshrt->next=0;
writeshrt(totshrt);
}


int crekeyword(char *ascend, char *nom)
{
int numascend, i, totkey, pdes;
char upnom[WIDTH_KS + 1];
static int first=1;

if(first) {
	first=0;
	init_liste_racines();
	}

padtosize(upnom,nom,WIDTH_KS);
for(i=0;i<WIDTH_KS;i++) upnom[i]=toupper(upnom[i]);
totkey=iknum(upnom,kkey);
while(totkey && pkey->syno >0) {
	totkey=pkey->syno;
	readkey(totkey);
	}
if(ascend==NULL) {
	numascend=2;
	}
else	{
	char upascend[WIDTH_KS + 1];
	padtosize(upascend,ascend,WIDTH_KS);
	for(i=0;i<WIDTH_KS;i++) upascend[i]=toupper(upascend[i]);
	numascend=iknum(upascend,kkey);
	if(numascend==0) 
		numascend=crekeyword(NULL,ascend);
	else	{
		readkey(numascend);
		while(pkey->syno >0) {
			numascend=pkey->syno;
			readkey(numascend);
			}
		}
	}
if(totkey) { /* nom existait deja */
	if(numascend!=2) { /* nom ne doit pas etre une racine */
		if( is_keyw_root(totkey) ) {
			/* enlever nom comme racine */
			unsetroot(totkey);
			/* mettre nom en descendant de ascend */
			setkeywbranch(numascend,totkey);
			}
		}
	return totkey;
	}
/* nom est un nouveau mot-cle */
if(numascend!=2) { /* nom ne doit pas etre une racine */
	/* recherche de la fin de la liste des descendants de ascend */
	pdes = fin_desc_list(numascend);
	}
else	{ /* nom doit etre une racine */
	pdes=end_of_root_chain;
	readshrt(pdes);
	while(pshrt->next!=0) {
		pdes=pshrt->next;
		readshrt(pdes);
		}
	end_of_root_chain=pdes;
	}
totkey= cre_key_spec_node(kkey,pdes,upnom);
if(numascend==2) mem_new_root(totkey);
return totkey;
}


void cre_new_division(char *name)
{
int totsmj, finsorted, tottxt, l;
unsigned taille, taille_previous;

if(divisions + 1 >= max_divisions) {
	fprintf(stderr,"too many divisions\n");
	exit(ERREUR);
	}

if( use_div_sizes && divisions>=0) {
	update_div_size(divisions);
	taille_previous = get_current_div_size(divisions);
	}
divisions++;

gcgname[divisions]=(char *)malloc(16);
strcpy(gcgname[divisions],name);
annotopened[divisions]=seqopened[divisions]=0;
if( use_div_sizes ) taille = get_current_div_size(divisions);

memset(psmj,0,lrsmj);
sprintf(psmj->name,"06FLT%-15.15s",name);
if(!flat_format) memcpy(psmj->name+2,"GCG",3);
if(big_annots)
	sprintf(ptxt,"rank:%d",divisions);
else
	sprintf(ptxt,"rank:%d size:%u",divisions,taille);
l=strlen(ptxt); memset(ptxt+l,' ',60-l);
tottxt=read_first_rec(ktxt,NULL)+1;
writetxt(tottxt);
write_first_rec(ktxt,tottxt,0);
psmj->libel=tottxt;
totsmj=read_first_rec(ksmj,&finsorted)+1;
writesmj(totsmj);
write_first_rec(ksmj,totsmj,finsorted);
if( use_div_sizes ) {
	if(divisions>=1) { 
		div_offset[divisions]= div_offset[divisions-1] + 
					(taille_previous/1024+1)*1024;
		}
	else	div_offset[divisions]=0;
	}
} /* end of cre_new_division */


void update_div_size(int numdiv)
/* calcule la taille actuelle d'une division, l'ecrit dans TEXT */
{
int num, totsmj, l;
unsigned taille;
char aux[21];

if( big_annots ) return;
/* calcul de la taille de la division */
taille=get_current_div_size(numdiv);

/* recherche de l'element ds SMJYT qui decrit la division */
sprintf(aux, "06FLT%-15.15s", gcgname[numdiv]);
if(!flat_format) memcpy(aux + 2, "GCG", 3);
totsmj=read_first_rec(ksmj,NULL);
for(num=totsmj; num>=2; num--) {
	readsmj(num);
	if(strncmp(aux,psmj->name,20)==0) break;
	}
if(num<2) {
	fprintf(stderr,"Last div not found.\n");
	exit(ERREUR);
	}

/* ecriture de la taille en libelle de cet element */
sprintf(ptxt,"rank:%d size:%u",numdiv,taille);
l=strlen(ptxt); memset(ptxt+l,' ',60-l);
writetxt(psmj->libel);
}


unsigned get_current_div_size(int numdiv)
{
unsigned taille, taille2;

if( big_annots ) return 0;
/* calcul de la taille de la division */
gcgors("inf",numdiv,FALSE);
if(annotopened[numdiv] == 0) return 0;
fseek(divannot[numdiv],0,SEEK_END);
taille = (unsigned)ftell(divannot[numdiv]);
taille2 = 0;
fclose(divannot[numdiv]);
annotopened[numdiv] = FALSE;
if(!flat_format) {
	gcgors("nuc",numdiv,TRUE);
	fseek(divseq[numdiv],0,SEEK_END);
	taille2 = (unsigned)ftell(divseq[numdiv]);
	fclose(divseq[numdiv]);
	seqopened[numdiv] = FALSE;
	}
if(taille2 > taille) taille = taille2;
return taille;
}


int read_addr_loc_qualif(long faddr, int div, 
	char *location, int maxlocat, char *type, char *qualif, int maxqualif)
/* lecture de la location a l'addr faddr + div et de ses qualifiers 
location != NULL: memoire pour mettre la location dedans
location == NULL: location ignoree
maxlocat: place memoire max disponible
type: (si != NULL) memoire pour mettre le type dedans
qualif != NULL: memoire pour metre les qualifiers
qualif == NULL: on ignore les qualifiers
maxqualif: place memoire max disponible
valeur rendue: FALSE si ok; TRUE si erreur pas assez de memoire
*/
{
int i, lloc=0, l, erreur, next_feat, lu;
char *ps, *p;

erreur = FALSE; next_feat = FALSE;
maxlocat--;
read_annots(faddr, div);
location[0] = 0;
if(type != NULL) { /* memorisation du type */
	l = 0; ps = pinfo->line + 4; 
	while(*(++ps) != ' ') type[l++] = toupper(*ps); type[l] = 0;
	}
while( (ps = strchr(pinfo->line,'/')) == NULL ) {
	if(location != NULL) {
		i = strlen(pinfo->line);
		while( pinfo->line[i - 1] == ' ') i--;
		i -= 21;
		if(lloc+i >= maxlocat) {
			erreur = TRUE;
			i = maxlocat - lloc;
			}
		if(i>0) {
			memcpy(location+lloc, pinfo->line + 21, i);
			lloc += i; location[lloc] = 0;
			}
		}
	next_annots(NULL);
	ps = NULL;
	if(strcmptrail(pinfo->line,20,NULL,0) != 0 && 
		strcmptrail(pinfo->line,20,"FT",2) !=0 ) {
		next_feat = TRUE;
		break;
		}
	}
if(location != NULL && ps != NULL) {
	i = ps - pinfo->line + 1;
	if( strcmptrail(pinfo->line+2, i-3, NULL, 0) != 0 && i > 22 ) {
		if(lloc+i-22 >= maxlocat) {
			erreur = TRUE;
			i = maxlocat - lloc + 22 - 1;
			}
		if(i > 22) {
			memcpy(location + lloc, pinfo->line + 21, i - 22);
			lloc += (i-22); location[lloc] = 0;
			}
		}
	}
if(location != NULL) {
	compact(location);
	majuscules(location);
	}

/* lecture des qualifiers */
if(qualif == NULL) return erreur;
qualif[0] = 0;
if (next_feat) return erreur;
maxqualif--;
l = 0;
do	{
	lu = strlen(ps);
	if( l + lu + 1>= maxqualif ) {
		erreur = TRUE;
		lu = maxqualif - l - 1;
		}
	if(lu > 1) {
		memcpy(qualif + l, ps, lu);
		l += lu + 1;
		qualif[l - 1] = ' '; qualif[l] = 0;
		}
	next_annots(NULL);
	ps = pinfo->line + 21;
	}
while (strcmptrail(pinfo->line, 21, NULL, 0) == 0 || 
			strcmptrail(pinfo->line, 21, "FT", 2) == 0);
majuscules(qualif);
return erreur;
}


char *get_qualif_value(char *p)
/* extraire la valeur apres le = d'un qualifier /mot=
p: pointe sur /
les " " ou ' ' encadrantes sont enlevees
rend la valeur dans un static char alloue ici
	ou NULL si la valeur est vide: /mot=""
*/
{
static char value[100];
char *q;
int l, use_quotes;

p = strchr(p, '=') + 1;
use_quotes = ( *p == '"' || *p == '\'' );
if( use_quotes ) p++;
q = strchr(p, '/');
if(q == NULL) q = p + strlen(p);
q--;
while( *q == ' ') q--;
if( use_quotes && *q == *(p - 1) ) q--;
l = q - p + 1;
if( l <= 0 ) return NULL;
if( l >= sizeof(value) ) l = sizeof(value) - 1;
memcpy(value, p, l); value[l] = 0;
return value;
}


int chgnam(int isub, char *qualif, int maxqual)
/* changer le nom d'une fille en fonction de ses qualifiers /GENE= etc...
rend TRUE si nom change', FALSE sinon
isub: numero de la fille
qualif: memoire pour charger les qualifiers dedans
maxqual: taille de cette place memoire
*/
{
char *p, *deb_qual, *fin_qual;
int lqual, l, div;
long faddr;
static char newname[L_MNEMO + 1];

readsub(isub);
if(psub->length == 0 || psub->pext <= 0) return FALSE;
p = strchr(psub->name, '.');
if(p == NULL || p >= psub->name + L_MNEMO ) return FALSE;
seq_to_annots(isub, &faddr, &div);
read_addr_loc_qualif(faddr, div, NULL, 0, NULL, qualif, maxqual);
lqual = strlen(qualif);
if(lqual == 0) return FALSE;
deb_qual = strstr(qualif, "/GENE=");
if(deb_qual == NULL) deb_qual = strstr(qualif, "/NOMGEN=");
if(deb_qual == NULL) deb_qual = strstr(qualif, "/STANDARD_NAME=");
if(deb_qual == NULL) return FALSE;
deb_qual = get_qualif_value(deb_qual);
if(deb_qual == NULL) return FALSE;
fin_qual = deb_qual + strlen(deb_qual) - 1;
if( (p - psub->name) + (fin_qual - deb_qual) + 2 > L_MNEMO )
	return FALSE;
l = p - psub->name + 1;
memcpy(newname, psub->name, l);
memcpy(newname + l, deb_qual, (fin_qual - deb_qual) + 1);
newname[l + (fin_qual - deb_qual) + 1] = 0;
/* replace internal spaces by '-' */
while ( ( p = strchr(newname, ' ') ) != NULL ) *p = '-';
/* remove @ characters */
while ( ( p = strchr(newname, '@') ) != NULL && p >= newname + l ) *p = ' ';
compact(newname);
padtosize(newname, newname, L_MNEMO);
if(isenum(newname) != 0) return FALSE;
suphsh(isub, ksub);
readsub(isub);
memcpy(psub->name, newname, L_MNEMO);
writesub(isub);
addhsh(isub, ksub);
return TRUE;
}


int loadtypes( int **type_ranks, char ***type_codes, char ***type_abbrevs)
/* lecture de tous les types (sauf LOCUS/ID) connus et chargement de leur
rang, nom, et abbreviation.
Tout en alloc dynamique de memoire.
Retour: nbre de types trouves (from 1), ou -1 si erreur de memoire
*/
{
#define FIXTYPES 6
#define MAXLOADTYPES 20
int num, totsmj, lcode, oldformat = FALSE;
int ntypes = -1;
char *p;
const int fixtypes = FIXTYPES;
static char fixtypec[FIXTYPES][10] = { 
	"CDS","RRNA","TRNA","SCRNA","SNRNA","MISC_RNA"
	};
static char fixtypeab[FIXTYPES][3] = { 
	"PE","RR","TR","SC","SN","RN"
	};
char aux[21];
int currmaxtypes = MAXLOADTYPES;

totsmj = read_first_rec(ksmj, NULL);

encore:
*type_ranks = (int *)malloc(currmaxtypes * sizeof(int));
if( *type_ranks == NULL ) return -1;
*type_codes = (char **)malloc(currmaxtypes * sizeof(char *));
if( *type_codes == NULL ) return -1;
*type_abbrevs = (char **)malloc(currmaxtypes * sizeof(char *));
if( *type_abbrevs == NULL ) return -1;

for(num = 2; num <= totsmj; num++) {
	readsmj(num);
	if(strncmp(psmj->name, "04", 2) != 0) continue;
	if(strcmptrail(psmj->name, 20, "04LOCUS", 7) == 0) continue;
	if(strcmptrail(psmj->name, 20, "04ID", 4) == 0) continue;
/* check that all type libels begin with a dot */
	if(psmj->libel != 0) readtxt(psmj->libel);
	if(psmj->libel == 0 || *ptxt != '.') {
		oldformat = TRUE;
		break;
		}
	ntypes++;
	if(ntypes >= MAXLOADTYPES) {
		currmaxtypes += 10;
		goto encore;
		}
	(*type_ranks)[ntypes] = num;
	lcode = sizeof(psmj->name) - 2; 
	while(psmj->name[2 + lcode - 1] == ' ') lcode--;
	(*type_codes)[ntypes] = malloc(lcode + 1);
	if( (*type_codes)[ntypes] == NULL) return -1;
	memcpy( (*type_codes)[ntypes], psmj->name + 2, lcode);
	(*type_codes)[ntypes][lcode] = 0;
	(*type_abbrevs)[ntypes] = malloc(3);
	if( (*type_abbrevs)[ntypes] == NULL) return -1;
	memcpy( (*type_abbrevs)[ntypes], ptxt + 1, 2);
	(*type_abbrevs)[ntypes][2] = 0;
	p = (*type_abbrevs)[ntypes] - 1; while( *(++p) != 0) *p = toupper(*p);
	}
if( !oldformat ) return ntypes + 1;

/* compatibility with old format where types are hard-coded in the routine */
ntypes = fixtypes;
for (num = 0; num < ntypes; num++) {
	(*type_codes)[num] = fixtypec[num];
	(*type_abbrevs)[num] = fixtypeab[num];
	sprintf(aux, "04%s", fixtypec[num]);
	(*type_ranks)[num] = fcode(ksmj, aux, 20);
	}
return ntypes;
#undef FIXTYPES
#undef MAXLOADTYPES
}


#ifdef unix 

void write_quick_meres(void)
/* creates index file MERES for fast db start
******** MUST BE USED WHEN DB IS CLOSED ******** 
*/
{
char *fname;
FILE *out;
int *blist;

dir_acnucopen("RO");
fname = prepare_env_var("acnuc");
strcat(fname, "MERES");
out = fopen(fname, "w");
blist = (int *)malloc(lenw * sizeof(int));
if(out != NULL && blist != NULL) {
	lngbit(2, blist);
	fwrite(blist, sizeof(int), lenw, out);
	fclose(out);
	free(blist);
	}
dir_acnucclose();
}


/* gestion de l'interruption d'un programme */
#include <signal.h>
/* private functions */
static void do_on_sigterm(int signum);

/* globals */
static int stop_next_point = FALSE;

int check_term(void)
{
static int first = TRUE;

if(first) {
	first = FALSE;
	signal(SIGTERM, do_on_sigterm);
	return FALSE;
	}
return stop_next_point;
}


static void do_on_sigterm(int signum)
{
if(signum == SIGTERM) stop_next_point = TRUE;
}

#endif
