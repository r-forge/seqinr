#include "dir_acnuc.h"
#include <math.h>
#include <ctype.h>
#define Boolean unsigned char

extern int *defbitlist;


/* prototypes des fonctions inclues */
void *mycalloc(int nbre, int size);
int notrail2(char *chaine, int len);
void prokey(int *slist, int *klist);
void prospe(int *sqlist, int *splist);
void proseq(int *slist, int *sqlist, char genre, int *aux);
void descen(DIR_FILE *kan, int recnum, int *blist);
void subloc(int *source, int *destin);
void locsub(int *source, int *destin);
int accseq(char *nom, int *list);
int bibseq(char *nom, int *list);
int smjseq(char *nom, int class, int *list);
int orgseq(char *nom,int *blist, int *aux);
int yeaseq(char *annee, char oper, int *blist, int *aux);
int autseq(char *nom, int *list);
void sel_seqs_1_node(DIR_FILE *kan, int num, int *seq_list, int hote);
int shkseq(char *nom, int *blist, int err);
int prepch(char *chaine, char **posmot);
int compch(char *cible, int lcible, char **posmot, int nbrmots);
int isenum_flou(char *nom, int *blist, int *locus);
void showannots(int iseq, Boolean *ctabchoix, char **record, int totrec, 
		int l_key, Boolean *nouveau_choix, void (*outoneline)(char *) );
int isfileokforwrite(char *fname, int *file_exists);
void memdelta(char *pos, int longueur, int delta);/* decalage en memoire */
int fileaccnums(FILE *fich,int *blist,int *locus);


void *mycalloc(int nbre, int size)
{
void *point;
point = calloc(nbre,size);
if(point == NULL) {
#ifdef __INTEL__
	FILE *tmp;
	tmp = fopen("queryerr.log", "w");
	fprintf(tmp, "Problem allocating memory.\n");
	fclose(tmp);
#else
	fprintf(stderr,"Problem allocating memory.\n");
#endif
	exit(ERREUR);
	}
return point;
}


int notrail2(char *chaine, int len)
{
len--;
while(len>=0 && chaine[len]==' ') len--;
return len+1;
}


void prokey(int *slist, int *klist)
{
int isub, next;
memset(klist, 0, longa * sizeof(int));
isub=1;
while( (isub=irbit(slist,isub,nseq))!=0) {
	readsub(isub);
	next=psub->plkey;
	while(next) {
		readshrt(next); next=pshrt->next;
		bit1(klist,pshrt->val);
		}
	}
}

void prospe(int *sqlist, int *splist)
{
int isub, next;
char nom[41];
nom[40]=0;
memset(splist, 0, longa * sizeof(int));
isub=1;
while( (isub=irbit(sqlist,isub,nseq))!=0) {
	readsub(isub);
	if(psub->pext>0) {
		readext(psub->pext);
		readsub(pext->mere);
		}
	readloc(psub->plinf);
	if(nbrf || swissprot) {
		next=ploc->spec;
		while(next) {
			readshrt(next); next=pshrt->next;
			bit1(splist,pshrt->val);
			}
		}
	else	bit1(splist,ploc->spec);
	}
}


void proseq(int *slist, int *sqlist, char genre, int *aux)
{
int num;
DIR_FILE *kan;
for (num = 0; num < lenw; num++) sqlist[num] = 0;
if(genre=='K') kan=kkey;
else kan=kspec;
num=1;
while( (num = irbit(slist, num, maxa))!=0 ) {
	dir_read(kan,num,1,pspec);
	if(pspec->plsub) {
		lngbit(pspec->plsub,aux);
		ou(sqlist,sqlist,aux,lenw);
		}
	}
}


static void descen_rec(int pdesc, int *blist)
{
int next;
readshrt(pdesc);
bit1(blist,abs(pshrt->val));
next=pshrt->next;
while(next) {
	readshrt(next); next=pshrt->next;
	descen_rec(pshrt->val,blist);
	}
}


void descen(DIR_FILE *kan, int recnum, int *blist)
{
int num;
memset(blist, 0, longa*sizeof(int));
dir_read(kan,recnum,1,pspec);
while(pspec->syno>0) {
	recnum=pspec->syno;
	dir_read(kan,recnum,1,pspec);
	}
descen_rec(pspec->desc,blist);
}


void subloc(int *source, int *destin)
{
int sn, point;
for(sn=0; sn<lenw; sn++) destin[sn]=0;
sn=1;
while( (sn=irbit(source,sn,nseq)) != 0 ) {
	readsub(sn); point=psub->pext;
	if(point<=0) bit1(destin,sn);
	else	{
		do	{
			readext(point); point=pext->next;
			bit1(destin,pext->mere);
			}
		while(point);
		}
	}
}


void locsub(int *source, int *destin)
{
int sn, point, val, i;
for(sn=0; sn<lenw; sn++) destin[sn]=source[sn];
sn=1;
while( (sn=irbit(source,sn,nseq)) != 0 ) {
	readsub(sn); point=psub->pext;
	if(point<0) {
		point= -point;
		do	{
			readlng(point); point=plng->next;
			for(i=0; i<SUBINLNG; i++)
				if( (val=plng->sub[i])!=0 ) bit1(destin,val);
			}
		while(point);
		}
	}
}


int accseq(char *nom, int *list)
{
int point, i;
point=fcode(kacc,nom,ACC_LENGTH);
if(!point)return 2;
for(i=0; i<lenw; i++) list[i]=0;
point=pacc->plsub;
if(!point)return 2;
do	{
	readshrt(point); point=pshrt->next;
	bit1(list,pshrt->val);
	}
while(point);
return 1;
}


int bibseq(char *nom, int *list)
{
int point, i;
point=fcode(kbib,nom,40);
if(!point)return 2;
for(i=0; i<lenw; i++) list[i]=0;
point=pbib->plsub;
if(!point)return 2;
do	{
	readshrt(point); point=pshrt->next;
	bit1(list,pshrt->val);
	}
while(point);
return 1;
}


int smjseq(char *nom, int class, int *list)
{
int point, i;
static char code[21];
sprintf(code,"%2.2d",class);
#ifdef vms
/* bug dans printf de vms */
if(*code==' ') *code='0';
#endif
strncpy(code+2,nom,18); code[20]=0;
point = fcode(ksmj,code,20);
if( point == 0 )return 2;
point = psmj->plong;
lngbit(point, list);
if( point == 0 )return 2;
return 1;
}


int orgseq(char *nom,int *blist, int *aux)
{
int last, num;
if( strncmp(nom,"NUCLEAR",7) ) {
	return smjseq(nom,5,blist);
	}
last=read_first_rec(ksmj,NULL);
memset(blist, 0, lenw*sizeof(int));
for (num=2; num<=last; num++) {
	readsmj(num);
	if(strncmp(psmj->name,"05",2) || !psmj->plong)continue;
	lngbit(psmj->plong,aux);
	ou(blist,blist,aux,lenw);
	}
non(blist,blist,lenw);
et(blist,blist,defbitlist + lenw,lenw);
return 1;
}


int yeaseq(char *annee, char oper, int *blist, int *aux)
{
int last, num, test, err=2;
last=read_first_rec(ksmj,NULL);
memset(blist, 0, lenw*sizeof(int));
for (num=2; num<=last; num++) {
	readsmj(num);
	if(strncmp(psmj->name,"03",2) || !psmj->plong)continue;
	test=strncmp(annee,psmj->name+2,4);
	if( (oper=='=' && !test) || (oper=='>' && test<=0) || 
			(oper=='<' && test>=0) ) {
		lngbit(psmj->plong,aux);
		ou(blist,blist,aux,lenw);
		err=1;
		if( oper=='=' && !test ) break;
		}
	}
return err;
}


static void numautseq(int plref,int *list)
{
int point2;
while(plref) {
	readshrt(plref); plref=pshrt->next;
	readbib(pshrt->val);
	point2=pbib->plsub;
	while(point2) {
		readshrt(point2); point2=pshrt->next;
		bit1(list,pshrt->val);
		}
	}
}


int autseq(char *nom, int *list)
{
int num;
for(num=0; num<lenw; num++) list[num]=0;

if( strchr(nom,'@') == NULL ) {
	num=fcode(kaut,nom,20);
	if(!num)return 2;
	num=paut->plref;
	if(!num)return 2;
	numautseq(num,list);
	}
else 	{
	int nbrmots, totaut;
	char *posmot[10];
	nbrmots=prepch(nom,posmot);
	totaut=read_first_rec(kaut,NULL);
	for(num=2; num<=totaut; num++) {
		readaut(num);
		if(paut->plref==0) continue;
		if(compch(paut->name,20,posmot,nbrmots))
			numautseq(paut->plref,list);
		}
	}
return 1;
}


void sel_seqs_1_node(DIR_FILE *kan, int num, int *seq_list, int hote)
{
int point;
do	{
	if(dir_read(kan,num,1,pspec)!=1) dir_readerr(kan,num);
	if(pspec->syno > 0) num=pspec->syno;
	}
while(pspec->syno > 0);
if(hote)
	point=pspec->plhost;
else
	point=pspec->plsub;
while(point !=0) {
	int i;
	readlng(point); point=plng->next;
	for(i=0; i<SUBINLNG;i++) if(plng->sub[i]) bit1(seq_list,plng->sub[i]);
	}
point=pspec->desc;
readshrt(point);
point=pshrt->next;
while(point!=0) {
	readshrt(point); point=pshrt->next;
	readshrt(pshrt->val);
	sel_seqs_1_node(kan,abs(pshrt->val),seq_list,hote);
	}
}


int shkseq(char *nom, int *blist, int err)
{
int num, ncoup, last, hote;
char *tab[20], nomc[61];
DIR_FILE *kan;

hote=0;
if(err==1) kan=kspec;
else if(err==2) { kan=kspec; hote=1; }
else if(err==3) kan=kkey;
else exit(1);

memset(blist, 0, lenw*sizeof(int));

strcpy(nomc,nom);
ncoup=prepch(nomc,tab);
if(ncoup>0) goto rech_floue;

num=iknum(nom,kan);
if(num==0) {
	return 2;
	}
sel_seqs_1_node(kan,num,blist,hote);
return 1;

rech_floue:
err=2;
last=read_first_rec(kan,NULL);
for(num=3; num<=last; num++) {
	dir_read(kan,num,1,pspec);
	if(pspec->desc == 0 && pspec->syno == 0)continue;
	if( compch(pspec->name,40,tab,ncoup) ) {
		err=1;
		sel_seqs_1_node(kan,num,blist,hote);
		}
	}
return err;
}


int prepch(char *chaine, char **posmot)
{
/*
chaine: template a rechercher qui contient des wildcard @
posmot: tableau de pointeurs vers char au retour rempli avec des pointeurs adequats qui pointent dans chaine qui ne doit plus etre modifiee
valeur rendue: nbre de pointeurs dans tableau posmot
*/
char *pos;
int nbrmots;
static char wildcard='@';

if(strchr(chaine,'@')==NULL) return 0;
nbrmots= -1;
pos=chaine+strlen(chaine)-1;
while( pos>=chaine && *pos==' ' ) pos--;
*(pos+1)=0;

pos=chaine;
while(*pos!=0) {
	if(*pos==wildcard) {
		posmot[++nbrmots]=NULL;
		*pos=0;
		while(*(pos+1)==wildcard) pos++;
		}
	else	{
		posmot[++nbrmots]=pos;
		while( *(pos+1)!=wildcard && *(pos+1) !=0 ) pos++;
		}
	pos++;
	}
return nbrmots+1;
}


int compch(char *cible, int lcible, char **posmot, int nbrmots)
{
/*
cible: chaine a tester pour presence du template
lcible: long. de cible qui n'est pas forcement finie par \0
	doit etre <=60
posmot: tableau fabrique par prepch
nbrmots: valeur rendue par prepch
valeur rendue: 1 ssi template present dans cible, 0 si absent
*/
int num= 0, l, total;
char *pos;
static char vcible[61];

pos=cible+lcible-1;
while( pos>=cible && *pos==' ' ) pos--;
lcible=pos-cible+1;
memcpy(vcible,cible,lcible);
vcible[lcible]=0;
cible=vcible;
if(posmot[nbrmots-1]==NULL)
	total=nbrmots-1;
else
	total=nbrmots-2;

if(posmot[0]!=NULL) { /* comparaison avec mot initial */
	l=strlen(posmot[0]);
	if(strncmp(cible,posmot[0],l)!=0) return 0;
	cible += l;
	num++;
	}
while(num<total) { /* recherche des mots internes */
	num++;
	pos=strstr(cible,posmot[num]);
	if(pos==NULL) return 0;
	l=strlen(posmot[num]);
	cible = pos+l;
	num++;
	}
if( total==nbrmots-1 ) return 1; /* template se termine par @ */
/* test si cible se termine par dernier mot du template */
l=strlen(posmot[nbrmots-1]); 
if( strcmp(vcible+lcible-l,posmot[nbrmots-1]) == 0 ) return 1;
return 0;
}


int isenum_flou(char *nom, int *blist, int *locus)
/* rend 0 si pas trouve
>0 si 1 seq trouvee
<0 si flou et trouve et alors les bits sont mis a 1
*locus rendu a FALSE si au moins 1 fille dedans
*/
{
int nbrpoints, debut, rapide, fin, num, trouve, l_partie_mere;
char *tab[20], *p;
if(strchr(nom,'@')==NULL) {
	num = isenum(nom);
	if(*locus && psub->pext>0) *locus=0;
	return num;
	}
fin=nseq;
nbrpoints=prepch(nom,tab);
read_first_rec(ksub,&num);
if(tab[0]!=NULL && num==nseq) {
	rapide=1;
	l_partie_mere=strlen(tab[0]);
	p=strchr(tab[0],'.');
	if( p!=NULL ) l_partie_mere=(p-tab[0]);
	debut=fcode(ksub,tab[0],l_partie_mere);
	if(debut==0) return 0;
	}
else	{ 
	rapide=0;
	debut=2;
	}
num=debut-1;
trouve=0;
do	{
	num++;
	readsub(num);
	if(compch(psub->name,16,tab,nbrpoints)) {
		bit1(blist,num);
		if(*locus && psub->pext>0) *locus=0;
		trouve=1;
		}
	}
while(num<fin && (!rapide || strncmp(psub->name,tab[0],l_partie_mere)==0 ));
if(rapide) {
	num=debut;
	do	{
		num--;
		readsub(num);
		if(compch(psub->name,16,tab,nbrpoints)) {
			bit1(blist,num);
			trouve=1;
			if(*locus && psub->pext>0) *locus=0;
			}
		}
	while(num>2 && strncmp(psub->name,tab[0],l_partie_mere)==0 );
	}
if(trouve) return -1;
else return 0;
}


static void showfeatures(int mere, long *pannot, int div, void (*outoneline)(char *) )
/* traitement table des features d'une mere GenBank ou EMBL
   mere: numero de la seq mere
   *pannot: pointeur vers debut FEATURES de mere 
   *outoneline:  fonction d'ecriture des lignes produites
*/
{
char *p, extens[50], nom_mere[16];
long *tabsub=NULL;
int pext, nft, nf, le, lk, block=0, point, fille, div_fille;
Boolean show;

readsub(mere); pext=psub->pext;
memcpy(nom_mere,psub->name,16);

(*outoneline)(pinfo->line);
if(embl && !strncmp(pinfo->line,"FH",2) ) {
		next_annots(pannot);
		(*outoneline)(pinfo->line);
		}
encore:
nft= -1;         /* memoriser la liste des filles */
block += 100;
if(tabsub!=NULL) free(tabsub);
tabsub=(long *)mycalloc(block, sizeof(long));
point=abs(pext);
while (point) {
	readlng(point); point=plng->next;
	for(nf=0; nf<SUBINLNG; nf++) {
		if(plng->sub[nf] == 0) break;
		if( (++nft) >= block ) goto encore; /*recommencer avec + de memoire*/
		tabsub[nft]=plng->sub[nf];
		}
	}
for(nf=0; nf<=nft; nf++) {  /* pour chaque fille */
	fille = tabsub[nf];
	readsub(fille);
	p=strchr(psub->name,'.');
/* si autre mere principale */
	if( strcmptrail(psub->name,p-psub->name,nom_mere,16) ) {
		sprintf(extens,"** used in      %16.16s",psub->name);
		(*outoneline)(extens);
		tabsub[nf]=0;
		continue;
		}
	/* lecture annots de la fille */
	seq_to_annots(fille, &tabsub[nf], &div_fille);
/*	if(tabsub[nf]==0) continue; */
	if( read_annots(tabsub[nf], div_fille) == NULL) continue;
	readsub(fille);
	memset(extens,' ',20); extens[20]=0;
	if(embl) {
		memcpy(extens,"FT ",3);
		memcpy(extens+3,p,psub->name-p+16);
		}
	else
		memcpy(extens,p,psub->name-p+16);
	p=extens+19;   /* calcul longueur de l'extension */
	while(*p==' ') p--;
	le=p-extens+1;
	if(le<4) le=4;
	p=pinfo->line+19;  /* calcul longueur de la cle de la feature */
	while(*p==' ') p--;
	lk=p-(pinfo->line+5)+1;
	if(le+lk>19) le=19-lk;
	extens[le]=' ';
	p = pinfo->line;
	memcpy(extens+le+1,p+5,lk);
	do	{        /* ecriture de la feature et de ses lignes associees */
		if(*p != 0) memcpy(p,extens,20);
		(*outoneline)(p);
		memset(extens,' ',20);
		if(embl) memcpy(extens,"FT",2);
		p = next_annots(NULL);
		}
	while( *p == 0 || !strcmptrail(p,10,NULL,0) || 
			!strcmptrail(p,10,"FT",2) );
	}
/* write other features */
read_annots(*pannot, div);
p = next_annots(pannot);
while( !strcmptrail(p,4,NULL,0) || !strncmp(p,"FT",2) ){
		show=TRUE;
		for(nf=0; nf<=nft; nf++) { /* est-ce feature deja traitee? */
			if(tabsub[nf] == *pannot) {
				show=FALSE;
				break;
				}
			}
		do	{ /* ecrire feature et ses lignes associees */
			if(show) (*outoneline)(p);
			p = next_annots(pannot);
			}
		while( *p == 0 || !strcmptrail(p,10,NULL,0) ||
			!strcmptrail(p,10,"FT",2) );
		}
free(tabsub);
return;
}


void showannots(int iseq, Boolean *ctabchoix, char **record, int totrec, 
		int l_key, Boolean *nouveau_choix, void (*outoneline)(char *) )
/* traitement des annotations
iseq: num de la seq
ctabchoix: tableau des options retenues (from 0)
record: tableau des noms des records existants (from 0)
totrec: nbre total d'options (from 1)
l_key: longueur des noms de records
*nouveau_choix: TRUE ssi ctabchoix a change
outoneline: fonction d'ecriture de la ligne traitee
*/
{
char line[81], *p;
int i, j, k, div, longueur;
static int *tabchoix=NULL, lastrec;
long pannot;

if(*nouveau_choix) { /* preparation des options */
	if(tabchoix==NULL) tabchoix=(int *)mycalloc(totrec,sizeof(int));
	for(i=0; i<totrec; i++) tabchoix[i]=ctabchoix[i];
	*nouveau_choix=FALSE;
/* add subfields */
	if(embl) {
		if(tabchoix[7]) tabchoix[6]=TRUE; /* OC==> OS */
		if(tabchoix[6]) tabchoix[8]=TRUE; /* OS==> OG */
		if(tabchoix[19]) tabchoix[18]=TRUE; /* FT==> FH */
		if(tabchoix[9]) { /* RN==> RC+RP+RX+RA+RT+RL */
			tabchoix[10]=tabchoix[11]=tabchoix[12]=
				tabchoix[13]=tabchoix[14]=tabchoix[15]=TRUE;
			}
		}
	else if(swissprot) {
		if(tabchoix[6]) tabchoix[4]=TRUE; /* OC==> OS */
		if(tabchoix[4]) tabchoix[5]=TRUE; /* OS==> OG */
		if(tabchoix[8]) { /* RN==> RP+RC+RX+RA+RL */
			tabchoix[9]=tabchoix[10]=tabchoix[11]=
				tabchoix[12]=tabchoix[13]=TRUE;
			}
		}
	else if(!nbrf) { /* GenBank */
		if(tabchoix[6]) tabchoix[7]=TRUE; /* SOURCE==> ORGANISM */
		if(tabchoix[8]) { /* REFERENCE==> AUT+TIT+JOU+MED+REM */
			tabchoix[9]=tabchoix[10]=tabchoix[11]=tabchoix[12]=
				tabchoix[13]=TRUE;
			}
		}
/* find lastrec the # of the last record type after those asked for */
	for(lastrec=totrec-2; lastrec>=1; lastrec--) {
		if(tabchoix[lastrec]) break;
		}
	lastrec++;
	if(embl) {
		if(lastrec < 5) lastrec = 5; /* DE apres AC+SV+NI+DT */
		if(lastrec>=10 && lastrec<=15) lastrec=16; /* for R- records */
		}
	else if(swissprot) {
		if(lastrec<3) lastrec=3; /* DE apres AC+DT */
		if(lastrec>=8 && lastrec<=12) lastrec=13; /* for R- records */
		}
	else if(!nbrf) { /* GenBank */
		if(lastrec>=9 && lastrec<=13) lastrec=14; /* for REFER records*/
		}
	}
/* traitement de la seq */
seq_to_annots(iseq, &pannot, &div);
readsub(iseq);
longueur=psub->length;
if( psub->pext>0 ) { /* subsequence */
	Boolean ignoretrans=FALSE;
	sprintf(line,"%16.16s     Location/Qualifiers    (length=%d bp)",
			psub->name,psub->length);
	(*outoneline)(line);
	if( read_annots(pannot, div) == NULL) return;
	while(TRUE) {
		if(!ignoretrans) (*outoneline)(pinfo->line);
		next_annots(&pannot);
		if( strcmptrail(pinfo->line+2,8,NULL,0) ) break;
		if( strstr(pinfo->line," /translation=") != NULL ) 
			ignoretrans=TRUE;
		else if(ignoretrans && strchr(pinfo->line,'/') != NULL) 
			ignoretrans=FALSE;
		}
	}
else { /* sequence mere */
	Boolean showit, embl_like;
	int choix, lcle;
	char cle[20];

	lcle=l_key;
	if(nbrf) lcle=2; /* ignore subfields for NBRF */
	if( read_annots(pannot, div) == NULL) return;
	(*outoneline)(pinfo->line);
	embl_like = (embl || swissprot);
/* recherche des cibles et de DE / DEFINITION / TITLE */
	next_annots(&pannot);
	do	{
		showit=FALSE;
		if( tabchoix[0] || !strncmp(pinfo->line,"DE",2) ||
			(nbrf && !strncmp(pinfo->line,"TITLE",5) ) ) 
			showit=TRUE;
		if( !showit ) {
			for(choix=1; choix<=totrec-2; choix++) {
			     if( !strncmp(pinfo->line,record[choix],
					notrail2(record[choix],l_key) ) ) {
				   if(tabchoix[choix]) { 
					showit=TRUE; break; }
				   if(choix>=lastrec) goto sequence;
				   }
			     }
			}
		memcpy(cle,pinfo->line,l_key);
		if( showit && !nbrf && ( !strncmp(pinfo->line,"FEAT",4) || 
			!strncmp(pinfo->line,"FH",2) )  )
			showfeatures(iseq,&pannot,div,outoneline);
		else 	{ /* normal treatment of records */
			do 	{
				if( showit ) (*outoneline)(pinfo->line);
				if( next_annots(&pannot) == NULL) break;
				}
			while( (embl_like && !strncmp(pinfo->line,cle,l_key) ) 
				|| 
			(!embl_like && !strcmptrail(pinfo->line,lcle,NULL,0)) );
			}
		}
	while( strncmp(cle,record[totrec-2],notrail2(record[totrec-2],l_key)));
	}

sequence:
if(tabchoix[totrec-1]) {      /* ecriture de la sequence */
    char outlig[81];
    if(!nbrf) {
	for(i=1; i<=longueur; i+=60) {
		k = gfrag(iseq,i,60,line);
		if(k<60) memset(line+k,' ',60-k);
		if(embl) { /* embl: numer des lignes */
			memset(outlig,' ',4);
			sprintf(outlig+70,"%10d",i+k-1);
			p=outlig+4;
			outlig[80]=0;
			}
		else if( swissprot ) { /* swissprot: pas de numer des lignes */
			memset(outlig,' ',4);
			p=outlig+4;
			outlig[70]=0;
			}
		else  { /* GenBank: numer des lignes */
			sprintf(outlig,"%9d",i);
			p=outlig+9;
			outlig[75]=0;
			}
		for(j=0; j<60; j+=10) { /* les residus */
			*p=' ';
			memcpy(++p,line+j,10);
			p+=10;
			}
		(*outoneline)(outlig);
		}
	}
    else {
	int vlong, l, nlu;
	(*outoneline)("SEQUENCE");
	strcpy(outlig,"               ");
	j=5;
	for(i=1;i<=6;i++) {
		sprintf(line,"%2d        ",j);
		strcat(outlig,line);
		j+=5;
		}
	(*outoneline)(outlig);
	vlong=2*longueur;
	if( longueur % 30 != 0 ) ++vlong;
	gfrag(iseq,1,1,line); /* charge tout avec ponctuation dans nucbuf */
	l=0; i=1;
	do	{
		nlu=l+60; if(nlu>vlong) nlu=vlong;
		sprintf(outlig,"%7d",i);
		p=outlig+7;
		for(j=l;j<nlu;j++) *(p++) = nucbuf[j];
		*p=0;
		(*outoneline)(outlig);
		l=nlu; i+=30;
		}
	while(l<vlong);
	}
    }
if(nbrf)
	(*outoneline)("///");
else
	(*outoneline)("//");
return;
} /* end of showannots */


int isfileokforwrite(char *fname, int *file_exists)
{
FILE *fich;
int retval=1;

fich=fopen(fname,"r");
if(fich==NULL)
	*file_exists=0;
else	{
	*file_exists=1;
	fclose(fich);
	}
fich=fopen(fname,"a");
if(fich==NULL)
	retval=0;
else
	fclose(fich);
return retval;
}


void memdelta(char *pos, int longueur, int delta)
{
/* pour decaler la chaine pos de longueur longueur de quantite delta
a droite si delta>0 a gauche si delta<0  et rien si delta=0 */
char *fin;
if(delta == 0) return;
if(delta > 0) {
	fin = pos;
	pos += longueur;
	while(--pos >= fin) 
		*(pos+delta) = *pos;
	}
else	{
	fin = pos + longueur - 1;
	pos--;
	while(++pos <= fin) 
		*(pos+delta) = *pos;
	}
}	


int fileaccnums(FILE *fich,int *blist,int *locus)
/* returns !=0 if no seq found */
{
char nom[42], *p, mnemo[17], *fin;
int nacc, next, nsub, val, retval=1;

*locus=TRUE;
while(fgets(nom,40,fich)!=NULL) {
	nom[strlen(nom)-1]=0;
	p=nom+strlen(nom)-1;
	while(*p==' ') p--;
	*(p+1)=0;
	p=nom-1;
	while(*(++p)!=0) *p=toupper(*p);
	p=strchr(nom,'.');
	if(p!=NULL) *p=0;
	nacc=fcode(kacc,nom,ACC_LENGTH);
	if(nacc==0) continue;
	readacc(nacc);
	next=pacc->plsub;
	while(next) {
		readshrt(next); next=pshrt->next;
		val=pshrt->val;
		if(p!=NULL) {
			readsub(val);
			memcpy(mnemo,psub->name,16);
			fin=strchr(mnemo,' ');
			*fin='.';
			strcpy(fin+1,p+1);
			nsub=isenum(mnemo);
			if(nsub!=0) {
				*locus=FALSE;
				val=nsub;
				}
			}
		bit1(blist,val);
		retval=0;
		}
	}
return retval;
}
