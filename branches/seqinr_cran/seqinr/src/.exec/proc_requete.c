#include "dir_acnuc.h"
#include <ctype.h>
#define Boolean unsigned char
#define boolean_differ(a,b) ((a && !b) || (!a && b))


/* extern variables */
extern int deflocus[], tlist;
extern int defoccup[];
extern char defgenre[];
extern int *defbitlist;
extern int defllen[];
extern char *deflnames[];

/* prototypes of included functions */
int proc_requete(char *requete, char *message, char *nomliste, int *numliste);
int getfreelist(int toplist);
static int syntax(char *requete, char *finreq, 
	int *numliste);
static char *nextquote(char *debut, char *fin);
static char *nextpar(char *debut, char *fin);
static int unary(char *oper, int numliste);
static int binary(char *oper, int r1, int r2);
static int operand(char *requete, char *finreq, int *numliste);
static void addquotes(char *requete, char *lcrequete);
static char *prepaddquotes(char *requete, char *lcrequete, char *pospar,
	 int *lreq);

/* used functions */
void *mycalloc(int nbre, int taille);
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
int isenum_flou(char *nom, int *blist, int *locus);
int fileaccnums(FILE *fich,int *blist,int *locus);
void memdelta(char *pos, int longueur, int delta);/* decalage en memoire */
int shkseq(char *nom, int *blist, int err);


/* global variables */
static char *debutsrequetes[2];
static char *poserr;
#define NUMOPER 10
static char operateur[]="ET/OU/NO/ME/FI/PS/PK/UN/SD/KD/";


int proc_requete(char *requetenoquotes, 
	char *message, char *nomliste, int *numliste)
{
static char *errmess[]={
	"ok", "syntax error", "not enough memory available, delete some list",
	"unknown species","unknown keyword","unknown operand",
	"unknown status, molecule, journal, year, type or organelle",
	"unknown reference","unknown author","unknown sequence name",
	"unknown accession number","file does not exist or is not adequate",
	"unknown list","operation invalid with its operands"
	};
int numerr, i, *oldoccup;
Boolean test;
char upnomliste[11], *p, *requete, *lcrequete;

strcpy(upnomliste,nomliste);
p=upnomliste;
while(*p!=0) {
	*p=toupper(*p);
	p++;
	}
i=strlen(requetenoquotes);
requete=mycalloc(i+30,1);
lcrequete=mycalloc(i+30,1);
strcpy(requete,requetenoquotes);
strcpy(lcrequete,requetenoquotes);
p=requete-1; while(*(++p)!=0) *p=toupper(*p);
if(strchr(requete,'"')==NULL) {
	addquotes(requete,lcrequete);
	}
debutsrequetes[0]=requete; debutsrequetes[1]=lcrequete;
if( getfreelist(1) == -1) { /* plus de liste vide */
	numerr=2;
	}
else {
/* memoriser l'etat d'occupation des listes */
    oldoccup=(int *)mycalloc(tlist,sizeof(int));
    memcpy(oldoccup,defoccup,tlist*sizeof(int));
    numerr=syntax(requete,requete-1+strlen(requete),numliste);
    if(numerr==0) { /* sortie sans erreur */
	i=1;
	do	{
		i++; 
		test=defoccup[i] && (i != *numliste) &&
		    !strcmp(upnomliste,deflnames[i]);
		}
	while( !test && i<tlist-1);
	if(test) {/* new list replaces a previous one */
		memcpy(defbitlist+i*lenw,defbitlist+(*numliste)*lenw,
				lenw*sizeof(int));
		defgenre[i]=defgenre[*numliste];
		deflocus[i]=deflocus[*numliste];
		defoccup[*numliste]=FALSE;
		*numliste=i;
		}
	else 	{ /* new list is new */
		deflnames[*numliste]=(char *)mycalloc(strlen(upnomliste)+1,1);
		strcpy(deflnames[*numliste],upnomliste);
		}
/* croiser la nouvelle liste avec liste des valides pour eliminer les seqs 
supprimees a l'installation par connectindex 
*/
	if(defgenre[*numliste] == 'S')
		et(defbitlist+(*numliste)*lenw, defbitlist+(*numliste)*lenw, 
			defbitlist, lenw);
	defllen[*numliste]=bcount(defbitlist+(*numliste)*lenw,lenbit);
	}
    else { /* sortie avec erreur:remettre ancien etat d'occupation des listes */
	memcpy(defoccup,oldoccup,tlist*sizeof(int));
	}
    free(oldoccup);
    }
if(numerr) {
	if(numerr==2) /* pas de position d'erreur pour memoire insuffisante */
		poserr=NULL;
	sprintf(message,"%s",errmess[numerr]);
	if(poserr!=NULL) {
		char keep;
		poserr= (poserr-requete)+lcrequete;
		keep= *poserr; *poserr=0;
		sprintf(message+strlen(message)," at (^): %s(^)",lcrequete);
		*poserr=keep;
		sprintf(message+strlen(message),"%s",poserr);
		}	
	}
free(requete); free(lcrequete);
return numerr;
} /* end of proc_requete */


int getfreelist(int toplist)
{
Boolean found=FALSE;
while(++toplist<tlist) {
	if(!defoccup[toplist]) {found=TRUE; break;}
	}
if(!found) return -1;
memset(defbitlist+toplist*lenw,0,lenw*sizeof(int));
return toplist;
}


static int syntax(char *requete, char *finreq, 
	int *numliste)
{
static int prio[NUMOPER];
static Boolean unaire[NUMOPER];
static Boolean first=TRUE;
int numerr, num, minp, goodnum;
Boolean test;
char *pos, oper[3], *p, *goodpos;
if(first) {
	int i;
	first=FALSE;
	unaire[0]=unaire[1]=FALSE;
	for(i=2; i<NUMOPER; i++) unaire[i]=TRUE;
	prio[0]=2;
	for(i=1; i<NUMOPER; i++) prio[i]=1;
	}
do	{
	while(*requete==' ' && requete<finreq) requete++;
	if(requete>finreq || *requete==' ') {
		numerr=1; poserr=requete; return numerr;
		}
	while(*finreq==' ') finreq--;
	test=FALSE;
	if(*requete=='(') {
		/* elimination des () exterieures */
		pos=nextpar(requete,finreq);
		if(pos==NULL) {
			numerr=1; poserr=requete; return numerr;
			}
		if(pos==finreq) {
			requete++; finreq--; test=TRUE;
			}
		}
	}
while(test);
if(*requete=='"') {
	/* cas "operande" */
	pos=nextquote(requete,finreq);
	if(pos==NULL) {
		numerr=1; poserr=requete; return numerr;
		}
	if(pos==finreq) {
		numerr=operand(requete,finreq,numliste);
		return numerr;
		}
	}
pos=requete; minp=10;
do	{
	/* recherche de l'operateur de prorite minimum le plus a droite */
	while(pos<finreq && *pos==' ') pos++;
	if(*pos=='(') {
		p=nextpar(pos,finreq);
		if(p==NULL) {
			numerr=1; poserr=pos; return numerr;
			}
		pos= ++p;
		}
	else if(*pos=='"') {
		p=nextquote(pos,finreq);
		if(p==NULL) {
			numerr=1; poserr=pos; return numerr;
			}
		pos = ++p;
		}
	while(pos<finreq && *pos==' ') pos++;
	if(pos<finreq) {
		memcpy(oper,pos,2); oper[2]=0;
		p=strstr(operateur,oper);
		if(p==NULL) {
			numerr=1; poserr=pos; return numerr;
			}
		num=(p-operateur)/3;
		if(prio[num]<=minp) {
			minp=prio[num]; goodpos=pos; goodnum=num;
			}
		}
	pos += 2;
	}
while(pos<=finreq);
while(goodpos>requete && unaire[goodnum]) {
	/* si l'operateur retenu est unaire et pas a l'extreme gauche 
	considerer l'operateur precedant */
	pos=goodpos-1;
	while(*pos==' ') pos--;
	pos--;
	memcpy(oper,pos,2); oper[2]=0;
	p=strstr(operateur,oper);
	if(p==NULL) {
		numerr=1; poserr=goodpos; return numerr;
		}
	goodnum=(p-operateur)/3;
	goodpos=pos;
	}
if(!unaire[goodnum]) {
	/* operateur binaire */
	char *d1, *d2, *f1, *f2;
	int totg=0, totd=0, r1, r2;
	for(pos=requete; pos<goodpos; pos++) if(*pos=='"') totg++;
	for(pos=goodpos+2; pos<=finreq; pos++) if(*pos=='"') totd++;
	if(totg<totd) {
		d1=goodpos+2; f1=finreq; d2=requete; f2=goodpos-1;
		}
	else {
		d1=requete; f1=goodpos-1; d2=goodpos+2; f2=finreq;
		}
	numerr=syntax(d1,f1,&r1);
	if(numerr) return numerr;
	numerr=syntax(d2,f2,&r2);
	if(numerr) return numerr;
	numerr=binary(goodpos,r1,r2);
	if(numerr) {
		poserr=goodpos; return  numerr;
		}
	defoccup[r2]=FALSE;
	if(deflnames[r2]!=NULL) {
		free(deflnames[r2]); deflnames[r2]=NULL;
		}
	*numliste=r1;
	}
else {
	/* operateur unaire */
	numerr=syntax(requete+2,finreq,numliste);
	if(numerr) return numerr;
	numerr=unary(goodpos,*numliste);
	if(numerr) {
		poserr=goodpos; return numerr;
		}
	}
return 0;
} /* end of syntax */


static char *nextquote(char *debut, char *fin)
{
debut++;
while(debut<fin && *debut!='"') debut++;
if(debut>fin || *debut!='"') return NULL;
return debut;
}

static char *nextpar(char *debut, char *fin)
{
debut++;
do	{
	if(*debut=='(') {
		debut=nextpar(debut,fin);
		if(debut==NULL) return NULL;
		debut++;
		}
	if(*debut=='"') {
		debut=nextquote(debut,fin);
		if(debut==NULL) return NULL;
		debut++;
		}
	if(*debut==')') {
		return debut;
		}
	else debut++;
	}
while(debut<=fin);
return NULL;
}


static int unary(char *oper, int numliste)
{
int r2, r3, *etliste;
r2=getfreelist(-1);
if(r2== -1) return 2;
if(!strncmp(oper,"ME",2)) {
	if(defgenre[numliste]!='S') return 13;
	subloc(defbitlist+numliste*lenw,defbitlist+r2*lenw);
	deflocus[r2]=TRUE;
	defgenre[r2]='S';
	}
else if(!strncmp(oper,"FI",2)) {
	if(defgenre[numliste]!='S') return 13;
	locsub(defbitlist+numliste*lenw,defbitlist+r2*lenw);
	deflocus[r2]=FALSE;
	defgenre[r2]='S';
	}
else if(!strncmp(oper,"NO",2)) {
	non(defbitlist+r2*lenw, defbitlist+numliste*lenw,lenw);
	if(defgenre[numliste]=='S') {
		if(deflocus[numliste]) etliste=defbitlist+lenw;
		else etliste=defbitlist;
		et(defbitlist+r2*lenw,defbitlist+r2*lenw,etliste,lenw);
		deflocus[r2]=deflocus[numliste];
		}
	else	{
		int i;
		bit0(defbitlist+r2*lenw,1);
		bit0(defbitlist+r2*lenw,2);
		if(lenw-longa>0)memset(defbitlist+r2*lenw+longa,0,
						(lenw-longa)*sizeof(int));
		for(i=maxa+1; i<=longa*lmot; i++) bit0(defbitlist+r2*lenw,i);
		}
	defgenre[r2]=defgenre[numliste];
	}
else if(!strncmp(oper,"PK",2)) {
	if(defgenre[numliste] != 'S') return 13;
	prokey(defbitlist + numliste * lenw, defbitlist + r2 * lenw);
	/* car prokey ne le fait que pour longa */
	if(lenw > longa) memset(defbitlist + r2 * lenw + longa, 0,
				(lenw - longa) * sizeof(int));
	defgenre[r2] = 'K';
	}
else if(!strncmp(oper,"PS",2)) {
	if(defgenre[numliste] != 'S') return 13;
	prospe(defbitlist + numliste * lenw, defbitlist + r2 * lenw);
	/* car prospe ne le fait que pour longa */
	if(lenw > longa) memset(defbitlist + r2 * lenw + longa, 0,
				(lenw - longa) * sizeof(int));
	defgenre[r2] = 'E';
	}
else if(!strncmp(oper,"UN",2)) {
	int i;
	if(defgenre[numliste]=='S') return 13;
	r3=getfreelist(r2);
	if(r3== -1) return 2;
	proseq(defbitlist+numliste*lenw,defbitlist+r2*lenw,defgenre[numliste],
		defbitlist+r3*lenw);
	defgenre[r2]='S';
/* y a-t-il des filles dans la liste r2 ? */
	et(defbitlist+r3*lenw, defbitlist+r2*lenw, defbitlist+lenw, lenw);
	deflocus[r2]=TRUE;
	for(i=0; i<lenw; i++) {
		if( *(defbitlist+r3*lenw+i) != *(defbitlist+r2*lenw+i) ) {
			deflocus[r2]=FALSE;
			break;
			}
		}
	}
else if( !strncmp(oper,"KD",2) || !strncmp(oper,"SD",2) ) {
	DIR_FILE *kan;
	int nkey;
	r3=getfreelist(r2);
	if(r3== -1) return 2;
	if(oper[0]=='K') {
		if(defgenre[numliste]!='K') return 13;
		kan=kkey; defgenre[r2]='K';
		}
	else {
		if(defgenre[numliste]!='E') return 13;
		kan=kspec; defgenre[r2]='E';
		}
	nkey=2;
	while ( (nkey=irbit(defbitlist+numliste*lenw,nkey,maxa)) != 0 ) {
		descen(kan,nkey,defbitlist+r3*lenw);
		ou(defbitlist+r2*lenw,defbitlist+r2*lenw,
					defbitlist+r3*lenw,lenw);
		}
	}
memcpy(defbitlist+numliste*lenw, defbitlist+r2*lenw,lenw*sizeof(int));
deflocus[numliste]=deflocus[r2];
defgenre[numliste]=defgenre[r2];
return 0;
} /* end of unary */


static int binary(char *oper, int r1, int r2)
{
int r3, *dr1, *dr2, *dr3;
dr1=defbitlist+r1*lenw; dr2=defbitlist+r2*lenw;
if(defgenre[r1]!=defgenre[r2]) return 13;
if(!strncmp(oper,"OU",2)) {
	if(defgenre[r1]=='S' && boolean_differ(deflocus[r1] , deflocus[r2]) ) {
		r3=getfreelist(-1);
		if(r3== -1) return 2;
		dr3=defbitlist+r3*lenw;
		if(deflocus[r1]) {
			locsub(dr1,dr3);
			ou(dr1,dr3,dr2,lenw);
			}
		else {
			locsub(dr2,dr3);
			ou(dr1,dr1,dr3,lenw);
			}
		}
	else ou(dr1,dr1,dr2,lenw);
	deflocus[r1]=deflocus[r1] && deflocus[r2];
	}
else if(!strncmp(oper,"ET",2)) {
	if(defgenre[r1]=='S' && boolean_differ(deflocus[r1] , deflocus[r2] ) ) {
		r3=getfreelist(-1);
		if(r3== -1) return 2;
		dr3=defbitlist+r3*lenw;
		if(deflocus[r1]) {
			locsub(dr1,dr3); et(dr1,dr2,dr3,lenw);
			}
		else {
			locsub(dr2,dr3); et(dr1,dr1,dr3,lenw);
			}
		}
	else et(dr1,dr1,dr2,lenw);
	deflocus[r1]=deflocus[r1] && deflocus[r2];
	}
return 0;
} /* end of binary */


static int operand(char *requete, char *finreq, int *numliste)
{
int *dr, m3, *dm3, numerr, lnom, i;
char *dnom, nom[82];

dnom=requete;
requete++; finreq--;
while(requete <finreq && *requete==' ') requete++;
while(finreq>requete && *requete==' ') finreq--;
if(finreq<requete || *requete==' ') {
	poserr=finreq;
	return 1;
	}
dnom=strchr(requete,'=');
if(dnom==NULL || dnom>=finreq)
	dnom=strchr(requete,'<');
if(dnom==NULL || dnom>=finreq)
	dnom=strchr(requete,'>');
if(dnom==NULL || dnom>finreq) dnom=requete;
else dnom++;
while(dnom<finreq && *dnom==' ') dnom++;
if(dnom>finreq) {
	poserr=requete; return 5;
	}
lnom=finreq-dnom+1;
memcpy(nom,dnom,lnom); 
i=lnom; while(i<40) nom[i++]=' '; nom[i]=0;
*numliste=getfreelist(-1);
if(*numliste== -1) {
	poserr=dnom; return 2;
	}
dr=defbitlist+(*numliste)*lenw;
if(!strncmp(requete,"SP=",3) || !strncmp(requete,"H=",2) ||
		!strncmp(requete,"K=",2) ) {
	/* species host or keyword */
	int option;
	if(*requete=='S') option=1;
	else if(*requete=='H') option=2;
	else option=3;
	numerr=shkseq(nom,dr,option);
	if(numerr!=1) {
		if(option==3) numerr=4;
		else numerr=3;
		poserr=dnom;
		return numerr;
		}
	if(option==3) {
		m3=getfreelist(*numliste);
		if(m3== -1) {  poserr=dnom; return 2; }
		dm3=defbitlist+m3*lenw;
		/* in list m3: subsequences */
		non(dm3,defbitlist+lenw,lenw);
		et(dm3,defbitlist,dm3,lenw);
		/* in m3: subsequences in list */
		et(dm3,dm3,dr,lenw);
		/* locus vrai ssi aucune fille dans liste */
		deflocus[*numliste]=( irbit(dm3,1,nseq)==0 );
		}
	else
		deflocus[*numliste]=TRUE;
	defgenre[*numliste]='S';
	}
else if(!strncmp(requete,"T=",2) ||
	!strncmp(requete,"ST=",3) ||
	!strncmp(requete,"M=",2) ||
	!strncmp(requete,"J=",2) ||
	!strncmp(requete,"Y=",2) ||
	!strncmp(requete,"Y>",2) ||
	!strncmp(requete,"Y<",2) ||
	!strncmp(requete,"O=",2) ) {
	if(*requete=='O') {
		m3=getfreelist(*numliste);
		if(m3== -1) { poserr=dnom; return 2; }
		dm3=defbitlist+m3*lenw;
		numerr=orgseq(nom,dr,dm3);
		}
	else if(*requete=='Y') {
		m3=getfreelist(*numliste);
		if(m3== -1) { poserr=dnom; return 2; }
		dm3=defbitlist+m3*lenw;
		numerr=yeaseq(nom,*(requete+1),dr,dm3);
		}
	else {
		static char choix[]="SMJYTO";
		m3= ( strchr(choix,*requete)-choix );
		numerr=smjseq(nom,m3,dr);
		}
	if(numerr!=1) { poserr=dnom; return 6; }
	deflocus[*numliste] = ( (*requete)!='T' || nbrf || swissprot);
	defgenre[*numliste]='S';
	}
else if(!strncmp(requete,"R=",2)) {
	numerr=bibseq(nom,dr);
	if(numerr!=1) {
		poserr=dnom; return 7;
		}
	deflocus[*numliste]=TRUE;
	defgenre[*numliste]='S';
	}
else if(!strncmp(requete,"AU=",3)) {
	numerr=autseq(nom,dr);
	if(numerr!=1) {
		poserr=dnom; return 8;
		}
	deflocus[*numliste]=TRUE;
	defgenre[*numliste]='S';
	}
else if(!strncmp(requete,"N=",2)) {
	deflocus[*numliste]=TRUE;
	numerr=isenum_flou(nom,dr,&deflocus[*numliste]);
	if(numerr==0) {
		poserr=dnom; return 9;
		}
	if(numerr>0) bit1(dr,numerr);
	defgenre[*numliste]='S';
	}
else if(!strncmp(requete,"AC=",3)) {
	numerr=accseq(nom,dr);
	if(numerr!=1) {
		poserr=dnom; return 10;
		}
	deflocus[*numliste]=TRUE;
	defgenre[*numliste]='S';
	}
else if(!strncmp(requete,"F=",2) || !strncmp(requete,"FA=",3) ||
	!strncmp(requete,"FS=",3) || !strncmp(requete,"FK=",3) ) {
	char *p;
	int num;
	FILE *fich;
	p=debutsrequetes[1]+(dnom-debutsrequetes[0]);
	memcpy(nom,p,lnom); nom[lnom]=0;
	fich=fopen(nom,"r");
	if(fich==NULL) {
		poserr=dnom; return 11;
		}
	memset(dr,0,lenw*sizeof(int));
	if(*(requete+1)=='A') {
		if( fileaccnums(fich,dr,&deflocus[*numliste]) != 0 ) return 11;
		}
	else 	{
		deflocus[*numliste]=TRUE;
		while( (p=fgets(nom,80,fich)) != NULL) {
			*(p+strlen(p)-1)=0;
			while(*p!=0) { *p=toupper(*p); p++; }
			if(*(requete+1)=='S') {
				if((num=iknum(nom,kspec))!=0) bit1(dr,num);
				}
			else if(*(requete+1)=='K') {
				if((num=iknum(nom,kkey))!=0) bit1(dr,num);
				}
			else if((num=isenum(nom))!=0) {
				bit1(dr,num);
				if(!testbit(defbitlist+lenw,num)) 
					deflocus[*numliste]=FALSE;
				}
			}
		}
	fclose(fich);
	if(*(requete+1)=='S')
		defgenre[*numliste]='E';
	else if(*(requete+1)=='K')
		defgenre[*numliste]='K';
	else
		defgenre[*numliste]='S';
	}
else 	{
/* previous list */
	int trouve=0;
	for(i=1;i<tlist;i++) 
		if(defoccup[i] && 
		   deflnames[i] != NULL &&
		   strcmptrail(nom,lnom,deflnames[i],strlen(deflnames[i]))==0) {
			trouve=1;
			break;
		}
	if(!trouve) {
		poserr=dnom; return 12;
		}
	memcpy(dr,defbitlist+i*lenw,lenw*sizeof(int));
	defgenre[*numliste]=defgenre[i];
	deflocus[*numliste]=deflocus[i];
	}
defoccup[*numliste]=TRUE;
return 0;
} /* end of operand */


static int is_operateur(char *requete, char *pos, size_t l_req)
{
int i, loper=strlen(operateur);
for(i=0; i<loper; i+=3) {
	if(strncmp(pos, operateur+i, 2) == 0) {
		if( (pos == requete || *(pos-1) == ' ' || *(pos-1) == ')' ||
			 *(pos-1) == '"') && 
		   	( pos + 2 >= requete + l_req || *(pos+2) == ' ' ||
			 	*(pos+2) == '(' ) || *(pos+2) == '"' )
			return TRUE;
		}
	}
return FALSE;
}

static void addquotes(char *requete , char *lcrequete )
/* 
pour ajouter des " autour des operandes de requete, et en parallele ds lcrequete
*/
{
int pos, lreq;
char *p, *done=NULL;
/* enlever les blancs initiaux et terminaux */
lreq=strlen(requete);
while(lreq>0 && requete[lreq-1]==' ') lreq--;
requete[lreq]=0;
p=requete;
while(*p==' ') p++;
if(p>requete) {
	memdelta(p,lreq+1,-(p-requete)); lreq-=(p-requete);
	}
while (*requete=='(' && nextpar(requete,requete+lreq-1)==requete+lreq-1) {
	/* supprimer les () englobantes totales */
	memdelta(requete+1,lreq-2,-1); lreq-=2; requete[lreq]=0;
	}
/* traiter la () initiale */
if(*requete=='(' && nextpar(requete,requete+lreq-1)!=NULL ) {
	p=prepaddquotes(requete,lcrequete,requete,&lreq);
	pos=(p-requete);
	done=p;
	}
else
	pos= -1;
while(++pos<lreq) {
	if(requete[pos]==' ') continue;
	if(pos==lreq-1) {
/* ajout de " en fin  */
		pos++; lreq++;
		requete[pos]='"';
		lcrequete[pos]='"';
		requete[lreq]=0;
		lcrequete[lreq]=0;
		}
/* rechercher les operateurs isoles ou entoures par ) ou ( */
	else if( is_operateur(requete, requete + pos, lreq) ) {
/* recherche a gauche de l'operateur */
		p=requete+pos-1;
		while(p>=requete && *p==' ') p--;
		if((done==NULL || p>done) && (p==requete || 
				(p-1>=requete && 
				!is_operateur(requete, p - 1, lreq)))){
			*(p+1)='"';
			*(lcrequete+(p-requete)+1)='"';
			done = p + 1;
			}
/* recherche a droite de l'operateur */
		p=requete+pos+2;
		while(p < requete+lreq && *p==' ') p++;
		if( p < requete+lreq && *p != '(' && ( p + 1 >= requete+lreq || 
				!is_operateur(requete, p, lreq) ) ) {
			*(p-1)='"';
			*(lcrequete+(p-requete)-1)='"';
			done = p - 1;
			}
		else if(p<requete+lreq && *p=='(') {
/* traiter recursivement le contenu de () apres un operateur */
			p=prepaddquotes(requete,lcrequete,p,&lreq);
			if(p==NULL) continue;
			pos=p-requete;
			done=p;
			}
		}
	}
/* ajout de " en debut si besoin: ni ( ni operateur */
pos=0;
while(requete[pos]==' ') pos++;
if( requete[pos]!='(' && !is_operateur(requete, requete + pos, lreq) ) {
	if(pos==0) {
		memdelta(requete,lreq+1,1);
		memdelta(lcrequete,lreq+1,1);
		pos++;
		}
	requete[pos-1]='"';
	lcrequete[pos-1]='"';
	}
}

static char *prepaddquotes(char *requete, char *lcrequete, char *pospar,
	 int *lreq)
{
/* traiter recursivement le contenu de () */
char *aux, *lcaux, *finpar;
int laux, laux2, delta, pos;
pos=pospar-requete;
finpar=nextpar(requete+pos,requete+ *lreq-1);
if(finpar==NULL) return NULL;
laux=finpar-(requete+pos)-1;
laux2=2*laux; if(laux2<laux+10)laux2=laux+10;
aux=mycalloc(laux2,1);
lcaux=mycalloc(laux2,1);
memcpy(aux,requete+pos+1,laux); aux[laux]=0;
memcpy(lcaux,lcrequete+pos+1,laux);
lcaux[laux]=0;
addquotes(aux,lcaux);
laux2=strlen(aux);
delta=laux2-laux;
if(delta!=0){
	memdelta(finpar,(requete+ *lreq)-finpar,
		delta);
	memdelta( lcrequete+(finpar-requete) ,
		(requete+ *lreq)-finpar,delta);
	}
memcpy(requete+pos+1,aux,laux2);
memcpy(lcrequete+pos+1,lcaux,laux2);
*lreq += delta;
free(aux); free(lcaux);
return finpar+delta;
}
