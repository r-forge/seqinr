#include "dir_acnuc.h"
#include <ctype.h>

#define plus_grand_int(a,b)  ( (a) >= (b) ? (a) : (b) )
#define plus_petit_int(a,b)  ( (a) <= (b) ? (a) : (b) )

/* prototypes of included functions */
DIR_FILE *open_access_file(char *fname, char *mode);
void dir_acnucopen(char *db_access);/*use RO, WP or WA to control file access*/
void acnucopen(void);
void dir_acnucclose(void);
void readacc(int recnum);
int read_first_rec(DIR_FILE *fp, int *endsort);
void padtosize(char *pname, char *name, int length);
int fcode(DIR_FILE *fp,char *search, int lcompar);
int hashmn(char *name);
int hasnum(char *name);
int isenum(char *name);
int iknum(char *name, DIR_FILE *fp);
void dir_readerr(DIR_FILE *fich, int recnum);
void gcgini(void);
void goffset(int point,int *div,int *offset);
int poffset(int div, int offset);
int is_rdp_residue(int c);
static int rdnuc(off_t offset, int div, int mlong);
int gfrag(int nsub,int first,int lfrag,char *dseq);
char complementer_base(char nucl);
void complementer_seq(char *deb_ch, int l);
static void gcgbinseq(unsigned char *seq, int length);
static void skipseq(FILE *fich, int length, int codage, off_t *position);
static void prepsq(unsigned *pnuc, int div);
void gcgors(char *type_fic, int div, int stop_if_error);
int gsnuml(char *name,  int *length, int *frame, int *gencode);
int strcmptrail(char *s1, int l1, char *s2, int l2);
char *short_descr(int seqnum, char *text, int maxlen);
char *short_descr_p(int seqnum, char *text, int maxlen);
void simpleopen(void); 
char *get_code_descr(int code);
int calc_codon_number(char *codon);
char codaa(char *codon, int code);
int get_ncbi_gc_number(int gc);
int get_acnuc_gc_number(int ncbi_gc);
char translate_init_codon(int numseq, int gc, int debut_codon /* 1, 2, or 3 */);
char init_codon_to_aa(char *codon, int gc);
void compact(char *chaine);
int trim_key(char *name); /* remove trailing spaces */
void seq_to_annots(int numseq, long *faddr, int *div);
char *read_annots(long faddr, int div);
char *next_annots(long *pfaddr);
void majuscules(char *name);
char *check_acnuc_gcgacnuc(void);
char *translate_cds(int seqnum);
int decode_locus_rec(char *line, char **pname, char **pmolec, 
	int *circular, char **pdivision, char **pdate);


/* allocation of global variables (most of them) */
DIR_FILE *kinf, *knuc;  /* pour acnuc ancienne structure */

#define MAXDIV 2000
int max_divisions=MAXDIV;
FILE *divannot[MAXDIV], *divseq[MAXDIV]; /* pour structures plat et gcg */
int annotopened[MAXDIV], seqopened[MAXDIV];
char *gcgname[MAXDIV];
unsigned div_offset[MAXDIV];


int lmot=(8*sizeof(int)),hoffst,hsub,hkwsp,nseq,nbrf,lenbit,lenw,maxa,longa;
int nbmrfa,flat_format,gcgcod,unixos,embl,genbank,swissprot,big_annots;
int use_div_sizes; /* flag TRUE iff size-based format of pointers to divisions*/
DIR_FILE *ksub,*kloc,*kkey,*kspec,*kbib,*kacc,*ktxt,*ksmj,
	*kext,*kaut,*kshrt,*klng;
struct rsub *psub;
struct rloc *ploc;
struct rkey *pkey;
struct rspec *pspec;
struct rbib *pbib;
struct racc *pacc;
struct rsmj *psmj;
struct rext *pext;
struct raut *paut;
struct rshrt *pshrt;
struct rlng *plng;
struct rinfo *pinfo;
char ptxt[lrtxt];

/* echange entre gfrag et gcgini */
int gfrag_premier;


/* data for handling sequences cut in pieces by gcg */
#define GCG_SLICE 100000  /* size of pieces to cut very long gcg sequences */
#define MAX_GCG_SLICES 50 /* max allowed number of 100 Kb pieces */
static unsigned int slice_number, slice_start[MAX_GCG_SLICES];

/* allocate here the nucbuf buffer to hold the currently read sequence */
/* larger piece ever read: GCG_SLICE */
char nucbuf[GCG_SLICE+1]; 

/* pointer to function to say what char is an allowed residue */
int (*is_residue)(int); 


static DIR_FILE *myopen(char *fname,char *mode,size_t rsize)
{
DIR_FILE *fich;
char context[40];
#ifdef vms
strcpy(context,"BANK");
#else
strcpy(context,"acnuc");
#endif
if( ( fich=dir_open(fname,context,mode,rsize,8*1024) )==NULL){
	fprintf(stderr,"Trouble opening file %s\n",fname);
	exit(ERREUR);
	}
return fich;
}

int acc_length_on_disk; /* taille chaine ACCESSION sur le disque */
int lracc_on_disk; /* taille records du fichier ACCESSION sur le disque */

DIR_FILE *open_access_file(char *fname, char *mode)
{
DIR_FILE *fich;
int totsmj, num, new_format = FALSE;
totsmj = read_first_rec(ksmj, NULL);
for(num = totsmj; num >= 2; num--) {
	readsmj(num);
	if(strncmp(psmj->name, "07ACCESSION ", 12) == 0) {
	new_format = TRUE;
		break;
		}
	}
if( !new_format ) { /* ancien format */
	acc_length_on_disk = 6;
	lracc_on_disk = 14;
	}
else	{ /* nouveau format */
	acc_length_on_disk = ACC_LENGTH;
/* taille sur disque = max( taille memoire , taille 1er record ) */
	lracc_on_disk = 6 + 2*sizeof(int); /* total + SORTED + fin_sort */
	if( lracc_on_disk < sizeof(struct racc) )
		lracc_on_disk = sizeof(struct racc);
	}

#if defined(vms) || defined(__alpha) 
/* machines avec records multiples de 4 octets en F77 */
lracc_on_disk = ( (lracc_on_disk + 3)/4 ) * 4;
#endif

fich = myopen(fname, mode, lracc_on_disk);
return fich;
}


void dir_acnucopen(char *db_access) /*use RO, WP or WA to control file access*/
{
int recn, tot, sorted;
char mode[4];

if(!strcmp(db_access,"RO") )
	strcpy(mode,"r");
else
	strcpy(mode,"r+");
#ifdef vms
	unixos=0;
#else
	unixos=1;
#endif
ksub=myopen("SUBSEQ",mode,lrsub);  psub = (struct rsub *)malloc(lrsub);
nseq=read_first_rec(ksub,NULL);
kloc=myopen("LOCUS",mode,lrloc);  ploc = (struct rloc *)malloc(lrloc);
kkey=myopen("KEYWORDS",mode,lrkey);  pkey = (struct rkey *)malloc(lrkey);
tot=read_first_rec(kkey,NULL);
kspec=myopen("SPECIES",mode,lrspec);  pspec = (struct rspec *)malloc(lrspec);
maxa=read_first_rec(kspec,NULL);
maxa=(tot > maxa ? tot : maxa);
kbib=myopen("BIBLIO",mode,lrbib);  pbib = (struct rbib *)malloc(lrbib);
ktxt=myopen("TEXT",mode,lrtxt); 
ksmj=myopen("SMJYT",mode,lrsmj);  psmj = (struct rsmj *)malloc(lrsmj);
kacc=open_access_file("ACCESS",mode);  
pacc = (struct racc *)calloc(1, lracc_on_disk);/*pour pouvoir faire dir_read*/
tot=read_first_rec(ksmj,&recn);  sorted = (recn==tot);
nbrf=swissprot=FALSE;
for(recn=2;recn<=tot;recn++) {
	readsmj(recn);
	if(!strncmp(psmj->name,"01PROTEIN ",10)) { nbrf = TRUE; break; }
	if(!strncmp(psmj->name,"01PRT ", 6)) { swissprot = TRUE; break; }
	if(sorted && strncmp(psmj->name,"01",2)>0 ) break;
	}
kext = NULL;
if( !( nbrf || swissprot) ) {
	kext=myopen("EXTRACT",mode,lrext);  pext = (struct rext *)malloc(lrext);
	}
kaut=myopen("AUTHOR",mode,lraut);  paut= (struct raut *)malloc(lraut);
kshrt=myopen("SHORTL",mode,lrshrt);  pshrt = (struct rshrt *)malloc(lrshrt);
readshrt(2);
if(pshrt->val < 0) { /* new hashing format */
	hoffst=2;
	hsub=abs(pshrt->val);  hkwsp=abs(pshrt->next);
	}
else	{ /* old hashing format */
	hoffst=1;
	hsub=9973; hkwsp=1999;
	}
klng=myopen("LONGL",mode,lrlng);  plng = (struct rlng*)malloc(lrlng);
pinfo=(struct rinfo *)malloc(lrinfo);

gcgini();
longa=(maxa-1)/lmot+1;
lenbit=(nseq>maxa? nseq : maxa);
lenw=(lenbit-1)/lmot+1;
return;
} /* end of dir_acnucopen */

void acnucopen(void)
{
dir_acnucopen("RO");
}


void dir_acnucclose(void)
{
int div;

dir_close(kshrt);
dir_close(ksub);
dir_close(kloc);
if( !( nbrf || swissprot) )dir_close(kext);
for(div=0; div<=divisions; div++) {
	if( annotopened[div] ) fclose(divannot[div]);
	if( seqopened[div] ) fclose(divseq[div]);
	}
if(ksmj != NULL) {
	dir_close(ksmj);
	dir_close(kkey);
	dir_close(kspec);
	dir_close(kbib);
	dir_close(kacc);
	dir_close(ktxt);
	dir_close(kaut);
	dir_close(klng);
	}
return;
}


void readacc(int recnum)
{
/* ici structure en memoire ne colle pas celle sur le disque, donc il faut
copier ce qui est lu vers la structure en memoire: pacc
*/
int numlu;
char *point;
point = (char *)dir_read_buff(kacc,recnum,&numlu);
if(numlu < 1) dir_readerr(kacc,recnum);
memcpy(pacc->name, point, acc_length_on_disk);
memcpy(&(pacc->plsub), point + acc_length_on_disk, sizeof(int));
if(acc_length_on_disk < ACC_LENGTH)
	memset(pacc->name + acc_length_on_disk, ' ',
		 ACC_LENGTH - acc_length_on_disk);
}


int read_first_rec(DIR_FILE *fp, int *endsort)
/* lecture complete du 1er rec d'un fichier 
au retour:
	valeur rendue =  nbre total de recs
	*endsort = total si SORTED; fin triee si 1/2SOR; 0 sinon
		rien n'est retourne si endsort est envoye NULL
*/
{
char *pos;
int total;
pos=(char *)dir_read_buff(fp,1,&total);
if(total <= 0)dir_readerr(fp,1);
total= *((int *)pos);
if(endsort != NULL) {
	if(!strncmp(pos+sizeof(int),"SORTED",6))
		*endsort= total;
	else if (!strncmp(pos+sizeof(int),"1/2SOR",6))
		memcpy( (char *)endsort, pos+sizeof(int)+6, sizeof(int) );
	else
		*endsort=0;
	}
return total;
}



void padtosize(char *pname, char *name, int length)
{
int i;
strncpy(pname,name,length);
pname[length]=0;
for(i=strlen(pname);i<length;i++) pname[i]=' ';
}


/* fin partie triee des fichiers kaut, kbib, kacc, ksmj, ksub */
int fcode_nsortd[5]={0,0,0,0,0};

int fcode(DIR_FILE *fp,char *search, int lcompar)
/*
marche avec fichiers kaut, kbib, kacc, ksmj, et ksub
en lecture et aussi en allongement de ces fichiers 
fp: fichier traite
search: cle recherchee dans le fichier
lcompar: longueur de la cle
*/
{
static int first[5]={1,1,1,1,0};
int lastrc;
static char target[41];
char *pbuff;
int numf,retval,dernier,premier,ord;
if( lcompar == 0)lcompar=strlen(search);
if(lcompar==0) return 0;
padtosize(target,search,lcompar);
dernier=0; while(dernier<lcompar) {
	target[dernier]= toupper(target[dernier]);
	dernier++;
	}
if (fp==kaut) { numf=0; pbuff=(char *)paut;  }
else if (fp==kbib) { numf=1; pbuff=(char *)pbib;  }
else if (fp==kacc)  { 
	numf=2; pbuff=(char *)pacc;  
	if(lcompar > acc_length_on_disk) lcompar = acc_length_on_disk;
	}
else if (fp==ksmj)  { numf=3; pbuff=(char *)psmj;  }
else if (fp==ksub) { numf=4; pbuff=(char *)psub; }

lastrc=read_first_rec(fp,fcode_nsortd+numf);
if(lastrc == 1) return 0;
if( first[numf] && fcode_nsortd[numf] == 0 ) {
	char pre[40];
	first[numf]=0;
	dir_read(fp,2,1,pbuff);
	strncpy(pre,pbuff,lcompar);
	retval=0;
	if(!strncmp(pre,target,lcompar)) retval = 2;
	for (dernier=2;dernier<=lastrc-1;dernier++) {
		dir_read(fp,dernier+1,1,pbuff);
		if(!strncmp(pbuff,target,lcompar)) retval = dernier+1;
		if(strncmp(pbuff,pre,lcompar) < 0) break;
		strncpy(pre,pbuff,lcompar);
		}
	fcode_nsortd[numf]=dernier;
	if(retval) goto lastop;
	goto sequential_search;
	}

/* recherche dichotomique dans partie triee */
premier=2;
dernier=fcode_nsortd[numf];
while(premier <= dernier) {
	retval=(dernier+premier)/2;
	dir_read(fp,retval,1,pbuff);
	if( (ord=strncmp(target,pbuff,lcompar)) == 0) goto lastop;
	else if(ord > 0) premier = retval+1;
	else dernier = retval-1;
	}

sequential_search:
for(retval=fcode_nsortd[numf]+1;retval<=lastrc;retval++) {
	dir_read(fp,retval,1,pbuff);
	if(!strncmp(target,pbuff,lcompar))goto lastop;
	}
return 0;
lastop:
if(fp==kacc) readacc(retval);
else dir_read(fp, retval, 1, pbuff); /* par securite, lire l'enreg trouve */
return retval;
} /* end of fcode */


int hashmn(char *name)
{
int retval=0;
char *pos;
pos=name+16;
while(name < pos) {
	retval += *((int *)name) % hsub; 
	name += sizeof(int); 
	}
return (abs(retval) % hsub ) +1;
}

int hasnum(char *name)
{
int retval=0;
char *pos;
pos=name + WIDTH_KS;
while(name < pos) {
	retval += *((int *)name) % hkwsp; 
	name += sizeof(int); 
	}
return (abs(retval) % hkwsp ) +1;
}

int isenum(char *name)
{
int retval, l;
static int copie[5];
char *pname=(char *)copie;

l=strlen(name); if(l>16) l=16;
memcpy(pname,name,l);
for(retval=0; retval<l; retval++) pname[retval]=toupper(pname[retval]);
while(l<16) { pname[l]=' '; l++; }
retval=hashmn(pname);
readshrt(hoffst+(retval+1)/2);
retval = (retval % 2 ? pshrt->val : pshrt->next);
while( retval ) {
	readsub(retval);
        if( !strncmp(psub->name,pname,16) ) return retval;
        retval = psub->h;
        }
return 0;
}

int iknum(char *name,DIR_FILE *fp)
{
int retval,point, l;
struct rspec *pbuff;
static int copie[11];
char *pname=(char *)copie;

l=strlen(name); if(l > WIDTH_KS) l = WIDTH_KS;
memcpy(pname,name,l);
for(retval=0; retval<l; retval++) pname[retval]=toupper(pname[retval]);
while(l < WIDTH_KS) { pname[l]=' '; l++; }
retval=hasnum(pname);
point=hoffst+(hsub+1)/2;
if(fp == kspec) {
	pbuff = pspec;
	point += (hkwsp+1)/2;
	}
else pbuff = (struct rspec *)pkey;
readshrt(point+(retval+1)/2);
retval = (retval % 2 ? pshrt->val : pshrt->next);
while( retval ) {
	dir_read(fp,retval,1,pbuff);
        if( !strncmp(pbuff->name, pname, WIDTH_KS) ) return retval;
        retval = pbuff->h;
        }
return 0;
}


void dir_readerr(DIR_FILE *fich, int recnum)
{
char text[150];
sprintf(text,"Error while reading in %s at record # %d",fich->filename,recnum);
perror(text);
exit(ERREUR);
}



static int calculer_nb_bit(int nb_div)
{
      int nb_bit = 0;

      do
      {
	nb_div >>= 1;
        ++nb_bit;

      } while(nb_div > 0);
      return(nb_bit);
}



static int rfquot(void)
{
int rfquota;

#if defined(vms)
#include JPIDEF
	int status;
	int lib$getjpi();

	status=lib$getjpi(&JPI$_FILCNT,0,0,&rfquota,0,0);
#elif defined(unix)
/* under unix, file descriptors are allocated serially after the last
currently used one, thus closing a file which was not the last one opened
does not really help!!
Here, we increase the max authorized file descriptors to its max possible value
and hope that it is enough to open all divisions. 
Empirically it seems reasonable.
*/
#ifdef sun
/* a strange warning on sun computers who do not like resource.h!
The necessary part is reproduced here without the cause of the warning.
*/
#define	RLIMIT_NOFILE	5		/* file descriptors */
	struct rlimit {
		unsigned long rlim_cur, rlim_max;
		};
	extern int setrlimit(int , void *);
	extern int getrlimit(int , void *);
#else
#include <sys/types.h>
#include <sys/time.h>
#include <sys/resource.h>
#endif
#ifdef RLIMIT_NOFILE  /* some systems (sgi) do not know RLIMIT_NOFILE limit */
	struct rlimit rlp;
	getrlimit(RLIMIT_NOFILE,&rlp);
	rlp.rlim_cur = rlp.rlim_max;
	setrlimit(RLIMIT_NOFILE,&rlp);
#endif
	rfquota = 200; /* semble que s'arrete � 240 sur sun */
	
#elif defined(__INTEL__)
	rfquota = FOPEN_MAX - 5;	
#else
	rfquota = 1000;
#endif
return rfquota;
}


/* allocate here: 
divisions	current total # of divisions (from 0) 
maxgcgfiles	max # of files that can be opened for divisions
*/
int divisions, maxgcgfiles;

void gcgini(void)
{
/* determine what format for the database: acnuc, flat or gcg
and identifies the names of divisions */

  int nb_enrtot ; /* nbre total d'enregistrements de smjyt */
  int num, suivant, debut, rank;
  unsigned taille, previous;
  char *p;
 
   divisions = -1 ; 
   nbmrfa = lmot; /* cas de l'ancienne structure: ABANDONNE! */
   genbank = embl = big_annots = FALSE;
   use_div_sizes=FALSE;  /* structure basee sur la taille des divisions */
   gfrag_premier=TRUE; /* pour que prochain appel a gfrag soit bien initialise*/

/* determiner format des addresses des annotations et sequences */
big_annots = ( fcode(ksmj, "07BIG_ANNOTS ", sizeof(psmj->name)) != 0 );
use_div_sizes = !big_annots;

if ( fcode(ksmj, "07ALLOW_PUNCTUATION ", sizeof(psmj->name)) != 0 ) 
	is_residue = is_rdp_residue;
else
	is_residue = isalpha;


/* PARSE FILE SMJYT FOR 06GCG or 06FLT GIVING FILE NAMES OF DIVISIONS */
nb_enrtot = read_first_rec(ksmj, NULL);
for(num = 2; num <= nb_enrtot; num++) {
	readsmj(num);
	if(!strncmp(psmj->name,"04ID  ",6) ) embl = TRUE;
	if(!strncmp(psmj->name,"04LOCUS ",8) ) genbank = TRUE;
	if(strncmp(psmj->name,"06",2) )continue;
	divisions++;
	if(divisions >= MAXDIV) {
		fprintf(stderr,"Error: too many divisions!\n");
		exit(ERREUR);
		}
	flat_format = ( strncmp(&(psmj->name[2]),"FLT",3) == 0 );
	p = NULL;
	if( psmj->libel != 0) {
		readtxt(psmj->libel);
		p = strstr(ptxt,"rank:");
		}
	if(p != NULL) {
		sscanf(p + 5, "%d", &rank);
		if(use_div_sizes) {
			p = strstr(ptxt,"size:")+5;
			sscanf(p,"%u",&taille);
			div_offset[rank] = taille;
			}
		}
	else	{
#ifdef unix
		fprintf(stderr, "ERROR: missing rank of division %.20s\n",
				psmj->name);
#endif
		exit(ERREUR);
		}
	gcgname[rank] = (char *)malloc(16);
	memcpy(gcgname[rank], &(psmj->name[5]), 15);
	gcgname[rank][15] = 0; trim_key(gcgname[rank]);
	}

if( use_div_sizes) {
/* calcul du tableau div_offset: taille cumulee des divisions precedentes */
	previous=0;
	for(num=0; num<=divisions; num++) {
		taille= previous + (div_offset[num]/1024+1)*1024;
		div_offset[num]=previous;
		previous=taille;
		}
	}

/* ALL .REF AND .SEQ FILES ARE INITIALLY CLOSED */
for(num=0; num<MAXDIV; num++ ) {
	annotopened[num]=seqopened[num]=0;
	}

/*  GET THE MAX NUMBER OF .DAT/.REF/.SEQ FILES THE PROGR MAY OPEN */
maxgcgfiles=rfquot() - 1;
if(maxgcgfiles<2) {
	fprintf(stderr,"I cannot open enough files to work.\n");
	exit(ERREUR);
	}

if(divisions == -1) {
#ifdef unix
	fprintf(stderr, "no division found\n");
#endif
	exit(ERREUR);
	}

nbmrfa = 0; /* plus vraiment utile */

return;
} /* end of gcgini */



void goffset(int point,int *div,int *offset)
{
/*
point	(input) mix of division + offset within division
div	(output) the division encoded in point
offset	(output) the offset encoded in point
*/
 
/* par taille cumulee des fichiers */
for(*div=divisions; *div>0; (*div)--)
		if( (unsigned) point >= div_offset[*div] ) break;
*(unsigned *)offset = (unsigned) point - div_offset[*div];
}


int poffset(int div, int offset)
/*
div	(input) a division number
offset	(input) an offset within the division
return value	a mix of division + offset in a single integer
*/
{
/* par taille cumulee des fichiers */
return (int) ( (unsigned) offset + div_offset[div] );
}


int is_rdp_residue(int c)
{
return isalpha(c) || ispunct(c);
}



static int rdnuc(off_t offset, int div, int mlong)
/* NEVER CALLED DIRECTLY, ONLY THROUGH GFRAG 		*/
/* 
reads mlong bases from pointer to sequence offset + div
and puts that in global variable nucbuf
return FALSE if ok, TRUE if error
*/
{
FILE *fichier;
char ligne[100], *p;
int tot,laux,deb,l;

if(flat_format) {
	if( !annotopened[div] ) gcgors("inf",div,TRUE);
	fichier= divannot[div];
	}
else	{
	if( !seqopened[div] ) gcgors("nuc",div,TRUE);
	fichier= divseq[div];
	}
fseeko(fichier, offset, SEEK_SET);
if(flat_format) {  /* format plat */
        l = 0;
        if(nbrf){ /* pour NBRF, on lit aussi la ponctuation */
		laux = 68;
		deb  = 9;
		do	{
			if( fgets(ligne,100,fichier) == NULL) return TRUE;
			tot=laux-deb+1;
			memcpy(nucbuf+l,ligne+deb-1,tot);
			l += tot;
			if(l >= mlong) return FALSE;
			}
		while(TRUE);
		}
	else	{
		do 	{
			if( fgets(ligne, sizeof(ligne), fichier) == NULL) 
				return TRUE;
			p=ligne;
			while(*p != 0 && l<mlong) {
				if(is_residue(*p)) nucbuf[l++]= *p;
				p++;
				}
			}
		while(l < mlong);
		}
      }
 else  /* cas format gcg */
      {
	l=mlong;
	if(gcgcod == 1)
		l=(mlong-1)/4 + 1; 
#ifdef vms
	{ int pos, lu;
	for(pos=0; pos<l; pos+=MAXRECLENGTH) {
                lu=l-pos;
                if(lu>MAXRECLENGTH) lu=MAXRECLENGTH;
		fread(nucbuf+pos,1,lu,fichier);
		if(pos+MAXRECLENGTH<l) {
			offset += lu+2; 
			fseeko(fichier, offset, SEEK_SET);
			}
		}
	}
#else
	fread(nucbuf,1,l,fichier);
#endif
	if(gcgcod == 1)
		gcgbinseq((unsigned char *)nucbuf,mlong);
      }
return FALSE;
 } /* ---- fin rdnuc ---- */


/* valeur de SIMEXT_NFMAX ici pas importante, valeur reelle ds extract.c */
#define SIMEXT_NFMAX 10 
struct simext {
	int nfmax;  /* maximum # of fragments in a virtual subseq */
	int newsim; /* TRUE when a new virtual subseq has just been created*/
	int pinfme; /* pointer to annots of parent of virtual subseq */
	int siminf; /* pointer to annots of virtual subseq */
	int div;    /* div of these annots */
	int nfrags; /* # of fragments in virtual subseq */
	int valext[SIMEXT_NFMAX][4];
/* valext(.,0)=beginning of fragment */
/* valext(.,1)=end of fragment	*/
/* valext(.,2)=pointer to nucleot for that fragment */
/* valext(.,3)=# of parent seq of that fragment	*/
	} *psimext=NULL;


/* ---------------------------------------------------------------- */
/*  fonction gfrag avec prise en compte des nouvelles structures    */
/*	recoit le numero d'enregistrement de la sequence ds subseq. */
/*	       le numero du premier nucleotide a extraire dans      */
/*	       la sequence.					    */
/*	       le nombre de nucleotides a extraire		    */
/*	       un pointeur sur un tableau de carateres		    */
/*								    */
/*	traitement						    */
/*		transfert lfrag nucleotides a partir de l'adresse   */
/*		pointee par dseq., ou moins si la sequence est plus */
/*              courte.        					    */
/*								    */
/* 	retourne lfrag >0 nombre de nucl. effectivement lus ds seq. */
/*		       =0 si first <= 0 ou lfrag <= 0	            */
/* ---------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
/* Le cas nsub = -1 correspond a l'extraction d'une sequence dont le decoupage*/
/* est decrit dans les features mais qui ne correspond pas a une sous-sequence*/
/* Une fonction appelee avant gfrag analyse les features et cree une pseudo   */
/* structure extract. L'ensemble des parametres necessaires sont enregistres  */
/* dans des variables globales.                  			      */
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------  */
/* NOM DES VARIABLES MODIFIEES PAR RAPPORT AU FORTRAN   		*/
/* long devient mlong (long est reserve en c)				*/
/* ploc devient mloc car ploc est une var globale de acnuc.h		*/
/* pnuc idem mnuc (idem)						*/
/*	(penser a des noms du type enr_nuc car ce sont des pointeurs) 	*/
/* pext devient mext							*/
/* le buffer AUX (et ses derives laux, vlaux ...) devient nucbuf        */
/* -------------------------------------------------------------------  */


int gfrag(int nsub, int first, int lfrag, char *dseq)
/*
 int nsub  ;  numero d'enr. ds subseq               
 int first ;  num. du premier nucleotide            
 int lfrag ;  nbre de nucleotides                   
 char *dseq;  pointeur sur un vecteur de caracteres 
 valeur rendue: nbre de nucl. lus effectivement
*/ 
{

    /* ------------------- DECLARATIONS ------------------------------  */
    /* EXPLANATIONS:							*/
    /* deb: 5'end of exon in mother coordinates				*/
    /* fin: 3'end ...							*/
    /* df : 5'end of exon in subseq coordinates				*/
    /* ff : 3'end ...							*/
    /* dm : 5'end of exon region to be read in mother coordinates	*/
    /* fm : 3'end ...							*/
    /* nucbuf: buffer							*/
    /* r  : Nucleot record # present in buffer				*/
    /* fr : 3'end of fragment present in buffer in mother coordinates	*/
    /* da : 5'end of buffer region to be read in buffer coordinates	*/
    /* fa : 3'end ...							*/
    /* daf: 5'end of exon frag. present in buffer in subseq coordinates */
    /* faf: 3'end ...							*/
    /* maxseq: # of residues maintained in nucbuf buffer 		*/
    /* lnucbuf: # of bytes read from sequence file to fill the buffer	*/
    /*       differs from MAXSEQ for GCG format because sequence may be */
    /*	     2-bit coded and for NBRF because sequence punctuation is   */
    /*       not kept in the buffer.					*/
    /* recoff: offset between the beginning in file of 2 buffer-full    */
    /*       residues							*/
    /* comprecoff: true until one full base line was read               */

    int ls,dm,fm, mloc,da,fa,la;
    static int oldsub,compl,oldpex,mext,df,ff,daf,faf,fr;
    static int deb,fin,mlong; 
    static int new_read,maxseq,recoff,lnucbuf,comprecoff, current_div;
    static unsigned int mnuc, r;


    /* nlles variables */
    int ii = 0; /*cpteur de boucle*/
    unsigned int mrfar;
    int vfr, vlnucbuf;


    if(gfrag_premier) 
    {
/* si c'est le prem. appel a la fct initialisations de diff var.*/
	gfrag_premier = FALSE;
	comprecoff = 0;
	oldsub  = 0;

	if(flat_format)  
	/* FOR FLAT FILE STRUCTURE */
	{
		if(nbrf)
                /* for flat nbrf, sequence records are of variable width,  */
                /* thus each sequence must be read completely by one rdnuc */
                /* call because file address computation does not work	   */
		{
			maxseq = 40000; 
			lnucbuf   = 2*maxseq; /* keep < sizeof(nucbuf) */
			recoff = 1;
		}

		else
		{
			lnucbuf   = 60 * ( (sizeof(nucbuf) - 1) / 60);
			maxseq = lnucbuf;
			comprecoff = 1;
	 	}
	}

	else /* FOR GCG STRUCTURE */
	{
		maxseq= GCG_SLICE;
      	}
	
     } /* fin cas premier appel */
	      
     ls = 0;

    if(first <= 0)
    {
	lfrag = 0; /* 1er cas d'erreur: first <= 0; on retourne lfrag = 0 */
	goto l1000;
    }


    else if( (nsub != oldsub) || ( (nsub== -1) && psimext->newsim)) 
     goto l10;  /* = nlle seq car nl enr de subseq ou seq. virtuelle nouvelle */

    lfrag = plus_petit_int(lfrag,mlong-first+1);
    if( lfrag <= 0 ) /* lg negative ou nulle */
    {
	lfrag = 0;
	goto l1000; /* on sort car erreur */
    }
    else if( (first >= daf) && (first <= faf) )
    {
		dm = ( plus_grand_int(first,df)-df ) * compl + deb;
		fm = ( plus_petit_int(first+lfrag-1,ff)-df) * compl + deb;
		goto l215;
    }

    else if( (first >= df) && (first <= ff) )
		goto l202;
	
    else if( first > ff)
		goto l201;

    else goto l200; /* attention aux else et else if */



l10:
	
	/* cas des sequences virtuelles */

	if(nsub == -1)  
	{
	   mlong = 0;
	   for(ii = 0; ii< psimext->nfrags;ii++)
	   {
		mlong  = mlong + 
		   abs(psimext->valext[ii][1] - psimext->valext[ii][0]) + 1;
		oldpex = 0;

	   }

	}
	else /*nl enr donc on lit subseq */
	{
	  readsub(nsub);
	  mlong  = psub->length;
	  oldpex = psub->pext;
	  mloc   = psub->plinf; /* si mere -> loc : si fille -> inf */
	}
        lfrag = plus_petit_int(lfrag,mlong - first +1); 

	if(lfrag <= 0)
	{
		lfrag = 0;
		goto l1000; /* on sort vers l'appelante */ 
	}
	oldsub = nsub;   /* on garde le numero subseq en memoire */

	/* nsub = -1 cas d'une sequence "virtuelle " */
       	if(nsub == -1)
	   psimext->newsim = 0; /*la sequence virtuelle n'est plus nouvelle  */
	if(oldpex > 0 || nsub == -1) /* seq fille on doit passer par extract */
		goto l200;

	readloc(mloc);    /* seq mere donc on lit locus a partir de subseq */
 	mnuc = ploc->pnuc; 			
	if(big_annots) current_div = ploc->div;
	else	{ int tmp;
		goffset( mnuc, &current_div, &tmp);
		mnuc = tmp;
		}
        df = deb  = 1;
	compl = 1;
	ff = fin = mlong;
	goto l203;
	
l200:
	ff = 0; 
	mext= oldpex;

l201: 
	if( nsub == -1)
	{
		deb  = psimext->valext[mext][0]; 
		fin  =  psimext->valext[mext][1];
		mnuc =  psimext->valext[mext][2];
		if(big_annots) {
			readsub(psimext->valext[mext][3]);
			readloc(psub->plinf);
			current_div = ploc->div;
			}
		else	{ int tmp;
			goffset( mnuc, &current_div, &tmp);
			mnuc = tmp;
			}
		mext =  mext + 1;
		if(mext >= psimext->nfrags)
			mext = 0;	
	}
	else /* on lit extract */
	{
		readext(mext);
		mnuc =  pext->pnuc;
		if(big_annots) {
			readsub(pext->mere); 
			readloc(psub->plinf); 
			current_div = ploc->div;
			}
		else	{ int tmp;
			goffset( mnuc, &current_div, &tmp);
			mnuc = tmp;
			}
		deb  =  pext->deb;
		fin  =  pext->fin;
		mext =  pext->next; /* enr suivant ds extract */
	}
	df = ff+ 1;
	ff = df + abs(fin - deb);
	if( ff < first ) 
		goto l201;
	if(deb <= fin) 
		compl =  1; 
	else 
		compl = -1; /* on doit complementer */



l203:
	if( !flat_format ) /* cas gcg */
		{
		prepsq(&mnuc, current_div);
		if(gcgcod == 1)
			lnucbuf = maxseq / 4;
		else
			lnucbuf = maxseq;
		if(unixos)
			recoff = lnucbuf;
		else
			recoff = lnucbuf + 2;			

		}
	else if(comprecoff) {
		int l,j,nbases = 0;
		pinfo->line[0] = 0;
		read_annots(mnuc, current_div);
		l = strlen(pinfo->line);
		for(j=0; j<l; j++) 
		  if(pinfo->line[j]!=' ' && !isdigit(pinfo->line[j]) ) nbases++;
		if(nbases == 60) {
			off_t fpos;
			fpos = ftello(divannot[current_div]);
			recoff = (lnucbuf/60) * (fpos - mnuc);
			comprecoff = 0;
			}
		else	recoff=0;
		}
l202:
	
	dm = ( plus_grand_int(first,df)-df ) * compl + deb;
	fm = ( plus_petit_int(first+lfrag-1,ff)-df) * compl + deb;

	r  = mnuc +   ( (dm-1) / maxseq );
	fr = (r - mnuc+1) * maxseq;


l210:
/* compute the file address mrfar of beginning of what we want to read */
	if(!flat_format) { /* format GCG is special: may be several pieces */
		int debut_buffer, current_slice, pos_in_slice;
		debut_buffer=(r-mnuc)*maxseq+1;
		if(slice_number > 0) { /* several pieces */
			current_slice=(debut_buffer-1)/GCG_SLICE;
			if(current_slice>slice_number) {
				fprintf(stderr,
				    "Error: part of sequence not found in gcg file\n");
				exit(ERREUR);
				}
			}
		else	current_slice = 0;
		pos_in_slice = debut_buffer - current_slice * GCG_SLICE;
		mrfar= ((pos_in_slice-1)/maxseq)*recoff +
					slice_start[current_slice];
		}
	else
		mrfar = (r - mnuc) * recoff + mnuc;
	if(nbrf && (r != mnuc)) {
		fprintf(stderr,
			"Error cannot read prot. seq. longer than %d\n",maxseq);
		exit(ERREUR);
        	}
/* compute the exact # of residues of the .SEQ record to read */
	vfr = plus_petit_int(fr,plus_grand_int(deb,fin));
	if(nbrf)
		vlnucbuf = 2 * (vfr - (r-mnuc) * maxseq);
	else
		vlnucbuf = vfr - (r-mnuc) * maxseq;

	if( rdnuc( (off_t)mrfar, current_div, vlnucbuf) ) {
		/* rdnuc return TRUE indicates error */
		dseq[0] = 0;
		oldsub = 0;
		return 0;
		}
	new_read = TRUE;

l215:

	da = (plus_grand_int
		(plus_petit_int(dm,fm),fr-maxseq+1) -1 ) % maxseq;	
	da = da + 1;

	fa = (plus_petit_int
		 (plus_grand_int(dm,fm),fr) -1 ) % maxseq;
	fa = fa+1;

	la = fa - da + 1;
	if(nbrf)
	{    
	   for (ii = 1; ii <= la; ii++) 
		dseq[ls+ii-1] = nucbuf[2*(da-1+ii)-1];
	}
	else
	{
	   /* transfert ds dseq */
		memcpy( &dseq[ls] , &nucbuf[da-1] , la );
        }

	if( compl == -1)
		complementer_seq(dseq+ls,la);
		ls +=  la;

	if( ((compl == 1) && (fr < fm)) ||
	    ((compl == -1) && ((fr+1-maxseq) > fm)) )
	{
			r += compl;
			fr += maxseq*compl;
			goto l210;
	}
	else if( ls<lfrag )
		   goto l201;

	else if (new_read) /* pour le prochain appel on change certaines variables */
	{
		if( compl == 1)
		{
			daf = plus_grand_int(deb,fr+1-maxseq)-deb+df;
			faf = plus_petit_int(fin,fr)- deb + df;
		}
		else
		{
			daf = deb - plus_petit_int(fr,deb)+df;
			faf = deb - plus_grand_int(fr+1-maxseq,fin)+df;
		}
		new_read = FALSE;
	}


l1000: /* on sort vers la fct appelante */
	*(dseq+lfrag)=0;
	return(lfrag);
										
} /* fin gfrag */


/* -------------- rend le complementaire d'une base ------------------- */
/* appele par complementer_seq					        */
/* -------------------------------------------------------------------- */

char complementer_base(char nucl)
{
    switch (nucl) {
        case 'a':
        case 'A': return('t');

        case 'c':
        case 'C': return('g');

        case 'g':
        case 'G': return('c');
  
        case 'u':
        case 'U':
        case 't':
        case 'T': return('a');

        case 'r':
        case 'R': return('y');

        case 'y':
        case 'Y': return('r');

	default : return('n');

	}
}
   

/*   ~~~~~~~~~~~~ retourne le complementaire d'une sequence ~~~~~~~~~~~ 
 * recoit l'adresse du debut d'un tableau de caractere et sa longueur
 * inverse et complemente cette sequence
 * prend en compte si c'est un adn ou un arn
 * -------------------------------------------------------------------- */

void complementer_seq(char *deb_ch, int l)
{
    int ii = 0;
    char compl1,compl2;

    for(ii = 0; ii <= (l-1)/2; ii++)
    {
	compl1 = complementer_base(*(deb_ch+ii));

	compl2 = complementer_base(*(deb_ch+l-ii-1));

	*(deb_ch+ii)     = compl2;
	*(deb_ch+l-ii-1) = compl1;
    }	

 
}


static void gcgbinseq(unsigned char *seq, int length)
/*
c to decode gcg encoding of sequences where each nucl is on 2 bits
c seq: (input) the encoded sequence
c      (output) the decoded sequence
c length: (readonly) the sequence length in nucleotides
*/
{
static char alphab[]="ctag";
static char table[4][256];
static int first=1;
int last, pos, pcode, i, j;
unsigned int code;
if(first) {
	int i,j,k,l, rank;
	first=0;
/*
 compute the corrispondence table:
 table(*,code)=the 4 nucleotides corresponding to the coded byte code
*/
	rank= -1;
	for(i=0; i<=3; i++)
		for(j=0; j<=3; j++)
			for(k=0; k<=3; k++)
				for(l=0; l<=3; l++) {
					rank++;
					table[0][rank] = alphab[i];
					table[1][rank] = alphab[j];
					table[2][rank] = alphab[k];
					table[3][rank] = alphab[l];
					}
	}
/* decode the sequence starting from its last coded byte */
last=(length-1)/4+1;
pos=4*last-3;
for(pcode=last-1;pcode>=0;pcode--) {
	code= seq[pcode];
	i= -1;
	for(j=pos-1; j<pos+3; j++) {
		seq[j]=table[++i][code];
		}
	pos -= 4;
	}
return;
} /* end of gcgbinseq */


static void skipseq(FILE *fich, int length, int codage, off_t *position)
{ /* skip one sequence in gcg format */
#ifdef vms
int nbr_lines, offset, reste;
#endif
off_t current;

if(codage == 1) length = (length+3)/4;
current = ftello(fich);
#ifndef vms
	current += length + 1;
	fseeko(fich, current, SEEK_SET);
#else
	nbr_lines=length/MAXRECLENGTH;
	reste=length % MAXRECLENGTH;
	if(reste%2==1) ++reste;
	offset=nbr_lines*(MAXRECLENGTH+2);
	if(reste!=0) offset += reste+2;
	current += offset;
	fseeko(fich, current, SEEK_SET);
	*position += offset;
#endif
}


static void prepsq(unsigned *pnuc, int div)
/*
 to read the >>>> line in file .seq pointed by *pnuc, update *pnuc to 
 beginning of the sequence,
 and return gcgcod=1 if sequence is 2BIT encoded
modified to handle very long seqs cut in pieces of 100000 bases
*/
{
char chaine[100], c;
off_t offset, position;
char name[20], slice_name[20], *p;
int lpiece, lname, cod;

/* goto the >>>> in the .seq file */
if(!seqopened[div]) gcgors("nuc",div,TRUE);
offset = *pnuc;
fseeko(divseq[div], offset, SEEK_SET);
/*  read the >>>> line */
#ifdef vms
fgets(chaine,MINRECSIZE+1,divseq[div]);
#else
fgets(chaine,sizeof(chaine),divseq[div]);
#endif
/* skip the next(title) line */
#ifdef vms
position = offset+MINRECSIZE+MAXDOCSTRING+4;
fseeko(divseq[div], position, SEEK_SET););
#else
do  
	c=getc(divseq[div]);
#ifdef __MWERKS__
while (c != '\r' && c != '\n');
#else
while (c != '\n');
#endif
#endif
/* update the pointer to file .seq */
offset = ftello(divseq[div]);
/* compute gcgcod */
if(strstr(chaine," 2BIT ") != NULL) 
	gcgcod=1;
else
	gcgcod=0;
*pnuc = (unsigned int)offset;

slice_number=0; slice_start[slice_number]= *pnuc;
/* get seq name */
sscanf(chaine+4,"%s",name); 
if( (p = strstr(name, "_0") ) == NULL) return;
lname = p - name + 1;
/* get slice length */
sscanf(strstr(chaine," Len: ")+5 , "%d", &lpiece);
position = offset;
while(TRUE) {
	skipseq(divseq[div],lpiece,gcgcod, &position);
#ifdef vms
	p=fgets(chaine,MINRECSIZE+1,divseq[div]);
#else
	p=fgets(chaine,sizeof(chaine),divseq[div]);
#endif
	if(p==NULL || strncmp(chaine,">>>>",4) != 0) break;
	sscanf(chaine+4,"%s",slice_name);
	if(strncmp(slice_name,name,lname) != 0) break;
#ifdef vms
	offset = position+MINRECSIZE+MAXDOCSTRING+4;
	fseeko(divseq[div], offset, SEEK_SET););
#else
	do  
		c=getc(divseq[div]);
#ifdef __MWERKS__
	while (c != '\r' && c != '\n');
#else
	while (c != '\n');
#endif
#endif
	position = ftello(divseq[div]);
	slice_number++;
	if(slice_number >= MAX_GCG_SLICES) {
		fprintf(stderr,"Error: too many gcg slices\n"); exit(ERREUR);
		}
	slice_start[slice_number] = (unsigned int)position;
	if(strstr(chaine," 2BIT ")!=NULL) cod=1;
	else cod=0;
	if(cod!=gcgcod) {
		fprintf(stderr,
			"Error: cannot process mixed 2BIT and ASCII slices\n"); 
		exit(ERREUR);
		}
	sscanf(strstr(chaine," Len: ")+5 , "%d", &lpiece);
	}
if(lpiece > GCG_SLICE) { /* cas dernier morceau entre 100000 et 110000 */
	slice_number++;
	if(slice_number >= MAX_GCG_SLICES) {
		fprintf(stderr,"Error: too many gcg slices\n"); exit(ERREUR);
		}
	slice_start[slice_number] = slice_start[slice_number - 1] + 
		(gcgcod == 1 ? GCG_SLICE / 4 : GCG_SLICE);
	}
} /* end of prepsq */


void gcgors(char *type_fic, int div, int stop_if_error) 
/*
to open a division file; div is the division number;
type_fic is "inf" for annotation files, or "nuc" for gcg .seq files
*/
{
  FILE *fichier;
  char fname[40], *fullname;
  char env_var[10];
  char ext[5];
  static int orderopen = 0;
  int info, current, i;

#ifdef vms
	strcpy(env_var,"GCGBANK");
#else
	strcpy(env_var,"gcgacnuc");
#endif



  if(flat_format) /*format flat */
  {

	if( embl || nbrf || swissprot ) strcpy(ext, ".dat");
	else strcpy(ext, ".seq");
	info = TRUE;
  }

  else if(strcmp(type_fic, "inf")==0)
  {
  	strcpy(ext, ".ref");
	info = TRUE;
  }

  else if(strcmp(type_fic, "nuc")==0)
  {
  	strcpy(ext, ".seq");
	info = FALSE;
  }

  else
  {
  	fprintf(stderr,"Erreur dans gcgors \n");
	exit(ERREUR);
  }

  strcpy(fname, gcgname[div]);
  strcat(fname, ext);

/* compute the current # of opened gcg files */
current=0;
for(i=0; i<=divisions; i++) {
	if(annotopened[i]) current++;
	if(!flat_format && seqopened[i] ) current++;
	}
while (current >= maxgcgfiles) {
	int sup, ndiv;
/* search for the least recently opened .SEQ or .REF file */
	sup=orderopen+1;
	for(i=0; i<=divisions; i++) {
		if(annotopened[i] && annotopened[i]<sup){
			sup=annotopened[i];
			ndiv=i+1;
			}
		if(!flat_format && seqopened[i] && seqopened[i]<sup) {
			sup=seqopened[i];
			ndiv= -(i+1);
			}
		}
	if(sup==orderopen+1) {
		fichier=NULL;
		goto erreur;
		}
/* close the least used file */
	if(ndiv > 0) {
		fclose(divannot[ndiv-1]);
		annotopened[ndiv-1]=0;
		}
	else	{
		fclose(divseq[-ndiv-1]);
		seqopened[-ndiv-1]=0;
		}
	current--;
	}

/* open the necessary file */
fullname = prepare_env_var(env_var);
strcat(fullname,fname);
fichier = fopen(fullname, 
#ifdef __INTEL__
	"rb"
#else
	"r"
#endif
	);
erreur:
  if(fichier==NULL) {
      if(stop_if_error) {
	fprintf(stderr,"Cannot open the necessary file %s\n",fname);
	exit(ERREUR);
	}
      else return;
      }
/* remember order of opening to then close the file opened most early */
  orderopen++;
  if(info)
  {
  	divannot[div] = fichier;
  	annotopened[div] = orderopen;
  }

  else
  {
  	divseq[div] = fichier;
  	seqopened[div] = orderopen;
  }
} /* end of gcgors */


int gsnuml(char *name,  int *length, int *frame, int *gencode)
/*
name	(input) sequence name, lower or upper case
return value	sequence number in acnuc, or 0 if does not exist
length	(output) sequence length in bases
frame	(output) nothing if NULL pointer or reading frame (0, 1, or 2)
gencode	(output) nothing if NULL pointer or genetic code (0 is universal code)

EXAMPLE:
char name="ecotgp";
int length, seq_num, gen_code;
seq_num=gsnuml(name,&length,NULL,&gen_code);
*/
{
int retval, frame_code;

retval=isenum(name);
if(retval==0)return 0;
*length = psub->length;
frame_code = psub->phase;
if(frame != NULL) *frame = frame_code % 100;
if(gencode != NULL) *gencode = frame_code/100;
return retval;
}


int strcmptrail(char *s1, int l1, char *s2, int l2)
/*
compare strings s1 and s2 of length l1 and l2 as done by strcmp
but ignores all trailing spaces
*/
{
char *fin;
int l, flag=1;

if(l1 > 0) {
	if( (fin = memchr(s1, 0, l1) ) != NULL) l1 = fin - s1;
	}
if(l2 > 0) {
	if( (fin = memchr(s2, 0, l2) ) != NULL) l2 = fin - s2;
	}

if(l2 > l1) {
	flag = -1;
	fin=s1; s1=s2; s2=fin;
	l=l1; l1=l2; l2=l;
	}
l = l2;
fin = s2 + l;
while(s2 < fin) {
	if( *s1 != *s2 ) return (*s1 - *s2)*flag;
	s1++; s2++;
	}
fin= s1+l1-l2;
while(s1 < fin)	{
	if( *s1 != ' ') return flag;
	s1++;
	}
return 0;
}


char *short_descr(int seqnum, char *text, int maxlen)
/*
to get a description of a sequence or of a subsequence
seqnum	the sequence number
text	the string to be loaded with description
maxlen	the max # of chars allowed in text (\0 is put but not counted in maxlen)
return value	a pointer to text 
*/
{
long pinf;
int div, parent, l, deb;
char *p;
text[maxlen]=0;
seq_to_annots(seqnum, &pinf, &div);
memcpy(text,psub->name,16);
p=text+15;
while(*p == ' ') p--;
*(p+1)=0;
l=strlen(text);
parent = (psub->pext <= 0);
if(!parent) { /* subsequence */
	if( read_annots(pinf, div) == NULL) return text;
	pinfo->line[20]=0;
	strcat(text,pinfo->line+4);
	l=strlen(text);
	while(text[l-1]==' ') l--;
	text[l]=0;
	if( ( p = strchr(pinfo->line + 21, '/') ) != NULL) {
		strncat(text, p, maxlen - l);
		l = strlen(text);
		if(l > 75) return text;
		}
	do	{
		next_annots(NULL);
		if( strcmptrail(pinfo->line,20,NULL,0) && 
			strcmptrail(pinfo->line,10,"FT",2) ) return text;
		}
	while(pinfo->line[21]!='/');
	do	{
		strncat(text,pinfo->line+20,maxlen-l);
		l=strlen(text);
		if(l>75) return text;
		next_annots(NULL);
		}
	while ( !strcmptrail(pinfo->line,20,NULL,0) || 
				!strcmptrail(pinfo->line,10,"FT",2) );
	}
else	{ /* parent sequence */
	if( read_annots(pinf, div) == NULL) return text;
	next_annots(NULL);
	if(nbrf) {
		deb=17;
		}
	else	{
		deb=13;
		if(embl || swissprot) {
			while (strncmp(pinfo->line,"DE",2)) {
				next_annots(NULL);
				}
			deb=6;
			}
		}
	do	{
		strncat(text,pinfo->line+deb-2,maxlen-l);
		l=strlen(text);
		if(l>=77) return text;
		next_annots(NULL);
		}
	while( !strncmp(pinfo->line,"  ",2) || !strncmp(pinfo->line,"DE",2) );
	}
return text;
}
		


char *short_descr_p(int seqnum, char *text, int maxlen)
/*
to get a description of a sequence derived from its parent's annotations
if it is a subsequence
seqnum	the sequence number
text	the string to be loaded with description
maxlen	the max # of chars allowed in text (\0 is put but not counted in maxlen)
return value	a pointer to text 
*/
{
char *p;
readsub(seqnum);
if(psub->pext > 0) {
	p=strchr(psub->name,'.');
	*p=0;
	seqnum=isenum(psub->name);
	}
return short_descr(seqnum, text, maxlen);
}

void simpleopen(void) 
{
int recn, tot, sorted;
char mode[5]="r";

#ifdef __INTEL__
strcat(mode, "b");
#endif
#ifdef vms
	unixos=0;
#else
	unixos=1;
#endif
ksub=myopen("SUBSEQ",mode,lrsub);  psub = (struct rsub *)malloc(lrsub);
nseq=read_first_rec(ksub,NULL);
kloc=myopen("LOCUS",mode,lrloc);  ploc = (struct rloc *)malloc(lrloc);
ksmj=myopen("SMJYT",mode,lrsmj);  psmj = (struct rsmj *)malloc(lrsmj);
tot=read_first_rec(ksmj,&recn);  sorted = (recn==tot);
nbrf=swissprot=FALSE;
for(recn=2;recn<=tot;recn++) {
	readsmj(recn);
	if(!strncmp(psmj->name,"01PROTEIN ",10)) { nbrf = TRUE; break; }
	if(!strncmp(psmj->name,"01PRT ", 6)) { swissprot = TRUE; break; }
	if(sorted && strncmp(psmj->name,"01",2)>0 ) break;
	}
kext = NULL;
if(!(nbrf || swissprot)) {
	kext=myopen("EXTRACT",mode,lrext);  pext = (struct rext *)malloc(lrext);
	}
kshrt=myopen("SHORTL",mode,lrshrt);  pshrt = (struct rshrt *)malloc(lrshrt);
readshrt(2);
if(pshrt->val < 0) { /* new hashing format */
	hoffst=2;
	hsub=abs(pshrt->val);  hkwsp=abs(pshrt->next);
	}
else	{ /* old hashing format */
	hoffst=1;
	hsub=9973; hkwsp=1999;
	}
hsub=abs(pshrt->val);  hkwsp=abs(pshrt->next);
pinfo=(struct rinfo *)malloc(lrinfo);

ktxt=myopen("TEXT",mode,lrtxt); 
gcgini();
dir_close(ksmj);
dir_close(ktxt);
ksmj=NULL;
maxa=1;
longa=(maxa-1)/lmot+1;
lenbit=(nseq>maxa? nseq : maxa);
lenw=(lenbit-1)/lmot+1;
klng=myopen("LONGL",mode,lrlng);  plng = (struct rlng*)malloc(lrlng);
return;
} /* end of simpleopen */



#define TOTCODES 17  /* nbre total de codes definis, 0 inclus */
int totcodes=TOTCODES;

char aminoacids[]="RLSTPAGVKNQHEDYCFIMW*X";

struct genetic_code_libel { /* definition d'un code genetique */
	char libel[61]; /* nom du code decrivant ses variants % code standard */
	int code[65]; /* tableau codon->acide amine */
	int ncbi_gc; /* numero NCBI du meme code */
	int codon_init[64]; /* tableau codon initiateur -> acide amine */
	};

/* 
les codons sont numerotes de 1 a 64 selon ordre alphabetique;
le numero 65 est attribue a tout codon avec base hors AcCcGgTtUu
les acides amines sont numerotes selon l'ordre de la variable aminoacids
de un a 20 + * pour stop et X pour inconnu
*/

/* initialisation de tous les codes genetiques */
struct genetic_code_libel genetic_code[TOTCODES] = 
{

{ /* 0: universel */
	{"Universal genetic code"},
/*ANN*/	{9,10,9,10,4,4,4,4,1,3,1,3,18,18,19,18,
/*CNN*/	11,12,11,12,5,5,5,5,1,1,1,1,2,2,2,2,
/*GNN*/	13,14,13,14,6,6,6,6,7,7,7,7,8,8,8,8,
/*TNN*/	21,15,21,15,3,3,3,3,21,16,20,16,2,17,2,17,22},
/*ncbi*/1,
/*init*/{0,0,0,0,0,0,0,0,0,0,0,0,0,0,19,0, /* AUG */
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,19,0, /* CUG */
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,19,0} /* UUG */
}
,
{ /* 1: yeast mt */
	{"CUN=T  AUA=M  UGA=W"},
	{9,10,9,10,4,4,4,4,1,3,1,3,19,18,19,18,
	11,12,11,12,5,5,5,5,1,1,1,1,4,4,4,4,
	13,14,13,14,6,6,6,6,7,7,7,7,8,8,8,8,
	21,15,21,15,3,3,3,3,20,16,20,16,2,17,2,17,22},
	3,
	{0,0,0,0,0,0,0,0,0,0,0,0,19,0,19,0, /* AUA, AUG */
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
}
,
{ /* 2: :    MITOCHONDRIAL CODE OF VERTEBRATES */
	{"AGR=*  AUA=M  UGA=W"},
     {9,10,9,10,4,4,4,4,21,3,21,3,19,18,19,18,11,12,11,12,
     5,5,5,5,1,1,1,1,2,2,2,2,13,14,13,14,6,6,6,6,7,7,7,7,8,8,8,8,21,15,
     21,15,3,3,3,3,20,16,20,16,2,17,2,17,22},
	2,
	{0,0,0,0,0,0,0,0,0,0,0,0,19,19,19,19, /* AUN */
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,19,0, /* GUG */
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
}
,
{ /* 3:   MITOCHONDRIAL CODE OF FILAMENTOUS FUNGI */
	{"UGA=W"},
     {9,10,9,10,4,4,4,4,1,3,1,3,18,18,19,18,11,12,11,12,5,5,5,
     5,1,1,1,1,2,2,2,2,13,14,13,14,6,6,6,6,7,7,7,7,8,8,8,8,21,15,21,
     15,3,3,3,3,20,16,20,16,2,17,2,17,22},
	4,
	{0,0,0,0,0,0,0,0,0,0,0,0,19,19,19,19, /* AUN */
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,19,0, /* CUG */
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,19,0, /* GUG */
	0,0,0,0,0,0,0,0,0,0,0,0,19,0,19,0} /* UUR */
}
,
{ /* 4:    MITOCHONDRIAL CODE OF INSECT AND PLATYHELMINTHES  */
	{"AUA=M  UGA=W  AGR=S"},
     {9,10,9,10,4,4,4,4,3,3,3,3,19,18,19,18,11,12,11,12,5,5,5,
     5,1,1,1,1,2,2,2,2,13,14,13,14,6,6,6,6,7,7,7,7,8,8,8,8,21,15,21,
     15,3,3,3,3,20,16,20,16,2,17,2,17,22},
	5,
	{0,0,0,0,0,0,0,0,0,0,0,0,19,19,19,19, /* AUN */
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,19,0, /* GUG */
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,19,0} /* UUG */
}
,
{ /* 5:    Nuclear code of Candida cylindracea (see nature 341:164) */
	{"CUG=S"},
     	{9,10,9,10,4,4,4,4,1,3,1,3,18,18,19,18,
	11,12,11,12,5,5,5,5,1,1,1,1,2,2,3,2,
	13,14,13,14,6,6,6,6,7,7,7,7,8,8,8,8,
	21,15,21,15,3,3,3,3,21,16,20,16,2,17,2,17,22},
	12,
	{0,0,0,0,0,0,0,0,0,0,0,0,0,0,19,0, /* AUG */
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,19,0, /* CUG */
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
}
,
{ /* 6:   NUCLEAR CODE OF CILIATA: UAR = Gln = Q */
	{"UAR=Q"},
     {9,10,9,10,4,4,4,4,1,3,1,3,18,18,19,18,11,12,11,12,5,5,5,
     5,1,1,1,1,2,2,2,2,13,14,13,14,6,6,6,6,7,7,7,7,8,8,8,8,11,15,11,
     15,3,3,3,3,21,16,20,16,2,17,2,17,22},
	6,
	{0,0,0,0,0,0,0,0,0,0,0,0,0,0,19,0, /* AUG */
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
}
,
{ /* 7:   NUCLEAR CODE OF EUPLOTES */
	{"UGA=C"},
     {9,10,9,10,4,4,4,4,1,3,1,3,18,18,19,18,11,12,11,12,5,5,5,
     5,1,1,1,1,2,2,2,2,13,14,13,14,6,6,6,6,7,7,7,7,8,8,8,8,21,15,21,
     15,3,3,3,3,16,16,20,16,2,17,2,17,22},
	10,
	{0,0,0,0,0,0,0,0,0,0,0,0,0,0,19,0, /* AUG */
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
}
,
{ /* 8:   MITOCHONDRIAL CODE OF ECHINODERMS */
	{"UGA=W  AGR=S  AAA=N"},
     	{10,10,9,10,4,4,4,4,3,3,3,3,18,18,19,18,
	11,12,11,12,5,5,5,5,1,1,1,1,2,2,2,2,
	13,14,13,14,6,6,6,6,7,7,7,7,8,8,8,8,
	21,15,21,15,3,3,3,3,20,16,20,16,2,17,2,17,22},
	9,
	{0,0,0,0,0,0,0,0,0,0,0,0,0,0,19,0, /* AUG */
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
}
,
{ /* 9:   MITOCHONDRIAL CODE OF ASCIDIACEA */
	{"UGA=W  AGR=G  AUA=M"},
     	{9,10,9,10,4,4,4,4,7,3,7,3,19,18,19,18,
	11,12,11,12,5,5,5,5,1,1,1,1,2,2,2,2,
	13,14,13,14,6,6,6,6,7,7,7,7,8,8,8,8,
	21,15,21,15,3,3,3,3,20,16,20,16,2,17,2,17,22},
	13,
	{0,0,0,0,0,0,0,0,0,0,0,0,0,0,19,0, /* AUG */
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
}
,
{ /* 10:   MITOCHONDRIAL CODE OF PLATYHELMINTHES */
	{"UGA=W  AGR=S  UAA=Y AAA=N"},
	{10,10,9,10,4,4,4,4,3,3,3,3,18,18,19,18,
	11,12,11,12,5,5,5,5,1,1,1,1,2,2,2,2,
	13,14,13,14,6,6,6,6,7,7,7,7,8,8,8,8,
	15,15,21,15,3,3,3,3,20,16,20,16,2,17,2,17,22},
	14,
	{0,0,0,0,0,0,0,0,0,0,0,0,0,0,19,0, /* AUG */
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
}
,
{ /* 11:   NUCLEAR CODE OF BLEPHARISMA */
	{"UAG=Q"},
/*ANN*/	{9,10,9,10,4,4,4,4,1,3,1,3,18,18,19,18,
/*CNN*/	11,12,11,12,5,5,5,5,1,1,1,1,2,2,2,2,
/*GNN*/	13,14,13,14,6,6,6,6,7,7,7,7,8,8,8,8,
/*TNN*/	21,15,11,15,3,3,3,3,21,16,20,16,2,17,2,17,22},
	15,
	{0,0,0,0,0,0,0,0,0,0,0,0,0,0,19,0, /* AUG */
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
}
,
{ /* 12:   NUCLEAR CODE OF BACTERIA: differs only for initiation codons */
	{"NUG=AUN=M when initiation codon"},
/*ANN*/	{9,10,9,10,4,4,4,4,1,3,1,3,18,18,19,18,
/*CNN*/	11,12,11,12,5,5,5,5,1,1,1,1,2,2,2,2,
/*GNN*/	13,14,13,14,6,6,6,6,7,7,7,7,8,8,8,8,
/*TNN*/	21,15,21,15,3,3,3,3,21,16,20,16,2,17,2,17,22},
	11,
	{0,0,0,0,0,0,0,0,0,0,0,0,19,19,19,19, /* AUN */
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,19,0, /* CUG */
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,19,0, /* GUG */
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,19,0} /* UUG */
}
,
{ /* 13: Chlorophycean Mitochondrial */
	{"UAG=Leu"},
/*ANN*/	{9,10,9,10,4,4,4,4,1,3,1,3,18,18,19,18,
/*CNN*/	11,12,11,12,5,5,5,5,1,1,1,1,2,2,2,2,
/*GNN*/	13,14,13,14,6,6,6,6,7,7,7,7,8,8,8,8,
/*TNN*/	21,15,2,15,3,3,3,3,21,16,20,16,2,17,2,17,22},
/*ncbi*/16,
/*init*/{0,0,0,0,0,0,0,0,0,0,0,0,0,0,19,0, /* AUG */
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
}
,
{ /* 14:    MITOCHONDRIAL CODE OF TREMATODE  */
	{"AUA=M  UGA=W  AGR=S AAA=N"},
     {10,10,9,10,4,4,4,4,3,3,3,3,19,18,19,18,11,12,11,12,5,5,5,
     5,1,1,1,1,2,2,2,2,13,14,13,14,6,6,6,6,7,7,7,7,8,8,8,8,21,15,21,
     15,3,3,3,3,20,16,20,16,2,17,2,17,22},
	21,
	{0,0,0,0,0,0,0,0,0,0,0,0,0,0,19,0, /* AUG */
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,19,0, /* GUG */
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0} 
}
,
{ /* 15: TAG-Leu,TCA-stop */
	{"UAG=L UCA=*"},
/*ANN*/	{9,10,9,10,4,4,4,4,1,3,1,3,18,18,19,18,
/*CNN*/	11,12,11,12,5,5,5,5,1,1,1,1,2,2,2,2,
/*GNN*/	13,14,13,14,6,6,6,6,7,7,7,7,8,8,8,8,
/*TNN*/	21,15,2,15,21,3,3,3,21,16,20,16,2,17,2,17,22},
/*ncbi*/22,
/*init*/{0,0,0,0,0,0,0,0,0,0,0,0,0,0,19,0, /* AUG */
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0} 
}
,
{ /* 16: Thraustochytrium-mt */
	{"UUA=*"},
/*ANN*/	{9,10,9,10,4,4,4,4,1,3,1,3,18,18,19,18,
/*CNN*/	11,12,11,12,5,5,5,5,1,1,1,1,2,2,2,2,
/*GNN*/	13,14,13,14,6,6,6,6,7,7,7,7,8,8,8,8,
/*TNN*/	21,15,21,15,3,3,3,3,21,16,20,16,21,17,2,17,22},
/*ncbi*/23,
/*init*/{0,0,0,0,0,0,0,0,0,0,0,0,0,0,19,19, /* AUG AUU */
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,19,0, /* GUG */
	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0} 
}


/*       1         2
1234567890123456789012
RLSTPAGVKNQHEDYCFIMW*X
*/

};


char *get_code_descr(int code)
/* 
get a 60-letter (or less) description of a variant genetic code
return value	pointer to the description, not to be altered!
*/
{
if(code>=0 && code<totcodes)
	return genetic_code[code].libel ;
else	
	return "Unknown genetic code. Standard code is used.";
}


int calc_codon_number(char *codon)
{
static char nucleotides[] = "AaCcGgTtUu";
static int nucnum[5] = {0,1,2,3,3};
int num, i, base;
char *p;

num = 0;
for(i = 1; i <= 3; i++) {
	p = strchr(nucleotides, *codon);
	if(p == NULL) {
		num = 64;
		break;
		}
	else
		base = (p-nucleotides)/2;
	num = num * 4 + nucnum[base];
	codon++;
	}
return num;
}


char codaa(char *codon, int code)
/*
amino acid translation:
codon	a 3-base string
code	the genetic code to be used
return value	the amino acid as 1 character
*/
{
struct genetic_code_libel *pdata;
int num;

num = calc_codon_number(codon);
if(code < 0 || code >= totcodes)code = 0;/*use regular code if unknown number */
pdata = &genetic_code[code]; /* ici ecriture plus compacte mal compilee sur PC*/
return aminoacids[ pdata->code[num] - 1 ];
}


int get_ncbi_gc_number(int gc)
{ /* from acnuc to ncbi genetic code number */
return genetic_code[gc].ncbi_gc;
}


int get_acnuc_gc_number(int ncbi_gc)
{ /* from ncbi to acnuc genetic code number (returns 0 if not found) */
int num;

for( num = 0; num < totcodes; num++ ) 
	if(genetic_code[num].ncbi_gc == ncbi_gc) return num;
return 0;
}


char translate_init_codon(int numseq, int gc, int debut_codon /* 1, 2, or 3 */)
{
char codon[4];
int point, special_init = TRUE;
static int num_5_partial = 0;

if(num_5_partial == 0) num_5_partial = iknum("5'-PARTIAL", kkey);
if(debut_codon != 1) special_init = FALSE;
else	{ /* la seq est-elle 5'-PARTIAL ? */
	readsub(numseq);
	point = psub->plkey;
	while(point != 0) {
		readshrt(point);
		if(pshrt->val == num_5_partial) {
			special_init = FALSE;
			break;
			}
		point = pshrt->next;
		}
	}
gfrag(numseq, debut_codon, 3, codon);
if(special_init)  /* traduction speciale du codon initiateur */
	return init_codon_to_aa(codon, gc);
else	return codaa(codon, gc);
}


char init_codon_to_aa(char *codon, int gc)
{
int num, aa;
struct genetic_code_libel *pdata;

num = calc_codon_number(codon);
if(num >= 64) return 'X';
/* use regular code if unknown number */
if(gc < 0 || gc >= totcodes) gc = 0; 
pdata = &genetic_code[gc];
aa = pdata->codon_init[num];
/* if not listed in expected init codons */
if(aa == 0) aa = pdata->code[num];
return aminoacids[aa - 1];
}


void compact(char *chaine)
{
int l;
char *p, *q;

l=strlen(chaine); p=chaine+l;
while( *(--p) == ' ' && p>=chaine) *p=0;
while((p=strchr(chaine,' '))!=NULL) {
	q=p+1;
	while(*q==' ') q++;
	l=q-p;
	while(*q!=0) {*(q-l) = *q; q++; }
	*(q-l)=0;
	}
}


int trim_key(char *name) /* remove trailing spaces */
{
char *p;
int l = strlen(name);
p = name + l - 1;
while( p >= name && *p == ' ' ) *(p--) = 0;
return (p + 1) - name;
}


void seq_to_annots(int numseq, long *faddr, int *pdiv)
{
unsigned pinf;
int int_faddr, mere;
char *p;
static char nom_mere[L_MNEMO + 1];

readsub(numseq);
if(psub->pext > 0) { /* seq fille */
	readshrt(psub->plinf);
	pinf = pshrt->val;
	if(big_annots) {
		p = strchr(psub->name, '.');
		*p = 0;
		strcpy(nom_mere, psub->name);
		mere = isenum(nom_mere);
		readsub(mere);
		readloc(psub->plinf);
		*pdiv = ploc->div;
		readsub(numseq); /* important pour programmeur */
		}
	}
else 	{ /* seq mere */
	readloc(psub->plinf);
	pinf = ploc->pinf;
	if(big_annots) *pdiv = ploc->div;
	}
if(big_annots) *faddr = (long) pinf;
else	{
	goffset(pinf, pdiv, &int_faddr);
	*faddr = (long) int_faddr;
	}
}

static int current_div_annots = -1;

char *read_annots(long faddr, int div)
/* returns NULL if error, or pinfo->line if ok */
{
int i;
off_t fpos;

if( div < 0 || div > divisions ) return NULL;
/* tester si le fichier a ouvrir n'est pas ouvert */
if(!annotopened[div]) gcgors("inf", div, TRUE); 
fpos = (off_t)(unsigned long)faddr;
if( fseeko(divannot[div], fpos, SEEK_SET) != 0 ) return NULL;
current_div_annots = div;
if( fgets(pinfo->line, sizeof(pinfo->line), divannot[div]) == NULL ) return NULL;
i = strlen(pinfo->line) - 1;
while( i >= 0 && isspace( pinfo->line[i] ) ) i--;
pinfo->line[i + 1] = 0;
return pinfo->line;
}


char *next_annots(long *pfaddr)
/* lit la ligne d'annotations qui suit la derniere lue
si pfaddr != NULL le charge avec faddr du debut de la ligne lue
rend la ligne lue, ou NULL si erreur 
*/
{
int l;
off_t fpos;

if( current_div_annots < 0 || current_div_annots > divisions || 
	!annotopened[current_div_annots] ) return NULL;
if(pfaddr != NULL) {
	fpos = ftello( divannot[current_div_annots] );
	*(unsigned long *)pfaddr = (unsigned long)fpos;
	}
if( fgets(pinfo->line, sizeof(pinfo->line), 
	divannot[current_div_annots]) == NULL ) return NULL;
l = strlen(pinfo->line) - 1;
while( l >= 0 && isspace( pinfo->line[l] ) ) l--;
pinfo->line[l + 1] = 0;
return pinfo->line;
}


void majuscules(char *name)
{
name--;
while(*(++name) != 0) *name = toupper(*name);
}


char *check_acnuc_gcgacnuc(void)
{
#ifdef unix
if(getenv("acnuc") == NULL || getenv("gcgacnuc") == NULL)
	return "Environment variables acnuc and gcgacnuc should be defined";
#elif defined(vms)
if(getenv("BANK") == NULL || getenv("GCGBANK") == NULL)
	return "Logical variables BANK and GCGBANK should be defined";
#else  /* Mac ou PC */
if(getenv("acnuc") == NULL || getenv("gcgacnuc") == NULL)
#ifdef __INTEL__ /* PC */
	return "File acnuc.ini locating bank is not found";
#else  /* Mac */
	return "File acnuc.preference locating bank is not found";
#endif
#endif
return NULL;
}


char *translate_cds(int seqnum)
/* traduction d'un cds avec codon initiateur traite et * internes ==> X
rendue dans memoire allouee ici qu'il ne faut pas modifier
retour NULL si pb lecture de la seq
*/
{
static char *buffer = NULL;
static int lbuffer = 0;
int debut_codon, longueur, pos, code;
char codon[4], *p;

readsub(seqnum);
debut_codon = (psub->phase % 100) + 1;
code = psub->phase / 100;
longueur = (psub->length - debut_codon + 1)/3;
if(longueur > lbuffer) {
	if(buffer != NULL) free(buffer);
	buffer = (char *)malloc(longueur + 1);
	lbuffer = longueur;
	}
if(buffer == NULL) {lbuffer = 0; return NULL; }
buffer[0] = translate_init_codon(seqnum, code, debut_codon);
debut_codon += 3;
for(pos = 1; pos < longueur; pos++) {
	if( gfrag(seqnum, debut_codon, 3, codon) == 0) return NULL;
	buffer[pos] = codaa(codon,code);
	debut_codon += 3;
	}
buffer[longueur] = 0;
while( (p = strchr(buffer, '*') ) != NULL && p - buffer < longueur - 1 )
	*p = 'X';
return buffer;
}


/*

New LOCUS Format:

---------+---------+---------+---------+---------+---------+---------+---------
1       10        20        30        40        50        60        70       79
LOCUS       16Char_LocusName 99999999999 bp ss-snoRNA  circular DIV DD-MMM-YYYY

Positions  Contents
---------  --------
01-05      LOCUS
06-12      spaces
13-28      Locus name
31-31      space
30-40      Length of sequence, right-justified
41-41      space
42-43      bp
44-44      space
45-47      spaces, ss- (single-stranded), ds- (double-stranded), or
           ms- (mixed-stranded)
48-53      NA, DNA, RNA, tRNA (transfer RNA), rRNA (ribosomal RNA), 
           mRNA (messenger RNA), uRNA (small nuclear RNA), snRNA,
           snoRNA. Left justified.
54-55      space
56-63      'linear' followed by two spaces, or 'circular'
64-64      space
65-67      The division code (see Section 3.3)
68-68      space
69-79      Date, in the form dd-MMM-yyyy (e.g., 15-MAR-1991)



Old LOCUS format :
1       10        20        30        40        50        60        70       79
LOCUS       ABCRRAA       118 bp ss-rRNA            RNA       15-SEP-1990
Positions   	Contents

1-12	LOCUS
13-22	Locus name
23-29	Length of sequence, right-justified
31-32	bp
34-36	Blank, ss- (single-stranded), ds- (double-stranded), or
	 ms- (mixed-stranded)
37-41	Blank, DNA, RNA, tRNA (transfer RNA), rRNA (ribosomal RNA), 
	mRNA (messenger RNA), uRNA (small nuclear RNA), snRNA, or scRNA
43-52	Blank (implies linear) or circular
53-55	The division code (see Section 3.3)
63-73	Date, in the form dd-MMM-yyyy (e.g., 15-MAR-1991)

*/

int decode_locus_rec(char *line, char **pname, char **pmolec, 
	int *circular, char **pdivision, char **pdate)
{
static char name[20], molec[10], division[5], date[12];
int newformat, length;
char *p, *q;

newformat = (line[70] == '-');
p = line + 12; 
q = p; while (*q != ' ') q++;
memcpy(name, p, q - p); name[q - p] = 0;
if(pname != NULL) *pname = name;
p = q + 1; while (*p == ' ') p++;
sscanf(p, "%d", &length);
if(newformat) {
	memcpy(molec, line + 47 , 6); 
	molec[6] = 0;
	}
else	{
	memcpy(molec, line +  36, 5); molec[5] = 0;
	}
if(pmolec != NULL) *pmolec = molec;
if(circular != NULL) *circular = 
	strncmp(line + (newformat ? 55 : 42), "circular", 8) == 0;
memcpy(division, line + (newformat ? 64 : 52),3); division[3] = 0;
if(pdivision != NULL) *pdivision = division;
memcpy(date, line + (newformat ? 68 : 62), 11); date[11] = 0;
if(pdate != NULL) *pdate = date;
return length;
}


void quick_list_meres(int *blist)
{
#ifdef unix
char *fname;
FILE *in;
int lu;

fname = prepare_env_var("acnuc");
strcat(fname, "MERES");
in = fopen(fname, "r");
if(in != NULL) {
	lu = fread(blist, sizeof(int), lenw, in);
	fclose(in);
	if(lu == lenw) return;
	}
#endif

lngbit(2, blist);

return;
}