#include "dir_acnuc.h"
#include <time.h>
#include <ctype.h>

void seq_to_annots(int numseq, long *faddr, int *div);
char *read_annots(long faddr, int div);
char *next_annots(long *pfaddr);
char *translate_cds(int seqnum);


#define NBRF_BACK_SLASHES "\\\\\\\n"

typedef enum { simple=1, translate, fragment, feature, region } extract_option;

typedef enum { acnuc=1, gcg, fasta, analseq, flat } out_option;

typedef struct {
	out_option option;
	int previous_len;
	int width;
	union 	{
		FILE *outfile;
		struct analseq_files {
			DIR_FILE *kin; int num_in;
			DIR_FILE *kse; int num_se;
			} *analseq;
		} fichier;
	} out_format_struct;
	
typedef struct {
	extract_option choix; /* donne le genre d'extraction */
	union	{ /* selon le genre d'extraction */
		char feature_name[20]; /* pour options feature + region */
		struct	{ /* pour options fragment + region */
			char *end_points;
			char *min_end_points;
			} extremites;
		} data;
	} struct_extract_data;
			
#define SIMEXT_NFMAX 500
typedef struct {
	int nfmax;  /* maximum # of fragments in a virtual subseq */
	int newsim; /* TRUE when a new virtual subseq has just been created*/
	int pinfme; /* pointer to annots of parent of virtual subseq */
	int siminf; /* pointer to annots of virtual subseq */
	int div;    /* div of these annots */
	int nfrags; /* # of fragments in virtual subseq (from 1)*/
	int valext[SIMEXT_NFMAX][4];
/* valext(.,0)=beginning of fragment */
/* valext(.,1)=end of fragment	*/
/* valext(.,2)=pointer to nucleot for that fragment */
/* valext(.,3)=# of parent seq of that fragment	*/
	} simext_struct;
	
typedef union  {
	struct {
		char name[20];
		int pseq;
		int length;
		int phase;
		} normal_rec;
	struct {
		int total;
		int nbseq;
		int espace;
		} first_rec;
	} kin_union;
#define kin_size sizeof(kin_union)
	
typedef union {
	char seq[90];
	int total;
	} kse_union;
#define kse_size 90
	
	
/* extern variables */
extern simext_struct *psimext; /* pointeur (pas structure) alloue dans gfrag */
	
/* global variables */
static simext_struct simext;  /* structure elle-meme allouee ici */
static out_format_struct out_struct; /* infos sur le format de sortie */
static struct_extract_data extract_data; /* infos sur le type d'extraction */
static struct analseq_files analseq_data;
static int gcg_checksum;
static int cds_type=0;
static kin_union kin_record;
static kse_union kse_record;
/* var concernant les bornes a extraire et les bornes minimales */
static int frag_begin_to_end; /* TRUE ssi bornes vont de debut a fin */
static int dlu, flu, xd, xf, dlum, flum, xdm, xfm;


/* prototypes */
static int read_location(long pannot, int div, char *location, int maxloclen);
int find_label(char *label, char *access, long *pannot, int *pdiv, int *isub);
static int write_one_seq_line(out_format_struct *out_format, char *seq, 
	int length);
static int comp_gcg_checksum(char *seq, int length, int previous_len);
static int write_seq_header(char *mnemo, int iseq, int debut, int fin,
	out_format_struct *out_format, char *descri, long deb_info_fille, 
	int div_fille);
int prep_extract(out_option option, char *fname, extract_option choix,
	 char *feature_name, char *bornes, char *min_bornes, char **message);
int extract_1_seq(int numseq, char *bornes);
static int write_seq_footer(out_format_struct *out_format, char *name, 
	char Nuc_ou_Prot);
static int extraire_bornes(char *bornes, int *debut, int *fin, 
	int *xdebut , int *xfin);
void fin_extract(void);
static int end_gcg(out_format_struct *out_format, char *mnemo, 
	char Nuc_ou_Prot);
static int proc_locat_fragment(char *location, char *debut, char *fin, 
	int tot_locat,
	int max_locat, long pinf, int div, int mere, int *part5, int *part3);
static int proc_location(long *pinf, int div, char *nom_feature, int mere, 
	int debut_info_mere);
static int proc_location_sw(long *pinf, int div, char *nom_feature, int mere, 
	int debut_info_mere);
static int mod_limites(void);
static int write_fragment(out_format_struct *out_format, char *mnemo, int iseq, 
	int deb, int fin, int debn, int finn, char *description, 
	long deb_info_fi, int div);

/* external prototypes */
char translate_init_codon(int numseq, int gc, int debut_codon);

static char *nextpar(char *debut)
{
debut++;
do	{
	if(*debut=='(') {
		debut=nextpar(debut);
		if(debut==NULL) return NULL;
		debut++;
		}
	if(*debut==')') {
		return debut;
		}
	else debut++;
	}
while(*debut!=0);
return NULL;
}


static int read_location(long pannot, int div, char *location, int maxloclen)
/* charge location avec la location lue a partir d'adresse pannot
rend la longueur chargee, ou 0 si location + longue que maxloclen */
{
int i, lloc=0;
char *ps;

if( read_annots(pannot, div) == NULL) return 0;
while((ps=strchr(pinfo->line,'/'))==NULL) {
	i=strlen(pinfo->line)-21;
	if(lloc+i>=maxloclen) return 0;
	if(i>0) {
		memcpy(location+lloc,pinfo->line+21,i);
		lloc+=i;
		}
	next_annots(NULL);
	if(strcmptrail(pinfo->line,20,NULL,0)!=0 && 
		strcmptrail(pinfo->line,20,"FT",2) !=0 ) break;
	}
if(ps!=NULL) {
	i=ps-pinfo->line+1;
	if(strcmptrail(pinfo->line+2,i-3,NULL,0)!=0 && i>22) {
		if(lloc+i-22>=maxloclen) return 0;
		memcpy(location+lloc,pinfo->line+21,i-22);
		lloc+=(i-22);
		}
	}
location[lloc]=0;
compact(location);
ps=location; while(*ps!=0) { *ps=toupper(*ps); ps++; }
return strlen(location);
}


int find_label(char *label, char *access, long *pannot, int *pdiv, int *isub)
/* recherche de la feature  access-label
si access==NULL, recherche locale ds la table de seq # *isub,
si trouve, *isub = # mere ou trouve et *pannot = pointeur info debut features
        et rendu 0
si pas trouve, rendu 1 si accession manque  et 2 si label manque 
*/
{
char *p;
int nacc, trouve;
long point;

if(access==NULL)
	nacc=0;
else	{
	nacc=fcode(kacc,access,ACC_LENGTH);
	if(nacc!=0)
		nacc=pacc->plsub;
	if(nacc==0) return 1;
	}
do	{
	if(nacc!=0) 	{
		readshrt(nacc);
		*isub=pshrt->val; nacc=pshrt->next;
		}
	seq_to_annots(*isub, pannot, pdiv);
	if( read_annots(*pannot, *pdiv) == NULL) return 2;
	trouve=FALSE;
	do	{
		next_annots(pannot);
		if(strncmp(pinfo->line,"BASE COUNT",10)==0 || 
			strncmp(pinfo->line,"SQ ",3)==0) break;
		trouve = !strncmp(pinfo->line,"FEATURES",8) || 
				!strncmp(pinfo->line,"FH ",3);
		}
	while(!trouve);
	if(!trouve) break;
	point= *pannot;
	do	next_annots(&point);
	while( embl && strncmp(pinfo->line,"FH",2) == 0);
	while( (genbank && strcmptrail(pinfo->line, 2, NULL, 0) == 0 ) || 
				(embl && strncmp(pinfo->line, "FT", 2) == 0) ) {
		*pannot = point;
		do	{
			p = pinfo->line - 1; while(*(++p) != 0) *p = toupper(*p);
			p=strstr(pinfo->line,"/LABEL=");
			if(p != NULL) {
				if(*(p+7)==0) {
					next_annots(&point);
					p = pinfo->line - 1; while(*(++p) != 0) *p = toupper(*p);
					p=pinfo->line+14;
					}
				if(strncmp(p+7,label,strlen(label)) == 0) return 0;
				}
			next_annots(&point);
			}
		while ( (genbank && strcmptrail(pinfo->line,10,NULL,0) == 0 ) ||
				(embl && strcmptrail(pinfo->line,10,"FT",2) == 0) );
		}
	}
while(nacc!=0);
return 2;
}


static int write_one_seq_line(out_format_struct *out_format, char *seq, int length)
/* pour ecrire une ligne de seq de length residus 
selon la structure de fichier define par out_format
retour: !=0 ssi erreur d'ecriture
*/
{
int retval;
char *p;

seq[length] = 0;
if( ( !(nbrf || swissprot) ) && out_format->option!=flat) {
	p=seq; while(*p!=0) {*p=toupper(*p); p++; }
	}
if(out_format->option==analseq) {
	if( !(nbrf || swissprot) && extract_data.choix != translate ) {
		p=seq; while(*p) {if(*p=='T') *p='U'; p++; }
		}
	memset(kse_record.seq,' ',kse_size);
	memcpy(kse_record.seq,seq,length);
	retval=dir_write(out_format->fichier.analseq->kse,
		++(out_format->fichier.analseq->num_se), 1, &kse_record);
	}
else if(out_format->option==gcg) {
	int k;
	fprintf(out_format->fichier.outfile,"%8d ",out_format->previous_len+1);
	for(k=0; k<length; k+=10)
		fprintf(out_format->fichier.outfile," %.10s",seq+k);
	fprintf(out_format->fichier.outfile,"\n");
	if(out_format->previous_len==0) gcg_checksum=0;
	gcg_checksum += comp_gcg_checksum(seq,length,out_format->previous_len);
	}
else if(out_format->option==acnuc || out_format->option==fasta) {
	fprintf(out_format->fichier.outfile,"%s\n",seq);
	}
else if(out_format->option==flat) {
	int k;
	if(nbrf) {
		fprintf(out_format->fichier.outfile,"%8d",
			out_format->previous_len+1);
		for(k=0; k<length && k<30; k++)
			fprintf(out_format->fichier.outfile," %c",seq[k]);
		fprintf(out_format->fichier.outfile,"\n");
		if(length>=31) {
			fprintf(out_format->fichier.outfile,"%8d",
				out_format->previous_len+31);
			for(k=30; k<length; k++)
			     fprintf(out_format->fichier.outfile," %c",seq[k]);
			fprintf(out_format->fichier.outfile,"\n");
			}
		}
	else if(embl || swissprot) {
		char outline[81];
		memset(outline,' ',80);
		p=outline+4;
		for(k=0; k<length; k+=10) {
			sprintf(p," %.10s",seq+k);
			p+=11;
			}
		if(embl) {
			outline[strlen(outline)]=' ';
			p=outline+70;
			sprintf(p,"%10d",out_format->previous_len+length);
			}
		fprintf(out_format->fichier.outfile,"%s\n",outline);
		}
	else	{
		fprintf(out_format->fichier.outfile, "%9d",
					out_format->previous_len+1);
		for(k=0; k<length; k+=10)
			fprintf(out_format->fichier.outfile," %.10s",seq+k);
		fprintf(out_format->fichier.outfile,"\n");
		}
	}
else	{
	retval=1;
	}
out_format->previous_len += length;
if( out_format->option != analseq) retval=ferror(out_format->fichier.outfile);
return retval;
}


static int comp_gcg_checksum(char *seq, int length, int previous_len)
{
int retval=0, i;
char *fin;

fin=seq+length;
i=previous_len;
while(seq<fin) {
	retval += ( (int)(*seq) ) * ( i % 57 + 1);
	seq++; i++;
	}
return retval;
}


static int write_seq_header(char *mnemo, int iseq, int debut, int fin,
	out_format_struct *out_format, char *descri, long deb_info_fille, int div_fille)
/* ecrire l'en tete d'une sequence ou fragment a extraire ds fichier et format 
decrit ds out_format 
rend 0 si ok, !=0 si erreur ecriture ds fichier
*/
{
int i, pha, sgc, length, erreur, div_mere;
long point, dime;
char mark;

if(iseq== -1) {
	length=0; pha=0;
	for(i=0; i<simext.nfrags; i++)
		length+=(abs(simext.valext[i][1]-simext.valext[i][0])+1);
	dime = simext.pinfme;
	deb_info_fille = simext.siminf;
	div_fille = div_mere = simext.div;
	}
else	{
	readsub(iseq); pha=psub->phase;
	length=(abs(fin-debut)+1);
	}
if(extract_data.choix==simple && psub->pext<=0 && out_format->option==flat) { 
	/* ecriture toutes annots */
	seq_to_annots(iseq, &point, &div_mere);
	if( read_annots(point, div_mere) != NULL) {
		fprintf(out_format->fichier.outfile,"%s\n",pinfo->line);
		do	{
			next_annots(NULL);
			fprintf(out_format->fichier.outfile,"%s\n",pinfo->line);
			}
		while((nbrf && strncmp(pinfo->line,"SUMMARY",7)!=0) || 
			(!nbrf && (strncmp(pinfo->line,"ORIGIN",6)!=0 && 
			strncmp(pinfo->line,"SQ ",3)!=0)));
		}
	if(nbrf) {
		fprintf(out_format->fichier.outfile,"SEQUENCE\n         ");
		for(i=5; i<=30; i+=5) fprintf(out_format->fichier.outfile,"%9d ",i);
		fprintf(out_format->fichier.outfile,"\n");
		}
	goto fin_header;
	}
sgc=pha/100; pha = pha % 100;
if(out_format->option == analseq) goto fin_header;
if(out_format->option == fasta)
	mark='>';
else if(out_format->option == acnuc)
	mark=';';
else
	mark=0;
if(out_format->option==flat) {
/* format plat pour autre que sequence mere */
	if(embl) {
		fprintf(out_format->fichier.outfile,
			"ID   %s  standard; ; ; %d BP.\n",mnemo,length);
		}
	else if(swissprot) {
		fprintf(out_format->fichier.outfile,
			"ID   %s  STANDARD;   PRT;  %d AA.\n",mnemo,length);
		}
	else	{ /* format GenBank */
		fprintf(out_format->fichier.outfile,
			"LOCUS       %-10.10s%7d BP\n",mnemo,length);
		}
	}
if(iseq != -1) {
	readsub(iseq);
	if(psub->pext > 0) {
		seq_to_annots( iseq, &deb_info_fille, &div_fille);
		readsub(iseq);
		readext(psub->pext);
		seq_to_annots(pext->mere, &dime, &div_mere);
		}
	else seq_to_annots(iseq, &dime, &div_mere);
	}
if(out_format->option == gcg) rewind(out_format->fichier.outfile);
else if(out_format->option != flat) 	{
	if(mark!=0) fputc(mark,out_format->fichier.outfile);
	fprintf(out_format->fichier.outfile,
			"%-20.20s %6d residues",mnemo,length);
	if(out_format->option == acnuc)
		fprintf(out_format->fichier.outfile,
				" Frame%2d Code%3d",pha,sgc);
	fprintf(out_format->fichier.outfile,"\n");
	}
if(out_format->option == fasta ) goto fin_header;
if( read_annots(dime, div_mere) == NULL) goto fin_header;
next_annots(NULL);
if( embl || swissprot ) {
/* format EMBL */
	while(strncmp(pinfo->line,"RN",2) !=0 && 
			strncmp(pinfo->line,"SQ",2) !=0) {
		if(strncmp(pinfo->line,"DE",2)==0 ||
				strncmp(pinfo->line,"AC",2)==0) {
			if(mark!=0) fputc(mark,out_format->fichier.outfile);
			fprintf(out_format->fichier.outfile,"%s\n",pinfo->line);
			}
		next_annots(NULL);
		}
	}
else if(nbrf) { /* format NBRF */
	while(strcmptrail(pinfo->line,10,NULL,0)==0 || 
		strncmp(pinfo->line,"TITLE",5)==0) {
		if(mark!=0) fputc(mark,out_format->fichier.outfile);
		fprintf(out_format->fichier.outfile,"%s\n",pinfo->line);
		next_annots(NULL);
		}
	}
else	{ /* format GenBank */
	while(strcmptrail(pinfo->line,10,NULL,0)==0 || 
		strncmp(pinfo->line,"DEFINITION",10)==0 ||
		strncmp(pinfo->line,"ACCESSION",9)==0) {
		if(mark!=0) fputc(mark,out_format->fichier.outfile);
		fprintf(out_format->fichier.outfile,"%s\n",pinfo->line);
		next_annots(NULL);
		}
	}
if(deb_info_fille != 0) { /* write subseq annotations skipping /translation= */
	int in_transl=FALSE;
	if(out_format->option==flat && !nbrf) {
		if(embl) {
			fprintf(out_format->fichier.outfile,
			"FH   Key             Location/Qualifiers\nFH\n");
			}
		else if(!swissprot)
			fprintf(out_format->fichier.outfile,
			"FEATURES             Location/Qualifiers\n");
		}
	if( read_annots(deb_info_fille, div_fille) != NULL) {
		do	{
			if(!in_transl) {
				if(mark!=0) 
				   fputc(mark, out_format->fichier.outfile);
				fprintf(out_format->fichier.outfile, "%s\n", 
					pinfo->line);
				}
			next_annots(NULL);
			if(strstr(pinfo->line,"/translation=")!=NULL)
				in_transl=TRUE;
			else if(in_transl && strchr(pinfo->line,'/')!=NULL)
				in_transl=FALSE;
			}
		while(strcmptrail(pinfo->line,10,NULL,0)==0 || 
			strcmptrail(pinfo->line,10,"FT",2)==0);
		}
	}
if(sgc!=0 && out_format->option != flat) {
	if(mark!=0) fputc(mark,out_format->fichier.outfile);
	fprintf(out_format->fichier.outfile,"Genetic code used: %s\n",
			get_code_descr(sgc));
	}
if(pha!=0 && out_format->option != flat) {
	if(mark!=0) fputc(mark,out_format->fichier.outfile);
	fprintf(out_format->fichier.outfile,
		"Reading frame shifted by %d bases.\n",pha);
	}
if(descri != NULL) {
	if(mark!=0) fputc(mark,out_format->fichier.outfile);
	else if(out_format->option == flat) {
		if(embl) fprintf(out_format->fichier.outfile,"CC   ");
		else fprintf(out_format->fichier.outfile,
				"COMMENT\n            ");
		}
	fprintf(out_format->fichier.outfile,
	    "Extracting fragment %s from above ",descri);
	if(deb_info_fille!=0) fprintf(out_format->fichier.outfile,"feature\n");
	else	fprintf(out_format->fichier.outfile,"sequence\n");
	}
	
if(out_format->option==flat) {
	if(embl) { /* cas EMBL */
		fprintf(out_format->fichier.outfile,"SQ   Sequence %d BP;\n",length);
		}
	else if(swissprot) { /* cas EMBL */
		fprintf(out_format->fichier.outfile,"SQ   SEQUENCE %d AA; ; ;\n",length);
		}
	else	{ /* cas GenBank */
		fprintf(out_format->fichier.outfile,"ORIGIN\n");
		}
	}
fin_header:
if(out_format->option == analseq) {
	int lastin;
	
	i=sgc;
	if(pha!=0 && debut==1) i=100*pha+i;
	padtosize(kin_record.normal_rec.name,mnemo,20);
	kin_record.normal_rec.pseq = out_format->fichier.analseq->num_se+1;
	kin_record.normal_rec.length=length;
	kin_record.normal_rec.phase=i;
	lastin = ++(out_format->fichier.analseq->num_in);
	erreur=dir_write(out_format->fichier.analseq->kin,lastin,1,&kin_record);
	if(erreur) return erreur;
	kin_record.first_rec.total=lastin;
	kin_record.first_rec.nbseq=lastin-1;
	kin_record.first_rec.espace=0;
	erreur=dir_write(out_format->fichier.analseq->kin,1,1,&kin_record);
	if(erreur) return erreur;
	}
else if(out_format->option == gcg) 
	fprintf(out_format->fichier.outfile,"..\n");
out_format->previous_len=0;
if(out_format->option != analseq) 
	erreur=ferror(out_format->fichier.outfile);
return erreur;
}


int prep_extract(out_option option, char *fname, extract_option choix,
	char *feature_name, char *bornes, char *min_bornes, char **message)
/* rend 0 ssi pas d'erreur 
min_bornes==NULL si pas besoin de description du fragment  minimal
*/
{
int exists;
FILE *test;
static char open_problem[]="File cannot be opened";
static char espace_problem[]="Ne pas utiliser de fichier de travail avec espace";
static char prot_DNA[]="Cannot extract protein sequence in DNA workfile";
static char DNA_prot[]="Cannot extract DNA sequence in protein workfile";

*message=NULL;
out_struct.option = option;
if(option == gcg ) {
	out_struct.width=50;
	out_struct.fichier.outfile=tmpfile();
	if(out_struct.fichier.outfile==NULL) goto Open_problem;
	}
else if(option == analseq ) {
	char tot_name[80], mode[3];
	out_struct.width=90;
	strcpy(tot_name,fname); strcat(tot_name,".in");
	test=fopen(tot_name,"r");
	if(test!=NULL) {
		exists=TRUE; fclose(test); strcpy(mode,"r+");
		}
	else	{
		exists=FALSE; strcpy(mode,"w+");
		}
	analseq_data.kin=dir_open(tot_name,NULL,mode,kin_size,kin_size);
	if(analseq_data.kin==NULL) goto Open_problem;
	strcpy(tot_name,fname); strcat(tot_name,".se");
	analseq_data.kse=dir_open(tot_name,NULL,mode,kse_size,kse_size);
	if(analseq_data.kse==NULL) { 
		dir_close(analseq_data.kin);
		goto Open_problem;
		}
	if(exists) {
		dir_read(analseq_data.kin,1,1,&kin_record);
		analseq_data.num_in = kin_record.first_rec.total;
		dir_read(analseq_data.kse,1,1,&kse_record);
		analseq_data.num_se = kse_record.total;
		if(kin_record.first_rec.espace != 0) {
			dir_close(analseq_data.kin);
			dir_close(analseq_data.kse);
			goto Espace_problem;
			}
		}
	else 	{
		kin_record.first_rec.total=1;
		kin_record.first_rec.nbseq=0;
		kin_record.first_rec.espace=0;
		dir_write(analseq_data.kin,1,1,&kin_record);
		kse_record.total=2;
		dir_write(analseq_data.kse,1,1,&kse_record);
		memset(&kse_record,' ',kse_size);
		if(nbrf || swissprot || choix==translate) 
			memcpy(kse_record.seq,"PROTEINES",9);
		else
			memcpy(kse_record.seq,"ACIDES NUCLEIQUES",17);
		dir_write(analseq_data.kse,2,1,&kse_record);
		analseq_data.num_in=1;
		analseq_data.num_se=2;
		}
	out_struct.fichier.analseq= &analseq_data;
	}
else	{
	out_struct.width=60;
	test=fopen(fname,"r");
	if(test != NULL) {
		exists=TRUE; fclose(test);
		}
	else	exists=FALSE;
	if(exists) {
		static char file_exists[]="File already exists";
		*message=file_exists;
		}
	out_struct.fichier.outfile=fopen(fname,"a");
	if(out_struct.fichier.outfile==NULL) goto Open_problem;
	if(nbrf && option == flat)
		fputs(NBRF_BACK_SLASHES, out_struct.fichier.outfile);
	}

if(choix==fragment || choix==region) { /* analyse de bornes et min_bornes */
	int err;
	err=extraire_bornes(bornes, &dlu, &flu, &xd, &xf);
	if(err) {
		static char pb_bornes[]="Incorrect endpoints";
		*message=pb_bornes;
		return 1;
		}
	frag_begin_to_end=(xd==0 && xf!=0);
	if( min_bornes==NULL ) {
		dlum=dlu; flum=flu; xdm=xd; xfm=xf;
		}
	else	{
		err=extraire_bornes(min_bornes, &dlum, &flum, &xdm, &xfm);
		if(err) {
			static char pb_min[]="Incorrect minimum endpoints";
			*message=pb_min;
			return 1;
			}
		}
	}

if(option==analseq && exists)  {
	dir_read(analseq_data.kse,2,1,&kse_record);
	if( (choix==translate) || swissprot) {
		if(strncmp(kse_record.seq,"PROTEINES",9) != 0) {
			dir_close(analseq_data.kin);
			dir_close(analseq_data.kse);
			goto Protseq_DNAfile;
			}
		}
	else	{
		if(strncmp(kse_record.seq,"PROTEINES",9) == 0 ) {
			dir_close(analseq_data.kin);
			dir_close(analseq_data.kse);
			goto DNAseq_Protfile;
			}
		}
	}
extract_data.choix=choix;
if(choix==feature || choix==region) {
	char *p;
	strcpy(extract_data.data.feature_name,feature_name);
	p=extract_data.data.feature_name - 1; while(*(++p)!=0) *p=toupper(*p);
	}
return 0;

/* traitement des erreurs */
Open_problem:
*message=open_problem;
return 1;

Protseq_DNAfile:
*message=prot_DNA;
return 1;

DNAseq_Protfile:
*message=DNA_prot;
return 1;

Espace_problem:
*message=espace_problem;
return 1;
} /* end of prep_extract */


int extract_1_seq(int numseq, char *bornes)
/* rend le nbre >=0 de seqs extraites par la fonction
ou -1 si erreur d'ecriture
*/
{
int erreur, longueur, pos, l, lmne, totfrags, debut_codon;
char mnemo[20], seq[101];
static const int echec= -1;

readsub(numseq); longueur = psub->length;
memcpy(mnemo,psub->name,16); mnemo[16]=0;
lmne=16; while(mnemo[--lmne]==' ') mnemo[lmne]=0; lmne++;
if(extract_data.choix == simple) {
	if(out_struct.option == fasta) {
		debut_codon = (psub->phase % 100) + 1;
		longueur = psub->length - debut_codon + 1;
		}
	else	debut_codon = 1;
	erreur = write_seq_header(mnemo, numseq, 1, longueur, &out_struct, 
			NULL, 0, 0);
	if(erreur) return echec;
	for(pos = debut_codon; pos <= longueur; pos += out_struct.width) {
		l = gfrag(numseq, pos, out_struct.width, seq);
		erreur = write_one_seq_line(&out_struct, seq, l);
		if(erreur) return echec;
		}
	if( write_seq_footer(&out_struct,mnemo,'N') ) return echec;
	return 1;
	}
else if(extract_data.choix == translate ) {
	char *protein, *p;
	int k;

	if(cds_type == 0) cds_type = fcode(ksmj, "04CDS", 20);
	if(psub->type != cds_type) return 0;
	protein = translate_cds(numseq);
	if(protein == NULL) return 0;
	longueur = strlen(protein);
	erreur = write_seq_header(mnemo,numseq,1,longueur,&out_struct,NULL,0,0);
	if(erreur) return echec;
	for(pos = 0; pos < longueur; pos += out_struct.width) {
		k = longueur - pos;
		if(k > out_struct.width) k = out_struct.width;
		memcpy(seq, protein + pos, k); seq[k] = 0;
		erreur = write_one_seq_line(&out_struct, seq, k);
		if(erreur) return echec;
		}
	if(out_struct.option == gcg) {
		p = strchr(mnemo, '.');
		if(p != NULL) *p = '-';
		strcat(mnemo, ".aa");
		}
	if( write_seq_footer(&out_struct, mnemo, 'P') ) return echec;
	return 1;
	}
else if(extract_data.choix == feature ) {
	long debut_info_mere, pft, previous;
	int div;
	
	if(psub->length==0 || psub->pext>0) return 0;
	seq_to_annots(numseq, &debut_info_mere, &div);
	totfrags=0;
	pft = debut_info_mere;
	if( read_annots(pft, div) == NULL) return totfrags;
	do	{
		previous = pft;
		next_annots(&pft);
		if( (genbank && strncmp(pinfo->line,"BASE COUNT",10) == 0 ) ||
			strncmp(pinfo->line,"SQ ",3)==0 ) return totfrags;
		}
	while( (genbank && strncmp(pinfo->line,"FEATURES",8) != 0 ) ||
		(embl && strncmp(pinfo->line,"FH ",3) != 0 ) ||
		(swissprot && strncmp(pinfo->line,"FT ",3) != 0) );
	if( swissprot ) pft = previous;
	while ( TRUE ) {
		if(embl || genbank)
			longueur = proc_location(&pft, div, 
					extract_data.data.feature_name, 
						numseq, debut_info_mere);
		else if (swissprot)
			longueur = proc_location_sw(&pft, div, 
					extract_data.data.feature_name, 
						numseq, debut_info_mere);
		else
			longueur = 0;
		if(longueur == 0) break;
		++totfrags;
		sprintf(mnemo+lmne,".F%d",totfrags);
		erreur =  write_seq_header(mnemo, -1, 1, longueur, &out_struct, 
			NULL, 0, 0);
		if(erreur) return echec;
		for(pos=1; pos<=longueur; pos += out_struct.width) {
			l=gfrag(-1,pos,out_struct.width,seq);
			erreur=write_one_seq_line(&out_struct,seq,l);
			if(erreur) return echec;
			}
		if( write_seq_footer(&out_struct,mnemo,'N') ) return echec;
		}
	return totfrags;
	}
else if(extract_data.choix == fragment) {
	int point, dm, fm, cpld, cplf, numd, numf, ld, lf, is,
		compl, deb, fin, debn, finn, div;
	long deb_info_fille;
		
	if(psub->pext>0) {
		seq_to_annots(numseq, &deb_info_fille, &div);
		readsub(numseq);
		readext(psub->pext); dm=pext->deb; fm=pext->fin;
		cpld= (dm > fm);
		numd=pext->mere; point=pext->next;
		readsub(numd); ld=psub->length;
		if( xd+xf != 0 && pext->next != 0 ) {
			do	{
				readext(point); fm=pext->fin;
				numf=pext->mere; point=pext->next;
				cplf= (pext->deb > pext->fin);
				}
			while(point);
			if(frag_begin_to_end && numd != numf) return 0;
			readsub(numf); lf=psub->length;
			}
		else if( xd+xf != 0 ) {
			cplf=cpld; numf=numd; lf=ld;
			}
		}
	else	{
		deb_info_fille=0;
		cpld=cplf=FALSE; dm=1; fm=ld=lf=longueur;
		numd=numf=numseq;
		}
	if( xd+xf == 0 ) {
		is=numd; compl=cpld;
		}
	else	{
		is=numf; compl=cplf;
		}
	deb=dlu; if(compl) deb = -dlu;
	if(xd==0) deb += dm;
	else deb += fm;
	fin=flu; if(compl) fin = -flu;
	if(xf==0) fin += dm;
	else fin += fm;
	debn=dlum-dlu; if(debn<0) debn=0;
	finn=abs(flum-flu);
	erreur=write_fragment(&out_struct,mnemo,is,deb,fin,debn,finn,bornes,
		deb_info_fille, div);
	if(erreur==2) return echec;
	else if(erreur==1) return 0;
	else return 1;
	}
else if(extract_data.choix == region ) {
	long debut_info_mere, pft, previous;
	int div;
	
	if(psub->length==0 || psub->pext>0) return 0;
	seq_to_annots(numseq, &debut_info_mere, &div);
	totfrags=0;
	pft = debut_info_mere;
	if( read_annots(pft, div) == NULL) return totfrags;
	do	{
		previous = pft;
		next_annots(&pft);
		if( (genbank && strncmp(pinfo->line,"BASE COUNT",10) == 0 ) ||
			strncmp(pinfo->line,"SQ ",3)==0 ) return totfrags;
		}
	while( (genbank && strncmp(pinfo->line,"FEATURES",8) != 0 ) ||
		(embl && strncmp(pinfo->line,"FH ",3) != 0 ) ||
		(swissprot && strncmp(pinfo->line,"FT ",3) != 0) );
	if( swissprot ) pft = previous;
	while ( TRUE ) {
		if(genbank || embl)
			longueur = proc_location(&pft, div, extract_data.data.feature_name, 
				numseq, debut_info_mere);
		else if(swissprot)
			longueur = proc_location_sw(&pft, div, extract_data.data.feature_name, 
				numseq, debut_info_mere);
		else
			longueur = 0;
		if(longueur == 0) break;
		longueur = mod_limites();
		if(longueur == 0) continue;
		++totfrags;
		sprintf(mnemo+lmne,".F%d",totfrags);
		erreur = write_fragment(&out_struct,mnemo,-1,1,longueur,0,0,
			bornes,0,0);
		if(erreur==1) continue;
		else if(erreur==2) return echec;
		}
	return totfrags;
	}
else	return 0;
} /* end of extract_1_seq */


static int write_seq_footer(out_format_struct *out_format, char *name, char Nuc_ou_Prot)
/* returns 0 for success, !=0 for write error */
{
int retval=0;

if(out_format->option==flat) {
	if(nbrf)
		fprintf(out_format->fichier.outfile,"///\n");
	else
		fprintf(out_format->fichier.outfile,"//\n");
	retval=ferror(out_format->fichier.outfile);
	}
else if(out_format->option == gcg) {
	retval=end_gcg(out_format,name, Nuc_ou_Prot);
	}
return retval;
} /* end of write_seq_footer */

	
static int extraire_bornes(char *bornes, int *debut, int *fin, int *xdebut , int *xfin)
/* interprete [e]+/-debut,[e]+/-fin
charge *debut, *fin
charge *xdebut!=0 ssi e a gauche
charge *xfin!=0 ssi e a droite
rend !=0 en cas d'erreur
*/
{
char *p, *px;
int tot;

p=bornes-1; while(*(++p)!=0) *p=toupper(*p);
compact(bornes);
p=strchr(bornes,',');
if(p==NULL) return 1;
px=strchr(bornes,'E'); if(px>p) px=NULL;
*xdebut = (px!=NULL);
if( px!=NULL && px==p-1 )
	*debut=0;
else	{
	if(px==NULL) px=bornes-1;
	tot=sscanf(px+1,"%d",debut);
	if(tot!=1) return 1;
	}
px=strchr(p+1,'E');
*xfin = (px!=NULL);
if( px!=NULL && *(px+1)==0) 
	*fin=0;
else	{
	if(px==NULL) px=p;
	tot=sscanf(px+1,"%d",fin);
	if(tot!=1) return 1;
	}
if( (*xdebut==0 && (*xfin==0 && *debut> *fin) ||
	(*xdebut!=0 && (*xfin==0 || (*xfin!=0 && *debut> *fin))))) return 1;
if(*xdebut==0 && *debut>0) (*debut)--;
if(*xfin==0 && *fin>0) (*fin)--;
return 0;
}


void fin_extract(void)
{
if(out_struct.option==analseq) {
	memset(&kin_record,0,kin_size);
	kin_record.first_rec.total=out_struct.fichier.analseq->num_in;
	kin_record.first_rec.nbseq=kin_record.first_rec.total-1;
	kin_record.first_rec.espace=0;
	dir_write(out_struct.fichier.analseq->kin,1,1,&kin_record);
	memset(&kse_record,0,kse_size);
	kse_record.total=out_struct.fichier.analseq->num_se;
	dir_write(out_struct.fichier.analseq->kse,1,1,&kse_record);
	dir_close(out_struct.fichier.analseq->kin);
	dir_close(out_struct.fichier.analseq->kse);
	}
else	{
	if( nbrf && out_struct.option == flat) 
		fputs(NBRF_BACK_SLASHES, out_struct.fichier.outfile);
	fclose(out_struct.fichier.outfile);
	}
}


static int end_gcg(out_format_struct *out_format, char *mnemo, char Nuc_ou_Prot)
/* rend !=0 si erreur */
{
time_t heure;
FILE *gcg_file;
int la, retval;
#define MAX_AUX 100
char *p, fname[40], line[MAX_AUX+1], aux[MAX_AUX+1], *fin, date[26];

strcpy(fname,mnemo);
p=fname-1; while(*(++p)!=0) *p=tolower(*p);
p=strchr(fname,'.');
if(p==NULL) strcat(fname,".seq");
p=strchr(fname,'.');
if(nbrf || swissprot || (p!=NULL && strcmp(p+1,"aa")==0)) {
	if(p!=NULL) strcpy(p,".pep");
	}
rewind(out_format->fichier.outfile);
gcg_file=fopen(fname,"w");
if(gcg_file==NULL) return 1;
while(strcmp(fgets(line,MAX_AUX,out_format->fichier.outfile),"..\n")!=0) {
/* ecriture du header
special attention needed for lines containing .. */
	la= -1;
	p=line;
	fin=line+strlen(line);
	do	{
		if(strncmp(p,"..",2)!=0)
			aux[++la]= *p;
		else	{
			memcpy(aux+la+1,". .",3);
			la+=3;
			p++;
			}
		p++;
		}
	while(p<fin);
	aux[++la]=0;
	fputs(aux,gcg_file);
	}
/* name + length + date + check line */
time(&heure);
strcpy(date,asctime(localtime(&heure)));
fprintf(gcg_file,
	"%s  Length: %d  %.6s, %.4s %.5s  Type: %c  Check: %d  ..\n",
	fname,out_format->previous_len,date+4,date+20,date+11,Nuc_ou_Prot,
	gcg_checksum % 10000 );
/* read & write the seq part */
while(fgets(line,MAX_AUX,out_format->fichier.outfile)!=NULL) {
	fputc('\n',gcg_file);
	fputs(line,gcg_file);
	}
retval=ferror(gcg_file);
fclose(gcg_file);
fclose(out_format->fichier.outfile);
out_format->fichier.outfile=tmpfile();
return retval;
#undef MAX_AUX
} /* end of end_gcg */


static int proc_locat_fragment(char *location, char *debut, char *fin, int tot_locat,
	int max_locat, long pinf, int div, int mere, int *part5, int *part3)
/* traitement d'un exon d'une location, avec cas envoi avec access # ou avec
label.
location, debut, fin: location traitee, extremites de sa partie traitee,
tot_locat: longueur utilisee de location,
max_locat: taille max de location,
pinf: debut annotation,
mere: # mere traitee,
part5, part3: mis a TRUE si partiel
rend 0 si ok, sinon code d'erreur
*/
{
int compl, point, next, nacc, a, b, err, expr_a_b, sv_mere, same_mere;
char *p, access[11], label[16];

sv_mere=mere;
boulot:
while( *debut==' ' && debut<=fin) debut++;
if(debut>fin) return 2;
if(strncmp(debut,"COMPLEMENT(",11)==0) {
	debut+=11;
	p=nextpar(debut-1)-1; if(p>fin) return 2;
	fin=p-1;
	compl=TRUE;
	}
else compl=FALSE;
/* recherche renvoi a autre banque */
if((p=strstr(debut,"::"))!=NULL && p<fin) return 3;
p=strchr(debut,':');
if(p!=NULL && p<fin) {
/* renvoi a autre seq par access # */
	memcpy(access,debut,p-debut);
	access[p-debut] = 0;
	if( (debut = strchr(access,'.')) != NULL) *debut = 0;
	debut = p + 1;
	nacc = fcode(kacc, access, ACC_LENGTH);
/* chercher une seq ayant cet acc # en primary ou seq unique */
	if(nacc == 0) return 5;
	readacc(nacc);
	point = pacc->plsub;
	if(point == 0) return 5;
	readshrt(point);  
	if(pshrt->next == 0) { /* une seule seq associee a access */
		mere = pshrt->val;
		}
	else	{ /* plusieurs seqs associees, chercher celle dont primary */
		while(TRUE) {
			if(point == 0) return 5;
			readshrt(point); next = pshrt->next; 
			mere = pshrt->val;
			readsub(mere); readloc(psub->plinf);
			readshrt(ploc->placc);
			if(pshrt->val == nacc) break;
			point = next;
			}
		}
	}
else	{
/* on reste dans la seq de la table de features */
	*access = 0;
	}
/* traitement partie location */
p=strstr(debut,".."); if(p!=NULL && p>=fin) p=NULL;
if(p!=NULL) 
	expr_a_b=TRUE;
else	{
	memcpy(label,debut,fin-debut+1); label[fin-debut+1]=0;
	if((p=strchr(label,'<'))!=NULL) *p=' ';
	if((p=strchr(label,'>'))!=NULL) *p=' ';
	err=sscanf(label,"%d",&a);
	if(err==1) { /* residu unique */
		*part5= *part3 = FALSE; /* ??? */
		b=a; expr_a_b=FALSE;
		}
	else	{ /* c'est un label */
		err = find_label(label,(*access==0 ? NULL : access), &pinf, &div,
				 &mere);
		if(err==1) return 5;
		else if(err==2) return 4;
		debut=location+tot_locat+1;
		err = read_location(pinf, div, debut, max_locat - tot_locat - 1);
		if(err==0) return 8;
		fin=debut+err-1;
		goto boulot;
		}
	}
if(expr_a_b) {
	if(*debut=='<' || *debut=='>') {
		if(compl) *part3=TRUE;
		else *part5=TRUE;
		debut++;
		}
	err=sscanf(debut,"%d",&a);
	if(err==0) return 2;
	p+=2;
	if(*p=='<' || *p=='>') {
		if(compl) *part5=TRUE;
		else *part3=TRUE;
		p++;
		}
	err=sscanf(p,"%d",&b);
	if(err==0) return 2;
	}
if(simext.nfrags>=simext.nfmax) return 2;
++simext.nfrags;
readsub(mere);
readloc(psub->plinf);
/* check fragment endpoints relatively to size of mother seq */
if(a<1 || b<1 || a>psub->length || b>psub->length) return 7;
if(compl) {
	next=a; a=b; b=next;
	}
simext.valext[simext.nfrags-1][0]=a;
simext.valext[simext.nfrags-1][1]=b;
simext.valext[simext.nfrags-1][2]=ploc->pnuc;
simext.valext[simext.nfrags-1][3]=mere;
if(mere==sv_mere) same_mere=TRUE;
return 0;
}
/* end of proc_locat_fragment */


static int proc_location(long *pinf, int div, char *nom_feature, int mere, int debut_info_mere)
/* parcours de la table des FT a partir de *pinf a la recherche de feature
nom_feature et cree la sous-seq virtuelle correspondante a la premiere rencontre
*pinf: debut ou lire annotations
nom_feature: nom de la feature cherchee
mere: # de la mere traitee
debut_info_mere: debut annot de la mere traitee
rend longueur de la sous-seq virtuelle fabriquee ou 0 si aucune trouvee
*/
{
int llocat, compl, i, j, err, cherche_v, length;
char location[5000], *debut, *fin, *p, *q, *v;


if( read_annots(*pinf, div) == NULL) return 0;
while(TRUE) {
	next_annots(pinf);
	if(strcmptrail(pinfo->line,10,NULL,0)==0 || 
		strcmptrail(pinfo->line,10,"FT",2)==0 ) continue;
	if(embl && strncmp(pinfo->line,"FH",2)==0) continue;
	if(pinfo->line[0]!=' ' && strncmp(pinfo->line,"FT",2)!=0) {
		length = 0;
		break;
		}
	p=pinfo->line+4;
	while(*(++p)!=' ') *p=toupper(*p);
	if(strcmptrail(pinfo->line+5,15,nom_feature,strlen(nom_feature))!=0) continue;
	llocat=read_location(*pinf,div,location,sizeof(location));
	if(llocat==0) continue;
	psimext= &simext;
	simext.nfrags=0;
	simext.nfmax=SIMEXT_NFMAX;
	simext.newsim=TRUE;
	simext.siminf= (int)(*pinf);
	simext.pinfme = debut_info_mere;
	simext.div = div;
	debut=location; fin=debut+llocat-1;
	compl=FALSE;
	if(strncmp(location,"COMPLEMENT(JOIN(",16)==0) {
		compl=TRUE;
		debut+=11;
		}
	if(strncmp(debut,"JOIN(",5)==0) {
		fin=nextpar(debut+4);
		debut+=5;
/* remove redundant join() within a previous join() */
		do	{
			p=strstr(debut,"JOIN("); if(p>fin) p=NULL;
			if(p==NULL) break;
			q=nextpar(p+4);
			memset(p,' ',5); *q=' ';
			}
		while(TRUE);
		cherche_v=TRUE;
		}
	else	{
		v=fin+1;
		cherche_v=FALSE;
		}
	do	{
		if(cherche_v) {
			v=debut;
			do	{
				v++;
				if(*v==',' || *v==')') break;
				if(*v=='(') v=nextpar(v);
				}
			while(TRUE);
			}
		err=proc_locat_fragment(location,debut,v-1,llocat,sizeof(location),
			*pinf,div,mere,&i,&i);
		debut=v+1;
		cherche_v=TRUE;
		}
	while(err==0 && v<fin);
	if(err!=0) continue;
	length=0;
	for(i=0; i<simext.nfrags; i++)
		length+=abs(simext.valext[i][1]-simext.valext[i][0])+1;
	if(compl) {
/* fin traitement du complement( externe */
		for(i=0; i<simext.nfrags/2; i++) 
			for(j=0; j<4; j++) {
				err=simext.valext[i][j];
				simext.valext[i][j]=simext.valext[simext.nfrags-i-1][j];
				simext.valext[simext.nfrags-i-1][j]=err;
				}
		for(i=0; i<simext.nfrags; i++) {
			err=simext.valext[i][0];
			simext.valext[i][0]=simext.valext[i][1];
			simext.valext[i][1]=err;
			}
		}
	break;
	}
return length;
} /* end of proc_location */


static int proc_location_sw(long *pinf, int div, char *nom_feature, int mere, int debut_info_mere)
/* parcours de la table des FT a partir de *pinf a la recherche de feature
nom_feature et cree la sous-seq virtuelle correspondante a la premiere rencontre
*pinf: debut ou lire annotations
nom_feature: nom de la feature cherchee
mere: # de la mere traitee
debut_info_mere: debut annot de la mere traitee
rend longueur de la sous-seq virtuelle fabriquee ou 0 si aucune trouvee
*/
{
int length, lnom, debut, fin;
char *p, *q;

lnom = strlen(nom_feature);
if( read_annots(*pinf, div) == NULL) return 0;
while(TRUE) {
	next_annots(pinf);
	if(strcmptrail(pinfo->line,10,"FT",2) == 0 ) continue;
	if(strncmp(pinfo->line,"FT",2) != 0) {
		length = 0;
		break;
		}
	p=pinfo->line+4;
	while( *(++p) != ' ') *p = toupper(*p);
	if(strcmptrail(pinfo->line + 5, lnom + 1, nom_feature, lnom) != 0) continue;
	psimext = &simext;
	simext.nfrags=1;
	simext.nfmax=SIMEXT_NFMAX;
	simext.newsim=TRUE;
	simext.siminf= (int)(*pinf);
	simext.pinfme = debut_info_mere;
	simext.div = div;
	pinfo->line[27] = 0; /* enlever les < ou > des FROM et TO */
	q = p - 1;
	while( *(++q) != 0 ) if( *q == '>' || *q == '<') *q = ' ';
	debut = fin = -1000;
	sscanf(p, "%d%d", &debut, &fin);
	if(debut == -1000 || fin == -1000) continue; /* lecture FROM ou TO impossible */
	simext.valext[0][0] = debut;
	simext.valext[0][1] = fin;
	length = simext.valext[0][1] - simext.valext[0][0] + 1;
	readsub(mere);
	readloc(psub->plinf);
	simext.valext[0][2]=ploc->pnuc;
	simext.valext[0][3]=mere;
	break;
	}
return length;
} /* end of proc_location_sw */


static int lim1(int b, int l, int lm, int lseq)
{
if(b+l>lseq || b+lm<=0) return 0;
if(b+l>0) return b+l;
else return 1;
}

static int lim2(int b, int l, int lm, int lseq)
{
if(b+l<=0 || b+lm>lseq) return 0;
if(b+l<=lseq) return b+l;
else return lseq;
}

static int mod_limites(void)
/* calcul des limites du fragment a extraire.
Rend 0 si fragment ne peut pas etre extrait;
sinon rend sa longueur
*/
{
int pmere, bd, point, bf, p2, d, f, compl;

if(xd==0) {
	pmere=simext.valext[0][2];
	bd=simext.valext[0][0];
	point=1;
	}
else	{
	pmere=simext.valext[simext.nfrags-1][2];
	bd=simext.valext[simext.nfrags-1][1];
	point=simext.nfrags;
	}
if(xf==0) {
	p2=simext.valext[0][2];
	bf=simext.valext[0][0];
	point=1;
	}
else	{
	p2=simext.valext[simext.nfrags-1][2];
	bf=simext.valext[simext.nfrags-1][1];
	point=simext.nfrags;
	}
if(pmere!=p2) return 0;
readsub(simext.valext[point-1][3]);
compl= ( simext.valext[point-1][0] > simext.valext[point-1][1] );
if(!compl) {
	d=lim1(bd,dlu,dlum,psub->length);
	f=lim2(bf,flu,flum,psub->length);
	}
else	{
	d=lim2(bd,-dlu,-dlum,psub->length);
	f=lim1(bf,-flu,-flum,psub->length);
	}
if( (d==0) || (f==0) || (!compl && d>f) || (compl && f>d) ) return 0;
simext.newsim=TRUE;
simext.nfrags=1;
simext.valext[0][0]=d;
simext.valext[0][1]=f;
simext.valext[0][2]=pmere;
return abs(f-d)+1;
}


static int write_fragment(out_format_struct *out_format, char *mnemo, int iseq, 
	int deb, int fin, int debn, int finn, char *description, long deb_info_fi, int div)
/* ecrit un fragment de deb a fin de seq # iseq prolonge par des N a debn->finn
rend 0 ssi ok; 1 si fragment non extractible, 2 si erreur ecriture dans fichier
*/
{
int dmod, fmod, dn, fn, compl, rb, i, l, r, width, doit, nbn, kd, erreur;
char lu[100];

dmod=deb; fmod=fin; dn=fn=0;
width=out_format->width;
compl= (deb>fin);
if(iseq != -1) {
	readsub(iseq);
	if(!compl) {
		if(deb<=0) {
			dn= -deb+1; dmod=1;
			if(dn>debn) return 1;
			}
			else dn=0;
		if(fin>psub->length) {
			fn=fin-psub->length; fmod=psub->length;
			if(fn>finn) return 1;
			}
		else fn=0;
		}
	else	{
		if(fin<=0) {
			fn= -fin+1; fmod=1;
			if(fn>finn) return 1;
			}
		else fn=0;
		if(deb>psub->length) {
			dmod=psub->length; dn=deb-psub->length;
			if(dn>debn) return 1;
			}
		else dn=0;
		}
	}
erreur=write_seq_header(mnemo,iseq,deb,fin,out_format,description,deb_info_fi,div);
if(erreur) return 2;
/* traitement du debut si existence de N */
doit=TRUE;
if(dn != 0) {
	rb = dn % width;
	memset(lu,'N',width);
	for(i=width; i<=dn; i+=width ) {
		erreur=write_one_seq_line(out_format,lu,width);
		if(erreur) return 2;
		}
	if(compl) {
		l=width-rb; i=dmod-fmod+1;
		if(l>i) l=i;
		l=gfrag(iseq,dmod-l+1,l,lu+rb);
		complementer_seq(lu+rb,l);
		r=rb+l;
		dmod -= l;
		if(dmod<fmod) doit=FALSE;
		else	{
			erreur=write_one_seq_line(out_format,lu,width);
			if(erreur) return 2;
			}
		}
	else	{
		l=width-rb; i=fmod-dmod+1;
		if(l>i) l=i;
		l=gfrag(iseq,dmod,l,lu+rb);
		r=rb+l;
		dmod += l;
		if(dmod>fmod) doit=FALSE;
		else	{
			erreur=write_one_seq_line(out_format,lu,width);
			if(erreur) return 2;
			}
		}
	}
/* ecriture du reste de la sequence */
if(doit) {
	l=width;
	if(compl) {
		for(i=dmod; i >= fmod+l-1; i -= l) {
			kd=i-width+1;
			r=gfrag(iseq,kd,width,lu);
			complementer_seq(lu,r);
			erreur=write_one_seq_line(out_format,lu,r);
			if(erreur) return 2;
			}
		}
	else	{
		for(i=dmod; i<= fmod-l+1; i+=l) {
			kd=i;
			r=gfrag(iseq,kd,width,lu);
			erreur=write_one_seq_line(out_format,lu,r);
			if(erreur) return 2;
			}
		}
	r= (abs(fmod-dmod)+1) % width;
	nbn=0;
	if(r==0) goto final;
	if(compl) {
		r = gfrag(iseq,fmod,r,lu);
		complementer_seq(lu,r);
		}
	else r = gfrag(iseq,fmod-r+1,r,lu);
	}
l=r+fn;
if(l>width) l=width;
memset(lu+r,'N',l-r);
if(r+fn < l) memset(lu+r+fn,' ',l-r-fn);
erreur=write_one_seq_line(out_format,lu,l);
if(erreur) return 2;
nbn=width-r;

final:
if(nbn<fn) {
	nbn=fn-nbn;
	memset(lu,'N',width);
	for(i=1; i<=nbn; i+=width) {
		r = width;
		if(i + r - 1 > nbn) r = nbn - i + 1;
		erreur=write_one_seq_line(out_format,lu,r);
		if(erreur) return 2;
		}
	}
erreur=write_seq_footer(out_format,mnemo,'N');
if(erreur) return 2;
return 0;
} /* end of write_fragment */
