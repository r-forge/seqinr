/* connectindex 

trois modes de fonctionnement:
1) update mode:
charger SMJYT avec les noms des divisions (au format avec libelle)
determiner les nouveautes dans une release par rapport a une banque acnuc
connecter ce qui est inchange aux nouveaux fichiers plats
produire fichiers:
	--.1  --.2  listes des noms des nouvelles et modifiees
	--.address  noms et adresses des nouvelles et modifiees
	--.lost     listes des noms des meres a passer dans processft
	disparu.mne listes des noms des disparues
	
2) install mode:
charger SMJYT avec les noms des divisions (au format avec libelle)
connecter les seqs trouvees dans les divisions avec les fichiers index,
mettre meres et filles non trouvees dans la liste des seqs non valides
mettre a jour en-tete de HELP et HELP_WIN

3) scan mode:
lire une ou plusieurs divisions, et remettre les pointeurs des meres
et filles inchangees vers leurs valeurs dans ces divisions.
Ce mode est prevu pour traiter une partie de la banque. Les 2 autres
traitent la banque en entier.

accepte format plat et format GCG

environnement: acnuc vers indexs
	       gcgacnuc: vers dir des nouveaux fichiers plats

les anciennes annotations (plats ou non) sont inutiles!

Le mode scan s'obtient par des arguments:
connectindex -<nbre de divisions> noms-des-divisions

Dialogue:
i  ou  u  (install / update)
f  ou  g  (flat / gcg)
si mode update:
	root-of-outfiles
si mode install:
	y ou n  (create new seqs division?)
	si y:
		new-seqs-div-name
# of divisions
division names

*/

#include "dir_acnuc.h"
#include <ctype.h>

FILE *newfile, *modfile, *addrfile, *dispfile;
int *blist, *invalidlist, *subpnuc, *suplist, *proclist, *auxlist,
	subseq_sorted, tloc, update_mode, count_scan = 0;
int *nom2loc, tot_nom2loc;
char *lismne, **newdivname;
#define MAXLOCAT 10000
#define MAXTYPES 100
char *list_types[MAXTYPES];
int tottypes;
typedef struct {
	int numero_fille;
	int newpinf;
	char *description;
	} fille_data;
#define maxfillespermere 10000
fille_data tab_fille_data[maxfillespermere+1];


/* prototypes */
void prepnewdiv(int totdiv);
void cre_new_if_needed(void);
int ask_div_names(void);
void prnumloc(void);
int numloc(char *name);
int proc_division(char *divname, int divnum, char *prefix, int *num_meres,
	int *num_filles, double *tot_long);
int read_all_types(void);
int compare_ft(FILE *flat, char *line, int lline, int numseq, int *changed);
void update_pinf_filles(int divnum, int totfilles);
char **prep_subdescr(char *acnuc_env);
static int flat_read_location(FILE *flat,char *line, int lline, 
	char *location,int maxlocat, long *newpos);
char *prepextract(int seqnum, int mere);
void preplocat(char *location, int maxlen);
static char *acctomere(char *access);
static char *nextpar(char *debut, char *fin);
void calc_gbk_date(char *gbk, char *newdate);
void calc_nbrf_date(FILE *flat, char *line, char *newdate);
char *fulldivname(char *divname, char *prefix);
void compute_pnuc(FILE *flat, char *line, int *newpnuc, size_t lline);
int fastfindacc(char *access);
static int n_r_compar(const void *p1, const void *p2);
void proc_help_file(char *name, int totmeres, double totlong, 
	int totfilles, int totbib, char *acnuc_env);
char *paqde3(double val);
void get_acnuc_dirs(char *acnuc_env, char *gcgacnuc_env);
#ifdef vms
void skipseq(char *nameline, FILE *fich, int cur_pos);
#endif
#ifdef unix
int proc_options(int argc, char **argv);
void write_quick_meres(void);
#endif


/* extern prototypes */
extern int poffset(int div, int offset);
extern void cre_new_division(char *fname);
extern FILE **divseq;
extern int *seqopened;
extern void *mycalloc(int n, size_t taille);


#ifdef __INTEL__
void attendre(void);
void attendre(void)
{
char rep[10];
printf("Hit <RETURN> key"); gets(rep);
}
char *mygets(char *ligne);
char *mygets(char *ligne)
{
char *p;
gets(ligne);
p = ligne + strlen(ligne) - 1;
while(*p == '\r' || *p == '\n') *(p--) = 0;
return ligne;
}
#define gets(p) mygets(p)
#endif


int main(int argc, char *argv[])
{
int lrootname, totlu, newtotdiv, newdivnum, tot, num, nfi, i,
	num_meres, num_filles, tot_bib;
double tot_long;
char rootname[100], *p, fname[100], name[17], mere[17], gcgacnuc_env[100],
	acnuc_env[100], new_seq_divname[40] = "";

#ifdef __INTEL__
atexit(attendre);
printf("Opening acnuc index files\n");
#endif
/* get acnuc + gcgacnuc environment */
get_acnuc_dirs(acnuc_env, gcgacnuc_env);
#ifdef __INTEL__
printf("from directory: %s\n", acnuc_env);
#endif
dir_acnucopen("WP");
count_scan = 0;
#ifdef unix
count_scan = proc_options(argc, argv);
#endif

tloc=read_first_rec(kloc,NULL);
printf("lenw=%d tloc=%d nseq=%d\n", lenw, tloc, nseq);

blist=(int *)mycalloc(lenw,sizeof(int));
invalidlist=(int *)mycalloc(lenw,sizeof(int));
auxlist=(int *)mycalloc(lenw,sizeof(int));
proclist=(int *)mycalloc(lenw,sizeof(int));
suplist=(int *)mycalloc(lenw,sizeof(int));
subpnuc=(int *)mycalloc(nseq+1,sizeof(int));
lismne = (char *)mycalloc(tloc, L_MNEMO);
nom2loc = (int *)mycalloc(tloc, sizeof(int));

update_mode = FALSE;
/* 
mode install = seqs trouvees sont traitees, autres seqs sont eliminees
mode update = install + cree liste nouvelles, disparues, modifiees, address
*/
if(count_scan == 0) {
#ifndef __INTEL__
	printf("Mode install(i) or update(u) ?");
	gets(fname);
	update_mode= (*fname=='u' || *fname=='U');
#else
	update_mode = FALSE;
#endif
	
	printf("Use flat(f) or gcg(g) format? ");
	gets(fname);
	flat_format= (*fname!='g' && *fname!='G');
	}

if(update_mode) {
	printf("Output filename (without extension)? ");
	gets(rootname);
	p=strchr(rootname,'.'); if(p!=NULL) *p=0;
	lrootname=strlen(rootname);
	strcpy(fname,rootname);
	strcat(fname,".1");
	printf("Names of new seqs go to %s\n",fname);
	newfile=fopen(fname,"w");
	if(newfile==NULL) {
		fprintf(stderr,"Error opening %s\n",fname);
		exit(ERREUR);
		}
	strcpy(fname,rootname);
	strcat(fname,".2");
	printf("Names of modified seqs go to %s\n",fname);
	modfile=fopen(fname,"w");
	if(modfile==NULL) {
		fprintf(stderr,"Error opening %s\n",fname);
		exit(ERREUR);
		}
	strcpy(fname,rootname);
	strcat(fname,".address");
	addrfile=fopen(fname,"w");
	if(addrfile==NULL) {
		fprintf(stderr,"Error opening %s\n",fname);
		exit(ERREUR);
		}
	}
else if(count_scan == 0)	{
#ifndef __INTEL__
	printf("Create a `new sequences' division? (y/n) ");
	gets(new_seq_divname);
#else
	*new_seq_divname = 'n';
#endif
	if(*new_seq_divname == 'y' || *new_seq_divname == 'Y') {
		printf("Name of this division? ");
		gets(new_seq_divname);
		}
	else	*new_seq_divname=0;
	}

/* read all division names */ 
if(count_scan == 0) newtotdiv = ask_div_names();

/* read all types */
tottypes=read_all_types();

printf("Reading all seq names in memory..."); fflush(stdout);
lngbit(2,blist); /* liste des meres */
/* invalidlist: all bits set to 1 */
memset(invalidlist, 0, lenw * sizeof(int) );
non(invalidlist, invalidlist, lenw);
for(i = nseq+1; i <= lenw*lmot; i++) bit0(invalidlist, i);

prnumloc();
printf("\n"); fflush(stdout);

/* create all divisions */
if(count_scan == 0) prepnewdiv(newtotdiv);
else	cre_new_if_needed();

totlu = num_meres = num_filles = 0;
tot_long = 0;
if(count_scan != 0) {
	for(i=0; i < count_scan; i++) {
		for(newdivnum = 0; newdivnum <= divisions; newdivnum++) {
			if(strcmp(gcgname[newdivnum], newdivname[i]) == 0) 
				break;
			}
		totlu += proc_division(newdivname[i], newdivnum, gcgacnuc_env,
				&num_meres, &num_filles, &tot_long);
		}
	}
else	{
	if(update_mode) printf("New or modified seqs found:\n");
	for(newdivnum=0; newdivnum<=newtotdiv; newdivnum++) {
		totlu += proc_division(newdivname[newdivnum], newdivnum, 
			gcgacnuc_env, &num_meres, &num_filles, &tot_long);
		}
	if(update_mode) {
		fclose(newfile); fclose(modfile); fclose(addrfile);
		}
	}

printf("Total # of seqs found in flat files=%d\n",totlu);

if(genbank || embl) {
	printf("Updating NUCLEOT pointers of EXTRACT\n"); fflush(stdout);
	tot = read_first_rec(kext,NULL);
	for(num=2; num<=tot; num++) {
		readext(num);
		if(pext->mere>=2 && pext->mere <= nseq && 
			subpnuc[pext->mere] != 0) {
				pext->pnuc = subpnuc[pext->mere];
			writeext(num);
			}
		}
	dir_flush(kext);
	}

if(update_mode) {
	printf("Building list of disappeared sequences\n"); fflush(stdout);
	dispfile=fopen("disparu.mne","w");
	num=1;
	while((num=irbit(blist,num,nseq))!=0) {
		readsub(num);
		fprintf(dispfile,"%.16s\n",psub->name);
		}
	fclose(dispfile);

	if(genbank || embl) {
		/* suplist = modif + disparues */
		ou(suplist,suplist,blist,lenw);
	/* 
	recherche des meres principales des filles a supprimer et qui ne sont
	pas supprimees
	*/
		num=1;
		while((num=irbit(suplist,num,nseq))!=0) {
			readsub(num);
			if(psub->pext>=0) continue;
			memcpy(name,psub->name,16); 
			p=name+16; *p=0; while(*(--p)==' ') *p=0;
			lngbit(-psub->pext,auxlist);
			nfi=1;
			while((nfi=irbit(auxlist,nfi,nseq))!=0) {
				readsub(nfi);
				p=strchr(psub->name,'.');
				if(p==NULL ||p>=psub->name+16) continue;
				memcpy(mere,psub->name,p-psub->name);
				mere[p-psub->name]=0;
				if(strcmp(mere,name)==0) continue;
				i=isenum(mere);
				bit1(proclist,i);
				}
			}
		/* proclist=meres avec filles perdues et non rechargees */
		non(suplist,suplist,lenw);
		et(proclist,proclist,suplist,lenw);
		printf("Building list of sequences having lost subsequences\n");
		fflush(stdout);
		strcpy(fname,rootname);
		strcat(fname,".lost");
		dispfile=fopen(fname,"w");
		num = 1;
		while((num=irbit(proclist,num,nseq))!=0) {
			readsub(num);
			fprintf(dispfile,"%.16s\n",psub->name);
			}
		fclose(dispfile);
		}
	} /* end of if update_mode */
else if(count_scan == 0) { /* cas install mode */
	int lastlng, point, avant, apres, difference;
/* write list of invalid seqs at rec=3 of LONGL */
	lngbit(3, auxlist);
	avant = bcount(auxlist, nseq);
	apres = bcount(invalidlist, nseq);
	difference = apres - avant;
	printf("Writing new list of invalid sequences\n");
	num=0; i=0; point=3;
	memset(plng,0,lrlng);
	lastlng = read_first_rec(klng,NULL);
	plng->next= ++lastlng;
	while( (num=irbit(invalidlist,num,nseq)) != 0 ) {
		if(i >= SUBINLNG ) {
			writelng(point); point=plng->next;
			memset(plng,0,lrlng);
			plng->next = ++lastlng;
			i=0;
			}
		plng->sub[i++]=num;
		}
	plng->next=0;
	writelng(point);
	write_first_rec(klng, lastlng - 1, 0);
	dir_flush(klng);
	if(difference != 0) printf("%d seqs & subseqs from index files "
		"were not found in data files\n", difference);
/* update files HELP and HELP_WIN */
	printf("Updating statistics in Help files\n");
	tot_bib=read_first_rec(kbib,NULL)-1;
	proc_help_file("HELP", num_meres, tot_long, num_filles, 
				tot_bib, acnuc_env);
	proc_help_file("HELP_WIN", num_meres, tot_long, num_filles, 
				tot_bib, acnuc_env);
	if( (int)strlen(new_seq_divname) > 0 ) {
/* create the division of new sequences */
		printf("Creating division %s for new sequences\n",
			new_seq_divname);
		cre_new_division(new_seq_divname);
		dir_flush(ksmj);
		dir_flush(ktxt);
		}
	}

dir_acnucclose();
#ifdef unix
if( (!update_mode) && count_scan == 0 ) write_quick_meres();
#endif
printf("Normal end\n");
return 0;
}


int read_all_types(void)
{
/* charger list_types avec noms de tous les types
return=nbre de types (from 0)
*/
int totrec, tot, num;
char *p;

totrec=read_first_rec(ksmj,NULL);
tot= -1;
for(num=2; num<=totrec; num++) {
	readsmj(num);
	if(strncmp(psmj->name,"04",2)!=0) continue;
	tot++;
	list_types[tot]=(char *)mycalloc(21,1);
	memcpy(list_types[tot],psmj->name+2,18);
	list_types[tot][18]=0;
	p=strchr(list_types[tot],' ');
	if(p!=NULL) *p=0;
	}
return tot;
}


void prepnewdiv(int totdiv)
/* effacer toutes les anciennes div et ecrire les nouvelles */
{
int i, totsmj, finsorted, div;

/* effacer toutes les anciennes divisions */
for(div=0; div<=divisions; div++) {
	if( annotopened[div] ) fclose(divannot[div]);
	if(gcgname[div] != NULL) free(gcgname[div]);
	}
free(gcgname); gcgname = NULL;
free(annotopened); annotopened = NULL; free(divannot); divannot = NULL;
if(!flat_format) {
	for(div=0; div<=divisions; div++) {
		if( seqopened[div] ) fclose(divseq[div]);
		}
	free(seqopened); seqopened = NULL; free(divseq); divseq = NULL;
	}
divisions= -1;

totsmj = read_first_rec(ksmj, &finsorted);
for(i=2; i<=totsmj; i++) {
	readsmj(i);
	if(strncmp(psmj->name,"06",2)!=0) continue;
	memset(psmj,0,lrsmj);
	memset(psmj->name,'x',20);
	writesmj(i);
	if(finsorted >= i) finsorted = i-1;
	}
write_first_rec(ksmj, totsmj, finsorted);

for(i=0; i<=totdiv; i++) {
	cre_new_division(newdivname[i]);
	}
dir_flush(ksmj);
dir_flush(ktxt);
}


void cre_new_if_needed(void)
{
int num, i;

for(num = 0; num < count_scan; num++) {
	for(i = 0; i <= divisions; i++) {
		if(strcmp(gcgname[i], newdivname[num]) == 0) break;
		}
	if(i > divisions) cre_new_division(newdivname[num]);
	}
}


int ask_div_names(void)
/* demander les nombre et noms des divisions, les memoriser.
return: # de div (from 0)
*/
{
char fname[20], *p;
int totdiv, i;

printf("# of divisions in database? "); gets(fname);
sscanf(fname,"%d",&totdiv);
printf("Enter the %d division names\n",totdiv);
totdiv--;
newdivname = (char **)mycalloc(totdiv+1,sizeof(char *));
for(i=0; i<=totdiv; i++) {
	printf("? "); gets(fname);
	p=strchr(fname,'.'); if(p!=NULL) *p=0;
	newdivname[i]=(char *)mycalloc(1,strlen(fname)+1);
	strcpy(newdivname[i],fname);
	}
return totdiv;
}


int proc_division(char *divname, int divnum, char *prefix, int *num_meres,
	int *num_filles, double *tot_long)
{
static char gbkdeb[]="LOCUS";
static char embldeb[]="ID ";
static char nbrfdeb[]="ENTRY ";
FILE *flat, *seq_fich;
char *fname, line[300], *target, *p, *q, flatname[20], newdate[9],
	gcgtarget[40], gcgtarget2[40], gcgline[2000], *gbk_date;
int newpinf, totlu=0, ltarget, nloc, new, modif, newlong, old1, old2, 
	newpnuc, totfilles, l, l2, i, mere, lnameline;
off_t fpos;

fname = fulldivname(divname, prefix);
flat = fopen(fname,
#ifdef __INTEL__
	"rb"
#else
	"r"
#endif
	);
if(flat == NULL) {
	printf("Warning division file %s not found\n",fname);
	return 0;
	}
setvbuf(flat, NULL, _IOFBF, 51200);
if(!flat_format) {
	l=strlen(fname);
	strcpy(fname+l-3,"seq");
	seq_fich=fopen(fname,
#ifdef __INTEL__
		"rb"
#else
		"r"
#endif
		);
	if(seq_fich==NULL) {
		fprintf(stderr,"Error: sequence file %s not found\n",fname);
		exit(ERREUR);
		}
	}
if(embl || swissprot) {
	target = embldeb;
	}
else if(genbank) {
	target=gbkdeb;
	}
else if(nbrf) {
	target=nbrfdeb;
	}
ltarget=strlen(target);
	
printf("Processing division %s\n",fname); fflush(stdout);
if(update_mode) fprintf(addrfile,"%s\n",divname);
#ifdef vms
lnameline=MINRECSIZE+1;
#else
lnameline=sizeof(gcgline);
#endif
while(TRUE) {
	fpos = ftello(flat);
	newpinf = fpos;
	p = fgets(line, sizeof(line), flat);
	if(p == NULL) {
		if(update_mode) fprintf(addrfile, "//\n");
		break;
		}
	if(strncmp(line, target, ltarget) != 0) continue;
	if(genbank) {
		char *gbk_name;
		newlong = decode_locus_rec(line, &gbk_name, NULL, NULL, NULL,
			&gbk_date);
		strcpy(flatname, gbk_name);
		}
	else {
		if(embl || swissprot) p = line + 5;
		else if(nbrf) {
			p=line+ltarget;
			while( *p==' ' ) p++;
			}
		q=strchr(p,' ');
		memcpy(flatname,p,q-p); flatname[q-p]=0;
		}
	if(!flat_format) {
/* recherche de la seq dans le fichier .seq pour format gcg 
On cherche soit   ">>>>nom "  soit   ">>>>nom_0" qui est utilise dans
le format GCG pour les longues sequences
*/
		sprintf(gcgtarget,">>>>%s ",flatname); l=strlen(gcgtarget);
		sprintf(gcgtarget2,">>>>%s_0",flatname); l2=strlen(gcgtarget2);
		do	{
			fpos = ftello(seq_fich);
			newpnuc = fpos;
/* read the name line */
			p=fgets(gcgline,lnameline,seq_fich);
			if(p==NULL) {
				fprintf(stderr,
"Synchro error .ref/.seq looking for %s\n",gcgtarget);
				exit(ERREUR);
				}
#ifdef vms
/* skip the rest of seq data only under vms+gcg */
			skipseq(gcgline,seq_fich,newpnuc);
#endif
			}
		while( strncmp(gcgtarget,gcgline,l)!=0 && 
				strncmp(gcgtarget2,gcgline,l2)!=0 );
		}
	p=flatname-1; while(*++p) *p=toupper(*p);
	nloc=numloc(flatname);
	if(swissprot) {
/* provisoire tant que bug swissprot avec plusieurs lignes ID */
		do	fgets(gcgline, sizeof(gcgline), flat);
		while(strncmp(gcgline, "ID ", 3) == 0);
		}
	totlu++;
	new=FALSE; modif=FALSE;
	if(nloc==0)
		new=TRUE;
	else	{
		readloc(nloc);
		mere=ploc->sub;
		if( ! testbit(blist, mere) ) { /* seq deja trouvee auparavant */
			if(update_mode)
				printf("Warning: %s re-appears.\n", flatname);
			continue; /* on l'ignore */
			}
		readsub(mere);
		bit0(blist, mere);
		/* comparaison des longueurs */
		if(embl) {
			p=strstr(line,"BP.");
			if(p!=NULL) {
				while(*(--p)==' '); /* recherche mot precedent*/
				while(*(--p)!=' ');
				sscanf(p+1,"%d",&newlong);
				}
			else newlong=0; /* pas de longeur ds ID */
			}
		else if(swissprot) {
			p=strstr(line,"AA.");
			if(p!=NULL) {
				while(*(--p)==' '); /* recherche mot precedent*/
				while(*(--p)!=' ');
				sscanf(p+1,"%d",&newlong);
				}
			else newlong=0; /* pas de longeur ds ID */
			}
		else if(nbrf) {
			newlong=psub->length;   /* longueur non testee */
			}
		if(update_mode)
			modif=(newlong!=psub->length);
		else if(count_scan == 0) { 
			/* pas de comparaison de longueur en mode install */
			if(newlong!=psub->length) printf(
"Warning %.16s has not the same length in index and in division files\n",
						psub->name);
			}
		/* comparaison des dates */
		if(genbank) {
			calc_gbk_date(gbk_date, newdate);
			}
		else if(embl) {
			do	{
				fgets(line,sizeof(line),flat);
				if(strncmp(line,"FH",2) == 0) break;
				if(strncmp(line,"SQ",2) == 0) break;
				}
			while(strncmp(line,"DT",2) != 0);
			if(strncmp(line,"DT",2)!=0) strcpy(newdate,"xxxx");
			else	{
				fgets(line,sizeof(line),flat);
				if(strncmp(line,"DT",2)!=0) 
					strcpy(newdate,"xxxx");
				else calc_gbk_date(line+5,newdate);
				}
			}
		else if(swissprot) {
			do	{
				fgets(line,sizeof(line),flat);
				if(strncmp(line,"FT",2) == 0) break;
				if(strncmp(line,"SQ",2) == 0) break;
				}
			while(strncmp(line,"DT",2) != 0);
			if(strncmp(line,"DT",2)!=0) strcpy(newdate,"xxxx");
			else	{
				fgets(line,sizeof(line),flat);
				fgets(line,sizeof(line),flat);
				if(strncmp(line,"DT",2)!=0) 
					strcpy(newdate,"xxxx");
				else calc_gbk_date(line+5,newdate);
				}
			}
		else if(nbrf) {
			do	{
				fgets(line,sizeof(line),flat);
				if(strncmp(line,"SUMMARY",7)==0) break;
				}
			while(strncmp(line,"DATE",4)!=0);
			if(strncmp(line,"SUMMARY",7)==0) strcpy(newdate,"xxxx");
			else calc_nbrf_date(flat,line,newdate);
			}
		if(!modif && update_mode) /* pas compar date en mode install */
			modif= (strncmp(ploc->date,newdate,8)!=0);
		/* comparaison des tables de features */
		totfilles= -1;
		if ( ( genbank || embl ) && !modif) {
			totfilles=compare_ft(flat,line,sizeof(line),
					ploc->sub,&modif);
			}
		}
	if(!new && ( count_scan > 0 || !update_mode || !modif) ) {
	/* unchanged seq, update INFO + NUCLEOT pointers */
		if(flat_format)compute_pnuc(flat, line, &newpnuc, sizeof(line));
		if( !big_annots ) {
			newpinf=poffset(divnum,newpinf);
			newpnuc=poffset(divnum,newpnuc);
			}
		readloc(nloc);
		ploc->pnuc=newpnuc; ploc->pinf=newpinf;
		if(big_annots)ploc->div = divnum;
		writeloc(nloc);
		subpnuc[ploc->sub] = newpnuc;
		bit0(invalidlist, mere);
		if(totfilles >= 0) {
			update_pinf_filles(divnum, totfilles);
			}
		++(*num_meres);
		*tot_long += (double)newlong;
		}
	if(update_mode) {
		if(new) {
			fprintf(newfile,"%s\n",flatname);
			printf("%s new \n",flatname);
			}
		else if(modif) {
			bit1(suplist,ploc->sub);
			printf("%s mod\n",flatname);
			fprintf(modfile,"%s\n",flatname);
			}
		if(new || modif) fprintf(addrfile,"%u %s\n",newpinf,flatname);
		}
	else if(count_scan == 0) { /* mode install */
		if(new)	{
			printf(
			"Warning: %s in bank but not in acnuc\n",flatname);
			}
		else 	{
		    for(i=0; i<=totfilles; i++) {
			if(tab_fille_data[i].newpinf != 0) {
				++(*num_filles);
				}
			}
		    }
		}
	} /* end of while(TRUE) */
dir_flush(kloc);
dir_flush(kshrt);
if(!flat_format)fclose(seq_fich);
fclose(flat);
printf("%d seqs found in this division\n",totlu); fflush(stdout);
return totlu;
} /* end of proc_division */


void prnumloc(void)
{
int finsort, isub, numloc, rank;

read_first_rec(ksub, &finsort);
subseq_sorted = (finsort == nseq);
if(!subseq_sorted) return;
rank = -1;
for(isub = 2; isub <= nseq; isub++) {
	readsub(isub);
	if(psub->length != 0 && psub->pext <= 0) {
		numloc=psub->plinf;
		rank++;
		memcpy(lismne + rank * L_MNEMO, psub->name, L_MNEMO);
		nom2loc[rank] = numloc;
		}
	}
tot_nom2loc = rank;
}


int numloc(char *mnemo)
{
int val, l, first, last, compar;

if(!subseq_sorted) {
	val=isenum(mnemo);
	if(val) {readsub(val); val=psub->plinf; }
	return val;
	}
first = 0; last = tot_nom2loc;
l=strlen(mnemo);
do	{
	val = (last+first)/2;
	compar = strcmptrail(mnemo, l, lismne + val * L_MNEMO, L_MNEMO);
	if(compar == 0) return nom2loc[val];
	if(compar > 0) first = val+1;
	else	last = val-1;
	}
while(first <= last);
return 0;
}


int compare_ft(FILE *flat, char *line, int lline, int mere, int *changed)
/* recherche ds les features a venir ds flat toutes les filles de mere.
*changed mis a TRUE si differences ds les location sont trouvees
(qualifiers et ordre n'est pas important)
return: nbre de filles de la mere from 0 (toutes testees)
ne detecte pas si nouvelles FT peut creer d'autres filles (cause cas ou filles
normalement non creees)
*/
{
int totfilles, point, i, nfille, l, err, num, trouve, ok=FALSE;
char nom_mere[17], key[20], *p;
static char location[MAXLOCAT+1];
char *oldlocat;
long newpos, current, debut_fi_annot;
off_t fpos;

readsub(mere); memcpy(nom_mere,psub->name,16); nom_mere[16]=0;
totfilles =  -1;
if(psub->pext < 0) { /* parcours de toutes les filles de mere */
	point = - psub->pext;
	while(point) {
		readlng(point); point = plng->next;
		for(i = 0; i < SUBINLNG; i++) {
			nfille = plng->sub[i];
			if(nfille == 0) break;
			readsub(nfille);
			p = strchr(psub->name,'.'); if(p == NULL) continue;
			if(strcmptrail(nom_mere,16,psub->name,p-psub->name)!= 0)
				continue;
			if(totfilles < maxfillespermere ) {
				totfilles++;
				tab_fille_data[totfilles].numero_fille = nfille;
				tab_fille_data[totfilles].newpinf = 0;
				tab_fille_data[totfilles].description = 
					prepextract(nfille, mere);
				}
			else if(update_mode) { /* abandon compar des features */
				fprintf(stderr,
			"Warning: increase parameter maxfillespermere\n");
				goto ret_changed;
				}
			else	{ /* compar limitee au debut des filles */
				point=0; break;
				}
			}
		}
	}
if(totfilles == -1) goto test_final;
/* recherche du debut de la nouvelle table des features */
do	{
	fgets(line,lline,flat);
	if( (genbank && strncmp(line,"ORIGIN",6)==0) || 
		(embl && strncmp(line,"SQ ",3)==0) ) goto test_final;
	}
while( (embl && strncmp(line,"FH",2)!=0) || 
	(genbank && strncmp(line,"FEATURES",8)!=0) );
if(embl) fgets(line,200,flat);
/* parcours de la nouvelle table des features */
do	{
	fpos = ftello(flat);
	newpos = fpos;
	fgets(line,lline,flat); l=strlen(line); line[--l]=0;
   suite:
	if( (embl && strncmp(line,"FT",2) != 0) || 
		(genbank && *line != 0 && strcmptrail(line,4,NULL,0)!=0)) break;
	if(l<=6 || line[5]==' ') continue;
	/* test de la cle de la feature */
	p=line+5; l=0;
	while(*p!=' ' && *p!=0) key[l++]= toupper(*(p++));
	key[l]=0;
	for(i=0; i<=tottypes; i++ )
		if(strcmp(key,list_types[i])==0) break;
	if(i>tottypes) continue;
	strcpy(location,key);
	debut_fi_annot=newpos;
	err = flat_read_location(flat,line,lline,location+l,MAXLOCAT-l,&newpos);
/* formatter la location lue */
	preplocat(location+l, MAXLOCAT-l);
	trouve=FALSE;
	for(num=0; num<=totfilles; num++) {
		if(tab_fille_data[num].newpinf != 0) continue;
		nfille = tab_fille_data[num].numero_fille;
		p = tab_fille_data[num].description;
		if(p != NULL && strcmp(p, location) == 0) {
			trouve=TRUE;
			break;
			}
		}
	if(trouve) tab_fille_data[num].newpinf= (int) debut_fi_annot;
	l=strlen(line);
	goto suite; /* pour travailler la derniere ligne lue par 
			flat_read_location */
	}
while(TRUE);
test_final:
for(num=0; num<=totfilles; num++) {
	if(tab_fille_data[num].description != NULL) 
		free(tab_fille_data[num].description);
	}
for(num=0; num<=totfilles; num++) {
	if(tab_fille_data[num].newpinf == 0) goto ret_changed;
	}
for(num=1; num<=totfilles; num++) { 
/* verifier que les filles sont en ordre croissant dans FT 
ceci pour detecter les changements d'ordre dans les FT 
*/
	if( (unsigned int)tab_fille_data[num].newpinf <= 
		(unsigned int)tab_fille_data[num-1].newpinf) 
			goto ret_changed;
	}
ok = TRUE;
ret_changed:
*changed = !ok;
return totfilles;
} /* end of compare_ft */


void update_pinf_filles(int divnum, int totfilles)
/* mettre a jour les pointeurs INFO des filles memorisees ds tab_fille_data
*/
{
int num, numfille, point, newpinf;

for(num=0; num<=totfilles; num++) {
	numfille=tab_fille_data[num].numero_fille;
	newpinf=tab_fille_data[num].newpinf;
	if(newpinf==0) continue;
	bit0(invalidlist, numfille);
	readsub(numfille);
	point=psub->plinf;
	readshrt(point);
	if( !big_annots ) newpinf=poffset(divnum,newpinf);
	pshrt->val=newpinf;
	writeshrt(point);
	}
}


static int flat_read_location(FILE *flat, char *line, int lline, 
	char *location, int maxloclen, long *newpos)
/* 
charge location avec la location ds fichier flat 
lue dans line transmise puis ensuite a partir d'adresse courante
*newpos mis a jour a chaque nouvelle lecture dans flat
rend la longueur chargee, ou 0 si location + longue que maxloclen 
*/
{
int i, lloc=0, l, nlu=0;
char *ps;
off_t fpos;

while((ps=strchr(line,'/'))==NULL) {
	i=strlen(line)-21;
	if(lloc+i>=maxloclen) return 0;
	if(i>0) {
		memcpy(location+lloc,line+21,i);
		lloc+=i;
		}
	fpos = ftello(flat);
	*newpos = fpos;
	fgets(line,lline,flat); l=strlen(line); line[l-1]=0; nlu++;
	if(strcmptrail(line,20,NULL,0)!=0 && 
		strcmptrail(line,20,"FT",2) !=0 ) break;
	}
if(ps!=NULL) {
	i=ps-line+1;
	if(strcmptrail(line+2,i-3,NULL,0)!=0 && i>22) {
		if(lloc+i-22>=maxloclen) return 0;
		memcpy(location+lloc,line+21,i-22);
		lloc+=(i-22);
		}
	}
location[lloc]=0;
compact(location);
ps=location; while(*ps!=0) { *ps=toupper(*ps); ps++; }
if(nlu==0) {
	fpos = ftello(flat);
	*newpos = fpos;
	fgets(line,lline,flat); l=strlen(line); line[l-1]=0;
	}
return strlen(location);
}


#define MAXELEMENTS 400
static char *elements[MAXELEMENTS];
int tri_chaines(const void *s1, const void *s2);

char *prepextract(int seqnum, int mere)
{
int point, l, ltot, totelements, i;
char *p, *result;
static char buffer[50];

readsub(seqnum);
point = psub->pext;
ltot = 0; totelements = -1;
while(point != 0) {
	readext(point); point = pext->next;
	l = 0;
	if(pext->mere != mere) {
		readsub(pext->mere);
		memcpy(buffer, psub->name, L_MNEMO);
		buffer[L_MNEMO] = 0;
		compact(buffer);
		strcat(buffer, ":");
		l = strlen(buffer);
		}
	if(pext->deb != pext->fin) {
		int a,b;
		if(pext->deb <= pext->fin) { a = pext->deb; b = pext->fin; }
		else { a = pext->fin; b = pext->deb; }
		sprintf(buffer + l, "%d..%d", a, b);
		}
	else
		sprintf(buffer + l, "%d", pext->deb);
	/* memoriser l'element */
	l = strlen(buffer);
	ltot += l;
	if(++totelements >= MAXELEMENTS ) {
		fprintf(stderr,"increase MAXELEMENTS\n");
		exit(ERREUR);
		}
	elements[totelements] = (char *)malloc(l+1);
	if(elements[totelements] == NULL) return NULL;
	memcpy(elements[totelements], buffer, l + 1);
	}
/* trier les elements */
if(totelements>0) {
	qsort(elements,totelements+1,sizeof(char *),tri_chaines);
	}
/* en tete le nom du type de la fille */
readsub(seqnum);
readsmj(psub->type);
memcpy(buffer, psmj->name + 2, 18); buffer[18] = 0; compact(buffer);
l = strlen(buffer);
/* fabriquer chaine resultat */
result = (char *)malloc(l + ltot + totelements + 1);
if(result == NULL) return NULL;
memcpy(result, buffer, l);
p = result + l;
for(i = 0; i <= totelements; i++) {
	l = strlen(elements[i]);
	memcpy(p, elements[i], l);
	p += l;
	*(p++) = ',';
	}
*(p-1) = 0;
/* liberer la memoire allouee ici */
for(i = 0; i <= totelements; i++)  free(elements[i]);
return result;
}


void preplocat(char *location, int maxlen)
/* 
formatter la location (sur place) en enlevant les operateurs et les < et >
et en triant par ordre alphab les a..b ou X1:a..b restants
*/
{
int lelement, totelements, i, llocat, newllocat;
char *p, *q, *debut_element, *virg, *fin_element;

newllocat=0;
totelements= -1;
/* supprimer les < et > */
p=location-1; while(*(++p)!=0) if(*p=='<' || *p=='>') *p=' ';
compact(location);
llocat=strlen(location);
debut_element=location;

suite:
/* isoler l'element en cherchant les , et sautant les () */
virg=debut_element;
while(*++virg!=0) {
	if(*virg==',')break;
	if(*virg=='('){
		virg=nextpar(virg,location+llocat-1);
		if(virg==NULL) { /* no matching ) */
			*location=0;
			return;
			}
		}
	}
fin_element=virg-1;
/* eliminer les () externes de l'element */
p=strchr(debut_element,'(');
if(p!=NULL && p<=fin_element) {
	q=nextpar(p,fin_element);
	if(q==NULL) { /* no matching ) */
		*location=0;
		return;
		}
	if(q==fin_element) {
		char *finp;
		debut_element=p+1; 
		finp = location+llocat;
		for(p=q; p<finp ; p++) *p = *(p+1);
		llocat--;
		goto suite;
		}
	}

/* remplacer pos..pos par pos */
p=strstr(debut_element,"..");
if(p!=NULL && p<fin_element) {
	q=strchr(debut_element,':');
	if(q==NULL || q>=p) q=debut_element-1;
	q++;
	if( (p-q == fin_element-p-1) && strncmp(q,p+2,p-q)==0) fin_element=p-1;
	}
/* memoriser l'element */
lelement=fin_element-debut_element+1;
newllocat+= (lelement+1);
if(++totelements>=MAXELEMENTS || newllocat>MAXLOCAT) {
	fprintf(stderr,"Location too long for preplocat: %s\n",location);
	exit(ERREUR);
	}
elements[totelements] = (char *)malloc(lelement+1);
memcpy(elements[totelements], debut_element, lelement);
elements[totelements][lelement]=0;
/* remplacer acc num par nom-mere:  AC000001.1:a..b ==> ECOTGP:a..b  */
p = strchr(elements[totelements], ':');
if(p != NULL ) {
	char access[40], *nmere, *tmp;
	i = p - elements[totelements];
	q = strchr(elements[totelements], '.');
	if(q != NULL && q < p) i = q - elements[totelements];
	memcpy(access, elements[totelements], i); access[i] = 0;
	nmere = acctomere(access);
	if(nmere != NULL) {
		tmp = (char *)malloc(strlen(nmere) + strlen(p) + 1);
		strcpy(tmp, nmere); strcat(tmp, p);
		free(elements[totelements]);
		elements[totelements] = tmp;
		}
	}
/* aller a l'element suivant */
debut_element=virg+1;
if(debut_element<location+llocat)goto suite;
/* trier les elements */
if(totelements>0) {
	qsort(elements,totelements+1,sizeof(char *),tri_chaines);
	}
/* refabriquer location comme suite d'elements */
p=location;
for(i=0; i<=totelements; i++) {
	lelement=strlen(elements[i]);
	memcpy(p,elements[i],lelement);
	p+=lelement;
	*(p++)=',';
	}
*(p-1)=0;
/* liberer la memoire allouee ici */
for(i=0; i<=totelements; i++)  free(elements[i]);
return;
#undef MAXELEMENTS
}
int tri_chaines(const void *s1, const void *s2)
{
return strcmp( *(char **)s1, *(char **)s2);
}


static char *acctomere(char *access)
{
int nacc, point, mere, next;
static char nom_mere[L_MNEMO+1];

nacc = fastfindacc(access);
if(nacc == 0) return NULL;
readacc(nacc);
point = pacc->plsub;
if(point == 0) return NULL;
readshrt(point);  
if(pshrt->next == 0) { /* une seule seq associee a access */
	mere = pshrt->val;
	}
else	{ /* plusieurs seqs associees, chercher celle dont primary */
	while(TRUE) {
		if(point == 0) return NULL;
		readshrt(point); next = pshrt->next; 
		mere = pshrt->val;
		readsub(mere); readloc(psub->plinf);
		readshrt(ploc->placc);
		if(pshrt->val == nacc) break;
		point = next;
		}
	}
readsub(mere);
memcpy(nom_mere, psub->name, L_MNEMO); nom_mere[L_MNEMO] = 0; compact(nom_mere);
return nom_mere;
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
	if(*debut==')') {
		return debut;
		}
	else debut++;
	}
while(debut<=fin);
return NULL;
}


static char mois[]="JANFEBMARAPRMAYJUNJULAUGSEPOCTNOVDEC";

void calc_gbk_date(char *gbk, char *newdate)
/* input gbk pointe sur 21-JUL-1994
retour: newdate rempli avec 07/21/94
*/
{
int jour, num_mois, an;
char *p;

sscanf(gbk,"%d",&jour);
sscanf(gbk+9,"%d",&an);
gbk[6]=0;
p=strstr(mois,gbk+3);
if(p!=NULL)  num_mois=(p-mois)/3+1;
else 	num_mois=0;
sprintf(newdate,"%2.2d/%2.2d/%2.2d",num_mois,jour,an);
return;
}


void calc_nbrf_date(FILE *flat, char *line, char *newdate)
{
char *p, *q;
int jour, num_mois, an, curr_date, max_date;

max_date=0;
do	{ /* boucle sur ligne line et sur les suivantes lues ds flat */
	/* remove DATE */
	memset(line,' ',4);
	/* remove all ; */
	while( (p=strchr(line,';')) !=NULL ) *p=' ';
	/* remove all \n */
	while( (p=strchr(line,'\n')) !=NULL ) *p=' ';
	/* remove all #qualif */
	while( (p=strchr(line,'#')) !=NULL ) {
		while(*p!=' ') {*p=' '; p++; }
		}
	/* to upper */
	p=line-1; while(*(++p)!=0) *p=toupper(*p);
	/* boucle sur tous les mots qui restent */
	p=line;
	while(TRUE) {
		while(*p==' ') p++;
		if(*p==0) break;
		sscanf(p,"%d",&jour);
		sscanf(p+7,"%d",&an);
		p[6]=0;
		q=strstr(mois,p+3);
		if(q!=NULL)  num_mois=(q-mois)/3+1;
		else 	num_mois=1;
		curr_date=an*10000 + num_mois*100 + jour;
		if(curr_date>max_date) max_date=curr_date;
		p += 11;
		}
	fgets(line,80,flat);
	}
while(*line==' ');
jour= max_date % 100;
num_mois= max_date/100;
num_mois= num_mois % 100;
an= max_date/10000;
an = an % 100;
sprintf(newdate,"%2.2d/%2.2d/%2.2d",num_mois,jour,an);
}


char *fulldivname(char *divname, char *prefix)
{
char *fname;

fname= (char *)malloc(strlen(prefix)+strlen(divname)+10);
strcpy(fname,prefix);
strcat(fname,divname);
if(!flat_format) {
	strcat(fname,".ref");
	}
else if(embl) {
	strcat(fname,".dat");
	}
else if(genbank) {
	strcat(fname,".seq");
	}
else if(nbrf) {
	strcat(fname,".dat");
	}
else if(swissprot) {
	strcat(fname,".dat");
	}
return fname;
}


void compute_pnuc(FILE *flat, char *line, int *newpnuc, size_t lline)
{
/* unchanged seq, update INFO + NUCLEOT pointers */
off_t fpos;

if(genbank) {
	while(strncmp(line,"ORIGIN",6)!=0)
		fgets(line,lline,flat);
	}
else if(embl || swissprot) {
	while(strncmp(line,"SQ ",3)!=0)
		fgets(line,lline,flat);
	}
else if(nbrf) {
	while(strncmp(line,"SEQUENCE",8)!=0)
		fgets(line,lline,flat);
/* sauter la ligne de numerotation des sequences */
	fgets(line,lline,flat);
	}
fpos = ftello(flat);
*newpnuc = fpos;
/* sauter la suite */
do	fgets(line, lline, flat);
while(strncmp(line, "//", 2) != 0);
}



int fastfindacc(char *access)
{
static int first = TRUE, tmem, sorted, size_n_r_struct;
int tacc, num, premier, dernier, ord, rank;
static char *p_access;
static char *new_accs;
char *p;

if(first) {
	first = FALSE;
	p_access = (char *)malloc(ACC_LENGTH + 1);
	size_n_r_struct = sizeof(int) + ACC_LENGTH;
	tacc = read_first_rec(kacc, &sorted);
	if(sorted == 0) sorted = 1;
	tmem = tacc - sorted;
	if(tmem > 0) {
		new_accs = (char *)malloc(tmem * size_n_r_struct);
		p = new_accs;
		for(num = sorted + 1; num <= tacc; num++, p += size_n_r_struct) {
			readacc(num);
			memcpy(p, pacc->name, ACC_LENGTH);
			memcpy(p + ACC_LENGTH, &num, sizeof(int));
			}
		qsort(new_accs, tmem, size_n_r_struct, n_r_compar);
		}
	}
padtosize(p_access, access, ACC_LENGTH);
premier = 2; dernier = sorted;
while(premier <= dernier) {
	num = (dernier + premier) / 2;
	readacc(num);
	ord = memcmp(p_access, pacc->name, ACC_LENGTH);
	if(ord == 0) 
		return num;
	else if(ord > 0) premier = num + 1;
	else dernier = num - 1;
	}

premier = 0; dernier = tmem - 1;
while(premier <= dernier) {
	num = (dernier + premier) / 2;
	ord = memcmp(p_access, new_accs + num * size_n_r_struct, ACC_LENGTH);
	if(ord == 0) {
		memcpy(&rank, new_accs + num * size_n_r_struct + ACC_LENGTH, sizeof(int) );
		return rank;
		}
	else if(ord > 0) premier = num + 1;
	else dernier = num - 1;
	}

return 0;
}


static int n_r_compar(const void *p1, const void *p2)
{
return memcmp(p1, p2, ACC_LENGTH);
}



void proc_help_file(char *name, int totmeres, double totlong, 
	int totfilles, int totbib, char *acnuc_env)
{
char fullname[100], ligne[100], fullname_new[100];
FILE *fichin, *fichout;
int i, ok=FALSE;

strcpy(fullname,acnuc_env);
strcat(fullname,name);
fichin = fopen(fullname,
#ifdef __INTEL__
	"rb"
#else
	"r"
#endif
	);
if(fichin == NULL) return;
strcpy(fullname_new,fullname);
strcat(fullname_new,".NEW");
fichout = fopen(fullname_new,"w");
if(fichout == NULL) {
	fclose(fichin);
	return;
	}
do	{
	for(i=1; i<=3; i++) {
		fgets(ligne,100,fichin);
		if( fputs(ligne,fichout) == EOF) break;
		}
	fgets(ligne,100,fichin);
	
	if(embl || genbank)
		sprintf(ligne,
		"%s bases; %s sequences; %s subsequences; %s references.\n",
		paqde3(totlong), paqde3(totmeres), paqde3(totfilles), 
		paqde3(totbib) );
	else
		sprintf(ligne,
		"%s aminoacids; %s sequences; %s references.\n",
		paqde3(totlong), paqde3(totmeres), paqde3(totbib) );
	if( fputs(ligne,fichout) == EOF) break;
	
	while( fgets(ligne,100,fichin)!=NULL ) {
		if( fputs(ligne,fichout) == EOF) break;
		}
	ok=TRUE;
	}
while(FALSE);
fclose(fichin); fclose(fichout);
if(ok)	{
	remove(fullname);
	rename(fullname_new,fullname);
	}
else remove(fullname_new);
}


char *paqde3(double val)
{
char *result, aux[40], *p, *q;
int l, n_virg, offset;

result = (char *)malloc(40);
sprintf(aux, "%.0f", val);
l = strlen(aux) ;
if(l <= 3) {
	strcpy(result, aux);
	return result;
	}
n_virg = (l - 1)/3;
offset = l % 3; if(offset == 0) offset = 3;
memcpy(result, aux, offset);
*(result + offset) = ',';
p = result + offset + 1;
q = aux + offset;
while( 1 ) {
	memcpy(p, q, 3);
	q += 3;
	p += 3;
	if(q >= aux + l) break;
	*p = ',';
	p++;
	}
*p = 0;
return result;
}


void get_acnuc_dirs(char *acnuc_env, char *gcgacnuc_env)
{
int erreur = FALSE;
char nom[20], *p;

#ifdef vms
#define ACNUC "BANK"
#define GCGACNUC "GCGBANK"
#else
#define ACNUC "acnuc"
#define GCGACNUC "gcgacnuc"
#endif

p = prepare_env_var(ACNUC);
if(p == NULL) {
	erreur = TRUE;
	strcpy(nom, ACNUC);
	}
else strcpy(acnuc_env, p);
if( !erreur ) {
	p = prepare_env_var(GCGACNUC);
	if(p == NULL) {
		erreur = TRUE;
		strcpy(nom, GCGACNUC);
		}
	else strcpy(gcgacnuc_env, p);
	}
if(erreur) {
	fprintf(stderr,"Error: environment variable %s is not defined\n", nom);
	exit(ERREUR);
	}
#undef ACNUC
#undef GCGACNUC
}


#ifdef vms
void skipseq(char *nameline, FILE *fich, int cur_pos)
{
char *p;
int longueur, binary, truelong, numfull, numreste, toskip;
off_t fpos;

p=strstr(nameline," Len:")+5;
sscanf(p,"%d",&longueur);
binary = (strstr(nameline," 2BIT ")!=NULL);
/* compute true seq data length */
if(binary)
	truelong= (longueur+3)/4;
else
	truelong=longueur;
/* compute # of complete seq records */
numfull = truelong/MAXRECLENGTH;
/* compute length of last record */
numreste = truelong % MAXRECLENGTH;
if(numreste%2==1) ++numreste; /* odd rec lengths are padded by one byte */
/* compute byte offset to skip all of seq data */
toskip = numfull*(MAXRECLENGTH+2);
if(numreste!=0) toskip += (numreste+2);
/* add header & title lines to skipped data */
toskip += (MINRECSIZE + MAXDOCSTRING+4);
/* skip all of seq data */
fpos = cur_pos+toskip;
fseeko(fich, fpos, SEEK_SET);
}
#endif


#ifdef unix
int proc_options(int argc, char **argv)
{
int num, count = 0, i;

num = 0;
while (++num < argc) { /* options -mmap xxx  */
	if(strcmp(argv[num], "-mmap") != 0) continue;
	num++;
	if(strcmp(argv[num], "kshrt") == 0) dir_set_mmap(kshrt);
	else if(strcmp(argv[num], "klng") == 0) dir_set_mmap(klng);
	else if(strcmp(argv[num], "ksub") == 0) dir_set_mmap(ksub);
	else if(strcmp(argv[num], "kloc") == 0) dir_set_mmap(kloc);
	}
num = 0; 
while (++num < argc) { /* option -<number> xxx yyy zzz  */
	if(argv[num][0] != '-') continue;
	if(strcmp(argv[num], "-mmap") == 0) continue;
	count = 0;
	sscanf(argv[num] + 1, "%d", &count);
	if(count <= 0 || count > 100) { count = 0; continue; }
	newdivname = (char **)malloc(count * sizeof(char *));
	for(i = 0; i < count; i++) {
		num++;
		newdivname[i] = (char *)malloc(strlen(argv[num]) + 1);
		strcpy(newdivname[i], argv[num]);
		}
	}
return count;
}
#endif

