#include "dir_acnuc.h"
#define Boolean char

extern void *mycalloc(int nbre, int size);
extern void seq_to_annots(int numseq, long *faddr, int *div);
extern char *read_annots(long faddr, int div);
extern char *next_annots(long *pfaddr);
extern char init_codon_to_aa(char *codon, int gc);


void pretty_seq(int numseq, int basesperline, Boolean translate, 
			void (*outline)(char *line));
static int is_5_partial(int numseq);


static void cplsq(char *seq, int longueur)
{
static char bases[]="aAcCgGtTuUrRyYnN";
static char compl[]="tTgGcCaAaAyYrRnN";
char *pos, *fin;
fin=seq+longueur;
while(seq<fin) {
	pos=strchr(bases,*seq);
	if(pos==NULL) pos=strchr(bases,'n');
	*seq=compl[pos-bases];
	seq++;
	}
}

typedef struct endpoint {
	int numseq;
	int position;
	int phase;
	int initiateur;
	struct endpoint *next;
	} endpoint;

static Boolean newendpoint(endpoint **debut_endpoints, int numseq, int phase, int extremite, 
	Boolean traduc, Boolean is_init)
{
/* pour chainer des struct endpoint par ordre croissant de abs(position) */
endpoint *pnew, *courant, *previous;
pnew=(endpoint *)mycalloc(1,sizeof(endpoint));
if(pnew==NULL) return TRUE;
courant= *debut_endpoints; previous=NULL;
while(courant != NULL) {
	if(abs(courant->position)>abs(extremite) ) break;
	previous=courant;
	courant=courant->next;
	}
if(previous!=NULL) {
	pnew->next=previous->next;
	previous->next=pnew;
	}
else	{
	pnew->next= *debut_endpoints;
	*debut_endpoints=pnew;
	}
pnew->numseq=numseq;
pnew->position=extremite;
if(traduc)
	pnew->phase=(phase%3)+1;
else
	pnew->phase=0;
pnew->initiateur = is_init;
return FALSE;
}


void pretty_seq(int numseq, int basesperline, Boolean translate, 
			void (*outline)(char *line))
{
endpoint *debut_endpoints=NULL, *current_endpoint, *first_endpoint_in_line, 
	*one_endpoint, *endpoint_in_next_line;
char ligne[3][150], transl[3][150], aux[150], oneline[150], seq[150],
	*carac;
static Boolean first=TRUE;
static int cdstype; 
static char codon[4];
int point, i, num, l, phase, j, codgen, pointext, numfille, phase_endpoint, 
	lseq, position, absposition, p, lm, f, larg, k, fin, flag, 
	pos, longueur, bpl2, np, is_init;
int trad[3], dtrad[3], ftrad[3], der[3], tradc[3], dtradc[3];
Boolean decoup[3], ntrad[3];
Boolean traduc, warn, init_in_lignes;

if(first) {
	first=FALSE;
	if(!nbrf) cdstype=fcode(ksmj,"04CDS",20);
	codon[3]=0;
	}

codgen = 0;
lseq=basesperline; bpl2=basesperline+2; 
larg=basesperline+20+basesperline/10-1; if(larg<80) larg=80;
readsub(numseq); longueur=psub->length;
sprintf(aux,"Name: %.*s Length:%d", L_MNEMO, psub->name,longueur);
(*outline)(aux);
if(translate && psub->type==cdstype) {
	phase=psub->phase%100; 
	is_init = (phase == 0 && !is_5_partial(numseq));
	codgen=psub->phase/100;
	newendpoint(&debut_endpoints,numseq,phase+1,1,TRUE,is_init);
	newendpoint(&debut_endpoints,numseq,phase+1,-longueur,TRUE,FALSE);
	}
/* parcourir la liste des filles et memoriser les points de decoupage */
if(psub->pext<0) (*outline)(
			"Subsequence names, lengths, types and qualifiers:");
point= -psub->pext;
while(point>0) {
	long faddr; int div, compl;
	readlng(point); point=plng->next;
	for(num=0; num<SUBINLNG; num++) {
		if(plng->sub[num]==0) break;
		numfille=plng->sub[num];
		warn=FALSE;
		readsub(numfille);
		memset(aux,' ',80);
		l=2;
		memcpy(aux+l,psub->name,16); l+=16;
		sprintf(aux+l,"%7d   ",psub->length); l+=10;
		seq_to_annots(numfille, &faddr, &div);
		if( read_annots(faddr, div) != NULL) {
			memcpy(aux+l,pinfo->line+5,15); l+=15;
			while(aux[l-1]==' ') l--;
			while(l<80) {
				next_annots(NULL);
				if(strcmptrail(pinfo->line,20,NULL,0) &&
				     strcmptrail(pinfo->line,20,"FT",2) ) break;
				if( (carac=strchr(pinfo->line,'/')) == NULL) 
					continue;
				k = strlen(pinfo->line) - (carac - pinfo->line);
				if(l+k>80) k = 80 - l;
				memcpy(aux+l,pinfo->line+21,k); l += k;
				while(aux[l-1] == ' ') l--;
				}
			aux[l]=0;  (*outline)(aux);
			}
		traduc=FALSE;
		if(translate) {
			traduc=(psub->type==cdstype);
			if(traduc) {
				phase=psub->phase%100; codgen=psub->phase/100;
				}
			}
		pointext=psub->pext;
		is_init = (traduc && phase == 0 && !is_5_partial(numfille));
		while(pointext) {
			readext(pointext); pointext=pext->next;
			if(traduc) {
				if(pext->deb<=pext->fin)
					phase_endpoint=pext->deb+phase;
				else
					phase_endpoint=pext->deb-phase;
				phase=(abs(pext->deb - pext->fin)+1-phase)%3;
				if(phase==1)phase=2;
				else if(phase==2)phase=1;
				}
			if(pext->mere!=numseq) warn=TRUE;
			else 	{
				if(pext->deb>pext->fin) {
					i=pext->deb; pext->deb=pext->fin; pext->fin=i;
					numfille= -abs(numfille);
					compl = TRUE;
					}
				else 	compl = FALSE;
				newendpoint(&debut_endpoints,numfille,
					phase_endpoint,pext->deb,traduc, (compl ? FALSE : is_init) );
				newendpoint(&debut_endpoints,numfille,
					phase_endpoint,-pext->fin,traduc, (compl ? is_init : FALSE) );
				}
			is_init = FALSE;
			}
		if(warn) {
			readsub(abs(numfille));
			sprintf(aux,
			"Subsequence: %.*s has also another parent sequence.",
				L_MNEMO, psub->name);
			(*outline)(aux);
			(*outline)("It will be partial only.");
			}
		}
	}
if(codgen != 0) {
	sprintf(aux, "Genetic code used: %s", get_code_descr(codgen));
	(*outline)(aux);
	}
(*outline)("");

/* edition de la seq et de son decoupage */
memset(trad,0,3*sizeof(int));
memset(tradc,0,3*sizeof(int));
memset(dtrad,0,3*sizeof(int));
memset(dtradc,0,3*sizeof(int));
l=0; current_endpoint=debut_endpoints; flag=0;
next_seq_piece:
i=gfrag(numseq,l+1,bpl2,seq);
if(i == 0) return;
if(i<lseq) lseq=i;
if(i<bpl2) memset(seq+i,' ',bpl2-i);
strcpy(aux,seq);
for(i=0; i<3; i++) {
	memset(ligne[i],' ',130);
	decoup[i]=FALSE;
	der[i]=0;
	memset(transl[i],' ',110);
	ntrad[i]=FALSE;
	ftrad[i]=0;
	}
init_in_lignes = FALSE;
pos=0; np=0; first_endpoint_in_line=current_endpoint;
while(current_endpoint!=NULL && abs(current_endpoint->position)<=l+lseq) {
	position=current_endpoint->position;
	absposition=abs(position);
	traduc=FALSE;
	numfille=current_endpoint->numseq;
	phase=current_endpoint->phase;
	if(numfille>0) {
		p=absposition-l+(absposition-l-1)/10+10;
		readsub(numfille);
		if(phase!=0) {
			traduc=TRUE; phase--;
			}
		lm=16; while(psub->name[lm-1]==' ') lm--;
		i=1;
		if(position>0) {
			f=p+lm; if(f>larg) f=larg;
			if(f==larg) i += (p+lm-larg);
			k=0;
			while(p<=der[k] && k<2) k++;
			if(p<=der[k]) k=0;
			if(current_endpoint->initiateur) {
				ligne[k][p-1] = '!'; /* presence d'un codon initiateur */
				init_in_lignes = TRUE;
				}
			else ligne[k][p-1]='>';
			fin=i;
			memcpy(ligne[k]+p,psub->name+fin-1,f-p);
			fin += (f-p);
			der[k]=f; decoup[k]=TRUE;
			if(traduc) {
				trad[phase]++;
				if(dtrad[phase]==0)
					dtrad[phase]=absposition+(phase-absposition+
						3*(absposition/3+1))%3;
				}
			}
		else	{
			f=p-lm; if(f<1) f=1;
			if(f==1) i += (1-p+lm);
			k=0;
			while(f<=der[k] && k<2) k++;
			if(f<=der[k]) k=0;
			fin=i;
			memcpy(ligne[k]+f-1,psub->name+fin-1,p-f);
			fin += (p-f);
			ligne[k][p-1]='<';
			der[k]=p; decoup[k]=TRUE;
			if(traduc && trad[phase]!=0) {
				if(dtrad[phase]<=absposition-2) {
					k=0;
					if(dtrad[phase]<=ftrad[k]) {
						k++;
						if(dtrad[phase]<=ftrad[k]) k++;
						}
					ntrad[k]=TRUE;
					i=dtrad[phase];
					while(i<=absposition-2) {
						j=i;
						transl[k][i-l] = codaa(aux+i-l-1, codgen);
						if(init_in_lignes) { /* traitement des codons initiateurs */	
							int tmp, kk; 
							tmp = i-l+(i-l-1)/10 + 10;
							for(kk = 0; kk < 3; kk++) 
								if(ligne[kk][tmp-1] == '!') {
									transl[k][i-l] = init_codon_to_aa(aux+i-l-1, codgen);
									ligne[kk][tmp-1] = '>';
									}
							}
						i+=3;
						}
			/* codon entre exons sur brin direct */
					if(j+2<absposition) {
						f=absposition-(j+2); traduc=FALSE;
						one_endpoint=current_endpoint;
						if(one_endpoint->next!=NULL) {
							do	{
								one_endpoint=one_endpoint->next;
								traduc=(abs(one_endpoint->numseq)==numfille);
								}
							while(!traduc && one_endpoint->next!=NULL);
							}
						if(traduc) {
							memcpy(codon,aux+2+j-l,f);
							p=3-f;
							gfrag(numseq,abs(one_endpoint->position),p,codon+f);
							j+=3;
							transl[k][j-l] = codaa(codon,codgen);
							}
						}
					dtrad[phase]=j+3; ftrad[k]=j+2;
					}
				trad[phase]--;
				if(trad[phase]==0) dtrad[phase]=0;
				}
			}
		}
	else	{
		np++;
		if(np==1) cplsq(seq,lseq+2);
		if(position>0) {
			flag++;
			if(flag==1 && position-l-1-pos>0) memset(seq+pos,' ',position-l-1-pos);
			}
		else	{
			flag--;
			if(flag==0) pos=absposition-l;
			}
		}
	current_endpoint=current_endpoint->next;
	}
for(phase=0; phase<=2; phase++) {
	if(trad[phase]!=0 && dtrad[phase]<=l+lseq) {
		k=0;
		if(dtrad[phase]<=ftrad[k]) {
			k++;
			if(dtrad[phase]<=ftrad[k]) k++;
			}
		ntrad[k]=TRUE;
		i=dtrad[phase];
		while(i<=l+lseq) {
			j=i;
			transl[k][i-l] = codaa(aux+i-l-1,codgen); /* sauf dernier bloc */
			if(init_in_lignes) { /* traitement des codons initiateurs */	
				int tmp, kk; 
				tmp = i-l+(i-l-1)/10 + 10;
				for(kk = 0; kk < 3; kk++) 
					if(ligne[kk][tmp-1] == '!') {
						transl[k][i-l] = init_codon_to_aa(aux+i-l-1, codgen);
						ligne[kk][tmp-1] = '>';
						}
				}
			i+=3;
			}
		dtrad[phase]=j+3; ftrad[k]=j+2;
		}
	}
/* toutes les editions du brin direct */
memset(oneline,' ',130); j=10;
for(i=1; i<=1+(lseq-1)/10; i++) {
	sprintf(oneline+j,"%10d ",l+i*10);
	j+=11;
	}
(*outline)(oneline); /*ligne de numerotation */
for(j=2; j>=0; j--) {
	if(ntrad[j]) {
		memset(oneline,' ',10); k=10;
		i=1;
		while(i<=lseq+1) {
			sprintf(oneline+k,"%10.10s ",transl[j]+i-1);
			k+=11; i+=10;
			}
		while(oneline[k-1]==' ') k--;
		oneline[k]=0; (*outline)(oneline); /* ligne de traduction */
		}
	}
if( (lseq%10) !=0) memset(aux+lseq,' ',102-lseq);
memset(oneline,' ',10); k=10;
i=1;
while(i<=lseq) {
	sprintf(oneline+k,"%10.10s ",aux+i-1); k+=11; i+=10;
	}
(*outline)(oneline); /* ligne de sequence */
for(i=0; i<=2; i++) {
	if(decoup[i]) {
		k=larg; while(k>0 && ligne[i][k-1]==' ') k--;
		ligne[i][k]=0; 
		(*outline)(ligne[i]); /* ligne de decoupage */
		}
	}
endpoint_in_next_line=current_endpoint;
if( flag==0 && np==0 )goto fin_seq_frag;

/* traitement des decoupages du brin complementaire */
current_endpoint=first_endpoint_in_line;
for(i=0; i<3; i++) {
	memset(ligne[i],' ',130);
	decoup[i]=FALSE;
	der[i]=0;
	memset(transl[i],' ',110);
	ntrad[i]=FALSE;
	ftrad[i]=0;
	}
init_in_lignes = FALSE;
if(np==0) cplsq(seq,lseq+2);
do	{
	if( np==0 || endpoint_in_next_line==first_endpoint_in_line ) break;
	if(current_endpoint->numseq > 0) continue;
	numfille=abs(current_endpoint->numseq);
	phase=abs(current_endpoint->phase);
	position=current_endpoint->position; absposition=abs(position);
	p=absposition-l+(absposition-l-1)/10+10;
	readsub(numfille);
	if(phase!=0) {
		traduc=TRUE; phase--;
		}
	else traduc=FALSE;
	lm=16; while(psub->name[lm-1]==' ') lm--;
	i=1;
	if(position>0) {
		f=p+lm; if(f>larg) f=larg;
		if(f==larg) i += (p+lm-larg);
		k=0;
		while(p<=der[k] && k<2) k++;
		if(p<=der[k]) k=0;	
		ligne[k][p-1]='>';
		fin=i;
		memcpy(ligne[k]+p,psub->name+fin-1,f-p);
		fin += (f-p);
		der[k]=f; decoup[k]=TRUE;
		if(traduc) {
			tradc[phase]++;
			if(dtradc[phase]==0)
				dtradc[phase]=absposition+(phase-absposition+
					3*(absposition/3+1))%3;
			}
		}
	else	{
		f=p-lm; if(f<1) f=1;
		if(f==1) i += (1-p+lm);
		k=0;
		while(f<=der[k] && k<2) k++;
		if(f<=der[k]) k=0;
		fin=i;
		memcpy(ligne[k]+f-1,psub->name+fin-1,p-f);
		fin += (p-f);
		if(current_endpoint->initiateur) {
			ligne[k][p-1] = '!';
			init_in_lignes = TRUE;
			}
		else 	ligne[k][p-1]='<';
		der[k]=p; decoup[k]=TRUE;
		if(traduc && tradc[phase]!=0) {
			if(dtradc[phase]<=absposition) {
				k=0;
				if(dtradc[phase]<=ftrad[k]) {
					k++;
					if(dtradc[phase]<=ftrad[k]) k++;
					}
				ntrad[k]=TRUE;
				i=dtradc[phase];
				while(i<=absposition) {
					j=i;
					strcpy(codon,"NNN");
					codon[0]=seq[i-l-1];
					if(i-l-2>=0) codon[1]=seq[i-l-2];
					if(i-l-3>=0) codon[2]=seq[i-l-3];
					if(i-l-3>=0) {
						transl[k][i-l-2] = codaa(codon,codgen);
						if(init_in_lignes) { /* traitement des codons initiateurs */	
							int tmp, kk; 
							tmp = i-l+(i-l-1)/10 + 10 - 1;
							for(kk = 0; kk < 3; kk++) 
								if(ligne[kk][tmp] == '!') {
									transl[k][i-l-2] = init_codon_to_aa(codon, codgen);
									ligne[kk][tmp] = '<';
									}
							}
						}
					i+=3;
					}
		/* codon entre exons sur brin complementaire */
				if(j<absposition) {
					f=absposition-j; traduc=FALSE;
					one_endpoint=current_endpoint;
					if(one_endpoint->next!=NULL) {
						do	{
							one_endpoint=one_endpoint->next;
							traduc=(abs(one_endpoint->numseq)==numfille);
							}
						while(!traduc && one_endpoint->next!=NULL);
						}
					if(traduc) {
						memcpy(aux,seq+j-l,f);
						p=3-f;
						gfrag(numseq,abs(one_endpoint->position),p,aux+f);
						cplsq(aux+f,p);
						aux[3]=aux[2]; aux[4]=aux[1]; aux[5]=aux[0];
						j+=3;
						transl[k][j-l-2] = codaa(aux+3,codgen);
						}
					}
				dtradc[phase]=j+3; ftrad[k]=j;
				}
			tradc[phase]--;
			if(tradc[phase]==0) dtradc[phase]=0;
			}
		}
	}
while( (current_endpoint=current_endpoint->next) != endpoint_in_next_line);
if(np!=0 && lseq>pos && flag==0) memset(seq+pos,' ',lseq-pos);
for(phase=0; phase<=2; phase++) {
	if(tradc[phase]!=0 && dtradc[phase]<=l+lseq+2) {
		k=0;
		if(dtradc[phase]<=ftrad[k]) {
			k++;
			if(dtradc[phase]<=ftrad[k]) k++;
			}
		ntrad[k]=TRUE;
		i=dtradc[phase];
		strcpy(codon,"NNN");
		while(i<=l+lseq+2) {
			j=i;
			codon[0]=seq[i-l-1];
			if(i-l-2>=0) codon[1]=seq[i-l-2];
			if(i-l-3>=0) codon[2]=seq[i-l-3];
			if(i-l-3>=0) transl[k][i-l-2] = codaa(codon,codgen);
			i+=3;
			}
		dtradc[phase]=j+3; ftrad[k]=j;
		}
	}
/* toutes les editions du brin complementaire */
if(lseq%10 != 0) memset(seq+lseq,' ',102-lseq);
memset(aux,' ',130);
j=1; i=10;
while(j<=lseq) {
	k=lseq-j+1; if(k>10) k=10;
	memcpy(aux+i,seq+j-1,k);
	i+=(k+1); j+=k;
	}
aux[i]=0;
(*outline)(aux); /* seq complementaire */
memset(aux,' ',130);
for(j=0; j<=2; j++) {
	if(ntrad[j]) {
		i=10; pos=1;
		while(pos<=lseq+1) {
			k=lseq-pos+2; if(k>10) k=10;
			memcpy(aux+i,transl[j]+pos-1,k);
			i+=(k+1); pos+=k;
			}
		while(aux[i-1]==' ') i--;
		aux[i]=0;
		(*outline)(aux); /* lignes de traduction */
		}
	}
for(i=0; i<=2; i++) {
	if(decoup[i]) {
		k=larg; while(k>0 && ligne[i][k-1]==' ') k--;
		ligne[i][k]=0; 
		if(init_in_lignes) { char *p; /* enlever les ! restant dans les decoupages */
			while( (p = strchr(ligne[i], '!') ) != NULL ) *p = '<';
			}
		(*outline)(ligne[i]); /* lignes de decoupage */
		}
	}
fin_seq_frag:
(*outline)("");
/* passage a la ligne suivante de la sequence */
current_endpoint=endpoint_in_next_line;
l += basesperline;
if(l<longueur) goto next_seq_piece;

/* liberer la memoire des points de decoupage */
current_endpoint=debut_endpoints;
while(current_endpoint!=NULL) {
	one_endpoint=current_endpoint->next;
	free(current_endpoint);
	current_endpoint=one_endpoint;
	}
return;
}


static int is_5_partial(int numseq)
{
static int num_5_partial = 0;
int point;

if(num_5_partial == 0) num_5_partial = iknum("5'-PARTIAL", kkey);
/* la seq est-elle 5'-PARTIAL ? */
readsub(numseq);
point = psub->plkey;
while(point != 0) {
	readshrt(point);
	if(pshrt->val == num_5_partial) return TRUE;
	point = pshrt->next;
	}
return FALSE;
}
