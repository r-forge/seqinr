#include "dir_acnuc.h"
#include <time.h>

/* prototypes */
char *paqde3( double val);
void proc_help_file(char *name);

/* global vars */
double totbib, totmeres=0, totfilles=0, totlong=0;
/* when TRUE date of operation is written in HELP files
run program with -noupdate not to write date in files
*/
int write_update; 

main(int argc, char **argv)
{
int i;
write_update =  (argc < 2 || strcmp(argv[1],"-noupdate") != 0);
	

acnucopen();

totbib=read_first_rec(kbib,NULL) - 1;

for(i=2; i<=nseq; i++) {
	readsub(i); if(psub->length==0) continue;
	if(psub->pext<=0) { /* mere */
		totmeres++;
		totlong += psub->length;
		}
	else	totfilles++;
	}

proc_help_file("HELP");
proc_help_file("HELP_WIN");
exit(0);
}


void proc_help_file(char *name)
{
char path[100], oldname[200], newname[200], ligne[100], ligne2[100],
	*debut_heure, *p;
FILE *fichin, *fichout;
int i;
time_t heure;

strcpy(oldname, prepare_env_var("acnuc"));
strcat(oldname, name);
fichin=fopen(oldname,"r");
if(fichin==NULL) return;
strcpy(newname, oldname);
strcat(newname, ".NEW");
fichout=fopen(newname,"w");
if(fichout==NULL) exit(ERREUR);
for(i=1; i<=2; i++) {
	fgets(ligne,100,fichin);
	fputs(ligne,fichout);
	}
/* process the release # line */
fgets(ligne,100,fichin);
if(write_update) {
	time(&heure);
	debut_heure = ctime(&heure);
	p = ligne; while(*p == ' ') p++;
	strcpy(ligne2,p);
	if( ( p = strstr(ligne2,"Last Updated") ) == NULL) {
		i=strlen(ligne2); p = ligne2+i-2; 
		}
	else	p--;
	while(*p == ' ') p--;
	sprintf(p+1," Last Updated: %.6s, %.4s\n",debut_heure+4, 
			debut_heure+20);
	i=strlen(ligne2) - 1; i = (80-i)/2; if(i<0) i=0;
	while(i-- > 0) fputc(' ',fichout);
	fputs(ligne2,fichout);
	}
else	fputs(ligne,fichout);

/* process the totals line */
fgets(ligne,100,fichin);
if(swissprot || nbrf) {
sprintf(ligne, "          %s amino acids; %s sequences; %s references.\n",
	paqde3(totlong), paqde3(totmeres), paqde3(totbib) );
	}
else	{
sprintf(ligne,"%s bases; %s sequences; %s subseqs; %s refers.\n",
	paqde3(totlong), paqde3(totmeres), paqde3(totfilles), 
		paqde3(totbib) );
	}
fputs(ligne,fichout);

while( fgets(ligne,100,fichin)!=NULL ) fputs(ligne,fichout);
fclose(fichin); fclose(fichout);
if( remove(oldname) != 0 ) exit(ERREUR);
if( rename(newname, oldname) != 0) exit(ERREUR);
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
