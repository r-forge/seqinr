#include "dir_acnuc.h"
#include <ctype.h>

void process_seq(int num, int translate);
void libsub_modif(int num, char *libel, int length);
char translate_init_codon(int numseq, int gc, int debut_codon /* 1, 2, or 3 */);

int main(int argc, char *argv[])
{
int tout, tloc, iloc, translate, num;
FILE *listfile;
char nom[20], *p;

if(argc < 2) {
	fprintf(stderr, 
	"Erreur: Donner le fichier liste des noms ou 'tout' en 1er arg\n");
	exit(1);
	}
tout = (strcmp(argv[1], "tout") == 0);
acnucopen();
#ifdef unix
dir_set_mmap(kshrt);
dir_set_mmap(ksub);
dir_set_mmap(kloc);
if(kext != NULL) dir_set_mmap(kext);
#endif
if(!tout)
	listfile = fopen(argv[1], "r");
else	{
	tloc = read_first_rec(kloc, NULL);
	}
translate = FALSE;
if( (!(nbrf || swissprot)) && argc >= 3 )
	translate = (*argv[2] == 't');
if( !tout ) {
	while(fgets(nom, sizeof(nom), listfile) != NULL) {
		p = nom + strlen(nom); if(*(p-1) == '\n') *(p-1) = 0;
		num = isenum(nom);
		if(num == 0) continue;
		process_seq(num, translate);
		}
	}
else	{
	for(iloc = 2; iloc <= tloc; iloc++) {
		readloc(iloc);
		if(ploc->sub == 0) continue;
		process_seq(ploc->sub, translate);
		}
	}
if(!tout) fclose(listfile);
return 0;
}


void process_seq(int num, int translate)
{
static char libel[2000], seq[200];
char *p, *q;
int step, phase, vlength, pos;

readsub(num);
if(translate) {
	phase = (psub->phase % 100 ) + 1;
	vlength = (psub->length - phase + 1) / 3;
	}
else	{
	step = 75;
	vlength = psub->length;
	}
libsub_modif(num, libel, vlength);
printf("%.100s\n", libel);
if(translate) {
	p = translate_cds(num);
	if(p == NULL) exit(ERREUR);
	vlength = strlen(p);
	for(q = p; q < p + vlength; q += 75)
		printf("%.75s\n", q);
	}
else	{
	for(pos = 1; pos <= vlength; pos += step) {
		if( gfrag(num, pos, step, seq) == 0) exit(ERREUR);
		if(!nbrf) majuscules(seq);
		puts(seq);
		}
	}
}


void libsub_modif(int num, char *libel, int length)
{
char *pos, *p;
int mere, div, deb, l;
long pflat;

readsub(num);
*libel = '>';
pos = libel + 1;
memcpy(pos, psub->name, sizeof(psub->name));
pos += sizeof(psub->name);
while( *(pos - 1) == ' ') pos--;
sprintf(pos, "  %d", length);
pos = libel + strlen(libel);
if(psub->pext > 0) { /* une fille */
	p = strchr(psub->name, '.');
	*p = 0;
	mere = isenum(psub->name);
	num = mere;
	}
seq_to_annots(num, &pflat, &div);
if( read_annots(pflat, div) == NULL) return;
next_annots(NULL);
if(nbrf)
	deb = 16;
else if (embl || swissprot) {
	while(strncmp(pinfo->line, "DE", 2) != 0) {
		next_annots(NULL);
		}
	deb = 5;
	}
else
	deb = 12;
do	{
	*(pos++) = ' ';
	l = strlen(pinfo->line + deb);
	memcpy(pos, pinfo->line + deb, l);
	pos += l; *pos = 0;
	next_annots(NULL);
	}
while(strncmp(pinfo->line, "DE", 2) == 0 || 
		strcmptrail(pinfo->line, 10, NULL, 0) == 0);
}



