/*
Pour compresser une ou plusieurs divisions en enlevant les seqs 
qui ne sont pas pointees.
Typiquement pour la division des new, mais aussi pour n'importe quelle autre.
Programme correct pour genbank, embl et swissprot.
Usage:
	compressnewdiv  div_name...
*/

#include "dir_acnuc.h"
#include <unistd.h>

typedef struct {
	off_t addr;
	off_t new_addr;
	int taille;
	int locus;
	int div;
	} extradata;

#ifdef unix
#include <sys/types.h>
#include <sys/stat.h>
#endif

extradata **calc_todo(int *targets, int n_targets, int *p_maxi, int *p_total);
void compress_div(extradata **liste, int total, int target_div, char *buffer,
	char *outbuf, int outsize, int *targets);
void tri(extradata **todo, int total);
int sort_by_addr(const void *p1, const void *p2);
void set_pointers(extradata **todo, unsigned *mere_delta, int total);
void arret_ok(char *message);
char *get_div_fname(int target_div);
#ifdef unix
int add_u_w(char *fname);
void mem_old_mode(char *fname);
DIR_FILE *reopenwr(DIR_FILE *kan);
void reset_old_modes(void *to_reset);

void *to_reset;
#endif

void *init_dynlist(void);
char *find_g_dynlist(char *name, void *tree, int create_if_not_found, 
	void **node);
void *dynlist_get_extra(void *node);
void dynlist_set_extra(void *node, void *extra);
char *next_dynlist(void *arbre, void **p_noeud);
int sizeof_dynlist(void *arbre);
char *prepare_env_var(char *env_var);
void gcgors(char *type_fic, int div, int stop_if_error);


int wid_line;


int main(int argc, char **argv)
{
extradata **todo;
int target_div, total;
unsigned *mere_delta;
int *targets, n_targets, num;
char *newdivname;
int maxi, outsize;
char *buffer, *outbuf;
time_t heure;

if(argc < 2) {
	printf("Usage: compressnewdiv  division_name...\n");
	exit(0);
	}
#ifdef unix
dir_acnucopen("RO");
to_reset = init_dynlist();
kloc = reopenwr(kloc);
if(!( nbrf || swissprot)) kext = reopenwr(kext);
else kext = (DIR_FILE *)1;
kshrt = reopenwr(kshrt);
if(kloc==NULL || kext==NULL || kshrt==NULL) 
	arret_ok("Ecriture impossible dans fichiers indexs");
#else
dir_acnucopen("WP");
#endif
if(!flat_format) arret_ok("Banque doit etre au format plat");
if(!big_annots) arret_ok("Banque doit etre au format BIG ANNOTS");
n_targets = argc - 1;
targets = (int *)malloc(n_targets * sizeof(int));
for(num = 1; num < argc; num++) {
	newdivname = argv[num];
	for(target_div=0; target_div<=divisions; target_div++) 
		if(strcmp(gcgname[target_div], newdivname) == 0) break;
	if(target_div > divisions) {
		fprintf(stderr, "bad division name %s\n", newdivname);
		exit(ERREUR);
		}
#ifdef unix
	{	char *fname; FILE *tmp;
	fname = get_div_fname(target_div);
	tmp = fopen(fname, "r+");
	if(tmp == NULL) {
		mem_old_mode(fname);
		add_u_w(fname);
		tmp = fopen(fname, "r+");
		}
	if(tmp == NULL)
		arret_ok("Ecriture impossible dans fichier plat");
	fclose(tmp);
	}
#endif
	targets[num - 1] = target_div;
	}

todo = calc_todo(targets, n_targets, &maxi, &total);
tri(todo, total);
mere_delta = (unsigned *)calloc(nseq + 1, sizeof(unsigned));
if(mere_delta == NULL) arret_ok("Pas assez de memoire");
buffer = (char *)malloc(maxi + 1);
if(buffer == NULL) arret_ok("Pas assez de memoire pour buffer");
outsize = 10000 * 1024;
if(2 * maxi > outsize) outsize = 2 * maxi;
outbuf = (char *)malloc(outsize + 1);
if(outbuf == NULL) arret_ok("Pas assez de memoire pour outbuf");
for(num = 0; num < n_targets; num++) 
	compress_div(todo, total, targets[num], buffer, outbuf, outsize, 
		targets);
free(buffer);
free(outbuf);
set_pointers(todo, mere_delta, total);
dir_acnucclose();
reset_old_modes(to_reset);
printf("Normal end\n");
return 0;
}


extradata **calc_todo(int *targets, int n_targets, int *p_maxi, int *p_total)
{
int totloc, numloc, comp_wid, rank, target_div;
off_t addr, deb_seq, fpos;
extradata **todo;
extradata *extra;
int length, taille, maxi;
int total;

totloc = read_first_rec(kloc, NULL);
todo = (extradata **)malloc(totloc * sizeof(extradata *));
if(todo == NULL) arret_ok("Memoire alloc todo");
maxi = 0;
wid_line = 100;  /* calcul exact plus loin */
comp_wid = TRUE;
total = 0;
for(numloc = 2; numloc <= totloc; numloc++) {
	readloc(numloc);
	if(ploc->sub == 0) continue;
	for(rank = 0; rank < n_targets; rank++)
		if(ploc->div == targets[rank]) break;
	if(rank >= n_targets) continue;
	target_div = ploc->div;
	addr = (unsigned)ploc->pinf;
	readsub(ploc->sub);
	extra = (extradata *)malloc(sizeof(extradata));
	if(extra == NULL) {
		arret_ok("Pas assez de memoire");
		}
	todo[total++] = extra;
	deb_seq = (unsigned)ploc->pnuc;
	length = psub->length;
	if( (length >= 60) && comp_wid) {
		read_annots((long)deb_seq, target_div);
		fpos = ftello(divannot[target_div]);
		wid_line = fpos - deb_seq;
		comp_wid = FALSE;
		}
	taille = ((length - 1)/60 + 1) * wid_line + 3;
	taille = (int)(deb_seq - addr) + taille;
	if(taille > maxi) maxi = taille;
	extra->addr = addr;
	extra->locus = numloc;
	extra->taille = taille;
	extra->div = rank;
	}
*p_maxi = maxi;
*p_total = total;
return todo;
}


void compress_div(extradata **liste, int total, int target_div, char *buffer,
	char *outbuf, int outsize, int *targets)
{
extradata *extra;
int fd, num, taille, decal;
off_t current_new, old_size;
char *newdiv_fname;
FILE *newdiv;
char *fin, *p;
int outlen, needwrite;
off_t current_out;

newdiv_fname = get_div_fname(target_div);
if( annotopened[target_div] ) {
	fclose(divannot[target_div]);
	annotopened[target_div] = 0;
	}
newdiv = fopen(newdiv_fname, "r+");
if(newdiv == NULL) {
	char mess[200];
	sprintf(mess, "Acces ecriture non autorise sur %s", newdiv_fname);
	arret_ok(mess);
	}
fd = fileno(newdiv);
old_size = lseek(fd, 0, SEEK_END);
current_new = outlen = current_out = 0;
needwrite = FALSE;
for(num = 0; num < total; num++) {
	extra = liste[num];
	if( targets[ extra->div ] != target_div ) continue;
	extra->new_addr = current_new;
	if(extra->new_addr != extra->addr) needwrite = TRUE;
	lseek(fd, extra->addr, SEEK_SET); 
	read(fd, buffer, extra->taille); 
	buffer[extra->taille] = 0;
	decal = extra->taille - 2 * wid_line;
	p = (decal > 0 ? buffer + decal : buffer);
	fin = strstr(p, "\n//");
	if(fin == NULL) {
		fin = buffer + extra->taille - 4;
		readloc(extra->locus);
		readsub(ploc->sub);
		fprintf(stderr,"Warning: // pas trouve pour seq %.16s\n", 
			psub->name);
		}
	taille = fin - buffer + 4;
	buffer[taille - 1] = '\n';
	if(outlen + taille >= outsize) {
		if(needwrite) {
			lseek(fd, current_out, SEEK_SET);
			write(fd, outbuf, outlen);
			}
		current_out += outlen;
		outlen = 0;
		}
	memcpy(outbuf + outlen, buffer, taille);
	outlen += taille;
	current_new += taille;
	}
if(outlen > 0) {
	if(needwrite) {
		lseek(fd, current_out, SEEK_SET);
		write(fd, outbuf, outlen);
		}
	current_out += outlen;
	outlen = 0;
	}
if(current_new != old_size) ftruncate(fd, current_new);
fclose(newdiv);
printf("File: %s\nold size: %lu new size: %lu difference: %lu\n", 
	newdiv_fname, 
	(unsigned long)old_size, (unsigned long)current_new, 
	(unsigned long)(old_size - current_new) );
fflush(stdout);
}


void tri(extradata **todo, int total)
{
qsort(todo, total, sizeof(extradata *), sort_by_addr);
}


int sort_by_addr(const void *p1, const void *p2)
{
extradata *q1, *q2;
int retval;

q1 = *(extradata **)p1;
q2 = *(extradata **)p2;

if(q1->div != q2->div) return q1->div - q2->div;

if(q1->addr > q2->addr) retval = 1;
else if (q1->addr == q2->addr) retval = 0;
else retval = -1;
return retval;
}


void set_pointers(extradata **todo, unsigned *mere_delta, int total)
{
int locus;
unsigned delta;
int mere, num, point, l, pre_mere;
char *q, nom_mere[L_MNEMO + 1], aux[L_MNEMO + 1];

/* process LOCUS */
#ifdef unix
dir_set_mmap(kloc);
#endif
for(num = 0; num < total; num++) {
	locus = todo[num]->locus;
	readloc(locus);
	mere = ploc->sub;
	ploc->pinf = (int)todo[num]->new_addr;
	delta = (unsigned)( todo[num]->addr - todo[num]->new_addr );
	mere_delta[mere] = delta;
	ploc->pnuc = (unsigned)ploc->pnuc - delta;
	if(delta != 0) writeloc(locus); 
	}
#ifdef unix
dir_set_normal(kloc);
#endif

/* process SHORTL & EXTRACT */
strcpy(nom_mere, "xxxx"); l = 4;
for(num = 2; num <= nseq; num++) {
	readsub(num);
	if(psub->length == 0 || psub->pext == 0) continue;
	if(psub->pext < 0) {
		memcpy(nom_mere, psub->name, L_MNEMO); nom_mere[L_MNEMO] = 0;
		l = trim_key(nom_mere); nom_mere[l] = '.';
		pre_mere = num;
		continue;
		}
	point = psub->pext;
	while(point != 0) {
		readext(point);
		mere = pext->mere;
		if( (delta = mere_delta[mere] ) != 0) {
			pext->pnuc = (unsigned)pext->pnuc - delta;
			writeext(point); 
			}
		point = pext->next;
		}
	if(strncmp(psub->name, nom_mere, l) == 0) 
		mere = pre_mere;
	else	{
		memcpy(aux, psub->name, L_MNEMO);
		q = strchr(aux, '.');
		*q = 0;
		mere = isenum(aux);
		readsub(num);
		}
	if( (delta = mere_delta[mere] ) != 0) {
		readshrt(psub->plinf);
		pshrt->val = (unsigned)pshrt->val - delta;
		writeshrt(psub->plinf); 
		}
	}

}


void arret_ok(char *message)
{
fprintf(stderr, "%s\nBanque inchangee\n", message);
exit(0);
}


char *get_div_fname(int target_div)
{
static char newdiv_fname[200];

strcpy(newdiv_fname, prepare_env_var("gcgacnuc") );
strcat(newdiv_fname, gcgname[target_div]);
if(genbank)
	strcat(newdiv_fname, ".seq");
else
	strcat(newdiv_fname, ".dat");
return newdiv_fname;
}


#ifdef unix

int add_u_w(char *fname)
{
return chmod(fname, S_IRUSR|S_IWUSR|S_IRGRP|S_IROTH);
}


void mem_old_mode(char *fname)
{
struct stat stat_data;
void *node;
mode_t *data;

stat(fname, &stat_data);
find_g_dynlist(fname, to_reset, TRUE, &node);
data = (mode_t *)malloc(sizeof(mode_t));
*data = stat_data.st_mode;
dynlist_set_extra(node, data);
}

DIR_FILE *reopenwr(DIR_FILE *kan)
{
static char fname[200];
size_t rl;

strcpy(fname, kan->filename);
rl = kan->record_length;
dir_close(kan);
kan = dir_open(fname, NULL, "r+", rl, 8*1024);
if(kan == NULL) {
	mem_old_mode(fname);
	add_u_w(fname);
	kan = dir_open(fname, NULL, "r+", rl, 8*1024);
	}
return kan;
}


void reset_old_modes(void *to_reset)
{
void *loop = NULL;
char *p;
mode_t *extra;

while( (p = next_dynlist(to_reset, &loop)) != NULL) {
	extra = (mode_t *)dynlist_get_extra(loop);
	chmod(p, *extra);
	}
}

#endif
