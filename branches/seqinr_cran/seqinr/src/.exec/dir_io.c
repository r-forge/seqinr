/*
	BIBLIO POUR ACCES DIRECT BUFFERISE A DES FICHIERS
	SOUS UNIX POSSIBILITE D'ACCES EN MODE MMAP'ed

marche tant que le nbre de records du fichier 
est representable par unsigned int (~ 4 milliards)
accepte donc des fichiers plus grands que 2 GB si nbre de rec est OK ET si
est compile tel que le type long soit sur 64 bits (gcc -m64)

*/
#include "dir_io.h"

#ifdef unix
#include <sys/types.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>

#ifndef MAP_FAILED  /* needed for SGI */
#define MAP_FAILED -1
#endif
#define USE_MAP
#endif

/* sometimes missing prototype from standard include files *
#ifndef fileno
extern int fileno(FILE *);
#endif
*/


/* prototypes */
DIR_FILE *dir_open(char *fname, char *context, char *mode, size_t record_length,
	size_t bytes_per_buffer);
int dir_read(DIR_FILE *fich, int firstrec, int recnum, void *record);
void *dir_read_buff(DIR_FILE *fich, int firstrec, int *num_lu);
int dir_write(DIR_FILE *fich, int firstrec, int recnum, void *record);
static int write_mod_part(DIR_FILE *fich);
int dir_flush(DIR_FILE *fich);
int dir_close(DIR_FILE *fich);
char *prepare_env_var(char *env_var);
int dir_set_normal(DIR_FILE *fich);
#ifdef USE_MAP
int dir_set_mmap(DIR_FILE *dfile);
static int dir_write_map(DIR_FILE *fich, off_t position, char *record, 
	size_t numbytes);
static int init_normal(DIR_FILE *fich);
static size_t check_mmap_arg(off_t file_size);
#endif


#ifdef USE_MAP
int dir_set_mmap(DIR_FILE *dfile)
/* pour mettre un fichier prealablement ouvert en mode mmap'ed
rend != 0 si erreur, 0 si OK
*/
{
caddr_t pointer;
size_t rsize, mmap_size;

if( dfile->max_buff_size == -1) return 0; /* deja fait */
dir_flush(dfile);
rsize = dfile->record_length;
if( dfile->curr_file_size < rsize ) {
	lseek(dfile->fd, 0, SEEK_SET);
	dfile->fich_pos->filepos = 0;
	if( write(dfile->fd, dfile->tot_rec, rsize) != rsize) return 1;
	dfile->curr_file_size = rsize;
	}
mmap_size = check_mmap_arg(dfile->curr_file_size);
if( mmap_size == 0 ) return 1;
pointer = mmap((caddr_t) 0, mmap_size, (PROT_READ | PROT_WRITE),
                         MAP_SHARED, dfile->fd, 0);
if(pointer == (caddr_t) MAP_FAILED) {
	if(errno == EACCES) { /* cas d'un fichier ouvert en lecture */
		pointer = mmap((caddr_t) 0, mmap_size, 
			PROT_READ, MAP_SHARED, dfile->fd, 0);
		}
	if(pointer == (caddr_t) MAP_FAILED) return 1;
	}
free(dfile->buffer);
free(dfile->tot_rec);
dfile->buffer = (char *)pointer;
dfile->max_buff_size = -1; /* indique mode mmap'ed */
/* taille en records de la partie mappee.
Peut etre > taille reelle du fichier car elle est agrandie par paquets
*/
dfile->curr_buff_size = dfile->curr_file_size / rsize;
return 0;
}


static size_t check_mmap_arg(off_t file_size)
{
size_t mmap_size;
off_t check;

mmap_size = (size_t)file_size;
check = (off_t)mmap_size;
return (check == file_size ? mmap_size : 0);
}

#endif


DIR_FILE *dir_open(char *fname, char *context, char *mode, size_t record_length,
	size_t bytes_per_buffer)
{
FILE *fich;
int fd, open_mode;
DIR_FILE *dir_file;
off_t total;
int max_recs_per_buff;
char *fullname;
#if defined(unix)
mode_t old_mask, file_mask;
#endif

#ifdef vms
char  mrs_option[15];
sprintf(mrs_option, "mrs=%d", record_length);
#endif

if(context != NULL) {
	fullname = prepare_env_var(context);
	if(fullname == NULL) return NULL;
	strcat(fullname,fname);
	}
else	fullname = fname;

#if defined(unix)  /* unix: ouverture par open */
	if( strchr(mode, '+') != NULL || strchr(mode, 'w') != NULL ) {
		open_mode = O_RDWR | O_CREAT;
		old_mask = umask(0); /* old_mask = mask courant */
		umask(old_mask); /* remettre mask courant */
		file_mask = 0666 - old_mask; /* mask pour fichier a creer */
		}
	else	{
		open_mode = O_RDONLY;
		file_mask = 0; /* en lecture mask est inutile */
		}
	fd = open(fullname, open_mode, file_mask); 
	if( fd == -1) return NULL;
	fich = NULL;
#elif defined(__INTEL__)  
	if( strchr(mode, '+') != NULL || strchr(mode, 'w') != NULL ) {
		open_mode = O_RDWR | O_CREAT | O_BINARY;
		}
	else	{
		open_mode = O_RDONLY | O_BINARY;
		}
	fd = open(fullname, open_mode); 
	if( fd == -1) return NULL;
	fich = NULL;
#else /* non unix: ouverture par fopen */
	fich = fopen(fullname, mode
#ifdef vms
		,"rfm=fix", mrs_option
#endif
		);
	if(fich == NULL) return NULL;
	setbuf(fich, NULL);
	fd = fileno(fich);
#endif

if( (dir_file = (DIR_FILE *)malloc(sizeof(DIR_FILE))) == NULL) return NULL;
if((dir_file->tot_rec = (char *)malloc(record_length)) == NULL) return NULL;
if((dir_file->filename = (char *)malloc(strlen(fullname)+1)) == NULL) return NULL;
if((dir_file->fich_pos = (DIR_FICH_POS *)malloc(sizeof(DIR_FICH_POS))) == NULL) return NULL;
strcpy(dir_file->filename, fullname);
dir_file->fd = fd;
dir_file->fich_pos->fich = fich;
max_recs_per_buff = bytes_per_buffer/record_length;
dir_file->record_length = record_length;
dir_file->curr_buff_size = 0;
dir_file->max_buff_size = max_recs_per_buff;
total = lseek(dir_file->fd, 0, SEEK_END);
if(total == -1)
	dir_file->curr_file_size = 0;
else
	dir_file->curr_file_size = total;
lseek(dir_file->fd, 0, SEEK_SET);
dir_file->fich_pos->filepos = 0;
if(dir_file->curr_file_size >= record_length) {
	if( read(dir_file->fd, dir_file->tot_rec, record_length) == -1) 
			return NULL;
	dir_file->fich_pos->filepos = record_length;
	}
dir_file->tot_rec_changed = 0;
dir_file->modif_recs_tot = 0;
dir_file->buffer = (char *)malloc(max_recs_per_buff * record_length);
if(dir_file->buffer == NULL) return NULL;
return dir_file;
}


int dir_read(DIR_FILE *fich, int firstrec_p, int recnum, void *record)
/* 
returns -1 in case of error
        # of records read else
*/
{
char *buff_pos;
int lu, totlu = 0;
size_t rsize;
unsigned firstrec = (unsigned) firstrec_p;

rsize = fich->record_length;
if(fich->max_buff_size == -1) { /* mode mmap'ed */
	int num_lu = 0;
	off_t position = (firstrec-1) * (off_t)rsize;
	unsigned last = fich->curr_file_size / rsize;
	unsigned total =  last - (firstrec - 1);
	if(total > recnum) total = recnum;
	if(total > 0) 
		memcpy(record, fich->buffer + position, total * rsize);
	else 	total = 0;
	return (int) total;
	}
else	{
	do 	{
		buff_pos=dir_read_buff(fich, (int)firstrec, &lu);
		if( lu == -1) return -1;
		if(lu>recnum)lu=recnum;
		if(lu>0)memcpy(record,buff_pos,lu*rsize);
		record= (char *)record + lu*rsize;
		recnum -= lu;
		firstrec += lu;
		totlu += lu;
		}
	while(recnum>0 && lu > 0);
	}
return totlu;
}


void *dir_read_buff(DIR_FILE *fich, int firstrec_p, int *num_lu)
/*
 *num_lu is returned at -1 in case of error and at # of records read else
*/
{
unsigned lastrec;
int lus;
size_t rsize;
off_t filepos;
unsigned firstrec = (unsigned) firstrec_p;

rsize=fich->record_length;

if(fich->max_buff_size == -1) { /* mode mmap'ed */
	long position = (firstrec - 1) * rsize;
	if(firstrec > fich->curr_file_size / rsize) goto erreur;
	*num_lu = (int)( (fich->curr_file_size - position) / rsize );
	return fich->buffer + position;
	}

lastrec = fich->first_buff_rec + fich->curr_buff_size;
if(firstrec == 1) {
	*num_lu = 1;
	return fich->tot_rec;
	}
else if(firstrec >= fich->first_buff_rec && firstrec < lastrec) {
	*num_lu = (int)( lastrec - firstrec );
	return fich->buffer + (firstrec - fich->first_buff_rec)*rsize;
	}
else 	{
	if( write_mod_part(fich) )goto erreur;
	filepos = (firstrec-1)*(off_t)rsize;
	if(filepos < fich->curr_file_size) {
		if(fich->fich_pos->filepos != filepos)
		    if(lseek(fich->fd,filepos,SEEK_SET) != filepos)goto erreur;
		lus = read(fich->fd, fich->buffer, rsize*fich->max_buff_size);
		if(lus == -1) goto erreur;
		fich->curr_buff_size = lus / rsize;
		fich->fich_pos->filepos = filepos + lus;
		}
	else	fich->curr_buff_size=0;
	fich->first_buff_rec=firstrec;
	*num_lu = (int)fich->curr_buff_size;
	return fich->buffer;
	}
erreur:
*num_lu= -1;
return NULL;
}


int dir_write(DIR_FILE *fich, int firstrec_p, int recnum, void *record)
/*
returns nonzero in case of any error
*/
{
size_t rsize;
unsigned lastrec, newsize;
int num_done, offset;
unsigned firstrec = (unsigned) firstrec_p;

rsize=fich->record_length;

if(fich->max_buff_size == -1) { /* mode mmap'ed */
	off_t position = (firstrec-1)*(off_t)rsize;
	off_t last = position + recnum*(off_t)rsize;
	if(firstrec - 1 + recnum > fich->curr_buff_size) {
#ifdef USE_MAP
		if( dir_write_map(fich, position, record, (size_t)(last-position) ) ) 
#endif
			return 1;
		}
	else 	{
		memcpy(fich->buffer+position, record, (last-position));
		if(fich->curr_file_size < last ) fich->curr_file_size = last;
		}
	return 0;
	}

if(fich->curr_buff_size == 0) fich->first_buff_rec=firstrec;
lastrec=fich->first_buff_rec + fich->max_buff_size;
if(firstrec==1) {
	memcpy(fich->tot_rec, record, rsize);
	fich->tot_rec_changed=1;
	num_done=1;
	}
else if (firstrec>=fich->first_buff_rec && 
				firstrec < lastrec) {
	num_done=recnum;
	if(firstrec+recnum > lastrec) num_done= (int)( lastrec-firstrec );
	memcpy(fich->buffer + (firstrec - fich->first_buff_rec)*rsize,
		record,rsize*num_done);
	newsize=firstrec+num_done - fich->first_buff_rec;
	if(newsize>fich->curr_buff_size) fich->curr_buff_size = newsize;
	offset = (int)( (firstrec - fich->first_buff_rec)*rsize );
	if(!fich->modif_recs_tot) {
		fich->buff_offset = offset;
		fich->modif_recs_tot=num_done;
		}
	else	{
		int oldlastmod, newlastmod;
		oldlastmod=fich->buff_offset/rsize + fich->modif_recs_tot;
		if(offset < fich->buff_offset) fich->buff_offset = offset;
		newlastmod= offset/rsize + num_done;
		if(newlastmod < oldlastmod) newlastmod = oldlastmod;
		fich->modif_recs_tot= newlastmod - fich->buff_offset/rsize;
		}
	}
else	{
	char *pos=NULL;
	if( dir_read(fich,(int)firstrec,0,pos) == -1) return 1;
	if( dir_write(fich,(int)firstrec,recnum,record) )return 1;
	num_done=recnum;
	}
if(num_done<recnum) 
	if( dir_write(fich, (int)( firstrec+num_done ), recnum-num_done,
	(char *)record + num_done*rsize) ) return 1;
return 0;
}


static int write_mod_part(DIR_FILE *fich)
{
int numbytes;
size_t rsize;
off_t newsize, filepos;

rsize = fich->record_length;
if(fich->modif_recs_tot) {
	if(fich->curr_file_size == 0) {
		if(fich->fich_pos->filepos != 0) lseek(fich->fd, 0, SEEK_SET);
		if( write(fich->fd, fich->tot_rec, rsize) != rsize) return 1;
		fich->tot_rec_changed = 0;
		fich->curr_file_size = rsize;
		fich->fich_pos->filepos = rsize;
	        }
	filepos = (fich->first_buff_rec-1) * (off_t)rsize + fich->buff_offset;
	if(fich->fich_pos->filepos != filepos)
		if(lseek(fich->fd, filepos, SEEK_SET) != filepos)return 1;
	numbytes = rsize * fich->modif_recs_tot;
	if( write(fich->fd, fich->buffer + fich->buff_offset, numbytes)
		!= numbytes) return 1;
	newsize = filepos + numbytes;
	fich->fich_pos->filepos = newsize;
	if(newsize > fich->curr_file_size) fich->curr_file_size = newsize;
	fich->modif_recs_tot = 0;
	}
return 0;
}


#ifdef USE_MAP
static int dir_write_map(DIR_FILE *fich, off_t position, char *record, 
	size_t numbytes)
{
/* nbre de bytes a ecrire d'avance en fin de fichier */
#define PREAMPTION 1000000 
#define CHUMPS 50000
static char dummy[CHUMPS];
char *old_map;
int rsize, num, new_recs;
caddr_t pointer;
off_t new_lmap;
size_t old_lmap, mmap_size;

rsize = fich->record_length;
old_map = fich->buffer;
old_lmap = fich->curr_buff_size * (size_t) rsize;
munmap((caddr_t) old_map, old_lmap);
lseek(fich->fd, position, SEEK_SET);
if( write(fich->fd, record, numbytes) != numbytes) return 1;
fich->curr_file_size = position + numbytes;
new_recs = PREAMPTION / rsize; if(new_recs < 1) new_recs = 1;
for(num = 0; num < PREAMPTION; num += CHUMPS)
	if( write(fich->fd, dummy, CHUMPS) != CHUMPS) return 1;
new_lmap = fich->curr_file_size + new_recs * rsize;
mmap_size = check_mmap_arg(new_lmap);
if(mmap_size == 0) pointer = (caddr_t) MAP_FAILED;
else	pointer = mmap((caddr_t) 0, mmap_size, (PROT_READ|PROT_WRITE),
		MAP_SHARED, fich->fd, 0);
if(pointer != (caddr_t) MAP_FAILED){
	fich->buffer = pointer;
	fich->curr_buff_size = new_lmap / rsize;
	return 0; /* sortie avec succes */
	}

/* pas assez de memoire pour mode mmap, revenir en mode simple */
return init_normal(fich);
}


static int init_normal(DIR_FILE *fich)
{
int rsize;
off_t length;

rsize = fich->record_length;
length = lseek(fich->fd, 0, SEEK_END);
if(length > fich->curr_file_size ) 
	ftruncate(fich->fd, fich->curr_file_size);
if((fich->tot_rec = malloc(rsize))==NULL) return 1;
fich->max_buff_size = (8*1024)/rsize;
fich->curr_buff_size = 0;
lseek(fich->fd, 0, SEEK_SET);
fich->fich_pos->filepos = 0;
if(fich->curr_file_size >= rsize) {
	if( read(fich->fd, fich->tot_rec, rsize) == -1) 
			return 1;
	fich->fich_pos->filepos = rsize;
	}
fich->tot_rec_changed = 0;
fich->modif_recs_tot = 0;
fich->buffer = (char *)malloc(fich->max_buff_size * rsize);
if(fich->buffer == NULL)return 1;
return 0;
}
#endif


int dir_set_normal(DIR_FILE *fich)
/* remettre fichier mmap'ed en mode normal */
{
#ifdef USE_MAP
char *old_map;
size_t old_lmap, new_lmap;
int rsize;

if(fich->max_buff_size == -1) {
	rsize = fich->record_length;
	old_map = fich->buffer;
	old_lmap = fich->curr_buff_size * (size_t) rsize;
	munmap((caddr_t) old_map, old_lmap);
	return init_normal(fich);
	}
#endif
return 0;
}


int dir_flush(DIR_FILE *fich)
{
size_t rsize;

if(fich->max_buff_size == -1) return 0;

rsize = fich->record_length;
if( write_mod_part(fich) ) return 1;
if (fich->tot_rec_changed) {
	lseek(fich->fd, 0, SEEK_SET);
	if( write(fich->fd, fich->tot_rec, rsize) != rsize) return 1;
	fich->tot_rec_changed = 0;
	if(fich->curr_file_size == 0) fich->curr_file_size = rsize;
	fich->fich_pos->filepos = rsize;
	}
return 0;
}


int dir_resize_buff(DIR_FILE *fich, void *buffer, size_t buff_size)
{
int rec_count;

if(fich->max_buff_size == -1) return 0; /* mode mmap'ed */
if( write_mod_part(fich) )return 1;
free(fich->buffer);
rec_count = buff_size / fich->record_length;
fich->buffer = buffer;
fich->curr_buff_size=0;
fich->max_buff_size=rec_count;
fich->modif_recs_tot=0;
return 0;
}


int dir_close(DIR_FILE *fich)
{
if(fich->max_buff_size == -1) { /* mode mmap'ed */
#ifdef USE_MAP
	size_t l_map;
	l_map = fich->curr_buff_size * (size_t) fich->record_length;
	munmap((caddr_t) fich->buffer, l_map);
	if( l_map > (size_t) fich->curr_file_size )
		ftruncate(fich->fd, fich->curr_file_size);
#endif
	;
	}
else	{
	if( dir_flush(fich) ) return 1;
	free(fich->buffer);
	free(fich->tot_rec);
	}
#if defined(unix) || defined(__INTEL__)
close(fich->fd);
#else
fclose(fich->fich_pos->fich);
#endif
free(fich->filename);
free(fich->fich_pos);
free(fich);
return 0;
}


char *prepare_env_var(char *env_var)
{
static char long_fname[200];
char *p;
size_t l;

if(env_var == NULL) return NULL;
p = getenv(env_var);
if(p == NULL) return NULL;
strcpy(long_fname, p); l = strlen(long_fname);
#ifdef unix
	if(long_fname[l-1] != '/')strcat(long_fname,"/");
#elif __INTEL__
	if(long_fname[l-1] != '\\')strcat(long_fname,"\\");
#else /* vms or Mac */
	if(long_fname[l-1] != ']' && long_fname[l-1] != ':')
						strcat(long_fname,":");
#endif
return long_fname;
}
