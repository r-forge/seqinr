/*
	BIBLIO POUR ACCES DIRECT BUFFERISE A DES FICHIERS
*/
#ifndef DIR_IO_H
#define DIR_IO_H

#if defined(unix) && ( defined(__sun) || defined(linux) )
	/* so that open lseek ftello fseeko use large files */ 
#define  _FILE_OFFSET_BITS  64  
#elif defined(__vms)  /* axp/vms */
#undef vms
#define vms
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef unix 
#include <unistd.h>
#if defined(__alpha) /* on alpha ftello == ftell && fseeko == fseek */
#define ftello ftell
#define fseeko fseek
#endif

#elif defined(macintosh)        /* Mac OS */
#include <unistd.h>
typedef long off_t;
#define ftello ftell
#define fseeko fseek


#elif defined(__INTEL__)
#include <fcntl.h>
#define ftello ftell
#define fseeko fseek

#elif defined(vms)
#include <unixio.h>
#endif

#if __INTEL__ || macintosh /* sur Mac et sur PC: fgets speciale */
char *my_fgets(char *s, int n, FILE *f);
#define fgets my_fgets
#endif

#ifndef SEEK_SET
#define SEEK_SET 0
#define SEEK_CUR 1
#define SEEK_END 2
#endif
#ifdef EXIT_FAILURE
#define ERREUR EXIT_FAILURE
#else
#define ERREUR 1
#endif

typedef struct {
	FILE *fich;
	off_t filepos; 
	} DIR_FICH_POS;
typedef struct {
	DIR_FICH_POS *fich_pos;
	int fd;   /* unix file descriptor associated to stream fich */
	size_t record_length;  /* en bytes */
	unsigned curr_buff_size; /* in records */
	int max_buff_size; /* in records (ou -1 si mode mmap) */
	unsigned first_buff_rec; /* numero du 1er record du fichier en memoire*/
	char *buffer; /* memoire */
	int buff_offset, modif_recs_tot; /* offset en bytes jusqu'a 1er record
			modifie et non ecrit; nbre de records modifies */
	char *tot_rec;       /* buffer pour contenu du 1er record */
	int tot_rec_changed; /* vrai ssi 1er record a ete modifie et non ecrit*/
	off_t curr_file_size; /* taille courante du fichier */
	char *filename;   /* le nom du fichier quand il a ete ouvert */
	} DIR_FILE;

/* prototypes pour interface c */

char *prepare_env_var(char *env_var);

DIR_FILE *dir_open(char *fname, char *context, char *mode, size_t record_length, size_t bytes_per_buffer);
/* rend struct ou NULL si erreur quelconque
	context: var d'environ pour unix ou logical pour vms;
		dir courante si NULL  !!! EN MAJUSCULES POUR VMS!!!
	mode: "rb" ou "wb+" etc
*/

int dir_read(DIR_FILE *fich, int firstrec, int recnum, void *record);
/* rend nbre de records lus (-1 indique erreur, 0 indique EOF)
	firstrec: num du 1er record a lire
	recnum: nbre de records a lire
	record: addresse ou copier ce qui est lu
*/

void *dir_read_buff(DIR_FILE *fich, int firstrec, int *num_lu);
/* lecture direct dans buffer: ok pour readonly, dangereux si modif
rend pointeur vers partie du buffer contenant donnees cherchees
	firstrec: num 1er record cherche
	*num_lu: (retour) nbre de records disponibles dans le buffer
		-1 indique erreur
		0 indique EOF
*/

int dir_write(DIR_FILE *fich, int firstrec, int recnum, void *record);
/* rend 0 ssi ok
	firstrec: num 1er record a ecrire
	recnum: nbre de records a ecrire
	record: adresse de memoire contenant donnees a ecrire
*/

int dir_close(DIR_FILE *fich);
/* rend 0 ssi ok
TOUJOURS APPELER SI ECRITURE FAITE! sinon le fichier n'est pas a jour
*/

int dir_flush(DIR_FILE *fich);
/* rend 0 ssi ok */

int dir_resize_buff(DIR_FILE *fich, void *buffer, size_t buff_size);
/* pour changer la taille du buffer d'un fichier 
(meme apres acces a ce fichier!)
rend 0 ssi ok
*/

int dir_set_mmap(DIR_FILE *fich);
/* pour acceder au fichier en mode mmap'ed (unix only) 
*/

int dir_set_normal(DIR_FILE *fich);
/* pour remettre un fichier mmap'ed en mode normal
*/
#endif /* DIR_IO_H */
