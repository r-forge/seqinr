#include "dir_acnuc.h"

#define MAX_NON_TRIE 1000
#define MAX_CHAR_WIDTH 100
#define COEFF 20

typedef struct bt_node { /* binary tree node */
	char *name;
	int recnum;
	struct bt_node *before;
	struct bt_node *after;
	} bt_node;
	
typedef struct bt_node_aux { /* binary tree node with int aux data */
	bt_node noeud_binaire;
	int aux_data;
	} bt_node_aux;

/* gestion complete d'un fichier index acnuc par binary tree */
typedef struct index_bt { 
	DIR_FILE *k;
	bt_node *racine_trie;
/* partie nouvelle de l'index qui croit jusqu'a max_non_trie
puis est transferee dans racine_trie */
	bt_node *racine_non_trie; 
	int total_fichier;
	int total_trie;
	int total_non_trie;
	int last_sorted;
	int char_width;
	int max_non_trie;
	char *buffer;
	size_t taille_node;
	} index_bt;
	

/* prototypes */
int find_key(char *name, index_bt *tree, int create_if_not_found,
	int **p_aux_data);
index_bt *load_index_bt(DIR_FILE *k, int char_width, int use_aux_data);
void free_index_bt(index_bt *tree);

/* fonctions privees */
static void init_aux_data(bt_node *racine);
static bt_node *find_node(char *name, bt_node *noeud, int create_if_not_found, 
	size_t taille_node, int *maxprof);
static int trim_key_max(char *record, int max_width);
static void reset_non_trie(index_bt *tree);
static bt_node **tri_non_trie(bt_node *racine, bt_node **liste_triee);
static void add_mid_node(int bas, int haut, bt_node **liste_triee,
	bt_node *tree);
static void add_node(bt_node *noeud, bt_node *racine);
static void free_bt_node(bt_node *racine);
static void resort_partie_non_triee(index_bt *tree);
static bt_node *opt_tree(int bas, int haut, bt_node **tableau);
int tri_par_nom(const void *elt_1, const void *elt_2);
bt_node **merge_liste_noeuds(bt_node **liste_noeuds, int total, int fin_trie);


int find_key(char *target, index_bt *tree, int create_if_not_found,
	int **p_aux_data)
/* chercher target dans l'arbre tree, rendre son numero de record dans le 
fichier index correspondant. 
Si p_aux_data != NULL, rendre adresse de aux_data
Si create_if_not_found==TRUE, le record est cree dans le fichier index et 
dans les arbres binaires associes si le nom n'existait pas. 
valeur rendue = 0 si pas assez de memoire disponible pour fonctionner
Sinon valeur rendue = 0 si pas trouve
*/
{
bt_node *noeud;
int total, l, maxprof;
static char name[MAX_CHAR_WIDTH + 1];

/* tronquer target a tree->char_width et enlever les espaces terminaux */
l = strlen(target);
if(l > tree->char_width) l = tree->char_width;
memcpy(name, target, l); name[l] = 0;
if(name[l - 1] == ' ') trim_key_max(name, tree->char_width);
noeud = find_node(name, tree->racine_trie, FALSE, tree->taille_node, NULL);
maxprof=0;
if(noeud == NULL) noeud = find_node(name, tree->racine_non_trie, 
			create_if_not_found, tree->taille_node, &maxprof);
if(noeud == NULL) return 0;
if(p_aux_data != NULL) *p_aux_data = &(((bt_node_aux *)noeud)->aux_data);
if(noeud->recnum != 0) return noeud->recnum; /* noeud trouve deja existant */
/* ajout du nouveau noeud dans le fichier */
if(tree->racine_non_trie == NULL) tree->racine_non_trie = noeud;
++(tree->total_non_trie);
total = ++(tree->total_fichier);
noeud->recnum = total;
if(tree->taille_node == sizeof(bt_node_aux)) 
	((bt_node_aux *)noeud)->aux_data = 0;
write_first_rec(tree->k, total, tree->last_sorted);
memset(tree->buffer, 0, tree->k->record_length);
padtosize(tree->buffer, name, tree->char_width);
if( dir_write(tree->k, total, 1, tree->buffer) ) dir_writeerr(tree->k, total);
if( maxprof > 50 )
	/* trier partie non triee */
	resort_partie_non_triee(tree);
if(tree->total_non_trie >= tree->max_non_trie) 
	/* transferer partie non triee dans partie triee */
	reset_non_trie(tree);
return total;
}


index_bt *load_index_bt(DIR_FILE *k, int char_width, int use_aux_data)
/* charger en memoire un fichier index 
use_aux_data TRUE pour place pour un entier en plus 
dans chaque noeud arbre binaire
en fin de chargement il est entierement en memoire 
sous forme d'arbre binaire trie
valeur rendue: pointeur vers struct index_bt cree
ou NULL si pas assez de memoire
*/
{
index_bt *tree;
bt_node *node, **liste_noeuds;
char croix[60];
int num, debut, croix_width, maxprof, 
	l, total, fin_trie;

if(char_width > MAX_CHAR_WIDTH) /* limite */
	return NULL;
croix_width = ( char_width < 60 ? char_width : 60 );
memset(croix, 'x', croix_width);
tree = (index_bt *)calloc(1, sizeof(index_bt));
if(tree == NULL) return NULL;
tree->buffer = (char *)malloc(k->record_length);
if(tree->buffer == NULL) return NULL;
tree->k = k;
tree->char_width = char_width;
tree->total_fichier = read_first_rec(k, &(tree->last_sorted) );
tree->max_non_trie = MAX_NON_TRIE;
tree->taille_node = (use_aux_data ?
	sizeof(bt_node_aux)
	: sizeof(bt_node) );

/* lire tous les noms du fichier sequentiellement */
liste_noeuds = (bt_node **) malloc(tree->total_fichier * sizeof(bt_node *));
total = fin_trie = 0;
if(liste_noeuds == NULL) return NULL;
for(num = 2; num <= tree->total_fichier; num++) {
	if( dir_read(k, num, 1, tree->buffer)!=1 ) dir_readerr(k,num);
	if(memcmp(tree->buffer, croix, croix_width) == 0) continue;
	l = trim_key_max(tree->buffer, char_width) + 1;
	node = (bt_node *)calloc(1, tree->taille_node);
	if(node == NULL) return NULL;
	node->name = (char *)malloc(l);
	memcpy(node->name, tree->buffer, l);
	node->recnum = num;
	liste_noeuds[total++] = node;
	if(num == tree->last_sorted) fin_trie = total;
	}
if(tree->last_sorted <= 2) fin_trie = 0;
/* trier la partie non triee du tableau liste_noeuds */
if(total > fin_trie) qsort(liste_noeuds + fin_trie, total - fin_trie, 
	sizeof(bt_node *), tri_par_nom);
/* trier toute la liste en mergeant ses 2 parties */
if( total > fin_trie && fin_trie > 0 )
	liste_noeuds = merge_liste_noeuds(liste_noeuds, total, fin_trie);
if(liste_noeuds == NULL) return NULL;
/* fabriquer arbre binaire optimal */
if(total > 0)  tree->racine_trie = opt_tree(0, total - 1, liste_noeuds);
free(liste_noeuds);
tree->total_trie = total;
tree->total_non_trie = 0;
tree->racine_non_trie = NULL;
if(use_aux_data) init_aux_data(tree->racine_trie);

return tree;
}


int tri_par_nom(const void *elt_1, const void *elt_2)
{
return strcmp( ((*(bt_node **)elt_1))->name, ((*(bt_node **)elt_2))->name);
}


bt_node **merge_liste_noeuds(bt_node **liste_noeuds, int total, int fin_trie)
{
bt_node **new_list;
int rank_1, rank_2, rank_new;

new_list = (bt_node **)malloc(total * sizeof(bt_node *));
if(new_list == NULL) {
	free(liste_noeuds);
	return NULL;
	}
rank_1 = 0; rank_2 = fin_trie; rank_new = 0;
while(rank_new < total) {
	while( rank_1 >= fin_trie || strcmp(liste_noeuds[rank_1]->name, 
			liste_noeuds[rank_2]->name) > 0) {
		new_list[rank_new++] = liste_noeuds[rank_2++];
		if(rank_2 >= total) 
			break;
		}
	if(rank_new >= total) break;
	while( rank_2 >= total || strcmp(liste_noeuds[rank_1]->name, 
			liste_noeuds[rank_2]->name) < 0) {
		new_list[rank_new++] = liste_noeuds[rank_1++];
		if(rank_1 >= fin_trie) 
			break;
		}
	}
free(liste_noeuds);
return new_list;
}


void free_index_bt(index_bt *tree)
/* liberer memoire d'un arbre tree + flush du fichier index associe */
{
dir_flush(tree->k);
free_bt_node(tree->racine_trie);
free_bt_node(tree->racine_non_trie);
free(tree->buffer);
free(tree);
}


static void init_aux_data(bt_node *racine)
{
if(racine == NULL) return;
((bt_node_aux *)racine)->aux_data = 0;
init_aux_data(racine->before);
init_aux_data(racine->after);
}


static void free_bt_node(bt_node *racine)
{
if(racine == NULL) return;
free_bt_node(racine->before);
free_bt_node(racine->after);
free(racine->name);
free(racine);
}


static bt_node *find_node(char *name, bt_node *noeud, int create_if_not_found,
	size_t taille_node, int *maxprof)
/* trouver noeud de nom name dans l'arbre binaire qui debute en noeud 
rend pointeur vers noeud correspondant, 
si create_if_not_found==TRUE celui-ci est cree si ce nom n'existait pas
	si le noeud a ete cree, on a ->recnum = 0
	et rend NULL si pas assez de memoire
sinon rend NULL si pas trouve
*/
{
int test;

if(noeud == NULL) {
	int l;
	if( !create_if_not_found ) return NULL;
	noeud = (bt_node *)calloc(1, taille_node);
	if(noeud == NULL) return NULL;
	l = strlen(name) + 1;
	noeud->name = (char *)malloc(l);
	if(noeud->name == NULL) return NULL;
	memcpy(noeud->name, name, l);
	return noeud;
	}
test = strcmp(name, noeud->name);
if(test == 0) return noeud;
if(maxprof != NULL) (*maxprof)++;
if(test > 0) {
	if(noeud->after == NULL) {
		if(create_if_not_found) noeud->after = 
			find_node(name, NULL, TRUE, taille_node, maxprof);
		return noeud->after;
		}
	else
		return find_node(name, noeud->after, create_if_not_found, 
				taille_node, maxprof);
	}
else	{
	if(noeud->before == NULL) {
		if(create_if_not_found) noeud->before = 
			find_node(name, NULL, TRUE, taille_node, maxprof);
		return noeud->before;
		}
	else
		return find_node(name, noeud->before, create_if_not_found, 
				taille_node, maxprof);
	}
}


static void reset_non_trie(index_bt *tree)
/* transferer toute la partie non triee de l'arbre tree vers sa partie triee
*/
{
bt_node **liste_triee;
int milieu;

liste_triee = (bt_node **)malloc(tree->total_non_trie * sizeof(bt_node *) );
if(liste_triee == NULL) return;
tri_non_trie(tree->racine_non_trie, liste_triee);
if( liste_triee[tree->total_non_trie - 1] == NULL ) {
	fprintf(stderr, "ERROR: repetition in file %s\n", tree->k->filename);
	exit(ERREUR);
	}
milieu = (tree->total_non_trie - 1)/2;
if(tree->racine_trie == NULL) tree->racine_trie = liste_triee[milieu];
add_mid_node(0, tree->total_non_trie - 1, liste_triee, tree->racine_trie);
tree->total_trie += tree->total_non_trie;
tree->total_non_trie = 0;
tree->racine_non_trie = NULL;
free(liste_triee);
}


static void resort_partie_non_triee(index_bt *tree)
/* retrier entierement la partie racine_non_trie de tree
*/
{
bt_node **liste_triee;

if(tree->total_non_trie <= 2) return;
liste_triee = (bt_node **)malloc(tree->total_non_trie * sizeof(bt_node *));
if(liste_triee == NULL) return;
tri_non_trie(tree->racine_non_trie, liste_triee);
if( liste_triee[tree->total_non_trie - 1] == NULL ) {
	fprintf(stderr, "ERROR: repetition in file %s\n", tree->k->filename);
	exit(ERREUR);
	}
tree->racine_non_trie = opt_tree(0, tree->total_non_trie - 1, liste_triee);
free(liste_triee);
}


static bt_node *opt_tree(int bas, int haut, bt_node **tableau)
/* creer les liens entre les noeuds tries dans tableau pour arbre optimal
*/
{
bt_node *noeud;
int milieu;

milieu = (bas + haut)/2;
noeud = tableau[milieu];
noeud->before = noeud->after = NULL;
if(bas <= milieu - 1) {
	noeud->before = opt_tree(bas, milieu-1, tableau);
	}
if(milieu + 1 <= haut) {
	noeud->after = opt_tree(milieu+1, haut, tableau);
	}
return noeud;
}


static bt_node **tri_non_trie(bt_node *racine, bt_node **liste_triee)
/* descendre arbre binaire debutant en racine selon ordre alphab
et mettre pointeur vers chaque noeud trouve dans tableau liste_triee
*/
{
if(racine == NULL) return liste_triee;
liste_triee = tri_non_trie(racine->before, liste_triee);
*(liste_triee++) = racine;
liste_triee = tri_non_trie(racine->after, liste_triee);
return liste_triee;
}


static void add_mid_node(int bas, int haut, bt_node **liste_triee, 
	bt_node *racine)
/* ajoute a l'arbre binaire racine les noeuds de rang compris entre bas et haut
dans le tableau liste_triee, ceci de maniere dichotomique recursive
*/
{
int milieu;
bt_node *new_racine;

milieu = (bas + haut)/2;
liste_triee[milieu]->before = liste_triee[milieu]->after = NULL;
add_node( liste_triee[milieu], racine);
if(bas <= milieu - 1) {
	new_racine = racine;
	while( new_racine->before != NULL && 
		strcmp(liste_triee[milieu - 1]->name, new_racine->name) < 0 )
		new_racine = new_racine->before;
	add_mid_node(bas, milieu - 1, liste_triee, new_racine );
	}
if(milieu + 1 <= haut) {
	new_racine = racine;
	while( new_racine->after != NULL && 
		strcmp(liste_triee[milieu + 1]->name, new_racine->name) > 0 )
		new_racine = new_racine->after;
	add_mid_node(milieu + 1, haut, liste_triee, new_racine );
	}
}


static void add_node(bt_node *noeud, bt_node *racine)
/* ajouter *noeud dans l'arbre binaire qui debute en racine 
*/
{
int test;

test = strcmp(noeud->name, racine->name);
if(test > 0) {
	if(racine->after == NULL) racine->after = noeud;
	else	add_node(noeud, racine->after);
	}
else if(test < 0)	{
	if(racine->before == NULL) racine->before = noeud;
	else	add_node(noeud, racine->before);
	}
}


static int trim_key_max(char *record, int max_width)
{
char *p;
p=record+max_width;
while(p>record && *(p-1)==' ') p--; 
*p=0;
return p-record;
}

/* pour le tester *
int find_key2(char *name, index_bt *tree, int *prof);
bt_node *find_node2(char *name, bt_node *noeud, int *prof);

static void test(index_bt *tree)
{
int tot, num, trouve, prof, maxprof, avprof, count;
char croix[60];

memset(croix, 'x', tree->char_width);
tot=read_first_rec(tree->k,NULL);
maxprof = 0; avprof = 0, count = 0;
printf("debut test\n");fflush(stdout);
for(num=2; num<=tot; num++) {
	if( dir_read(tree->k,num,1,tree->buffer) != 1) dir_readerr(tree->k,num);
	if(memcmp(tree->buffer, croix, tree->char_width)==0) continue;
	trim_key_max(tree->buffer, tree->char_width);
	prof = 0 ;
	trouve = find_key2(tree->buffer, tree, &prof);
	avprof += prof; count++;
	if(prof>maxprof) maxprof = prof;
	if(trouve != num) {
		printf("pb %s\n", tree->buffer);
		}
	}
printf("fin test\n");fflush(stdout);
printf("maxprof=%d   aver. prof=%d\n", maxprof, avprof/count );
}


static void test_aux_val(bt_node *racine, index_bt *tree)
{
static int trouve, *paux;

if(racine == NULL) return;
trouve = find_key(racine->name, tree, FALSE, &paux);
if(trouve != racine->recnum) 
	printf("pas trouve %s\n", racine->name);
else if(paux != &( ((bt_node_aux *)racine)->aux_data ) )
	printf("pb addr aux_data %s\n", racine->name);
if( ((bt_node_aux *)racine)->aux_data != 0) {
	printf("pb aux_data %d  %s\n", 
		((bt_node_aux *)racine)->aux_data, racine->name);
	}
test_aux_val(racine->after, tree);
test_aux_val(racine->before, tree);
}


static void test_aux(index_bt *tree)
{
printf("debut test_aux\n");fflush(stdout);
test_aux_val(tree->racine_trie, tree);
test_aux_val(tree->racine_non_trie, tree);
printf("fin test_aux\n");fflush(stdout);
}

#include <math.h>
main()
{
index_bt *tree;
char rep[10];

dir_acnucopen("RO");

printf("\ntotal acc=%d\n",read_first_rec(kacc,NULL));
tree = load_index_bt(kacc, 8, FALSE);
printf("arbre acc=%p\n",tree); fflush(stdout);
if(tree != NULL) 
	printf("total trie=%d non trie=%d log2=%f\n", tree->total_trie, 
		tree->total_non_trie, log(tree->total_trie)/log(2) ); 
test(tree);
free_index_bt(tree);

printf("\ntotal bib=%d\n",read_first_rec(kbib,NULL));
tree = load_index_bt(kbib, 40, FALSE);
printf("arbre bib=%p\n",tree); fflush(stdout);
test(tree);
free_index_bt(tree);

printf("\ntotal aut=%d\n",read_first_rec(kaut,NULL));
tree = load_index_bt(kaut, 20, FALSE);
printf("arbre aut=%p\n",tree); fflush(stdout);
test(tree);
free_index_bt(tree);

printf("\ntotal smj=%d\n",read_first_rec(ksmj,NULL));
tree = load_index_bt(ksmj, 20, FALSE);
printf("arbre smj=%p\n",tree); fflush(stdout);
test(tree);
free_index_bt(tree);

printf("\ntotal spec=%d\n",read_first_rec(kspec,NULL)); fflush(stdout);
tree = load_index_bt(kspec, 40, TRUE);
printf("arbre spec=%p\n",tree); fflush(stdout);
if(tree != NULL) {
	printf("total trie=%d non trie=%d log2=%f\n", tree->total_trie, 
		tree->total_non_trie, log(tree->total_trie)/log(2) ); 
	fflush(stdout);
	test(tree);
	test_aux(tree);
	free_index_bt(tree);
	}

printf("attendre..."); gets(rep);
}

bt_node *find_node2(char *name, bt_node *noeud, int *prof)
{
int test;

if(noeud == NULL) {
	return NULL;
	}
test = strcmp(name, noeud->name);
if(test == 0) return noeud;
if(test > 0) {
	if(noeud->after == NULL) {
		return NULL;
		}
	else	{
		++(*prof);
		return find_node2(name, noeud->after, prof);
		}
	}
else	{
	if(noeud->before == NULL) {
		return NULL;
		}
	else	{
		++(*prof);
		return find_node2(name, noeud->before, prof);
		}
	}
}

int find_key2(char *name, index_bt *tree, int *prof)
{
bt_node *noeud;

noeud = find_node2(name, tree->racine_trie, prof);
if(noeud != NULL) return noeud->recnum;
noeud = find_node2(name, tree->racine_non_trie, prof);
if(noeud == NULL) return 0;
else return noeud->recnum;
}


* fin du test */
