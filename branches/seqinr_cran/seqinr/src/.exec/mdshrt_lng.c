#include "dir_acnuc.h"

/* fonctions contenues */
int addshrt(int point, int val);
int supshrt(int point, int val);
int mdshrt(DIR_FILE *kan, int nrec, int offset, int val, int *newplist);
int addlng(int point, int val);
int suplng(int point, int val);
int mdlng(DIR_FILE *kan, int nrec, int offset, int val, int *newplist);


int addshrt(int point, int val)
{
unsigned last, next;

readshrt(1);
last=pshrt->val;
next=point;
if((unsigned)point <= last) {
	while(next) {
		point=next;
		readshrt(point);
		if(pshrt->val==val)return 2;
		next=pshrt->next;
		}
	pshrt->next=last+1;
	writeshrt(point);
	}
pshrt->val=val;
pshrt->next=0;
writeshrt(++last);
pshrt->val=last; pshrt->next=0;
writeshrt(1);
return 1;
}

int supshrt(int point, int val)
{
int pre, next;
if(!point) return 3;
pre=point;
readshrt(point);
if(pshrt->val != val) {
   	do	{
		if(!pshrt->next)return 3;
		pre=point;
		point=pshrt->next;
		readshrt(point);
		next=pshrt->next;
		}
   	while(pshrt->val != val);
	readshrt(pre);
	pshrt->next=next;
   	}
else 	{
	if( pshrt->next==0) return 2;
	readshrt(pshrt->next);
	}
writeshrt(pre);
return 1;
}



int mdshrt(DIR_FILE *kan, int nrec, int offset, int val, int *newplist)
{
int retval=1, tot, *point, ier, new=0, next;
char *buffer;
if(kan==kloc) {
	point=(int *)ploc;
	tot=13;
	buffer=(char *)ploc;
	}
else if(kan==ksub) {
	tot=7;
	point= &(psub->length);
	buffer=(char *)psub;
	}
else if(kan==kbib) {
	tot=4;
	point= &(pbib->plsub);
	buffer=(char *)pbib;
	}
else if(kan==kacc) {
	tot=1;
	point= &(pacc->plsub);
	buffer=(char *)pacc;
	}
else if(kan==kaut) {
	tot=2;
	point= &(paut->plref);
	buffer=(char *)paut;
	}
else if(kan==kspec) {
	tot=6;
	point= &(pspec->libel);
	buffer=(char *)pspec;
	}
else if(kan==kkey) {
	tot=5;
	point= &(pkey->libel);
	buffer=(char *)pkey;
	}
else
	return 3;
if(offset==0 || abs(offset)>tot) return 3;
if(kan!=kacc)
	dir_read(kan,nrec,1,buffer);
else
	readacc(nrec);
point += abs(offset)-1;
next= *point;
if(offset>0) {
	if(next==0) {
		readshrt(1);
		next= (unsigned)pshrt->val + 1;
		new=1;
		}
	if(addshrt(next,val)==2)retval=2;
	}
else 	{
	ier=supshrt(next,val);
	if(ier==2) {
		next=0;
		new=1;
		}
	else if(ier==3) {
		retval=2;
		}
	}
if(new)	{
	*point=next;
	if(kan!=kacc) {
		if(dir_write(kan,nrec,1,buffer))dir_writeerr(kan,nrec);
		}
	else
		writeacc(nrec);
	}
if(newplist!=NULL) *newplist=next;
return retval;
}

int addlng(int point, int val)
{
int last, next, loop;
if(val==0)return 2;
readlng(1);
last=plng->sub[0];
next=point;
loop=0;
if(next<=last) {
	while(next) {
		point=next;
		readlng(point);
		for(loop=0; loop<SUBINLNG; loop++) {
			if(plng->sub[loop]==val) return 2;
			}
		next=plng->next;
		}
	if(plng->sub[SUBINLNG-1]!=0) {
		plng->next= ++last;
		write_first_rec(klng,last,0);
		writelng(point);
		loop=0;
		next=last;
		}
	else 	{
		for(loop=SUBINLNG-2; loop>=0; loop--) if(plng->sub[loop]) break;
		++loop;
		next=point;
		}
	}
else	{
	write_first_rec(klng,next,0);
	}
plng->sub[loop]=val;
while(++loop<SUBINLNG) plng->sub[loop]=0;
plng->next=0;
writelng(next);
return 1;
}

int suplng(int point, int val)
{
int pre=0, last, found=0, loop1, loop2, new=0;
if(!point) return 3;
do	{
	readlng(point);
	loop1= -1;
	while(++loop1<SUBINLNG) {
		if(plng->sub[loop1]==val) {
			found=1; 
			break; 
			}
		}
	if(!found) {
		if(plng->next==0)return 3;
		pre=point;
		point=plng->next;
		}
	}
while(!found);
last=point;
while(plng->next) {
	pre=last;
	last=plng->next;
	readlng(last);
	}	
loop2= 0;
while(++loop2<SUBINLNG) if(plng->sub[loop2]==0) break;
loop2--;
new=plng->sub[loop2];
if(loop2>0) {
	plng->sub[loop2]=0;
	writelng(last);
	}
else	{
	if(pre==0) return 2;
	readlng(pre);
	plng->next=0;
	writelng(pre);
	}
if(new!=val) {
	readlng(point);
	plng->sub[loop1]=new;
	writelng(point);
	}
return 1;
}



int mdlng(DIR_FILE *kan, int nrec, int offset, int val, int *newplist)
{
int retval=1, tot, *point, ier, new=FALSE, next, cas_sp_kw=FALSE, desc, premier;
char *buffer;
if(kan==ksub) {
	tot=7;
	point= &(psub->length);
	buffer=(char *)psub;
	}
else if(kan==ksmj) {
	tot=2;
	point= &(psmj->plong);
	buffer=(char *)psmj;
	}
else if(kan==kspec) {
	tot=6;
	point= &(pspec->libel);
	buffer=(char *)pspec;
	cas_sp_kw=TRUE;
	}
else if(kan==kkey) {
	tot=5;
	point= &(pkey->libel);
	buffer=(char *)pkey;
	cas_sp_kw=TRUE;
	}
else
	return 3;
if(offset==0 || abs(offset)>tot) return 3;
dir_read(kan,nrec,1,buffer);
point += abs(offset)-1;
next= *point;
if(kan==ksub && abs(offset)==3) next=abs(next);
if(offset>0) {
	if(next==0) {
		readlng(1);
		next = plng->sub[0] + 1;
		new = TRUE;
		}
	if(addlng(next,val)==2)retval=2;
	}
else 	{
	ier=suplng(next,val);
	if(ier==2) {
		next = 0;
		new = TRUE;
		}
	else if(ier==3) {
		retval=2;
		}
	}
if(new)	{
	if(kan==ksub && abs(offset)==3)next = -next;
	*point=next;
	if(dir_write(kan,nrec,1,buffer)) dir_writeerr(kan,nrec);
	}
if(newplist != NULL) *newplist = next;
if(cas_sp_kw && new) { /* controle du signe debut liste desc */
	if(kan == kspec) {
		ier = (pspec->plsub != 0 || pspec->plhost != 0); 
		desc = pspec->desc;
		}
	else	{
		ier = (pkey->plsub != 0);
		desc = pkey->desc;
		}
	readshrt(desc);
	premier = pshrt->val;
	if( ier && premier > 0 ) { /* liste non vide: mettre signe < 0 */
		pshrt->val = - premier;
		writeshrt(desc);
		}
	else if( (!ier) && premier < 0 ) { /* liste vide: mettre signe > 0 */
		pshrt->val = abs(premier);
		writeshrt(desc);
		}
	}
return retval;
}
