CC = gcc
OPT = -O
CFLAGS  = $(OPT) -Dunix


#useful for SGI !!!
#SHELL = /usr/bin/sh

F77 = f77
FFLAGS = -O -w

LIBNAME = cacnuc
LIBFULLNAME = lib$(LIBNAME).a

LIB_OBJS = utilquery.o utilgener.o utilgener2.o arbrebin.o gestion_acnuc.o \
	use_acnuc.o dynlist.o lngbit.o mdshrt_lng.o bit.o dir_io.o hashacc.o \
	prep_acnuc_requete.o two_banks.o proc_requete.o misc_acnuc.o


#each of these is a possible target of this makefile
CPROGS = gbemgener swgener updatehelp newordalphab supold \
        connectindex conv_to_bigannots testmatchindex acnuc2fasta \
        initf sortsubseq compressnewdiv renamediv smjytload

QUERYFILES = query.o extract.o pretty_seq.o


all :  query cprogs

cprogs : $(CPROGS)

redoall : clean all

$(LIBFULLNAME) : $(LIB_OBJS)
	ar rv $(LIBFULLNAME) $(LIB_OBJS)
	@if [ `uname` = Darwin ] ; \
	then set -x; ranlib $(LIBFULLNAME) ; \
	fi

query : $(LIBFULLNAME) $(QUERYFILES)
	$(CC) -o $@ $(QUERYFILES) -L. -l$(LIBNAME)

$(CPROGS) : $(LIBFULLNAME) 
	$(CC) -o $@ $(CFLAGS) $@.c -L. -l$(LIBNAME)

querydiv : querydiv.o coperations.o
	$(F77) -o $@ $@.o coperations.o

clean :
	 rm -f $(LIB_OBJS) $(CPROGS) $(LIBFULLNAME) \
		$(QUERYFILES) query querydiv
	@set -x; for objname in $(CPROGS); do rm -f $$objname.o; done


## Implicit actions
.c.o :	
	@if [ `uname` = OSF1 ] ; \
	then set -x; c89 -c $(CFLAGS) $< ; \
	elif [ `uname` = HP-UX ] ; \
	then set -x; cc -Ae -c $(CFLAGS) $< ; \
	else set -x; $(CC) -c $(CFLAGS) $< ; \
	fi
.f.o :  
	$(F77) -c $(FFLAGS) $<
