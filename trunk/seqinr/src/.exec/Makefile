CC = gcc
OPT = -O
CFLAGS  = $(OPT) -Dunix
SRCDIR = /bge/charif/seqinr/src/.exec/
pkg = seqinr


#useful for SGI !!!
#SHELL = /usr/bin/sh

F77 = f77
FFLAGS = -O -w

LIBNAME = cacnuc
LIBFULLNAME = lib$(LIBNAME).a

LIB_OBJS = utilquery.o utilgener.o utilgener2.o arbrebin.o gestion_acnuc.o \
	use_acnuc.o dynlist.o lngbit.o mdshrt_lng.o bit.o dir_io.o hashacc.o \
	prep_acnuc_requete.o two_banks.o timing.o proc_requete.o

LIB_OBJS_EXEC = $(LIB_OBJS:%=$(SRCDIR)%)

SHARED_OBJ = R_acnuc.c

#each of these is a possible target of this makefile
CPROGS = gbemgener swgener updatehelp newordalphab supold \
        connectindex conv_to_bigannots testmatchindex acnuc2fasta \
        initf sortsubseq compressnewdiv renamediv smjytload

QUERYFILES = query.o extract.o pretty_seq.o

QUERYFILES_EXEC= $(QUERYFILES:%=$(SRCDIR)%)


all :  query cprogs

cprogs : $(CPROGS)

redoall : clean all

$(LIBFULLNAME) : $(LIB_OBJS_EXEC)
	cd $(SRCDIR) ; 	ar rv $(LIBFULLNAME) $(LIB_OBJS)
	@if [ `uname` = Darwin ] ; \
	then set -x; ranlib $(LIBFULLNAME) ; \
	fi

query : $(LIBFULLNAME) $(QUERYFILES_EXEC)
	cd $(SRCDIR); $(CC) -o $@ $(QUERYFILES) -L. -l$(LIBNAME)

$(CPROGS) : $(LIBFULLNAME) 
	cd $(SRCDIR); $(CC) -o $@ $(CFLAGS) $@.c -L. -l$(LIBNAME)

querydiv : querydiv.o coperations.o
	$(F77) -o $@ $@.o coperations.o

clean :
	 rm -f $(LIB_OBJS) $(CPROGS) $(LIBFULLNAME) \
		$(QUERYFILES) query querydiv
	@set -x; for objname in $(CPROGS); do rm -f $$objname.o; done


## Implicit actions
.c.o :	
	@if [ `uname` = OSF1 ] ; \
	then cd $(SRCDIR) ; set -x; c89 -c $(CFLAGS) $(<F) ; \
	elif [ `uname` = HP-UX ] ; \
	then cd $(SRCDIR) ; set -x; cc -Ae -c $(CFLAGS) $(<F) ; \
	else cd $(SRCDIR) ; set -x; $(CC) -c $(CFLAGS) $(<F) ; \
	fi
.f.o :  
	$(F77) -c $(FFLAGS) $<
