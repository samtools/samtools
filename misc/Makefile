CC=			gcc
CXX=		g++
CFLAGS=		-g -Wall -O2 #-m64 #-arch ppc
CXXFLAGS=	$(CFLAGS)
DFLAGS=		-D_FILE_OFFSET_BITS=64
OBJS=		
PROG=		md5sum-lite md5fa maq2sam-short maq2sam-long wgsim seqtk
INCLUDES=	-I..
SUBDIRS=	.

.SUFFIXES:.c .o

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

lib-recur all-recur clean-recur cleanlocal-recur install-recur:
		@target=`echo $@ | sed s/-recur//`; \
		wdir=`pwd`; \
		list='$(SUBDIRS)'; for subdir in $$list; do \
			cd $$subdir; \
			$(MAKE) CC="$(CC)" DFLAGS="$(DFLAGS)" CFLAGS="$(CFLAGS)" \
				INCLUDES="$(INCLUDES)" $$target || exit 1; \
			cd $$wdir; \
		done;

lib:

seqtk:seqtk.o
		$(CC) $(CFLAGS) -o $@ seqtk.o -lm -lz

wgsim:wgsim.o
		$(CC) $(CFLAGS) -o $@ wgsim.o -lm -lz

md5fa:md5.o md5fa.o md5.h ../kseq.h
		$(CC) $(CFLAGS) -o $@ md5.o md5fa.o -lz

md5sum-lite:md5sum-lite.o
		$(CC) $(CFLAGS) -o $@ md5sum-lite.o

md5sum-lite.o:md5.c md5.h
		$(CC) -c $(CFLAGS) -DMD5SUM_MAIN -o $@ md5.c

maq2sam-short:maq2sam.c
		$(CC) $(CFLAGS) -o $@ maq2sam.c -lz

maq2sam-long:maq2sam.c
		$(CC) $(CFLAGS) -DMAQ_LONGREADS -o $@ maq2sam.c -lz

md5fa.o:md5.h md5fa.c
		$(CC) $(CFLAGS) -c -I.. -o $@ md5fa.c

seqtk.o:seqtk.c ../khash.h ../kseq.h
		$(CC) $(CFLAGS) -c -I.. -o $@ seqtk.c

wgsim.o:wgsim.c ../kseq.h
		$(CC) $(CFLAGS) -c -I.. -o $@ wgsim.c

cleanlocal:
		rm -fr gmon.out *.o a.out *.exe *.dSYM $(PROG) *~ *.a

clean:cleanlocal-recur
