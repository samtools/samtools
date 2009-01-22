CC=			gcc
CXX=		g++
CFLAGS=		-g -Wall -O2 -m64 #-arch ppc
CXXFLAGS=	$(CFLAGS)
DFLAGS=		#-D_FILE_OFFSET_BITS=64
OBJS=		
PROG=		faidx md5sum-lite md5fa maq2sam-short maq2sam-long wgsim
INCLUDES=
LIBS=		-lm -lz
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
			$(MAKE) CC="$(CC)" CXX="$(CXX)" DFLAGS="$(DFLAGS)" CFLAGS="$(CFLAGS)" \
				INCLUDES="$(INCLUDES)" $$target || exit 1; \
			cd $$wdir; \
		done;

lib:

wgsim:wgsim.o
		$(CC) $(CFLAGS) -o $@ wgsim.o -lm

faidx:../faidx.c ../faidx.h
		$(CC) $(CFLAGS) -D_NO_RAZF -DFAIDX_MAIN -o $@ ../faidx.c

md5fa:md5.o md5fa.o md5.h ../kseq.h
		$(CC) $(CFLAGS) -o $@ md5.o md5fa.o -lz

md5sum-lite:md5.c md5.h
		$(CC) $(CFLAGS) -DMD5SUM_MAIN -o $@ md5.c

maq2sam-short:maq2sam.c
		$(CC) $(CFLAGS) -o $@ maq2sam.c -lz

maq2sam-long:maq2sam.c
		$(CC) $(CFLAGS) -DMAQ_LONGREADS -o $@ maq2sam.c -lz

md5fa.o:md5.h md5fa.c
		$(CC) $(CFLAGS) -c -I.. -o $@ md5fa.c

cleanlocal:
		rm -fr gmon.out *.o a.out *.dSYM $(PROG) *~ *.a

clean:cleanlocal-recur
