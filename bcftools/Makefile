CC=			gcc
CFLAGS=		-g -Wall -O2 #-m64 #-arch ppc
DFLAGS=		-D_FILE_OFFSET_BITS=64 -D_USE_KNETFILE
LOBJS=		bcf.o vcf.o bcfutils.o prob1.o em.o kfunc.o kmin.o index.o fet.o mut.o bcf2qcall.o
OMISC=		..
AOBJS=		call1.o main.o $(OMISC)/kstring.o $(OMISC)/bgzf.o $(OMISC)/knetfile.o $(OMISC)/bedidx.o
PROG=		bcftools
INCLUDES=	
SUBDIRS=	.

.SUFFIXES:.c .o

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) -I.. $(INCLUDES) $< -o $@

all-recur lib-recur clean-recur cleanlocal-recur install-recur:
		@target=`echo $@ | sed s/-recur//`; \
		wdir=`pwd`; \
		list='$(SUBDIRS)'; for subdir in $$list; do \
			cd $$subdir; \
			$(MAKE) CC="$(CC)" DFLAGS="$(DFLAGS)" CFLAGS="$(CFLAGS)" \
				INCLUDES="$(INCLUDES)" LIBPATH="$(LIBPATH)" $$target || exit 1; \
			cd $$wdir; \
		done;

all:$(PROG)

lib:libbcf.a

libbcf.a:$(LOBJS)
		$(AR) -csru $@ $(LOBJS)

bcftools:lib $(AOBJS)
		$(CC) $(CFLAGS) -o $@ $(AOBJS) -L. $(LIBPATH) -lbcf -lm -lz

bcf.o:bcf.h
vcf.o:bcf.h
index.o:bcf.h
bcfutils.o:bcf.h
prob1.o:prob1.h bcf.h
call1.o:prob1.h bcf.h
bcf2qcall.o:bcf.h
main.o:bcf.h

bcf.pdf:bcf.tex
		pdflatex bcf

cleanlocal:
		rm -fr gmon.out *.o a.out *.dSYM $(PROG) *~ *.a bcf.aux bcf.log bcf.pdf *.class libbcf.*.dylib libbcf.so*

clean:cleanlocal-recur
