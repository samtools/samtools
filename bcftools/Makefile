CC=			gcc
CFLAGS=		-g -Wall -O2 #-m64 #-arch ppc
DFLAGS=		-D_FILE_OFFSET_BITS=64 -D_USE_KNETFILE -D_BCF_QUAD
LOBJS=		bcf.o vcf.o bcfutils.o prob1.o kfunc.o index.o fet.o
OMISC=		..
AOBJS=		call1.o main.o $(OMISC)/kstring.o $(OMISC)/bgzf.o $(OMISC)/knetfile.o
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
		$(AR) -cru $@ $(LOBJS)

bcftools:lib $(AOBJS)
		$(CC) $(CFLAGS) -o $@ $(AOBJS) -lm $(LIBPATH) -lz -L. -lbcf

bcf.o:bcf.h
prob1.o:prob1.h bcf.h
vcfout.o:bcf.h prob1.h
main.o:bcf.h

bcf.pdf:bcf.tex
		pdflatex bcf

cleanlocal:
		rm -fr gmon.out *.o a.out *.dSYM $(PROG) *~ *.a bcf.aux bcf.log bcf.pdf *.class libbcf.*.dylib libbcf.so*

clean:cleanlocal-recur
