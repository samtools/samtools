CC=			gcc
CFLAGS=		-g -Wall -O2 #-m64 #-arch ppc
DFLAGS=		-D_FILE_OFFSET_BITS=64 -D_USE_KNETFILE -D_CURSES_LIB=1
KNETFILE_O=	knetfile.o
LOBJS=		bgzf.o kstring.o bam_aux.o bam.o bam_import.o sam.o bam_index.o	\
			bam_pileup.o bam_lpileup.o bam_md.o glf.o razf.o faidx.o \
			$(KNETFILE_O) bam_sort.o sam_header.o bam_reheader.o
AOBJS=		bam_tview.o bam_maqcns.o bam_plcmd.o sam_view.o	\
			bam_rmdup.o bam_rmdupse.o bam_mate.o bam_stat.o bam_color.o	\
			bamtk.o kaln.o bam2bcf.o errmod.o sample.o
PROG=		samtools
INCLUDES=	-I.
SUBDIRS=	. bcftools misc
LIBPATH=
LIBCURSES=	-lcurses # -lXCurses

.SUFFIXES:.c .o

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

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

.PHONY:all lib clean cleanlocal
.PHONY:all-recur lib-recur clean-recur cleanlocal-recur install-recur

lib:libbam.a

libbam.a:$(LOBJS)
		$(AR) -cru $@ $(LOBJS)

samtools:lib-recur $(AOBJS)
		$(CC) $(CFLAGS) -o $@ $(AOBJS) libbam.a -lm $(LIBPATH) $(LIBCURSES) -lz -Lbcftools -lbcf

razip:razip.o razf.o $(KNETFILE_O)
		$(CC) $(CFLAGS) -o $@ razf.o razip.o $(KNETFILE_O) -lz

bgzip:bgzip.o bgzf.o $(KNETFILE_O)
		$(CC) $(CFLAGS) -o $@ bgzf.o bgzip.o $(KNETFILE_O) -lz

razip.o:razf.h
bam.o:bam.h razf.h bam_endian.h kstring.h sam_header.h
sam.o:sam.h bam.h
bam_import.o:bam.h kseq.h khash.h razf.h
bam_pileup.o:bam.h razf.h ksort.h
bam_plcmd.o:bam.h faidx.h bam_maqcns.h glf.h bcftools/bcf.h bam2bcf.h
bam_index.o:bam.h khash.h ksort.h razf.h bam_endian.h
bam_lpileup.o:bam.h ksort.h
bam_tview.o:bam.h faidx.h bam_maqcns.h
bam_maqcns.o:bam.h ksort.h bam_maqcns.h kaln.h
bam_sort.o:bam.h ksort.h razf.h
bam_md.o:bam.h faidx.h
glf.o:glf.h
sam_header.o:sam_header.h khash.h
bcf.o:bcftools/bcf.h
bam2bcf.o:bam2bcf.h errmod.h bcftools/bcf.h
errmod.o:errmod.h

faidx.o:faidx.h razf.h khash.h
faidx_main.o:faidx.h razf.h


libbam.1.dylib-local:$(LOBJS)
		libtool -dynamic $(LOBJS) -o libbam.1.dylib -lc -lz

libbam.so.1-local:$(LOBJS)
		$(CC) -shared -Wl,-soname,libbam.so -o libbam.so.1 $(LOBJS) -lc -lz

dylib:
		@$(MAKE) cleanlocal; \
		case `uname` in \
			Linux) $(MAKE) CFLAGS="$(CFLAGS) -fPIC" libbam.so.1-local;; \
			Darwin) $(MAKE) CFLAGS="$(CFLAGS) -fPIC" libbam.1.dylib-local;; \
			*) echo 'Unknown OS';; \
		esac


cleanlocal:
		rm -fr gmon.out *.o a.out *.exe *.dSYM razip bgzip $(PROG) *~ *.a *.so.* *.so *.dylib

clean:cleanlocal-recur
