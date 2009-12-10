CC=			gcc.exe
AR=			ar.exe
CFLAGS=		-g -Wall -O2
DFLAGS=		-D_CURSES_LIB=2 -D_USE_KNETFILE
KNETFILE_O=	knetfile.o
LOBJS=		bgzf.o kstring.o bam_aux.o bam.o bam_import.o sam.o bam_index.o	\
			bam_pileup.o bam_lpileup.o bam_md.o glf.o razf.o faidx.o bam_sort.o \
			$(KNETFILE_O)
AOBJS=		bam_tview.o bam_maqcns.o bam_plcmd.o sam_view.o	\
			bam_rmdup.o bam_rmdupse.o bam_mate.o bam_stat.o bam_color.o	\
			bamtk.o kaln.o sam_header.o
PROG=		samtools
INCLUDES=	-Iwin32
SUBDIRS=	.
LIBPATH=	

.SUFFIXES:.c .o

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

.PHONY:all lib clean cleanlocal
.PHONY:all-recur lib-recur clean-recur cleanlocal-recur install-recur

lib:libbam.a

libbam.a:$(LOBJS)
		$(AR) -cru $@ $(LOBJS)

samtools:$(AOBJS) libbam.a
		$(CC) $(CFLAGS) -o $@ $(AOBJS) $(LIBPATH) -lm -L. -lbam -Lwin32 -lz -lcurses -lws2_32

razip:razip.o razf.o $(KNETFILE_O)
		$(CC) $(CFLAGS) -o $@ razf.o razip.o $(KNETFILE_O) -lz

bgzip:bgzip.o bgzf.o $(KNETFILE_O)
		$(CC) $(CFLAGS) -o $@ bgzf.o bgzip.o $(KNETFILE_O) -lz

razip.o:razf.h
bam.o:bam.h razf.h bam_endian.h kstring.h
sam.o:sam.h bam.h
bam_import.o:bam.h kseq.h khash.h razf.h
bam_pileup.o:bam.h razf.h ksort.h
bam_plcmd.o:bam.h faidx.h bam_maqcns.h glf.h
bam_index.o:bam.h khash.h ksort.h razf.h bam_endian.h
bam_lpileup.o:bam.h ksort.h
bam_tview.o:bam.h faidx.h bam_maqcns.h
bam_maqcns.o:bam.h ksort.h bam_maqcns.h
bam_sort.o:bam.h ksort.h razf.h
bam_md.o:bam.h faidx.h
glf.o:glf.h

faidx.o:faidx.h razf.h khash.h
faidx_main.o:faidx.h razf.h

clean:
		rm -fr gmon.out *.o *.exe *.dSYM razip bgzip $(PROG) *~ *.a
