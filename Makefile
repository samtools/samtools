# The default version string in bam.h and bcftools/bcf.h can be overriden directly
#   make VERSION="-DVERSION='\\\"my-version\\\"'"
# or using the git-stamp rule
#   make git-stamp
VERSION=

# Adjust $(HTSDIR) to point to your top-level htslib directory
HTSDIR = ../htslib
HTSLIB = $(HTSDIR)/htslib/libhts.a

CC=			gcc
CFLAGS=		-g -Wall $(VERSION) -O2
#LDFLAGS=		-Wl,-rpath,\$$ORIGIN/../lib
DFLAGS=		-D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -D_USE_KNETFILE -D_CURSES_LIB=1
LOBJS=		bam_aux.o bam.o bam_import.o sam.o bam_index.o	\
			bam_pileup.o bam_lpileup.o sam_header.o
AOBJS=		bam_tview.o bam_plcmd.o sam_view.o \
			bam_cat.o bam_md.o bam_reheader.o bam_sort.o bedidx.o kprobaln.o \
			bam_rmdup.o bam_rmdupse.o bam_mate.o bam_stat.o bam_color.o \
			bamtk.o kaln.o bam2bcf.o bam2bcf_indel.o errmod.o sample.o \
			cut_target.o phase.o bam2depth.o padding.o bedcov.o bamshuf.o \
			bam_tview_curses.o bam_tview_html.o
PROG=		samtools
INCLUDES=	-I. -I$(HTSDIR)
SUBDIRS=	. bcftools misc
LIBPATH=
LIBCURSES=	-lcurses # -lXCurses


.SUFFIXES:.c .o
.PHONY: all lib

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

all-recur lib-recur clean-recur cleanlocal-recur install-recur:
		@target=`echo $@ | sed s/-recur//`; \
		wdir=`pwd`; \
		list='$(SUBDIRS)'; for subdir in $$list; do \
			cd $$subdir; \
			$(MAKE) CC="$(CC)" DFLAGS="$(DFLAGS)" CFLAGS="$(CFLAGS)" \
				HTSDIR="$(HTSDIR)" HTSLIB="$(HTSLIB)" \
				INCLUDES="$(INCLUDES)" LIBPATH="$(LIBPATH)" $$target || exit 1; \
			cd $$wdir; \
		done;

all:$(PROG)

git-stamp:
		make VERSION="-DVERSION='\\\"`git describe --always --dirty`\\\"'"

.PHONY:all lib clean cleanlocal
.PHONY:all-recur lib-recur clean-recur cleanlocal-recur install-recur

lib:libbam.a

libbam.a:$(LOBJS)
		$(AR) -csru $@ $(LOBJS)

samtools:lib-recur $(AOBJS) $(HTSLIB)
		$(CC) $(CFLAGS) -o $@ $(AOBJS) $(LDFLAGS) libbam.a -Lbcftools -lbcf $(LIBPATH) $(HTSLIB) $(LIBCURSES) -lm -lz -lpthread

razip:razip.o $(HTSLIB)
		$(CC) $(CFLAGS) -o $@ razip.o $(HTSLIB) -lz

bgzip:bgzip.o $(HTSLIB)
		$(CC) $(CFLAGS) -o $@ bgzip.o $(HTSLIB) -lz -lpthread

razip.o:$(HTSDIR)/htslib/razf.h
bam.o:bam.h bam_endian.h $(HTSDIR)/htslib/kstring.h sam_header.h
sam.o:sam.h bam.h
bam_import.o:bam.h $(HTSDIR)/htslib/kseq.h $(HTSDIR)/htslib/khash.h
bam_pileup.o:bam.h $(HTSDIR)/htslib/ksort.h
bam_plcmd.o:bam.h $(HTSDIR)/htslib/faidx.h bcftools/bcf.h bam2bcf.h
bam_index.o:bam.h $(HTSDIR)/htslib/khash.h $(HTSDIR)/htslib/ksort.h bam_endian.h
bam_lpileup.o:bam.h $(HTSDIR)/htslib/ksort.h
bam_tview.o:bam.h $(HTSDIR)/htslib/faidx.h bam_tview.h
bam_tview_curses.o:bam.h $(HTSDIR)/htslib/faidx.h bam_tview.h
bam_tview_html.o:bam.h $(HTSDIR)/htslib/faidx.h bam_tview.h
bam_sort.o:bam.h $(HTSDIR)/htslib/ksort.h
bam_md.o:bam.h $(HTSDIR)/htslib/faidx.h
sam_header.o:sam_header.h $(HTSDIR)/htslib/khash.h
bcf.o:bcftools/bcf.h
bam2bcf.o:bam2bcf.h errmod.h bcftools/bcf.h
bam2bcf_indel.o:bam2bcf.h
errmod.o:errmod.h
phase.o:bam.h $(HTSDIR)/htslib/khash.h $(HTSDIR)/htslib/ksort.h
bamtk.o:bam.h



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
