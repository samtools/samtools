# The default version string in bam.h and bcftools/bcf.h can be overriden directly
#   make VERSION="-DVERSION='\\\"my-version\\\"'"
# or using the git-stamp rule
#   make git-stamp
VERSION=

all: samtools bgzip razip

# Adjust $(HTSDIR) to point to your top-level htslib directory
HTSDIR = ../htslib
include $(HTSDIR)/htslib.mk
HTSLIB = $(HTSDIR)/libhts.a

CC=			gcc
CFLAGS=		-g -Wall $(VERSION) -O2
#LDFLAGS=		-Wl,-rpath,\$$ORIGIN/../lib
DFLAGS=		-D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -D_USE_KNETFILE -D_CURSES_LIB=1
LOBJS=		bam_aux.o bam.o bam_import.o sam.o \
			bam_pileup.o bam_lpileup.o sam_header.o
AOBJS=		bam_index.o bam_tview.o bam_plcmd.o sam_view.o \
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

bam_h = bam.h $(htslib_bgzf_h) $(htslib_sam_h)
bam2bcf_h = bam2bcf.h errmod.h bcftools/bcf.h
bam_tview_h = bam_tview.h $(bam_h) $(htslib_faidx_h) $(bam2bcf_h) sam_header.h $(HTSDIR)/htslib/khash.h
bcftools_bcf_h = bcftools/bcf.h $(htslib_bgzf_h)
sam_h = sam.h $(htslib_sam_h) $(bam_h)
sample_h = sample.h $(HTSDIR)/htslib/kstring.h

bam.o: bam.c $(bam_h) sam_header.h
bam2bcf.o: bam2bcf.c
bam2bcf_indel.o: bam2bcf_indel.c
bam2depth.o: bam2depth.c
bam_aux.o: bam_aux.c
bam_cat.o: bam_cat.c $(HTSDIR)/htslib/knetfile.h $(htslib_bgzf_h) $(bam_h)
bam_color.o: bam_color.c $(bam_h)
bam_import.o: bam_import.c $(bam_h) $(HTSDIR)/htslib/kseq.h
bam_index.o: bam_index.c $(bam_h)
bam_lpileup.o: bam_lpileup.c $(bam_h) $(HTSDIR)/htslib/ksort.h
bam_mate.o: bam_mate.c $(bam_h)
bam_md.o: bam_md.c $(htslib_faidx_h) $(sam_h) kaln.h kprobaln.h
bam_pileup.o: bam_pileup.c $(sam_h)
bam_plcmd.o: bam_plcmd.c $(sam_h) $(htslib_faidx_h) sam_header.h $(bam2bcf_h) $(sample_h)
bam_reheader.o: bam_reheader.c $(HTSDIR)/htslib/knetfile.h $(htslib_bgzf_h) $(bam_h)
bam_rmdup.o: bam_rmdup.c $(sam_h) $(HTSDIR)/htslib/khash.h
bam_rmdupse.o: bam_rmdupse.c $(sam_h) $(HTSDIR)/htslib/khash.h $(HTSDIR)/htslib/klist.h
bam_sort.o: bam_sort.c $(bam_h) $(HTSDIR)/htslib/ksort.h
bam_stat.o: bam_stat.c $(bam_h)
bam_tview.o: bam_tview.c $(bam_tview_h)
bam_tview_curses.o: bam_tview_curses.c $(bam_tview_h)
bam_tview_html.o: bam_tview_html.c $(bam_tview_h)
bamshuf.o: bamshuf.c $(htslib_sam_h) $(HTSDIR)/htslib/ksort.h
bamtk.o: bamtk.c $(bam_h) $(HTSDIR)/htslib/knetfile.h
bedcov.o: bedcov.c $(htslib_bgzf_h) $(bam_h) $(HTSDIR)/htslib/kseq.h
bedidx.o: bedidx.c $(HTSDIR)/htslib/ksort.h $(HTSDIR)/htslib/kseq.h $(HTSDIR)/htslib/khash.h
bgzip.o: bgzip.c $(htslib_bgzf_h)
cut_target.o: cut_target.c $(bam_h) errmod.h $(htslib_faidx_h)
errmod.o: errmod.c errmod.h $(HTSDIR)/htslib/ksort.h
kaln.o: kaln.c kaln.h
kprobaln.o: kprobaln.c kprobaln.h
padding.o: padding.c sam_header.h $(sam_h) $(bam_h) $(htslib_faidx_h)
phase.o: phase.c $(bam_h) errmod.h $(HTSDIR)/htslib/kseq.h $(HTSDIR)/htslib/khash.h $(HTSDIR)/htslib/ksort.h
razip.o: razip.c $(htslib_razf_h)
sam.o: sam.c $(htslib_faidx_h) $(sam_h)
sam_header.o: sam_header.c sam_header.h $(HTSDIR)/htslib/khash.h
sam_view.o: sam_view.c sam_header.h $(sam_h) $(htslib_faidx_h) $(HTSDIR)/htslib/khash.h
sample.o: sample.c $(sample_h) $(HTSDIR)/htslib/khash.h


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

clean:cleanlocal-recur clean-htslib
