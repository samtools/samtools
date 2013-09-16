
PROG=samtools bgzip
all:$(PROG)

# Adjust $(HTSDIR) to point to your top-level htslib directory
HTSDIR = ../htslib
include $(HTSDIR)/htslib.mk
HTSLIB = $(HTSDIR)/libhts.a

CC=			gcc
CFLAGS=		-g -Wall $(VERSION) -O2
#LDFLAGS=		-Wl,-rpath,\$$ORIGIN/../lib
DFLAGS=		-D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -D_CURSES_LIB=1
LOBJS=		bam_aux.o bam.o bam_import.o sam.o \
			sam_header.o
AOBJS=		bam_index.o bam_plcmd.o sam_view.o \
			bam_cat.o bam_md.o bam_reheader.o bam_sort.o bedidx.o kprobaln.o \
			bam_rmdup.o bam_rmdupse.o bam_mate.o bam_stat.o bam_color.o \
			bamtk.o kaln.o bam2bcf.o bam2bcf_indel.o errmod.o sample.o \
			cut_target.o phase.o bam2depth.o padding.o bedcov.o bamshuf.o \
            faidx.o stats.o 
            # tview todo: bam_tview.o bam_tview_curses.o bam_tview_html.o bam_lpileup.o
INCLUDES=	-I. -I$(HTSDIR)
SUBDIRS=	. misc
LIBPATH=
LIBCURSES=	-lcurses # -lXCurses


.SUFFIXES:.c .o
.PHONY: all lib test force

force:

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

test:
		test/test.pl

.PHONY:all lib clean cleanlocal
.PHONY:all-recur lib-recur clean-recur cleanlocal-recur install-recur

# See htslib/Makefile
PACKAGE_VERSION  = 0.0.1
LIBHTS_SOVERSION = 0
NUMERIC_VERSION  = $(PACKAGE_VERSION)
ifneq "$(wildcard .git)" ""
original_version := $(PACKAGE_VERSION)
PACKAGE_VERSION := $(shell git describe --always --dirty)
ifneq "$(subst ..,.,$(subst 0,,$(subst 1,,$(subst 2,,$(subst 3,,$(subst 4,,$(subst 5,,$(subst 6,,$(subst 7,,$(subst 8,,$(subst 9,,$(PACKAGE_VERSION))))))))))))" "."
empty :=
NUMERIC_VERSION := $(subst $(empty) ,.,$(wordlist 1,2,$(subst ., ,$(original_version))) 255)
endif
version.h: $(if $(wildcard version.h),$(if $(findstring "$(PACKAGE_VERSION)",$(shell cat version.h)),,force))
endif
version.h:
	echo '#define SAMTOOLS_VERSION "$(PACKAGE_VERSION)"' > $@

lib:libbam.a

libbam.a:$(LOBJS)
		$(AR) -csru $@ $(LOBJS)

samtools:lib-recur $(AOBJS) $(HTSLIB)
		$(CC) $(CFLAGS) -o $@ $(AOBJS) $(LDFLAGS) libbam.a $(LIBPATH) $(HTSLIB) $(LIBCURSES) -lm -lz -lpthread

bgzip:bgzip.o $(HTSLIB)
		$(CC) $(CFLAGS) -o $@ bgzip.o $(HTSLIB) -lz -lpthread

bam_h = bam.h $(htslib_bgzf_h) $(htslib_sam_h)
bam2bcf_h = bam2bcf.h errmod.h
bam_tview_h = bam_tview.h $(bam_h) $(htslib_faidx_h) $(bam2bcf_h) sam_header.h $(HTSDIR)/htslib/khash.h
sam_h = sam.h $(htslib_sam_h) $(bam_h)
sample_h = sample.h $(HTSDIR)/htslib/kstring.h

bam.o: bam.c $(bam_h) sam_header.h
bam2bcf.o: bam2bcf.c bam2bcf.h $(HTSDIR)/htslib/kfunc.h
bam2bcf_indel.o: bam2bcf_indel.c bam2bcf.h
bam2depth.o: bam2depth.c
bam_aux.o: bam_aux.c
bam_cat.o: bam_cat.c $(htslib_bgzf_h) $(bam_h)
bam_color.o: bam_color.c $(bam_h)
bam_import.o: bam_import.c $(HTSDIR)/htslib/kstring.h $(bam_h) $(HTSDIR)/htslib/kseq.h
bam_index.o: bam_index.c $(bam_h)
bam_lpileup.o: bam_lpileup.c $(bam_h) $(HTSDIR)/htslib/ksort.h
bam_mate.o: bam_mate.c $(bam_h)
bam_md.o: bam_md.c $(htslib_faidx_h) $(sam_h) kaln.h kprobaln.h
bam_pileup.o: bam_pileup.c $(sam_h)
bam_plcmd.o: bam_plcmd.c $(sam_h) $(htslib_faidx_h) sam_header.h $(bam2bcf_h) $(sample_h)
bam_reheader.o: bam_reheader.c $(htslib_bgzf_h) $(bam_h)
bam_rmdup.o: bam_rmdup.c $(sam_h) $(HTSDIR)/htslib/khash.h
bam_rmdupse.o: bam_rmdupse.c $(sam_h) $(HTSDIR)/htslib/khash.h $(HTSDIR)/htslib/klist.h
bam_sort.o: bam_sort.c $(bam_h) $(HTSDIR)/htslib/ksort.h
bam_stat.o: bam_stat.c $(bam_h)
bam_tview.o: bam_tview.c $(bam_tview_h)
bam_tview_curses.o: bam_tview_curses.c $(bam_tview_h)
bam_tview_html.o: bam_tview_html.c $(bam_tview_h)
bamshuf.o: bamshuf.c $(htslib_sam_h) $(HTSDIR)/htslib/ksort.h
bamtk.o: bamtk.c $(bam_h) version.h samtools.h
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
stats.o: sam.c $(bam_h) $(HTSDIR)/htslib/khash.h $(htslib_faidx_h)

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
