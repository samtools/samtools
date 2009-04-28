EXTRA_DIST = autogen.sh 
bin_PROGRAMS = bgzip samtools razip
SUBDIRS = . misc

bgzip_SOURCES = bgzip.c bgzf.c
bgzip_LDADD = -lm -lz

samtools_SOURCES = bamtk.c bam.c bam_aux.c bam_import.c bam_index.c \
				   bam_lpileup.c bam_maqcns.c bam_mate.c bam_md.c bam_pileup.c \
				   bam_plcmd.c bam_rmdup.c bam_sort.c bam_stat.c \
				   bam_tview.c bgzf.c faidx.c glf.c kstring.c razf.c \
				   bam.h bam_endian.h bam_maqcns.h bgzf.h faidx.h glf.h \
				   khash.h kseq.h kstring.h ksort.h razf.h zutil.h
samtools_LDADD = -lm -lz -lcurses

razip_SOURCES = razip.c razf.c
razip_LDADD = -lm -lz

