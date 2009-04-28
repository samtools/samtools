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
dist_bin_SCRIPTS = bowtie2sam.pl export2sam.pl novo2sam.pl samtools.pl soap2sam.pl wgsim_eval.pl zoom2sam.pl
bin_PROGRAMS = faidx maq2sam-short maq2sam-long md5fa md5sum-lite wgsim

faidx_SOURCES = ../faidx.c ../faidx.c
#faidx_LDFLAGS = -D_NO_RAZF -DFAIDX_MAIN 

md5fa_SOURCES = md5.c md5.h md5fa.c ../kseq.h
#md5fa_LDADD = -lz

md5sum_lite_SOURCES = md5.c md5.h
#md5sum_lite_LDFLAGS = -DMD5SUM_MAIN

maq2sam_short_SOURCES = maq2sam.c
#maq2sam_short_LDADD = -lz

maq2sam_long_SOURCES = maq2sam.c
#maq2sam_long_LDADD = -lz

wgsim_SOURCES = wgsim.c 
wgsim_LDADD = -lm 

faidx:../faidx.c ../faidx.h
	$(CC) $(CFLAGS) -D_NO_RAZF -DFAIDX_MAIN -o $@ ../faidx.c
md5fa:md5.o md5fa.o md5.h ../kseq.h
	$(CC) $(CFLAGS) -o $@ md5.o md5fa.o -lz
md5fa.o:md5.h md5fa.c
	$(CC) $(CFLAGS) -c -I.. -o $@ md5fa.c
md5sum-lite:md5sum-lite.o
	$(CC) $(CFLAGS) -o $@ md5sum-lite.o
md5sum-lite.o:md5.c md5.h
	$(CC) -c $(CFLAGS) -DMD5SUM_MAIN -o $@ md5.c
maq2sam-short:maq2sam-short.o
	$(CC) $(CFLAGS) -o $@ maq2sam-short.o -lz 
maq2sam-short.o:maq2sam.c
	$(CC) -c $(CFLAGS) -o $@ maq2sam.c 
maq2sam-long:maq2sam-long.o
	$(CC) $(CFLAGS) -o $@ maq2sam-long.o -lz 
maq2sam-long.o:maq2sam.c
	$(CC) -c $(CFLAGS) -DMAQ_LONGREADS -D_MAQ_LONGREADS -o $@ maq2sam.c 
