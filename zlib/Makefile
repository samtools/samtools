CC=			gcc
CFLAGS=		-g -Wall -O2
DFLAGS=		-D_FILE_OFFSET_BITS=64
ZOBJS=		adler32.o compress.o crc32.o deflate.o gzio.o inffast.o inflate.o \
			infback.o inftrees.o trees.o uncompr.o zutil.o
LIBS=		-lm

.SUFFIXES:.c .o

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $< -o $@

all:lib

lib:libz.a

libz.a:$(ZOBJS)
		$(AR) -cru $@ $(ZOBJS)

adler32.o: zlib.h zconf.h
compress.o: zlib.h zconf.h
crc32.o: crc32.h zlib.h zconf.h
deflate.o: deflate.h zutil.h zlib.h zconf.h
gzio.o: zutil.h zlib.h zconf.h
inffast.o: zutil.h zlib.h zconf.h inftrees.h inflate.h inffast.h
inflate.o: zutil.h zlib.h zconf.h inftrees.h inflate.h inffast.h
infback.o: zutil.h zlib.h zconf.h inftrees.h inflate.h inffast.h
inftrees.o: zutil.h zlib.h zconf.h inftrees.h
trees.o: deflate.h zutil.h zlib.h zconf.h trees.h
uncompr.o: zlib.h zconf.h
zutil.o: zutil.h zlib.h zconf.h

cleanlocal:
		rm -f *.o a.out *~ *.a libzr.a gmon.out

clean:cleanlocal

