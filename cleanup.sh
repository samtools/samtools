if [ -f "Makefile" ]; then
	make distclean
fi
rm -fr *~ .in .gdb_history Makefile.in aclocal.m4 configure autom4*.cache config.guess config.h.in config.sub depcomp install-sh missing mkinstalldirs misc/Makefile.in
