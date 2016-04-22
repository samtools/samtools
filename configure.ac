# Configure script for samtools, utilities for Sequence Alignment/Map files.
#
#    Copyright (C) 2015 Genome Research Ltd.
#
#    Author: John Marshall <jm18@sanger.ac.uk>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.

dnl Process this file with autoconf to produce a configure script
AC_INIT([Samtools], m4_esyscmd_s([make print-version]),
        [samtools-help@lists.sourceforge.net], [], [http://www.htslib.org/])
AC_PREREQ([2.63])  dnl This version introduced 4-argument AC_CHECK_HEADER
AC_CONFIG_SRCDIR([bamtk.c])
AC_CONFIG_HEADERS([config.h])

dnl Copyright notice to be copied into the generated configure script
AC_COPYRIGHT([Portions copyright (C) 2015 Genome Research Ltd.

This configure script is free software: you are free to change and
redistribute it.  There is NO WARRANTY, to the extent permitted by law.])

AC_PROG_CC

AC_SYS_LARGEFILE

AX_WITH_HTSLIB
if test "$ax_cv_htslib" != yes; then
  AC_MSG_ERROR([HTSlib development files not found

Samtools uses HTSlib to parse bioinformatics file formats etc.  Building it
requires an unpacked HTSlib source tree (which will be built in conjunction
with samtools) or a previously-installed HTSlib.  In either case you may
need to configure --with-htslib=DIR to locate the appropriate HTSlib.

FAILED.  You must supply an HTSlib in order to build samtools successfully.])
fi

if test "$ax_cv_htslib_which" = source; then
  Hsource=
  Hinstall='#'
else
  Hsource='#'
  Hinstall=
fi
AC_SUBST([Hsource])
AC_SUBST([Hinstall])

AC_ARG_WITH([curses],
  [AS_HELP_STRING([--without-curses],
                  [omit curses support, so no curses library needed])])

if test "$with_curses" != no; then
  AX_WITH_CURSES
  if test "$ax_cv_curses" != yes; then
    AC_MSG_ERROR([curses development files not found

The 'samtools tview' command uses the curses text user interface library.
Building samtools with tview requires curses/ncurses/etc development files
to be installed on the build machine; you may need to ensure a package such
as libncurses5-dev (on Debian or Ubuntu Linux) or ncurses-devel (on RPM-based
Linux distributions) is installed.

FAILED.  Either configure --without-curses or resolve this error to build
samtools successfully.])
  fi
fi

save_LIBS=$LIBS
zlib_devel=ok
dnl Set a trivial non-empty INCLUDES to avoid excess default includes tests
AC_CHECK_HEADER([zlib.h], [], [zlib_devel=missing], [;])
AC_CHECK_LIB(z, inflate,  [], [zlib_devel=missing])
LIBS=$save_LIBS

if test $zlib_devel != ok; then
  AC_MSG_ERROR([zlib development files not found

Samtools uses compression routines from the zlib library <http://zlib.net>.
Building samtools requires zlib development files to be installed on the build
machine; you may need to ensure a package such as zlib1g-dev (on Debian or
Ubuntu Linux) or zlib-devel (on RPM-based Linux distributions) is installed.

FAILED.  This error must be resolved in order to build samtools successfully.])
fi

AC_CONFIG_FILES([config.mk])
AC_OUTPUT
