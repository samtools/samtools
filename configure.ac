AC_INIT(SAMTOOLS, 0.1.3)
AM_CONFIG_HEADER(config.h)
AM_INIT_AUTOMAKE([no-dependencies])
AC_CANONICAL_HOST

AC_PROG_INSTALL

AC_PROG_CC

# set CFLAGS and CXXFLAGS
user_CFLAGS="${CFLAGS} -D_IOLIB=2 -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE"
generic_CFLAGS="-Wall -L."

AC_ARG_ENABLE(norazf, [  --disable-razf	  disable razf support],
			  [user_DFLAGS="${user_DFLAGS} -D_NO_RAZF"], [])
AC_ARG_ENABLE(nocurses, [  --disable-curses	  disable curses support],
			  [user_DFLAGS="${user_DFLAGS} -D_NO_CURSES"], [])
CFLAGS="${generic_CFLAGS} ${user_CFLAGS}"
CXXFLAGS=$CFLAGS

AC_STDC_HEADERS
AC_CHECK_FUNCS([bzero memset pow sqrt strcasecmp strdup strstr strtol])
AC_CHECK_HEADERS([zlib.h fcntl.h malloc.h stddef.h])
AC_CHECK_TYPES([ptrdiff_t])
AC_C_CONST
AC_C_INLINE
AC_FUNC_FSEEKO
AC_FUNC_MALLOC
AC_FUNC_MEMCMP
AC_FUNC_REALLOC
AC_HEADER_STDBOOL
AC_TYPE_OFF_T
AC_TYPE_SIZE_T

AC_OUTPUT([Makefile misc/Makefile])