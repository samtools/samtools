dnl @synopsis HTS_PROG_CC_WARNINGS([ANSI])
dnl
dnl Derived from
dnl     http://ac-archive.sourceforge.net/ac-archive/vl_prog_cc_warnings.html
dnl
dnl Enables a reasonable set of warnings for the C compiler.
dnl Optionally, if the first argument is nonempty, turns on flags which
dnl enforce and/or enable proper ANSI C if such are known with the
dnl compiler used.
dnl
dnl Currently this macro knows about GCC, Solaris C compiler, Digital
dnl Unix C compiler, C for AIX Compiler, HP-UX C compiler, IRIX C
dnl compiler, NEC SX-5 (Super-UX 10) C compiler, and Cray J90 (Unicos
dnl 10.0.0.8) C compiler.
dnl
dnl @category C
dnl @author Ville Laurikari <vl@iki.fi>
dnl Updated by Rob Davies <rmd@sanger.ac.uk> for HTSlib
dnl @license AllPermissive
dnl Copying and distribution of this file, with or without modification,
dnl are permitted in any medium without royalty provided the copyright notice
dnl and this notice are preserved. Users of this software should generally
dnl follow the principles of the MIT License including its disclaimer.
dnl Original Copyright (c) Ville Laurikari 2002
dnl Modifications Copyright (c) Genome Research Limited 2015,2017

AC_DEFUN([HTS_PROG_CC_WARNINGS],[
  AC_ARG_ENABLE([warnings],
    [AS_HELP_STRING([--disable-warnings], [turn off compiler warnings])],
    [],
    [enable_warnings=yes])

  AS_IF([test "x$enable_warnings" != xno],[
    AC_REQUIRE([AC_PROG_GREP])

    ansi="$1"
    AS_IF([test "x$ansi" = "x"],
          [msg="for C compiler warning flags"],
          [msg="for C compiler warning and ANSI conformance flags"])

    AC_MSG_CHECKING($msg)
    AC_CACHE_VAL(hts_cv_prog_cc_warnings,[dnl
      hts_cv_prog_cc_warnings=""
      AS_IF([test "x$CC" != "x"],[
        cat > conftest.c <<EOF
int main(int argc, char **argv) { return 0; }
EOF

dnl Most compilers print some kind of a version string with some command
dnl line options (often "-V").  The version string should be checked
dnl before doing a test compilation run with compiler-specific flags.
dnl This is because some compilers (like the Cray compiler) only
dnl produce a warning message for unknown flags instead of returning
dnl an error, resulting in a false positive.  Also, compilers may do
dnl erratic things when invoked with flags meant for a different
dnl compiler.

dnl We attempt to strip out any flags that are already on CFLAGS.
dnl If an option needs more than one word (e.g. see Cray below) then
dnl they should be separated by hash signs (#), which will be converted
dnl to spaces before comparing and possibly adding to CFLAGS.
dnl This separator will need to be changed if a new compiler ever needs
dnl an option that includes a hash sign...

        # Tests for flags to enable C compiler warnings
        # GCC compatible
        AS_IF([test "x$GCC" = "xyes" &&
               "$CC" -c -Wall conftest.c > /dev/null 2>&1 &&
               test -f conftest.o],[dnl
          AS_IF([test "x$ansi" = "x"],
                [hts_cv_prog_cc_warnings="-Wall"],
                [hts_cv_prog_cc_warnings="-Wall -ansi -pedantic"])
        ],
        # Sun Studio or Solaris C compiler
        ["$CC" -V 2>&1 | $GREP -i -E "WorkShop|Sun C" > /dev/null 2>&1 &&
         "$CC" -c -v -Xc conftest.c > /dev/null 2>&1 &&
         test -f conftest.o],[dnl
          AS_IF([test "x$ansi" = "x"],
                [hts_cv_prog_cc_warnings="-v"],
                [hts_cv_prog_cc_warnings="-v -Xc"])
        ],
        # Digital Unix C compiler
        ["$CC" -V 2>&1 | $GREP -i "Digital UNIX Compiler" > /dev/null 2>&1 &&
         "$CC" -c -verbose -w0 -warnprotos -std1 conftest.c > /dev/null 2>&1 &&
         test -f conftest.o],[dnl
          AS_IF([test "x$ansi" = "x"],
                [hts_cv_prog_cc_warnings="-verbose -w0 -warnprotos"],
                [hts_cv_prog_cc_warnings="-verbose -w0 -warnprotos -std1"])
        ],
        # C for AIX Compiler
        ["$CC" 2>&1 | $GREP -i "C for AIX Compiler" > /dev/null 2>&1 &&
         "$CC" -c -qlanglvl=ansi -qinfo=all conftest.c > /dev/null 2>&1 &&
         test -f conftest.o],[dnl
          AS_IF([test "x$ansi" = "x"],
                [hts_cv_prog_cc_warnings="-qsrcmsg -qinfo=all:noppt:noppc:noobs:nocnd"],
                [hts_cv_prog_cc_warnings="-qsrcmsg -qinfo=all:noppt:noppc:noobs:nocnd -qlanglvl=ansi"])
        ],
        # IRIX C compiler
        ["$CC" -version 2>&1 | $GREP -i "MIPSpro Compilers" > /dev/null 2>&1 &&
         "$CC" -c -fullwarn -ansi -ansiE conftest.c > /dev/null 2>&1 &&
         test -f conftest.o],[dnl
          AS_IF([test "x$ansi" = "x"],
                [hts_cv_prog_cc_warnings="-fullwarn"],
                [hts_cv_prog_cc_warnings="-fullwarn -ansi -ansiE"])
        ],
        # HP-UX C compiler
        [what "$CC" 2>&1 | $GREP -i "HP C Compiler" > /dev/null 2>&1 &&
         "$CC" -c -Aa +w1 conftest.c > /dev/null 2>&1 &&
         test -f conftest.o],[dnl
          AS_IF([test "x$ansi" = "x"],
                [hts_cv_prog_cc_warnings="+w1"],
                [hts_cv_prog_cc_warnings="+w1 -Aa"])
        ],
        # The NEC SX series (Super-UX 10) C compiler
        ["$CC" -V 2>&1 | $GREP "/SX" > /dev/null 2>&1 &&
         "$CC" -c -pvctl[,]fullmsg -Xc conftest.c > /dev/null 2>&1 &&
         test -f conftest.o],[dnl
          AS_IF([test "x$ansi" = "x"],
                [hts_cv_prog_cc_warnings="-pvctl[,]fullmsg"],
                [hts_cv_prog_cc_warnings="-pvctl[,]fullmsg -Xc"])
        ],
        # The Cray C compiler (Unicos)
        ["$CC" -V 2>&1 | $GREP -i "Cray" > /dev/null 2>&1 &&
         "$CC" -c -h msglevel_2 conftest.c > /dev/null 2>&1 &&
         test -f conftest.o],[dnl
          AS_IF([test "x$ansi" = "x"],
                [hts_cv_prog_cc_warnings="-h#msglevel_2"],
                [hts_cv_prog_cc_warnings="-h#msglevel_2,conform"])
        ],
        # The Tiny C Compiler
        ["$CC" -v 2>&1 | $GREP "tcc version" > /dev/null &&
         "$CC" -Wall -c conftest.c > /dev/null 2>&1 &&
         test -f conftest.o],[dnl
         hts_cv_prog_cc_warnings="-Wall"
        ])
        rm -f conftest.*
      ])
    ])

    AS_IF([test "x$hts_cv_prog_cc_warnings" != "x"],[
dnl Print result, with underscores as spaces
ac_arg_result=`echo "$hts_cv_prog_cc_warnings" | tr '#' ' '`
AC_MSG_RESULT($ac_arg_result)

dnl Add options to CFLAGS only if they are not already present
ac_arg_needed=""
for ac_arg in $hts_cv_prog_cc_warnings
do
  ac_arg_sp=`echo "$ac_arg" | tr '#' ' '`
  AS_CASE([" $CFLAGS "],
[*" $ac_arg_sp "*], [],
[ac_arg_needed="$ac_arg_all $ac_arg_sp"])
done
CFLAGS="$ac_arg_needed $CFLAGS"],[dnl
      AC_MSG_RESULT(unknown)
    ])
  ])
])dnl HTS_PROG_CC_WARNINGS

# SYNOPSIS
#
# HTS_PROG_CC_WERROR(FLAGS_VAR)
#
# Set FLAGS_VAR to the flags needed to make the C compiler treat warnings
# as errors.

AC_DEFUN([HTS_PROG_CC_WERROR], [
  AC_ARG_ENABLE([werror],
    [AS_HELP_STRING([--enable-werror], [change warnings into errors, where supported])],
    [],
    [enable_werror=no])

  AS_IF([test "x$enable_werror" != xno],[
    AC_MSG_CHECKING([for C compiler flags to error on warnings])
    AC_CACHE_VAL(hts_cv_prog_cc_werror,[dnl
      hts_cv_prog_cc_werror=""
      AS_IF([test "x$CC" != "x"],[
        cat > conftest.c <<EOF
int main(int argc, char **argv) { return 0; }
EOF

        AS_IF(dnl
         # Tests for flags to make the C compiler treat warnings as errors
         # GCC compatible
         [test "x$GCC" = "xyes" &&
          "$CC" -c -Werror conftest.c > /dev/null 2>&1 &&
          test -f conftest.o],[hts_cv_prog_cc_werror="-Werror"],
         # Sun Studio or Solaris C compiler
         ["$CC" -V 2>&1 | $GREP -i -E "WorkShop|Sun C" > /dev/null 2>&1 &&
          "$CC" -c -errwarn=%all conftest.c > /dev/null 2>&1 &&
          test -f conftest.o],[hts_cv_prog_cc_werror="-errwarn=%all"],
         # The Tiny C Compiler
         ["$CC" -v 2>&1 | $GREP "tcc version" > /dev/null &&
          "$CC" -Wall -c conftest.c > /dev/null 2>&1 &&
          test -f conftest.o],[hts_cv_prog_cc_werror="-Werror"]
         dnl TODO: Add more compilers
        )
        rm -f conftest.*
      ])
    ])
    AS_IF([test "x$hts_cv_prog_cc_werror" != x],[dnl
      AC_MSG_RESULT($hts_cv_prog_cc_werror)
      AS_IF([test "x$1" != x],[eval AS_TR_SH([$1])="$hts_cv_prog_cc_werror"])
    ],[dnl
      AC_MSG_RESULT(unknown)
    ])
  ])
])dnl HTS_PROG_CC_WERROR
