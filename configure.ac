# Configure script for samtools, utilities for Sequence Alignment/Map files.
#
#    Copyright (C) 2015, 2017, 2019 Genome Research Ltd.
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
AC_INIT([Samtools], m4_esyscmd_s([./version.sh 2>/dev/null]),
        [samtools-help@lists.sourceforge.net], [], [http://www.htslib.org/])
AC_PREREQ([2.63])  dnl This version introduced 4-argument AC_CHECK_HEADER
AC_CONFIG_SRCDIR([bamtk.c])
AC_CONFIG_HEADERS([config.h])

m4_include([m4/hts_prog_cc_warnings.m4])
m4_include([m4/hts_with_curses.m4])
m4_include([m4/ax_with_htslib.m4])

dnl Copyright notice to be copied into the generated configure script
AC_COPYRIGHT([Portions copyright (C) 2015,2017, 2019 Genome Research Ltd.

This configure script is free software: you are free to change and
redistribute it.  There is NO WARRANTY, to the extent permitted by law.])

dnl Notes to be copied (by autoheader) into the generated config.h.in
AH_TOP([/* If you use configure, this file provides @%:@defines reflecting your
   configuration choices.  If you have not run configure, suitable
   conservative defaults will be used.

   Autoheader adds a number of items to this template file that are not
   used by samtools: STDC_HEADERS and most HAVE_*_H header file defines
   are immaterial, as we assume standard ISO C headers and facilities;
   the PACKAGE_* defines are unused and are overridden by the more
   accurate PACKAGE_VERSION as computed by the Makefile.  */])

AC_PROG_AWK
AC_PROG_CC

dnl Turn on compiler warnings, if possible
HTS_PROG_CC_WARNINGS
dnl Flags to treat warnings as errors.  These need to be applied to CFLAGS
dnl later as they can interfere with some of the tests (notably AC_SEARCH_LIBS)
HTS_PROG_CC_WERROR(hts_late_cflags)

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
  HTS_WITH_CURSES
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
else
  CURSES_LIB=
  AC_SUBST([CURSES_LIB])
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

dnl CUDA support
AC_ARG_ENABLE([cuda],
  [AS_HELP_STRING([--enable-cuda], [Enable CUDA GPU acceleration (requires NVIDIA GPU)])
AS_HELP_STRING([--disable-cuda], [Disable CUDA GPU acceleration])],
  [enable_cuda=$enableval], [enable_cuda=auto])

dnl Intelligent CUDA auto-detection with hardware checks
if test "x$enable_cuda" = "xauto"; then
  AC_MSG_CHECKING([for NVIDIA GPU hardware])
  
  # Step 1: Check for NVIDIA GPU hardware first
  gpu_detected=no
  gpu_info=""
  
  # Method 1: Try nvidia-smi for GPU detection
  AC_PATH_PROG([NVIDIA_SMI], [nvidia-smi], [no])
  if test "x$NVIDIA_SMI" != "xno"; then
    # Check if nvidia-smi can detect GPUs
    gpu_count=`$NVIDIA_SMI -L 2>/dev/null | wc -l`
    if test "$gpu_count" -gt 0; then
      gpu_detected=yes
      gpu_info=`$NVIDIA_SMI -L 2>/dev/null | head -1`
      AC_MSG_RESULT([yes - found $gpu_count NVIDIA GPU(s)])
      AC_MSG_NOTICE([GPU detected: $gpu_info])
    fi
  fi
  
  # Method 2: Check for NVIDIA devices in /proc/driver/nvidia if nvidia-smi failed
  if test "x$gpu_detected" = "xno"; then
    if test -d "/proc/driver/nvidia" && test -r "/proc/driver/nvidia/version"; then
      gpu_detected=yes
      nvidia_version=`cat /proc/driver/nvidia/version 2>/dev/null | head -1`
      AC_MSG_RESULT([yes - NVIDIA driver detected])
      AC_MSG_NOTICE([NVIDIA driver: $nvidia_version])
    fi
  fi
  
  # Method 3: Check lspci for NVIDIA VGA/3D controllers
  if test "x$gpu_detected" = "xno"; then
    AC_PATH_PROG([LSPCI], [lspci], [no])
    if test "x$LSPCI" != "xno"; then
      nvidia_devices=`$LSPCI 2>/dev/null | grep -i "nvidia.*\(vga\|3d\|display\)" | wc -l`
      if test "$nvidia_devices" -gt 0; then
        gpu_detected=yes
        AC_MSG_RESULT([yes - found $nvidia_devices NVIDIA device(s) via lspci])
        nvidia_device_info=`$LSPCI 2>/dev/null | grep -i "nvidia.*\(vga\|3d\|display\)" | head -1`
        AC_MSG_NOTICE([NVIDIA device: $nvidia_device_info])
      fi
    fi
  fi
  
  # If no GPU detected, inform user and skip CUDA
  if test "x$gpu_detected" = "xno"; then
    AC_MSG_RESULT([no])
    AC_MSG_NOTICE([])
    AC_MSG_NOTICE([=================================================================])
    AC_MSG_NOTICE([No NVIDIA GPU hardware detected on this system.])
    AC_MSG_NOTICE([CUDA GPU acceleration will be disabled.])
    AC_MSG_NOTICE([])
    AC_MSG_NOTICE([If you have an NVIDIA GPU:])
    AC_MSG_NOTICE([  1. Install NVIDIA drivers: https://www.nvidia.com/drivers])
    AC_MSG_NOTICE([  2. Install CUDA toolkit: https://developer.nvidia.com/cuda-downloads])
    AC_MSG_NOTICE([  3. Reconfigure with: ./configure --enable-cuda])
    AC_MSG_NOTICE([])
    AC_MSG_NOTICE([If you want to use CPU-only mode, this is normal.])
    AC_MSG_NOTICE([=================================================================])
    AC_MSG_NOTICE([])
    enable_cuda=no
  else
    # Step 2: GPU detected, now check for CUDA software stack
    AC_MSG_CHECKING([for CUDA software stack])
    
    AC_PATH_PROG([NVCC], [nvcc], [no])
    
    if test "x$NVCC" = "xno"; then
      AC_MSG_RESULT([no - nvcc not found])
      AC_MSG_NOTICE([])
      AC_MSG_NOTICE([=================================================================])
      AC_MSG_NOTICE([NVIDIA GPU detected but CUDA toolkit not installed.])
      AC_MSG_NOTICE([])
      AC_MSG_NOTICE([To enable GPU acceleration:])
      AC_MSG_NOTICE([  1. Install CUDA toolkit: https://developer.nvidia.com/cuda-downloads])
      AC_MSG_NOTICE([  2. Make sure 'nvcc' is in your PATH])
      AC_MSG_NOTICE([  3. Reconfigure with: ./configure --enable-cuda])
      AC_MSG_NOTICE([])
      AC_MSG_NOTICE([Building without GPU acceleration for now.])
      AC_MSG_NOTICE([=================================================================])
      AC_MSG_NOTICE([])
      enable_cuda=no
    else
      # Step 3: Test CUDA functionality
      AC_MSG_RESULT([found nvcc])
      AC_MSG_CHECKING([if CUDA is functional])
      
      # Create test program
      cat > conftest.cu << 'EOF'
#include <cuda_runtime.h>
#include <stdio.h>
int main() { 
  int count = 0;
  cudaError_t err = cudaGetDeviceCount(&count);
  if (err != cudaSuccess) {
    printf("CUDA Error: %s\n", cudaGetErrorString(err));
    return 1;
  }
  if (count == 0) {
    printf("No CUDA devices found\n");
    return 1;
  }
  printf("Found %d CUDA device(s)\n", count);
  return 0;
}
EOF
      
      # Get CUDA path from nvcc
      CUDA_PATH_AUTO=`dirname $(dirname $(which nvcc 2>/dev/null) 2>/dev/null) 2>/dev/null`
      if test -z "$CUDA_PATH_AUTO" || test ! -d "$CUDA_PATH_AUTO/include"; then
        # Try common locations
        for cuda_path in /usr/local/cuda /opt/cuda /usr/cuda; do
          if test -d "$cuda_path/include"; then
            CUDA_PATH_AUTO="$cuda_path"
            break
          fi
        done
      fi
      
      if test -d "$CUDA_PATH_AUTO/include"; then
        # Try lib64 first, then lib
        if test -d "$CUDA_PATH_AUTO/lib64"; then
          CUDA_LIB_PATH="$CUDA_PATH_AUTO/lib64"
        elif test -d "$CUDA_PATH_AUTO/lib"; then
          CUDA_LIB_PATH="$CUDA_PATH_AUTO/lib"
        else
          CUDA_LIB_PATH=""
        fi
        
        if test -n "$CUDA_LIB_PATH"; then
          TEST_CPPFLAGS="-I$CUDA_PATH_AUTO/include"
          TEST_LDFLAGS="-L$CUDA_LIB_PATH"
          TEST_LIBS="-lcudart"
          
          # Try to compile and run test
          if $NVCC $TEST_CPPFLAGS conftest.cu -o conftest $TEST_LDFLAGS $TEST_LIBS >/dev/null 2>&1; then
            test_output=`./conftest 2>&1`
            test_result=$?
            
            if test $test_result -eq 0; then
              enable_cuda=yes
              AC_MSG_RESULT([yes])
              AC_MSG_NOTICE([CUDA test result: $test_output])
            else
              enable_cuda=no
              AC_MSG_RESULT([no - runtime test failed])
              AC_MSG_NOTICE([CUDA test failed: $test_output])
              AC_MSG_NOTICE([])
              AC_MSG_NOTICE([Possible issues:])
              AC_MSG_NOTICE([  - NVIDIA driver not loaded or incompatible])
              AC_MSG_NOTICE([  - CUDA driver version mismatch])
              AC_MSG_NOTICE([  - GPU not accessible to current user])
              AC_MSG_NOTICE([  - Try: sudo modprobe nvidia])
              AC_MSG_NOTICE([])
            fi
          else
            enable_cuda=no
            AC_MSG_RESULT([no - compilation failed])
            AC_MSG_NOTICE([])
            AC_MSG_NOTICE([CUDA compilation test failed.])
            AC_MSG_NOTICE([Check that CUDA toolkit is properly installed.])
            AC_MSG_NOTICE([])
          fi
        else
          enable_cuda=no
          AC_MSG_RESULT([no - CUDA libraries not found])
        fi
      else
        enable_cuda=no
        AC_MSG_RESULT([no - CUDA headers not found])
      fi
      
      rm -f conftest.cu conftest
    fi
  fi
fi

if test "x$enable_cuda" = "xyes"; then
  AC_PATH_PROG([NVCC], [nvcc], [no])
  AC_PATH_PROG([NVIDIA_SMI], [nvidia-smi], [no])
  
  if test "x$NVCC" = "xno"; then
    AC_MSG_ERROR([NVCC compiler not found. Please install CUDA toolkit or use --disable-cuda.])
  fi
  
  if test "x$NVIDIA_SMI" = "xno"; then
    AC_MSG_WARN([nvidia-smi not found. GPU detection may be limited.])
  fi
  
  # Auto-detect CUDA installation path
  AC_MSG_CHECKING([for CUDA installation path])
  
  # Try multiple methods to find CUDA path
  if test -n "$CUDA_PATH" && test -d "$CUDA_PATH"; then
    CUDA_ROOT="$CUDA_PATH"
  elif test -n "$CUDA_HOME" && test -d "$CUDA_HOME"; then
    CUDA_ROOT="$CUDA_HOME"
  else
    # Try to get path from nvcc
    CUDA_ROOT=`dirname $(dirname $(which nvcc 2>/dev/null) 2>/dev/null) 2>/dev/null`
    if test ! -d "$CUDA_ROOT/include"; then
      # Common installation paths
      for path in /usr/local/cuda /opt/cuda /usr/cuda; do
        if test -d "$path/include"; then
          CUDA_ROOT="$path"
          break
        fi
      done
    fi
  fi
  
  if test -d "$CUDA_ROOT/include" && test -d "$CUDA_ROOT/lib64"; then
    AC_MSG_RESULT([$CUDA_ROOT])
    CUDA_CPPFLAGS="-I$CUDA_ROOT/include"
    CUDA_LDFLAGS="-L$CUDA_ROOT/lib64"
  elif test -d "$CUDA_ROOT/include" && test -d "$CUDA_ROOT/lib"; then
    AC_MSG_RESULT([$CUDA_ROOT])
    CUDA_CPPFLAGS="-I$CUDA_ROOT/include"
    CUDA_LDFLAGS="-L$CUDA_ROOT/lib"
  else
    AC_MSG_ERROR([CUDA installation not found. Please set CUDA_PATH or install CUDA toolkit.])
  fi
  
  # Check for CUDA headers and libraries
  save_CPPFLAGS=$CPPFLAGS
  save_LDFLAGS=$LDFLAGS
  save_LIBS=$LIBS
  
  CPPFLAGS="$CPPFLAGS $CUDA_CPPFLAGS"
  LDFLAGS="$LDFLAGS $CUDA_LDFLAGS"
  
  AC_CHECK_HEADER([cuda_runtime.h], [], [
    AC_MSG_ERROR([CUDA headers not found in $CUDA_ROOT/include.])
  ])
  
  AC_CHECK_LIB([cudart], [cudaMalloc], [], [
    AC_MSG_ERROR([CUDA runtime library not found. Please check CUDA installation.])
  ])
  
  # Auto-detect GPU architecture
  AC_MSG_CHECKING([for optimal CUDA architecture])
  
  if test "x$NVIDIA_SMI" != "xno"; then
    # Get GPU compute capability from nvidia-smi
    GPU_ARCH=`$NVIDIA_SMI --query-gpu=compute_cap --format=csv,noheader,nounits | head -1 | tr -d '.' 2>/dev/null`
    
    if test -n "$GPU_ARCH" && test "$GPU_ARCH" -ge 50; then
      CUDA_ARCH="sm_$GPU_ARCH"
      AC_MSG_RESULT([$CUDA_ARCH (detected)])
    else
      # Fallback: try to compile test programs for different architectures
      for arch in 80 75 70 61 60 52 50; do
        cat > conftest_arch.cu << 'EOF'
__global__ void test_kernel() { }
int main() { test_kernel<<<1,1>>>(); return 0; }
EOF
        if $NVCC -arch=sm_$arch conftest_arch.cu -o conftest_arch $CUDA_CPPFLAGS $CUDA_LDFLAGS -lcudart >/dev/null 2>&1; then
          CUDA_ARCH="sm_$arch"
          break
        fi
      done
      rm -f conftest_arch.cu conftest_arch
      
      if test -z "$CUDA_ARCH"; then
        CUDA_ARCH="sm_50"  # Safe fallback
      fi
      AC_MSG_RESULT([$CUDA_ARCH (fallback)])
    fi
  else
    CUDA_ARCH="sm_50"  # Default fallback
    AC_MSG_RESULT([$CUDA_ARCH (default)])
  fi
  
  # Restore original flags
  CPPFLAGS=$save_CPPFLAGS
  LDFLAGS=$save_LDFLAGS
  LIBS=$save_LIBS
  
  # Set CUDA variables
  CUDA_CFLAGS="-DENABLE_CUDA"
  CUDA_LIBS="-lcudart -lcuda"
  CUDA_OBJS="cuda/cuda_config.o cuda/cuda_sort.o cuda/cuda_pileup.o cuda/cuda_stats.o"
  
  AC_SUBST([NVCC])
  AC_SUBST([CUDA_CPPFLAGS])
  AC_SUBST([CUDA_CFLAGS])
  AC_SUBST([CUDA_LDFLAGS])
  AC_SUBST([CUDA_LIBS])
  AC_SUBST([CUDA_ARCH])
  AC_SUBST([CUDA_OBJS])
  
  AC_DEFINE([ENABLE_CUDA], [1], [Enable CUDA GPU acceleration])
  enable_cuda_result="yes"
else
  CUDA_CPPFLAGS=""
  CUDA_CFLAGS=""
  CUDA_LDFLAGS=""
  CUDA_LIBS=""
  CUDA_ARCH=""
  CUDA_OBJS=""
  NVCC=""
  
  AC_SUBST([NVCC])
  AC_SUBST([CUDA_CPPFLAGS])
  AC_SUBST([CUDA_CFLAGS])
  AC_SUBST([CUDA_LDFLAGS])
  AC_SUBST([CUDA_LIBS])
  AC_SUBST([CUDA_ARCH])
  AC_SUBST([CUDA_OBJS])
  
  enable_cuda_result="no"
fi

dnl Look for regcomp in various libraries (needed on windows/mingw).
AC_SEARCH_LIBS(regcomp, regex, [libregex=needed], [])

dnl Force POSIX mode on Windows/Mingw
test -n "$host_alias" || host_alias=unknown-`uname -s`
case $host_alias in
  *-msys* | *-MSYS* | *-mingw* | *-MINGW*)
    host_result="MSYS dll"
    PLATFORM=MSYS
    PLUGIN_EXT=.dll
    # This also sets __USE_MINGW_ANSI_STDIO which in turn makes PRId64,
    # %lld and %z printf formats work.  It also enforces the snprintf to
    # be C99 compliant so it returns the correct values (in kstring.c).
    CPPFLAGS="$CPPFLAGS -D_XOPEN_SOURCE=700"
    ;;
esac

dnl Apply value from HTS_PROG_CC_WERROR (if set)
AS_IF([test "x$hts_late_cflags" != x],[CFLAGS="$CFLAGS $hts_late_cflags"])

AC_CONFIG_FILES([config.mk])
AC_OUTPUT
