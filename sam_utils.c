/*  sam_utils.c -- various utilities internal to samtools.

    Copyright (C) 2014-2016 Genome Research Ltd.

    Author: John Marshall <jm18@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#include <config.h>

#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <errno.h>

#include "samtools.h"

static void vprint_error_core(const char *subcommand, const char *format, va_list args, const char *extra)
{
    fflush(stdout);
    if (subcommand && *subcommand) fprintf(stderr, "samtools %s: ", subcommand);
    else fprintf(stderr, "samtools: ");
    vfprintf(stderr, format, args);
    if (extra) fprintf(stderr, ": %s\n", extra);
    else fprintf(stderr, "\n");
    fflush(stderr);
}

void print_error(const char *subcommand, const char *format, ...)
{
    va_list args;
    va_start(args, format);
    vprint_error_core(subcommand, format, args, NULL);
    va_end(args);
}

void print_error_errno(const char *subcommand, const char *format, ...)
{
    int err = errno;
    va_list args;
    va_start(args, format);
    vprint_error_core(subcommand, format, args, err? strerror(err) : NULL);
    va_end(args);
}
