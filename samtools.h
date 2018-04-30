/*  samtools.h -- utility routines.

    Copyright (C) 2013-2015 Genome Research Ltd.

    Author: Petr Danecek <pd3@sanger.ac.uk>

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

#ifndef SAMTOOLS_H
#define SAMTOOLS_H

const char *samtools_version(void);
const char *samtools_version_short(void);

#if defined __GNUC__ && __GNUC__ >= 2
#define CHECK_PRINTF(fmt,args) __attribute__ ((format (printf, fmt, args)))
#else
#define CHECK_PRINTF(fmt,args)
#endif

void print_error(const char *subcommand, const char *format, ...) CHECK_PRINTF(2, 3);
void print_error_errno(const char *subcommand, const char *format, ...) CHECK_PRINTF(2, 3);

#endif
