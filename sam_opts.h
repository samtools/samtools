/*  sam_opts.h -- utilities to aid parsing common command line options.

    Copyright (C) 2015, 2019 Genome Research Ltd.

    Author: James Bonfield <jkb@sanger.ac.uk>

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

#ifndef SAM_OPTS_H
#define SAM_OPTS_H

#include <stdio.h>
#include <limits.h>
#include <getopt.h>
#include <htslib/hts.h>

typedef struct sam_global_args {
    htsFormat in;
    htsFormat out;
    char *reference;
    int nthreads;
    int write_index;
} sam_global_args;

#define SAM_GLOBAL_ARGS_INIT {{0},{0}}

enum {
    SAM_OPT_INPUT_FMT = CHAR_MAX+1,
    SAM_OPT_INPUT_FMT_OPTION,
    SAM_OPT_OUTPUT_FMT,
    SAM_OPT_OUTPUT_FMT_OPTION,
    SAM_OPT_REFERENCE,
    SAM_OPT_NTHREADS,
    SAM_OPT_WRITE_INDEX,
    SAM_OPT_VERBOSITY,
};

#define SAM_OPT_VAL(val, defval) ((val) == '-')? '?' : (val)? (val) : (defval)

// Use this within struct option lopts[] = {...} to add the standard global
// options.  The arguments determine whether the corresponding option is
// enabled and, if so, whether it has a short option equivalent:
// 0      No short option has been assigned. Use --long-opt only.
// '-'    Both long and short options are disabled.
// <c>    Otherwise the equivalent short option is character <c>.
#define SAM_OPT_GLOBAL_OPTIONS(o1, o2, o3, o4, o5, o6) \
    {"input-fmt",         required_argument, NULL, SAM_OPT_VAL(o1, SAM_OPT_INPUT_FMT)}, \
    {"input-fmt-option",  required_argument, NULL, SAM_OPT_VAL(o2, SAM_OPT_INPUT_FMT_OPTION)}, \
    {"output-fmt",        required_argument, NULL, SAM_OPT_VAL(o3, SAM_OPT_OUTPUT_FMT)}, \
    {"output-fmt-option", required_argument, NULL, SAM_OPT_VAL(o4, SAM_OPT_OUTPUT_FMT_OPTION)}, \
    {"reference",         required_argument, NULL, SAM_OPT_VAL(o5, SAM_OPT_REFERENCE)}, \
    {"threads",           required_argument, NULL, SAM_OPT_VAL(o6, SAM_OPT_NTHREADS)}, \
    {"write-index",       no_argument,       NULL, SAM_OPT_WRITE_INDEX}, \
    {"verbosity",         required_argument, NULL, SAM_OPT_VERBOSITY}

/*
 * Processes a standard "global" samtools long option.
 *
 * The 'c' value is the return value from a getopt_long() call.  It is checked
 * against the lopt[] array to find the corresponding value as this may have
 * been reassigned by the individual subcommand.
 *
 * Having found the entry, the corresponding long form is used to apply the
 * option, storing the setting in sam_global_args *ga.
 *
 * Returns 0 on success,
 *        -1 on failure.
 */
int parse_sam_global_opt(int c, const char *optarg, const struct option *lopt,
                         sam_global_args *ga);

/*
 * Report the usage for global options.
 *
 * This accepts a string with one character per SAM_OPT_GLOBAL_OPTIONS option
 * to determine which options need to be printed and how.
 * Each character should be one of:
 * '.'    No short option has been assigned. Use --long-opt only.
 * '-'    The long (and short) option has been disabled.
 * <c>    Otherwise the short option is character <c>.
 */
void sam_global_opt_help(FILE *fp, const char *shortopts);


void sam_global_args_init(sam_global_args *ga);
void sam_global_args_free(sam_global_args *ga);

#endif
