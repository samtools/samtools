/*  sam_opts.h -- utilities to aid parsing common command line options.

    Copyright (C) 2015 Genome Research Ltd.

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

#ifndef BAM_OPTS_H
#define BAM_OPTS_H

#include <limits.h>
#include <getopt.h>
#include <htslib/hts.h>

typedef struct sam_global_args {
    htsFileOpts in;
    htsFileOpts out;
    char *reference; // TO DO
    int verbosity;
} sam_global_args;

#define SAM_GLOBAL_ARGS_INIT {HTS_FILE_OPTS_INIT, HTS_FILE_OPTS_INIT, NULL, 0}

enum {
    SAM_OPT_INPUT_FMT = CHAR_MAX+1,
    SAM_OPT_INPUT_FMT_OPTION,
    SAM_OPT_OUTPUT_FMT,
    SAM_OPT_OUTPUT_FMT_OPTION,
    SAM_OPT_VERBOSE
};

// Use this within an existing struct option lopts[] = {...}, used where
// the subcommmand is already using some long options.
#define SAM_GLOBAL_LOPTS \
    {"input-fmt",         required_argument, NULL, SAM_OPT_INPUT_FMT}, \
    {"input-fmt-option",  required_argument, NULL, SAM_OPT_INPUT_FMT_OPTION}, \
    {"output-fmt",        required_argument, NULL, SAM_OPT_OUTPUT_FMT}, \
    {"output-fmt-option", required_argument, NULL, SAM_OPT_OUTPUT_FMT_OPTION}, \
    {"verbose",           no_argument,       NULL, SAM_OPT_VERBOSE}

// Use this to completely assign all long options, used where the subcommand
// doesn't have any long options of its own.
//
// Eg: struct option lopts[] = SAM_GLOBAL_LOPTS_INIT;
#define SAM_GLOBAL_LOPTS_INIT {SAM_GLOBAL_LOPTS, {NULL,0,NULL,0}}


/*
 * Assign a short option to each of the long options listed above.
 *
 * There should be one character per option. This will be one of:
 * '.'    No short option has been assigned. Use --long-opt only.
 * '-'    The long (and short) option has been disabled.
 * <c>    Otherwise the short option is character <c>.
 */
void assign_short_opts(struct option lopts[], const char *shortopts);


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
int parse_sam_global_opt(int c, char *optarg, struct option *lopt, 
			 sam_global_args *ga);

/*
 * Report the usage for global options.
 *
 * This accepts the same shortopts string as used by assign_short_opts()
 * to determine which options need to be printed and how.
 */
void sam_global_opt_help(FILE *fp, char *shortopts);


void sam_global_args_init(sam_global_args *ga);
void sam_global_args_free(sam_global_args *ga);

#endif
