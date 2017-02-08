/*  sam_opts.c -- utilities to aid parsing common command line options.

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

#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sam_opts.h"

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
                         sam_global_args *ga) {
    int r = 0;

    while (lopt->name) {
        if (c != lopt->val) {
            lopt++;
            continue;
        }

        if (strcmp(lopt->name, "input-fmt") == 0) {
            r = hts_parse_format(&ga->in, optarg);
            break;
        } else if (strcmp(lopt->name, "input-fmt-option") == 0) {
            r = hts_opt_add((hts_opt **)&ga->in.specific, optarg);
            break;
        } else if (strcmp(lopt->name, "output-fmt") == 0) {
            r = hts_parse_format(&ga->out, optarg);
            break;
        } else if (strcmp(lopt->name, "output-fmt-option") == 0) {
            r = hts_opt_add((hts_opt **)&ga->out.specific, optarg);
            break;
        } else if (strcmp(lopt->name, "reference") == 0) {
            char *ref = malloc(10 + strlen(optarg) + 1);
            sprintf(ref, "reference=%s", optarg);
            ga->reference = strdup(optarg);
            r  = hts_opt_add((hts_opt **)&ga->in.specific, ref);
            r |= hts_opt_add((hts_opt **)&ga->out.specific, ref);
            free(ref);
            break;
        } else if (strcmp(lopt->name, "threads") == 0) {
            ga->nthreads = atoi(optarg);
            break;
//      } else if (strcmp(lopt->name, "verbose") == 0) {
//          ga->verbosity++;
//          break;
        }
    }

    if (!lopt->name) {
        fprintf(stderr, "Unexpected global option: %s\n", lopt->name);
        return -1;
    }

    return r;
}

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
void sam_global_opt_help(FILE *fp, const char *shortopts) {
    int i = 0;

    static const struct option lopts[] = {
        SAM_OPT_GLOBAL_OPTIONS(0,0,0,0,0,0),
        { NULL, 0, NULL, 0 }
    };

    for (i = 0; shortopts && shortopts[i] && lopts[i].name; i++) {
        if (shortopts[i] == '-')
            continue;

        if (shortopts[i] == '.')
            fprintf(fp, "      --");
        else
            fprintf(fp, "  -%c, --", shortopts[i]);

        if (strcmp(lopts[i].name, "input-fmt") == 0)
            fprintf(fp,"input-fmt FORMAT[,OPT[=VAL]]...\n"
                    "               Specify input format (SAM, BAM, CRAM)\n");
        else if (strcmp(lopts[i].name, "input-fmt-option") == 0)
            fprintf(fp,"input-fmt-option OPT[=VAL]\n"
                    "               Specify a single input file format option in the form\n"
                    "               of OPTION or OPTION=VALUE\n");
        else if (strcmp(lopts[i].name, "output-fmt") == 0)
            fprintf(fp,"output-fmt FORMAT[,OPT[=VAL]]...\n"
                    "               Specify output format (SAM, BAM, CRAM)\n");
        else if (strcmp(lopts[i].name, "output-fmt-option") == 0)
            fprintf(fp,"output-fmt-option OPT[=VAL]\n"
                    "               Specify a single output file format option in the form\n"
                    "               of OPTION or OPTION=VALUE\n");
        else if (strcmp(lopts[i].name, "reference") == 0)
            fprintf(fp,"reference FILE\n"
                    "               Reference sequence FASTA FILE [null]\n");
        else if (strcmp(lopts[i].name, "threads") == 0)
            fprintf(fp,"threads INT\n"
                    "               Number of additional threads to use [0]\n");
//      else if (strcmp(lopts[i].name, "verbose") == 0)
//          fprintf(fp,"verbose\n"
//                  "               Increment level of verbosity\n");
    }
}

void sam_global_args_init(sam_global_args *ga) {
    if (!ga)
        return;

    memset(ga, 0, sizeof(*ga));
}

void sam_global_args_free(sam_global_args *ga) {
    if (ga->in.specific)
        hts_opt_free(ga->in.specific);

    if (ga->out.specific)
        hts_opt_free(ga->out.specific);

    if (ga->reference)
        free(ga->reference);
}
