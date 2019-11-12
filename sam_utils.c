/*  sam_utils.c -- various utilities internal to samtools.

    Copyright (C) 2014-2016, 2018, 2019 Genome Research Ltd.

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
#include <stdlib.h>

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

void check_sam_close(const char *subcmd, samFile *fp, const char *fname, const char *null_fname, int *retp)
{
    int r = sam_close(fp);
    if (r >= 0) return;

    // TODO Need error infrastructure so we can print a message instead of r
    if (fname) print_error(subcmd, "error closing \"%s\": %d", fname, r);
    else print_error(subcmd, "error closing %s: %d", null_fname, r);

    *retp = EXIT_FAILURE;
}

/* Pick an index suffix based on the output file descriptor type. */
static char *idx_suffix(htsFile *fp) {
    switch (fp->format.format) {
    case sam:
    case bam:
        // Tough cheese if you wanted bai!
        // New feature => mandatory new index too, for simplicity of CLI.
        return "csi";

    case cram:
        return "crai";

    default:
        return NULL;
    }
}

/*
 * Utility function to add an index to a file we've opened for write.
 * NB: Call this after writing the header and before writing sequences.
 *
 * The returned index filename should be freed by the caller, but only
 * after sam_idx_save has been called.
 *
 * Returns index filename on success,
 *         NULL on failure.
 */
char *auto_index(htsFile *fp, const char *fn, bam_hdr_t *header) {
    char *fn_idx;
    int min_shift = 14; /* CSI */
    if (!fn || !*fn || strcmp(fn, "-") == 0)
        return NULL;

    char *delim = strstr(fn, HTS_IDX_DELIM);
    if (delim != NULL) {
        delim += strlen(HTS_IDX_DELIM);

        fn_idx = strdup(delim);
        if (!fn_idx)
            return NULL;

        size_t l = strlen(fn_idx);
        if (l >= 4 && strcmp(fn_idx + l - 4, ".bai") == 0)
            min_shift = 0;
    } else {
        char *suffix = idx_suffix(fp);
        if (!suffix)
            return NULL;

        fn_idx = malloc(strlen(fn)+6);
        if (!fn_idx)
            return NULL;

        sprintf(fn_idx, "%s.%s", fn, suffix);
    }

    if (sam_idx_init(fp, header, min_shift, fn_idx) < 0) {
        print_error_errno("auto_index", "failed to open index \"%s\" for writing", fn_idx);
        free(fn_idx);
        return NULL;
    }

    return fn_idx;
}
