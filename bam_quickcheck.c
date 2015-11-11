/*  bam_quickcheck.c -- quickcheck subcommand.

    Copyright (C) 2015 Genome Research Ltd.

    Author: Joshua C. Randall <jcrandall@alum.mit.edu>

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

#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/bgzf.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

static void usage_quickcheck(FILE *write_to)
{
    fprintf(write_to,
"Usage: samtools quickcheck [options] <input> [...]\n"
"Options:\n"
"  -v              verbose output (repeat for more verbosity)\n"
"\n"
    );
}

int main_quickcheck(int argc, char** argv)
{
    int verbose = 0;
    hts_verbose = 0;

    const char* optstring = "v";
    int opt;
    while ((opt = getopt(argc, argv, optstring)) != -1) {
        switch (opt) {
        case 'v':
            verbose++;
            break;
        default:
            usage_quickcheck(stderr);
            return 1;
        }
    }

    argc -= optind;
    argv += optind;

    if (argc < 1) {
        usage_quickcheck(stdout);
        return 1;
    }

    if (verbose >= 2) {
        fprintf(stderr, "verbosity set to %d\n", verbose);
    }

    if (verbose >= 4) {
        hts_verbose = 3;
    }

    int ret = 0;
    int i;

    for (i = 0; i < argc; i++) {
        char* fn = argv[i];
        int file_state = 0;

        if (verbose >= 3) fprintf(stderr, "checking %s\n", fn);

        // attempt to open
        htsFile *hts_fp = hts_open(fn, "r");
        if (hts_fp == NULL) {
            if (verbose >= 2) fprintf(stderr, "%s could not be opened for reading\n", fn);
            file_state |= 2;
        }
        else {
            if (verbose >= 3) fprintf(stderr, "opened %s\n", fn);
            // make sure we have sequence data
            const htsFormat *fmt = hts_get_format(hts_fp);
            if (fmt->category != sequence_data ) {
                if (verbose >= 2) fprintf(stderr, "%s was not identified as sequence data\n", fn);
                file_state |= 4;
            }
            else {
                if (verbose >= 3) fprintf(stderr, "%s is sequence data\n", fn);
                // check header
                bam_hdr_t *header = sam_hdr_read(hts_fp);
                if (header->n_targets <= 0) {
                    if (verbose >= 2) fprintf(stderr, "%s had no targets in header\n", fn);
                    file_state |= 8;
                }
                else {
                    if (verbose >= 3) fprintf(stderr, "%s has %d targets in header\n", fn, header->n_targets);
                }

                // only check EOF on BAM for now
                // TODO implement and use hts_check_EOF() to include CRAM support
                if (fmt->format == bam) {
                    if (bgzf_check_EOF(hts_fp->fp.bgzf) <= 0) {
                        if (verbose >= 2) fprintf(stderr, "%s was missing EOF block\n", fn);
                        file_state |= 16;
                    }
                    else {
                        if (verbose >= 3) fprintf(stderr, "%s has good EOF block\n", fn);
                    }
                }
            }

            hts_close(hts_fp);
        }

        if (file_state > 0 && verbose >= 1) {
            fprintf(stdout, "%s\n", fn);
        }
        ret |= file_state;
    }

    return ret;
}
