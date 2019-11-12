/*  bam_quickcheck.c -- quickcheck subcommand.

    Copyright (C) 2015-2017 Genome Research Ltd.

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

#include <config.h>

#include <htslib/hts.h>
#include <htslib/sam.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

/* File status flags (zero means OK). It's possible for more than one to be
 * set on a single file.   The final exit status is the bitwise-or of the
 * status of all the files. */
#define QC_FAIL_OPEN     2
#define QC_NOT_SEQUENCE  4
#define QC_BAD_HEADER    8
#define QC_NO_EOF_BLOCK 16
#define QC_FAIL_CLOSE   32

static void usage_quickcheck(FILE *write_to)
{
    fprintf(write_to,
"Usage: samtools quickcheck [options] <input> [...]\n"
"Options:\n"
"  -v              verbose output (repeat for more verbosity)\n"
"  -q              suppress warning messages\n"
"  -u              unmapped input (do not require targets in header)\n"
"\n"
"Notes:\n"
"\n"
"1. By default quickcheck will emit a warning message if and only if a file\n"
"   fails the checks, in which case the exit status is non-zero.  Under normal\n"
"   behaviour with valid data it will be silent and has a zero exit status.\n"
"   The warning messages are purely for manual inspection and should not be \n"
"   parsed by scripts.\n"
"\n"
"2. In order to use this command programmatically, you should check its exit\n"
"   status.  One way to use quickcheck might be as a check that all BAM files in\n"
"   a directory are okay:\n"
"\n"
"\tsamtools quickcheck *.bam && echo 'all ok' \\\n"
"\t   || echo 'fail!'\n"
"\n"
"   The first level of verbosity lists only files that fail to stdout.\n"
"   To obtain a parsable list of files that have failed, use this option:\n"
"\n"
"\tsamtools quickcheck -qv *.bam > bad_bams.fofn \\\n"
"\t   && echo 'all ok' \\\n"
"\t   || echo 'some files failed check, see bad_bams.fofn'\n"
    );
}

#define QC_ERR(state, v, msg, arg1)                                     \
    file_state |= (state);                                              \
    if (!quiet || verbose >= (v)) fprintf(stderr, (msg), (arg1))

int main_quickcheck(int argc, char** argv)
{
    int verbose = 0, quiet = 0, unmapped = 0;
    hts_verbose = 0;

    const char* optstring = "vqu";
    int opt;
    while ((opt = getopt(argc, argv, optstring)) != -1) {
        switch (opt) {
        case 'u':
            unmapped = 1;
            break;
        case 'v':
            verbose++;
            break;
        case 'q':
            quiet = 1;
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
            QC_ERR(QC_FAIL_OPEN, 2, "%s could not be opened for reading.\n", fn);
        }
        else {
            if (verbose >= 3) fprintf(stderr, "opened %s\n", fn);
            // make sure we have sequence data
            const htsFormat *fmt = hts_get_format(hts_fp);
            if (fmt->category != sequence_data ) {
                QC_ERR(QC_NOT_SEQUENCE, 2, "%s was not identified as sequence data.\n", fn);
            }
            else {
                if (verbose >= 3) fprintf(stderr, "%s is sequence data\n", fn);
                // check header
                sam_hdr_t *header = sam_hdr_read(hts_fp);
                if (header == NULL) {
                    QC_ERR(QC_BAD_HEADER, 2, "%s caused an error whilst reading its header.\n", fn);
                } else {
                    if (!unmapped && sam_hdr_nref(header) <= 0) {
                        QC_ERR(QC_BAD_HEADER, 2, "%s had no targets in header.\n", fn);
                    }
                    else {
                        if (verbose >= 3) fprintf(stderr, "%s has %d targets in header.\n", fn, sam_hdr_nref(header));
                    }
                    sam_hdr_destroy(header);
                }
            }
            // check EOF on formats that support this
            int ret;
            if ((ret = hts_check_EOF(hts_fp)) < 0) {
                QC_ERR(QC_NO_EOF_BLOCK, 2, "%s caused an error whilst checking for EOF block.\n", fn);
           }
            else {
                switch (ret) {
                    case 0:
                        QC_ERR(QC_NO_EOF_BLOCK, 2, "%s was missing EOF block when one should be present.\n", fn);
                        break;
                    case 1:
                        if (verbose >= 3) fprintf(stderr, "%s has good EOF block.\n", fn);
                        break;
                    case 2:
                        if (verbose >= 3) fprintf(stderr, "%s cannot be checked for EOF block as it is not seekable.\n", fn);
                        break;
                    case 3:
                        if (verbose >= 3) fprintf(stderr, "%s cannot be checked for EOF block because its filetype does not contain one.\n", fn);
                        break;
                }
            }

            if (hts_close(hts_fp) < 0) {
                QC_ERR(QC_FAIL_CLOSE, 2, "%s did not close cleanly.\n", fn);
            }
        }

        if (file_state > 0 && verbose >= 1) {
            fprintf(stdout, "%s\n", fn);
        }
        ret |= file_state;
    }

    return ret;
}
