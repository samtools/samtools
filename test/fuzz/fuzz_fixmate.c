/*  test/fuzz/fuzz_fixmate.c -- Fuzz driver for samtools fixmate

    Copyright (C) 2026 Genome Research Ltd.

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

#include "fuzz_common.h"

int bam_mating(int argc, char *argv[]);

/*
Fuzzes samtools fixmate

See fuzz_common.h for build and run details, including FUZZ_EXTRACT for
reproducing a crash outside the harness.

The leading control bytes of each input select the command line options.
The remainder is the alignment data, passed to samtools as a data: URL.

fixmate requires an output filename; we always give "-" (stdout), which
-close_fd_mask=1 discards.  The fuzzer writes no files and needs none.
CRAM output uses no_ref, so no reference is required.

The corpus should be name-collated, with pairs, singletons, secondary and
supplementary records, unmapped mates, records with and without MC/ms
tags, and base modification (MM/ML/MN) tags including hard-clipped
secondaries, which is what -M is for.

No validity of result is assumed.  We're simply checking for crashes or
illegal program behaviour.
*/

/*
 * Control bytes:
 *
 * c[0] bit  0    -r
 *      bit  1    -p
 *      bit  2    -c
 *      bit  3    -m
 *      bit  4    -M
 *      bit  5    -u
 *      bits 6-7  -O
 * c[1] bits 0-1  -z
 *      bit  2    --no-PG
 *      bits 3-7  unused
 */
static void args(fuzz_args *a, const uint8_t *c, char *in) {
    ARG("fixmate");

    if (c[0] & 0x01)   ARG("-r");
    if (c[0] & 0x02)   ARG("-p");
    if (c[0] & 0x04)   ARG("-c");
    if (c[0] & 0x08)   ARG("-m");
    if (c[0] & 0x10)   ARG("-M");
    if (c[0] & 0x20)   ARG("-u");
    CHOOSE("-O", (c[0] >> 6) & 3, NULL, "sam", "bam", "cram,no_ref");

    CHOOSE("-z", (c[1]) & 3,    NULL, "off", "on", "cigar,cigdup,cigarx,aux");
    if (c[1] & 0x04)   ARG("--no-PG");
    // Threads make runs non-deterministic.
    //if (c[1] & 0x08) ARG2("-@", "2");

    ENDOPTS();
    ARG(in);
    ARG("-");   // output to stdout; never write a file

}

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
    return fuzz_one_input(data, size, 2, FUZZ_PARSE, args, bam_mating);
}
