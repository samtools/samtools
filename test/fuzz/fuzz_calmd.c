/*  test/fuzz/fuzz_calmd.c -- Fuzz driver for samtools calmd

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

int bam_fillmd(int argc, char *argv[]);

/*
Fuzzes samtools calmd

See fuzz_common.h for build and run details, including FUZZ_EXTRACT for
reproducing a crash outside the harness.

The leading control bytes of each input select the command line options.
The remainder is the alignment data, passed to samtools as a data: URL.

calmd has no -o option and always writes to stdout; -close_fd_mask=1
discards it.  The fuzzer writes no files, provided the fixed reference
files below have their indices already built (otherwise htslib will
write a .fai / .gzi on first use).

Fixed reference files, which must exist in the working directory:
    ref.fa,       ref.fa.fai                  the reference for the corpus
    ref.fa.gz,    ref.fa.gz.fai, ref.fa.gz.gzi    bgzip -i of ref.fa
    ref_short.fa, ref_short.fa.fai            same names as ref.fa, but
                                              sequences truncated (and
                                              ideally one name missing),
                                              so reads run off the end
The corpus should be aligned against ref.fa, and include records with
existing MD/NM tags (correct and incorrect), "=" bases, indels, clips,
"*" SEQ/QUAL, and unsorted orderings, which use calmd's caching mode.

No validity of result is assumed.  We're simply checking for crashes or
illegal program behaviour.
*/

/*
 * Control bytes:
 *
 * c[0] bit  0    -A
 *      bit  1    -e
 *      bit  2    -u
 *      bit  3    -b
 *      bits 4-5  -C
 *      bit  6    -r
 *      bit  7    -E
 * c[1] bits 0-1  reference file
 *      bit  2    --no-PG
 *      bits 3-7  unused
 */
static void args(fuzz_args *a, const uint8_t *c, char *in) {
    ARG("calmd");

    if (c[0] & 0x01)   ARG("-A");
    if (c[0] & 0x02)   ARG("-e");
    if (c[0] & 0x04)   ARG("-u");
    if (c[0] & 0x08)   ARG("-b");
    CHOOSE("-C", (c[0] >> 4) & 3, NULL, "0", "50", "255");
    if (c[0] & 0x40)   ARG("-r");
    if (c[0] & 0x80)   ARG("-E");

    if (c[1] & 0x04)   ARG("--no-PG");
    // Threads make runs non-deterministic.
    //if (c[1] & 0x08) ARG2("-@", "2");

    ENDOPTS();
    ARG(in);

    // ref.fa is listed twice so the matching reference is used half the time
    static const char *ref[4] = {"ref.fa", "ref.fa.gz", "ref_short.fa", "ref.fa"};
    ARG(ref[c[1] & 3]);

}

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
    return fuzz_one_input(data, size, 2, FUZZ_PARSE, args, bam_fillmd);
}
