/*  test/fuzz/fuzz_depth.c -- Fuzz driver for samtools depth

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

int main_depth(int argc, char *argv[]);

/*
Fuzzes samtools depth.  See fuzz_common.h for build and run details.

Link bam2depth.o.  Output goes to stdout.  Runtime file: regions.bed for -b.

-r and -X need an index, and -f a file of filenames, so aren't fuzzed.
One control bit passes the same input twice, for multi-file depth.

Inputs must be coordinate sorted with short @SQ lines (for -a/-aa).

Control bytes:
 c[0] bit  0    -a
      bit  1    -a (again; -aa)
      bit  2    -b regions.bed
      bit  3    -H
      bits 4-5  -l
      bit  6    -J
      bit  7    -s
 c[1] bits 0-1  -q
      bits 2-3  -Q
      bits 4-5  -g
      bits 6-7  -G
 c[2] bits 0-1  --incl-flags
      bits 2-3  --require-flags
      bit  4    -m 10 (deprecated, ignored)
      bit  5    a second copy of the input
      bits 6-7  unused
*/

static void args(fuzz_args *a, const uint8_t *c, char *in) {
    ARG("depth");

    if (c[0] & 0x01)   ARG("-a");
    if (c[0] & 0x02)   ARG("-a");
    if (c[0] & 0x04)   ARG2("-b", "regions.bed");
    if (c[0] & 0x08)   ARG("-H");
    CHOOSE("-l", (c[0] >> 4) & 3, NULL, "0", "20", "1000");
    if (c[0] & 0x40)   ARG("-J");
    if (c[0] & 0x80)   ARG("-s");

    CHOOSE("-q", c[1] & 3,        NULL, "0", "20", "99");
    CHOOSE("-Q", (c[1] >> 2) & 3, NULL, "0", "20", "255");
    CHOOSE("-g", (c[1] >> 4) & 3, NULL, "0x900", "UNMAP,DUP", "0xfff");
    CHOOSE("-G", (c[1] >> 6) & 3, NULL, "0", "0x400", "0xfff");

    CHOOSE("--incl-flags",    c[2] & 3,        NULL, "0x1", "READ2,REVERSE", "0");
    CHOOSE("--require-flags", (c[2] >> 2) & 3, NULL, "0x1", "PAIRED,PROPER_PAIR", "0xfff");
    if (c[2] & 0x10)   ARG2("-m", "10");
    // Threads make runs non-deterministic.
    //ARG2("-@", "2");

    ENDOPTS();
    ARG(in);
    if (c[2] & 0x20)   ARG(in);
}

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
    return fuzz_one_input(data, size, 3, FUZZ_SORTED | FUZZ_SHORT_REFS,
                          args, main_depth);
}
