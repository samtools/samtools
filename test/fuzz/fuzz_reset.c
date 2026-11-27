/*  test/fuzz/fuzz_reset.c -- Fuzz driver for samtools reset

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

int main_reset(int argc, char *argv[]);

/*
Fuzzes samtools reset.  See fuzz_common.h for build and run details.

Link reset.o.  Output goes to stdout.  No runtime files.

Control bytes:
 c[0] bits 0-1  -x: none, MD,NM, ^ (remove all), ^RG,BC (keep only)
      bits 2-3  --keep-tag: none, RG,BC, empty, XA
      bits 4-5  --reject-PG: none, bwa, absent, empty
      bit  6    --no-RG
      bit  7    --no-PG
 c[1] bit  0    --dupflag
      bits 1-2  -O
      bit  3    -x given twice (-x XA -x XB)
      bits 4-7  unused
*/

static void args(fuzz_args *a, const uint8_t *c, char *in) {
    ARG("reset");
    CHOOSE("-x",          c[0] & 3,        NULL, "MD,NM", "^", "^RG,BC");
    CHOOSE("--keep-tag",  (c[0] >> 2) & 3, NULL, "RG,BC", "", "XA");
    CHOOSE("--reject-PG", (c[0] >> 4) & 3, NULL, "bwa", "absent", "");
    if (c[0] & 0x40)   ARG("--no-RG");
    if (c[0] & 0x80)   ARG("--no-PG");

    if (c[1] & 0x01)   ARG("--dupflag");
    CHOOSE("-O", (c[1] >> 1) & 3, NULL, "sam", "bam", "cram,no_ref");
    if (c[1] & 0x08) { ARG2("-x", "XA"); ARG2("-x", "XB"); }

    ENDOPTS();
    ARG(in);
}

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
    return fuzz_one_input(data, size, 2, FUZZ_PARSE, args, main_reset);
}
