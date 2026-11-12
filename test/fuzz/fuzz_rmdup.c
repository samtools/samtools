/*  test/fuzz/fuzz_rmdup.c -- Fuzz driver for samtools rmdup

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

int bam_rmdup(int argc, char *argv[]);

/*
Fuzzes samtools rmdup (deprecated, but still shipped).  See fuzz_common.h
for build and run details.

Link bam_rmdup.o, bam_rmdupse.o.  Output is "-" (stdout).  No runtime
files.

Inputs must be coordinate sorted.

Control byte:
 c[0] bit  0    -s
      bit  1    -S
      bits 2-7  unused
*/

static void args(fuzz_args *a, const uint8_t *c, char *in) {
    ARG("rmdup");
    if (c[0] & 0x01)   ARG("-s");
    if (c[0] & 0x02)   ARG("-S");
    ENDOPTS();
    ARG(in);
    ARG("-");
}

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
    return fuzz_one_input(data, size, 1, FUZZ_SORTED, args, bam_rmdup);
}
