/*  test/fuzz/fuzz_phase.c -- Fuzz driver for samtools phase

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

int main_phase(int argc, char *argv[]);

/*
Fuzzes samtools phase.  See fuzz_common.h for build and run details.

-b writes <prefix>.0.bam, <prefix>.1.bam and <prefix>.chimera.bam; the
prefix is always fuzz_prefix, so these are fixed names removed after each
run.

Link phase.o.  Output goes to stdout.  No runtime files.

Inputs must be coordinate sorted.

Control bytes:
 c[0] bit  0    -A
      bit  1    -F
      bits 2-3  -k
      bits 4-5  -q
      bits 6-7  -Q
 c[1] bit  0    -b fuzz_prefix
      bits 1-7  unused
*/

static void args(fuzz_args *a, const uint8_t *c, char *in) {
    ARG("phase");
    if (c[0] & 0x01)   ARG("-A");
    if (c[0] & 0x02)   ARG("-F");
    CHOOSE("-k", (c[0] >> 2) & 3, NULL, "1", "13", "64");
    CHOOSE("-q", (c[0] >> 4) & 3, NULL, "0", "37", "1000");
    CHOOSE("-Q", (c[0] >> 6) & 3, NULL, "0", "13", "60");
    if (c[1] & 0x01)   ARG2("-b", fuzz_prefix);
    ENDOPTS();
    ARG(in);
}

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
    return fuzz_one_input(data, size, 2, FUZZ_SORTED, args, main_phase);
}
