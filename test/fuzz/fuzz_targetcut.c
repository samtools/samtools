/*  test/fuzz/fuzz_targetcut.c -- Fuzz driver for samtools targetcut

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

int main_cut_target(int argc, char *argv[]);

/*
Fuzzes samtools targetcut.  See fuzz_common.h for build and run details.

Link cut_target.o.  Output goes to stdout.  Runtime file: ref.fa (-f).

Inputs must be coordinate sorted.

Control bytes:
 c[0] bits 0-1  -Q
      bits 2-3  -i
      bits 4-5  -0
      bits 6-7  -1
 c[1] bits 0-1  -2
      bit  2    -f ref.fa
      bits 3-7  unused
*/

static void args(fuzz_args *a, const uint8_t *c, char *in) {
    ARG("targetcut");
    CHOOSE("-Q", c[0] & 3,        NULL, "0", "13", "60");
    CHOOSE("-i", (c[0] >> 2) & 3, NULL, "0", "5", "100");
    CHOOSE("-0", (c[0] >> 4) & 3, NULL, "0", "-6", "20");
    CHOOSE("-1", (c[0] >> 6) & 3, NULL, "0", "-6", "40");
    CHOOSE("-2", c[1] & 3,        NULL, "0", "-6", "60");
    if (c[1] & 0x04)   ARG2("-f", "ref.fa");
    ENDOPTS();
    ARG(in);
}

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
    return fuzz_one_input(data, size, 2, FUZZ_SORTED, args, main_cut_target);
}
