/*  test/fuzz/fuzz_depad.c -- Fuzz driver for samtools depad

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

int main_pad2unpad(int argc, char *argv[]);

/*
Fuzzes samtools depad (pad2unpad).  See fuzz_common.h for build and run
details.

Seeds should include padded alignments (P CIGAR operations).

Link padding.o.  Output goes to stdout.  Runtime file: ref.fa (for -T).

Control byte:
 c[0] bit  0    -s (SAM output)
      bit  1    -S
      bit  2    -C (CRAM output)
      bit  3    -u
      bit  4    -1
      bit  5    -T ref.fa
      bit  6    --no-PG
      bit  7    unused
*/

static void args(fuzz_args *a, const uint8_t *c, char *in) {
    ARG("depad");
    if (c[0] & 0x01)   ARG("-s");
    if (c[0] & 0x02)   ARG("-S");
    if (c[0] & 0x04)   ARG("-C");
    if (c[0] & 0x08)   ARG("-u");
    if (c[0] & 0x10)   ARG("-1");
    if (c[0] & 0x20)   ARG2("-T", "ref.fa");
    if (c[0] & 0x40)   ARG("--no-PG");
    ENDOPTS();
    ARG(in);
}

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
    return fuzz_one_input(data, size, 1, FUZZ_PARSE, args, main_pad2unpad);
}
