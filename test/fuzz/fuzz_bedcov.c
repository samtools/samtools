/*  test/fuzz/fuzz_bedcov.c -- Fuzz driver for samtools bedcov

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

int main_bedcov(int argc, char *argv[]);

/*
Fuzzes samtools bedcov.  See fuzz_common.h for the input layout (file
byte, region bytes) and build/run details.

The BED file is the driver's generated one (regions derived from the
fuzzed header and region bytes), so every input is self-contained.

Link bedcov.o.  Output goes to stdout.  No runtime files.

Inputs must be coordinate sorted.

Tool bytes:
 c[0] bits 0-1  -Q
      bits 2-3  -g
      bits 4-5  -G
      bit  6    -j
      bit  7    -c
 c[1] bits 0-1  -d
      bit  2    -H
      bit  3    unused
      bit  4    a second copy of the input
      bit  5    -X, giving the index path explicitly
      bits 6-7  unused
*/

static void args(fuzz_args *a, const uint8_t *c, fuzz_file *f) {
    ARG("bedcov");

    CHOOSE("-Q", c[0] & 3,        NULL, "0", "20", "255");
    CHOOSE("-g", (c[0] >> 2) & 3, NULL, "0x900", "UNMAP,DUP", "0xfff");
    CHOOSE("-G", (c[0] >> 4) & 3, NULL, "0", "0x400", "0xfff");
    if (c[0] & 0x40)   ARG("-j");
    if (c[0] & 0x80)   ARG("-c");

    CHOOSE("-d", c[1] & 3, NULL, "0", "1", "100");
    if (c[1] & 0x04)   ARG("-H");
    int x = c[1] & 0x20, two = c[1] & 0x10;
    if (x)             ARG("-X");

    ENDOPTS();
    ARG(f->bed_fn);
    ARG(x ? f->fn : f->in);
    if (two)           ARG(x ? f->fn : f->in);
    if (x) {
        ARG(f->idx_path);
        if (two)       ARG(f->idx_path);
    }
}

static const fuzz_step steps[] = {{args, main_bedcov}};

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
    return fuzz_one_file_input(data, size, 2, FUZZ_SORTED, 1, steps, 1);
}
