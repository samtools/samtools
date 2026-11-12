/*  test/fuzz/fuzz_samples.c -- Fuzz driver for samtools samples

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

int main_samples(int argc, char *argv[]);

/*
Fuzzes samtools samples.  See fuzz_common.h for the input layout (file
byte, region bytes) and build/run details.

Link bam_samples.o.  Output goes to stdout.  Runtime file: ref.fa (-f).

Tool byte:
 c[0] bit  0    -h
      bit  1    -i
      bits 2-3  -T: none, SM, LB, ID
      bit  4    -f ref.fa
      bit  5    -X (alignment file then index file)
      bit  6    a second copy of the input
      bit  7    unused
*/

static void args(fuzz_args *a, const uint8_t *c, fuzz_file *f) {
    int x = c[0] & 0x20, two = c[0] & 0x40;
    ARG("samples");
    if (c[0] & 0x01)   ARG("-h");
    if (c[0] & 0x02)   ARG("-i");
    CHOOSE("-T", (c[0] >> 2) & 3, NULL, "SM", "LB", "ID");
    if (c[0] & 0x10)   ARG2("-f", "ref.fa");
    if (x)             ARG("-X");
    ENDOPTS();
    for (int i = 0; i < (two ? 2 : 1); i++) {
        ARG(x ? f->fn : f->in);
        if (x)
            ARG(f->idx_path);
    }
}

static const fuzz_step steps[] = {{args, main_samples, NULL}};

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
    return fuzz_one_file_input(data, size, 1, FUZZ_PARSE, FUZZ_INDEX,
                               steps, 1);
}
