/*  test/fuzz/fuzz_cram_size.c -- Fuzz driver for samtools cram-size

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

int main_cram_size(int argc, char *argv[]);

/*
Fuzzes samtools cram-size.  See fuzz_common.h for the input layout (file
byte, region bytes) and build/run details.

cram-size reads CRAM only: file-byte form 1 (CRAM re-encoded by htslib) or
3 (raw, with CRAM seeds) reach the interesting code; the other forms test
its rejection of other formats.

Link cram_size.o.  Output goes to stdout or /dev/null.  No runtime files.

Tool byte:
 c[0] bit  0    -v
      bit  1    -e
      bit  2    -o /dev/null
      bits 3-7  unused
*/

static void args(fuzz_args *a, const uint8_t *c, fuzz_file *f) {
    ARG("cram-size");
    if (c[0] & 0x01)   ARG("-v");
    if (c[0] & 0x02)   ARG("-e");
    if (c[0] & 0x04)   ARG2("-o", "/dev/null");
    ENDOPTS();
    ARG(f->fn);
}

static const fuzz_step steps[] = {{args, main_cram_size, NULL}};

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
    return fuzz_one_file_input(data, size, 1, FUZZ_PARSE, 0, steps, 1);
}
