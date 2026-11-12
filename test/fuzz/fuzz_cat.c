/*  test/fuzz/fuzz_cat.c -- Fuzz driver for samtools cat

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

int main_cat(int argc, char *argv[]);

/*
Fuzzes samtools cat, which concatenates BAM or CRAM files without decoding
them.  See fuzz_common.h for the input layout (file byte, region bytes,
extra payload) and build/run details.

The extra payload is a second alignment file, written in the same format
as the first (BGZF SAM, from file byte form 2, exercises cat's rejection
of unsupported formats), and may also serve as the -h header file.

Link bam_cat.o.  Output is "-o -".

Tool byte:
 c[0] bits 0-1  inputs: main + extra; main twice; main + extra + main;
                -b <file of filenames> (main only)
      bit  2    -h <extra payload as a header file>
      bit  3    --no-PG
      bits 4-7  unused
*/

static void args(fuzz_args *a, const uint8_t *c, fuzz_file *f) {
    ARG("cat");
    ARG2("-o", "-");
    if (c[0] & 0x04)   ARG2("-h", f->extra_raw_fn);
    if (c[0] & 0x08)   ARG("--no-PG");
    if ((c[0] & 3) == 3)
        ARG2("-b", f->list_fn);

    ENDOPTS();
    switch (c[0] & 3) {
    case 0: ARG(f->fn); ARG(f->extra_fn); break;
    case 1: ARG(f->fn); ARG(f->fn); break;
    case 2: ARG(f->fn); ARG(f->extra_fn); ARG(f->fn); break;
    }
}

static const fuzz_step steps[] = {{args, main_cat}};

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
    return fuzz_one_file_input(data, size, 1, FUZZ_PARSE, FUZZ_EXTRA,
                               steps, 1);
}
