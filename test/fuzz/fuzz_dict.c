/*  test/fuzz/fuzz_dict.c -- Fuzz driver for samtools dict

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

int dict_main(int argc, char *argv[]);

/*
Fuzzes samtools dict.  See fuzz_common.h for the input layout (file byte,
region bytes, extra payload) and build/run details.

The input is FASTA (plain or BGZF), written raw.  The extra payload is the
bwa-style alt-contig file for --alt.

Link dict.o.  Output goes to stdout.  No runtime files.

Tool bytes:
 c[0] bits 0-1  -a: none, GRCh38, empty, with a space
      bit  2    -A
      bit  3    -H
      bits 4-5  -s: none, "Homo sapiens", empty, x
      bits 6-7  -u: none, file URI, empty, http URI
 c[1] bit  0    --alt <extra payload>
      bits 1-7  unused
*/

static void args(fuzz_args *a, const uint8_t *c, fuzz_file *f) {
    ARG("dict");
    CHOOSE("-a", c[0] & 3,        NULL, "GRCh38", "", "a b");
    if (c[0] & 0x04)   ARG("-A");
    if (c[0] & 0x08)   ARG("-H");
    CHOOSE("-s", (c[0] >> 4) & 3, NULL, "Homo sapiens", "", "x");
    CHOOSE("-u", (c[0] >> 6) & 3, NULL, "file:///ref.fa", "", "http://x/y");
    if (c[1] & 0x01)   ARG2("--alt", f->extra_raw_fn);
    ENDOPTS();
    ARG(f->fn);
}

static const fuzz_step steps[] = {{args, dict_main, NULL}};

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
    return fuzz_one_file_input(data, size, 2, FUZZ_ANYFMT, FUZZ_EXTRA,
                               steps, 1);
}
