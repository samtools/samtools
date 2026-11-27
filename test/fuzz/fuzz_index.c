/*  test/fuzz/fuzz_index.c -- Fuzz driver for samtools index / idxstats

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

int bam_index(int argc, char *argv[]);
int bam_idxstats(int argc, char *argv[]);
int main_samview(int argc, char *argv[]);

/*
Fuzzes samtools index, then uses the index it built with samtools idxstats
and a region query (samtools view -c).  See fuzz_common.h for the input
layout (file byte, region bytes) and build/run details.

The driver writes the input to a file but does not index it; "samtools
index" does that, with fuzzed -b/-c/-m.  The index goes to the default name,
or to a fixed alternate name (via -o or the old positional form), in which
case later steps find it through "file##idx##index".

Link bam_index.o and sam_view.o (and their dependencies).  No runtime files.

Inputs must be coordinate sorted.

Tool byte:
 c[0] bits 0-1  index type: default, -b, -c, -c with -m
      bits 2-3  -m value (only with the last of those)
      bit  4    -o alternate-name
      bit  5    -M, indexing the file twice (overrides -o)
      bit  6    positional output name instead of -o (with bit 4)
      bit  7    unused
*/

static int alt_name(const uint8_t *c) {
    return (c[0] & 0x10) && !(c[0] & 0x20);
}

static void index_args(fuzz_args *a, const uint8_t *c, fuzz_file *f) {
    ARG("index");
    switch (c[0] & 3) {
    case 1: ARG("-b"); break;
    case 2: ARG("-c"); break;
    case 3: ARG("-c"); CHOOSE("-m", (c[0] >> 2) & 3, "14", "0", "1", "31"); break;
    }
    if (c[0] & 0x20) {
        ARG("-M");
        ENDOPTS();
        ARG(f->fn);
        ARG(f->fn);
        return;
    }
    if (alt_name(c) && !(c[0] & 0x40))
        ARG2("-o", f->idx_fn);
    ENDOPTS();
    ARG(f->fn);
    if (alt_name(c) && (c[0] & 0x40))
        ARG(f->idx_fn);
}

static void idxstats_args(fuzz_args *a, const uint8_t *c, fuzz_file *f) {
    ARG("idxstats");
    ENDOPTS();
    ARG(alt_name(c) ? f->idx_url : f->fn);
}

static void view_args(fuzz_args *a, const uint8_t *c, fuzz_file *f) {
    ARG("view");
    ARG("-c");
    ENDOPTS();
    ARG(alt_name(c) ? f->idx_url : f->fn);
    for (int i = 0; i < f->nreg; i++)
        ARG_DATA(f->reg[i].s);
}

static const fuzz_step steps[] = {
    {index_args,    bam_index},
    {idxstats_args, bam_idxstats},
    {view_args,     main_samview},
};

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
    return fuzz_one_file_input(data, size, 1, FUZZ_SORTED, 0, steps, 3);
}
