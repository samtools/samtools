/*  test/fuzz/fuzz_merge.c -- Fuzz driver for samtools merge

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

int bam_merge(int argc, char *argv[]);

/*
Fuzzes samtools merge.  See fuzz_common.h for the input layout (file byte,
region bytes, extra payload) and build/run details.

The extra payload is a second alignment file (merged with the first, which
exercises header merging: @SQ mismatches, @RG/@PG ID collisions, ...), and
may also serve as the -h header file.  Both inputs are indexed, for -R.

Link bam_sort.o (which hosts merge) and its dependencies.  Output is "-o -".

Tool bytes:
 c[0] bits 0-2  sort order: coordinate, -n, -N, -t RG, -n -t XB,
                --template-coordinate, -t NM, coordinate
      bit  3    -r
      bit  4    -c
      bit  5    -p
      bit  6    -u
      bit  7    -1
 c[1] bits 0-1  inputs: main + extra file; main + extra data: URL;
                main twice; main + extra + main
      bit  2    -h <extra payload as a header file>
      bit  3    -R <first region>
      bit  4    -L <generated BED>
      bit  5    -b <file of filenames> as well
      bit  6    -s 42
      bit  7    --no-PG
 c[2] bits 0-1  -O
      bits 2-3  -l
      bits 4-7  unused
*/

static void args(fuzz_args *a, const uint8_t *c, fuzz_file *f) {
    ARG("merge");
    ARG2("-o", "-");

    switch (c[0] & 7) {
    case 1: ARG("-n"); break;
    case 2: ARG("-N"); break;
    case 3: ARG2("-t", "RG"); break;
    case 4: ARG("-n"); ARG2("-t", "XB"); break;
    case 5: ARG("--template-coordinate"); break;
    case 6: ARG2("-t", "NM"); break;
    }
    if (c[0] & 0x08)   ARG("-r");
    if (c[0] & 0x10)   ARG("-c");
    if (c[0] & 0x20)   ARG("-p");
    if (c[0] & 0x40)   ARG("-u");
    if (c[0] & 0x80)   ARG("-1");

    if (c[1] & 0x04)   ARG2("-h", f->extra_raw_fn);
    if (c[1] & 0x08)   ARG_OPTDATA("-R", f->reg[0].s);
    if (c[1] & 0x10)   ARG2("-L", f->bed_fn);
    if (c[1] & 0x20)   ARG2("-b", f->list_fn);
    if (c[1] & 0x40)   ARG2("-s", "42");
    if (c[1] & 0x80)   ARG("--no-PG");

    CHOOSE("-O", c[2] & 3,        NULL, "sam", "bam", "cram,no_ref");
    CHOOSE("-l", (c[2] >> 2) & 3, NULL, "0", "1", "9");

    ENDOPTS();
    ARG(f->in);
    switch (c[1] & 3) {
    case 0: ARG(f->extra_fn); break;   // missing if extra didn't parse
    case 1: ARG_DATA(f->extra_url); break;
    case 2: ARG(f->in); break;
    case 3: ARG(f->extra_fn); ARG(f->in); break;
    }
}

static const fuzz_step steps[] = {{args, bam_merge}};

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
    return fuzz_one_file_input(data, size, 3, FUZZ_PARSE,
                               FUZZ_INDEX | FUZZ_EXTRA, steps, 1);
}
