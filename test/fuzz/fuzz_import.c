/*  test/fuzz/fuzz_import.c -- Fuzz driver for samtools import

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

int main_import(int argc, char *argv[]);

/*
Fuzzes samtools import.  See fuzz_common.h for the input layout (file
byte, region bytes, extra payload) and build/run details.

The main input and the extra payload are both FASTQ (or FASTA), written
raw; together they provide READ1/READ2 and index-read files.

Link bam_import.o.  Output goes to stdout.  No runtime files.

Tool bytes:
 c[0] bits 0-2  inputs:
                  0 main
                  1 main extra (positional pair)
                  2 -s main
                  3 -0 main
                  4 -1 main -2 extra
                  5 -1 main -2 extra --i1 main --i2 extra
                  6 -1 main
                  7 main extra main (three positional)
      bit  3    -i
      bit  4    -U
      bit  5    -N
      bit  6    -u
      bit  7    --no-PG
 c[1] bits 0-2  read group:
                  0 none           1 -r ID:x
                  2 -r "@RG<TAB>ID:x<TAB>PL:ILLUMINA"
                  3 -r ID:x -r SM:y
                  4 -R rg9         5 -r SM:noid (no ID)
                  6 -r '@RG\tID:esc' (backslash escape)
                  7 -R "" (empty)
      bits 3-4  -T: none, *, empty, BC,QT
      bits 5-6  --order: none, ro, rn:12, rn:0
      bit  7    --UMI-tag OX
 c[2] bits 0-1  -O
      bits 2-3  --barcode-tag: none, BC, XZ
      bits 4-5  --quality-tag: none, QT, QX
      bits 6-7  unused
*/

static void args(fuzz_args *a, const uint8_t *c, fuzz_file *f) {
    const char *m = f->fn, *x = f->extra_raw_fn;
    int mode = c[0] & 7;

    ARG("import");
    switch (mode) {
    case 2: ARG2("-s", m); break;
    case 3: ARG2("-0", m); break;
    case 4: ARG2("-1", m); ARG2("-2", x); break;
    case 5: ARG2("-1", m); ARG2("-2", x); ARG2("--i1", m); ARG2("--i2", x); break;
    case 6: ARG2("-1", m); break;
    }
    if (c[0] & 0x08)   ARG("-i");
    if (c[0] & 0x10)   ARG("-U");
    if (c[0] & 0x20)   ARG("-N");
    if (c[0] & 0x40)   ARG("-u");
    if (c[0] & 0x80)   ARG("--no-PG");

    switch (c[1] & 7) {
    case 1: ARG2("-r", "ID:x"); break;
    case 2: ARG2("-r", "@RG\tID:x\tPL:ILLUMINA"); break;
    case 3: ARG2("-r", "ID:x"); ARG2("-r", "SM:y"); break;
    case 4: ARG2("-R", "rg9"); break;
    case 5: ARG2("-r", "SM:noid"); break;
    case 6: ARG2("-r", "@RG\\tID:esc"); break;
    case 7: ARG2("-R", ""); break;
    }
    CHOOSE("-T",      (c[1] >> 3) & 3, NULL, "*", "", "BC,QT");
    CHOOSE("--order", (c[1] >> 5) & 3, NULL, "ro", "rn:12", "rn:0");
    if (c[1] & 0x80)   ARG2("--UMI-tag", "OX");

    CHOOSE("-O", c[2] & 3, NULL, "sam", "bam", "cram,no_ref");
    CHOOSE("--barcode-tag", (c[2] >> 2) & 3, NULL, "BC", "XZ");
    CHOOSE("--quality-tag", (c[2] >> 4) & 3, NULL, "QT", "QX");

    ENDOPTS();
    switch (mode) {
    case 0: ARG(m); break;
    case 1: ARG(m); ARG(x); break;
    case 7: ARG(m); ARG(x); ARG(m); break;
    }
}

static const fuzz_step steps[] = {{args, main_import, NULL}};

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
    return fuzz_one_file_input(data, size, 3, FUZZ_ANYFMT, FUZZ_EXTRA,
                               steps, 1);
}
