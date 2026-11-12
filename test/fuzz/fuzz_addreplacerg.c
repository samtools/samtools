/*  test/fuzz/fuzz_addreplacerg.c -- Fuzz driver for samtools addreplacerg

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

int main_addreplacerg(int argc, char *argv[]);

/*
Fuzzes samtools addreplacerg.  See fuzz_common.h for build and run details.

Link bam_addrprg.o.  Output is "-o -".  No runtime files.

Control bytes:
 c[0] bits 0-2  read group:
                  0 -r ID:fish
                  1 -r "@RG<TAB>ID:fish<TAB>SM:x" (literal tabs)
                  2 -r '@RG\tID:fish\tSM:x' (backslash escapes)
                  3 -r ID:fish -r LB:1334 -r SM:alpha
                  4 -R rg1 (exists in the seed headers)
                  5 -R absent
                  6 -r SM:x (no ID)
                  7 -r ID:rg1 (duplicates an existing ID)
      bits 3-4  -m: none, overwrite_all, orphan_only, invalid
      bit  5    -w
      bit  6    --no-PG
      bit  7    unused
 c[1] bits 0-1  -O
      bits 2-3  -l
      bits 4-7  unused
*/

static void args(fuzz_args *a, const uint8_t *c, char *in) {
    ARG("addreplacerg");
    ARG2("-o", "-");

    switch (c[0] & 7) {
    case 0: ARG2("-r", "ID:fish"); break;
    case 1: ARG2("-r", "@RG\tID:fish\tSM:x"); break;
    case 2: ARG2("-r", "@RG\\tID:fish\\tSM:x"); break;
    case 3: ARG2("-r", "ID:fish"); ARG2("-r", "LB:1334"); ARG2("-r", "SM:alpha"); break;
    case 4: ARG2("-R", "rg1"); break;
    case 5: ARG2("-R", "absent"); break;
    case 6: ARG2("-r", "SM:x"); break;
    case 7: ARG2("-r", "ID:rg1"); break;
    }
    CHOOSE("-m", (c[0] >> 3) & 3, NULL, "overwrite_all", "orphan_only", "bogus");
    if (c[0] & 0x20)   ARG("-w");
    if (c[0] & 0x40)   ARG("--no-PG");

    CHOOSE("-O", c[1] & 3,        NULL, "sam", "bam", "cram,no_ref");
    CHOOSE("-l", (c[1] >> 2) & 3, NULL, "0", "1", "9");

    ENDOPTS();
    ARG(in);
}

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
    return fuzz_one_input(data, size, 2, FUZZ_PARSE, args, main_addreplacerg);
}
