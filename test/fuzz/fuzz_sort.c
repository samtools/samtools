/*  test/fuzz/fuzz_sort.c -- Fuzz driver for samtools sort

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

int bam_sort(int argc, char *argv[]);

/*
Fuzzes samtools sort.  See fuzz_common.h for build and run details.

Link bam_sort.o and the objects it needs (bam_sort.c also hosts merge).

Output goes to stdout.  Temporary files use -T fuzz_prefix and are removed
after each run.  Runtime file: ref.fa (for -I).

Note sort enforces a 1M minimum for -m, so fuzz-sized inputs never spill
to temporary files and the merge-from-disk path isn't reached.  Building
the fuzzer with that floor lowered would fix this, if bam_sort.c gains a
compile-time override.

Control bytes:
 c[0] bits 0-2  sort order: coordinate, -n, -N, -t RG, -t NM, -n -t XB,
                --template-coordinate, -M
      bit  3    -R
      bit  4    -H
      bits 5-6  -K
      bit  7    -I ref.fa (implies -M)
 c[1] bits 0-1  -w
      bits 2-3  -O
      bit  4    -u
      bits 5-6  -l
      bit  7    --no-PG
 c[2] bits 0-1  -m
      bits 2-7  unused
*/

static void args(fuzz_args *a, const uint8_t *c, char *in) {
    ARG("sort");
    ARG2("-T", fuzz_prefix);

    switch (c[0] & 7) {
    case 1: ARG("-n"); break;
    case 2: ARG("-N"); break;
    case 3: ARG2("-t", "RG"); break;
    case 4: ARG2("-t", "NM"); break;
    case 5: ARG("-n"); ARG2("-t", "XB"); break; // B-array tag
    case 6: ARG("--template-coordinate"); break;
    case 7: ARG("-M"); break;
    }
    if (c[0] & 0x08)   ARG("-R");
    if (c[0] & 0x10)   ARG("-H");
    CHOOSE("-K", (c[0] >> 5) & 3, NULL, "1", "15", "64");
    if (c[0] & 0x80)   ARG2("-I", "ref.fa");

    CHOOSE("-w", c[1] & 3,        NULL, "1", "10", "1000");
    CHOOSE("-O", (c[1] >> 2) & 3, NULL, "sam", "bam", "cram,no_ref");
    if (c[1] & 0x10)   ARG("-u");
    CHOOSE("-l", (c[1] >> 5) & 3, NULL, "0", "1", "9");
    if (c[1] & 0x80)   ARG("--no-PG");

    CHOOSE("-m", c[2] & 3, NULL, "1M", "1K", "1G");
    // Threads make runs non-deterministic.
    //ARG2("-@", "2");

    ENDOPTS();
    ARG(in);
}

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
    return fuzz_one_input(data, size, 3, FUZZ_PARSE, args, bam_sort);
}
