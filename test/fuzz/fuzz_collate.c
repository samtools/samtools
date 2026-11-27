/*  test/fuzz/fuzz_collate.c -- Fuzz driver for samtools collate

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

int main_bamshuf(int argc, char *argv[]);

/*
Fuzzes samtools collate.  See fuzz_common.h for build and run details.

Temporary files use -T fuzz_prefix, and are removed after each run.  Output
is "-o -", or with the old positional form the fixed name fuzz_prefix.bam
(also removed).

Link bamshuf.o.  No runtime files.

Control bytes:
 c[0] bit  0    -u
      bits 1-2  -O
      bits 3-4  -l
      bits 5-6  -n (number of temporary files)
      bit  7    -f (fast mode)
 c[1] bits 0-1  -r (reads held in memory, with -f)
      bit  2    --no-PG
      bit  3    positional output prefix instead of -o -
      bits 4-7  unused
*/

static void args(fuzz_args *a, const uint8_t *c, char *in) {
    ARG("collate");
    ARG2("-T", fuzz_prefix);

    if (c[0] & 0x01)   ARG("-u");
    CHOOSE("-O", (c[0] >> 1) & 3, NULL, "sam", "bam", "cram,no_ref");
    CHOOSE("-l", (c[0] >> 3) & 3, NULL, "0", "1", "9");
    CHOOSE("-n", (c[0] >> 5) & 3, NULL, "1", "2", "64");
    if (c[0] & 0x80)   ARG("-f");

    CHOOSE("-r", c[1] & 3, NULL, "1", "10", "10000");
    if (c[1] & 0x04)   ARG("--no-PG");
    if (!(c[1] & 0x08))
        ARG2("-o", "-");

    ENDOPTS();
    ARG(in);
    if (c[1] & 0x08)
        ARG(fuzz_prefix);   // writes fuzz_prefix.bam
}

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
    return fuzz_one_input(data, size, 2, FUZZ_PARSE, args, main_bamshuf);
}
