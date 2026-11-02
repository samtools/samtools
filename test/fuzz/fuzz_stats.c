/*  test/fuzz/fuzz_stats.c -- Fuzz driver for samtools stats

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

int main_stats(int argc, char *argv[]);

/*
Fuzzes samtools stats.  See fuzz_common.h for build and run details.

Link stats.o, stats_isize.o and bam_sort.o's dependencies as needed.

Output goes to stdout.  Runtime files: ref.fa (+ .fai) for -r, and
targets.txt (chr, from, to; 1-based) for -t.

-S/--split is deliberately not used: it writes <prefix>_<value>.bamstat
where <value> comes from tag data in the input, making the output path
fuzzable.

Region arguments and -X need an index, so aren't fuzzed here.

Control bytes:
 c[0] bit  0    -d
      bits 1-2  -f
      bits 3-4  -F
      bit  5    -x
      bit  6    -p
      bit  7    -r ref.fa
 c[1] bits 0-1  -c
      bits 2-3  -i
      bits 4-5  -l
      bits 6-7  -q
 c[2] bits 0-1  -m
      bits 2-3  -g
      bits 4-5  --GC-depth
      bit  6    -t targets.txt
      bit  7    -I rg1
 c[3] bits 0-1  unused (formerly -S; see above)
      bit  2    --ref-stats
      bits 3-4  --ref-stats-chunk
      bits 5-7  unused
*/

static void args(fuzz_args *a, const uint8_t *c, char *in) {
    ARG("stats");

    if (c[0] & 0x01)   ARG("-d");
    CHOOSE("-f", (c[0] >> 1) & 3, NULL, "1", "PAIRED,READ1", "0xfff");
    CHOOSE("-F", (c[0] >> 3) & 3, NULL, "0", "0x900", "UNMAP,DUP");
    if (c[0] & 0x20)   ARG("-x");
    if (c[0] & 0x40)   ARG("-p");
    if (c[0] & 0x80)   ARG2("-r", "ref.fa");

    // Coverage min,max,step: normal, all-zero, and min > max
    CHOOSE("-c", c[1] & 3,        NULL, "1,10,1", "0,0,0", "5,1,1");
    CHOOSE("-i", (c[1] >> 2) & 3, NULL, "0", "10", "100000");
    CHOOSE("-l", (c[1] >> 4) & 3, NULL, "0", "30", "1");
    CHOOSE("-q", (c[1] >> 6) & 3, NULL, "0", "20", "255");

    CHOOSE("-m", c[2] & 3,        NULL, "0", "0.5", "1.5");
    CHOOSE("-g", (c[2] >> 2) & 3, NULL, "0", "1", "1000");
    CHOOSE("--GC-depth", (c[2] >> 4) & 3, NULL, "1", "100", "0");
    if (c[2] & 0x40)   ARG2("-t", "targets.txt");
    if (c[2] & 0x80)   ARG2("-I", "rg1");

    // -S/-P: never; output filenames would derive from the data.
    if (c[3] & 0x04)   ARG("--ref-stats");
    CHOOSE("--ref-stats-chunk", (c[3] >> 3) & 3, NULL, "1", "0", "1000");
    // Threads make runs non-deterministic.
    //ARG2("-@", "2");

    ENDOPTS();
    ARG(in);
}

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
    return fuzz_one_input(data, size, 4, FUZZ_PARSE, args, main_stats);
}
