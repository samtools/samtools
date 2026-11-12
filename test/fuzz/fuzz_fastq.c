/*  test/fuzz/fuzz_fastq.c -- Fuzz driver for samtools fastq / fasta

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

int main_bam2fq(int argc, char *argv[]);

/*
Fuzzes samtools fastq and samtools fasta.  See fuzz_common.h for build and
run details.

Link bam2fq.o.  No runtime files.

The default output goes to stdout.  The split outputs (-0 -1 -2 -s -o
--i1 --i2) all go to /dev/null, so that code runs without creating files.
As a result -c (compression level) only matters if a filename ends .gz,
which none do.

Control bytes:
 c[0] bit  0    fasta instead of fastq
      bit  1    -n
      bit  2    -N
      bit  3    -O
      bit  4    -t
      bits 5-6  -T
      bit  7    -i
 c[1] bits 0-1  -f
      bits 2-3  -F
      bits 4-5  --rf
      bits 6-7  -G
 c[2] bit  0    -1 /dev/null -2 /dev/null
      bit  1    -0 /dev/null
      bit  2    -s /dev/null
      bit  3    -o /dev/null
      bits 4-5  -v
      bits 6-7  --index-format
 c[3] bit  0    --i1 /dev/null --i2 /dev/null
      bits 1-2  --barcode-tag
      bits 3-4  --quality-tag
      bits 5-6  -c
      bit  7    unused
*/

static void args(fuzz_args *a, const uint8_t *c, char *in) {
    ARG((c[0] & 0x01) ? "fasta" : "fastq");

    if (c[0] & 0x02)   ARG("-n");
    if (c[0] & 0x04)   ARG("-N");
    if (c[0] & 0x08)   ARG("-O");
    if (c[0] & 0x10)   ARG("-t");
    CHOOSE("-T", (c[0] >> 5) & 3, NULL, "*", "RG,BC,MM,ML", "XA,XB,XH");
    if (c[0] & 0x80)   ARG("-i");

    CHOOSE("-f",   c[1] & 3,        NULL, "0x1", "READ1", "0xfff");
    CHOOSE("-F",   (c[1] >> 2) & 3, NULL, "0", "0x900", "0x4");
    CHOOSE("--rf", (c[1] >> 4) & 3, NULL, "0x40", "READ2", "0");
    CHOOSE("-G",   (c[1] >> 6) & 3, NULL, "0x80", "0x900", "0xfff");

    if (c[2] & 0x01) { ARG2("-1", "/dev/null"); ARG2("-2", "/dev/null"); }
    if (c[2] & 0x02)   ARG2("-0", "/dev/null");
    if (c[2] & 0x04)   ARG2("-s", "/dev/null");
    if (c[2] & 0x08)   ARG2("-o", "/dev/null");
    CHOOSE("-v", (c[2] >> 4) & 3, NULL, "0", "1", "93");
    CHOOSE("--index-format", (c[2] >> 6) & 3, NULL, "i*", "n2i*", "i6i6");

    if (c[3] & 0x01) { ARG2("--i1", "/dev/null"); ARG2("--i2", "/dev/null"); }
    CHOOSE("--barcode-tag", (c[3] >> 1) & 3, NULL, "BC", "RX", "XZ");
    CHOOSE("--quality-tag", (c[3] >> 3) & 3, NULL, "QT", "QX", "XZ");
    CHOOSE("-c",            (c[3] >> 5) & 3, NULL, "0", "1", "9");
    // Threads make runs non-deterministic.
    //ARG2("-@", "2");

    ENDOPTS();
    ARG(in);
}

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
    return fuzz_one_input(data, size, 4, FUZZ_PARSE, args, main_bam2fq);
}
