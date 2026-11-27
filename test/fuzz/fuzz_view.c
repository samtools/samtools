/*  test/fuzz/fuzz_view.c -- Fuzz driver for samtools view

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

int main_samview(int argc, char *argv[]);

/*
Fuzzes samtools view.  See fuzz_common.h for build and run details.

We construct valid argument combinations in addition to the input data.
Both are controlled by the fuzzer.

Runtime files, which should be viable given the initial test corpus:
ref.fa.fai (-t), qname.txt (-N), RG.txt (-R), TAG.txt (--tag-file),
skip.bed (-L).  Modifying the corpus will permit testing of invalid data
relative to the filtering files.

All outputs (-o, -U, --save-counts) go to /dev/null, so their code runs
without creating files.

Options needing an index (-P, region arguments) are not fuzzed here.

The control-byte layout matches the original fuzz_view.c, except that
-s now has its own bit (c[4] bit 0) instead of sharing c[0] bit 0 with -b.

Control bytes:
 c[0] bit  0    -b
      bit  1    -c
      bit  2    -C
      bit  3    -F 0xf00
      bit  4    -n
      bit  5    -t ref.fa.fai
      bit  6    unused (-P needs an index)
      bit  7    -H
 c[1] bit  0    --rf 0xc0
      bit  1    --remove-tag MD,NM
      bit  2    --keep-tag MD,NM
      bit  3    --library L
      bit  4    -q 10
      bit  5    unused (region needs an index)
      bit  6    -m 10
      bit  7    --no-header
 c[2] bit  0    --no-PG
      bit  1    -o /dev/null
      bit  2    -U /dev/null
      bit  3    -N qname.txt
      bit  4    -r rg1
      bit  5    -R RG.txt
      bit  6    --remove-flags 0xfff
      bit  7    --save-counts /dev/null
 c[3] bit  0    --tag NM
      bit  1    --tag NM:0
      bit  2    --tag-file NM:TAG.txt
      bit  3    -L skip.bed
      bit  4    --unmap
      bit  5    -M (needs an index, so exercises the error path)
      bit  6    -z unmap,pos,mqual,cigar,cigdup,cigarx,aux
      bit  7    --remove-ur
 c[4] bit  0    -s auto.1
      bits 1-7  unused

TODO: fuzz the expression language (-e), perhaps with a dedicated fuzzer.
*/

static void args(fuzz_args *a, const uint8_t *c, char *in) {
    ARG("view");

    if (c[0] & 0x01)   ARG("-b");
    if (c[0] & 0x02)   ARG("-c");
    if (c[0] & 0x04)   ARG("-C");
    if (c[0] & 0x08)   ARG2("-F", "0xf00");
    if (c[0] & 0x10)   ARG("-n");
    // -t takes a fai-style list of names and lengths, not a FASTA file
    if (c[0] & 0x20)   ARG2("-t", "ref.fa.fai");
    //if (c[0] & 0x40) ARG("-P"); // needs index
    if (c[0] & 0x80)   ARG("-H");

    if (c[1] & 0x01)   ARG2("--rf", "0xc0");
    if (c[1] & 0x02)   ARG2("--remove-tag", "MD,NM");
    if (c[1] & 0x04)   ARG2("--keep-tag", "MD,NM");
    if (c[1] & 0x08)   ARG2("--library", "L");
    if (c[1] & 0x10)   ARG2("-q", "10");
    //if (c[1] & 0x20) region "1:10-20" after the input; needs index
    if (c[1] & 0x40)   ARG2("-m", "10");
    if (c[1] & 0x80)   ARG("--no-header");

    if (c[2] & 0x01)   ARG("--no-PG");
    if (c[2] & 0x02)   ARG2("-o", "/dev/null");
    if (c[2] & 0x04)   ARG2("-U", "/dev/null");
    if (c[2] & 0x08)   ARG2("-N", "qname.txt");
    if (c[2] & 0x10)   ARG2("-r", "rg1");
    if (c[2] & 0x20)   ARG2("-R", "RG.txt");
    if (c[2] & 0x40)   ARG2("--remove-flags", "0xfff");
    if (c[2] & 0x80)   ARG2("--save-counts", "/dev/null");

    if (c[3] & 0x01)   ARG2("--tag", "NM");
    if (c[3] & 0x02)   ARG2("--tag", "NM:0");
    if (c[3] & 0x04)   ARG2("--tag-file", "NM:TAG.txt");
    if (c[3] & 0x08)   ARG2("-L", "skip.bed");
    if (c[3] & 0x10)   ARG("--unmap");
    if (c[3] & 0x20)   ARG("-M");
    if (c[3] & 0x40)   ARG2("-z", "unmap,pos,mqual,cigar,cigdup,cigarx,aux");
    if (c[3] & 0x80)   ARG("--remove-ur");

    if (c[4] & 0x01)   ARG2("-s", "auto.1");

    ENDOPTS();
    ARG(in);
}

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
    // The original rejected @SQ lengths over 1000; kept for continuity
    return fuzz_one_input(data, size, 5, FUZZ_SHORT_REFS, args, main_samview);
}
