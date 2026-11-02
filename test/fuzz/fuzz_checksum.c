/*  test/fuzz/fuzz_checksum.c -- Fuzz driver for samtools checksum

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

int main_checksum(int argc, char *argv[]);

/*
Fuzzes samtools checksum.  See fuzz_common.h for the input layout (file
byte, region bytes, extra payload) and build/run details.

The main input may be SAM/BAM/CRAM (re-encoded per the file byte) or
FASTQ (passed raw).  In -m (merge) mode the inputs are checksum report
files: the extra payload, then optionally the main input.

Link bam_checksum.o.  Output goes to stdout.  No runtime files.

Tool bytes:
 c[0] bit  0    -a
      bit  1    -P
      bit  2    -C
      bit  3    -M
      bit  4    -O
      bit  5    -c
      bit  6    -T
      bit  7    -B
 c[1] bits 0-1  -F
      bits 2-3  -f
      bits 4-5  -b
      bits 6-7  -t
 c[2] bits 0-1  -z
      bit  2    -q
      bit  3    -v
      bit  4    -m, merging checksum files
      bit  5    a second input (the extra payload)
      bits 6-7  unused
*/

static void args(fuzz_args *a, const uint8_t *c, fuzz_file *f) {
    ARG("checksum");
    if (c[0] & 0x01)   ARG("-a");
    if (c[0] & 0x02)   ARG("-P");
    if (c[0] & 0x04)   ARG("-C");
    if (c[0] & 0x08)   ARG("-M");
    if (c[0] & 0x10)   ARG("-O");
    if (c[0] & 0x20)   ARG("-c");
    if (c[0] & 0x40)   ARG("-T");
    if (c[0] & 0x80)   ARG("-B");

    CHOOSE("-F", c[1] & 3,        NULL, "0", "0x900", "UNMAP");
    CHOOSE("-f", (c[1] >> 2) & 3, NULL, "0", "QCFAIL", "0x1");
    CHOOSE("-b", (c[1] >> 4) & 3, NULL, "0xfff", "0", "PAIRED,READ1,READ2");
    CHOOSE("-t", (c[1] >> 6) & 3, NULL, "*", "*,MD,NM", "BC,QT,RG");

    CHOOSE("-z", c[2] & 3, NULL, "all", "off", "cigarx");
    if (c[2] & 0x04)   ARG("-q");
    if (c[2] & 0x08)   ARG("-v");
    int merge = c[2] & 0x10;
    if (merge)         ARG("-m");

    ENDOPTS();
    if (merge) {
        ARG(f->extra_raw_fn);
        if (c[2] & 0x20)   ARG(f->fn);
    } else {
        ARG(f->in);
        if (c[2] & 0x20)   ARG(f->extra_raw_fn);
    }
}

static const fuzz_step steps[] = {{args, main_checksum, NULL}};

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
    return fuzz_one_file_input(data, size, 3, FUZZ_ANYFMT, FUZZ_EXTRA,
                               steps, 1);
}
