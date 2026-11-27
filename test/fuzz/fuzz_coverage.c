/*  test/fuzz/fuzz_coverage.c -- Fuzz driver for samtools coverage

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

int main_coverage(int argc, char *argv[]);

/*
Fuzzes samtools coverage, with and without -r regions.  See fuzz_common.h
for the input layout (file byte, region bytes) and build/run details.

Link coverage.o.  Output goes to stdout or /dev/null.  No runtime files
(-b uses the driver's generated file of filenames).

Inputs must be coordinate sorted.

Tool bytes:
 c[0] bits 0-1  -l
      bits 2-3  -q
      bits 4-5  -Q
      bit  6    -m
      bit  7    -D
 c[1] bit  0    -A
      bit  1    -H
      bits 2-3  -w
      bits 4-5  -d
      bits 6-7  --min-depth
 c[2] bits 0-1  --rf
      bits 2-3  --ff
      bit  4    -r <first region>
      bit  5    a second copy of the input
      bit  6    -b <file of filenames> instead of positional input
      bit  7    -o /dev/null
*/

static void args(fuzz_args *a, const uint8_t *c, fuzz_file *f) {
    ARG("coverage");

    CHOOSE("-l", c[0] & 3,        NULL, "0", "20", "1000");
    CHOOSE("-q", (c[0] >> 2) & 3, NULL, "0", "20", "255");
    CHOOSE("-Q", (c[0] >> 4) & 3, NULL, "0", "20", "99");
    if (c[0] & 0x40)   ARG("-m");
    if (c[0] & 0x80)   ARG("-D");

    if (c[1] & 0x01)   ARG("-A");
    if (c[1] & 0x02)   ARG("-H");
    CHOOSE("-w",          (c[1] >> 2) & 3, NULL, "0", "1", "80");
    CHOOSE("-d",          (c[1] >> 4) & 3, NULL, "0", "1", "1000000");
    CHOOSE("--min-depth", (c[1] >> 6) & 3, NULL, "0", "2", "1000");

    CHOOSE("--rf", c[2] & 3,        NULL, "0x1", "READ2,REVERSE", "0");
    CHOOSE("--ff", (c[2] >> 2) & 3, NULL, "0", "0x900", "0xfff");
    if (c[2] & 0x10)   ARG_OPTDATA("-r", f->reg[0].s);
    if (c[2] & 0x80)   ARG2("-o", "/dev/null");
    if (c[2] & 0x40)   ARG2("-b", f->list_fn);

    ENDOPTS();
    if (!(c[2] & 0x40)) {
        ARG(f->in);
        if (c[2] & 0x20)   ARG(f->in);
    }
}

static const fuzz_step steps[] = {{args, main_coverage}};

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
    return fuzz_one_file_input(data, size, 3, FUZZ_SORTED | FUZZ_SHORT_REFS,
                               1, steps, 1);
}
