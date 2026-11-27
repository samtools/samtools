/*  test/fuzz/fuzz_quickcheck.c -- Fuzz driver for samtools quickcheck

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

int main_quickcheck(int argc, char *argv[]);

/*
Fuzzes samtools quickcheck, including on damaged files.  See fuzz_common.h
for the input layout (file byte, region bytes) and build/run details.

The input is written in the chosen form (or raw, if it isn't alignment
data), then a harness-run "truncate" step optionally cuts bytes off the
end.  Cutting 28 bytes removes exactly the BGZF EOF block; other amounts
leave partial blocks or containers.  repro.sh records the same truncate,
so it is not idempotent: re-extract before running it a second time.

Link bam_quickcheck.o.  Output goes to stdout.  No runtime files.

Tool byte:
 c[0] bit  0    -v
      bit  1    -v (again; -vv)
      bit  2    -q
      bit  3    -u
      bit  4    a second copy of the input
      bits 5-7  bytes to cut from the end: 0, 1, 2, 27, 28, 29, 100, all
*/

static void truncate_args(fuzz_args *a, const uint8_t *c, fuzz_file *f) {
    static const char *cut[8] = {
        "+0", "-1", "-2", "-27", "-28", "-29", "-100", "0"
    };
    ARG2("-s", cut[(c[0] >> 5) & 7]);
    ENDOPTS();
    ARG(f->fn);
}

static void args(fuzz_args *a, const uint8_t *c, fuzz_file *f) {
    ARG("quickcheck");
    if (c[0] & 0x01)   ARG("-v");
    if (c[0] & 0x02)   ARG("-v");
    if (c[0] & 0x04)   ARG("-q");
    if (c[0] & 0x08)   ARG("-u");
    ENDOPTS();
    ARG(f->fn);
    if (c[0] & 0x10)   ARG(f->fn);
}

static const fuzz_step steps[] = {
    {truncate_args, NULL,            "truncate"},
    {args,          main_quickcheck, NULL},
};

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
    return fuzz_one_file_input(data, size, 1, FUZZ_ANYFMT, 0, steps, 2);
}
