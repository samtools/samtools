/*  test/fuzz/fuzz_reheader.c -- Fuzz driver for samtools reheader

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

int main_reheader(int argc, char *argv[]);
int main_samview(int argc, char *argv[]);

/*
Fuzzes samtools reheader, then reads the result back with samtools view.
See fuzz_common.h for the input layout (file byte, region bytes, extra
payload) and build/run details.

The extra payload is the new header (SAM header text).  With -i the
alignment file (a fixed name under fuzz_prefix) is rewritten in place;
this is only supported for CRAM, so other forms test the error path.  The
second step then decodes the patched file.  Without -i the output goes to
stdout and the second step reads the unmodified file.

-c CMD is never used: it runs a shell command.

Link bam_reheader.o and sam_view.o.  No runtime files.

Tool byte:
 c[0] bit  0    -P
      bit  1    -i
      bits 2-7  unused
*/

static void reheader_args(fuzz_args *a, const uint8_t *c, fuzz_file *f) {
    ARG("reheader");
    if (c[0] & 0x01)   ARG("-P");
    if (c[0] & 0x02)   ARG("-i");
    ENDOPTS();
    ARG(f->extra_raw_fn);
    ARG(f->fn);
}

static void view_args(fuzz_args *a, const uint8_t *c, fuzz_file *f) {
    (void) c;
    ARG("view");
    ENDOPTS();
    ARG(f->fn);
}

static const fuzz_step steps[] = {
    {reheader_args, main_reheader},
    {view_args,     main_samview},
};

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
    return fuzz_one_file_input(data, size, 1, FUZZ_PARSE, FUZZ_EXTRA,
                               steps, 2);
}
