/*  test/fuzz/fuzz_reference.c -- Fuzz driver for samtools reference

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
int main_reference(int argc, char *argv[]);

/*
Fuzzes samtools reference.  See fuzz_common.h for the input layout (file
byte, region bytes) and build/run details.

Without -e the reference is built from MD tags in the input.  -e needs a
CRAM with an embedded reference, so a first step converts the input with
"samtools view -O cram,embed_ref=2" (htslib derives the embedded
reference from the data, as no reference is available) to a fixed name.

Link bam_reference.o (or reference.c) and sam_view.o.  Output goes to
stdout or /dev/null.  No runtime files.

Tool byte:
 c[0] bit  0    -e, on the embedded-reference CRAM
      bit  1    -q
      bit  2    -r <first region> (needs the index; only the main input
                has one)
      bit  3    -o /dev/null
      bits 4-7  unused
*/

static const char *emb_cram(void) {
    static char fn[340];
    snprintf(fn, sizeof(fn), "%s.emb.cram", fuzz_prefix);
    return fn;
}

static void embed_args(fuzz_args *a, const uint8_t *c, fuzz_file *f) {
    (void) c;
    ARG("view");
    ARG2("-O", "cram,embed_ref=2");
    ARG2("-o", emb_cram());
    ENDOPTS();
    ARG(f->fn);
}

static void ref_args(fuzz_args *a, const uint8_t *c, fuzz_file *f) {
    ARG("reference");
    if (c[0] & 0x01)   ARG("-e");
    if (c[0] & 0x02)   ARG("-q");
    if (c[0] & 0x04)   ARG_OPTDATA("-r", f->reg[0].s);
    if (c[0] & 0x08)   ARG2("-o", "/dev/null");
    ENDOPTS();
    ARG((c[0] & 0x01) ? emb_cram() : f->in);
}

static const fuzz_step steps[] = {
    {embed_args, main_samview,   NULL},
    {ref_args,   main_reference, NULL},
};

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
    return fuzz_one_file_input(data, size, 1, FUZZ_PARSE, FUZZ_INDEX,
                               steps, 2);
}
