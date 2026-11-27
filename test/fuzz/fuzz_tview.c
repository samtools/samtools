/*  test/fuzz/fuzz_tview.c -- Fuzz driver for samtools tview (text and HTML modes)

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

int bam_tview_main(int argc, char *argv[]);

/*
Fuzzes samtools tview in its non-interactive modes: -d T (text) and -d H
(HTML) render one screen to stdout and exit.  The curses mode (-d C, the
default on a terminal) is never used.  See fuzz_common.h for the input
layout (file byte, region bytes) and build/run details.

Link bam_tview.o, bam_tview_html.o, bam_tview_curses.o and bam_lpileup.o.
Output goes to stdout.  Runtime file: ref.fa (optional reference).

Tool byte:
 c[0] bit  0    -d H (else -d T)
      bit  1    -p <first region>
      bit  2    ref.fa as the reference
      bits 3-4  -w: none, 1, 40, 1000
      bit  5    -s rg1
      bits 6-7  unused
*/

static void args(fuzz_args *a, const uint8_t *c, fuzz_file *f) {
    ARG("tview");
    ARG2("-d", (c[0] & 0x01) ? "H" : "T");
    if (c[0] & 0x02)   ARG_OPTDATA("-p", f->reg[0].s);
    CHOOSE("-w", (c[0] >> 3) & 3, NULL, "1", "40", "1000");
    if (c[0] & 0x20)   ARG2("-s", "rg1");
    ENDOPTS();
    ARG(f->in);
    if (c[0] & 0x04)   ARG("ref.fa");
}

static const fuzz_step steps[] = {{args, bam_tview_main, NULL}};

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
    return fuzz_one_file_input(data, size, 1, FUZZ_SORTED, FUZZ_INDEX,
                               steps, 1);
}
