/*  test/fuzz/fuzz_consensus.c -- Fuzz driver for samtools consensus

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

int main_consensus(int argc, char *argv[]);

/*
Fuzzes samtools consensus

See fuzz_common.h for build and run details, including FUZZ_EXTRACT for
reproducing a crash outside the harness.

The leading control bytes of each input select the command line options.
The remainder is the alignment data, passed to samtools as a data: URL.
Both are controlled by the fuzzer.

-o is never used, so consensus writes to stdout.  Use -close_fd_mask=1
to discard it.

Options taking filenames use fixed names which must exist in the working
directory:
    ref.fa, ref.fa.fai   for -T  (create the .fai up front, otherwise
                                  htslib will attempt to write one)
    qcal.txt             for -t  (QUAL value sub undercall overcall lines)
These should be viable given the initial test corpus, e.g. ref.fa holding
the @SQ names used.  Modifying the corpus then tests mismatches between
the data and these files.

Options which need an index (-r, --regions-file) are commented out, as the
data: URL has none.  Threading (-@, -Z) is also disabled as it makes runs
non-reproducible.

No validity of result is assumed.  We're simply checking for crashes or
illegal program behaviour.
*/

/*
 * Control bytes:
 *
 * c[0] bits 0-1  -f
 *      bits 2-3  -m
 *      bit  4    -a
 *      bit  5    -a (again; -aa)
 *      bit  6    -A
 *      bit  7    --mark-ins
 * c[1] bit  0    --show-del yes
 *      bit  1    --show-ins no
 *      bits 2-3  -l
 *      bit  4    --rf 0x2
 *      bit  5    --ff 0
 *      bit  6    --min-MQ 10
 *      bit  7    --min-BQ 10
 * c[2] bits 0-1  -d
 *      bit  2    -T ref.fa
 *      bit  3    --ref-qual 20
 *      bit  4    -q
 *      bits 5-6  -c
 *      bit  7    -H 0.1
 * c[3] bits 0-1  -C
 *      bit  2    --no-use-MQ
 *      bit  3    --no-adj-MQ
 *      bit  4    --no-adj-qual
 *      bit  5    --NM-halo 1
 *      bit  6    --low-MQ 20
 *      bit  7    --high-MQ 10 (with bit 6, low > high)
 * c[4] bit  0    --scale-MQ 1.5
 *      bit  1    --P-het 0.1
 *      bit  2    --P-indel 0.01
 *      bit  3    --het-scale 0.37
 *      bit  4    -p
 *      bit  5    --homopoly-score 0.3
 *      bits 6-7  -t
 * c[5] bits 0-2  -X config
 *      bit  6    -X first rather than last (later options override it)
 *      bits 3-5, 7 unused
 */
static void args(fuzz_args *a, const uint8_t *c, char *in) {
    ARG("consensus");

    // -X sets several parameters at once, and later options override
    // earlier ones, so c[5] bit 6 decides whether it goes first or last.
    static const char *config[8] = {
        NULL, "hiseq", "hifi", "r10.4_sup", "r10.4_dup", "ultima", NULL, NULL
    };
    const char *X = config[c[5] & 7];
    if (X && (c[5] & 0x40)) ARG2("-X", X);

    // c[0]: general options
    CHOOSE("-f", (c[0]) & 3,    NULL, "fasta",  "fastq",        "pileup");
    CHOOSE("-m", (c[0] >> 2) & 3, NULL, "simple", "bayesian_116", "Bayesian");
    if (c[0] & 0x10)   ARG("-a");
    if (c[0] & 0x20)   ARG("-a"); // -a -a: refs with no data too
    if (c[0] & 0x40)   ARG("-A");
    if (c[0] & 0x80)   ARG("--mark-ins");

    // c[1]: general options, mostly filtering
    if (c[1] & 0x01) ARG2("--show-del", "yes");
    if (c[1] & 0x02) ARG2("--show-ins", "no");
    CHOOSE("-l", (c[1] >> 2) & 3, NULL, "0", "1", "15");
    if (c[1] & 0x10) ARG2("--rf", "0x2");
    if (c[1] & 0x20) ARG2("--ff", "0");
    if (c[1] & 0x40) ARG2("--min-MQ", "10");
    if (c[1] & 0x80) ARG2("--min-BQ", "10");

    // c[2]: general and simple-mode options
    CHOOSE("-d", (c[2]) & 3,    NULL, "0", "2", "10");
    if (c[2] & 0x04) ARG2("-T", "ref.fa");
    if (c[2] & 0x08) ARG2("--ref-qual", "20");
    if (c[2] & 0x10)   ARG("-q");
    CHOOSE("-c", (c[2] >> 5) & 3, NULL, "0", "0.5", "1");
    if (c[2] & 0x80) ARG2("-H", "0.1");

    // c[3]: Bayesian-mode options
    CHOOSE("-C", (c[3]) & 3,    NULL, "0", "30", "99");
    if (c[3] & 0x04)   ARG("--no-use-MQ");
    if (c[3] & 0x08)   ARG("--no-adj-MQ");
    if (c[3] & 0x10)   ARG("--no-adj-qual");
    if (c[3] & 0x20) ARG2("--NM-halo", "1");
    if (c[3] & 0x40) { ARG("--low-MQ");  ARG("20"); }
    if (c[3] & 0x80) ARG2("--high-MQ", "10");
    // (both together gives low > high)

    // c[4]: more Bayesian-mode options
    if (c[4] & 0x01) ARG2("--scale-MQ", "1.5");
    if (c[4] & 0x02) { ARG("--P-het");    ARG("0.1"); }
    if (c[4] & 0x04) { ARG("--P-indel");  ARG("0.01"); }
    if (c[4] & 0x08) ARG2("--het-scale", "0.37");
    if (c[4] & 0x10)   ARG("-p");
    if (c[4] & 0x20) ARG2("--homopoly-score", "0.3");
    CHOOSE("-t", (c[4] >> 6) & 3, NULL, ":hiseq", ":r10.4_sup", "qcal.txt");

    // c[5]: bits 0-2 and 6 are -X, handled above and below.
    // These need an index.  We would need to first save the data and
    // index it before fuzzing these, so they're commented out for now.
    //if (c[5] & 0x08) ARG2("-r", "1:10-20");
    //if (c[5] & 0x10) ARG2("--regions-file", "regions.bed");
    // Threads make runs non-deterministic.
    //if (c[5] & 0x20) { ARG("-@"); ARG("2");
    //                   ARG("-Z"); ARG("100"); }
    // c[5] bit 7 unused for now

    if (X && !(c[5] & 0x40)) ARG2("-X", X);

    // Create a data file to read in.

    ENDOPTS();
    ARG(in);

}

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
    return fuzz_one_input(data, size, 6, FUZZ_SORTED | FUZZ_SHORT_REFS, args, main_consensus);
}
