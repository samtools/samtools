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

#include <config.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <setjmp.h>
#include <getopt.h>
#include <htslib/sam.h>
#include <htslib/hfile.h>
#include <htslib/hts.h>

int main_samview(int argc, char *argv[]);

/*
Fuzzes samtools view

Build with ./build_view.sh

Run with: ./run_view.sh
export ASAN_SYMBOLIZER_PATH=/software/badger/opt/llvm/21.1.8/bin/llvm-symbolizer 
ASAN_OPTIONS=allow_addr2line=1:abort_on_error=1:detect_leaks=0 ./a.out  -rss_limit_mb=2560 -timeout=10 

We construct valid argument combinations in addition to the input data.
Both are controlled by the fuzzer.

Arguments may attempt to load files exclude,file, skip.bed, and RG.txt.
These should be constant that are viable given the initial test corpus.
Modifying the corpus will permit testing of invalid data relative to the
filtering files.

No validity of result is assumed.  We're simply checking for crashes or
illegal program behaviour.
*/

// Taken from bamtk.c which we exclude from our build line
const char *samtools_version(void) {
    return "fuzzer";
}

// Returns "data:,<base64-data>" as a malloced buffer
char *base64_encode(const uint8_t *data, size_t len) {
    static const char table[] =
        "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

    size_t out_len = 4 * ((len + 2) / 3);

    char *output = malloc(out_len + 20);
    if (!output)
        return NULL;

    size_t i = 0, j = 0;
    // data: URI
    strcpy(output, "data:;base64,");
    j = strlen(output);

    // base-64
    while (i < len) {
        unsigned int a =             data[i++];
        unsigned int b = (i < len) ? data[i++] : 0;
        unsigned int c = (i < len) ? data[i++] : 0;

        unsigned int triple = (a << 16) | (b << 8) | c;

        output[j++] = table[(triple >> 18) & 0x3F];
        output[j++] = table[(triple >> 12) & 0x3F];
        output[j++] = (i - 1 < len) ? table[(triple >> 6) & 0x3F] : '=';
        output[j++] = (i < len)     ? table[triple & 0x3F]        : '=';
    }

    output[j] = '\0';
    return output;
}

void dump(int argc, char **argv) {
    fprintf(stderr, "Calling: ");
    for (int i = 0; i < argc; i++)
	fprintf(stderr, " %s", argv[i]);
    fprintf(stderr, "\n");
}

// Mpileup calls exit() directly.  This is problamtic for the fuzzer, but
// we can compile it with -Dexit=fuzz_exit and write our own fuzz_exit
// function that unwinds back to LLVMFuzzerTestOneInput via setjmp
static jmp_buf fuzz_jmp;
void fuzz_exit(int status) {
    fprintf(stderr, "Caught exit of status %d\n", status);
    longjmp(fuzz_jmp, status + 1);
}

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
    int jret = setjmp(fuzz_jmp);
    if (jret != 0) {
	// if non-zero don't add it to the test corpus
	return jret ? -1 : 0;
    }

    if (size < 10)
	return 0;

    // TODO: validate the input data is parseable first, so we don't spend
    // a lot of time on argument processing only to instantly reject it.
    // We can do this via a mem: and a read loop.
    hFILE *memfile;
    uint8_t *copy = malloc(size);
    if (copy == NULL) {
        abort();
    }

    memcpy(copy, data+5, size-5);
    if (!(memfile = hopen("mem:", "rb:", copy, size-5))) {
	free(copy);
	return -1;
    }
    htsFile *in = hts_hopen(memfile, "data", "rb");
    if (in == NULL) {
        if (hclose(memfile) != 0) {
            abort();
        }
        return -1;
    }

    sam_hdr_t *hdr = sam_hdr_read(in);
    if (!hdr) {
	hts_close(in);
	return -1;
    }

    // Reject long SQ lines, so we have faster mpileup -a
    int n_targets = sam_hdr_nref(hdr);
    for (int i = 0; i < n_targets; i++) {
        hts_pos_t len = sam_hdr_tid2len(hdr, i);
	if (len > 1000) {
	    hts_close(in);
	    sam_hdr_destroy(hdr);
	    return -1;
	}
    }

    // Reject inputs that don't parse, as we don't need to fuzz that
    // (htslib does already), and our tests are slow so reserve it for
    // correct data.
    //
    // We could do the same with unsorted data, but don't yet.
    // TODO: maybe we need to truncate the input somehow?  So we get
    // the valid data only?
    bam1_t *b = bam_init1();
    if (!b)
	abort();
    int r;
    while ((r = sam_read1(in, hdr, b)) >= 0)
	;
    bam_destroy1(b);
    sam_hdr_destroy(hdr);
    hts_close(in);
    if (r < -1)
	return -1; // error


    // Use the first 5 bytes as argument control so we can explore
    // the parameter space as well as the data space.
    char *argv[100] = {0}, argc = 0;

    argv[argc++] = "samtools";
    argv[argc++] = "view";

    if (data[0] & 0x01) { argv[argc++] = "-s"; argv[argc++] = "auto.1"; }
    // TODO: fuzz expression language.  Use a dedicated fuzzer maybe

    // Some options need an index.  We would need to first save the data and
    // index it before fuzzing these, so they're commented out for now.

    if (data[0] & 0x01)   argv[argc++] = "-b";
    if (data[0] & 0x02)   argv[argc++] = "-c";
    if (data[0] & 0x04)   argv[argc++] = "-C";
    if (data[0] & 0x08) { argv[argc++] = "-F"; argv[argc++] = "0xf00"; }
    if (data[0] & 0x10)   argv[argc++] = "-n";
    if (data[0] & 0x20) { argv[argc++] = "-t"; argv[argc++] = "ref.fa"; }
    //if (data[0] & 0x40)   argv[argc++] = "-P"; // needs index
    if (data[0] & 0x80)   argv[argc++] = "-H";

    if (data[1] & 0x01) { argv[argc++] = "--rf"; argv[argc++] = "0xc0"; }
    if (data[1] & 0x02) { argv[argc++] = "--remove-tag"; argv[argc++] = "MD,NM"; }
    if (data[1] & 0x04) { argv[argc++] = "--keep-tag"; argv[argc++] = "MD,NM"; }
    if (data[1] & 0x08) { argv[argc++] = "--library"; argv[argc++] = "L"; }
    if (data[1] & 0x10) { argv[argc++] = "-q"; argv[argc++] = "10"; }
    //if (data[1] & 0x20) { argv[argc++] = "-r"; argv[argc++] = "1:10-20"; }
    if (data[1] & 0x40) { argv[argc++] = "-m"; argv[argc++] = "10"; }
    if (data[1] & 0x80)   argv[argc++] = "--no-header";
	       
    if (data[2] & 0x01)   argv[argc++] = "--no-PG";
    if (data[2] & 0x02) { argv[argc++] = "-o"; argv[argc++] = "_out.sam"; }
    if (data[2] & 0x04) { argv[argc++] = "-U"; argv[argc++] = "_unout.sam"; }
    if (data[2] & 0x08) { argv[argc++] = "-N"; argv[argc++] = "qname.txt"; }
    if (data[2] & 0x10) { argv[argc++] = "-r"; argv[argc++] = "rg1"; }
    if (data[2] & 0x20) { argv[argc++] = "-R"; argv[argc++] = "RG.txt"; }
    if (data[2] & 0x40) { argv[argc++] = "--remove-flags"; argv[argc++] = "0xfff"; }
    if (data[2] & 0x80) { argv[argc++] = "--save-counts"; argv[argc++] = "_counts"; }
        
    if (data[3] & 0x01) { argv[argc++] = "--tag"; argv[argc++]="NM"; }
    if (data[3] & 0x02) { argv[argc++] = "--tag"; argv[argc++]="NM:0"; }
    if (data[3] & 0x04) { argv[argc++] = "--tag-file"; argv[argc++]="NM:TAG.txt"; }
    if (data[3] & 0x08) { argv[argc++] = "-L"; argv[argc++] = "skip.bed"; }
    if (data[3] & 0x10)   argv[argc++] = "--unmap";
    if (data[3] & 0x20)   argv[argc++] = "-M"; // even though it won't work
    if (data[3] & 0x40) { argv[argc++] = "-z"; argv[argc++] = "unmap,pos,mqual,cigar,cigdup,cigarx,aux"; }
    if (data[3] & 0x80)   argv[argc++] = "--remove-ur";

    // data[4] unused for now

    // Create a data file to read in.
    char *b64 = base64_encode(data+5, size-5);
    argv[argc++] = b64;

    dump(argc, argv);

    optind = 1;
    opterr = 0;

    int ret = main_samview(argc-1, argv+1);
    free(b64);

    return 0;
}
