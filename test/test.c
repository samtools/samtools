/*  test/test.c -- test harness utility routines.

    Copyright (C) 2014, 2016, 2019 Genome Research Ltd.

    Author: Martin O. Pollard <mp15@sanger.ac.uk>

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

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <errno.h>
#include <inttypes.h>
#include <htslib/sam.h>

#include "test.h"

void xfreopen(const char *path, const char *mode, FILE *stream)
{
    if (freopen(path, mode, stream) == NULL) {
        fprintf(stderr, __FILE__": error reopening %s: %s\n",
                path, strerror(errno));
        exit(2);
    }
}

int redirect_stderr(const char *path) {
    int fd = open(path, O_WRONLY|O_TRUNC|O_CREAT, 0666);
    if (!fd) {
        fprintf(stderr, "Couldn't open \"%s\" : %s\n", path, strerror(errno));
        exit(2);
    }
    fflush(stderr);
    dup2(fd, STDERR_FILENO);
    return fd;
}

void flush_and_restore_stderr(int orig_stderr, int redirect_fd) {
    fflush(stderr);
    dup2(orig_stderr, STDERR_FILENO);
    close(redirect_fd);
}

void dump_hdr(const sam_hdr_t* hdr)
{
    printf("n_targets: %d\n", sam_hdr_nref(hdr));
    printf("ignore_sam_err: %d\n", hdr->ignore_sam_err);
    printf("l_text: %zu\n", (size_t) sam_hdr_length((sam_hdr_t*)hdr));
    printf("idx\ttarget_len\ttarget_name:\n");
    int32_t target;
    for (target = 0; target < sam_hdr_nref(hdr); ++target) {
        printf("%d\t%"PRId64"\t\"%s\"\n", target, (int64_t) sam_hdr_tid2len(hdr, target), sam_hdr_tid2name(hdr, target));
    }
    printf("text: \"%s\"\n", sam_hdr_str((sam_hdr_t*)hdr));
}

// For tests, just return a constant that can be embedded in expected output.
const char *samtools_version(void)
{
    return "x.y.test";
}
