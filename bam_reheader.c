/*  bam_reheader.c -- reheader subcommand.

    Copyright (C) 2010 Broad Institute.
    Copyright (C) 2012-2015 Genome Research Ltd.

    Author: Heng Li <lh3@sanger.ac.uk>

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

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <getopt.h>

#include "htslib/bgzf.h"
#include "htslib/sam.h"
#include "htslib/hfile.h"
#include "htslib/cram.h"
#include "samtools.h"

#define BUF_SIZE 0x10000

/*
 * Reads a file and outputs a new BAM file to fd with 'h' replaced as
 * the header.    No checks are made to the validity.
 */
int bam_reheader(BGZF *in, bam_hdr_t *h, int fd,
                 const char *arg_list, int add_PG)
{
    BGZF *fp = NULL;
    ssize_t len;
    uint8_t *buf = NULL;
    SAM_hdr *sh = NULL;
    if (in->is_write) return -1;
    buf = malloc(BUF_SIZE);
    if (!buf) {
        fprintf(stderr, "Out of memory\n");
        return -1;
    }
    if (bam_hdr_read(in) == NULL) {
        fprintf(stderr, "Couldn't read header\n");
        goto fail;
    }
    fp = bgzf_fdopen(fd, "w");
    if (!fp) {
        print_error_errno("reheader", "Couldn't open output file");
        goto fail;
    }

    if (add_PG) {
        // Around the houses, but it'll do until we can manipulate bam_hdr_t natively.
        sh = sam_hdr_parse_(h->text, h->l_text);
        if (!sh)
            goto fail;
        if (sam_hdr_add_PG(sh, "samtools",
                           "VN", samtools_version(),
                           arg_list ? "CL": NULL,
                           arg_list ? arg_list : NULL,
                           NULL) != 0)
            goto fail;

        free(h->text);
        h->text = strdup(sam_hdr_str(sh));
        h->l_text = sam_hdr_length(sh);
        if (!h->text)
            goto fail;
        sam_hdr_free(sh);
        sh = NULL;
    }

    if (bam_hdr_write(fp, h) < 0) {
        print_error_errno("reheader", "Couldn't write header");
        goto fail;
    }
    if (in->block_offset < in->block_length) {
        if (bgzf_write(fp, in->uncompressed_block + in->block_offset, in->block_length - in->block_offset) < 0) goto write_fail;
        if (bgzf_flush(fp) < 0) goto write_fail;
    }
    while ((len = bgzf_raw_read(in, buf, BUF_SIZE)) > 0) {
        if (bgzf_raw_write(fp, buf, len) < 0) goto write_fail;
    }
    if (len < 0) {
        fprintf(stderr, "[%s] Error reading input file\n", __func__);
        goto fail;
    }
    free(buf);
    fp->block_offset = in->block_offset = 0;
    if (bgzf_close(fp) < 0) {
        fprintf(stderr, "[%s] Error closing output file\n", __func__);
        return -1;
    }
    return 0;

 write_fail:
    print_error_errno("reheader", "Error writing to output file");
 fail:
    bgzf_close(fp);
    free(buf);
    sam_hdr_free(sh);
    return -1;
}

/*
 * Reads a file and outputs a new CRAM file to stdout with 'h'
 * replaced as the header.  No checks are made to the validity.
 *
 * FIXME: error checking
 */
int cram_reheader(cram_fd *in, bam_hdr_t *h, const char *arg_list, int add_PG)
{
    htsFile *h_out = hts_open("-", "wc");
    cram_fd *out = h_out->fp.cram;
    cram_container *c = NULL;
    int ret = -1;

    // Attempt to fill out a cram->refs[] array from @SQ headers
    cram_fd_set_header(out, sam_hdr_parse_(h->text, h->l_text));
    if (add_PG) {
        if (sam_hdr_add_PG(cram_fd_get_header(out), "samtools",
                           "VN", samtools_version(),
                           arg_list ? "CL": NULL,
                           arg_list ? arg_list : NULL,
                           NULL) != 0)
            goto err;

        // Covert back to bam_hdr_t struct
        free(h->text);
        h->text = strdup(sam_hdr_str(cram_fd_get_header(out)));
        h->l_text = sam_hdr_length(cram_fd_get_header(out));
        if (!h->text)
            goto err;
    }

    if (sam_hdr_write(h_out, h) != 0)
        goto err;
    cram_set_option(out, CRAM_OPT_REFERENCE, NULL);

    while ((c = cram_read_container(in))) {
        int32_t i, num_blocks = cram_container_get_num_blocks(c);
        if (cram_write_container(out, c) != 0)
            goto err;

        for (i = 0; i < num_blocks; i++) {
            cram_block *blk = cram_read_block(in);
            if (!blk || cram_write_block(out, blk) != 0) {
                if (blk) cram_free_block(blk);
                goto err;
            }
            cram_free_block(blk);
        }
        cram_free_container(c);
    }

    ret = 0;

 err:
    if (hts_close(h_out) != 0)
        ret = -1;

    return ret;
}



/*
 * Reads a version 2 CRAM file and replaces the header in-place,
 * provided the header is small enough to fit without growing the
 * entire file.
 *
 * Version 2 format has an uncompressed SAM header with multiple nul
 * termination bytes to permit inline header editing.
 *
 * Returns 0 on success;
 *        -1 on general failure;
 *        -2 on failure due to insufficient size
 */
int cram_reheader_inplace2(cram_fd *fd, const bam_hdr_t *h, const char *arg_list,
                          int add_PG)
{
    cram_container *c = NULL;
    cram_block *b = NULL;
    SAM_hdr *hdr = NULL;
    off_t start;
    int ret = -1;

    if (cram_major_vers(fd) < 2 ||
        cram_major_vers(fd) > 3) {
        fprintf(stderr, "[%s] unsupported CRAM version %d\n", __func__,
                cram_major_vers(fd));
        goto err;
    }

    if (!(hdr = sam_hdr_parse_(h->text, h->l_text)))
        goto err;

    if (add_PG && sam_hdr_add_PG(hdr, "samtools", "VN", samtools_version(),
                                 arg_list ? "CL": NULL,
                                 arg_list ? arg_list : NULL,
                                 NULL))
        goto err;

    int header_len = sam_hdr_length(hdr);
    /* Fix M5 strings? Maybe out of scope for this tool */

    // Load the existing header
    if ((start = hseek(cram_fd_get_fp(fd), 26, SEEK_SET)) != 26)
        goto err;

    if (!(c = cram_read_container(fd)))
        goto err;

    // Version 2.1 has a single uncompressed block which is nul
    // terminated with many nuls to permit growth.
    //
    // So load old block and keep all contents identical bar the
    // header text itself
    if (!(b = cram_read_block(fd)))
        goto err;

    if (cram_block_get_uncomp_size(b) < header_len+4) {
        fprintf(stderr, "New header will not fit. Use non-inplace version (%d > %d)\n",
                header_len+4, cram_block_get_uncomp_size(b));
        ret = -2;
        goto err;
    }

    cram_block_set_offset(b, 0);   // rewind block
    int32_put_blk(b, header_len);
    cram_block_append(b, sam_hdr_str(hdr), header_len);
    // Zero the remaining block
    memset(cram_block_get_data(b)+cram_block_get_offset(b), 0,
           cram_block_get_uncomp_size(b) - cram_block_get_offset(b));
    // Make sure all sizes and byte-offsets are consistent after memset
    cram_block_set_offset(b, cram_block_get_uncomp_size(b));
    cram_block_set_comp_size(b, cram_block_get_uncomp_size(b));

    if (hseek(cram_fd_get_fp(fd), start, SEEK_SET) != start)
        goto err;

    if (cram_write_container(fd, c) == -1)
        goto err;

    if (cram_write_block(fd, b) == -1)
        goto err;

    ret = 0;
 err:
    if (c) cram_free_container(c);
    if (b) cram_free_block(b);
    if (hdr) sam_hdr_free(hdr);

    return ret;
}


/*
 * Reads a version 3 CRAM file and replaces the header in-place,
 * provided the header is small enough to fit without growing the
 * entire file.
 *
 * Version 3 format has a SAM header held as an (optionally)
 * compressed block within the header container.  Additional
 * uncompressed blocks or simply unallocated space (the difference
 * between total block sizes and the container size) are used to
 * provide room for growth or contraction of the compressed header.
 *
 * Returns 0 on success;
 *        -1 on general failure;
 *        -2 on failure due to insufficient size
 */
int cram_reheader_inplace3(cram_fd *fd, const bam_hdr_t *h, const char *arg_list,
                          int add_PG)
{
    cram_container *c = NULL;
    cram_block *b = NULL;
    SAM_hdr *hdr = NULL;
    off_t start, sz, end;
    int container_sz, max_container_sz;
    char *buf = NULL;
    int ret = -1;

    if (cram_major_vers(fd) < 2 ||
        cram_major_vers(fd) > 3) {
        fprintf(stderr, "[%s] unsupported CRAM version %d\n", __func__,
                cram_major_vers(fd));
        goto err;
    }

    if (!(hdr = sam_hdr_parse_(h->text, h->l_text)))
        goto err;

    if (add_PG && sam_hdr_add_PG(hdr, "samtools", "VN", samtools_version(),
                                 arg_list ? "CL": NULL,
                                 arg_list ? arg_list : NULL,
                                 NULL))
        goto err;

    int header_len = sam_hdr_length(hdr);
    /* Fix M5 strings? Maybe out of scope for this tool */

    // Find current size of SAM header block
    if ((start = hseek(cram_fd_get_fp(fd), 26, SEEK_SET)) != 26)
        goto err;

    if (!(c = cram_read_container(fd)))
        goto err;

    // +5 allows num_landmarks to increase from 0 to 1 (Cramtools)
    max_container_sz = cram_container_size(c)+5;

    sz = htell(cram_fd_get_fp(fd)) + cram_container_get_length(c) - start;
    end = htell(cram_fd_get_fp(fd)) + cram_container_get_length(c);

    // We force 1 block instead of (optionally) 2.  C CRAM
    // implementations for v3 were writing 1 compressed block followed
    // by 1 uncompressed block.  However this is tricky to deal with
    // as changing block sizes can mean the block header also changes
    // size due to itf8 and variable size integers.
    //
    // If we had 1 block, this doesn't change anything.
    // If we had 2 blocks, the new container header will be smaller by
    // 1+ bytes, requiring the cram_container_get_length(c) to be larger in value.
    // However this is an int32 instead of itf8 so the container
    // header structure stays the same size.  This means we can always
    // reduce the number of blocks without running into size problems.
    cram_container_set_num_blocks(c, 1);
    int32_t *landmark;
    int32_t num_landmarks;
    landmark = cram_container_get_landmarks(c, &num_landmarks);
    if (num_landmarks && landmark) {
        num_landmarks = 1;
        landmark[0] = 0;
    } else {
        num_landmarks = 0;
    }
    cram_container_set_landmarks(c, num_landmarks, landmark);

    buf = malloc(max_container_sz);
    container_sz = max_container_sz;
    if (cram_store_container(fd, c, buf, &container_sz) != 0)
        goto err;

    if (!buf)
        goto err;

    // Proposed new length, but changing cram_container_get_length(c) may change the
    // container_sz and thus the remainder (cram_container_get_length(c) itself).
    cram_container_set_length(c, sz - container_sz);

    int old_container_sz = container_sz;
    container_sz = max_container_sz;
    if (cram_store_container(fd, c, buf, &container_sz) != 0)
        goto err;

    if (old_container_sz != container_sz) {
        fprintf(stderr, "Quirk of fate makes this troublesome! "
                "Please use non-inplace version.\n");
        goto err;
    }



    // Version 3.0 supports compressed header
    b = cram_new_block(FILE_HEADER, 0);
    int32_put_blk(b, header_len);
    cram_block_append(b, sam_hdr_str(hdr), header_len);
    cram_block_update_size(b);

    cram_compress_block(fd, b, NULL, -1, -1);

    if (hseek(cram_fd_get_fp(fd), 26, SEEK_SET) != 26)
        goto err;

    if (cram_block_size(b) > cram_container_get_length(c)) {
        fprintf(stderr, "New header will not fit. Use non-inplace version"
                " (%d > %d)\n",
                (int)cram_block_size(b), cram_container_get_length(c));
        ret = -2;
        goto err;
    }

    if (cram_write_container(fd, c) == -1)
        goto err;

    if (cram_write_block(fd, b) == -1)
        goto err;

    // Blank out the remainder
    int rsz = end - htell(cram_fd_get_fp(fd));
    assert(rsz >= 0);
    if (rsz) {
        char *rem = calloc(1, rsz);
        ret = hwrite(cram_fd_get_fp(fd), rem, rsz) == rsz ? 0 : -1;
        free(rem);
    }

 err:
    if (c) cram_free_container(c);
    if (buf) free(buf);
    if (b) cram_free_block(b);
    if (hdr) sam_hdr_free(hdr);

    return ret;
}

int cram_reheader_inplace(cram_fd *fd, const bam_hdr_t *h, const char *arg_list,
                         int add_PG)
{
    switch (cram_major_vers(fd)) {
    case 2: return cram_reheader_inplace2(fd, h, arg_list, add_PG);
    case 3: return cram_reheader_inplace3(fd, h, arg_list, add_PG);
    default:
        fprintf(stderr, "[%s] unsupported CRAM version %d\n", __func__,
                cram_major_vers(fd));
        return -1;
    }
}

static void usage(FILE *fp, int ret) {
    fprintf(fp,
           "Usage: samtools reheader [-P] in.header.sam in.bam > out.bam\n"
           "   or  samtools reheader [-P] -i in.header.sam file.bam\n"
           "\n"
           "Options:\n"
           "    -P, --no-PG      Do not generate an @PG header line.\n"
           "    -i, --in-place   Modify the bam/cram file directly.\n"
           "                     (Defaults to outputting to stdout.)\n");
    exit(ret);
}

int main_reheader(int argc, char *argv[])
{
    int inplace = 0, r, add_PG = 1, c;
    bam_hdr_t *h;
    samFile *in;
    char *arg_list = stringify_argv(argc+1, argv-1);

    static const struct option lopts[] = {
        {"help",     no_argument, NULL, 'h'},
        {"in-place", no_argument, NULL, 'i'},
        {"no-PG",    no_argument, NULL, 'P'},
        {NULL, 0, NULL, 0}
    };

    while ((c = getopt_long(argc, argv, "hiP", lopts, NULL)) >= 0) {
        switch (c) {
        case 'P': add_PG = 0; break;
        case 'i': inplace = 1; break;
        case 'h': usage(stdout, 0); break;
        default:
            fprintf(stderr, "Invalid option '%c'\n", c);
            usage(stderr, 1);
        }
    }

    if (argc - optind != 2)
        usage(stderr, 1);

    { // read the header
        samFile *fph = sam_open(argv[optind], "r");
        if (fph == 0) {
            print_error_errno("reheader", "fail to read the header from '%s'", argv[optind]);
            return 1;
        }
        h = sam_hdr_read(fph);
        sam_close(fph);
        if (h == NULL) {
            fprintf(stderr, "[%s] failed to read the header for '%s'.\n",
                    __func__, argv[1]);
            return 1;
        }
    }
    in = sam_open(argv[optind+1], inplace?"r+":"r");
    if (in == 0) {
        print_error_errno("reheader", "fail to open file '%s'", argv[optind+1]);
        return 1;
    }
    if (hts_get_format(in)->format == bam) {
        r = bam_reheader(in->fp.bgzf, h, fileno(stdout), arg_list, add_PG);
    } else {
        if (inplace)
            r = cram_reheader_inplace(in->fp.cram, h, arg_list, add_PG);
        else
            r = cram_reheader(in->fp.cram, h, arg_list, add_PG);
    }

    if (sam_close(in) != 0)
        r = -1;

    bam_hdr_destroy(h);

    if (arg_list)
        free(arg_list);

    return -r;
}
