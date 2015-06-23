/*  bam_reheader.c -- reheader subcommand.

    Copyright (C) 2010 Broad Institute.
    Copyright (C) 2012, 2013 Genome Research Ltd.

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

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <getopt.h>

#include "htslib/bgzf.h"
#include "htslib/sam.h"
#include "htslib/hfile.h"
#include "cram/cram.h"
#include "version.h"

#define BUF_SIZE 0x10000

/*
 * Reads a file and outputs a new BAM file to fd with 'h' replaced as
 * the header.    No checks are made to the validity.
 */
int bam_reheader(BGZF *in, bam_hdr_t *h, int fd,
                 const char *arg_list, int add_PG)
{
    BGZF *fp;
    ssize_t len;
    uint8_t *buf;
    if (in->is_write) return -1;
    buf = malloc(BUF_SIZE);
    (void)bam_hdr_read(in);
    fp = bgzf_fdopen(fd, "w");

    if (add_PG) {
        // Around the houses, but it'll do until we can manipulate bam_hdr_t natively.
        SAM_hdr *sh = sam_hdr_parse_(h->text, h->l_text);
        if (sam_hdr_add_PG(sh, "samtools",
                           "VN", SAMTOOLS_VERSION,
                           arg_list ? "CL": NULL,
                           arg_list ? arg_list : NULL,
                           NULL) != 0)
            return -1;

        free(h->text);
        h->text = strdup(sam_hdr_str(sh));
        h->l_text = sam_hdr_length(sh);
        if (!h->text)
            return -1;
        sam_hdr_free(sh);
    }

    bam_hdr_write(fp, h);
    if (in->block_offset < in->block_length) {
        bgzf_write(fp, in->uncompressed_block + in->block_offset, in->block_length - in->block_offset);
        bgzf_flush(fp);
    }
    while ((len = bgzf_raw_read(in, buf, BUF_SIZE)) > 0)
        bgzf_raw_write(fp, buf, len);
    free(buf);
    fp->block_offset = in->block_offset = 0;
    bgzf_close(fp);
    return 0;
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
    out->header = sam_hdr_parse_(h->text, h->l_text);
    if (add_PG) {
        if (sam_hdr_add_PG(out->header, "samtools",
                           "VN", SAMTOOLS_VERSION,
                           arg_list ? "CL": NULL,
                           arg_list ? arg_list : NULL,
                           NULL) != 0)
            goto err;

        // Covert back to bam_hdr_t struct
        free(h->text);
        h->text = strdup(sam_hdr_str(out->header));
        h->l_text = sam_hdr_length(out->header);
        if (!h->text)
            goto err;
    }

    if (sam_hdr_write(h_out, h) != 0)
        goto err;
    cram_set_option(out, CRAM_OPT_REFERENCE, NULL);

    while ((c = cram_read_container(in))) {
        int i;
        if (cram_write_container(out, c) != 0)
            goto err;

        for (i = 0; i < c->num_blocks; i++) {
            cram_block *blk = cram_read_block(in);
            if (!blk || cram_write_block(out, blk) != 0) {
                if (blk) cram_free_block(blk);
                goto err;
            }
            cram_free_block(blk);
        }
    }

    ret = 0;

 err:
    if (hts_close(h_out) != 0)
        ret = -1;
    if (c) cram_free_container(c);

    return ret;
}


/*
 * CRAM manipulation functions.  A case could be made for moving these
 * to htslib, especially as these expose internal knowledge of the
 * CRAM encodings.
 *
 * However for simplicity and the fact they're needed by precisely one
 * tool (currently) we'll keep them here for now.  Shout if you need
 * them elsewhere.
 */

/* MAXIMUM storage size needed for the container. */
static inline int cram_container_size(cram_container *c) {
    return 55 + 5*c->num_landmarks;
}

/*
 * Stores the container structure in dat and returns *size as the
 * number of bytes written to dat[].  The input size of dat is also
 * held in *size and should be initialised to cram_container_size(c).
 *
 * Returns 0 on success;
 *        -1 on failure
 */
static int cram_store_container(cram_fd *fd, cram_container *c, char *dat,
                                int *size) {
    char *cp = dat;
    int i;
    
    // Check the input buffer is large enough according to our stated
    // requirements. (NOTE: it may actually take less.)
    if (cram_container_size(c) > *size)
        return -1;

    if (CRAM_MAJOR_VERS(fd->version) == 1) {
	cp += itf8_put(cp, c->length);
    } else {
	*(int32_t *)cp = le_int4(c->length);
	cp += 4;
    }
    if (c->multi_seq) {
	cp += itf8_put(cp, -2);
	cp += itf8_put(cp, 0);
	cp += itf8_put(cp, 0);
    } else {
	cp += itf8_put(cp, c->ref_seq_id);
	cp += itf8_put(cp, c->ref_seq_start);
	cp += itf8_put(cp, c->ref_seq_span);
    }
    cp += itf8_put(cp, c->num_records);
    if (CRAM_MAJOR_VERS(fd->version) == 2) {
	cp += itf8_put(cp, c->record_counter);
	cp += ltf8_put(cp, c->num_bases);
    } else if (CRAM_MAJOR_VERS(fd->version) >= 3) {
	cp += ltf8_put(cp, c->record_counter);
	cp += ltf8_put(cp, c->num_bases);
    }

    cp += itf8_put(cp, c->num_blocks);
    cp += itf8_put(cp, c->num_landmarks);
    for (i = 0; i < c->num_landmarks; i++)
	cp += itf8_put(cp, c->landmark[i]);

    if (CRAM_MAJOR_VERS(fd->version) >= 3) {
	c->crc32 = crc32(0L, (uc *)dat, cp-dat);
	cp[0] =  c->crc32        & 0xff;
	cp[1] = (c->crc32 >>  8) & 0xff;
	cp[2] = (c->crc32 >> 16) & 0xff;
	cp[3] = (c->crc32 >> 24) & 0xff;
	cp += 4;
    }

    *size = cp-dat; // actual used size

    return 0;
}


/*
 * Computes the size of a cram block, including the block
 * header itself.
 */
uint32_t cram_block_size(cram_block *b) {
    unsigned char dat[100], *cp = dat;;
    uint32_t sz;

    *cp++ = b->method;
    *cp++ = b->content_type;
    cp += itf8_put(cp, b->content_id);
    cp += itf8_put(cp, b->comp_size);
    cp += itf8_put(cp, b->uncomp_size);

    sz = cp-dat + 4;
    sz += b->method == RAW ? b->uncomp_size : b->comp_size;

    return sz;
}


/*
 * Reads a version 2 CRAM file and replaces the header in situ,
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
int cram_reheader_insitu2(cram_fd *fd, const bam_hdr_t *h, const char *arg_list,
                          int add_PG)
{
    cram_container *c = NULL;
    cram_block *b = NULL;
    SAM_hdr *hdr = NULL;
    off_t start;
    int ret = -1;

    if (CRAM_MAJOR_VERS(fd->version) < 2 ||
        CRAM_MAJOR_VERS(fd->version) > 3) {
        fprintf(stderr, "[%s] unsupported CRAM version %d\n", __func__,
                CRAM_MAJOR_VERS(fd->version));
        goto err;
    }

    if (!(hdr = sam_hdr_parse_(h->text, h->l_text)))
        goto err;

    if (add_PG && sam_hdr_add_PG(hdr, "samtools", "VN", SAMTOOLS_VERSION,
                                 arg_list ? "CL": NULL,
                                 arg_list ? arg_list : NULL,
                                 NULL))
        goto err;

    int header_len = sam_hdr_length(hdr);
    /* Fix M5 strings? Maybe out of scope for this tool */

    // Load the existing header
    if ((start = hseek(fd->fp, 26, SEEK_SET)) != 26)
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

    if (b->uncomp_size < header_len+4) {
        fprintf(stderr, "New header will not fit. Use non-insitu version (%d > %d)\n",
                header_len+4, b->uncomp_size);
        ret = -2;
        goto err;
    }

    BLOCK_SIZE(b) = 0;
    int32_put_blk(b, header_len);
    BLOCK_APPEND(b, sam_hdr_str(hdr), header_len);
    memset(BLOCK_DATA(b)+BLOCK_SIZE(b), 0, b->uncomp_size - BLOCK_SIZE(b));
    BLOCK_SIZE(b) = b->uncomp_size;
    BLOCK_UPLEN(b);

    if (hseek(fd->fp, start, SEEK_SET) != start)
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
 * Reads a version 3 CRAM file and replaces the header in situ,
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
int cram_reheader_insitu3(cram_fd *fd, const bam_hdr_t *h, const char *arg_list,
                          int add_PG)
{
    cram_container *c = NULL;
    cram_block *b = NULL;
    SAM_hdr *hdr;
    off_t start, sz, end;
    int container_sz, max_container_sz;
    char *buf = NULL;
    int ret = -1;

    if (CRAM_MAJOR_VERS(fd->version) < 2 ||
        CRAM_MAJOR_VERS(fd->version) > 3) {
        fprintf(stderr, "[%s] unsupported CRAM version %d\n", __func__,
                CRAM_MAJOR_VERS(fd->version));
        goto err;
    }

    if (!(hdr = sam_hdr_parse_(h->text, h->l_text)))
        goto err;

    if (add_PG && sam_hdr_add_PG(hdr, "samtools", "VN", SAMTOOLS_VERSION,
                                 arg_list ? "CL": NULL,
                                 arg_list ? arg_list : NULL,
                                 NULL))
        goto err;

    int header_len = sam_hdr_length(hdr);
    /* Fix M5 strings? Maybe out of scope for this tool */

    // Find current size of SAM header block
    if ((start = hseek(fd->fp, 26, SEEK_SET)) != 26)
        goto err;

    if (!(c = cram_read_container(fd)))
        goto err;

    // +5 allows num_landmarks to increase from 0 to 1 (Cramtools)
    max_container_sz = cram_container_size(c)+5; 

    sz = htell(fd->fp) + c->length - start;
    end = htell(fd->fp) + c->length;

    // We force 1 block instead of (optionally) 2.  C CRAM
    // implementations for v3 were writing 1 compressed block followed
    // by 1 uncompressed block.  However this is tricky to deal with
    // as changing block sizes can mean the block header also changes
    // size due to itf8 and variable size integers.
    //
    // If we had 1 block, this doesn't change anything.
    // If we had 2 blocks, the new container header will be smaller by
    // 1+ bytes, requiring the c->length to be larger in value.
    // However this is an int32 instead of itf8 so the container
    // header structure stays the same size.  This means we can always
    // reduce the number of blocks without running into size problems.
    c->num_blocks = 1;
    if (c->num_landmarks && c->landmark) {
        c->num_landmarks = 1;
        c->landmark[0] = 0;
    } else {
        c->num_landmarks = 0;
    }

    buf = malloc(max_container_sz);
    container_sz = max_container_sz;
    if (cram_store_container(fd, c, buf, &container_sz) != 0)
        goto err;

    if (!buf)
        goto err;

    // Proposed new length, but changing c->length may change the
    // container_sz and thus the remainder (c->length itself).
    c->length = sz - container_sz;
    
    int old_container_sz = container_sz;
    container_sz = max_container_sz;
    if (cram_store_container(fd, c, buf, &container_sz) != 0)
        goto err;

    if (old_container_sz != container_sz) {
        fprintf(stderr, "Quirk of fate makes this troublesome! "
                "Please use non-insitu version.\n");
        goto err;
    }
    
    

    // Version 3.0 supports compressed header
    b = cram_new_block(FILE_HEADER, 0);
    int32_put_blk(b, header_len);
    BLOCK_APPEND(b, sam_hdr_str(hdr), header_len);
    BLOCK_UPLEN(b);

    if (fd->level > 0) {
        int method = 1<<GZIP;
        if (fd->use_bz2)
            method |= 1<<BZIP2;
        if (fd->use_lzma)
            method |= 1<<LZMA;
        cram_compress_block(fd, b, NULL, method, fd->level);
    }

    if (hseek(fd->fp, 26, SEEK_SET) != 26)
        goto err;

    if (cram_block_size(b) > c->length) {
        fprintf(stderr, "New header will not fit. Use non-insitu version"
                " (%d > %d)\n",
                cram_block_size(b), c->length);
        ret = -2;
        goto err;
    }

    if (cram_write_container(fd, c) == -1)
        goto err;

    if (cram_write_block(fd, b) == -1)
        goto err;

    // Blank out the remainder
    int rsz = end - htell(fd->fp);
    assert(rsz >= 0);
    if (rsz) {
        char *rem = calloc(1, rsz);
        ret = hwrite(fd->fp, rem, rsz) == rsz ? 0 : -1;
        free(rem);
    }

 err:
    if (c) cram_free_container(c);
    if (buf) free(buf);
    if (b) cram_free_block(b);
    if (hdr) sam_hdr_free(hdr);

    return ret;
}

int cram_reheader_insitu(cram_fd *fd, const bam_hdr_t *h, const char *arg_list,
                         int add_PG)
{
    switch (CRAM_MAJOR_VERS(fd->version)) {
    case 2: return cram_reheader_insitu2(fd, h, arg_list, add_PG);
    case 3: return cram_reheader_insitu3(fd, h, arg_list, add_PG);
    default:
        fprintf(stderr, "[%s] unsupported CRAM version %d\n", __func__,
                CRAM_MAJOR_VERS(fd->version));
        return -1;
    }
}

static void usage(int ret) {
    printf("Usage: samtools reheader [-P] in.header.sam in.bam > out.bam\n"
           "   or  samtools reheader [-P] -I in.header.sam file.bam\n"
           "\n"
           "Options:\n"
           "    -P, --no-PG      Do not generate an @PG header line.\n"
           "    -I, --insitu     Modify the bam/cram file directly.\n"
           "                     (Defaults to outputting to stdout.)\n");
    exit(ret);
}

int main_reheader(int argc, char *argv[])
{
    int insitu = 0, r, add_PG = 1, c;
    bam_hdr_t *h;
    samFile *in;
    char *arg_list = stringify_argv(argc+1, argv-1);

    static const struct option lopts[] = {
        {"insitu", no_argument, NULL, 'I'},
        {"no-PG",  no_argument, NULL, 'P'},
        {NULL, 0, NULL, 0}
    };

    while ((c = getopt_long(argc, argv, "hIP", lopts, NULL)) >= 0) {
        switch (c) {
        case 'P': add_PG = 0; break;
        case 'I': insitu = 1; break;
        case 'h': usage(0);
        default:
            fprintf(stderr, "Invalid option '%c'\n", c);
            usage(1);
        }
    }

    if (argc - optind != 2)
        usage(1);

    { // read the header
        samFile *fph = sam_open(argv[optind], "r");
        if (fph == 0) {
            fprintf(stderr, "[%s] fail to read the header from %s.\n", __func__, argv[optind]);
            return 1;
        }
        h = sam_hdr_read(fph);
        sam_close(fph);
    }
    in = sam_open(argv[optind+1], insitu?"r+":"r");
    if (in == 0) {
        fprintf(stderr, "[%s] fail to open file %s.\n", __func__, argv[optind+1]);
        return 1;
    }
    if (in->format.format == bam) {
        r = bam_reheader(in->fp.bgzf, h, fileno(stdout), arg_list, add_PG);
    } else {
        if (insitu)
            r = cram_reheader_insitu(in->fp.cram, h, arg_list, add_PG);
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
