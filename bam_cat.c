/*  bam_cat.c -- efficiently concatenates bam files.

    Copyright (C) 2008-2009, 2011-2013 Genome Research Ltd.
    Modified SAMtools work copyright (C) 2010 Illumina, Inc.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notices and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

/*
bam_cat can be used to concatenate BAM files. Under special
circumstances, it can be used as an alternative to 'samtools merge' to
concatenate multiple sorted files into a single sorted file. For this
to work each file must be sorted, and the sorted files must be given
as command line arguments in order such that the final read in file i
is less than or equal to the first read in file i+1.

This code is derived from the bam_reheader function in samtools 0.1.8
and modified to perform concatenation by Chris Saunders on behalf of
Illumina.
*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

#include "htslib/bgzf.h"
#include "htslib/sam.h"
#include "htslib/cram.h"
#include "htslib/khash.h"

KHASH_MAP_INIT_STR(s2i, int)

// Bi-directional lookup.
// We can go from name to ID or ID to name.
typedef struct khash_s2i {
    khash_t(s2i) *h;
    int n_id, a_id;
    const char **id; // map Nth entry back to key
    const char **line;
} khash_s2i;

static int hash_s2i_inc(khash_s2i *hash, const char *str, const char *line, int *added) {
    // loosly based on khash_str2int_inc
    khint_t k;
    int n;

    if ( !hash ) return -1;
    // inefficient, but works
    char *my_str = strdup(str);
    k = kh_put(s2i, hash->h, my_str, added);
    if (*added == 0) {
        free(my_str);
        return kh_val(hash->h, k);
    }
    n = hash->n_id++;
    kh_val(hash->h, k) = n;
    if (hash->a_id <= n) {
        const char **id;
        hash->a_id = (n+1)*2;
        if (!(id = realloc(hash->id, hash->a_id*sizeof(*hash->id))))
            return -1;
        hash->id = id;
        if (!(id = realloc(hash->line, hash->a_id*sizeof(*hash->line))))
            return -1;
        hash->line = id;
    }
    hash->id[n] = my_str; // reverse map
    if (line)
        hash->line[n] = line;

    return n;
}

khash_s2i *hash_s2i_create(void) {
    khash_s2i *h = calloc(1, sizeof(*h));
    if (!h)
        return NULL;

    h->h = kh_init(s2i);
    if (!h->h) {
        free(h);
        return NULL;
    }
    return h;
}

static void hash_s2i_free(khash_s2i *hash) {
    // based on khash_str2int_destroy_free
    khint_t k;
    if (!hash) return;
    if (hash->h) {
        for (k = 0; k < kh_end(hash->h); ++k)
            if (kh_exist(hash->h, k)) free((char*)kh_key(hash->h, k));
        kh_destroy(s2i, hash->h);
    }
    if (hash->id)
        free(hash->id);
    if (hash->line)
        free(hash->line);

    free(hash);
}

static khash_s2i *hash_rg(const bam_hdr_t *h) {
    khash_s2i *rg2id = hash_s2i_create();
    char *cp, *line;
    int j, l;

    if (!h)
        return rg2id;

    if (!rg2id)
        return NULL;

    cp = h->text;

    for (l = 0; l+3 < h->l_text; l++) {
        line = &cp[l];
        if (!(cp[l] == '@' && cp[l+1] == 'R' && cp[l+2] == 'G')) {
            while (l < h->l_text && cp[l] != '\n')
                l++;
            continue;
        }

        // Found an @RG line; add to hash
        while (cp[l] != '\n') {
            while (l < h->l_text && cp[l] != '\n' && cp[l] != '\t')
                l++;
            if (l+4 < h->l_text && cp[l+1] == 'I' && cp[l+2] == 'D')
                break;
        }
        if (cp[l] == '\n')
            continue;
        l = (j = l+4);
        while (l < h->l_text && cp[l] != '\n' && cp[l] != '\t')
            l++;

        // To do: save id and keep realloc as needed, as hash_s2i_inc strdups.
        char *id = malloc(l-j+1);
        strncpy(id, &cp[j], l-j);
        id[l-j] = 0;

        int added;
        hash_s2i_inc(rg2id, id, line, &added);
        free(id);

        while (l < h->l_text && cp[l] != '\n')
            l++;
    }

    return rg2id;
}

/*
 * Check the files are consistent and capable of being concatenated.
 * Also fills out the rg2id read-group hash and the version numbers
 * and produces a new bam_hdr_t structure with merged RG lines.
 * Note it is only a simple merge, as we lack the niceties of a proper
 * header API.
 *
 * Returns updated header on success;
 *        NULL on failure.
 */
static bam_hdr_t *cram_cat_check_hdr(int nfn, char * const *fn, const bam_hdr_t *h,
                                     khash_s2i **rg2id, int *vers_maj_p, int *vers_min_p) {
    int i, vers_maj = -1, vers_min = -1;
    bam_hdr_t *new_h = NULL;

    if (h) {
        new_h = bam_hdr_dup(h);
        *rg2id = hash_rg(new_h);
    }

    for (i = 0; i < nfn; ++i) {
        samFile *in;
        cram_fd *in_c;
        khint_t ki;
        int new_rg = -1;

        in = sam_open(fn[i], "rc");
        if (in == 0) {
            fprintf(stderr, "[%s] ERROR: fail to open file '%s'.\n", __func__, fn[i]);
            return NULL;
        }
        in_c = in->fp.cram;

        int vmaj = cram_major_vers(in_c);
        int vmin = cram_minor_vers(in_c);
        if ((vers_maj != -1 && vers_maj != vmaj) ||
            (vers_min != -1 && vers_min != vmin)) {
            fprintf(stderr, "[%s] ERROR: input files have differing version numbers.\n",
                    __func__);
            return NULL;
        }
        vers_maj = vmaj;
        vers_min = vmin;

        bam_hdr_t *old = sam_hdr_read(in);
        khash_s2i *rg2id_in = hash_rg(old);

        if (!new_h) {
            new_h = bam_hdr_dup(old);
            *rg2id = hash_rg(new_h);
        }

        // Add any existing @RG entries to our global @RG hash.
        for (ki = 0; ki < rg2id_in->n_id; ki++) {
            int added;

            new_rg = hash_s2i_inc(*rg2id, rg2id_in->id[ki], rg2id_in->line[ki], &added);
            //fprintf(stderr, "RG %s: #%d -> #%d\n",
            //        rg2id_in->id[ki], ki, new_rg);

            if (added) {
                // Also add to new_h
                const char *line = rg2id_in->line[ki];
                const char *line_end = line;
                while (*line && *line_end++ != '\n')
                    ;
                new_h->l_text += line_end - line;
                new_h->text = realloc(new_h->text, new_h->l_text+1);
                strncat(&new_h->text[new_h->l_text - (line_end - line)],
                        line, line_end - line);
            }

            if (new_rg != ki && rg2id_in->n_id > 1) {
                fprintf(stderr, "[%s] ERROR: Same size @RG lists but differing order / contents\n",
                        __func__);
                return NULL;
            }
        }

        hash_s2i_free(rg2id_in);
        bam_hdr_destroy(old);
        sam_close(in);
    }

    *vers_maj_p = vers_maj;
    *vers_min_p = vers_min;

    return new_h;
}


/*
 * CRAM files don't store the RG:Z:ID per read in the aux field.
 * Instead they have a numerical data series (RG) to point each read
 * back to the Nth @RG line in the file.  This means that we may need
 * to edit the RG data series (if the files were produced from
 * "samtools split" for example).
 *
 * The encoding method is stored in the compression header. Typical
 * examples:
 *
 * RG => EXTERNAL {18}           # Block content-id 18 holds RG values
 *                               # as a series of ITF8 encoded values
 *
 * RG => HUFFMAN {1, 255, 255, 255, 255, 255, 1, 0}
 *                               # One RG value #-1.  (No RG)
 *
 * RG => HUFFMAN {1, 0, 1, 0}    # One RG value #0 (always first RG)
 *
 * RG => HUFFMAN {2, 0, 1, 2, 1, 1}
 *                               # Two RG values, #0 and #1, written
 *                               # to the CORE block and possibly
 *                               # mixed with other data series.
 *
 * A single value can (but may not be) implemented as a zero bit
 * huffman code.  In this situation we can change the meta-data in the
 * compression header to renumber an RG value..
 */
int cram_cat(int nfn, char * const *fn, const bam_hdr_t *h, const char* outcram)
{
    samFile *out;
    cram_fd *out_c;
    int i, vers_maj, vers_min;
    khash_s2i *rg2id = NULL;
    bam_hdr_t *new_h = NULL;

    /* Check consistent versioning and compatible headers */
    if (!(new_h = cram_cat_check_hdr(nfn, fn, h, &rg2id, &vers_maj, &vers_min)))
        return -1;

    /* Open the file with cram_vers */
    char vers[100];
    sprintf(vers, "%d.%d", vers_maj, vers_min);
    out = sam_open(outcram, "wc");
    if (out == 0) {
        fprintf(stderr, "[%s] ERROR: fail to open output file '%s'.\n", __func__, outcram);
        return 1;
    }
    out_c = out->fp.cram;
    cram_set_option(out_c, CRAM_OPT_VERSION, vers);
    //fprintf(stderr, "Creating cram vers %s\n", vers);

    cram_fd_set_header(out_c, sam_hdr_parse_(new_h->text,  new_h->l_text)); // needed?
    sam_hdr_write(out, new_h);

    for (i = 0; i < nfn; ++i) {
        samFile *in;
        cram_fd *in_c;
        cram_container *c;
        bam_hdr_t *old;
        int new_rg = -1;

        in = sam_open(fn[i], "rc");
        if (in == 0) {
            fprintf(stderr, "[%s] ERROR: fail to open file '%s'.\n", __func__, fn[i]);
            return -1;
        }
        in_c = in->fp.cram;

        old = sam_hdr_read(in);
        khash_s2i *rg2id_in = hash_rg(old);

        // Compute RG mapping if suitable for changing.
        if (rg2id_in->n_id == 1) {
            int _;
            new_rg = hash_s2i_inc(rg2id, rg2id_in->id[0], NULL, &_);
        } else {
            new_rg = 0;
        }

        hash_s2i_free(rg2id_in);


        // Copy contains and blocks within them
        while ((c = cram_read_container(in_c))) {
            cram_block *blk;

           if (cram_container_is_empty(in_c)) {
                if (cram_write_container(out_c, c) != 0)
                    return -1;

                // Container compression header
                if (!(blk = cram_read_block(in_c)))
                    return -1;
                if (cram_write_block(out_c, blk) != 0) {
                    cram_free_block(blk);
                    return -1;
                }
                cram_free_block(blk);
                cram_free_container(c);

                continue;
            }

            // If we have just one RG key and new_rg != 0 then
            // we need to edit the compression header. IF WE CAN.
            if (new_rg) {
                int zero = 0;
                //fprintf(stderr, "Transcode RG %d to %d\n", 0, new_rg);
                cram_transcode_rg(in_c, out_c, c, 1, &zero, &new_rg);
            } else {
                int32_t num_slices;

                // Not switching rg so do the usual read/write loop
                if (cram_write_container(out_c, c) != 0)
                    return -1;

                // Container compression header
                if (!(blk = cram_read_block(in_c)))
                    return -1;
                if (cram_write_block(out_c, blk) != 0) {
                    cram_free_block(blk);
                    return -1;
                }
                cram_free_block(blk);


                // Container num_blocks can be invalid, due to a bug.
                // Instead we iterate in slice context instead.
                (void)cram_container_get_landmarks(c, &num_slices);
                cram_copy_slice(in_c, out_c, num_slices);
            }

            cram_free_container(c);
        }

        bam_hdr_destroy(old);
        sam_close(in);
    }
    sam_close(out);

    hash_s2i_free(rg2id);
    bam_hdr_destroy(new_h);

    return 0;
}


#define BUF_SIZE 0x10000

#define GZIPID1 31
#define GZIPID2 139

#define BGZF_EMPTY_BLOCK_SIZE 28

int bam_cat(int nfn, char * const *fn, const bam_hdr_t *h, const char* outbam)
{
    BGZF *fp;
    uint8_t *buf;
    uint8_t ebuf[BGZF_EMPTY_BLOCK_SIZE];
    const int es=BGZF_EMPTY_BLOCK_SIZE;
    int i;

    fp = strcmp(outbam, "-")? bgzf_open(outbam, "w") : bgzf_fdopen(fileno(stdout), "w");
    if (fp == 0) {
        fprintf(stderr, "[%s] ERROR: fail to open output file '%s'.\n", __func__, outbam);
        return 1;
    }
    if (h) bam_hdr_write(fp, h);

    buf = (uint8_t*) malloc(BUF_SIZE);
    for(i = 0; i < nfn; ++i){
        BGZF *in;
        bam_hdr_t *old;
        int len,j;

        in = strcmp(fn[i], "-")? bgzf_open(fn[i], "r") : bgzf_fdopen(fileno(stdin), "r");
        if (in == 0) {
            fprintf(stderr, "[%s] ERROR: fail to open file '%s'.\n", __func__, fn[i]);
            return -1;
        }
        if (in->is_write) return -1;

        old = bam_hdr_read(in);
        if (old == NULL) {
            fprintf(stderr, "[%s] ERROR: couldn't read header for '%s'.\n",
                    __func__, fn[i]);
            bgzf_close(in);
            return -1;
        }
        if (h == 0 && i == 0) bam_hdr_write(fp, old);

        if (in->block_offset < in->block_length) {
            bgzf_write(fp, in->uncompressed_block + in->block_offset, in->block_length - in->block_offset);
            bgzf_flush(fp);
        }

        j=0;
        while ((len = bgzf_raw_read(in, buf, BUF_SIZE)) > 0) {
            if(len<es){
                int diff=es-len;
                if(j==0) {
                    fprintf(stderr, "[%s] ERROR: truncated file?: '%s'.\n", __func__, fn[i]);
                    return -1;
                }
                bgzf_raw_write(fp, ebuf, len);
                memcpy(ebuf,ebuf+len,diff);
                memcpy(ebuf+diff,buf,len);
            } else {
                if(j!=0) bgzf_raw_write(fp, ebuf, es);
                len-= es;
                memcpy(ebuf,buf+len,es);
                bgzf_raw_write(fp, buf, len);
            }
            j=1;
        }

        /* check final gzip block */
        {
            const uint8_t gzip1=ebuf[0];
            const uint8_t gzip2=ebuf[1];
            const uint32_t isize=*((uint32_t*)(ebuf+es-4));
            if(((gzip1!=GZIPID1) || (gzip2!=GZIPID2)) || (isize!=0)) {
                fprintf(stderr, "[%s] WARNING: Unexpected block structure in file '%s'.", __func__, fn[i]);
                fprintf(stderr, " Possible output corruption.\n");
                bgzf_raw_write(fp, ebuf, es);
            }
        }
        bam_hdr_destroy(old);
        bgzf_close(in);
    }
    free(buf);
    bgzf_close(fp);
    return 0;
}


int main_cat(int argc, char *argv[])
{
    bam_hdr_t *h = 0;
    char *outfn = 0;
    int c, ret;
    samFile *in;

    while ((c = getopt(argc, argv, "h:o:")) >= 0) {
        switch (c) {
            case 'h': {
                samFile *fph = sam_open(optarg, "r");
                if (fph == 0) {
                    fprintf(stderr, "[%s] ERROR: fail to read the header from '%s'.\n", __func__, argv[1]);
                    return 1;
                }
                h = sam_hdr_read(fph);
                if (h == NULL) {
                    fprintf(stderr,
                            "[%s] ERROR: failed to read the header for '%s'.\n",
                            __func__, argv[1]);
                    return 1;
                }
                sam_close(fph);
                break;
            }
            case 'o': outfn = strdup(optarg); break;
        }
    }
    if (argc - optind < 1) {
        fprintf(stderr, "Usage: samtools cat [-h header.sam] [-o out.bam] <in1.bam> [...]\n");
        return 1;
    }

    in = sam_open(argv[optind], "r");
    if (!in) {
        fprintf(stderr, "[%s] ERROR: failed to open file '%s'.\n", __func__, argv[optind]);
        return 1;
    }

    switch (hts_get_format(in)->format) {
    case bam:
        sam_close(in);
        ret = bam_cat(argc - optind, argv + optind, h, outfn? outfn : "-");
        break;

    case cram:
        sam_close(in);
        ret = cram_cat(argc - optind, argv + optind, h, outfn? outfn : "-");
        break;

    default:
        sam_close(in);
        fprintf(stderr, "[%s] ERROR: input is not BAM or CRAM\n", __func__);
        return 1;
    }
    free(outfn);

    if (h)
        bam_hdr_destroy(h);

    return ret;
}
