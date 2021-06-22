/*  bam_cat.c -- efficiently concatenates bam files.

    Copyright (C) 2008-2009, 2011-2013, 2015-2017, 2019, 2021 Genome Research Ltd.
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

#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <strings.h>

#include "htslib/bgzf.h"
#include "htslib/sam.h"
#include "htslib/cram.h"
#include "htslib/kstring.h"
#include "samtools.h"
#include "sam_opts.h"

/*
 * Check the files are consistent and capable of being concatenated.
 * Also fills out the version numbers and produces a new sam_hdr_t
 * structure with merged RG lines.
 * Note it is only a simple merge.
 *
 * Returns updated header on success;
 *        NULL on failure.
 */
static sam_hdr_t *cram_cat_check_hdr(int nfn, char * const *fn, const sam_hdr_t *h,
                                     int *vers_maj_p, int *vers_min_p) {
    int i, vers_maj = -1, vers_min = -1;
    sam_hdr_t *new_h = NULL, *old_h = NULL;
    samFile *in = NULL;
    kstring_t ks = KS_INITIALIZE;

    if (h) {
        new_h = sam_hdr_dup(h);
        if (!new_h) {
            fprintf(stderr, "[%s] ERROR: header duplication failed.\n",
                    __func__);
            goto fail;
        }
    }

    for (i = 0; i < nfn; ++i) {
        cram_fd *in_c;
        int ki;

        in = sam_open(fn[i], "rc");
        if (in == 0) {
            print_error_errno("cat", "fail to open file '%s'", fn[i]);
            goto fail;
        }
        in_c = in->fp.cram;

        int vmaj = cram_major_vers(in_c);
        int vmin = cram_minor_vers(in_c);
        if ((vers_maj != -1 && vers_maj != vmaj) ||
            (vers_min != -1 && vers_min != vmin)) {
            fprintf(stderr, "[%s] ERROR: input files have differing version numbers.\n",
                    __func__);
            goto fail;
        }
        vers_maj = vmaj;
        vers_min = vmin;

        old_h = sam_hdr_read(in);
        if (!old_h) {
            fprintf(stderr, "[%s] ERROR: header reading for file '%s' filed.\n",
                    __func__, fn[i]);
            goto fail;
        }

        if (!new_h) {
            new_h = sam_hdr_dup(old_h);
            if (!new_h) {
                fprintf(stderr, "[%s] ERROR: header duplication for file '%s' failed.\n",
                        __func__, fn[i]);
                goto fail;
            }
            sam_hdr_destroy(old_h);
            sam_close(in);
            continue;
        }

        int old_count = sam_hdr_count_lines(old_h, "RG");
        for (ki = 0; ki < old_count; ki++) {
            const char *old_name = sam_hdr_line_name(old_h, "RG", ki);
            if (old_name) {
                int new_i = sam_hdr_line_index(new_h, "RG", old_name);
                if (-1 == new_i) { // line does not exist in the new header
                    if (sam_hdr_find_line_pos(old_h, "RG", ki, &ks) ||
                        !ks.s || sam_hdr_add_lines(new_h, ks.s, ks.l)) {
                        fprintf(stderr, "[%s] ERROR: failed to add @RG line 'ID:%s' from file '%s'\n",
                                __func__, old_name, fn[i]);
                        goto fail;
                    }
                    ks_free(&ks);
                }
            } else {
                fprintf(stderr, "[%s] ERROR: failed to read %d @RG line from file '%s'\n",
                        __func__, ki, fn[i]);
                goto fail;
            }
        }

        if (old_count > 1 && sam_hdr_count_lines(new_h, "RG") == old_count) {
            for (ki = 0; ki < old_count; ki++) {
                const char *old_name = sam_hdr_line_name(old_h, "RG", ki);
                const char *new_name = sam_hdr_line_name(new_h, "RG", ki);
                if (!old_name || !new_name || strcmp(old_name, new_name)) {
                    fprintf(stderr, "[%s] ERROR: Same size @RG lists but differing order / contents\n",
                            __func__);
                    goto fail;
                }
            }
        }

        sam_hdr_destroy(old_h);
        sam_close(in);
    }

    ks_free(&ks);

    *vers_maj_p = vers_maj;
    *vers_min_p = vers_min;

    return new_h;

fail:
    ks_free(&ks);
    if (old_h) sam_hdr_destroy(old_h);
    if (new_h) sam_hdr_destroy(new_h);
    if (in) sam_close(in);

    return NULL;
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
int cram_cat(int nfn, char * const *fn, const sam_hdr_t *h, const char* outcram, sam_global_args *ga, char *arg_list, int no_pg)
{
    samFile *out;
    cram_fd *out_c;
    int i, vers_maj, vers_min;
    sam_hdr_t *new_h = NULL;

    /* Check consistent versioning and compatible headers */
    if (!(new_h = cram_cat_check_hdr(nfn, fn, h, &vers_maj, &vers_min)))
        return -1;

    /* Open the file with cram_vers */
    char vers[100];
    sprintf(vers, "%d.%d", vers_maj, vers_min);
    out = sam_open_format(outcram, "wc", &ga->out);
    if (out == 0) {
        print_error_errno("cat", "fail to open output file '%s'", outcram);
        return -1;
    }
    out_c = out->fp.cram;
    cram_set_option(out_c, CRAM_OPT_VERSION, vers);
    //fprintf(stderr, "Creating cram vers %s\n", vers);

    if (!no_pg && sam_hdr_add_pg(new_h, "samtools",
                                 "VN", samtools_version(),
                                 arg_list ? "CL": NULL,
                                 arg_list ? arg_list : NULL,
                                 NULL))
        return -1;

    if (sam_hdr_write(out, new_h) < 0) {
        print_error_errno("cat", "Couldn't write header");
        return -1;
    }

    for (i = 0; i < nfn; ++i) {
        samFile *in;
        cram_fd *in_c;
        cram_container *c;
        sam_hdr_t *old_h;
        int new_rg = -1;

        in = sam_open(fn[i], "rc");
        if (in == 0) {
            print_error_errno("cat", "fail to open file '%s'", fn[i]);
            return -1;
        }
        in_c = in->fp.cram;

        old_h = sam_hdr_read(in);
        if (!old_h) {
            print_error("cat", "fail to read the header of file '%s'", fn[i]);
            return -1;
        }

        // Compute RG mapping if suitable for changing.
        if (sam_hdr_count_lines(old_h, "RG") == 1) {
            const char *old_name = sam_hdr_line_name(old_h, "RG", 0);
            if (old_name) {
                new_rg = sam_hdr_line_index(new_h, "RG", old_name);
                if (new_rg < 0) {
                    print_error("cat", "fail to find @RG line '%s' in the new header", old_name);
                    return -1;
                }
            } else {
                print_error("cat", "fail to find @RG line in file '%s'", fn[i]);
                return -1;
            }
        } else {
            new_rg = 0;
        }

        // Copy contains and blocks within them
        while ((c = cram_read_container(in_c))) {
            if (cram_container_is_empty(in_c)) {
                cram_block *blk;
                // Container compression header
                if (!(blk = cram_read_block(in_c)))
                    return -1;
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
                cram_block *blk;

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

        sam_hdr_destroy(old_h);
        sam_close(in);
    }
    sam_close(out);
    sam_hdr_destroy(new_h);

    return 0;
}


#define BUF_SIZE 0x10000

#define GZIPID1 31
#define GZIPID2 139

#define BGZF_EMPTY_BLOCK_SIZE 28

int bam_cat(int nfn, char * const *fn, sam_hdr_t *h, const char* outbam, char *arg_list, int no_pg)
{
    BGZF *fp, *in = NULL;
    uint8_t *buf = NULL;
    uint8_t ebuf[BGZF_EMPTY_BLOCK_SIZE];
    const int es=BGZF_EMPTY_BLOCK_SIZE;
    int i;

    fp = strcmp(outbam, "-")? bgzf_open(outbam, "w") : bgzf_fdopen(fileno(stdout), "w");
    if (fp == 0) {
        print_error_errno("cat", "fail to open output file '%s'", outbam);
        return -1;
    }
    if (h) {
        if (!no_pg && sam_hdr_add_pg(h, "samtools",
                                     "VN", samtools_version(),
                                     arg_list ? "CL": NULL,
                                     arg_list ? arg_list : NULL,
                                     NULL))
            goto fail;

        if (bam_hdr_write(fp, h) < 0) {
            print_error_errno("cat", "Couldn't write header");
            goto fail;
        }
    }

    buf = (uint8_t*) malloc(BUF_SIZE);
    if (!buf) {
        fprintf(stderr, "[%s] Couldn't allocate buffer\n", __func__);
        goto fail;
    }
    for(i = 0; i < nfn; ++i){
        sam_hdr_t *old;
        int len,j;

        in = strcmp(fn[i], "-")? bgzf_open(fn[i], "r") : bgzf_fdopen(fileno(stdin), "r");
        if (in == 0) {
            print_error_errno("cat", "fail to open file '%s'", fn[i]);
            goto fail;
        }
        if (in->is_write) return -1;

        old = bam_hdr_read(in);
        if (old == NULL) {
            fprintf(stderr, "[%s] ERROR: couldn't read header for '%s'.\n",
                    __func__, fn[i]);
            goto fail;
        }
        if (h == 0 && i == 0) {
            if (!no_pg && sam_hdr_add_pg(old, "samtools",
                                         "VN", samtools_version(),
                                         arg_list ? "CL": NULL,
                                         arg_list ? arg_list : NULL,
                                         NULL))
                goto fail;

            if (bam_hdr_write(fp, old) < 0) {
                print_error_errno("cat", "Couldn't write header");
                goto fail;
            }
        }

        if (in->block_offset < in->block_length) {
            if (bgzf_write(fp, (char *)in->uncompressed_block + in->block_offset, in->block_length - in->block_offset) < 0) goto write_fail;
            if (bgzf_flush(fp) != 0) goto write_fail;
        }

        j=0;
        while ((len = bgzf_raw_read(in, buf, BUF_SIZE)) > 0) {
            if(len<es){
                int diff=es-len;
                if(j==0) {
                    fprintf(stderr, "[%s] ERROR: truncated file?: '%s'.\n", __func__, fn[i]);
                    goto fail;
                }
                if (bgzf_raw_write(fp, ebuf, len) < 0) goto write_fail;

                memcpy(ebuf,ebuf+len,diff);
                memcpy(ebuf+diff,buf,len);
            } else {
                if(j!=0) {
                    if (bgzf_raw_write(fp, ebuf, es) < 0) goto write_fail;
                }
                len-= es;
                memcpy(ebuf,buf+len,es);
                if (bgzf_raw_write(fp, buf, len) < 0) goto write_fail;
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
                if (bgzf_raw_write(fp, ebuf, es) < 0) goto write_fail;
            }
        }
        sam_hdr_destroy(old);
        bgzf_close(in);
        in = NULL;
    }
    free(buf);
    if (bgzf_close(fp) < 0) {
        fprintf(stderr, "[%s] Error on closing '%s'.\n", __func__, outbam);
        return -1;
    }
    return 0;

 write_fail:
    fprintf(stderr, "[%s] Error writing to '%s'.\n", __func__, outbam);
 fail:
    if (in) bgzf_close(in);
    if (fp) bgzf_close(fp);
    free(buf);
    return -1;
}


int main_cat(int argc, char *argv[])
{
    sam_hdr_t *h = 0;
    char *outfn = 0;
    char **infns = NULL; // files to concatenate
    int infns_size = 0;
    int c, ret = 0, no_pg = 0, usage = 0;
    samFile *in;
    sam_global_args ga;

    static const struct option lopts[] = {
        SAM_OPT_GLOBAL_OPTIONS('-', '-', '-', 0, '-', '@'),
        {"no-PG", no_argument, NULL, 1},
        { NULL, 0, NULL, 0 }
    };

    char *arg_list = NULL;

    sam_global_args_init(&ga);

    while ((c = getopt_long(argc, argv, "h:o:b:@:", lopts, NULL)) >= 0) {
        switch (c) {
            case 'h': {
                samFile *fph = sam_open(optarg, "r");
                if (fph == 0) {
                    fprintf(stderr, "[%s] ERROR: fail to read the header from '%s'.\n", __func__, optarg);
                    return 1;
                }
                h = sam_hdr_read(fph);
                if (h == NULL) {
                    fprintf(stderr,
                            "[%s] ERROR: failed to read the header from '%s'.\n",
                            __func__, optarg);
                    return 1;
                }
                sam_close(fph);
                break;
            }
            case 'o': outfn = strdup(optarg); break;
            case 'b': {
                // add file names in "optarg" to the list
                // of files to concatenate
                int nfns;
                char **fns_read = hts_readlines(optarg, &nfns);
                if (fns_read) {
                    infns = realloc(infns, (infns_size + nfns) * sizeof(char*));
                    if (infns == NULL) { ret = 1; goto end; }
                    memcpy(infns+infns_size, fns_read, nfns * sizeof(char*));
                    infns_size += nfns;
                    free(fns_read);
                } else {
                    print_error("cat", "Invalid file list \"%s\"", optarg);
                    ret = 1;
                }
                break;
            }
            case 1:
                no_pg = 1;
                break;
            default:
                if (parse_sam_global_opt(c, optarg, lopts, &ga) == 0) break;
                /* else fall-through */
            case '?': usage=1; break;
        }
    }

    if (!no_pg && !(arg_list = stringify_argv(argc+1, argv-1))) {
        print_error("cat", "failed to create arg_list");
        return 1;
    }

    // Append files specified in argv to the list.
    int nargv_fns = argc - optind;
    if (nargv_fns > 0) {
        infns = realloc(infns, (infns_size + nargv_fns) * sizeof(char*));
        if (infns == NULL) { ret = 1; goto end; }
        memcpy(infns + infns_size, argv + optind, nargv_fns * sizeof(char*));
    }

    // Require at least one input file
    if (infns_size + nargv_fns == 0 || usage) {
        fprintf(stderr, "Usage: samtools cat [options] <in1.bam>  [... <inN.bam>]\n");
        fprintf(stderr, "       samtools cat [options] <in1.cram> [... <inN.cram>]\n\n");
        fprintf(stderr, "Concatenate BAM or CRAM files, first those in <bamlist.fofn>, then those\non the command line.\n\n");
        fprintf(stderr, "Options: -b FILE  list of input BAM/CRAM file names, one per line\n");
        fprintf(stderr, "         -h FILE  copy the header from FILE [default is 1st input file]\n");
        fprintf(stderr, "         -o FILE  output BAM/CRAM\n");
        fprintf(stderr, "         --no-PG  do not add a PG line\n");
        sam_global_opt_help(stderr, "--..-@-.");
        return 1;
    }

    in = sam_open(infns[0], "r");
    if (!in) {
        print_error_errno("cat", "failed to open file '%s'", infns[0]);
        return 1;
    }

    switch (hts_get_format(in)->format) {
    case bam:
        sam_close(in);
        if (bam_cat(infns_size+nargv_fns, infns, h, outfn? outfn : "-", arg_list, no_pg) < 0)
            ret = 1;
        break;

    case cram:
        sam_close(in);
        if (cram_cat(infns_size+nargv_fns, infns, h, outfn? outfn : "-", &ga, arg_list, no_pg) < 0)
            ret = 1;
        break;

    default:
        sam_close(in);
        fprintf(stderr, "[%s] ERROR: input is not BAM or CRAM\n", __func__);
        return 1;
    }

 end:
    if (infns_size > 0) {
        int i;
        for (i=0; i<infns_size; i++)
            free(infns[i]);
    }

    free(outfn);
    free(infns);
    free(arg_list);
    if (h)
        sam_hdr_destroy(h);

    return ret;
}
