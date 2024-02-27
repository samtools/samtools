/*  bam_cat.c -- efficiently concatenates bam files.

    Copyright (C) 2008-2009, 2011-2013, 2015-2017, 2019, 2021,
                  2023-2024 Genome Research Ltd.
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
#include "htslib/hfile.h"
#include "samtools.h"
#include "sam_opts.h"

/// cat_check_merge_hdr - check compatibility and merge RG hearders merges RGon both CRAM and BAM.
/** @param firstfile - pointer to the 1sr file opened in caller
 *  @param nfn - number of files to be processed, including the firstfile
 *  @param fn - array of file paths to be processed
 *  @param h - sam header pointer which contains explicitly given header
 *  @param vers_maj_p - cram major version set and send out for output creation
 *  @param vers_min_p - cram min version set and send out for output creation
 *  @param out_h - pointer to sam header pointer, outputs the merged header
 * returns array of opened samFile pointers on success and NULL on failure
 * This method has the merged header processing for cram and bam.
 * RG lines are merged for both cram and bam. For cram, version match for each
 * file and order match of RG lines are compared as well.
 * Note: it is a simple merge of RG lines alone.
*/
static samFile** cat_check_merge_hdr(samFile * const firstfile, int nfn, char * const *fn, const sam_hdr_t *h,
                                     int *vers_maj_p, int *vers_min_p, sam_hdr_t **out_h) {
    int i, vers_maj = -1, vers_min = -1;
    sam_hdr_t *new_h = NULL, *old_h = NULL;
    samFile *in = NULL;
    kstring_t ks = KS_INITIALIZE;
    samFile **files = calloc(nfn, sizeof(samFile *));
    if(!files) {
        fprintf(stderr, "[%s] ERROR: failed to allocate space for file handles.\n", __func__);
        return NULL;
    }
    if (!out_h || !firstfile) {
        fprintf(stderr, "[%s] ERROR: header check failed.\n", __func__);
        goto fail;
    }
    if (*out_h) {           //use header if one is already present
        new_h = *out_h;
    }
    else {
        if (h) {            //use the explicit header given
            new_h = sam_hdr_dup(h);
            if (!new_h) {
                fprintf(stderr, "[%s] ERROR: header duplication failed.\n",
                        __func__);
                goto fail;
            }
        }
    }

    for (i = 0; i < nfn; ++i) {
        int ki;
        //1st file is already open and passed, rest open locally
        files[i] = in = i ? sam_open(fn[i], "r") : firstfile;
        if (in == 0) {
            print_error_errno("cat", "fail to open file '%s'", fn[i]);
            goto fail;
        }
        if (firstfile->format.format != in->format.format) {
            print_error("cat", "File %s is of different format!", fn[i]);
            goto fail;
        }
        if (firstfile->format.format == cram) {     //version check for cram
            cram_fd *in_c;
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
        }

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
            old_h = NULL;
            continue;
        }
        //merge RG lines
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

        if (firstfile->format.format == cram && old_count > 1 && sam_hdr_count_lines(new_h, "RG") == old_count) {
            //RG order check for cram
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

        sam_hdr_destroy(old_h); old_h = NULL;
    }

    ks_free(&ks);

    if (vers_maj_p) {
        *vers_maj_p = vers_maj;
    }
    if (vers_min_p) {
        *vers_min_p = vers_min;
    }
    *out_h = new_h;
    return files;

fail:
    ks_free(&ks);
    if (old_h) sam_hdr_destroy(old_h);
    if (new_h) sam_hdr_destroy(new_h);
    *out_h = NULL;
    for (i = 1; i < nfn; ++i) {         //close files other than the firstfile
        if (files[i]) {
            sam_close(files[i]);
        }
    }
    free(files);

    return NULL;
}


/* ----------------------------------------------------------------------
 * CRAM cat
 */

// Reports the number of CRAM containers spanning a specified region if
// specified, or the entire file if not.
// This is the implements the "samtools cat -q [-r reg]" functionality.
//
// Returns 0 on success, <0 on error.
static int cram_query_ncont(int nfn, char * const *fn, char *reg) {
    int i;
    hts_idx_t *idx = NULL;
    sam_hdr_t *hdr = NULL;
    hts_itr_t *iter = NULL;
    samFile *in = NULL;

    for (i = 0; i < nfn; i++) {
        in = sam_open(fn[i], "r");
        if (!in) {
            print_error_errno("cat", "Couldn't open file %s", fn[i]);
            return -1;
        }
        idx = sam_index_load(in, fn[i]);
        if (!idx) {
            print_error("cat", "No index found");
            goto err;
        }

        off_t cstart = 0, cend = 0;
        if (reg) {
            sam_hdr_t *hdr = sam_hdr_read(in);
            hts_itr_t *iter = sam_itr_querys(idx, hdr, reg);
            if (!hdr) {
                print_error("cat", "Unable to read header");
                goto err;
            }
            if (!iter) {
                print_error("cat", "Unable to parse region");
                goto err;
            }

            if (cram_index_extents(in->fp.cram, iter->tid, iter->beg,
                                   iter->end, &cstart, &cend) < 0) {
                print_error("cat", "Failed to query index");
                goto err;
            }

            hts_itr_destroy(iter);
            sam_hdr_destroy(hdr);
        }

        int64_t first, last;
        int64_t nc = cram_num_containers_between(in->fp.cram, cstart, cend,
                                                 &first, &last);
        printf("%s\t%"PRId64"\t%"PRId64"\t%"PRId64"\n",
               fn[i], nc, first, last);
        sam_close(in);
        hts_idx_destroy(idx);
    }

    return 0;

 err:
    if (idx)
        hts_idx_destroy(idx);
    if (hdr)
        sam_hdr_destroy(hdr);
    if (iter)
        hts_itr_destroy(iter);

    sam_close(in);

    return -1;
}

// Container range #:A-B or #:A.
// Returns 0 on success, -1 on failure.
static int cram_handle_cnum_region(cram_fd *fd, hts_idx_t *idx,
                                   char *reg, off_t *cstart, off_t *cend) {
    int cnum_start, cnum_end; // container versions
    int n = sscanf(reg, "#:%d-%d", &cnum_start, &cnum_end);

    if (n == 0) {
        print_error("cat", "ERROR: Malformed region: %s", reg);
        return -1;
    } else if (n == 1) {
        cnum_end = cnum_start;
    }

    int64_t nc = cram_num_containers(fd);
    if (cnum_end >= nc) {
        print_error("cat", "Too many containers.  "
                    "The end range should be < %"PRId64, nc);
        return -1;
    }

    // Container number to offset
    *cstart = cram_container_num2offset(fd, cnum_start);
    *cend   = cram_container_num2offset(fd, cnum_end);
    if (*cstart < 0 || *cend < 0)
        return -1;

    // Seek manually
    return cram_seek(fd, *cstart, SEEK_SET);
}

// Normal range chr:start-end.
// Returns an hts iterator on success, NULL on failure.
static hts_itr_t *cram_handle_region(cram_fd *fd, hts_idx_t *idx, sam_hdr_t *h,
                                     char *reg, off_t *cstart, off_t *cend) {
    hts_itr_t *iter;

    if (!idx) {
        fprintf(stderr, "[%s] ERROR: No index found.\n", __func__);
        return NULL;
    }

    // This does an implicit seek and modifies the cram_fd.
    if (!(iter = sam_itr_querys(idx, h, reg))) {
        print_error("cat", "Unable to parse region %s", reg);
        return NULL;
    }
    if (cram_index_extents(fd, iter->tid, iter->beg,
                           iter->end, cstart, cend) < 0) {
        print_error("cat", "Failed to query index");
        return NULL;
    }

    return iter;
}

// Handle the -p A/B option to subdivide our region or file into portions.
// Updates cstart/cend.
// returns 0 on success (do something),
//         1 on success (but nothing to do),
//        -1 on failure.
static int cram_subdivide_part(cram_fd *fd, hts_idx_t *idx, char *part,
                               off_t *cstart, off_t *cend) {
    int a, b;
    // Part N of M
    if (sscanf(part, "%d/%d", &a, &b) != 2) {
        print_error("cat", "malformed region %s. Should be e.g. '1/10'", part);
        return -1;
    }

    // Inclusive container numbers for range, 0 to NC-1
    // Our part N/M is in container percentages as we can't have
    // partial containers, so convert to that first.
    int64_t cnum1, cnum2;
    if (cstart) {
        cnum1 = cram_container_offset2num(fd, *cstart);
        cnum2 = cram_container_offset2num(fd, *cend);
    } else {
        cnum1 = 0;
        cnum2 = cram_num_containers(fd)-1;
    }

    // Subdivide cnum1/cnum2 container numbers to new range cnum_start/end
    int64_t nc = cnum2 - cnum1 + 1;
    if (b > nc)
        b = nc;

    int cnum_start = (a-1)*(double)nc/b;
    int cnum_end   = a*(double)nc/b - 1;
    if (cnum_start < 0 || cnum_end >= nc)
        return 1;

    // Then convert back to file offsets so we can seek and do htell
    // to detect EOR/EOF.
    *cstart = cram_container_num2offset(fd, cnum_start + cnum1);
    *cend   = cram_container_num2offset(fd, cnum_end   + cnum1);

    return 0;
}

// The main cram_cat interface.
// Returns 0 on success, < 0 on error.
int cram_cat(samFile * const firstfile, int nfn, char * const *fn,
             const sam_hdr_t *h, const char* outcram, sam_global_args *ga,
             char *arg_list, int no_pg, char *reg, char *part, int fast_reg)
{
    samFile *out = NULL;
    cram_fd *out_c;
    int i, vers_maj, vers_min, ret = -1;
    sam_hdr_t *new_h = NULL;
    samFile **files = NULL;
    hts_idx_t *idx = NULL;
    hts_itr_t *iter = NULL;
    sam_hdr_t *old_h = NULL;

    // Check consistent versioning and compatible headers;
    // merges RG lines, opens all files and returns them that multiple
    // non-seekable stream inputs can be handled
    if (!(files = cat_check_merge_hdr(firstfile, nfn, fn, h, &vers_maj,
                                      &vers_min, &new_h)))
        return -1;

    if (!new_h) {
        print_error_errno("cat", "failed to make output header");
        goto closefiles;
    }

    /* Open the file with cram_vers */
    char vers[100];
    snprintf(vers, sizeof(vers), "%d.%d", vers_maj, vers_min);

    // Can override level=1 with e.g. "--output-fmt-option level=9"
    out = sam_open_format(outcram, "wc1", &ga->out);
    if (out == 0) {
        print_error_errno("cat", "fail to open output file '%s'", outcram);
        goto closefiles;
    }
    out_c = out->fp.cram;
    cram_set_option(out_c, CRAM_OPT_VERSION, vers);

    if (!no_pg && sam_hdr_add_pg(new_h, "samtools",
                                 "VN", samtools_version(),
                                 arg_list ? "CL": NULL,
                                 arg_list ? arg_list : NULL,
                                 NULL))
        goto closefiles;

    if (sam_hdr_write(out, new_h) < 0) {
        print_error_errno("cat", "Couldn't write header");
        goto closefiles;
    }
    out_c = out->fp.cram;

    for (i = 0; i < nfn; ++i) {
        samFile *in;
        cram_fd *in_c;
        cram_container *c;
        int new_rg = -1;

        in = files[i];
        if (in == 0) {
            print_error_errno("cat", "fail to open file '%s'", fn[i]);
            goto closefiles;
        }
        in_c = in->fp.cram;

        old_h = sam_hdr_read(in);
        if (!old_h) {
            print_error("cat", "fail to read the header of file '%s'", fn[i]);
            goto closefiles;
        }

        // Compute RG mapping if suitable for changing.
        if (sam_hdr_count_lines(old_h, "RG") == 1) {
            const char *old_name = sam_hdr_line_name(old_h, "RG", 0);
            if (old_name) {
                new_rg = sam_hdr_line_index(new_h, "RG", old_name);
                if (new_rg < 0) {
                    print_error("cat", "fail to find @RG line '%s' in the new header", old_name);
                    goto closefiles;
                }
            } else {
                print_error("cat", "fail to find @RG line in file '%s'", fn[i]);
                goto closefiles;
            }
        } else {
            new_rg = 0;
        }

        // We have multiple region syntax.  Either the standard chr:start-end
        // or a cat-specific #:num-num for explicit container numbers.
        // Both of these seek to a specific file offset and also have the
        // end offset known so we can use htell to detect when we're done.
        // Those offsets are in cstart and cend.
        //
        // However for the e.g. -p 1/10 syntax we need to know the container
        // numbers corresponding to the cstart/cend offsets as we can't
        // start half way through a container when doing fractions.
        // These are cnum1 and cnum2.  (We could short cut this for the
        // #:num-num, but it's simpler to just treat all regions identically.)
        off_t cstart = 0, cend = 0;
        int filter_by_cnum = 0;

        if (reg || part) {
            idx = sam_index_load(in, fn[i]);
            if (!idx) {
                print_error("cat", "failed to load index");
                goto closefiles;
            }
        }

        if (reg) {
            if (strncmp(reg, "#:", 2) == 0) {
                // Region as container numbers
                if (cram_handle_cnum_region(in_c, idx, reg, &cstart,&cend) < 0)
                    goto closefiles;

                filter_by_cnum = 1;
            } else {
                // Normal range chr:start-end
                if (!(iter = cram_handle_region(in_c, idx, old_h, reg,
                                                &cstart, &cend)))
                    goto closefiles;
            }
        }

        // We can also take a range above and subdivide it into parts.
        // Eg -r chr1 -p 1/10 (... to -p 10/10).  Part only implies
        // portions of the entire file.
        if (part) {
            int r = cram_subdivide_part(in_c, idx, part, &cstart, &cend);
            if (r != 0) {
                if (r > 0) // Not an error, just nothing to do
                    ret = 0;
                goto closefiles;
            }
        }

        if (cstart) // reg or part
            if (0 != cram_seek(in_c, cstart, SEEK_SET))
                goto closefiles;


        // Make refid -2 ("*") come after other chromosomes, for easy sort
        int itid = iter
            ? (iter->tid == HTS_IDX_NOCOOR ? INT_MAX : iter->tid)
            : 0;
        int last_ref_id = -99;

        off_t before_hdr = htell(cram_fd_get_fp(in_c));

        // Copy contains and blocks within them
        while ((c = cram_read_container(in_c))) {
            if (cram_container_is_empty(in_c)) {
                cram_block *blk;
                // Container compression header
                if (!(blk = cram_read_block(in_c)))
                    goto closefiles;
                cram_free_block(blk);
                cram_free_container(c);
                continue;
            }

            int filter = 0;

            // If we have just one RG key and new_rg != 0 then
            // we need to edit the compression header. IF WE CAN.
            if (new_rg) {
                if (reg) {
                    print_error("cat", "Cannot specify a region while "
                                "transcoding RG lines");
                    goto closefiles;
                }
                int zero = 0;
                //fprintf(stderr, "Transcode RG %d to %d\n", 0, new_rg);
                cram_transcode_rg(in_c, out_c, c, 1, &zero, &new_rg);
            } else {
                int32_t num_slices;
                int refid;
                hts_pos_t start, span, end;

                if (reg) {
                    cram_container_get_coords(c, &refid, &start, &span);
                    end = start+span;
                    if (before_hdr > cend) {
                        cram_free_container(c);
                        break;
                    }
                }

                // For chr:start-end regions, do we need to filter or skip?
                if (iter && reg) {
                    // Make refid -1 ("*") come after other chromosomes
                    if (refid == -1)
                        refid = INT_MAX;

                    if (refid > itid || start > iter->end) {
                        // Beyond the requested range
                        break;
                    } else if (refid == -2) {
                        // Multi-ref containers.  We only support this if
                        // the RI data series is in a container by itself.
                        filter = 3;
                    } else if (refid < itid || end < iter->beg) {
                        // Skip, in case of mixed size containers.
                        // Eg: Use "=", skip "-".
                        // ==========|======|=
                        //   ----  ==|==  ==|== ----
                        //     ----- | ====   ----
                        filter = 2;
                    } else if (start < iter->beg || end > iter->end) {
                        // Container overlaps region.
                        // Fast mode just copies overlapping containers.
                        // Slow mode does a precise filtering by reading and
                        // filtering each record in turn so it's the same as
                        // a samtools view command.
                        filter = fast_reg ? 0 : 1;
                    }
                    // else we're in an "internal" container, so just copy
                }

                if (filter && last_ref_id == -1 && itid == INT_MAX)
                    // Multi-ref containers consisting solely of ref "*" are
                    // common, but if it's sorted then we know it's "*" from
                    // here on so we don't need to filter despite multi-ref.
                    filter = 0;

                if (filter) {
                    // Filter or skip
                    cram_filter_container(in_c, out_c, c, &last_ref_id);
                } else {
                    // Copy. Consider adding a cram_copy_container API instead.

                    // Container compression header
                    cram_block *blk;
                    if (!(blk = cram_read_block(in_c)))
                        goto closefiles;

                    // Not switching rg so do the usual read/write loop
                    if (cram_write_container(out_c, c) != 0)
                        goto closefiles;

                    // Contatiner compression header
                    if (cram_write_block(out_c, blk) != 0) {
                        cram_free_block(blk);
                        goto closefiles;
                    }

                    // Container num_blocks can be invalid, due to a bug.
                    // Instead we iterate in slice context instead.
                    (void)cram_container_get_landmarks(c, &num_slices);
                    if (cram_copy_slice(in_c, out_c, num_slices) < 0) {
                        cram_free_block(blk);
                        goto closefiles;
                    }
                    cram_free_block(blk);
                }
            }
            cram_free_container(c);

            // Location of next container start
            before_hdr = htell(cram_fd_get_fp(in_c));
            if (filter_by_cnum && before_hdr > cend)
                break;
        }
        sam_hdr_destroy(old_h);
        old_h = NULL;

        if (idx) {
            hts_idx_destroy(idx);
            idx = NULL;
        }

        if (iter) {
            hts_itr_destroy(iter);
            iter = NULL;
        }
    }
    ret = 0;

closefiles:
    if (old_h)
        sam_hdr_destroy(old_h);

    if (idx)
        hts_idx_destroy(idx);

    if (iter)
        hts_itr_destroy(iter);

    if (out)
        sam_close(out);

    if (new_h)
        sam_hdr_destroy(new_h);

    for (i = 1; i < nfn; ++i) {     //skip firstfile and close rest
        if (files[i]) {
            sam_close(files[i]);
        }
    }
    free(files);
    return ret;
}

/* ----------------------------------------------------------------------
 * BAM cat
 */
#define BUF_SIZE 0x10000

#define GZIPID1 31
#define GZIPID2 139

#define BGZF_EMPTY_BLOCK_SIZE 28

int bam_cat(samFile * const firstfile, int nfn, char * const *fn, sam_hdr_t *h, const char* outbam, char *arg_list, int no_pg)
{
    BGZF *fp = NULL, *in = NULL;
    uint8_t *buf = NULL;
    uint8_t ebuf[BGZF_EMPTY_BLOCK_SIZE];
    const int es=BGZF_EMPTY_BLOCK_SIZE;
    int i;
    samFile **files = NULL;
    sam_hdr_t *new_h = NULL;

    /* merges RG lines, opens all files and returns them that multiple non-seekable
    stream inputs can be handled */
    if (!(files = cat_check_merge_hdr(firstfile, nfn, fn, h, NULL, NULL, &new_h)))
        return -1;
    if (!new_h) {
        print_error_errno("cat", "failed to make output header");
        goto fail;
    }
    fp = strcmp(outbam, "-")? bgzf_open(outbam, "w") : bgzf_fdopen(fileno(stdout), "w");
    if (fp == 0) {
        print_error_errno("cat", "fail to open output file '%s'", outbam);
        goto fail;
    }

    if (!no_pg && sam_hdr_add_pg(new_h, "samtools",
                                    "VN", samtools_version(),
                                    arg_list ? "CL": NULL,
                                    arg_list ? arg_list : NULL,
                                    NULL))
        goto fail;

    if (bam_hdr_write(fp, new_h) < 0) {
        print_error_errno("cat", "Couldn't write header");
        goto fail;
    }

    buf = (uint8_t*) malloc(BUF_SIZE);
    if (!buf) {
        fprintf(stderr, "[%s] Couldn't allocate buffer\n", __func__);
        goto fail;
    }
    for(i = 0; i < nfn; ++i){
        int len,j;
        in = files[i]->fp.bgzf;
        if (in == 0) {
            print_error_errno("cat", "fail to open file '%s'", fn[i]);
            goto fail;
        }
        if (in->is_write) goto fail;

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
        in = NULL;
    }
    free(buf);
    if (bgzf_close(fp) < 0) {
        fprintf(stderr, "[%s] Error on closing '%s'.\n", __func__, outbam);
        goto fail;
    }
    for (i = 1; i < nfn; ++i) {     //skip firstfile and close rest
        if (files[i]) {
            sam_close(files[i]);
        }
    }
    free(files);
    sam_hdr_destroy(new_h);
    return 0;

 write_fail:
    fprintf(stderr, "[%s] Error writing to '%s'.\n", __func__, outbam);
 fail:
    if (new_h) {
        sam_hdr_destroy(new_h);
    }
    if (fp) bgzf_close(fp);
    free(buf);

    if (files) {
        for(i = 1; i < nfn; ++i) {  //except the firstfile
            if(files[i]) {
                sam_close(files[i]);
            }
        }
        free(files);
    }
    return -1;
}

int main_cat(int argc, char *argv[])
{
    sam_hdr_t *h = 0;
    char *outfn = 0;
    char **infns = NULL; // files to concatenate
    int infns_size = 0;
    int c, ret = 0, no_pg = 0, usage = 0, query_ncont = 0, fast_mode = 0;
    samFile *in;
    sam_global_args ga;
    char *reg = NULL, *part = NULL;

    static const struct option lopts[] = {
        SAM_OPT_GLOBAL_OPTIONS('-', '-', '-', 0, '-', '-'),
        {"no-PG", no_argument, NULL, 1},
        { NULL, 0, NULL, 0 }
    };

    char *arg_list = NULL;

    sam_global_args_init(&ga);
    while ((c = getopt_long(argc, argv, "h:o:b:r:p:qf", lopts, NULL)) >= 0) {
        switch (c) {
            case 'h': {
                samFile *fph = sam_open(optarg, "r");
                if (fph == 0) {
                    fprintf(stderr, "[%s] ERROR: fail to read the header from '%s'.\n", __func__, optarg);
                    ret = 1;
                    goto end;
                }
                h = sam_hdr_read(fph);
                if (h == NULL) {
                    fprintf(stderr,
                            "[%s] ERROR: failed to read the header from '%s'.\n",
                            __func__, optarg);
                    ret = 1;
                    goto end;
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
            case 'r':
                reg = optarg;
                break;
            case 'p':
                part = optarg;
                break;
            case 'f':
                fast_mode = 1;
                break;
            case 'q':
                query_ncont=1;
                break;

            default:
                if (parse_sam_global_opt(c, optarg, lopts, &ga) == 0) break;
                /* else fall-through */
            case '?': usage=1; break;
        }
    }

    if (!no_pg && !(arg_list = stringify_argv(argc+1, argv-1))) {
        print_error("cat", "failed to create arg_list");
        ret = 1;
        goto end;
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
        fprintf(stderr, "\nCRAM only options for filtering:\n");
        fprintf(stderr, "         -r REG   filter to region REG.\n");
        fprintf(stderr, "                  REG can also be #:cstart-cend for specific container numbers\n");
        fprintf(stderr, "         -p N/M   Specify part N of M (where N is 1 to M inclusive)\n");
        fprintf(stderr, "         -f       Fast mode: don't filter containers to exactly match region\n");
        fprintf(stderr, "         -q       Query the total number of indexed containers\n");
        fprintf(stderr, "\nStandard options:\n");
        sam_global_opt_help(stderr, "---.---.");
        ret = 1;
        goto end;
    }

    in = sam_open(infns[0], "r");
    if (!in) {
        print_error_errno("cat", "failed to open file '%s'", infns[0]);
        ret = 1;
        goto end;
    }

    switch (hts_get_format(in)->format) {
    case bam:
        if (bam_cat(in, infns_size+nargv_fns, infns, h, outfn? outfn : "-", arg_list, no_pg) < 0)
            ret = 1;
        break;

    case cram:
        if (query_ncont) {
            if (cram_query_ncont(infns_size+nargv_fns, infns, reg) < 0)
                ret = 1;
        } else {
            if (cram_cat(in, infns_size+nargv_fns, infns, h,
                         outfn? outfn : "-", &ga, arg_list, no_pg, reg,
                         part, fast_mode) < 0)
            ret = 1;
        }
        break;

    default:
        fprintf(stderr, "[%s] ERROR: input is not BAM or CRAM\n", __func__);
        ret = 1;
    }
    sam_close(in);

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
    sam_global_args_free(&ga);

    return ret;
}
