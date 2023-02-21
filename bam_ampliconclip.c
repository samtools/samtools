/*  bam_ampliconclip.c -- loads amplicon primers from a BED file and cuts reads
                          from the 5' end.

    Copyright (C) 2020-2022 Genome Research Ltd.

    Authors: Andrew Whitwham <aw7@sanger.ac.uk>
             Rob Davies <rmd+git@sanger.ac.uk>

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
DEALINGS IN THE SOFTWARE
*/

#include <config.h>

#include <assert.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>
#include "htslib/thread_pool.h"
#include "sam_opts.h"
#include <htslib/hts.h>
#include "htslib/hfile.h"
#include "htslib/kstring.h"
#include "htslib/sam.h"
#include "samtools.h"
#include "bam_ampliconclip.h"

typedef enum {
    soft_clip,
    hard_clip
} clipping_type;

typedef struct {
    int add_pg;
    int use_strand;
    int write_clipped;
    int mark_fail;
    int both;
    int fail_len;
    int filter_len;
    int unmapped;
    int oa_tag;
    int del_tag;
    int tol;
    char *arg_list;
    char *stats_file;
    char *rejects_file;
} cl_param_t;


static int bed_entry_sort(const void *av, const void *bv) {
    bed_entry_t *a = (bed_entry_t *) av;
    bed_entry_t *b = (bed_entry_t *) bv;
    return a->right < b->right ? -1 : (a->right == b->right ? 0 : 1);
}


int load_bed_file_multi_ref(char *infile, int get_strand, int sort_by_pos, khash_t(bed_list_hash) *bed_lists) {
    hFILE *fp;
    int line_count = 0, ret;
    int64_t left, right;
    kstring_t line = KS_INITIALIZE;
    bed_entry_list_t *list;
    khiter_t bed_itr;

    if ((fp = hopen(infile, "r")) == NULL) {
        print_error_errno("amplicon", "unable to open file %s.", infile);
        return 1;
    }

    char ref[1024];

    while (line.l = 0, kgetline(&line, (kgets_func *)hgets, fp) >= 0) {
        line_count++;
        int hret;
        char strand;

        if (line.l == 0 || *line.s == '#') continue;
        if (strncmp(line.s, "track ", 6) == 0) continue;
        if (strncmp(line.s, "browser ", 8) == 0) continue;

        if (get_strand) {
            if (sscanf(line.s, "%1023s %"SCNd64" %"SCNd64" %*s %*s %c",
                       ref, &left, &right, &strand) != 4) {
                fprintf(stderr, "[amplicon] error: bad bed file format in line %d of %s.\n"
                                "(N.B. ref/chrom name limited to 1023 characters.)\n",
                                    line_count, infile);
                ret = 1;
                goto error;
            }
        } else {
            if (sscanf(line.s, "%1023s %"SCNd64" %"SCNd64,
                       ref, &left, &right) != 3) {
                fprintf(stderr, "[amplicon] error: bad bed file format in line %d of %s\n"
                                "(N.B. ref/chrom name limited to 1023 characters.)\n",
                                    line_count, infile);
                ret = 1;
                goto error;
            }
        }

        bed_itr = kh_get(bed_list_hash, bed_lists, ref);

        if (bed_itr == kh_end(bed_lists)) { // new ref entry
            char *ref_name = strdup(ref); // need a copy for the hash key

            if (!ref_name) {
                fprintf(stderr, "[amplicon] error: unable to allocate memory for ref name.\n");
                ret = 1;
                goto error;
            }

            bed_itr = kh_put(bed_list_hash, bed_lists, ref_name, &hret);

            if (hret > 0) {
                list = &kh_val(bed_lists, bed_itr);

                // initialise the new hash entry
                list->longest = 0;
                list->size = 0;
                list->length = 0;
                list->bp = NULL;
            } else {
                fprintf(stderr, "[amplicon] error: ref hashing failure.\n");
                ret = 1;
                goto error;
            }
        } else { // existing ref
            list = &kh_val(bed_lists, bed_itr);
        }

        if (list->length == list->size) {
           bed_entry_t *tmp;

           list->size += list->size / 2 + 256;

           if ((tmp = realloc(list->bp, list->size * sizeof(bed_entry_t))) == NULL) {
               fprintf(stderr, "[amplicon] error: unable to allocate more memory for bed data.\n");
               ret = 1;
               goto error;
           }

           list->bp = tmp;
        }

        list->bp[list->length].left  = left;
        list->bp[list->length].right = right;

        if (get_strand) {
            if (strand == '+') {
                list->bp[list->length].rev = 0;
            } else if (strand == '-') {
                list->bp[list->length].rev = 1;
            } else {
                fprintf(stderr, "[amplicon] error: bad strand value in line %d, expecting '+' or '-', found '%c'.\n",
                            line_count, strand);
                ret = 1;
                goto error;
            }
        }

        if (right - left > list->longest)
            list->longest = right - left;

        list->length++;
    }

    if (sort_by_pos) {
        for (bed_itr = kh_begin(bed_lists); bed_itr != kh_end(bed_lists); ++bed_itr) {
            if (kh_exist(bed_lists, bed_itr)) {
                list = &kh_val(bed_lists, bed_itr);
                qsort(list->bp, list->length, sizeof(list->bp[0]), bed_entry_sort);
            }
        }
    }

    if (kh_size(bed_lists) > 0) {// any entries
        ret = 0;
    } else {
        ret = 1;
    }

error:
    ks_free(&line);

    if (hclose(fp) != 0) {
        fprintf(stderr, "[amplicon] warning: failed to close %s", infile);
    }

    return ret;
}


void destroy_bed_hash(khash_t(bed_list_hash) *hash) {
    khiter_t itr;

    for (itr = kh_begin(hash); itr != kh_end(hash); ++itr) {
       if (kh_exist(hash, itr)) {
           free(kh_val(hash, itr).bp);
           free((char *)kh_key(hash, itr));
           kh_key(hash, itr) = NULL;
        }
    }

    kh_destroy(bed_list_hash, hash);
}


static int matching_clip_site(bed_entry_list_t *sites, hts_pos_t pos,
                              int is_rev, int use_strand, int64_t longest,
                              cl_param_t *param) {
    int i, size;  // may need this to be variable
    int tol = param->tol;
    int l = 0, mid = sites->length / 2, r = sites->length;
    int pos_tol = is_rev ? (pos > tol ? pos - tol : 0) : pos;

    while (r - l > 1) {
        if (sites->bp[mid].right <= pos_tol) {
            l = mid;
        } else {
            r = mid;
        }
        mid = (l + r) / 2;
    }

    size = 0;

    for (i = l; i < sites->length; i++) {
        hts_pos_t mod_left, mod_right;

        if (use_strand && is_rev != sites->bp[i].rev)
            continue;

        if (is_rev) {
            mod_left = sites->bp[i].left;
            mod_right = sites->bp[i].right + tol;
        } else {
            if (sites->bp[i].left > tol) {
                mod_left = sites->bp[i].left - tol;
            } else {
                mod_left = 0;
            }
            mod_right = sites->bp[i].right;
        }

        if (pos + longest + tol < mod_right)
            break;

        if (pos >= mod_left && pos <= mod_right) {
            if (is_rev) {
                if (size < pos - sites->bp[i].left) {
                    size = pos - sites->bp[i].left;
                }
            } else {
                if (size < sites->bp[i].right - pos) {
                    size = sites->bp[i].right - pos;
                }
            }
        }
    }

    return size;
}


static int bam_trim_left(bam1_t *rec, bam1_t *rec_out, uint32_t bases,
                         clipping_type clipping) {
    uint32_t *orig_cigar = bam_get_cigar(rec);
    uint8_t *orig_seq = bam_get_seq(rec);
    uint8_t *orig_qual = bam_get_qual(rec);
    uint8_t *orig_aux = bam_get_aux(rec);
    uint32_t *new_cigar;
    uint8_t *new_qual;
    size_t orig_l_aux = bam_get_l_aux(rec);
    uint32_t i, j;
    uint32_t ref_remove = bases, qry_removed = 0, hardclip = 0;
    hts_pos_t new_pos = rec->core.pos;
    uint32_t cig_type, cig_op;

    if (rec->l_data + 8 > rec_out->m_data) {
        uint8_t *new_data = realloc(rec_out->data, rec->l_data + 8);
        if (!new_data) {
            fprintf(stderr, "[ampliconclip] error: could not allocate memoy for new bam record\n");
            return 1;
        }
        rec_out->data = new_data;
        rec_out->m_data = rec->l_data + 8;
    }

    // Copy core data & name
    memcpy(&rec_out->core, &rec->core, sizeof(rec->core));
    memcpy(rec_out->data, rec->data, rec->core.l_qname);

    // Modify CIGAR
    new_cigar = bam_get_cigar(rec_out);

    for (i = 0;  i < rec->core.n_cigar; i++) {
        cig_op = bam_cigar_op(orig_cigar[i]);
        cig_type = bam_cigar_type(cig_op);

        if (cig_op == BAM_CHARD_CLIP) {
            hardclip += bam_cigar_oplen(orig_cigar[i]);
        } else {
            if (cig_type & 2) {
                if (bam_cigar_oplen(orig_cigar[i]) <= ref_remove) {
                    ref_remove -= bam_cigar_oplen(orig_cigar[i]);
                } else {
                    break;
                }
                new_pos += bam_cigar_oplen(orig_cigar[i]);
            }
            if (cig_type & 1) {
                qry_removed += bam_cigar_oplen(orig_cigar[i]);
            }
        }
    }

    if (i < rec->core.n_cigar) {
        cig_type = bam_cigar_type(bam_cigar_op(orig_cigar[i]));

        // account for the last operation
        if (cig_type & 2) {
            new_pos += ref_remove;
        }
        if (cig_type & 1) {
            qry_removed += ref_remove;
        }
    } else {
        if (clipping == hard_clip) {

            rec_out->core.l_qseq = 0;
            rec_out->core.n_cigar = 0;

            if (orig_l_aux)
                memcpy(bam_get_aux(rec_out), orig_aux, orig_l_aux);

            rec_out->l_data = bam_get_aux(rec_out) - rec_out->data + orig_l_aux;

            return 0;
        }

        qry_removed = rec->core.l_qseq;
    }

    j = 0;
    if (clipping == hard_clip && hardclip + qry_removed > 0) {
        new_cigar[j++] = bam_cigar_gen(hardclip + qry_removed, BAM_CHARD_CLIP);
    }
    if (clipping == soft_clip) {
        if (hardclip > 0) {
            new_cigar[j++] = bam_cigar_gen(hardclip, BAM_CHARD_CLIP);
        }
        if (qry_removed > 0) {
            new_cigar[j++] = bam_cigar_gen(qry_removed, BAM_CSOFT_CLIP);
        }
    }

    if (i < rec->core.n_cigar
        && bam_cigar_oplen(orig_cigar[i]) > ref_remove) {
        new_cigar[j++] = bam_cigar_gen(bam_cigar_oplen(orig_cigar[i]) - ref_remove, bam_cigar_op(orig_cigar[i]));

        // fill in the rest of the cigar
        i++;

        for (; i < rec->core.n_cigar; i++) {
            new_cigar[j++] = orig_cigar[i];
        }
    }

    rec_out->core.n_cigar = j;

    if (clipping == soft_clip) {
        qry_removed = 0; // Copy all the sequence and confidence values
    }

    new_qual = bam_get_seq(rec_out) + (rec->core.l_qseq - qry_removed + 1) / 2;
    // Copy remaining SEQ
    if ((qry_removed & 1) == 0) {
        memcpy(bam_get_seq(rec_out), orig_seq + (qry_removed / 2),
                (rec->core.l_qseq - qry_removed + 1) / 2); // +1 to account for odd numbers
    } else {
        uint8_t *in = orig_seq + qry_removed / 2;
        uint8_t *out = bam_get_seq(rec_out);
        uint32_t i;
        for (i = qry_removed; i < rec->core.l_qseq - 1; i += 2) {
            *out++ = ((in[0] & 0x0f) << 4) | ((in[1] & 0xf0) >> 4);
            in++;
        }
        if (i < rec->core.l_qseq) {
            *out++ = (in[0] & 0x0f) << 4;
        }
        assert(out == new_qual);
    }

    // Copy remaining QUAL
    memmove(new_qual, orig_qual, rec->core.l_qseq - qry_removed);

    // Set new l_qseq
    rec_out->core.l_qseq -= qry_removed;

    // Move AUX
    if (orig_l_aux)
        memcpy(bam_get_aux(rec_out), orig_aux, orig_l_aux);

    // Set new l_data
    rec_out->l_data = bam_get_aux(rec_out) - rec_out->data + orig_l_aux;

    // put in new pos
    rec_out->core.pos = new_pos;

    return 0;
}


static int bam_trim_right(bam1_t *rec, bam1_t *rec_out, uint32_t bases,
                          clipping_type clipping) {
    uint32_t *orig_cigar = bam_get_cigar(rec);
    uint8_t *orig_seq = bam_get_seq(rec);
    uint8_t *orig_qual = bam_get_qual(rec);
    uint8_t *orig_aux = bam_get_aux(rec);
    uint32_t *new_cigar;
    uint32_t new_n_cigar = 0;
    uint8_t *new_qual;
    size_t orig_l_aux = bam_get_l_aux(rec);
    int32_t i;
    int32_t j;
    uint32_t ref_remove = bases, qry_removed = 0, hardclip = 0;
    uint32_t cig_type, cig_op;

    if (rec->l_data + 8 > rec_out->m_data) {
        uint8_t *new_data = realloc(rec_out->data, rec->l_data + 8);
        if (!new_data) {
            fprintf(stderr, "[ampliconclip] error: could not allocate memoy for new bam record\n");
            return 1;
        }
        rec_out->data = new_data;
        rec_out->m_data = rec->l_data + 8;
    }

    // Copy core data & name
    memcpy(&rec_out->core, &rec->core, sizeof(rec->core));
    memcpy(rec_out->data, rec->data, rec->core.l_qname);

    // Modify CIGAR here
    new_cigar = bam_get_cigar(rec_out);

    for (i = rec->core.n_cigar - 1;  i >= 0; --i) {
        cig_op = bam_cigar_op(orig_cigar[i]);
        cig_type = bam_cigar_type(cig_op);

        if (cig_op == BAM_CHARD_CLIP) {
            hardclip += bam_cigar_oplen(orig_cigar[i]);
        } else {
            if (cig_type & 2) {
                if (bam_cigar_oplen(orig_cigar[i]) <= ref_remove) {
                    ref_remove -= bam_cigar_oplen(orig_cigar[i]);
                } else {
                    break;
                }
            }
            if (cig_type & 1) {
                qry_removed += bam_cigar_oplen(orig_cigar[i]);
            }
        }
    }

    if (i >= 0) {
        cig_type = bam_cigar_type(bam_cigar_op(orig_cigar[i]));
        if (cig_type & 1) {
            qry_removed += ref_remove;
        }
        j = i;
        if (qry_removed > 0) j++;
        if (hardclip > 0 && (clipping == soft_clip || qry_removed == 0)) j++;
    } else {
        if (clipping == hard_clip) {

            rec_out->core.l_qseq = 0;
            rec_out->core.n_cigar = 0;

            if (orig_l_aux)
                memcpy(bam_get_aux(rec_out), orig_aux, orig_l_aux);

            rec_out->l_data = bam_get_aux(rec_out) - rec_out->data + orig_l_aux;

            return 0;
        }

        qry_removed = rec->core.l_qseq;
        j = 0;
        if (hardclip > 0 && clipping == soft_clip) j++;
    }

    if (clipping == hard_clip && hardclip + qry_removed > 0) {
        new_cigar[j] = bam_cigar_gen(hardclip + qry_removed, BAM_CHARD_CLIP);
        new_n_cigar++;
    }
    if (clipping == soft_clip) {
        if (hardclip > 0) {
            new_cigar[j] = bam_cigar_gen(hardclip, BAM_CHARD_CLIP);
            new_n_cigar++;
            if (qry_removed > 0) --j;
        }
        if (qry_removed > 0) {
            new_cigar[j] = bam_cigar_gen(qry_removed, BAM_CSOFT_CLIP);
            new_n_cigar++;
        }
    }

    if (j > 0) {
        new_cigar[--j] = bam_cigar_gen(bam_cigar_oplen(orig_cigar[i]) - ref_remove, bam_cigar_op(orig_cigar[i]));
        new_n_cigar++;
    }

    // fill in the rest of the cigar
    while (j > 0) {
        new_cigar[--j] = orig_cigar[--i];
        new_n_cigar++;
    }

    rec_out->core.n_cigar = new_n_cigar;

    if (clipping == soft_clip)
        qry_removed = 0; // Copy all the sequence and confidence values

    new_qual = bam_get_seq(rec_out) + (rec->core.l_qseq - qry_removed + 1) / 2;
    // Copy remaining SEQ
    memcpy(bam_get_seq(rec_out), orig_seq, (rec->core.l_qseq - qry_removed + 1) / 2);

    // Copy remaining QUAL
    memcpy(new_qual, orig_qual, rec->core.l_qseq - qry_removed);

    // Set new l_qseq
    rec_out->core.l_qseq -= qry_removed;

    // Copy AUX
    if (orig_l_aux)
        memcpy(bam_get_aux(rec_out), orig_aux, orig_l_aux);

    // Set new l_data
    rec_out->l_data = bam_get_aux(rec_out) - rec_out->data + orig_l_aux;

    return 0;
}


static hts_pos_t active_query_len(bam1_t *b) {
    uint32_t *cigar = bam_get_cigar(b);
    uint32_t cig_type, cig_op;
    hts_pos_t len = 0;
    int i;

    for (i = 0; i < b->core.n_cigar; i++) {
        cig_op =  bam_cigar_op(cigar[i]);
        cig_type = bam_cigar_type(cig_op);

        if ((cig_type & 1) && (cig_op != BAM_CSOFT_CLIP)) {
            len += bam_cigar_oplen(cigar[i]);
        }
    }

    return len;
}


static inline void swap_bams(bam1_t **a, bam1_t **b) {
    bam1_t *tmp = *a;
    *a = *b;
    *b = tmp;
}


// Format OA:Z:(RNAME,POS,strand,CIGAR,MAPQ,NM;
static inline int tag_original_data(bam1_t *orig, kstring_t *oa_tag) {
    char strand;
    uint8_t *nm_tag, *old_oa_tag;
    uint32_t *cigar;
    int64_t nm = 0;
    int i, res = 0;

    ks_clear(oa_tag);

    // if there is an existing OA tag the new one gets appended to it
    if ((old_oa_tag = bam_aux_get(orig, "OA"))) {
        res |= ksprintf(oa_tag, "%s", bam_aux2Z(old_oa_tag)) < 0;
    }

    if (orig->core.flag & BAM_FREVERSE)
        strand = '-';
    else
        strand = '+';

    if ((nm_tag = bam_aux_get(orig, "NM"))) {
        nm = bam_aux2i(nm_tag);
    }

    res |= ksprintf(oa_tag, "%s,%"PRIhts_pos",%c,", bam_get_qname(orig), orig->core.pos + 1, strand) < 0;

    for (i = 0, cigar = bam_get_cigar(orig); i < orig->core.n_cigar && res == 0; ++i) {
        res |= kputw(bam_cigar_oplen(cigar[i]), oa_tag) < 0;
        res |= kputc(bam_cigar_opchr(cigar[i]), oa_tag) < 0;
    }

    if (nm_tag) {
        res |= ksprintf(oa_tag, ",%d,%"PRId64";", orig->core.qual, nm) < 0;
    } else {
        res |= ksprintf(oa_tag, "%d,;", orig->core.qual) < 0;
    }

    return res;
}


static int bam_clip(samFile *in, samFile *out, samFile *reject, char *bedfile,
                    clipping_type clipping, cl_param_t *param) {
    int ret = 1, r, file_open = 0;

    bam_hdr_t *header = NULL;
    bam1_t *b = NULL, *b_tmp = NULL;
    long f_count = 0, r_count = 0, n_count = 0, l_count = 0, l_exclude = 0, b_count = 0;
    long filtered = 0, written = 0, failed = 0;
    kstring_t str = KS_INITIALIZE;
    kstring_t oat = KS_INITIALIZE;
    bed_entry_list_t *sites;
    FILE *stats_fp = stderr;
    khash_t(bed_list_hash) *bed_hash = kh_init(bed_list_hash);

    if (load_bed_file_multi_ref(bedfile, param->use_strand, 1, bed_hash)) {
        fprintf(stderr, "[ampliconclip] error: unable to load bed file.\n");
        goto fail;
    }

    if ((header = sam_hdr_read(in)) == NULL) {
        fprintf(stderr, "[ampliconclip] error: could not read header\n");
        goto fail;
    }

    // changing pos can ruin coordinate sort order
    if (sam_hdr_find_tag_hd(header, "SO", &str) == 0 && str.s && strcmp(str.s, "coordinate") == 0) {
        const char *new_order = "unknown";

        if (sam_hdr_update_hd(header, "SO", new_order) == -1) {
            fprintf(stderr, "[ampliconclip] error: unable to change sort order to 'SO:%s'\n", new_order);
            goto fail;
        }
    }

    ks_free(&str);

    if (param->add_pg && sam_hdr_add_pg(header, "samtools", "VN", samtools_version(),
                        param->arg_list ? "CL" : NULL,
                        param->arg_list ? param->arg_list : NULL,
                        NULL) != 0) {
        fprintf(stderr, "[ampliconclip] warning: unable to add @PG line to header.\n");
    }
    if (sam_hdr_write(out, header) < 0) {
        fprintf(stderr, "[ampliconclip] error: could not write header.\n");
        goto fail;
    }

    if (reject) {
       if (sam_hdr_write(reject, header) < 0) {
           fprintf(stderr, "[ampliconclip] error: could not write header to rejects file.\n");
           goto fail;
       }
    }

    b = bam_init1();
    b_tmp = bam_init1();
    if (!b || !b_tmp) {
        fprintf(stderr, "[ampliconclip] error: out of memory when trying to create record.\n");
        goto fail;
    }

    int32_t last_tid = -1;
    int ref_found = 0;

    while ((r = sam_read1(in, header, b)) >= 0) {
        hts_pos_t pos;
        int is_rev;
        int p_size;
        int been_clipped  = 0, filter = 0;
        int exclude = (BAM_FUNMAP | BAM_FQCFAIL);
        khiter_t itr;

        l_count++;

        if (b->core.tid != last_tid) {
            const char *ref_name;

            ref_found = 0;
            last_tid = b->core.tid;

            if ((ref_name = sam_hdr_tid2name(header, b->core.tid)) != NULL) {
                itr = kh_get(bed_list_hash, bed_hash, ref_name);

                if (itr != kh_end(bed_hash)) {
                    sites = &kh_val(bed_hash, itr);
                    ref_found = 1;
                }
            }
        }

        if (!(b->core.flag & exclude) && ref_found) {
            if (param->oa_tag)
                if (tag_original_data(b, &oat))
                    goto fail;

            if (!param->both) {
                if (bam_is_rev(b)) {
                    pos = bam_endpos(b);
                    is_rev = 1;
                } else {
                    pos = b->core.pos;
                    is_rev = 0;
                }

                if ((p_size = matching_clip_site(sites, pos, is_rev, param->use_strand, sites->longest, param))) {
                    if (is_rev) {
                        if (bam_trim_right(b, b_tmp, p_size, clipping) != 0)
                            goto fail;

                        swap_bams(&b, &b_tmp);
                        r_count++;
                    } else {
                        if (bam_trim_left(b, b_tmp, p_size, clipping) != 0)
                            goto fail;

                        swap_bams(&b, &b_tmp);
                        f_count++;
                    }

                    if (param->oa_tag) {
                        if (bam_aux_update_str(b, "OA", oat.l + 1, (const char *)oat.s))
                            goto fail;
                    }

                    if (param->del_tag) {
                        uint8_t *tag;

                        if ((tag = bam_aux_get(b, "NM")))
                            bam_aux_del(b, tag);

                        if ((tag = bam_aux_get(b, "MD")))
                            bam_aux_del(b, tag);
                    }

                    been_clipped = 1;
                } else {
                    if (param->mark_fail) {
                        b->core.flag |= BAM_FQCFAIL;
                    }

                    n_count++;
                }
            } else {
                int left = 0, right = 0;

                // left first
                pos = b->core.pos;
                is_rev = 0;

                if ((p_size = matching_clip_site(sites, pos, is_rev, param->use_strand, sites->longest, param))) {
                    if (bam_trim_left(b, b_tmp, p_size, clipping) != 0)
                        goto fail;

                    swap_bams(&b, &b_tmp);
                    f_count++;
                    left = 1;
                    been_clipped = 1;
                }

                // the right
                pos = bam_endpos(b);
                is_rev = 1;

                if ((p_size = matching_clip_site(sites, pos, is_rev, param->use_strand, sites->longest, param))) {
                    if (bam_trim_right(b, b_tmp, p_size, clipping) != 0)
                        goto fail;

                    swap_bams(&b, &b_tmp);
                    r_count++;
                    right = 1;
                    been_clipped = 1;
                }

                if (left || right) {
                    uint8_t *tag;

                    if (param->oa_tag) {
                        if (bam_aux_update_str(b, "OA", oat.l + 1, (const char *)oat.s))
                            goto fail;
                    }

                    if (param->del_tag) {
                        if ((tag = bam_aux_get(b, "NM")))
                            bam_aux_del(b, tag);

                        if ((tag = bam_aux_get(b, "MD")))
                            bam_aux_del(b, tag);
                    }
                }

                if (left && right) {
                    b_count++;
                } else if (!left && !right) {
                    if (param->mark_fail) {
                        b->core.flag |= BAM_FQCFAIL;
                    }

                    n_count++;
                }
            }

            if (param->fail_len >= 0 || param->filter_len >= 0) {
               hts_pos_t aql = active_query_len(b);

               if (param->fail_len >= 0 && aql <= param->fail_len) {
                   b->core.flag |= BAM_FQCFAIL;
               }

               if (param->filter_len >= 0 && aql <= param->filter_len) {
                   filter = 1;
               }
           }

           if (b->core.flag & BAM_FQCFAIL) {
               failed++;
           }

           if (param->write_clipped && !been_clipped) {
               filter = 1;
           }

        } else {
            l_exclude++;

            if (param->unmapped) {
                filter = 1;
            }
        }

        if (!filter) {
            if (sam_write1(out, header, b) < 0) {
                fprintf(stderr, "[ampliconclip] error: could not write line %ld.\n", l_count);
                goto fail;
            }

            written++;
        } else {
            if (reject) {
                if (sam_write1(reject, header, b) < 0) {
                    fprintf(stderr, "[ampliconclip] error: could not write to reject file %s\n",
                            param->rejects_file);
                    goto fail;
                }
            }

            filtered++;
        }
    }

    if (r < -1) {
        fprintf(stderr, "[ampliconclip] error: failed to read input.\n");
        goto fail;
    }

    if (param->stats_file) {
        if ((stats_fp = fopen(param->stats_file, "w")) == NULL) {
            fprintf(stderr, "[ampliconclip] warning: cannot write stats to %s.\n", param->stats_file);
        } else {
            file_open = 1;
        }
    }

    fprintf(stats_fp, "COMMAND: %s\n"
                    "TOTAL READS: %ld\n"
                    "TOTAL CLIPPED: %ld\n"
                    "FORWARD CLIPPED: %ld\n"
                    "REVERSE CLIPPED: %ld\n"
                    "BOTH CLIPPED: %ld\n"
                    "NOT CLIPPED: %ld\n"
                    "EXCLUDED: %ld\n"
                    "FILTERED: %ld\n"
                    "FAILED: %ld\n"
                    "WRITTEN: %ld\n", param->arg_list, l_count, f_count + r_count,
                                    f_count, r_count, b_count, n_count, l_exclude,
                                    filtered, failed, written);

    if (file_open) {
        fclose(stats_fp);
    }

    ret = 0;

fail:
    destroy_bed_hash(bed_hash);
    ks_free(&oat);
    sam_hdr_destroy(header);
    bam_destroy1(b);
    bam_destroy1(b_tmp);
    return ret;
}


static void usage(void) {
    fprintf(stderr, "Usage: samtools ampliconclip -b BED file <input.bam> -o <output.bam>\n\n");
    fprintf(stderr, "Option: \n");
    fprintf(stderr, " -b  FILE            BED file of regions (eg amplicon primers) to be removed.\n");
    fprintf(stderr, " -o  FILE            output file name (default stdout).\n");
    fprintf(stderr, " -f  FILE            write stats to file name (default stderr)\n");
    fprintf(stderr, " -u                  Output uncompressed data\n");
    fprintf(stderr, " --soft-clip         soft clip amplicon primers from reads (default)\n");
    fprintf(stderr, " --hard-clip         hard clip amplicon primers from reads.\n");
    fprintf(stderr, " --both-ends         clip on both 5' and 3' ends.\n");
    fprintf(stderr, " --strand            use strand data from BED file to match read direction.\n");
    fprintf(stderr, " --clipped           only output clipped reads.\n");
    fprintf(stderr, " --fail              mark unclipped, mapped reads as QCFAIL.\n");
    fprintf(stderr, " --filter-len INT    do not output reads INT size or shorter.\n");
    fprintf(stderr, " --fail-len   INT    mark as QCFAIL reads INT size or shorter.\n");
    fprintf(stderr, " --no-excluded       do not write excluded reads (unmapped or QCFAIL).\n");
    fprintf(stderr, " --rejects-file FILE file to write filtered reads.\n");
    fprintf(stderr, " --original          for clipped entries add an OA tag with original data.\n");
    fprintf(stderr, " --keep-tag          for clipped entries keep the old NM and MD tags.\n");
    fprintf(stderr, " --tolerance         match region within this number of bases, default 5.\n");
    fprintf(stderr, " --no-PG             do not add an @PG line.\n");
    sam_global_opt_help(stderr, "-.O..@-.");
    fprintf(stderr, "\nAbout: Soft clips read alignments where they match BED file defined regions.\n"
                    "Default clipping is only on the 5' end.\n\n");
}


int amplicon_clip_main(int argc, char **argv) {
    int c, ret;
    char wmode[4] = {'w', 'b', 0, 0};
    char *bedfile = NULL, *fnout = "-";
    sam_global_args ga = SAM_GLOBAL_ARGS_INIT;
    htsThreadPool p = {NULL, 0};
    samFile *in = NULL, *out = NULL, *reject = NULL;
    clipping_type clipping = soft_clip;
    cl_param_t param = {1, 0, 0, 0, 0, -1, -1, 0, 0, 1, 5, NULL, NULL, NULL};

    static const struct option lopts[] = {
        SAM_OPT_GLOBAL_OPTIONS('-', 0, 'O', 0, 0, '@'),
        {"no-PG", no_argument, NULL, 1002},
        {"soft-clip", no_argument, NULL, 1003},
        {"hard-clip", no_argument, NULL, 1004},
        {"strand", no_argument, NULL, 1005},
        {"clipped", no_argument, NULL, 1006},
        {"fail", no_argument, NULL, 1007},
        {"both-ends", no_argument, NULL, 1008},
        {"filter-len", required_argument, NULL, 1009},
        {"fail-len", required_argument, NULL, 1010},
        {"no-excluded", no_argument, NULL, 1011},
        {"rejects-file", required_argument, NULL, 1012},
        {"original", no_argument, NULL, 1013},
        {"keep-tag", no_argument, NULL, 1014},
        {"tolerance", required_argument, NULL, 1015},
        {NULL, 0, NULL, 0}
    };

    while ((c = getopt_long(argc, argv, "b:@:o:O:f:u", lopts, NULL)) >= 0) {
        switch (c) {
            case 'b': bedfile = optarg; break;
            case 'o': fnout = optarg; break;
            case 'f': param.stats_file = optarg; break;
            case 'u': wmode[2] = '0'; break;
            case 1002: param.add_pg = 0; break;
            case 1003: clipping = soft_clip; break;
            case 1004: clipping = hard_clip; break;
            case 1005: param.use_strand = 1; break;
            case 1006: param.write_clipped = 1; break;
            case 1007: param.mark_fail = 1; break;
            case 1008: param.both = 1; break;
            case 1009: param.filter_len = atoi(optarg); break;
            case 1010: param.fail_len = atoi(optarg); break;
            case 1011: param.unmapped = 1; break;
            case 1012: param.rejects_file = optarg; break;
            case 1013: param.oa_tag = 1; break;
            case 1014: param.del_tag = 0; break;
            case 1015: param.tol = atoi(optarg); break;
            default:  if (parse_sam_global_opt(c, optarg, lopts, &ga) == 0) break;
                      /* else fall-through */
            case '?': usage(); exit(1);
        }
    }

    if (!bedfile) {
        usage();
        return 1;
    }

    if (optind + 1 > argc) {
        usage();
        return 1;
    }

    if (param.tol < 0) {
        fprintf(stderr, "[ampliconclip] warning: invalid tolerance of %d,"
                        " reseting tolerance to default of 5.\n", param.tol);
        param.tol = 5;
    }

    if ((in = sam_open_format(argv[optind], "rb", &ga.in)) == NULL) {
        print_error_errno("ampliconclip", "cannot open input file");
        return 1;
    }

    sam_open_mode(wmode+1, fnout, NULL);

    if ((out = sam_open_format(fnout, wmode, &ga.out)) == NULL) {
        print_error_errno("ampliconclip", "cannot open output file");
        return 1;
    }

    if (param.rejects_file) {
        sam_open_mode(wmode+1, param.rejects_file, NULL);

        if ((reject = sam_open_format(param.rejects_file, wmode, &ga.out)) == NULL) {
            print_error_errno("ampliconclip", "cannot open rejects file");
            return 1;
        }
    }

    if (ga.nthreads > 0) {
        if (!(p.pool = hts_tpool_init(ga.nthreads))) {
            fprintf(stderr, "[ampliconclip] error: cannot create thread pool.\n");
            return 1;
        }
        hts_set_opt(in,  HTS_OPT_THREAD_POOL, &p);
        hts_set_opt(out, HTS_OPT_THREAD_POOL, &p);

        if (reject) {
           hts_set_opt(reject,  HTS_OPT_THREAD_POOL, &p);
        }
    }

    param.arg_list = stringify_argv(argc + 1, argv - 1);

    ret = bam_clip(in, out, reject, bedfile, clipping, &param);

    // cleanup
    sam_close(in);

    if (sam_close(out) < 0) {
        fprintf(stderr, "[ampliconclip] error: error while closing output file %s.\n", argv[optind+1]);
        ret = 1;
    }

    if (reject) {
        if (sam_close(reject) < 0) {
            fprintf(stderr, "[ampliconclip] error: error while closing reject file %s.\n", param.rejects_file);
            ret = 1;
        }
    }

    if (p.pool) hts_tpool_destroy(p.pool);

    sam_global_args_free(&ga);
    free(param.arg_list);

    return ret;
}

