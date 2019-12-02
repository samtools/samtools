/*  bam_sort.c -- sorting and merging.

    Copyright (C) 2008-2019 Genome Research Ltd.
    Portions copyright (C) 2009-2012 Broad Institute.

    Author: Heng Li <lh3@sanger.ac.uk>
    Author: Martin Pollard <mp15@sanger.ac.uk>

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

#include <stdbool.h>
#include <stdlib.h>
#include <ctype.h>
#include <errno.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <sys/stat.h>
#include <unistd.h>
#include <getopt.h>
#include <assert.h>
#include <pthread.h>
#include "htslib/ksort.h"
#include "htslib/hts_os.h"
#include "htslib/khash.h"
#include "htslib/klist.h"
#include "htslib/kstring.h"
#include "htslib/sam.h"
#include "htslib/hts_endian.h"
#include "sam_opts.h"
#include "samtools.h"


// Struct which contains the a record, and the pointer to the sort tag (if any) or
// a combined ref / position / strand.
// Used to speed up tag and position sorts.
typedef struct bam1_tag {
    bam1_t *bam_record;
    union {
        const uint8_t *tag;
        uint8_t pos_tid[12];
    } u;
} bam1_tag;

/* Minimum memory required in megabytes before sort will attempt to run. This
   is to prevent accidents where failing to use the -m option correctly results
   in the creation of a temporary file for each read in the input file.
   Don't forget to update the man page if you change this. */
const size_t SORT_MIN_MEGS_PER_THREAD = 1;

/* Default per-thread memory for sort. Must be >= SORT_MIN_MEGS_PER_THREAD.
   Don't forget to update the man page if you change this. */
const size_t SORT_DEFAULT_MEGS_PER_THREAD = 768;

#if !defined(__DARWIN_C_LEVEL) || __DARWIN_C_LEVEL < 900000L
#define NEED_MEMSET_PATTERN4
#endif

#ifdef NEED_MEMSET_PATTERN4
void memset_pattern4(void *target, const void *pattern, size_t size) {
    uint32_t* target_iter = target;
    size_t loops = size/4;
    size_t i;
    for (i = 0; i < loops; ++i) {
        memcpy(target_iter, pattern, 4);
        ++target_iter;
    }
    if (size%4 != 0)
        memcpy(target_iter, pattern, size%4);
}
#endif

KHASH_INIT(c2c, char*, char*, 1, kh_str_hash_func, kh_str_hash_equal)
KHASH_INIT(cset, char*, char, 0, kh_str_hash_func, kh_str_hash_equal)
KHASH_MAP_INIT_STR(c2i, int)

#define hdrln_free_char(p)
KLIST_INIT(hdrln, char*, hdrln_free_char)

static int g_is_by_qname = 0;
static int g_is_by_tag = 0;
static char g_sort_tag[2] = {0,0};

static int strnum_cmp(const char *_a, const char *_b)
{
    const unsigned char *a = (const unsigned char*)_a, *b = (const unsigned char*)_b;
    const unsigned char *pa = a, *pb = b;
    while (*pa && *pb) {
        if (isdigit(*pa) && isdigit(*pb)) {
            while (*pa == '0') ++pa;
            while (*pb == '0') ++pb;
            while (isdigit(*pa) && isdigit(*pb) && *pa == *pb) ++pa, ++pb;
            if (isdigit(*pa) && isdigit(*pb)) {
                int i = 0;
                while (isdigit(pa[i]) && isdigit(pb[i])) ++i;
                return isdigit(pa[i])? 1 : isdigit(pb[i])? -1 : (int)*pa - (int)*pb;
            } else if (isdigit(*pa)) return 1;
            else if (isdigit(*pb)) return -1;
            else if (pa - a != pb - b) return pa - a < pb - b? 1 : -1;
        } else {
            if (*pa != *pb) return (int)*pa - (int)*pb;
            ++pa; ++pb;
        }
    }
    return *pa? 1 : *pb? -1 : 0;
}

#define HEAP_EMPTY (UINT64_MAX >> 1)

typedef struct {
    int i;
    uint32_t tid;
    uint64_t pos:63, rev:1, idx;
    bam1_tag entry;
} heap1_t;

static inline int bam1_cmp_by_tag(const bam1_tag a, const bam1_tag b);

// Function to compare reads in the heap and determine which one is < the other
static inline int heap_lt(const heap1_t a, const heap1_t b)
{
    if (!a.entry.bam_record)
        return 1;
    if (!b.entry.bam_record)
        return 0;

    if (g_is_by_tag) {
        int t;
        t = bam1_cmp_by_tag(a.entry, b.entry);
        if (t != 0) return t > 0;
    } else if (g_is_by_qname) {
        int t, fa, fb;
        t = strnum_cmp(bam_get_qname(a.entry.bam_record), bam_get_qname(b.entry.bam_record));
        if (t != 0) return t > 0;
        fa = a.entry.bam_record->core.flag & 0xc0;
        fb = b.entry.bam_record->core.flag & 0xc0;
        if (fa != fb) return fa > fb;
    } else {
        if (a.tid != b.tid) return a.tid > b.tid;
        if (a.pos != b.pos) return a.pos > b.pos;
        if (a.rev != b.rev) return a.rev > b.rev;
    }
    // This compares by position in the input file(s)
    if (a.i != b.i) return a.i > b.i;
    return a.idx > b.idx;
}

KSORT_INIT(heap, heap1_t, heap_lt)

typedef struct merged_header {
    sam_hdr_t    *hdr;
    kstring_t     out_rg;
    kstring_t     out_pg;
    kstring_t     out_co;
    char        **target_name;
    uint32_t     *target_len;
    size_t        n_targets;
    size_t        targets_sz;
    khash_t(c2i) *sq_tids;
    khash_t(cset) *rg_ids;
    khash_t(cset) *pg_ids;
    bool          have_hd;
} merged_header_t;

typedef struct trans_tbl {
    int32_t n_targets;
    int* tid_trans;
    kh_c2c_t* rg_trans;
    kh_c2c_t* pg_trans;
    bool lost_coord_sort;
} trans_tbl_t;

static void trans_tbl_destroy(trans_tbl_t *tbl) {
    khiter_t iter;

    free(tbl->tid_trans);

    /*
     * The values for the tbl->rg_trans and tbl->pg_trans hashes are pointers
     * to keys in the rg_ids and pg_ids sets of the merged_header_t, so
     * they should not be freed here.
     *
     * The keys are unique to each hash entry, so they do have to go.
     */

    for (iter = kh_begin(tbl->rg_trans); iter != kh_end(tbl->rg_trans); ++iter) {
        if (kh_exist(tbl->rg_trans, iter)) {
            free(kh_key(tbl->rg_trans, iter));
        }
    }
    for (iter = kh_begin(tbl->pg_trans); iter != kh_end(tbl->pg_trans); ++iter) {
        if (kh_exist(tbl->pg_trans, iter)) {
            free(kh_key(tbl->pg_trans, iter));
        }
    }

    kh_destroy(c2c,tbl->rg_trans);
    kh_destroy(c2c,tbl->pg_trans);
}

/*
 *  Create a merged_header_t struct.
 */

static merged_header_t * init_merged_header() {
    merged_header_t *merged_hdr;

    merged_hdr = calloc(1, sizeof(*merged_hdr));
    if (merged_hdr == NULL) return NULL;

    merged_hdr->hdr = sam_hdr_init();
    if (!merged_hdr->hdr) goto fail;

    merged_hdr->targets_sz   = 16;
    merged_hdr->target_name = malloc(merged_hdr->targets_sz
                                     * sizeof(*merged_hdr->target_name));
    if (NULL == merged_hdr->target_name) goto fail;

    merged_hdr->target_len = malloc(merged_hdr->targets_sz
                                    * sizeof(*merged_hdr->target_len));
    if (NULL == merged_hdr->target_len) goto fail;

    merged_hdr->sq_tids = kh_init(c2i);
    if (merged_hdr->sq_tids == NULL) goto fail;

    merged_hdr->rg_ids = kh_init(cset);
    if (merged_hdr->rg_ids == NULL) goto fail;

    merged_hdr->pg_ids = kh_init(cset);
    if (merged_hdr->pg_ids == NULL) goto fail;

    return merged_hdr;

 fail:
    perror("[init_merged_header]");
    kh_destroy(cset, merged_hdr->pg_ids);
    kh_destroy(cset, merged_hdr->rg_ids);
    kh_destroy(c2i, merged_hdr->sq_tids);
    free(merged_hdr->target_name);
    free(merged_hdr->target_len);
    sam_hdr_destroy(merged_hdr->hdr);
    free(merged_hdr);
    return NULL;
}

/* Some handy kstring manipulating functions */

// Append char range to kstring
static inline int range_to_ks(const char *src, int from, int to,
                              kstring_t *dest) {
    return kputsn(src + from, to - from, dest) != to - from;
}

// Append a kstring to a kstring
static inline int ks_to_ks(kstring_t *src, kstring_t *dest) {
    return kputsn(ks_str(src), ks_len(src), dest) != ks_len(src);
}

/*
 * Generate a unique ID by appending a random suffix to a given prefix.
 * existing_ids is the set of IDs that are already in use.
 * If always_add_suffix is true, the suffix will always be included.
 * If false, prefix will be returned unchanged if it isn't in existing_ids.
 */

static int gen_unique_id(char *prefix, khash_t(cset) *existing_ids,
                         bool always_add_suffix, kstring_t *dest) {
    khiter_t iter;

    if (!always_add_suffix) {
        // Try prefix on its own first
        iter = kh_get(cset, existing_ids, prefix);
        if (iter == kh_end(existing_ids)) { // prefix isn't used yet
            dest->l = 0;
            if (kputs(prefix, dest) == EOF) return -1;
            return 0;
        }
    }

    do {
        dest->l = 0;
        ksprintf(dest, "%s-%0lX", prefix, lrand48());
        iter = kh_get(cset, existing_ids, ks_str(dest));
    } while (iter != kh_end(existing_ids));

    return 0;
}

/*
 * Add the @HD line to the new header
 * In practice the @HD line will come from the first input header.
 */

static int trans_tbl_add_hd(merged_header_t* merged_hdr,
                            sam_hdr_t *translate) {
    kstring_t hd_line = { 0, 0, NULL };
    int res;

    // TODO: handle case when @HD needs merging.
    if (merged_hdr->have_hd) return 0;

    res = sam_hdr_find_hd(translate, &hd_line);
    if (res < -1) {
        print_error("merge", "failed to get @HD line from header");
        return -1;
    }

    if (res < 0) // Not found
        return 0;

    if (sam_hdr_add_lines(merged_hdr->hdr, hd_line.s, hd_line.l) < 0) {
        print_error("merge", "failed to add @HD line to new header");
        free(hd_line.s);
        return -1;
    }

    free(hd_line.s);
    merged_hdr->have_hd = true;

    return 0;
}

/*
 * Add @SQ records to the translation table.
 *
 * Go through the target list for the input header.  Any new targets found
 * are added to the output header target list.  At the same time, a mapping
 * from the input to output target ids is stored in tbl.
 *
 * If any new targets are found, the header text is scanned to find the
 * corresponding @SQ records.  They are then copied into the
 * merged_hdr->out_text kstring (which will eventually become the
 * output header text).
 *
 * Returns 0 on success, -1 on failure.
 */

static int trans_tbl_add_sq(merged_header_t* merged_hdr, sam_hdr_t *translate,
                            trans_tbl_t* tbl) {
    int32_t i;
    int min_tid = -1, res;
    kstring_t sq_line = { 0, 0, NULL }, sq_sn = { 0, 0, NULL };

    // Fill in the tid part of the translation table, adding new targets
    // to the merged header as we go.

    for (i = 0; i < sam_hdr_nref(translate); ++i) {
        int trans_tid;
        sq_sn.l = 0;
        res = sam_hdr_find_tag_pos(translate, "SQ", i, "SN", &sq_sn);
        if (res < 0) {
            print_error("merge", "failed to get @SQ SN #%d from header", i + 1);
            goto fail;
        }

        trans_tid = sam_hdr_name2tid(merged_hdr->hdr, sq_sn.s);
        if (trans_tid < -1) {
            print_error("merge", "failed to lookup ref");
            goto fail;
        }

        if (trans_tid < 0) {
            // Append missing entries to out_hdr
            sq_line.l = 0;
            res = sam_hdr_find_line_id(translate, "SQ", "SN", sq_sn.s, &sq_line);
            if (res < 0) {
                print_error("merge", "failed to get @SQ SN:%s from header", sq_sn.s);
                goto fail;
            }

            trans_tid = sam_hdr_nref(merged_hdr->hdr);

            res = sam_hdr_add_lines(merged_hdr->hdr, sq_line.s, sq_line.l);
            if (res < 0) {
                print_error("merge", "failed to add @SQ SN:%s to new header", sq_sn.s);
                goto fail;
            }
        }
        tbl->tid_trans[i] = trans_tid;

        if (tbl->tid_trans[i] > min_tid) {
            min_tid = tbl->tid_trans[i];
        } else {
            tbl->lost_coord_sort = true;
        }
    }

    free(sq_line.s);
    free(sq_sn.s);

    return 0;

 fail:
    free(sq_line.s);
    free(sq_sn.s);
    return -1;
}

/*
 * Common code for setting up RG and PG record ID tag translation.
 *
 * is_rg is true for RG translation, false for PG.
 * translate is the input bam header
 * merge is true if tags with the same ID are to be merged.
 * known_ids is the set of IDs already in the output header.
 * id_map is the translation map from input header IDs to output header IDs
 * If override is set, it will be used to replace the existing ID (RG only)
 *
 * known_ids and id_map have entries for the new IDs added to them.
 *
 * Return value is a linked list of header lines with the translated IDs,
 * or NULL if something went wrong (probably out of memory).
 *
 */

static klist_t(hdrln) * trans_rg_pg(bool is_rg, sam_hdr_t *translate,
                                    bool merge, khash_t(cset)* known_ids,
                                    khash_t(c2c)* id_map, char *override) {
    khiter_t iter;
    int num_ids, i;
    const char *rec_type = is_rg ? "RG" : "PG";
    klist_t(hdrln) *hdr_lines;

    hdr_lines = kl_init(hdrln);

    // Search through translate's header
    num_ids = sam_hdr_count_lines(translate, rec_type);
    if (num_ids < 0)
        goto fail;

    for (i = 0; i < num_ids; i++) {
        kstring_t orig_id = { 0, 0, NULL };        // ID in original header
        kstring_t transformed_id = { 0, 0, NULL }; // ID in output header
        char *map_value;    // Value to store in id_map
        bool id_changed;    // Have we changed the ID?
        bool not_found_in_output; // ID isn't in the output header (yet)

        if (sam_hdr_find_tag_pos(translate, rec_type, i, "ID", &orig_id) < 0)
            goto fail;

        // is our matched ID in our output ID set already?
        iter = kh_get(cset, known_ids, ks_str(&orig_id));
        not_found_in_output = (iter == kh_end(known_ids));

        if (override) {
            // Override original ID (RG only)
#ifdef OVERRIDE_DOES_NOT_MERGE
            if (gen_unique_id(override, known_ids, false, &transformed_id))
                goto memfail;
            not_found_in_output = true;  // As ID now unique
#else
            if (kputs(override, &transformed_id) == EOF) goto memfail;
            // Know about override already?
            iter = kh_get(cset, known_ids, ks_str(&transformed_id));
            not_found_in_output = (iter == kh_end(known_ids));
#endif
            id_changed = true;
        } else {
            if ( not_found_in_output || merge) {
                // Not in there or merging so can add it as 1-1 mapping
                if (ks_to_ks(&orig_id, &transformed_id)) goto memfail;
                id_changed = false;
            } else {
                // It's in there so we need to transform it by appending
                // a random number to the id
                if (gen_unique_id(ks_str(&orig_id), known_ids,
                                  true, &transformed_id))
                    goto memfail;
                id_changed = true;
                not_found_in_output = true;  // As ID now unique
            }
        }

        // Does this line need to go into our output header?
        if (not_found_in_output) {
            // Take matched line and replace ID with transformed_id
            kstring_t new_hdr_line = { 0, 0, NULL };
            if (sam_hdr_find_line_id(translate, rec_type,
                                     "ID", ks_str(&orig_id), &new_hdr_line) < 0){
                goto fail;
            }

            if (id_changed) {
                char *idp = strstr(ks_str(&new_hdr_line), "\tID:"), *id_end;
                ptrdiff_t id_offset, id_len;
                if (!idp) {
                    print_error("merge", "failed to find ID in \"%s\"\n",
                                ks_str(&new_hdr_line));
                    goto fail;
                }
                idp += 4;
                for (id_end = idp; *id_end >= '\n'; id_end++) {}

                id_offset = idp - new_hdr_line.s;
                id_len = id_end - idp;

                if (id_len < transformed_id.l) {
                    if (ks_resize(&new_hdr_line, new_hdr_line.l + transformed_id.l - id_len))
                        goto fail;
                }
                if (id_len != transformed_id.l) {
                    memmove(new_hdr_line.s + id_offset + transformed_id.l,
                            new_hdr_line.s + id_offset + id_len,
                            new_hdr_line.l - id_offset - id_len + 1);
                }
                memcpy(new_hdr_line.s + id_offset, transformed_id.s,
                       transformed_id.l);
            }

            // append line to output linked list
            char** ln = kl_pushp(hdrln, hdr_lines);
            *ln = ks_release(&new_hdr_line);  // Give away to linked list

            // Need to add it to known_ids set
            int in_there = 0;
            iter = kh_put(cset, known_ids, ks_str(&transformed_id), &in_there);
            if (in_there < 0) goto memfail;
            assert(in_there > 0);  // Should not already be in the map
            map_value = ks_release(&transformed_id);
        } else {
            // Use existing string in id_map
            assert(kh_exist(known_ids, iter));
            map_value = kh_key(known_ids, iter);
            free(ks_release(&transformed_id));
        }

        // Insert it into our translation map
        int in_there = 0;
        iter = kh_put(c2c, id_map, ks_release(&orig_id), &in_there);
        kh_value(id_map, iter) = map_value;
    }

    // If there are no RG lines in the file and we are overriding add one
    if (is_rg && override && kl_begin(hdr_lines) == NULL) {
        kstring_t new_id = {0, 0, NULL};
        kstring_t line = {0, 0, NULL};
        kstring_t empty = {0, 0, NULL};
        int in_there = 0;
        char** ln;

        // Get the new ID
        if (gen_unique_id(override, known_ids, false, &new_id))
            goto memfail;

        // Make into a header line and add to linked list
        ksprintf(&line, "@RG\tID:%s", ks_str(&new_id));
        ln = kl_pushp(hdrln, hdr_lines);
        *ln = ks_release(&line);

        // Put into known_ids set
        iter = kh_put(cset, known_ids, ks_str(&new_id), &in_there);
        if (in_there < 0) goto memfail;
        assert(in_there > 0);  // Should be a new entry

        // Put into translation map (key is empty string)
        if (kputs("", &empty) == EOF) goto memfail;
        iter = kh_put(c2c, id_map, ks_release(&empty), &in_there);
        if (in_there < 0) goto memfail;
        assert(in_there > 0);  // Should be a new entry
        kh_value(id_map, iter) = ks_release(&new_id);
    }

    return hdr_lines;

 memfail:
    perror(__func__);
 fail:
    if (hdr_lines) kl_destroy(hdrln, hdr_lines);
    return NULL;
}

/*
 * Common code for completing RG and PG record translation.
 *
 * Input is a list of header lines, and the mapping from input to
 * output @PG record IDs.
 *
 * RG and PG records can contain tags that cross-reference to other @PG
 * records.  This fixes the tags to contain the new IDs before adding
 * them to the output header text.
 */

static int finish_rg_pg(bool is_rg, klist_t(hdrln) *hdr_lines,
                        khash_t(c2c)* pg_map, kstring_t *out_text) {
    const char *search = is_rg ? "\tPG:" : "\tPP:";
    khiter_t idx;
    char *line = NULL;

    while ((kl_shift(hdrln, hdr_lines, &line)) == 0) {
        char *id = strstr(line, search); // Look for tag to fix
        int pos1 = 0, pos2 = 0;
        char *new_id = NULL;

        if (id) {
            // Found a tag.  Look up the value in the translation map
            // to see what it should be changed to in the output file.
            char *end, tmp;

            id += 4; // Point to value
            end = strchr(id, '\t');  // Find end of tag
            if (!end) end = id + strlen(id);

            tmp = *end;
            *end = '\0'; // Temporarily get the value on its own.

            // Look-up in translation table
            idx = kh_get(c2c, pg_map, id);
            if (idx == kh_end(pg_map)) {
                // Not found, warn.
                fprintf(stderr, "[W::%s] Tag %s%s not found in @PG records\n",
                        __func__, search + 1, id);
            } else {
                // Remember new id and splice points on original string
                new_id = kh_value(pg_map, idx);
                pos1 = id - line;
                pos2 = end - line;
            }

            *end = tmp; // Restore string
        }

        // Copy line to output:
        // line[0..pos1), new_id (if not NULL), line[pos2..end), '\n'

        if (pos1 && range_to_ks(line, 0, pos1, out_text)) goto memfail;
        if (new_id && kputs(new_id, out_text) == EOF) goto memfail;
        if (kputs(line + pos2, out_text) == EOF) goto memfail;
        if (kputc('\n', out_text) == EOF) goto memfail;
        free(line);   // No longer needed
        line = NULL;
    }

    return 0;

 memfail:
    perror(__func__);
    free(line);  // Prevent leakage as no longer on list
    return -1;
}

/*
 * Build the translation table for an input *am file.  This stores mappings
 * which allow IDs to be converted from those used in the input file
 * to the ones which will be used in the output.  The mappings are for:
 *   Reference sequence IDs (for @SQ records)
 *   @RG record ID tags
 *   @PG record ID tags
 *
 * At the same time, new header text is built up by copying records
 * from the input bam file.  This will eventually become the header for
 * the output file.  When copied, the ID tags for @RG and @PG records
 * are replaced with their values.  The @PG PP: and @RG PG: tags
 * are also modified if necessary.
 *
 * merged_hdr holds state on the output header (which IDs are present, etc.)
 * translate is the input header
 * tbl is the translation table that gets filled in.
 * merge_rg controls merging of @RG records
 * merge_pg controls merging of @PG records
 * If rg_override is not NULL, it will be used to replace the existing @RG ID
 *
 * Returns 0 on success, -1 on failure.
 */

static int trans_tbl_init(merged_header_t* merged_hdr, sam_hdr_t* translate,
                          trans_tbl_t* tbl, bool merge_rg, bool merge_pg,
                          bool copy_co, char* rg_override)
{
    kstring_t lines = { 0, 0, NULL };
    klist_t(hdrln) *rg_list = NULL;
    klist_t(hdrln) *pg_list = NULL;

    tbl->n_targets = sam_hdr_nref(translate);
    tbl->rg_trans = tbl->pg_trans = NULL;
    tbl->tid_trans = (int*)calloc(tbl->n_targets ? tbl->n_targets : 1,
                                  sizeof(int));
    if (tbl->tid_trans == NULL) goto memfail;
    tbl->rg_trans = kh_init(c2c);
    if (tbl->rg_trans == NULL) goto memfail;
    tbl->pg_trans = kh_init(c2c);
    if (tbl->pg_trans == NULL) goto memfail;

    tbl->lost_coord_sort = false;

    // Get the @HD record (if not there already).
    if (trans_tbl_add_hd(merged_hdr, translate)) goto fail;

    // Fill in map and add header lines for @SQ records
    if (trans_tbl_add_sq(merged_hdr, translate, tbl)) goto fail;

    // Get translated header lines and fill in map for @RG records
    rg_list = trans_rg_pg(true, translate, merge_rg, merged_hdr->rg_ids,
                          tbl->rg_trans, rg_override);
    if (!rg_list) goto fail;

    // Get translated header lines and fill in map for @PG records
    pg_list = trans_rg_pg(false, translate, merge_pg, merged_hdr->pg_ids,
                          tbl->pg_trans, NULL);

    // Fix-up PG: tags in the new @RG records and add to output
    if (finish_rg_pg(true, rg_list, tbl->pg_trans, &merged_hdr->out_rg))
        goto fail;

    // Fix-up PP: tags in the new @PG records and add to output
    lines.l = 0;
    if (finish_rg_pg(false, pg_list, tbl->pg_trans, &merged_hdr->out_pg))
        goto fail;

    kl_destroy(hdrln, rg_list); rg_list = NULL;
    kl_destroy(hdrln, pg_list); pg_list = NULL;

    if (copy_co) {
        // Just append @CO headers without translation
        int num_co = sam_hdr_count_lines(translate, "CO"), i;
        if (num_co < 0)
            goto fail;

        for (i = 0; i < num_co; i++) {
            if (sam_hdr_find_line_pos(translate, "CO", i, &lines) < 0)
                goto fail;
            if (ks_to_ks(&lines, &merged_hdr->out_co))
                goto fail;
            if (kputc('\n', &merged_hdr->out_co) < 0)
                goto fail;
        }
    }

    free(lines.s);

    return 0;

 memfail:
    perror(__func__);
 fail:
    trans_tbl_destroy(tbl);
    if (rg_list) kl_destroy(hdrln, rg_list);
    if (pg_list) kl_destroy(hdrln, pg_list);
    free(lines.s);
    return -1;
}

static int finish_merged_header(merged_header_t *merged_hdr) {
    if (sam_hdr_add_lines(merged_hdr->hdr, ks_c_str(&merged_hdr->out_rg),
                          ks_len(&merged_hdr->out_rg)) < 0)
        return -1;
    if (sam_hdr_add_lines(merged_hdr->hdr, ks_c_str(&merged_hdr->out_pg),
                          ks_len(&merged_hdr->out_pg)) < 0)
        return -1;
    if (sam_hdr_add_lines(merged_hdr->hdr, ks_c_str(&merged_hdr->out_co),
                          ks_len(&merged_hdr->out_co)) < 0)
        return -1;

    return 0;
}

/*
 * Free a merged_header_t struct and all associated data.
 *
 * Note that the keys to the rg_ids and pg_ids sets are also used as
 * values in the translation tables.  This function should therefore not
 * be called until the translation tables are no longer needed.
 */

static void free_merged_header(merged_header_t *merged_hdr) {
    size_t i;
    khiter_t iter;
    if (!merged_hdr) return;
    free(ks_release(&merged_hdr->out_rg));
    free(ks_release(&merged_hdr->out_pg));
    free(ks_release(&merged_hdr->out_co));
    if (merged_hdr->target_name) {
        for (i = 0; i < merged_hdr->n_targets; i++) {
            free(merged_hdr->target_name[i]);
        }
        free(merged_hdr->target_name);
    }
    free(merged_hdr->target_len);
    kh_destroy(c2i, merged_hdr->sq_tids);

    if (merged_hdr->rg_ids) {
        for (iter = kh_begin(merged_hdr->rg_ids);
             iter != kh_end(merged_hdr->rg_ids); ++iter) {
            if (kh_exist(merged_hdr->rg_ids, iter))
                free(kh_key(merged_hdr->rg_ids, iter));
        }
        kh_destroy(cset, merged_hdr->rg_ids);
    }

    if (merged_hdr->pg_ids) {
        for (iter = kh_begin(merged_hdr->pg_ids);
             iter != kh_end(merged_hdr->pg_ids); ++iter) {
            if (kh_exist(merged_hdr->pg_ids, iter))
                free(kh_key(merged_hdr->pg_ids, iter));
        }
        kh_destroy(cset, merged_hdr->pg_ids);
    }

    free(merged_hdr);
}

static void bam_translate(bam1_t* b, trans_tbl_t* tbl)
{
    // Update target id if not unmapped tid
    if ( b->core.tid >= 0 ) { b->core.tid = tbl->tid_trans[b->core.tid]; }
    if ( b->core.mtid >= 0 ) { b->core.mtid = tbl->tid_trans[b->core.mtid]; }

    // If we have a RG update it
    uint8_t *rg = bam_aux_get(b, "RG");
    if (rg) {
        char* decoded_rg = bam_aux2Z(rg);
        khiter_t k = kh_get(c2c, tbl->rg_trans, decoded_rg);
        if (k != kh_end(tbl->rg_trans)) {
            char* translate_rg = kh_value(tbl->rg_trans,k);
            bam_aux_del(b, rg);
            if (translate_rg) {
                bam_aux_append(b, "RG", 'Z', strlen(translate_rg) + 1,
                               (uint8_t*)translate_rg);
            }
        } else {
            char *tmp = strdup(decoded_rg);
            fprintf(stderr,
                    "[bam_translate] RG tag \"%s\" on read \"%s\" encountered "
                    "with no corresponding entry in header, tag lost. "
                    "Unknown tags are only reported once per input file for "
                    "each tag ID.\n",
                    decoded_rg, bam_get_qname(b));
            bam_aux_del(b, rg);
            // Prevent future whinges
            if (tmp) {
                int in_there = 0;
                k = kh_put(c2c, tbl->rg_trans, tmp, &in_there);
                if (in_there > 0) kh_value(tbl->rg_trans, k) = NULL;
            }
        }
    }

    // If we have a PG update it
    uint8_t *pg = bam_aux_get(b, "PG");
    if (pg) {
        char* decoded_pg = bam_aux2Z(pg);
        khiter_t k = kh_get(c2c, tbl->pg_trans, decoded_pg);
        if (k != kh_end(tbl->pg_trans)) {
            char* translate_pg = kh_value(tbl->pg_trans,k);
            bam_aux_del(b, pg);
            if (translate_pg) {
                bam_aux_append(b, "PG", 'Z', strlen(translate_pg) + 1,
                               (uint8_t*)translate_pg);
            }
        } else {
            char *tmp = strdup(decoded_pg);
            fprintf(stderr,
                    "[bam_translate] PG tag \"%s\" on read \"%s\" encountered "
                    "with no corresponding entry in header, tag lost. "
                    "Unknown tags are only reported once per input file for "
                    "each tag ID.\n",
                    decoded_pg, bam_get_qname(b));
            bam_aux_del(b, pg);
            // Prevent future whinges
            if (tmp) {
                int in_there = 0;
                k = kh_put(c2c, tbl->pg_trans, tmp, &in_there);
                if (in_there > 0) kh_value(tbl->pg_trans, k) = NULL;
            }
        }
    }
}

int* rtrans_build(int n, int n_targets, trans_tbl_t* translation_tbl)
{
    // Create reverse translation table for tids
    int* rtrans = (int*)malloc(sizeof(int32_t)*n*n_targets);
    const int32_t NOTID = INT32_MIN;
    if (!rtrans) return NULL;
    memset_pattern4((void*)rtrans, &NOTID, sizeof(int32_t)*n*n_targets);
    int i;
    for (i = 0; i < n; ++i) {
        int j;
        for (j = 0; j < (translation_tbl+i)->n_targets; ++j) {
            if ((translation_tbl+i)->tid_trans[j] != -1) {
                rtrans[i*n_targets + (translation_tbl+i)->tid_trans[j]] = j;
            }
        }
    }

    return rtrans;
}

#define MERGE_RG          1 // Attach RG tag based on filename
#define MERGE_UNCOMP      2 // Generate uncompressed BAM
#define MERGE_LEVEL1      4 // Compress the BAM at level 1 (fast) mode
#define MERGE_FORCE       8 // Overwrite output BAM if it exists
#define MERGE_COMBINE_RG 16 // Combine RG tags frather than redefining them
#define MERGE_COMBINE_PG 32 // Combine PG tags frather than redefining them
#define MERGE_FIRST_CO   64 // Use only first file's @CO headers (sort cmd only)

/*
 * How merging is handled
 *
 * If a hheader is defined use we will use that as our output header
 * otherwise we use the first header from the first input file.
 *
 * Now go through each file and create a translation table for that file for:
 * -RG
 * -tid
 * -PG tags
 *
 * Then whenever we read a record from a bam we translate that read before
 * stashing it in the hash.
 *
 * In the actual merge, a read is read from each input file, translated and
 * stashed in the hash. This assumes that all input files are sorted in the
 * same way.  Next we just extract the next position ordered read from the
 * hash, and replace it if there are still reads left in it's source input
 * file. Finally we write our chosen read it to the output file.
 */

/*!
  @abstract    Merge multiple sorted BAM.
  @param  by_qname    whether to sort by query name
  @param  sort_tag    if non-null, sort by the given tag
  @param  out         output BAM file name
  @param  mode        sam_open() mode to be used to create the final output file
                      (overrides level settings from UNCOMP and LEVEL1 flags)
  @param  headers     name of SAM file from which to copy '@' header lines,
                      or NULL to copy them from the first file to be merged
  @param  n           number of files to be merged
  @param  fn          names of files to be merged
  @param  flag        flags that control how the merge is undertaken
  @param  reg         region to merge
  @param  n_threads   number of threads to use (passed to htslib)
  @param  cmd         command name (used in print_error() etc)
  @param  in_fmt      format options for input files
  @param  out_fmt     output file format and options
  @param  write_index create the index, together with the output file
  @param  arg_list    command string for PG line
  @param  no_pg       if 1, do not add a new PG line
  @discussion Padding information may NOT correctly maintained. This
  function is NOT thread safe.
 */
int bam_merge_core2(int by_qname, char* sort_tag, const char *out, const char *mode,
                    const char *headers, int n, char * const *fn, char * const *fn_idx,
                    int flag, const char *reg, int n_threads, const char *cmd,
                    const htsFormat *in_fmt, const htsFormat *out_fmt, int write_index,
                    char *arg_list, int no_pg)
{
    samFile *fpout, **fp = NULL;
    heap1_t *heap = NULL;
    sam_hdr_t *hout = NULL;
    sam_hdr_t *hin  = NULL;
    int i, j, *RG_len = NULL;
    uint64_t idx = 0;
    char **RG = NULL;
    hts_itr_t **iter = NULL;
    sam_hdr_t **hdr = NULL;
    trans_tbl_t *translation_tbl = NULL;
    int *rtrans = NULL;
    char *out_idx_fn = NULL;
    merged_header_t *merged_hdr = init_merged_header();
    if (!merged_hdr) return -1;

    // Is there a specified pre-prepared header to use for output?
    if (headers) {
        samFile* fpheaders = sam_open(headers, "r");
        if (fpheaders == NULL) {
            print_error_errno(cmd, "cannot open \"%s\"", headers);
            return -1;
        }
        hin = sam_hdr_read(fpheaders);
        sam_close(fpheaders);
        if (hin == NULL) {
            print_error(cmd, "couldn't read headers from \"%s\"", headers);
            goto mem_fail;
        }
    }

    g_is_by_qname = by_qname;
    if (sort_tag) {
        g_is_by_tag = 1;
        g_sort_tag[0] = sort_tag[0];
        g_sort_tag[1] = sort_tag[0] ? sort_tag[1] : '\0';
    }

    fp = (samFile**)calloc(n, sizeof(samFile*));
    if (!fp) goto mem_fail;
    heap = (heap1_t*)calloc(n, sizeof(heap1_t));
    if (!heap) goto mem_fail;
    iter = (hts_itr_t**)calloc(n, sizeof(hts_itr_t*));
    if (!iter) goto mem_fail;
    hdr = (sam_hdr_t**)calloc(n, sizeof(sam_hdr_t*));
    if (!hdr) goto mem_fail;
    translation_tbl = (trans_tbl_t*)calloc(n, sizeof(trans_tbl_t));
    if (!translation_tbl) goto mem_fail;
    RG = (char**)calloc(n, sizeof(char*));
    if (!RG) goto mem_fail;

    // prepare RG tag from file names
    if (flag & MERGE_RG) {
        RG_len = (int*)calloc(n, sizeof(int));
        if (!RG_len) goto mem_fail;
        for (i = 0; i != n; ++i) {
            int l = strlen(fn[i]);
            const char *s = fn[i];
            if (l > 4 && (strcmp(s + l - 4, ".bam") == 0 || strcmp(s + l - 4, ".sam") == 0)) l -= 4;
            if (l > 5 && strcmp(s + l - 5, ".cram") == 0) l -= 5;
            for (j = l - 1; j >= 0; --j) if (s[j] == '/') break;
            ++j; l -= j;
            RG[i] = (char*)calloc(l + 1, 1);
            if (!RG[i]) goto mem_fail;
            RG_len[i] = l;
            strncpy(RG[i], s + j, l);
        }
    }

    if (hin) {
        // Popluate merged_hdr from the pre-prepared header
        trans_tbl_t dummy;
        int res;
        res = trans_tbl_init(merged_hdr, hin, &dummy, flag & MERGE_COMBINE_RG,
                             flag & MERGE_COMBINE_PG, true, NULL);
        trans_tbl_destroy(&dummy);
        if (res) return -1; // FIXME: memory leak
    }

    // open and read the header from each file
    for (i = 0; i < n; ++i) {
        sam_hdr_t *hin;
        fp[i] = sam_open_format(fn[i], "r", in_fmt);
        if (fp[i] == NULL) {
            print_error_errno(cmd, "fail to open \"%s\"", fn[i]);
            goto fail;
        }
        hin = sam_hdr_read(fp[i]);
        if (hin == NULL) {
            print_error(cmd, "failed to read header from \"%s\"", fn[i]);
            goto fail;
        }

        if (trans_tbl_init(merged_hdr, hin, translation_tbl+i,
                           flag & MERGE_COMBINE_RG, flag & MERGE_COMBINE_PG,
                           (flag & MERGE_FIRST_CO)? (i == 0) : true,
                           RG[i]))
            return -1; // FIXME: memory leak

        // TODO sam_itr_next() doesn't yet work for SAM files,
        // so for those keep the headers around for use with sam_read1()
        if (hts_get_format(fp[i])->format == sam) hdr[i] = hin;
        else { sam_hdr_destroy(hin); hdr[i] = NULL; }

        if ((translation_tbl+i)->lost_coord_sort && !by_qname) {
            fprintf(stderr, "[bam_merge_core] Order of targets in file %s caused coordinate sort to be lost\n", fn[i]);
        }

        // Potential future improvement is to share headers between CRAM files for
        // samtools sort (where all headers are identical.
        // Eg:
        //
        // if (i > 1) {
        //     sam_hdr_free(cram_fd_get_header(fp[i]->fp.cram));
        //     cram_fd_set_header(fp[i]->fp.cram, cram_fd_get_header(fp[0]->fp.cram));
        //     sam_hdr_incr_ref(cram_fd_get_header(fp[0]->fp.cram));
        // }
    }

    // Did we get an @HD line?
    if (!merged_hdr->have_hd) {
        fprintf(stderr, "[W::%s] No @HD tag found.\n", __func__);
        /* FIXME:  Should we add an @HD line here, and if so what should
           we put in it? Ideally we want a way of getting htslib to tell
           us the SAM version number to assume given no @HD line.  Is
           it also safe to assume that the output is coordinate sorted?
           SO: is optional so we don't have to have it.*/
        /* ksprintf(&merged_hdr->out_hd, "@HD\tVN:1.5\tSO:coordinate\n"); */
    }

    // Transform the header into standard form
    if (finish_merged_header(merged_hdr) < 0)
        goto fail;

    hout = merged_hdr->hdr;
    if (!hout) return -1;  // FIXME: memory leak

    // If we're only merging a specified region move our iters to start at that point
    if (reg) {
        int tid;
        hts_pos_t beg, end;

        rtrans = rtrans_build(n, sam_hdr_nref(hout), translation_tbl);
        if (!rtrans) goto mem_fail;

        if (!sam_parse_region(hout, reg, &tid, &beg, &end, 0)) {
            fprintf(stderr, "[%s] Badly formatted region or unknown reference name: \"%s\"\n", __func__, reg);
            goto fail;
        }
        for (i = 0; i < n; ++i) {
            hts_idx_t *idx = NULL;
            // If index filename has not been specfied, look in BAM folder
            if (fn_idx != NULL) {
                idx = sam_index_load2(fp[i], fn[i], fn_idx[i]);
            } else {
                idx = sam_index_load(fp[i], fn[i]);
            }
            // (rtrans[i*n+tid]) Look up what hout tid translates to in input tid space
            int mapped_tid = rtrans[i*sam_hdr_nref(hout)+tid];
            if (idx == NULL) {
                fprintf(stderr, "[%s] failed to load index for %s.  Random alignment retrieval only works for indexed BAM or CRAM files.\n",
                        __func__, fn[i]);
                goto fail;
            }
            if (mapped_tid != INT32_MIN) {
                iter[i] = sam_itr_queryi(idx, mapped_tid, beg, end);
            } else {
                iter[i] = sam_itr_queryi(idx, HTS_IDX_NONE, 0, 0);
            }
            hts_idx_destroy(idx);
            if (iter[i] == NULL) {
                if (mapped_tid != INT32_MIN) {
                    fprintf(stderr,
                            "[%s] failed to get iterator over "
                            "{%s, %d, %"PRIhts_pos", %"PRIhts_pos"}\n",
                            __func__, fn[i], mapped_tid, beg, end);
                } else {
                    fprintf(stderr,
                            "[%s] failed to get iterator over "
                            "{%s, HTS_IDX_NONE, 0, 0}\n",
                            __func__, fn[i]);
                }
                goto fail;
            }
        }
        free(rtrans);
        rtrans = NULL;
    } else {
        for (i = 0; i < n; ++i) {
            if (hdr[i] == NULL) {
                iter[i] = sam_itr_queryi(NULL, HTS_IDX_REST, 0, 0);
                if (iter[i] == NULL) {
                    fprintf(stderr, "[%s] failed to get iterator\n", __func__);
                    goto fail;
                }
            }
            else iter[i] = NULL;
        }
    }

    // Load the first read from each file into the heap
    for (i = 0; i < n; ++i) {
        heap1_t *h = heap + i;
        int res;
        h->i = i;
        h->entry.bam_record = bam_init1();
        h->entry.u.tag = NULL;
        if (!h->entry.bam_record) goto mem_fail;
        res = iter[i] ? sam_itr_next(fp[i], iter[i], h->entry.bam_record) : sam_read1(fp[i], hdr[i], h->entry.bam_record);
        if (res >= 0) {
            bam_translate(h->entry.bam_record, translation_tbl + i);
            h->tid = h->entry.bam_record->core.tid;
            h->pos = (uint64_t)(h->entry.bam_record->core.pos + 1);
            h->rev = bam_is_rev(h->entry.bam_record);
            h->idx = idx++;
            if (g_is_by_tag) {
                h->entry.u.tag = bam_aux_get(h->entry.bam_record, g_sort_tag);
            } else {
                h->entry.u.tag = NULL;
            }
        }
        else if (res == -1 && (!iter[i] || iter[i]->finished)) {
            h->pos = HEAP_EMPTY;
            bam_destroy1(h->entry.bam_record);
            h->entry.bam_record = NULL;
            h->entry.u.tag = NULL;
        } else {
            print_error(cmd, "failed to read first record from \"%s\"", fn[i]);
            goto fail;
        }
    }

    // Open output file and write header
    if ((fpout = sam_open_format(out, mode, out_fmt)) == 0) {
        print_error_errno(cmd, "failed to create \"%s\"", out);
        return -1;
    }
    if (!no_pg && sam_hdr_add_pg(hout, "samtools",
                                 "VN", samtools_version(),
                                 arg_list ? "CL": NULL,
                                 arg_list ? arg_list : NULL,
                                 NULL)) {
        print_error(cmd, "failed to add PG line to the header of \"%s\"", out);
        sam_close(fpout);
        return -1;
    }
    if (sam_hdr_write(fpout, hout) != 0) {
        print_error_errno(cmd, "failed to write header to \"%s\"", out);
        sam_close(fpout);
        return -1;
    }
    if (write_index) {
        if (!(out_idx_fn = auto_index(fpout, out, hout))){
            sam_close(fpout);
            return -1;
        }
    }
    if (!(flag & MERGE_UNCOMP)) hts_set_threads(fpout, n_threads);

    // Begin the actual merge
    ks_heapmake(heap, n, heap);
    while (heap->pos != HEAP_EMPTY) {
        bam1_t *b = heap->entry.bam_record;
        if (flag & MERGE_RG) {
            uint8_t *rg = bam_aux_get(b, "RG");
            if (rg) bam_aux_del(b, rg);
            bam_aux_append(b, "RG", 'Z', RG_len[heap->i] + 1, (uint8_t*)RG[heap->i]);
        }
        if (sam_write1(fpout, hout, b) < 0) {
            print_error_errno(cmd, "failed writing to \"%s\"", out);
            sam_close(fpout);
            free(out_idx_fn);
            return -1;
        }
        if ((j = (iter[heap->i]? sam_itr_next(fp[heap->i], iter[heap->i], b) : sam_read1(fp[heap->i], hdr[heap->i], b))) >= 0) {
            bam_translate(b, translation_tbl + heap->i);
            heap->tid = b->core.tid;
            heap->pos = (uint64_t)(b->core.pos + 1);
            heap->rev = bam_is_rev(b);
            heap->idx = idx++;
            if (g_is_by_tag) {
                heap->entry.u.tag = bam_aux_get(heap->entry.bam_record, g_sort_tag);
            } else {
                heap->entry.u.tag = NULL;
            }
        } else if (j == -1 && (!iter[heap->i] || iter[heap->i]->finished)) {
            heap->pos = HEAP_EMPTY;
            bam_destroy1(heap->entry.bam_record);
            heap->entry.bam_record = NULL;
            heap->entry.u.tag = NULL;
        } else {
            print_error(cmd, "\"%s\" is truncated", fn[heap->i]);
            goto fail;
        }
        ks_heapadjust(heap, 0, n, heap);
    }

    if (write_index) {
        if (sam_idx_save(fpout) < 0) {
            print_error_errno("merge", "writing index failed");
            goto fail;
        }
    }
    free(out_idx_fn);

    // Clean up and close
    if (flag & MERGE_RG) {
        for (i = 0; i != n; ++i) free(RG[i]);
        free(RG_len);
    }
    for (i = 0; i < n; ++i) {
        trans_tbl_destroy(translation_tbl + i);
        hts_itr_destroy(iter[i]);
        sam_hdr_destroy(hdr[i]);
        sam_close(fp[i]);
    }
    sam_hdr_destroy(hin);
    sam_hdr_destroy(hout);
    free_merged_header(merged_hdr);
    free(RG); free(translation_tbl); free(fp); free(heap); free(iter); free(hdr);
    if (sam_close(fpout) < 0) {
        print_error(cmd, "error closing output file");
        return -1;
    }
    return 0;

 mem_fail:
    print_error(cmd, "Out of memory");

 fail:
    if (flag & MERGE_RG) {
        if (RG) {
            for (i = 0; i != n; ++i) free(RG[i]);
        }
        free(RG_len);
    }
    for (i = 0; i < n; ++i) {
        if (translation_tbl && translation_tbl[i].tid_trans) trans_tbl_destroy(translation_tbl + i);
        if (iter && iter[i]) hts_itr_destroy(iter[i]);
        if (hdr && hdr[i]) sam_hdr_destroy(hdr[i]);
        if (fp && fp[i]) sam_close(fp[i]);
        if (heap && heap[i].entry.bam_record) bam_destroy1(heap[i].entry.bam_record);
    }
    if (hout) sam_hdr_destroy(hout);
    free(RG);
    free(translation_tbl);
    free(hdr);
    free(iter);
    free(heap);
    free(fp);
    free(rtrans);
    free(out_idx_fn);
    return -1;
}

// Unused here but may be used by legacy samtools-using third-party code
int bam_merge_core(int by_qname, const char *out, const char *headers, int n, char * const *fn, int flag, const char *reg)
{
    char mode[12];
    strcpy(mode, "wb");
    if (flag & MERGE_UNCOMP) strcat(mode, "0");
    else if (flag & MERGE_LEVEL1) strcat(mode, "1");
    return bam_merge_core2(by_qname, NULL, out, mode, headers, n, fn, NULL, flag, reg, 0, "merge", NULL, NULL, 0, NULL, 1);
}

static void merge_usage(FILE *to)
{
    fprintf(to,
"Usage: samtools merge [-nurlf] [-h inh.sam] [-b <bamlist.fofn>] <out.bam> <in1.bam> [<in2.bam> ... <inN.bam>]\n"
"\n"
"Options:\n"
"  -n         Input files are sorted by read name\n"
"  -t TAG     Input files are sorted by TAG value\n"
"  -r         Attach RG tag (inferred from file names)\n"
"  -u         Uncompressed BAM output\n"
"  -f         Overwrite the output BAM if exist\n"
"  -1         Compress level 1\n"
"  -l INT     Compression level, from 0 to 9 [-1]\n"
"  -R STR     Merge file in the specified region STR [all]\n"
"  -h FILE    Copy the header in FILE to <out.bam> [in1.bam]\n"
"  -c         Combine @RG headers with colliding IDs [alter IDs to be distinct]\n"
"  -p         Combine @PG headers with colliding IDs [alter IDs to be distinct]\n"
"  -s VALUE   Override random seed\n"
"  -b FILE    List of input BAM filenames, one per line [null]\n"
"  -X         Use customized index files\n"
"  --no-PG    do not add a PG line\n");
    sam_global_opt_help(to, "-.O..@..");
}

int bam_merge(int argc, char *argv[])
{
    int c, is_by_qname = 0, flag = 0, ret = 0, level = -1, has_index_file = 0;
    char *fn_headers = NULL, *reg = NULL, mode[12];
    char *sort_tag = NULL, *arg_list = NULL;
    long random_seed = (long)time(NULL);
    char** fn = NULL;
    char** fn_idx = NULL;
    int fn_size = 0, no_pg = 0;

    sam_global_args ga = SAM_GLOBAL_ARGS_INIT;
    static const struct option lopts[] = {
        SAM_OPT_GLOBAL_OPTIONS('-', 0, 'O', 0, 0, '@'),
        { "threads", required_argument, NULL, '@' },
        {"no-PG", no_argument, NULL, 1},
        { NULL, 0, NULL, 0 }
    };

    if (argc == 1) {
        merge_usage(stdout);
        return 0;
    }

    while ((c = getopt_long(argc, argv, "h:nru1R:f@:l:cps:b:O:t:X", lopts, NULL)) >= 0) {
        switch (c) {
        case 'r': flag |= MERGE_RG; break;
        case 'f': flag |= MERGE_FORCE; break;
        case 'h': fn_headers = optarg; break;
        case 'n': is_by_qname = 1; break;
        case 't': sort_tag = optarg; break;
        case '1': flag |= MERGE_LEVEL1; level = 1; break;
        case 'u': flag |= MERGE_UNCOMP; level = 0; break;
        case 'R': reg = strdup(optarg); break;
        case 'l': level = atoi(optarg); break;
        case 'c': flag |= MERGE_COMBINE_RG; break;
        case 'p': flag |= MERGE_COMBINE_PG; break;
        case 's': random_seed = atol(optarg); break;
        case 'X': has_index_file = 1; break; // -X flag for index filename
        case 'b': {
            // load the list of files to read
            if (has_index_file) {
                fprintf(stderr,"Error: The -b option cannot be combined with -X\n");
                ret = 1; goto end;
            }
            int nfiles;
            char **fn_read = hts_readlines(optarg, &nfiles);
            if (fn_read) {
                // Append to end of array
                fn = realloc(fn, (fn_size+nfiles) * sizeof(char*));
                if (fn == NULL) { ret = 1; goto end; }
                memcpy(fn+fn_size, fn_read, nfiles * sizeof(char*));
                fn_size += nfiles;
                free(fn_read);
            }
            else {
                print_error("merge", "Invalid file list \"%s\"", optarg);
                ret = 1;
            }
            break;
        }
        case 1: no_pg = 1; break;
        default:  if (parse_sam_global_opt(c, optarg, lopts, &ga) == 0) break;
                  /* else fall-through */
        case '?': merge_usage(stderr); return 1;
        }
    }
    if ( argc - optind < 1 ) {
        print_error("merge", "You must at least specify the output file");
        merge_usage(stderr);
        return 1;
    }

    if (!no_pg && !(arg_list = stringify_argv(argc+1, argv-1))) {
        print_error("merge", "failed to create arg_list");
        return 1;
    }

    srand48(random_seed);
    if (!(flag & MERGE_FORCE) && strcmp(argv[optind], "-")) {
        FILE *fp = fopen(argv[optind], "rb");
        if (fp != NULL) {
            fclose(fp);
            fprintf(stderr, "[%s] File '%s' exists. Please apply '-f' to overwrite. Abort.\n", __func__, argv[optind]);
            return 1;
        }
    }

    int nargcfiles = 0;
    if (has_index_file) { // Calculate # of input BAM files
        if ((argc - optind - 1) % 2 != 0) {
            fprintf(stderr, "Odd number of filenames detected! Each BAM file should have an index file\n");
            return 1;
        }
        nargcfiles = (argc - optind - 1) / 2;
    } else {
        nargcfiles = argc - optind - 1;
    }

    if (nargcfiles > 0) {
        // Add argc files to end of array
        fn = realloc(fn, (fn_size+nargcfiles) * sizeof(char*));
        if (fn == NULL) { ret = 1; goto end; }
        memcpy(fn+fn_size, argv + (optind+1), nargcfiles * sizeof(char*));

        if(has_index_file) {
            fn_idx = realloc(fn_idx, nargcfiles * sizeof(char*));
            if (fn_idx == NULL) { ret = 1; goto end; }
            memcpy(fn_idx+fn_size, argv + nargcfiles + (optind+1), nargcfiles * sizeof(char*));
        }
    }
    if (fn_size+nargcfiles < 1) {
        print_error("merge", "You must specify at least one (and usually two or more) input files");
        merge_usage(stderr);
        free(fn_idx);
        return 1;
    }
    strcpy(mode, "wb");
    sam_open_mode(mode+1, argv[optind], NULL);
    if (level >= 0) sprintf(strchr(mode, '\0'), "%d", level < 9? level : 9);
    if (bam_merge_core2(is_by_qname, sort_tag, argv[optind], mode, fn_headers,
                        fn_size+nargcfiles, fn, fn_idx, flag, reg, ga.nthreads,
                        "merge", &ga.in, &ga.out, ga.write_index, arg_list, no_pg) < 0)
        ret = 1;

end:
    if (fn_size > 0) {
        int i;
        for (i=0; i<fn_size; i++) free(fn[i]);
    }
    free(fn);
    free(fn_idx);
    free(reg);
    free(arg_list);
    sam_global_args_free(&ga);
    return ret;
}

/***************
 * BAM sorting *
 ***************/

typedef struct {
    size_t from;
    size_t to;
} buf_region;

/* Simplified version of bam_merge_core2() for merging part-sorted
   temporary files.  No need for header merging or translation,
   it just needs to read data into the heap and push it out again. */

static inline int heap_add_read(heap1_t *heap, int nfiles, samFile **fp,
                                int num_in_mem, buf_region *in_mem,
                                bam1_tag *buf, uint64_t *idx, sam_hdr_t *hout) {
    int i = heap->i, res;
    if (i < nfiles) { // read from file
        res = sam_read1(fp[i], hout, heap->entry.bam_record);
    } else { // read from memory
        if (in_mem[i - nfiles].from < in_mem[i - nfiles].to) {
            heap->entry.bam_record = buf[in_mem[i - nfiles].from++].bam_record;
            res = 0;
        } else {
            res = -1;
        }
    }
    if (res >= 0) {
        heap->tid = heap->entry.bam_record->core.tid;
        heap->pos = (uint64_t)(heap->entry.bam_record->core.pos + 1);
        heap->rev = bam_is_rev(heap->entry.bam_record);
        heap->idx = (*idx)++;
        if (g_is_by_tag) {
            heap->entry.u.tag = bam_aux_get(heap->entry.bam_record, g_sort_tag);
        } else {
            heap->entry.u.tag = NULL;
        }
    } else if (res == -1) {
        heap->pos = HEAP_EMPTY;
        if (i < nfiles) bam_destroy1(heap->entry.bam_record);
        heap->entry.bam_record = NULL;
        heap->entry.u.tag = NULL;
    } else {
        return -1;
    }
    return 0;
}

static int bam_merge_simple(int by_qname, char *sort_tag, const char *out,
                            const char *mode, sam_hdr_t *hout,
                            int n, char * const *fn, int num_in_mem,
                            buf_region *in_mem, bam1_tag *buf, int n_threads,
                            const char *cmd, const htsFormat *in_fmt,
                            const htsFormat *out_fmt, char *arg_list, int no_pg,
                            int write_index) {
    samFile *fpout = NULL, **fp = NULL;
    heap1_t *heap = NULL;
    uint64_t idx = 0;
    int i, heap_size = n + num_in_mem;
    char *out_idx_fn = NULL;

    g_is_by_qname = by_qname;
    if (sort_tag) {
        g_is_by_tag = 1;
        g_sort_tag[0] = sort_tag[0];
        g_sort_tag[1] = sort_tag[0] ? sort_tag[1] : '\0';
    }
    if (n > 0) {
        fp = (samFile**)calloc(n, sizeof(samFile*));
        if (!fp) goto mem_fail;
    }
    heap = (heap1_t*)calloc(heap_size, sizeof(heap1_t));
    if (!heap) goto mem_fail;

    // Open each file, read the header and put the first read into the heap
    for (i = 0; i < heap_size; i++) {
        sam_hdr_t *hin;
        heap1_t *h = &heap[i];

        if (i < n) {
            fp[i] = sam_open_format(fn[i], "r", in_fmt);
            if (fp[i] == NULL) {
                print_error_errno(cmd, "fail to open \"%s\"", fn[i]);
                goto fail;
            }

            // Read header ...
            hin = sam_hdr_read(fp[i]);
            if (hin == NULL) {
                print_error(cmd, "failed to read header from \"%s\"", fn[i]);
                goto fail;
            }
            // ... and throw it away as we don't really need it
            sam_hdr_destroy(hin);
        }

        // Get a read into the heap
        h->i = i;
        h->entry.u.tag = NULL;
        if (i < n) {
            h->entry.bam_record = bam_init1();
            if (!h->entry.bam_record) goto mem_fail;
        }
        if (heap_add_read(h, n, fp, num_in_mem, in_mem, buf, &idx, hout) < 0) {
            assert(i < n);
            print_error(cmd, "failed to read first record from \"%s\"", fn[i]);
            goto fail;
        }
    }

    // Open output file and write header
    if ((fpout = sam_open_format(out, mode, out_fmt)) == 0) {
        print_error_errno(cmd, "failed to create \"%s\"", out);
        return -1;
    }

    if (!no_pg && sam_hdr_add_pg(hout, "samtools",
                                 "VN", samtools_version(),
                                 arg_list ? "CL": NULL,
                                 arg_list ? arg_list : NULL,
                                 NULL)) {
        print_error(cmd, "failed to add PG line to the header of \"%s\"", out);
        sam_close(fpout);
        return -1;
    }

    if (n_threads > 1) hts_set_threads(fpout, n_threads);

    if (sam_hdr_write(fpout, hout) != 0) {
        print_error_errno(cmd, "failed to write header to \"%s\"", out);
        sam_close(fpout);
        return -1;
    }

    if (write_index) {
        if (!(out_idx_fn = auto_index(fpout, out, hout))){
            sam_close(fpout);
            return -1;
        }
    }

    // Now do the merge
    ks_heapmake(heap, heap_size, heap);
    while (heap->pos != HEAP_EMPTY) {
        bam1_t *b = heap->entry.bam_record;
        if (sam_write1(fpout, hout, b) < 0) {
            print_error_errno(cmd, "failed writing to \"%s\"", out);
            goto fail;
        }
        if (heap_add_read(heap, n, fp, num_in_mem, in_mem, buf, &idx, hout) < 0) {
            assert(heap->i < n);
            print_error(cmd, "Error reading \"%s\" : %s",
                        fn[heap->i], strerror(errno));
            goto fail;
        }
        ks_heapadjust(heap, 0, heap_size, heap);
    }
    // Clean up and close
    for (i = 0; i < n; i++) {
        if (sam_close(fp[i]) != 0) {
            print_error(cmd, "Error on closing \"%s\" : %s",
                        fn[i], strerror(errno));
        }
    }
    free(fp);
    free(heap);

    if (write_index) {
        if (sam_idx_save(fpout) < 0) {
            print_error_errno("merge", "writing index failed");
            goto fail;
        }
        free(out_idx_fn);
    }

    if (sam_close(fpout) < 0) {
        print_error(cmd, "error closing output file");
        return -1;
    }
    return 0;
 mem_fail:
    print_error(cmd, "Out of memory");

 fail:
    for (i = 0; i < n; i++) {
        if (fp && fp[i]) sam_close(fp[i]);
    }
    for (i = 0; i < heap_size; i++) {
        if (heap && heap[i].i < n && heap[i].entry.bam_record)
            bam_destroy1(heap[i].entry.bam_record);
    }
    free(fp);
    free(heap);
    if (fpout) sam_close(fpout);
    free(out_idx_fn);
    return -1;
}

// Function to compare reads and determine which one is < or > the other
// Handle sort-by-pos and sort-by-name. Used as the secondary sort in bam1_lt_by_tag, if reads are equivalent by tag.
// Returns a value less than, equal to or greater than zero if a is less than,
// equal to or greater than b, respectively.
static inline int bam1_cmp_core(const bam1_tag a, const bam1_tag b)
{
    uint64_t pa, pb;
    if (!a.bam_record)
        return 1;
    if (!b.bam_record)
        return 0;

    if (g_is_by_qname) {
        int t = strnum_cmp(bam_get_qname(a.bam_record), bam_get_qname(b.bam_record));
        if (t != 0) return t;
        return (int) (a.bam_record->core.flag&0xc0) - (int) (b.bam_record->core.flag&0xc0);
    } else {
        pa = a.bam_record->core.tid;
        pb = b.bam_record->core.tid;

        if (pa == pb) {
            pa = (uint64_t)(a.bam_record->core.pos+1);
            pb = (uint64_t)(b.bam_record->core.pos+1);
        }

        if (pa == pb) {
            pa = bam_is_rev(a.bam_record);
            pb = bam_is_rev(b.bam_record);
        }

        return pa < pb ? -1 : (pa > pb ? 1 : 0);
    }
}

uint8_t normalize_type(const uint8_t* aux) {
    if (*aux == 'c' || *aux == 'C' || *aux == 's' || *aux == 'S' || *aux == 'i' || *aux == 'I') {
        return 'c';
    } else if (*aux == 'f' || *aux == 'd') {
        return 'f';
    } else if (*aux == 'H' || *aux == 'Z') {
         return 'H';
    } else {
        return *aux;
    }
}

// Sort record by tag, using pos or read name as a secondary key if tags are identical. Reads not carrying the tag sort first.
// Tags are first sorted by the type character (in case the types differ), or by the appropriate comparator for that type if they agree.
// Returns a value less than, equal to or greater than zero if a is less than,
// equal to or greater than b, respectively.
static inline int bam1_cmp_by_tag(const bam1_tag a, const bam1_tag b)
{
    const uint8_t* aux_a = a.u.tag;
    const uint8_t* aux_b = b.u.tag;

    if (aux_a == NULL && aux_b != NULL) {
        return -1;
    } else if (aux_a != NULL && aux_b == NULL) {
        return 1;
    } else if (aux_a == NULL && aux_b == NULL) {
        return bam1_cmp_core(a,b);
    }

    // 'Normalize' the letters of the datatypes to a canonical letter,
    // so that comparison of different types
    // forms a correct total ordering.
    uint8_t a_type = normalize_type(aux_a);
    uint8_t b_type = normalize_type(aux_b);

    if (a_type != b_type) {
        // Fix int to float comparisons by using bam_aux2f() to read the int
        if (a_type == 'c' && b_type == 'f') {
            a_type = 'f';
        } else if (a_type == 'f' && b_type == 'c') {
            b_type = 'f';
        } else {
            // Unfixable mismatched types
            return a_type < b_type ? -1 : 1;
        }
    }

    if (a_type == 'c') {
        int64_t va = bam_aux2i(aux_a);
        int64_t vb = bam_aux2i(aux_b);
        if (va != vb) return va < vb ? -1 : 1;
        return bam1_cmp_core(a, b);
    } else if (a_type == 'f') {
        double va = bam_aux2f(aux_a);
        double vb = bam_aux2f(aux_b);
        if (va != vb) return va < vb ? -1 : 1;
        return bam1_cmp_core(a, b);
    } else if (a_type == 'A') {
        unsigned char va = bam_aux2A(aux_a);
        unsigned char vb = bam_aux2A(aux_b);
        if (va != vb) return va < vb ? -1 : 1;
        return bam1_cmp_core(a, b);
    } else if (a_type == 'H') {
        int t = strcmp(bam_aux2Z(aux_a), bam_aux2Z(aux_b));
        if (t) return t;
        return bam1_cmp_core(a, b);
    } else {
        return bam1_cmp_core(a,b);
    }
}

// Function to compare reads and determine which one is < the other
// Handle sort-by-pos, sort-by-name, or sort-by-tag
static inline int bam1_lt(const bam1_tag a, const bam1_tag b)
{
    if (g_is_by_tag) {
        return bam1_cmp_by_tag(a, b) < 0;
    } else {
        return bam1_cmp_core(a,b) < 0;
    }
}



KSORT_INIT(sort, bam1_tag, bam1_lt)

typedef struct {
    size_t buf_len;
    const char *prefix;
    bam1_tag *buf;
    const sam_hdr_t *h;
    int index;
    int error;
    int no_save;
} worker_t;

// Returns 0 for success
//        -1 for failure
static int write_buffer(const char *fn, const char *mode, size_t l, bam1_tag *buf,
                        const sam_hdr_t *h, int n_threads, const htsFormat *fmt,
                        char *arg_list, int no_pg, int write_index)
{
    size_t i;
    samFile* fp;
    char *out_idx_fn = NULL;

    fp = sam_open_format(fn, mode, fmt);
    if (fp == NULL) return -1;
    if (!no_pg && sam_hdr_add_pg((sam_hdr_t *)h, "samtools",
                                 "VN", samtools_version(),
                                 arg_list ? "CL": NULL,
                                 arg_list ? arg_list : NULL,
                                 NULL)) {
        goto fail;
    }
    if (sam_hdr_write(fp, (sam_hdr_t *)h) != 0) goto fail;

    if (write_index) {
        if (!(out_idx_fn = auto_index(fp, fn, (sam_hdr_t *)h))) goto fail;
    }

    if (n_threads > 1) hts_set_threads(fp, n_threads);
    for (i = 0; i < l; ++i) {
        if (sam_write1(fp, (sam_hdr_t *)h, buf[i].bam_record) < 0) goto fail;
    }

    if (write_index) {
        if (sam_idx_save(fp) < 0) {
            print_error_errno("merge", "writing index failed");
            goto fail;
        }
        free(out_idx_fn);
    }


    if (sam_close(fp) < 0) return -1;
    return 0;
 fail:
    sam_close(fp);
    free(out_idx_fn);
    return -1;
}

#define NUMBASE 256

static int ks_radixsort(size_t n, bam1_tag *buf, const sam_hdr_t *h)
{
    int curr = 0, ret = -1;
    ssize_t i;
    bam1_tag *buf_ar2[2], *bam_a, *bam_b;
    uint64_t max_pos = 1;
    uint32_t max_tid = 1, tid_bytes = 0, pos_bytes = 0, byte = 0;
    uint32_t tid_shift_l, tid_shift_r;
    int nref = sam_hdr_nref(h);

    // Count number of bytes needed for biggest tid and pos
    //  Notes: Add 1 to core.pos so always positive.
    //         Convert unmapped tid (-1) to number of references so unmapped
    //         sort to the end.
    for (i = 0; i < n; i++) {
        bam1_t *b = buf[i].bam_record;
        uint32_t tid = b->core.tid == -1 ? nref : b->core.tid;
        uint64_t pos = ((uint64_t)(b->core.pos + 1) << 1) | bam_is_rev(b);
        if (max_tid < tid)
            max_tid = tid;
        if (max_pos < pos)
            max_pos = pos;
    }

    for (; max_pos > 0; max_pos >>= 8) pos_bytes++;
    for (; max_tid > 0; max_tid >>= 8) tid_bytes++;
    assert(pos_bytes + tid_bytes < sizeof(buf[0].u.pos_tid));

    tid_shift_l = pos_bytes * 8;
    tid_shift_r = 64 - tid_shift_l;

    // Write position and tid into bam1_tag::u::pos_tid using minimum number
    // of bytes required.  Values are stored little-endian so that we
    // get a least-significant digit (byte) radix sort.
    for (i = 0; i < n; i++) {
        bam1_t *b = buf[i].bam_record;
        uint32_t tid = b->core.tid == -1 ? nref : b->core.tid;
        // 'pos' here includes as many bytes of tid as will fit
        // in the space remaining above pos_bytes.  The rest of tid
        // is written out separately.
        uint64_t pos = (bam_is_rev(b) |
                        ((uint64_t)(b->core.pos + 1) << 1) |
                        (tid_shift_l < 64 ? (uint64_t) tid << tid_shift_l : 0));
        u64_to_le(pos, buf[i].u.pos_tid);
        u32_to_le(tid_shift_r < 32 ? tid >> tid_shift_r : 0,
                  &buf[i].u.pos_tid[8]);
    }

    buf_ar2[0] = buf;
    buf_ar2[1] = (bam1_tag *)malloc(sizeof(bam1_tag) * n);
    if (buf_ar2[1] == NULL) {
        print_error("sort", "couldn't allocate memory for temporary buf");
        goto err;
    }

    // Least-significant digit radix sort (where "digits" are bytes)
    for (byte = 0; byte < pos_bytes + tid_bytes; byte++) {
        size_t remainders[NUMBASE] = { 0 };
        bam_a = buf_ar2[curr]; bam_b = buf_ar2[1-curr];
        for (i = 0; i < n; ++i)
            remainders[bam_a[i].u.pos_tid[byte]]++;
        for (i = 1; i < NUMBASE; ++i)
            remainders[i] += remainders[i - 1];
        for (i = n - 1; i >= 0; i--) {
            size_t j = --remainders[bam_a[i].u.pos_tid[byte]];
            bam_b[j] = bam_a[i];
        }
        curr = 1 - curr;
    }
    if (curr == 1) {
        bam1_tag *end = buf + n;
        bam_a = buf_ar2[0]; bam_b = buf_ar2[1];
        while (bam_a < end) *bam_a++ = *bam_b++;
    }

    ret = 0;
err:
    free(buf_ar2[1]);
    return ret;
}

static void *worker(void *data)
{
    worker_t *w = (worker_t*)data;
    char *name;
    w->error = 0;

    if (!g_is_by_qname && !g_is_by_tag) {
        if (ks_radixsort(w->buf_len, w->buf, w->h) < 0) {
            w->error = errno;
            return NULL;
        }
    } else {
        ks_mergesort(sort, w->buf_len, w->buf, 0);
    }

    if (w->no_save)
        return 0;

    name = (char*)calloc(strlen(w->prefix) + 20, 1);
    if (!name) { w->error = errno; return 0; }
    sprintf(name, "%s.%.4d.bam", w->prefix, w->index);

    uint32_t max_ncigar = 0;
    int i;
    for (i = 0; i < w->buf_len; i++) {
        uint32_t nc = w->buf[i].bam_record->core.n_cigar;
        if (max_ncigar < nc)
            max_ncigar = nc;
    }

    if (max_ncigar > 65535) {
        htsFormat fmt;
        memset(&fmt, 0, sizeof(fmt));
        if (hts_parse_format(&fmt, "cram,version=3.0,no_ref,seqs_per_slice=1000") < 0) {
            w->error = errno;
            free(name);
            return 0;
        }

        if (write_buffer(name, "wcx1", w->buf_len, w->buf, w->h, 0, &fmt, NULL, 1, 0) < 0)
            w->error = errno;
    } else {
        if (write_buffer(name, "wbx1", w->buf_len, w->buf, w->h, 0, NULL, NULL, 1, 0) < 0)
            w->error = errno;
    }

    free(name);
    return 0;
}

static int sort_blocks(int n_files, size_t k, bam1_tag *buf, const char *prefix,
                       const sam_hdr_t *h, int n_threads, buf_region *in_mem)
{
    int i;
    size_t pos, rest;
    pthread_t *tid;
    pthread_attr_t attr;
    worker_t *w;
    int n_failed = 0;

    if (n_threads < 1) n_threads = 1;
    if (k < n_threads * 64) n_threads = 1; // use a single thread if we only sort a small batch of records
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    w = (worker_t*)calloc(n_threads, sizeof(worker_t));
    if (!w) return -1;
    tid = (pthread_t*)calloc(n_threads, sizeof(pthread_t));
    if (!tid) { free(w); return -1; }
    pos = 0; rest = k;
    for (i = 0; i < n_threads; ++i) {
        w[i].buf_len = rest / (n_threads - i);
        w[i].buf = &buf[pos];
        w[i].prefix = prefix;
        w[i].h = h;
        w[i].index = n_files + i;
        if (in_mem) {
            w[i].no_save = 1;
            in_mem[i].from = pos;
            in_mem[i].to = pos + w[i].buf_len;
        } else {
            w[i].no_save = 0;
        }
        pos += w[i].buf_len; rest -= w[i].buf_len;
        pthread_create(&tid[i], &attr, worker, &w[i]);
    }
    for (i = 0; i < n_threads; ++i) {
        pthread_join(tid[i], 0);
        if (w[i].error != 0) {
            errno = w[i].error;
            print_error_errno("sort", "failed to create temporary file \"%s.%.4d.bam\"", prefix, w[i].index);
            n_failed++;
        }
    }
    free(tid); free(w);
    if (n_failed) return -1;
    if (in_mem) return n_threads;
    return n_files + n_threads;
}

/*!
  @abstract Sort an unsorted BAM file based on the chromosome order
  and the leftmost position of an alignment

  @param  is_by_qname whether to sort by query name
  @param  sort_by_tag if non-null, sort by the given tag
  @param  fn       name of the file to be sorted
  @param  prefix   prefix of the temporary files (prefix.NNNN.bam are written)
  @param  fnout    name of the final output file to be written
  @param  modeout  sam_open() mode to be used to create the final output file
  @param  max_mem  approxiate maximum memory (very inaccurate)
  @param  in_fmt   input file format options
  @param  out_fmt  output file format and options
  @param  arg_list    command string for PG line
  @param  no_pg       if 1, do not add a new PG line
  @paran  write_index create index for the output file
  @return 0 for successful sorting, negative on errors

  @discussion It may create multiple temporary subalignment files
  and then merge them by calling bam_merge_simple(). This function is
  NOT thread safe.
 */
int bam_sort_core_ext(int is_by_qname, char* sort_by_tag, const char *fn, const char *prefix,
                      const char *fnout, const char *modeout,
                      size_t _max_mem, int n_threads,
                      const htsFormat *in_fmt, const htsFormat *out_fmt,
                      char *arg_list, int no_pg, int write_index)
{
    int ret = -1, res, i, n_files = 0;
    size_t max_k, k, max_mem, bam_mem_offset;
    sam_hdr_t *header = NULL;
    samFile *fp;
    bam1_tag *buf = NULL;
    bam1_t *b = bam_init1();
    uint8_t *bam_mem = NULL;
    char **fns = NULL;
    const char *new_so;
    buf_region *in_mem = NULL;
    int num_in_mem = 0;

    if (!b) {
        print_error("sort", "couldn't allocate memory for bam record");
        return -1;
    }

    if (n_threads < 2) n_threads = 1;
    g_is_by_qname = is_by_qname;
    if (sort_by_tag) {
        g_is_by_tag = 1;
        g_sort_tag[0] = sort_by_tag[0];
        g_sort_tag[1] = sort_by_tag[0] ? sort_by_tag[1] : '\0';
    }

    max_mem = _max_mem * n_threads;
    buf = NULL;
    fp = sam_open_format(fn, "r", in_fmt);
    if (fp == NULL) {
        print_error_errno("sort", "can't open \"%s\"", fn);
        goto err;
    }
    header = sam_hdr_read(fp);
    if (header == NULL) {
        print_error("sort", "failed to read header from \"%s\"", fn);
        goto err;
    }

    if (sort_by_tag != NULL)
        new_so = "unknown";
    else if (is_by_qname)
        new_so = "queryname";
    else
        new_so = "coordinate";

    if ((-1 == sam_hdr_update_hd(header, "SO", new_so))
     && (-1 == sam_hdr_add_line(header, "HD", "VN", SAM_FORMAT_VERSION, "SO", new_so, NULL))
     ) {
        print_error("sort", "failed to change sort order header to '%s'\n", new_so);
        goto err;
    }

    if (-1 == sam_hdr_remove_tag_hd(header, "GO")) {
        print_error("sort", "failed to delete group order header\n");
        goto err;
    }

    // No gain to using the thread pool here as the flow of this code
    // is such that we are *either* reading *or* sorting.  Hence a shared
    // pool makes no real difference except to reduce the thread count a little.
    if (n_threads > 1)
        hts_set_threads(fp, n_threads);

    if ((bam_mem = malloc(max_mem)) == NULL) {
        print_error("sort", "couldn't allocate memory for bam_mem");
        goto err;
    }

    // write sub files
    k = max_k = bam_mem_offset = 0;
    while ((res = sam_read1(fp, header, b)) >= 0) {
        int mem_full = 0;

        if (k == max_k) {
            bam1_tag *new_buf;
            max_k = max_k? max_k<<1 : 0x10000;
            if ((new_buf = realloc(buf, max_k * sizeof(bam1_tag))) == NULL) {
                print_error("sort", "couldn't allocate memory for buf");
                goto err;
            }
            buf = new_buf;
        }

        // Check if the BAM record will fit in the memory limit
        if (bam_mem_offset + sizeof(*b) + b->l_data < max_mem) {
            // Copy record into the memory block
            buf[k].bam_record = (bam1_t *)(bam_mem + bam_mem_offset);
            *buf[k].bam_record = *b;
            buf[k].bam_record->data = (uint8_t *)((char *)buf[k].bam_record + sizeof(bam1_t));
            memcpy(buf[k].bam_record->data, b->data, b->l_data);
            // store next BAM record in next 8-byte-aligned address after
            // current one
            bam_mem_offset = (bam_mem_offset + sizeof(*b) + b->l_data + 8 - 1) & ~((size_t)(8 - 1));
        } else {
            // Add a pointer to the remaining record
            buf[k].bam_record = b;
            mem_full = 1;
        }

        // Pull out the value of the position
        // or the pointer to the sort tag if applicable
        if (g_is_by_tag) {
            buf[k].u.tag = bam_aux_get(buf[k].bam_record, g_sort_tag);
        } else {
            buf[k].u.tag = NULL;
        }
        ++k;

        if (mem_full) {
            n_files = sort_blocks(n_files, k, buf, prefix, header, n_threads,
                                  NULL);
            if (n_files < 0) {
                goto err;
            }
            k = 0;
            bam_mem_offset = 0;
        }
    }
    if (res != -1) {
        print_error("sort", "truncated file. Aborting");
        goto err;
    }

    // Sort last records
    if (k > 0) {
        in_mem = calloc(n_threads > 0 ? n_threads : 1, sizeof(in_mem[0]));
        if (!in_mem) goto err;
        num_in_mem = sort_blocks(n_files, k, buf, prefix, header, n_threads,
                                 in_mem);
        if (num_in_mem < 0) goto err;
    } else {
        num_in_mem = 0;
    }

    // write the final output
    if (n_files == 0 && num_in_mem < 2) { // a single block
        if (write_buffer(fnout, modeout, k, buf, header, n_threads, out_fmt, arg_list, no_pg, write_index) != 0) {
            print_error_errno("sort", "failed to create \"%s\"", fnout);
            goto err;
        }
    } else { // then merge
        fprintf(stderr,
                "[bam_sort_core] merging from %d files and %d in-memory blocks...\n",
                n_files, num_in_mem);
        fns = (char**)calloc(n_files, sizeof(char*));
        if (!fns) goto err;
        for (i = 0; i < n_files; ++i) {
            fns[i] = (char*)calloc(strlen(prefix) + 20, 1);
            if (!fns[i]) goto err;
            sprintf(fns[i], "%s.%.4d.bam", prefix, i);
        }
        if (bam_merge_simple(is_by_qname, sort_by_tag, fnout, modeout, header,
                             n_files, fns, num_in_mem, in_mem, buf,
                             n_threads, "sort", in_fmt, out_fmt, arg_list,
                             no_pg, write_index) < 0) {
            // Propagate bam_merge_simple() failure; it has already emitted a
            // message explaining the failure, so no further message is needed.
            goto err;
        }
    }

    ret = 0;

 err:
    // free
    if (fns) {
        for (i = 0; i < n_files; ++i) {
            if (fns[i]) {
                unlink(fns[i]);
                free(fns[i]);
            }
        }
        free(fns);
    }
    bam_destroy1(b);
    free(buf);
    free(bam_mem);
    free(in_mem);
    sam_hdr_destroy(header);
    if (fp) sam_close(fp);
    return ret;
}

// Unused here but may be used by legacy samtools-using third-party code
int bam_sort_core(int is_by_qname, const char *fn, const char *prefix, size_t max_mem)
{
    int ret;
    char *fnout = calloc(strlen(prefix) + 4 + 1, 1);
    if (!fnout) return -1;
    sprintf(fnout, "%s.bam", prefix);
    ret = bam_sort_core_ext(is_by_qname, NULL, fn, prefix, fnout, "wb", max_mem, 0, NULL, NULL, NULL, 1, 0);
    free(fnout);
    return ret;
}

static void sort_usage(FILE *fp)
{
    fprintf(fp,
"Usage: samtools sort [options...] [in.bam]\n"
"Options:\n"
"  -l INT     Set compression level, from 0 (uncompressed) to 9 (best)\n"
"  -m INT     Set maximum memory per thread; suffix K/M/G recognized [768M]\n"
"  -n         Sort by read name\n"
"  -t TAG     Sort by value of TAG. Uses position as secondary index (or read name if -n is set)\n"
"  -o FILE    Write final output to FILE rather than standard output\n"
"  -T PREFIX  Write temporary files to PREFIX.nnnn.bam\n"
"  --no-PG    do not add a PG line\n");
    sam_global_opt_help(fp, "-.O..@-.");
}

static void complain_about_memory_setting(size_t max_mem) {
    char  *suffix = "";
    const size_t nine_k = 9<<10;
    if (max_mem > nine_k) { max_mem >>= 10; suffix = "K"; }
    if (max_mem > nine_k) { max_mem >>= 10; suffix = "M"; }

    fprintf(stderr,
"[bam_sort] -m setting (%zu%s bytes) is less than the minimum required (%zuM).\n\n"
"Trying to run with -m too small can lead to the creation of a very large number\n"
"of temporary files.  This may make sort fail due to it exceeding limits on the\n"
"number of files it can have open at the same time.\n\n"
"Please check your -m parameter.  It should be an integer followed by one of the\n"
"letters K (for kilobytes), M (megabytes) or G (gigabytes).  You should ensure it\n"
"is at least the minimum above, and much higher if you are sorting a large file.\n",
            max_mem, suffix, SORT_MIN_MEGS_PER_THREAD);
}

int bam_sort(int argc, char *argv[])
{
    size_t max_mem = SORT_DEFAULT_MEGS_PER_THREAD << 20;
    int c, nargs, is_by_qname = 0, ret, o_seen = 0, level = -1, no_pg = 0;
    char* sort_tag = NULL, *arg_list = NULL;
    char *fnout = "-", modeout[12];
    kstring_t tmpprefix = { 0, 0, NULL };
    struct stat st;
    sam_global_args ga = SAM_GLOBAL_ARGS_INIT;

    static const struct option lopts[] = {
        SAM_OPT_GLOBAL_OPTIONS('-', 0, 'O', 0, 0, '@'),
        { "threads", required_argument, NULL, '@' },
        {"no-PG", no_argument, NULL, 1},
        { NULL, 0, NULL, 0 }
    };

    while ((c = getopt_long(argc, argv, "l:m:no:O:T:@:t:", lopts, NULL)) >= 0) {
        switch (c) {
        case 'o': fnout = optarg; o_seen = 1; break;
        case 'n': is_by_qname = 1; break;
        case 't': sort_tag = optarg; break;
        case 'm': {
                char *q;
                max_mem = strtol(optarg, &q, 0);
                if (*q == 'k' || *q == 'K') max_mem <<= 10;
                else if (*q == 'm' || *q == 'M') max_mem <<= 20;
                else if (*q == 'g' || *q == 'G') max_mem <<= 30;
                break;
            }
        case 'T': kputs(optarg, &tmpprefix); break;
        case 'l': level = atoi(optarg); break;
        case 1: no_pg = 1; break;

        default:  if (parse_sam_global_opt(c, optarg, lopts, &ga) == 0) break;
                  /* else fall-through */
        case '?': sort_usage(stderr); ret = EXIT_FAILURE; goto sort_end;
        }
    }

    nargs = argc - optind;
    if (nargs == 0 && isatty(STDIN_FILENO)) {
        sort_usage(stdout);
        ret = EXIT_SUCCESS;
        goto sort_end;
    }
    else if (nargs >= 2) {
        // If exactly two, user probably tried to specify legacy <out.prefix>
        if (nargs == 2)
            fprintf(stderr, "[bam_sort] Use -T PREFIX / -o FILE to specify temporary and final output files\n");

        sort_usage(stderr);
        ret = EXIT_FAILURE;
        goto sort_end;
    }

    if (ga.write_index && (is_by_qname || sort_tag)) {
        fprintf(stderr, "[W::bam_sort] Ignoring --write-index as it only works for position sorted files.\n");
        ga.write_index = 0;
    }

    if (!no_pg && !(arg_list = stringify_argv(argc+1, argv-1))) {
        print_error("sort", "failed to create arg_list");
        return 1;
    }

    if (max_mem < (SORT_MIN_MEGS_PER_THREAD << 20)) {
        complain_about_memory_setting(max_mem);
        ret = EXIT_FAILURE;
        goto sort_end;
    }

    strcpy(modeout, "wb");
    sam_open_mode(modeout+1, fnout, NULL);
    if (level >= 0) sprintf(strchr(modeout, '\0'), "%d", level < 9? level : 9);

    if (tmpprefix.l == 0) {
        if (strcmp(fnout, "-") != 0) ksprintf(&tmpprefix, "%s.tmp", fnout);
        else kputc('.', &tmpprefix);
    }
    if (stat(tmpprefix.s, &st) == 0 && S_ISDIR(st.st_mode)) {
        unsigned t = ((unsigned) time(NULL)) ^ ((unsigned) clock());
        if (tmpprefix.s[tmpprefix.l-1] != '/') kputc('/', &tmpprefix);
        ksprintf(&tmpprefix, "samtools.%d.%u.tmp", (int) getpid(), t % 10000);
    }

    ret = bam_sort_core_ext(is_by_qname, sort_tag, (nargs > 0)? argv[optind] : "-",
                            tmpprefix.s, fnout, modeout, max_mem, ga.nthreads,
                            &ga.in, &ga.out, arg_list, no_pg, ga.write_index);
    if (ret >= 0)
        ret = EXIT_SUCCESS;
    else {
        char dummy[4];
        // If we failed on opening the input file & it has no .bam/.cram/etc
        // extension, the user probably tried legacy -o <infile> <out.prefix>
        if (ret == -2 && o_seen && nargs > 0 && sam_open_mode(dummy, argv[optind], NULL) < 0)
            fprintf(stderr, "[bam_sort] Note the <out.prefix> argument has been replaced by -T/-o options\n");

        ret = EXIT_FAILURE;
    }

sort_end:
    free(tmpprefix.s);
    free(arg_list);
    sam_global_args_free(&ga);

    return ret;
}
