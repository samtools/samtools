/*  bam_sort.c -- sorting and merging.

    Copyright (C) 2008-2016 Genome Research Ltd.
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
#include "sam_opts.h"
#include "samtools.h"


// Struct which contains the a record, and the pointer to the sort tag (if any) or
// a combined ref / position / strand.
// Used to speed up tag and position sorts.
typedef struct bam1_tag {
    bam1_t *bam_record;
    union {
        const uint8_t *tag;
        uint64_t pos;
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

#define HEAP_EMPTY UINT64_MAX

typedef struct {
    int i;
    uint32_t rev;
    uint64_t pos, idx;
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
        if (a.pos != b.pos) return a.pos > b.pos;
        if (a.rev != b.rev) return a.rev > b.rev;
    }
    // This compares by position in the input file(s)
    if (a.i != b.i) return a.i > b.i;
    return a.idx > b.idx;
}

KSORT_INIT(heap, heap1_t, heap_lt)

typedef struct merged_header {
    kstring_t     out_hd;
    kstring_t     out_sq;
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

/* Something to look like a regmatch_t */
typedef struct hdr_match {
    ptrdiff_t rm_so;
    ptrdiff_t rm_eo;
} hdr_match_t;

/*
 * Search for header lines of a particular record type.
 *
 * This replaces a regex search for something like /^@SQ.*\tSN:([^\t]+).*$/
 * but is much quicker.  The locations found are returned in *matches,
 * which has a signature the same as that of a regmatch_t.
 *
 * rec is the record type to match (i.e. @HD, @SQ, @PG or @RG)
 * tag is a tag type in the record to match (SN for @SQ, ID for @PG or @RG)
 *
 * The location of the record (if found) is returned in matches[0]
 * If tag is not NULL, the record is searched for the presence of the
 * given tag.  If found, the location of the value is returned in matches[1].
 * If the tag isn't found then the record is ignored and the search resumes
 * on the next header line.
 *
 * For simplicity, some assumptions are made about rec and tag:
 *   rec should include the leading '@' sign and be three characters long.
 *   tag should be exactly two characters long.
 * These are always string constants when this is called below, so we don't
 * bother to check here.
 *
 * Returns 0 if a match was found, -1 if not.
 */


static int hdr_line_match(const char *text, const char *rec,
                          const char *tag,  hdr_match_t *matches) {
    const char *line_start, *line_end = text;
    const char *tag_start, *tag_end;

    for (;;) {
        // Find record, ensure either at start of text or follows '\n'
        line_start = strstr(line_end, rec);
        while (line_start && line_start > text && *(line_start - 1) != '\n') {
            line_start = strstr(line_start + 3, rec);
        }
        if (!line_start) return -1;

        // Find end of header line
        line_end = strchr(line_start, '\n');
        if (!line_end) line_end = line_start + strlen(line_start);

        matches[0].rm_so = line_start - text;
        matches[0].rm_eo = line_end - text;
        if (!tag) return 0;  // Match found if not looking for tag.

        for (tag_start = line_start + 3; tag_start < line_end; tag_start++) {
            // Find possible tag start.  Hacky but quick.
            while (*tag_start > '\n') tag_start++;

            // Check it
            if (tag_start[0] == '\t'
                && strncmp(tag_start + 1, tag, 2) == 0
                && tag_start[3] == ':') {
                // Found tag, record location and return.
                tag_end = tag_start + 4;
                while (*tag_end && *tag_end != '\t' && *tag_end != '\n')
                    ++tag_end;
                matches[1].rm_so = tag_start - text + 4;
                matches[1].rm_eo = tag_end - text;
                return 0;
            }
        }
        // Couldn't find tag, try again from end of current record.
    }
}

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
    free(merged_hdr);
    return NULL;
}

/* Some handy kstring manipulating functions */

// Append char range to kstring
static inline int range_to_ks(const char *src, int from, int to,
                              kstring_t *dest) {
    return kputsn(src + from, to - from, dest) != to - from;
}

// Append a header line match to kstring
static inline int match_to_ks(const char *src, const hdr_match_t *match,
                              kstring_t *dest) {
    return range_to_ks(src, match->rm_so, match->rm_eo, dest);
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
                            bam_hdr_t *translate) {
    hdr_match_t match = {0, 0};

    // TODO: handle case when @HD needs merging.
    if (merged_hdr->have_hd) return 0;

    if (hdr_line_match(translate->text, "@HD", NULL, &match) != 0) {
        return 0;
    }

    if (match_to_ks(translate->text, &match, &merged_hdr->out_hd)) goto memfail;
    if (kputc('\n', &merged_hdr->out_hd) == EOF) goto memfail;
    merged_hdr->have_hd = true;

    return 0;

 memfail:
    perror(__func__);
    return -1;
}

static inline int grow_target_list(merged_header_t* merged_hdr) {
    size_t     new_size;
    char     **new_names;
    uint32_t  *new_len;

    new_size = merged_hdr->targets_sz * 2;
    new_names = realloc(merged_hdr->target_name, sizeof(*new_names) * new_size);
    if (!new_names) goto fail;
    merged_hdr->target_name = new_names;

    new_len = realloc(merged_hdr->target_len, sizeof(*new_len) * new_size);
    if (!new_len) goto fail;
    merged_hdr->target_len = new_len;

    merged_hdr->targets_sz = new_size;

    return 0;

 fail:
    perror(__func__);
    return -1;
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

static int trans_tbl_add_sq(merged_header_t* merged_hdr, bam_hdr_t *translate,
                            trans_tbl_t* tbl) {

    kstring_t *out_text = &merged_hdr->out_sq;
    khash_t(c2i)* sq_tids = merged_hdr->sq_tids;
    hdr_match_t *new_sq_matches = NULL;
    char *text;
    hdr_match_t matches[2];
    int32_t i;
    int32_t old_n_targets = merged_hdr->n_targets;
    khiter_t iter;
    int min_tid = -1;

    // Fill in the tid part of the translation table, adding new targets
    // to the merged header as we go.

    for (i = 0; i < translate->n_targets; ++i) {

        // Check if it's a new target.
        iter = kh_get(c2i, sq_tids, translate->target_name[i]);

        if (iter == kh_end(sq_tids)) {
            int ret;
            // Append missing entries to out_hdr

            if (merged_hdr->n_targets == merged_hdr->targets_sz) {
                if (grow_target_list(merged_hdr)) goto fail;
            }

            merged_hdr->target_name[merged_hdr->n_targets] = strdup(translate->target_name[i]);
            if (merged_hdr->target_name[merged_hdr->n_targets] == NULL) goto memfail;
            merged_hdr->target_len[merged_hdr->n_targets] = translate->target_len[i];

            // Record the new identifier for reference below,
            // and when building the ttable for other inputs.
            iter = kh_put(c2i, sq_tids,
                          merged_hdr->target_name[merged_hdr->n_targets], &ret);
            if (ret < 0) {
                free(merged_hdr->target_name[merged_hdr->n_targets]);
                goto memfail;
            }
            assert(ret > 0);  // Should not be in hash already.

            kh_value(sq_tids, iter) = merged_hdr->n_targets;
            tbl->tid_trans[i] = merged_hdr->n_targets++;
        } else {
            tbl->tid_trans[i] = kh_value(sq_tids, iter);
        }

        if (tbl->tid_trans[i] > min_tid) {
            min_tid = tbl->tid_trans[i];
        } else {
            tbl->lost_coord_sort = true;
        }
    }

    if (merged_hdr->n_targets == old_n_targets)
        return 0;  // Everything done if no new targets.

    // Otherwise, find @SQ lines in translate->text for all newly added targets.

    new_sq_matches = malloc((merged_hdr->n_targets - old_n_targets)
                            * sizeof(*new_sq_matches));
    if (new_sq_matches == NULL) goto memfail;

    for (i = 0; i < merged_hdr->n_targets - old_n_targets; i++) {
        new_sq_matches[i].rm_so = new_sq_matches[i].rm_eo = -1;
    }

    text = translate->text;
    while (hdr_line_match(text, "@SQ", "SN", matches) == 0) {
        // matches[0] is whole line, matches[1] is SN value.

        // This is a bit disgusting, but avoids a copy...
        char c = text[matches[1].rm_eo];
        int idx;

        text[matches[1].rm_eo] = '\0';

        // Look up the SN value in the sq_tids hash.
        iter = kh_get(c2i, sq_tids, text + matches[1].rm_so);
        text[matches[1].rm_eo] = c; // restore text

        if (iter == kh_end(sq_tids)) {
            // Warn about this, but it's not really fatal.
            fprintf(stderr, "[W::%s] @SQ SN (%.*s) found in text header but not binary header.\n",
                    __func__,
                    (int) (matches[1].rm_eo - matches[1].rm_so),
                    text + matches[1].rm_so);
            text += matches[0].rm_eo;
            continue;  // Skip to next
        }

        idx = kh_value(sq_tids, iter);
        if (idx >= old_n_targets) {
            // is a new SQ, so record position so we can add it to out_text.
            assert(idx < merged_hdr->n_targets);
            ptrdiff_t off = text - translate->text;
            new_sq_matches[idx - old_n_targets].rm_so = matches[0].rm_so + off;
            new_sq_matches[idx - old_n_targets].rm_eo = matches[0].rm_eo + off;
        }

        // Carry on searching from end of current match
        text += matches[0].rm_eo;
    }

    // Copy the @SQ headers found and recreate any missing from binary header.
    for (i = 0; i < merged_hdr->n_targets - old_n_targets; i++) {
        if (new_sq_matches[i].rm_so >= 0) {
            if (match_to_ks(translate->text, &new_sq_matches[i], out_text))
                goto memfail;
            if (kputc('\n', out_text) == EOF) goto memfail;
        } else {
            if (kputs("@SQ\tSN:", out_text) == EOF ||
                kputs(merged_hdr->target_name[i + old_n_targets], out_text) == EOF ||
                kputs("\tLN:", out_text) == EOF ||
                kputuw(merged_hdr->target_len[i + old_n_targets], out_text) == EOF ||
                kputc('\n', out_text) == EOF) goto memfail;
        }
    }

    free(new_sq_matches);
    return 0;

 memfail:
    perror(__func__);
 fail:
    free(new_sq_matches);
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

static klist_t(hdrln) * trans_rg_pg(bool is_rg, bam_hdr_t *translate,
                                    bool merge, khash_t(cset)* known_ids,
                                    khash_t(c2c)* id_map, char *override) {
    hdr_match_t matches[2];
    khiter_t iter;
    const char *text = translate->text;
    const char *rec_type = is_rg ? "@RG" : "@PG";
    klist_t(hdrln) *hdr_lines;

    hdr_lines = kl_init(hdrln);

    // Search through translate's header
    while (hdr_line_match(text, rec_type, "ID", matches) == 0) {
        // matches[0] is the whole @RG/PG line; matches[1] is the ID field value

        kstring_t orig_id = { 0, 0, NULL };        // ID in original header
        kstring_t transformed_id = { 0, 0, NULL }; // ID in output header
        char *map_value;    // Value to store in id_map
        bool id_changed;    // Have we changed the ID?
        bool not_found_in_output; // ID isn't in the output header (yet)

        // Take a copy of the ID as we'll need it for a hash key.
        if (match_to_ks(text, &matches[1], &orig_id)) goto memfail;

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

            if (!id_changed) { // Can just copy
                if (match_to_ks(text, &matches[0], &new_hdr_line)) goto memfail;
            } else { // Substitute new name for original
                if (range_to_ks(text, matches[0].rm_so, matches[1].rm_so,
                                &new_hdr_line)) goto memfail;
                if (ks_to_ks(&transformed_id, &new_hdr_line)) goto memfail;
                if (range_to_ks(text, matches[1].rm_eo, matches[0].rm_eo,
                                &new_hdr_line)) goto memfail;
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

        text += matches[0].rm_eo; // next!
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

static int trans_tbl_init(merged_header_t* merged_hdr, bam_hdr_t* translate,
                          trans_tbl_t* tbl, bool merge_rg, bool merge_pg,
                          bool copy_co, char* rg_override)
{
    klist_t(hdrln) *rg_list = NULL;
    klist_t(hdrln) *pg_list = NULL;

    tbl->n_targets = translate->n_targets;
    tbl->rg_trans = tbl->pg_trans = NULL;
    tbl->tid_trans = (int*)calloc(translate->n_targets, sizeof(int));
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
    if (finish_rg_pg(false, pg_list, tbl->pg_trans, &merged_hdr->out_pg))
        goto fail;

    kl_destroy(hdrln, rg_list); rg_list = NULL;
    kl_destroy(hdrln, pg_list); pg_list = NULL;

    if (copy_co) {
        // Just append @CO headers without translation
        const char *line, *end_pointer;
        for (line = translate->text; *line; line = end_pointer + 1) {
            end_pointer = strchr(line, '\n');
            if (strncmp(line, "@CO", 3) == 0) {
                if (end_pointer) {
                    if (kputsn(line, end_pointer - line + 1, &merged_hdr->out_co) == EOF)
                        goto memfail;
                } else { // Last line with no trailing '\n'
                    if (kputs(line, &merged_hdr->out_co) == EOF) goto memfail;
                    if (kputc('\n', &merged_hdr->out_co) == EOF) goto memfail;
                }
            }
            if (end_pointer == NULL) break;
        }
    }

    return 0;

 memfail:
    perror(__func__);
 fail:
    trans_tbl_destroy(tbl);
    if (rg_list) kl_destroy(hdrln, rg_list);
    if (pg_list) kl_destroy(hdrln, pg_list);
    return -1;
}

static inline void move_kstr_to_text(char **text, kstring_t *ks) {
    memcpy(*text, ks_str(ks), ks_len(ks));
    *text += ks_len(ks);
    **text = '\0';
    free(ks_release(ks));
}

/*
 * Populate a bam_hdr_t struct from data in a merged_header_t.
 */

static bam_hdr_t * finish_merged_header(merged_header_t *merged_hdr) {
    size_t     txt_sz;
    char      *text;
    bam_hdr_t *hdr;

    // Check output text size
    txt_sz = (ks_len(&merged_hdr->out_hd)
              + ks_len(&merged_hdr->out_sq)
              + ks_len(&merged_hdr->out_rg)
              + ks_len(&merged_hdr->out_pg)
              + ks_len(&merged_hdr->out_co));
    if (txt_sz >= INT32_MAX) {
        fprintf(stderr, "[%s] Output header text too long\n", __func__);
        return NULL;
    }

    // Allocate new header
    hdr = bam_hdr_init();
    if (hdr == NULL) goto memfail;

    // Transfer targets arrays to new header
    hdr->n_targets = merged_hdr->n_targets;
    if (hdr->n_targets > 0) {
        // Try to shrink targets arrays to correct size
        hdr->target_name = realloc(merged_hdr->target_name,
                                   hdr->n_targets * sizeof(char*));
        if (!hdr->target_name) hdr->target_name = merged_hdr->target_name;

        hdr->target_len = realloc(merged_hdr->target_len,
                                  hdr->n_targets * sizeof(uint32_t));
        if (!hdr->target_len) hdr->target_len = merged_hdr->target_len;

        // These have either been freed by realloc() or, in the unlikely
        // event that failed, have had their ownership transferred to hdr
        merged_hdr->target_name = NULL;
        merged_hdr->target_len  = NULL;
    }
    else {
        hdr->target_name = NULL;
        hdr->target_len  = NULL;
    }

    // Allocate text
    text = hdr->text = malloc(txt_sz + 1);
    if (!text) goto memfail;

    // Put header text in order @HD, @SQ, @RG, @PG, @CO
    move_kstr_to_text(&text, &merged_hdr->out_hd);
    move_kstr_to_text(&text, &merged_hdr->out_sq);
    move_kstr_to_text(&text, &merged_hdr->out_rg);
    move_kstr_to_text(&text, &merged_hdr->out_pg);
    move_kstr_to_text(&text, &merged_hdr->out_co);
    hdr->l_text = txt_sz;

    return hdr;

 memfail:
    perror(__func__);
    bam_hdr_destroy(hdr);
    return NULL;
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
    free(ks_release(&merged_hdr->out_hd));
    free(ks_release(&merged_hdr->out_sq));
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
  @discussion Padding information may NOT correctly maintained. This
  function is NOT thread safe.
 */
int bam_merge_core2(int by_qname, char* sort_tag, const char *out, const char *mode,
                    const char *headers, int n, char * const *fn, int flag,
                    const char *reg, int n_threads, const char *cmd,
                    const htsFormat *in_fmt, const htsFormat *out_fmt)
{
    samFile *fpout, **fp = NULL;
    heap1_t *heap = NULL;
    bam_hdr_t *hout = NULL;
    bam_hdr_t *hin  = NULL;
    int i, j, *RG_len = NULL;
    uint64_t idx = 0;
    char **RG = NULL;
    hts_itr_t **iter = NULL;
    bam_hdr_t **hdr = NULL;
    trans_tbl_t *translation_tbl = NULL;
    int *rtrans = NULL;
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
        g_sort_tag[1] = sort_tag[1];
    }

    fp = (samFile**)calloc(n, sizeof(samFile*));
    if (!fp) goto mem_fail;
    heap = (heap1_t*)calloc(n, sizeof(heap1_t));
    if (!heap) goto mem_fail;
    iter = (hts_itr_t**)calloc(n, sizeof(hts_itr_t*));
    if (!iter) goto mem_fail;
    hdr = (bam_hdr_t**)calloc(n, sizeof(bam_hdr_t*));
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
        bam_hdr_t *hin;
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
        else { bam_hdr_destroy(hin); hdr[i] = NULL; }

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
    hout = finish_merged_header(merged_hdr);
    if (!hout) return -1;  // FIXME: memory leak

    // If we're only merging a specified region move our iters to start at that point
    if (reg) {
        int tid, beg, end;
        const char *name_lim;

        rtrans = rtrans_build(n, hout->n_targets, translation_tbl);
        if (!rtrans) goto mem_fail;

        name_lim = hts_parse_reg(reg, &beg, &end);
        if (name_lim) {
            char *name = malloc(name_lim - reg + 1);
            if (!name) goto mem_fail;
            memcpy(name, reg, name_lim - reg);
            name[name_lim - reg] = '\0';
            tid = bam_name2id(hout, name);
            free(name);
        }
        else {
            // not parsable as a region, but possibly a sequence named "foo:a"
            tid = bam_name2id(hout, reg);
            beg = 0;
            end = INT_MAX;
        }
        if (tid < 0) {
            if (name_lim) fprintf(stderr, "[%s] Region \"%s\" specifies an unknown reference name\n", __func__, reg);
            else fprintf(stderr, "[%s] Badly formatted region: \"%s\"\n", __func__, reg);
            goto fail;
        }
        for (i = 0; i < n; ++i) {
            hts_idx_t *idx = sam_index_load(fp[i], fn[i]);
            // (rtrans[i*n+tid]) Look up what hout tid translates to in input tid space
            int mapped_tid = rtrans[i*hout->n_targets+tid];
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
                            "{%s, %d, %d, %d}\n",
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
            h->pos = ((uint64_t)h->entry.bam_record->core.tid<<32) | (uint32_t)((int32_t)h->entry.bam_record->core.pos+1);
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
    if (sam_hdr_write(fpout, hout) != 0) {
        print_error_errno(cmd, "failed to write header to \"%s\"", out);
        sam_close(fpout);
        return -1;
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
            return -1;
        }
        if ((j = (iter[heap->i]? sam_itr_next(fp[heap->i], iter[heap->i], b) : sam_read1(fp[heap->i], hdr[heap->i], b))) >= 0) {
            bam_translate(b, translation_tbl + heap->i);
            heap->pos = ((uint64_t)b->core.tid<<32) | (uint32_t)((int)b->core.pos+1);
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

    // Clean up and close
    if (flag & MERGE_RG) {
        for (i = 0; i != n; ++i) free(RG[i]);
        free(RG_len);
    }
    for (i = 0; i < n; ++i) {
        trans_tbl_destroy(translation_tbl + i);
        hts_itr_destroy(iter[i]);
        bam_hdr_destroy(hdr[i]);
        sam_close(fp[i]);
    }
    bam_hdr_destroy(hin);
    bam_hdr_destroy(hout);
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
        if (hdr && hdr[i]) bam_hdr_destroy(hdr[i]);
        if (fp && fp[i]) sam_close(fp[i]);
        if (heap && heap[i].entry.bam_record) bam_destroy1(heap[i].entry.bam_record);
    }
    if (hout) bam_hdr_destroy(hout);
    free(RG);
    free(translation_tbl);
    free(hdr);
    free(iter);
    free(heap);
    free(fp);
    free(rtrans);
    return -1;
}

// Unused here but may be used by legacy samtools-using third-party code
int bam_merge_core(int by_qname, const char *out, const char *headers, int n, char * const *fn, int flag, const char *reg)
{
    char mode[12];
    strcpy(mode, "wb");
    if (flag & MERGE_UNCOMP) strcat(mode, "0");
    else if (flag & MERGE_LEVEL1) strcat(mode, "1");
    return bam_merge_core2(by_qname, NULL, out, mode, headers, n, fn, flag, reg, 0, "merge", NULL, NULL);
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
"  -b FILE    List of input BAM filenames, one per line [null]\n");
    sam_global_opt_help(to, "-.O..@");
}

int bam_merge(int argc, char *argv[])
{
    int c, is_by_qname = 0, flag = 0, ret = 0, level = -1;
    char *fn_headers = NULL, *reg = NULL, mode[12];
    char *sort_tag = NULL;
    long random_seed = (long)time(NULL);
    char** fn = NULL;
    int fn_size = 0;

    sam_global_args ga = SAM_GLOBAL_ARGS_INIT;
    static const struct option lopts[] = {
        SAM_OPT_GLOBAL_OPTIONS('-', 0, 'O', 0, 0, '@'),
        { "threads", required_argument, NULL, '@' },
        { NULL, 0, NULL, 0 }
    };

    if (argc == 1) {
        merge_usage(stdout);
        return 0;
    }

    while ((c = getopt_long(argc, argv, "h:nru1R:f@:l:cps:b:O:t:", lopts, NULL)) >= 0) {
        switch (c) {
        case 'r': flag |= MERGE_RG; break;
        case 'f': flag |= MERGE_FORCE; break;
        case 'h': fn_headers = strdup(optarg); break;
        case 'n': is_by_qname = 1; break;
        case 't': sort_tag = strdup(optarg); break;
        case '1': flag |= MERGE_LEVEL1; level = 1; break;
        case 'u': flag |= MERGE_UNCOMP; level = 0; break;
        case 'R': reg = strdup(optarg); break;
        case 'l': level = atoi(optarg); break;
        case 'c': flag |= MERGE_COMBINE_RG; break;
        case 'p': flag |= MERGE_COMBINE_PG; break;
        case 's': random_seed = atol(optarg); break;
        case 'b': {
            // load the list of files to read
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

    srand48(random_seed);
    if (!(flag & MERGE_FORCE) && strcmp(argv[optind], "-")) {
        FILE *fp = fopen(argv[optind], "rb");
        if (fp != NULL) {
            fclose(fp);
            fprintf(stderr, "[%s] File '%s' exists. Please apply '-f' to overwrite. Abort.\n", __func__, argv[optind]);
            return 1;
        }
    }

    int nargcfiles = argc - (optind+1);
    if (nargcfiles > 0) {
        // Add argc files to end of array
        fn = realloc(fn, (fn_size+nargcfiles) * sizeof(char*));
        if (fn == NULL) { ret = 1; goto end; }
        memcpy(fn+fn_size, argv + (optind+1), nargcfiles * sizeof(char*));
    }
    if (fn_size+nargcfiles < 1) {
        print_error("merge", "You must specify at least one (and usually two or more) input files");
        merge_usage(stderr);
        return 1;
    }
    strcpy(mode, "wb");
    sam_open_mode(mode+1, argv[optind], NULL);
    if (level >= 0) sprintf(strchr(mode, '\0'), "%d", level < 9? level : 9);
    if (bam_merge_core2(is_by_qname, sort_tag, argv[optind], mode, fn_headers,
                        fn_size+nargcfiles, fn, flag, reg, ga.nthreads,
                        "merge", &ga.in, &ga.out) < 0)
        ret = 1;

end:
    if (fn_size > 0) {
        int i;
        for (i=0; i<fn_size; i++) free(fn[i]);
    }
    free(fn);
    free(reg);
    free(fn_headers);
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
                                bam1_tag *buf, uint64_t *idx, bam_hdr_t *hout) {
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
        heap->pos = (((uint64_t)heap->entry.bam_record->core.tid<<32)
                     | (uint32_t)((int32_t)heap->entry.bam_record->core.pos+1));
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
                            const char *mode, bam_hdr_t *hout,
                            int n, char * const *fn, int num_in_mem,
                            buf_region *in_mem, bam1_tag *buf, int n_threads,
                            const char *cmd, const htsFormat *in_fmt,
                            const htsFormat *out_fmt) {
    samFile *fpout = NULL, **fp = NULL;
    heap1_t *heap = NULL;
    uint64_t idx = 0;
    int i, heap_size = n + num_in_mem;

    g_is_by_qname = by_qname;
    if (sort_tag) {
        g_is_by_tag = 1;
        g_sort_tag[0] = sort_tag[0];
        g_sort_tag[1] = sort_tag[1];
    }
    if (n > 0) {
        fp = (samFile**)calloc(n, sizeof(samFile*));
        if (!fp) goto mem_fail;
    }
    heap = (heap1_t*)calloc(heap_size, sizeof(heap1_t));
    if (!heap) goto mem_fail;

    // Open each file, read the header and put the first read into the heap
    for (i = 0; i < heap_size; i++) {
        bam_hdr_t *hin;
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
            bam_hdr_destroy(hin);
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

    if (n_threads > 1) hts_set_threads(fpout, n_threads);

    if (sam_hdr_write(fpout, hout) != 0) {
        print_error_errno(cmd, "failed to write header to \"%s\"", out);
        sam_close(fpout);
        return -1;
    }

    // Now do the merge
    ks_heapmake(heap, heap_size, heap);
    while (heap->pos != HEAP_EMPTY) {
        bam1_t *b = heap->entry.bam_record;
        if (sam_write1(fpout, hout, b) < 0) {
            print_error_errno(cmd, "failed writing to \"%s\"", out);
            sam_close(fpout);
            return -1;
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
        if (heap && heap[i].entry.bam_record) bam_destroy1(heap[i].entry.bam_record);
    }
    free(fp);
    free(heap);
    if (fpout) sam_close(fpout);
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
        pa = (uint64_t)a.bam_record->core.tid<<32|(a.bam_record->core.pos+1);
        pb = (uint64_t)b.bam_record->core.tid<<32|(b.bam_record->core.pos+1);

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
    const bam_hdr_t *h;
    int index;
    int error;
    int no_save;
} worker_t;

// Returns 0 for success
//        -1 for failure
static int write_buffer(const char *fn, const char *mode, size_t l, bam1_tag *buf, const bam_hdr_t *h, int n_threads, const htsFormat *fmt)
{
    size_t i;
    samFile* fp;
    fp = sam_open_format(fn, mode, fmt);
    if (fp == NULL) return -1;
    if (sam_hdr_write(fp, h) != 0) goto fail;
    if (n_threads > 1) hts_set_threads(fp, n_threads);
    for (i = 0; i < l; ++i) {
        if (sam_write1(fp, h, buf[i].bam_record) < 0) goto fail;
    }
    if (sam_close(fp) < 0) return -1;
    return 0;
 fail:
    sam_close(fp);
    return -1;
}

#define NUMBASE 256
#define STEP 8

static int ks_radixsort(size_t n, bam1_tag *buf, const bam_hdr_t *h)
{
    int curr = 0, ret = -1;
    ssize_t i;
    bam1_tag *buf_ar2[2], *bam_a, *bam_b;
    uint64_t max_pos = 0, max_digit = 0, shift = 0;

    for (i = 0; i < n; i++) {
        bam1_t *b = buf[i].bam_record;
        int32_t tid = b->core.tid == -1 ? h->n_targets : b->core.tid;
        buf[i].u.pos = (uint64_t)tid<<32 | (b->core.pos+1)<<1 | bam_is_rev(b);
        if (max_pos < buf[i].u.pos)
            max_pos = buf[i].u.pos;
    }

    while (max_pos) {
        ++max_digit;
        max_pos = max_pos >> 1;
    }

    buf_ar2[0] = buf;
    buf_ar2[1] = (bam1_tag *)malloc(sizeof(bam1_tag) * n);
    if (buf_ar2[1] == NULL) {
        print_error("sort", "couldn't allocate memory for temporary buf");
        goto err;
    }

    while (shift < max_digit){
        size_t remainders[NUMBASE] = { 0 };
        bam_a = buf_ar2[curr]; bam_b = buf_ar2[1-curr];
        for (i = 0; i < n; ++i)
            remainders[(bam_a[i].u.pos >> shift) % NUMBASE]++;
        for (i = 1; i < NUMBASE; ++i)
            remainders[i] += remainders[i - 1];
        for (i = n - 1; i >= 0; i--) {
            size_t j = --remainders[(bam_a[i].u.pos >> shift) % NUMBASE];
            bam_b[j] = bam_a[i];
        }
        shift += STEP;
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

        if (write_buffer(name, "wcx1", w->buf_len, w->buf, w->h, 0, &fmt) < 0)
            w->error = errno;
    } else {
        if (write_buffer(name, "wbx1", w->buf_len, w->buf, w->h, 0, NULL) < 0)
            w->error = errno;
    }

    free(name);
    return 0;
}

static int sort_blocks(int n_files, size_t k, bam1_tag *buf, const char *prefix,
                       const bam_hdr_t *h, int n_threads, buf_region *in_mem)
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
  @return 0 for successful sorting, negative on errors

  @discussion It may create multiple temporary subalignment files
  and then merge them by calling bam_merge_simple(). This function is
  NOT thread safe.
 */
int bam_sort_core_ext(int is_by_qname, char* sort_by_tag, const char *fn, const char *prefix,
                      const char *fnout, const char *modeout,
                      size_t _max_mem, int n_threads,
                      const htsFormat *in_fmt, const htsFormat *out_fmt)
{
    int ret = -1, res, i, n_files = 0;
    size_t max_k, k, max_mem, bam_mem_offset;
    bam_hdr_t *header = NULL;
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
        strncpy(g_sort_tag, sort_by_tag, 2);
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

    if (sam_hdr_change_HD(header, "SO", new_so) != 0) {
        print_error("sort",
                    "failed to change sort order header to '%s'\n", new_so);
        goto err;
    }
    if (sam_hdr_change_HD(header, "GO", NULL) != 0) {
        print_error("sort",
                    "failed to delete group order header\n");
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
        if (write_buffer(fnout, modeout, k, buf, header, n_threads, out_fmt) != 0) {
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
                             n_threads, "sort", in_fmt, out_fmt) < 0) {
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
    bam_hdr_destroy(header);
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
    ret = bam_sort_core_ext(is_by_qname, NULL, fn, prefix, fnout, "wb", max_mem, 0, NULL, NULL);
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
"  -T PREFIX  Write temporary files to PREFIX.nnnn.bam\n");
    sam_global_opt_help(fp, "-.O..@");
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
    int c, nargs, is_by_qname = 0, ret, o_seen = 0, level = -1;
    char* sort_tag = NULL;
    char *fnout = "-", modeout[12];
    kstring_t tmpprefix = { 0, 0, NULL };
    struct stat st;
    sam_global_args ga = SAM_GLOBAL_ARGS_INIT;

    static const struct option lopts[] = {
        SAM_OPT_GLOBAL_OPTIONS('-', 0, 'O', 0, 0, '@'),
        { "threads", required_argument, NULL, '@' },
        { NULL, 0, NULL, 0 }
    };

    while ((c = getopt_long(argc, argv, "l:m:no:O:T:@:t:", lopts, NULL)) >= 0) {
        switch (c) {
        case 'o': fnout = optarg; o_seen = 1; break;
        case 'n': is_by_qname = 1; break;
        case 't': sort_tag = strdup(optarg); break;
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
                            &ga.in, &ga.out);
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
    sam_global_args_free(&ga);

    return ret;
}
