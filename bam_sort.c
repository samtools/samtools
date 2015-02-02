/*  bam_sort.c -- sorting and merging.

    Copyright (C) 2008-2014 Genome Research Ltd.
    Portions copyright (C) 2009-2012 Broad Institute.

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

#include <stdbool.h>
#include <stdlib.h>
#include <ctype.h>
#include <errno.h>
#include <stdio.h>
#include <string.h>
#include <regex.h>
#include <time.h>
#include <unistd.h>
#include "htslib/ksort.h"
#include "htslib/khash.h"
#include "htslib/klist.h"
#include "htslib/kstring.h"
#include "htslib/sam.h"

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
KHASH_MAP_INIT_STR(c2i, int)

#define __free_char(p)
KLIST_INIT(hdrln, char*, __free_char)

static int g_is_by_qname = 0;

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
    uint64_t pos, idx;
    bam1_t *b;
} heap1_t;

#define __pos_cmp(a, b) ((a).pos > (b).pos || ((a).pos == (b).pos && ((a).i > (b).i || ((a).i == (b).i && (a).idx > (b).idx))))

// Function to compare reads in the heap and determine which one is < the other
static inline int heap_lt(const heap1_t a, const heap1_t b)
{
    if (g_is_by_qname) {
        int t;
        if (a.b == NULL || b.b == NULL) return a.b == NULL? 1 : 0;
        t = strnum_cmp(bam_get_qname(a.b), bam_get_qname(b.b));
        return (t > 0 || (t == 0 && (a.b->core.flag&0xc0) > (b.b->core.flag&0xc0)));
    } else return __pos_cmp(a, b);
}

KSORT_INIT(heap, heap1_t, heap_lt)

typedef struct trans_tbl {
    int32_t n_targets;
    int* tid_trans;
    kh_c2c_t* rg_trans;
    kh_c2c_t* pg_trans;
    bool lost_coord_sort;
} trans_tbl_t;

static void trans_tbl_destroy(trans_tbl_t *tbl) {
    free(tbl->tid_trans);
    khiter_t iter;
    for (iter = kh_begin(tbl->rg_trans); iter != kh_end(tbl->rg_trans); ++iter) {
        if (kh_exist(tbl->rg_trans, iter)) {
            free(kh_value(tbl->rg_trans, iter));
            free(kh_key(tbl->rg_trans, iter));
        }
    }
    for (iter = kh_begin(tbl->pg_trans); iter != kh_end(tbl->pg_trans); ++iter) {
        if (kh_exist(tbl->pg_trans, iter)) {
            free(kh_value(tbl->pg_trans, iter));
            free(kh_key(tbl->pg_trans, iter));
        }
    }

    kh_destroy(c2c,tbl->rg_trans);
    kh_destroy(c2c,tbl->pg_trans);
}

// Takes in existing header and rewrites it in the usual order HD, SQ, RG, PG CO, other
static void pretty_header(char** text_in_out, int32_t text_len)
{
    char* output, *output_pointer;
    output = output_pointer = (char*)calloc(1,text_len+1);
    output[text_len] = '\0';

    // Read @HD and write
    regex_t hd_regex, sq_regex, pg_regex, rg_regex, co_regex, other_regex;
    regmatch_t matches[1];
    if (regcomp( &hd_regex, "^@HD.*$", REG_EXTENDED|REG_NEWLINE ))
        abort();
    if (regexec( &hd_regex, *text_in_out, 1, &matches[0], 0 ) == 0) {
        size_t match_size = matches[0].rm_eo - matches[0].rm_so;
        memcpy(output_pointer, *text_in_out+matches[0].rm_so, match_size);
        output_pointer[match_size] = '\n';
        output_pointer += match_size + 1;
    }
    regfree(&hd_regex);

    // Read @SQ's and write
    if (regcomp( &sq_regex, "^@SQ.*$", REG_EXTENDED|REG_NEWLINE )) abort();
    char* sq_pointer = *text_in_out;
    while (*text_in_out+text_len > sq_pointer && regexec( &sq_regex, sq_pointer, 1, &matches[0], 0) == 0) {
        size_t match_size = matches[0].rm_eo - matches[0].rm_so;
        memcpy(output_pointer, sq_pointer+matches[0].rm_so, match_size);
        output_pointer[match_size] = '\n';
        output_pointer += match_size + 1;
        sq_pointer += matches[0].rm_eo + 1;
    }
    regfree(&sq_regex);

    // Read @RG's and write
    if (regcomp( &rg_regex, "^@RG.*$", REG_EXTENDED|REG_NEWLINE )) abort();
    char* rg_pointer = *text_in_out;
    while (*text_in_out+text_len > rg_pointer && regexec( &rg_regex, rg_pointer, 1, &matches[0], 0) == 0) {
        size_t match_size = matches[0].rm_eo - matches[0].rm_so;
        memcpy(output_pointer, rg_pointer+matches[0].rm_so, match_size);
        output_pointer[match_size] = '\n';
        output_pointer += match_size + 1;
        rg_pointer += matches[0].rm_eo + 1;
    }
    regfree(&rg_regex);

    // Read @PG's and write
    if (regcomp( &pg_regex, "^@PG.*$", REG_EXTENDED|REG_NEWLINE )) abort();
    char* pg_pointer = *text_in_out;
    while (*text_in_out+text_len > pg_pointer && regexec( &pg_regex, pg_pointer, 1, &matches[0], 0) == 0) {
        size_t match_size = matches[0].rm_eo - matches[0].rm_so;
        memcpy(output_pointer, pg_pointer+matches[0].rm_so, match_size);
        output_pointer[match_size] = '\n';
        output_pointer += match_size + 1;
        pg_pointer += matches[0].rm_eo + 1;
    }
    regfree(&pg_regex);

    // Read @CO's and write
    if (regcomp( &co_regex, "^@CO.*$", REG_EXTENDED|REG_NEWLINE )) abort();
    char* co_pointer = *text_in_out;
    while (*text_in_out+text_len > co_pointer && regexec( &co_regex, co_pointer, 1, &matches[0], 0) == 0) {
        size_t match_size = matches[0].rm_eo - matches[0].rm_so;
        memcpy(output_pointer, co_pointer+matches[0].rm_so, match_size);
        output_pointer[match_size] = '\n';
        output_pointer += match_size + 1;
        co_pointer += matches[0].rm_eo + 1;
    }
    regfree(&co_regex);

    // Read any other not HD,SQ,RG,PG,CO tags and write
    if (regcomp( &other_regex, "^@([^HSCPR]|H[^D]|S[^Q]|[PR][^G]|C[^O]).*$", REG_EXTENDED|REG_NEWLINE )) abort();
    char* other_pointer = *text_in_out;
    while (*text_in_out+text_len > other_pointer && regexec( &other_regex, other_pointer, 1, &matches[0], 0) == 0) {
        size_t match_size = matches[0].rm_eo - matches[0].rm_so;
        memcpy(output_pointer, other_pointer+matches[0].rm_so, match_size);
        output_pointer[match_size] = '\n';
        output_pointer += match_size + 1;
        other_pointer += matches[0].rm_eo + 1;
    }
    regfree(&other_regex);

    // Safety check, make sure we copied it all, if we didn't something is wrong with the header
    if ( output+text_len != output_pointer ) {
        fprintf(stderr, "[pretty_header] invalid header\n");
        exit(1);
    }
    free(*text_in_out);
    *text_in_out = output;
}

static void trans_tbl_init(bam_hdr_t* out, bam_hdr_t* translate, trans_tbl_t* tbl, bool merge_rg, bool merge_pg)
{
    tbl->n_targets = translate->n_targets;
    tbl->tid_trans = (int*)calloc(translate->n_targets, sizeof(int));
    tbl->rg_trans = kh_init(c2c);
    tbl->pg_trans = kh_init(c2c);
    if (!tbl->tid_trans || !tbl->rg_trans || !tbl->pg_trans) { perror("out of memory"); exit(-1); }

    int32_t out_len = out->l_text;
    while (out_len > 0 && out->text[out_len-1] == '\n') {--out_len; } // strip trailing \n's
    kstring_t out_text = { 0, 0, NULL };
    kputsn(out->text, out_len, &out_text);

    int i, min_tid = -1;
    tbl->lost_coord_sort = false;

    khash_t(c2i) *out_tid = kh_init(c2i);
    for (i = 0; i < out->n_targets; ++i) {
        int ret;
        khiter_t iter = kh_put(c2i, out_tid, out->target_name[i], &ret);
        if (ret <= 0) abort();
        kh_value(out_tid, iter) = i;
    }

    for (i = 0; i < translate->n_targets; ++i) {
        khiter_t iter = kh_get(c2i, out_tid, translate->target_name[i]);

        if (iter == kh_end(out_tid)) { // Append missing entries to out
            tbl->tid_trans[i] = out->n_targets++;
            out->target_name = (char**)realloc(out->target_name, sizeof(char*)*out->n_targets);
            out->target_name[out->n_targets-1] = strdup(translate->target_name[i]);
            out->target_len = (uint32_t*)realloc(out->target_len, sizeof(uint32_t)*out->n_targets);
            out->target_len[out->n_targets-1] = translate->target_len[i];
            // grep line with regex '^@SQ.*\tSN:%s(\t.*$|$)', translate->target_name[i]
            // from translate->text
            regex_t sq_id;
            regmatch_t* matches = (regmatch_t*)calloc(2, sizeof(regmatch_t));
            if (matches == NULL) { perror("out of memory"); exit(-1); }
            kstring_t seq_regex = { 0, 0, NULL };
            ksprintf(&seq_regex, "^@SQ.*\tSN:%s(\t.*$|$)", translate->target_name[i]);
            regcomp(&sq_id, seq_regex.s, REG_EXTENDED|REG_NEWLINE);
            free(seq_regex.s);
            if (regexec(&sq_id, translate->text, 1, matches, 0) != 0)
            {
                fprintf(stderr, "[trans_tbl_init] @SQ SN (%s) found in binary header but not text header.\n",translate->target_name[i]);
                exit(1);
            }
            regfree(&sq_id);

            // Produce our output line and append it to out_text
            kputc('\n', &out_text);
            kputsn(translate->text+matches[0].rm_so, matches[0].rm_eo-matches[0].rm_so, &out_text);

            free(matches);
        } else {
            tbl->tid_trans[i] = kh_value(out_tid, iter);
        }
        if (tbl->tid_trans[i] > min_tid) {
            min_tid = tbl->tid_trans[i];
        } else {
            tbl->lost_coord_sort = true;
        }
    }
    kh_destroy(c2i, out_tid);

    // grep @RG id's
    regex_t rg_id;
    regmatch_t* matches = (regmatch_t*)calloc(2, sizeof(regmatch_t));
    if (matches == NULL) { perror("out of memory"); exit(-1); }
    regcomp(&rg_id, "^@RG.*\tID:([!-)+-<>-~][ !-~]*)(\t.*$|$)", REG_EXTENDED|REG_NEWLINE);
    char* text = translate->text;
    klist_t(hdrln) *rg_list = kl_init(hdrln);
    while(1) { //   foreach rg id in translate's header
        if (regexec(&rg_id, text, 2, matches, 0) != 0) break;
        // matches[0] is the whole @RG line; matches[1] is the ID field value
        kstring_t match_id = { 0, 0, NULL };
        kputsn(text+matches[1].rm_so, matches[1].rm_eo-matches[1].rm_so, &match_id);

        // is our matched ID in our output list already
        regex_t rg_id_search;
        kstring_t rg_regex = { 0, 0, NULL };
        ksprintf(&rg_regex, "^@RG.*\tID:%s(\t.*$|$)", match_id.s);
        regcomp(&rg_id_search, rg_regex.s, REG_EXTENDED|REG_NEWLINE|REG_NOSUB);
        free(rg_regex.s);
        kstring_t transformed_id = { 0, 0, NULL };
        bool transformed_equals_match;
        if (regexec(&rg_id_search, out->text, 0, NULL, 0) != 0  || merge_rg) {
            // Not in there so can add it as 1-1 mapping
            kputs(match_id.s, &transformed_id);
            transformed_equals_match = true;
        } else {
            // It's in there so we need to transform it by appending random number to id
            ksprintf(&transformed_id, "%s-%0lX", match_id.s, lrand48());
            transformed_equals_match = false;
        }
        regfree(&rg_id_search);

        // Insert it into our translation map
        int in_there = 0;
        khiter_t iter = kh_put(c2c, tbl->rg_trans, ks_release(&match_id), &in_there);
        char *transformed_id_s = ks_release(&transformed_id);
        kh_value(tbl->rg_trans,iter) = transformed_id_s;
        // take matched line and replace ID with transformed_id
        kstring_t transformed_line = { 0, 0, NULL };
        if (transformed_equals_match) {
            kputsn(text+matches[0].rm_so, matches[0].rm_eo-matches[0].rm_so, &transformed_line);
        } else {
            kputsn(text+matches[0].rm_so, matches[1].rm_so-matches[0].rm_so, &transformed_line);
            kputs(transformed_id_s, &transformed_line);
            kputsn(text+matches[1].rm_eo, matches[0].rm_eo-matches[1].rm_eo, &transformed_line);
        }

        if (!(transformed_equals_match && merge_rg)) {
            // append line to linked list for PG processing
            char** ln = kl_pushp(hdrln, rg_list);
            *ln = ks_release(&transformed_line);  // Give away to linked list
        }
        else free(transformed_line.s);

        text += matches[0].rm_eo; // next!
    }
    regfree(&rg_id);

    // Do same for PG id's
    regex_t pg_id;
    regcomp(&pg_id, "^@PG.*\tID:([!-)+-<>-~][ !-~]*)(\t.*$|$)", REG_EXTENDED|REG_NEWLINE);
    text = translate->text;
    klist_t(hdrln) *pg_list = kl_init(hdrln);
    while(1) { //   foreach pg id in translate's header
        if (regexec(&pg_id, text, 2, matches, 0) != 0) break;
        kstring_t match_id = { 0, 0, NULL };
        kputsn(text+matches[1].rm_so, matches[1].rm_eo-matches[1].rm_so, &match_id);

        // is our matched ID in our output list already
        regex_t pg_id_search;
        kstring_t pg_regex = { 0, 0, NULL };
        ksprintf(&pg_regex, "^@PG.*\tID:%s(\t.*$|$)", match_id.s);
        regcomp(&pg_id_search, pg_regex.s, REG_EXTENDED|REG_NEWLINE|REG_NOSUB);
        free(pg_regex.s);
        kstring_t transformed_id = { 0, 0, NULL };
        bool transformed_equals_match;
        if (regexec(&pg_id_search, out->text, 0, NULL, 0) != 0 || merge_pg) {
            // Not in there so can add it as 1-1 mapping
            kputs(match_id.s, &transformed_id);
            transformed_equals_match = true;
        } else {
            // It's in there so we need to transform it by appending random number to id
            ksprintf(&transformed_id, "%s-%0lX", match_id.s, lrand48());
            transformed_equals_match = false;
        }
        regfree(&pg_id_search);

        // Insert it into our translation map
        int in_there = 0;
        khiter_t iter = kh_put(c2c, tbl->pg_trans, ks_release(&match_id), &in_there);
        char *transformed_id_s = ks_release(&transformed_id);
        kh_value(tbl->pg_trans,iter) = transformed_id_s;
        // take matched line and replace ID with transformed_id
        kstring_t transformed_line = { 0, 0, NULL };
        if (transformed_equals_match) {
            kputsn(text+matches[0].rm_so, matches[0].rm_eo-matches[0].rm_so, &transformed_line);
        } else {
            kputsn(text+matches[0].rm_so, matches[1].rm_so-matches[0].rm_so, &transformed_line);
            kputs(transformed_id_s, &transformed_line);
            kputsn(text+matches[1].rm_eo, matches[0].rm_eo-matches[1].rm_eo, &transformed_line);
        }

        if (!(transformed_equals_match && merge_pg)) {
            // append line to linked list for PP processing
            char** ln = kl_pushp(hdrln, pg_list);
            *ln = ks_release(&transformed_line);  // Give away to linked list
        }
        else free(transformed_line.s);
        text += matches[0].rm_eo; // next!
    }
    regfree(&pg_id);
    // need to translate PP's on the fly in second pass because they may not be in correct order and need complete tbl->pg_trans to do this
    // for each line {
    // with ID replaced with tranformed_id and PP's transformed using the translation table
    // }
    regex_t pg_pp;
    regcomp(&pg_pp, "^@PG.*\tPP:([!-)+-<>-~][!-~]*)(\t.*$|$)", REG_EXTENDED|REG_NEWLINE);
    kliter_t(hdrln) *iter = kl_begin(pg_list);
    while (iter != kl_end(pg_list)) {
        char* data = kl_val(iter);

        kstring_t transformed_line = { 0, 0, NULL };
        // Find PP tag
        if (regexec(&pg_pp, data, 2, matches, 0) == 0) {
            // Lookup in hash table
            kstring_t pp_id = { 0, 0, NULL };
            kputsn(data+matches[1].rm_so, matches[1].rm_eo-matches[1].rm_so, &pp_id);

            khiter_t k = kh_get(c2c, tbl->pg_trans, pp_id.s);
            free(pp_id.s);
            char* transformed_id = kh_value(tbl->pg_trans,k);
            // Replace
            kputsn(data, matches[1].rm_so-matches[0].rm_so, &transformed_line);
            kputs(transformed_id, &transformed_line);
            kputsn(data+matches[1].rm_eo, matches[0].rm_eo-matches[1].rm_eo, &transformed_line);
        } else { kputs(data, &transformed_line); }
        // Produce our output line and append it to out_text
        kputc('\n', &out_text);
        kputsn(transformed_line.s, transformed_line.l, &out_text);

        free(transformed_line.s);
        free(data);
        iter = kl_next(iter);
    }
    regfree(&pg_pp);

    // Need to also translate @RG PG's on the fly too
    regex_t rg_pg;
    regcomp(&rg_pg, "^@RG.*\tPG:([!-)+-<>-~][!-~]*)(\t.*$|$)", REG_EXTENDED|REG_NEWLINE);
    kliter_t(hdrln) *rg_iter = kl_begin(rg_list);
    while (rg_iter != kl_end(rg_list)) {
        char* data = kl_val(rg_iter);

        kstring_t transformed_line = { 0, 0, NULL };
        // Find PG tag
        if (regexec(&rg_pg, data, 2, matches, 0) == 0) {
            // Lookup in hash table
            kstring_t pg_id = { 0, 0, NULL };
            kputsn(data+matches[1].rm_so, matches[1].rm_eo-matches[1].rm_so, &pg_id);

            khiter_t k = kh_get(c2c, tbl->pg_trans, pg_id.s);
            free(pg_id.s);
            char* transformed_id = kh_value(tbl->pg_trans,k);
            // Replace
            kputsn(data, matches[1].rm_so-matches[0].rm_so, &transformed_line);
            kputs(transformed_id, &transformed_line);
            kputsn(data+matches[1].rm_eo, matches[0].rm_eo-matches[1].rm_eo, &transformed_line);
        } else { kputs(data, &transformed_line); }
        // Produce our output line and append it to out_text
        kputc('\n', &out_text);
        kputsn(transformed_line.s, transformed_line.l, &out_text);

        free(transformed_line.s);
        free(data);
        rg_iter = kl_next(rg_iter);
    }

    regfree(&rg_pg);
    kl_destroy(hdrln,pg_list);
    kl_destroy(hdrln,rg_list);
    free(matches);

    // Add trailing \n and write back to header
    free(out->text);
    kputc('\n', &out_text);
    out->l_text = out_text.l;
    out->text = ks_release(&out_text);
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
            bam_aux_append(b, "RG", 'Z', strlen(translate_rg) + 1, (uint8_t*)translate_rg);
        } else {
            fprintf(stderr, "[bam_translate] RG tag \"%s\" on read \"%s\" encountered with no corresponding entry in header, tag lost\n",decoded_rg, bam_get_qname(b));
            bam_aux_del(b, rg);
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
            bam_aux_append(b, "PG", 'Z', strlen(translate_pg) + 1, (uint8_t*)translate_pg);
        } else {
            fprintf(stderr, "[bam_translate] PG tag \"%s\" on read \"%s\" encountered with no corresponding entry in header, tag lost\n",decoded_pg, bam_get_qname(b));
            bam_aux_del(b, pg);
        }
    }
}

int* rtrans_build(int n, int n_targets, trans_tbl_t* translation_tbl)
{
    // Create reverse translation table for tids
    int* rtrans = (int*)malloc(sizeof(int32_t)*n*n_targets);
    const int32_t NOTID = INT32_MIN;
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
  @param  is_by_qname whether to sort by query name
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
  @discussion Padding information may NOT correctly maintained. This
  function is NOT thread safe.
 */
int bam_merge_core2(int by_qname, const char *out, const char *mode, const char *headers, int n, char * const *fn, int flag, const char *reg, int n_threads)
{
    samFile *fpout, **fp;
    heap1_t *heap;
    bam_hdr_t *hout = NULL;
    int i, j, *RG_len = NULL;
    uint64_t idx = 0;
    char **RG = NULL;
    hts_itr_t **iter = NULL;
    bam_hdr_t **hdr = NULL;
    trans_tbl_t *translation_tbl = NULL;

    // Is there a specified pre-prepared header to use for output?
    if (headers) {
        samFile* fpheaders = sam_open(headers, "r");
        if (fpheaders == NULL) {
            const char *message = strerror(errno);
            fprintf(stderr, "[bam_merge_core] cannot open '%s': %s\n", headers, message);
            return -1;
        }
        hout = sam_hdr_read(fpheaders);
        sam_close(fpheaders);
    }

    g_is_by_qname = by_qname;
    fp = (samFile**)calloc(n, sizeof(samFile*));
    heap = (heap1_t*)calloc(n, sizeof(heap1_t));
    iter = (hts_itr_t**)calloc(n, sizeof(hts_itr_t*));
    hdr = (bam_hdr_t**)calloc(n, sizeof(bam_hdr_t*));
    translation_tbl = (trans_tbl_t*)calloc(n, sizeof(trans_tbl_t));
    // prepare RG tag from file names
    if (flag & MERGE_RG) {
        RG = (char**)calloc(n, sizeof(char*));
        RG_len = (int*)calloc(n, sizeof(int));
        for (i = 0; i != n; ++i) {
            int l = strlen(fn[i]);
            const char *s = fn[i];
            if (l > 4 && strcmp(s + l - 4, ".bam") == 0) l -= 4;
            for (j = l - 1; j >= 0; --j) if (s[j] == '/') break;
            ++j; l -= j;
            RG[i] = (char*)calloc(l + 1, 1);
            RG_len[i] = l;
            strncpy(RG[i], s + j, l);
        }
    }
    // open and read the header from each file
    for (i = 0; i < n; ++i) {
        bam_hdr_t *hin;
        fp[i] = sam_open(fn[i], "r");
        if (fp[i] == NULL) {
            int j;
            fprintf(stderr, "[bam_merge_core] fail to open file %s\n", fn[i]);
            for (j = 0; j < i; ++j) sam_close(fp[j]);
            free(fp); free(heap);
            // FIXME: possible memory leak
            return -1;
        }
        hin = sam_hdr_read(fp[i]);
        if (hout)
            trans_tbl_init(hout, hin, translation_tbl+i, flag & MERGE_COMBINE_RG, flag & MERGE_COMBINE_PG);
        else {
            // As yet, no headers to merge into...
            hout = bam_hdr_dup(hin);
            // ...so no need to translate header into itself
            trans_tbl_init(hout, hin, translation_tbl+i, true, true);
        }

        // TODO sam_itr_next() doesn't yet work for SAM files,
        // so for those keep the headers around for use with sam_read1()
        if (hts_get_format(fp[i])->format == sam) hdr[i] = hin;
        else { bam_hdr_destroy(hin); hdr[i] = NULL; }

        if ((translation_tbl+i)->lost_coord_sort && !by_qname) {
            fprintf(stderr, "[bam_merge_core] Order of targets in file %s caused coordinate sort to be lost\n", fn[i]);
        }
    }

    // Transform the header into standard form
    pretty_header(&hout->text,hout->l_text);

    // If we're only merging a specified region move our iters to start at that point
    if (reg) {
        int* rtrans = rtrans_build(n, hout->n_targets, translation_tbl);

        int tid, beg, end;
        const char *name_lim = hts_parse_reg(reg, &beg, &end);
        char *name = malloc(name_lim - reg + 1);
        memcpy(name, reg, name_lim - reg);
        name[name_lim - reg] = '\0';
        tid = bam_name2id(hout, name);
        free(name);
        if (tid < 0) {
            fprintf(stderr, "[%s] Malformated region string or undefined reference name\n", __func__);
            return -1;
        }
        for (i = 0; i < n; ++i) {
            hts_idx_t *idx = sam_index_load(fp[i], fn[i]);
            // (rtrans[i*n+tid]) Look up what hout tid translates to in input tid space
            int mapped_tid = rtrans[i*hout->n_targets+tid];
            if (mapped_tid != INT32_MIN) {
                iter[i] = sam_itr_queryi(idx, mapped_tid, beg, end);
            } else {
                iter[i] = sam_itr_queryi(idx, HTS_IDX_NONE, 0, 0);
            }
            hts_idx_destroy(idx);
            if (iter[i] == NULL) break;
        }
        free(rtrans);
    } else {
        for (i = 0; i < n; ++i) {
            if (hdr[i] == NULL) {
                iter[i] = sam_itr_queryi(NULL, HTS_IDX_REST, 0, 0);
                if (iter[i] == NULL) break;
            }
            else iter[i] = NULL;
        }
    }

    if (i < n) {
        fprintf(stderr, "[%s] Memory allocation failed\n", __func__);
        return -1;
    }

    // Load the first read from each file into the heap
    for (i = 0; i < n; ++i) {
        heap1_t *h = heap + i;
        h->i = i;
        h->b = bam_init1();
        if ((iter[i]? sam_itr_next(fp[i], iter[i], h->b) : sam_read1(fp[i], hdr[i], h->b)) >= 0) {
            bam_translate(h->b, translation_tbl + i);
            h->pos = ((uint64_t)h->b->core.tid<<32) | (uint32_t)((int32_t)h->b->core.pos+1)<<1 | bam_is_rev(h->b);
            h->idx = idx++;
        }
        else {
            h->pos = HEAP_EMPTY;
            bam_destroy1(h->b);
            h->b = NULL;
        }
    }

    // Open output file and write header
    if ((fpout = sam_open(out, mode)) == 0) {
        fprintf(stderr, "[%s] fail to create the output file.\n", __func__);
        return -1;
    }
    sam_hdr_write(fpout, hout);
    if (!(flag & MERGE_UNCOMP)) hts_set_threads(fpout, n_threads);

    // Begin the actual merge
    ks_heapmake(heap, n, heap);
    while (heap->pos != HEAP_EMPTY) {
        bam1_t *b = heap->b;
        if (flag & MERGE_RG) {
            uint8_t *rg = bam_aux_get(b, "RG");
            if (rg) bam_aux_del(b, rg);
            bam_aux_append(b, "RG", 'Z', RG_len[heap->i] + 1, (uint8_t*)RG[heap->i]);
        }
        sam_write1(fpout, hout, b);
        if ((j = (iter[heap->i]? sam_itr_next(fp[heap->i], iter[heap->i], b) : sam_read1(fp[heap->i], hdr[heap->i], b))) >= 0) {
            bam_translate(b, translation_tbl + heap->i);
            heap->pos = ((uint64_t)b->core.tid<<32) | (uint32_t)((int)b->core.pos+1)<<1 | bam_is_rev(b);
            heap->idx = idx++;
        } else if (j == -1) {
            heap->pos = HEAP_EMPTY;
            bam_destroy1(heap->b);
            heap->b = NULL;
        } else fprintf(stderr, "[bam_merge_core] '%s' is truncated. Continue anyway.\n", fn[heap->i]);
        ks_heapadjust(heap, 0, n, heap);
    }

    // Clean up and close
    if (flag & MERGE_RG) {
        for (i = 0; i != n; ++i) free(RG[i]);
        free(RG); free(RG_len);
    }
    for (i = 0; i < n; ++i) {
        trans_tbl_destroy(translation_tbl + i);
        hts_itr_destroy(iter[i]);
        bam_hdr_destroy(hdr[i]);
        sam_close(fp[i]);
    }
    bam_hdr_destroy(hout);
    sam_close(fpout);
    free(translation_tbl); free(fp); free(heap); free(iter); free(hdr);
    return 0;
}

int bam_merge_core(int by_qname, const char *out, const char *headers, int n, char * const *fn, int flag, const char *reg)
{
    char mode[12];
    strcpy(mode, "wb");
    if (flag & MERGE_UNCOMP) strcat(mode, "0");
    else if (flag & MERGE_LEVEL1) strcat(mode, "1");
    return bam_merge_core2(by_qname, out, mode, headers, n, fn, flag, reg, 0);
}

static void merge_usage(FILE *to)
{
    fprintf(to, "Usage:   samtools merge [-nurlf] [-h inh.sam] [-b <bamlist.fofn>] <out.bam> <in1.bam> <in2.bam> [<in3.bam> ... <inN.bam>]\n\n");
    fprintf(to, "Options: -n       sort by read names\n");
    fprintf(to, "         -r       attach RG tag (inferred from file names)\n");
    fprintf(to, "         -u       uncompressed BAM output\n");
    fprintf(to, "         -f       overwrite the output BAM if exist\n");
    fprintf(to, "         -1       compress level 1\n");
    fprintf(to, "         -l INT   compression level, from 0 to 9 [-1]\n");
    fprintf(to, "         -@ INT   number of BAM compression threads [0]\n");
    fprintf(to, "         -R STR   merge file in the specified region STR [all]\n");
    fprintf(to, "         -h FILE  copy the header in FILE to <out.bam> [in1.bam]\n");
    fprintf(to, "         -c       combine RG tags with colliding IDs rather than amending them\n");
    fprintf(to, "         -p       combine PG tags with colliding IDs rather than amending them\n");
    fprintf(to, "         -s VALUE override random seed\n");
    fprintf(to, "         -b FILE  list of input BAM filenames, one per line [null]\n\n");
}

int bam_merge(int argc, char *argv[])
{
    int c, is_by_qname = 0, flag = 0, ret = 0, n_threads = 0, level = -1;
    char *fn_headers = NULL, *reg = NULL, mode[12];
    long random_seed = (long)time(NULL);
    char** fn = NULL;
    int fn_size = 0;

    if (argc == 1) {
        merge_usage(stdout);
        return 0;
    }

    while ((c = getopt(argc, argv, "h:nru1R:f@:l:cps:b:")) >= 0) {
        switch (c) {
        case 'r': flag |= MERGE_RG; break;
        case 'f': flag |= MERGE_FORCE; break;
        case 'h': fn_headers = strdup(optarg); break;
        case 'n': is_by_qname = 1; break;
        case '1': flag |= MERGE_LEVEL1; level = 1; break;
        case 'u': flag |= MERGE_UNCOMP; level = 0; break;
        case 'R': reg = strdup(optarg); break;
        case 'l': level = atoi(optarg); break;
        case '@': n_threads = atoi(optarg); break;
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
            }
            else {
                fprintf(stderr, "[%s] Invalid file list \"%s\"\n", __func__, optarg);
                ret = 1;
            }
            break;
        }
        }
    }
    if ( argc - optind < 1 ) {
        fprintf(stderr, "You must at least specify the output file.\n");
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
    if (fn_size+nargcfiles < 2) {
        fprintf(stderr, "You must specify at least 2 input files.\n");
        merge_usage(stderr);
        return 1;
    }
    strcpy(mode, "wb");
    if (level >= 0) sprintf(strchr(mode, '\0'), "%d", level < 9? level : 9);
    if (bam_merge_core2(is_by_qname, argv[optind], mode, fn_headers, fn_size+nargcfiles, fn, flag, reg, n_threads) < 0) ret = 1;
end:
    if (fn_size > 0) {
        int i;
        for (i=0; i<fn_size; i++) free(fn[i]);
        free(fn);
    }
    free(reg);
    free(fn_headers);
    return ret;
}

/***************
 * BAM sorting *
 ***************/

#include <pthread.h>

typedef bam1_t *bam1_p;

static int change_SO(bam_hdr_t *h, const char *so)
{
    char *p, *q, *beg = NULL, *end = NULL, *newtext;
    if (h->l_text > 3) {
        if (strncmp(h->text, "@HD", 3) == 0) {
            if ((p = strchr(h->text, '\n')) == 0) return -1;
            *p = '\0';
            if ((q = strstr(h->text, "\tSO:")) != 0) {
                *p = '\n'; // change back
                if (strncmp(q + 4, so, p - q - 4) != 0) {
                    beg = q;
                    for (q += 4; *q != '\n' && *q != '\t'; ++q);
                    end = q;
                } else return 0; // no need to change
            } else beg = end = p, *p = '\n';
        }
    }
    if (beg == NULL) { // no @HD
        h->l_text += strlen(so) + 15;
        newtext = (char*)malloc(h->l_text + 1);
        sprintf(newtext, "@HD\tVN:1.3\tSO:%s\n", so);
        strcat(newtext, h->text);
    } else { // has @HD but different or no SO
        h->l_text = (beg - h->text) + (4 + strlen(so)) + (h->text + h->l_text - end);
        newtext = (char*)malloc(h->l_text + 1);
        strncpy(newtext, h->text, beg - h->text);
        sprintf(newtext + (beg - h->text), "\tSO:%s", so);
        strcat(newtext, end);
    }
    free(h->text);
    h->text = newtext;
    return 0;
}

// Function to compare reads and determine which one is < the other
static inline int bam1_lt(const bam1_p a, const bam1_p b)
{
    if (g_is_by_qname) {
        int t = strnum_cmp(bam_get_qname(a), bam_get_qname(b));
        return (t < 0 || (t == 0 && (a->core.flag&0xc0) < (b->core.flag&0xc0)));
    } else return (((uint64_t)a->core.tid<<32|(a->core.pos+1)<<1|bam_is_rev(a)) < ((uint64_t)b->core.tid<<32|(b->core.pos+1)<<1|bam_is_rev(b)));
}
KSORT_INIT(sort, bam1_p, bam1_lt)

typedef struct {
    size_t buf_len;
    const char *prefix;
    bam1_p *buf;
    const bam_hdr_t *h;
    int index;
} worker_t;

static void write_buffer(const char *fn, const char *mode, size_t l, bam1_p *buf, const bam_hdr_t *h, int n_threads)
{
    size_t i;
    samFile* fp;
    fp = sam_open(fn, mode);
    if (fp == NULL) return;
    sam_hdr_write(fp, h);
    if (n_threads > 1) hts_set_threads(fp, n_threads);
    for (i = 0; i < l; ++i)
        sam_write1(fp, h, buf[i]);
    sam_close(fp);
}

static void *worker(void *data)
{
    worker_t *w = (worker_t*)data;
    char *name;
    ks_mergesort(sort, w->buf_len, w->buf, 0);
    name = (char*)calloc(strlen(w->prefix) + 20, 1);
    sprintf(name, "%s.%.4d.bam", w->prefix, w->index);
    write_buffer(name, "wb1", w->buf_len, w->buf, w->h, 0);
    free(name);
    return 0;
}

static int sort_blocks(int n_files, size_t k, bam1_p *buf, const char *prefix, const bam_hdr_t *h, int n_threads)
{
    int i;
    size_t rest;
    bam1_p *b;
    pthread_t *tid;
    pthread_attr_t attr;
    worker_t *w;

    if (n_threads < 1) n_threads = 1;
    if (k < n_threads * 64) n_threads = 1; // use a single thread if we only sort a small batch of records
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    w = (worker_t*)calloc(n_threads, sizeof(worker_t));
    tid = (pthread_t*)calloc(n_threads, sizeof(pthread_t));
    b = buf; rest = k;
    for (i = 0; i < n_threads; ++i) {
        w[i].buf_len = rest / (n_threads - i);
        w[i].buf = b;
        w[i].prefix = prefix;
        w[i].h = h;
        w[i].index = n_files + i;
        b += w[i].buf_len; rest -= w[i].buf_len;
        pthread_create(&tid[i], &attr, worker, &w[i]);
    }
    for (i = 0; i < n_threads; ++i) pthread_join(tid[i], 0);
    free(tid); free(w);
    return n_files + n_threads;
}

/*!
  @abstract Sort an unsorted BAM file based on the chromosome order
  and the leftmost position of an alignment

  @param  is_by_qname whether to sort by query name
  @param  fn       name of the file to be sorted
  @param  prefix   prefix of the temporary files (prefix.NNNN.bam are written)
  @param  fnout    name of the final output file to be written
  @param  modeout  sam_open() mode to be used to create the final output file
  @param  max_mem  approxiate maximum memory (very inaccurate)
  @return 0 for successful sorting, negative on errors

  @discussion It may create multiple temporary subalignment files
  and then merge them by calling bam_merge_core(). This function is
  NOT thread safe.
 */
int bam_sort_core_ext(int is_by_qname, const char *fn, const char *prefix, const char *fnout, const char *modeout, size_t _max_mem, int n_threads)
{
    int ret, i, n_files = 0;
    size_t mem, max_k, k, max_mem;
    bam_hdr_t *header;
    samFile *fp;
    bam1_t *b, **buf;

    if (n_threads < 2) n_threads = 1;
    g_is_by_qname = is_by_qname;
    max_k = k = 0; mem = 0;
    max_mem = _max_mem * n_threads;
    buf = NULL;
    fp = sam_open(fn, "r");
    if (fp == NULL) {
        fprintf(stderr, "[bam_sort_core] fail to open file %s\n", fn);
        return -1;
    }
    header = sam_hdr_read(fp);
    if (is_by_qname) change_SO(header, "queryname");
    else change_SO(header, "coordinate");
    // write sub files
    for (;;) {
        if (k == max_k) {
            size_t kk, old_max = max_k;
            max_k = max_k? max_k<<1 : 0x10000;
            buf = (bam1_t**)realloc(buf, max_k * sizeof(bam1_t*));
            for (kk = old_max; kk < max_k; ++kk) buf[kk] = NULL;
        }
        if (buf[k] == NULL) buf[k] = bam_init1();
        b = buf[k];
        if ((ret = sam_read1(fp, header, b)) < 0) break;
        if (b->l_data < b->m_data>>2) { // shrink
            b->m_data = b->l_data;
            kroundup32(b->m_data);
            b->data = (uint8_t*)realloc(b->data, b->m_data);
        }
        mem += sizeof(bam1_t) + b->m_data + sizeof(void*) + sizeof(void*); // two sizeof(void*) for the data allocated to pointer arrays
        ++k;
        if (mem >= max_mem) {
            n_files = sort_blocks(n_files, k, buf, prefix, header, n_threads);
            mem = k = 0;
        }
    }
    if (ret != -1)
        fprintf(stderr, "[bam_sort_core] truncated file. Continue anyway.\n");
    // write the final output
    if (n_files == 0) { // a single block
        ks_mergesort(sort, k, buf, 0);
        write_buffer(fnout, modeout, k, buf, header, n_threads);
    } else { // then merge
        char **fns;
        n_files = sort_blocks(n_files, k, buf, prefix, header, n_threads);
        fprintf(stderr, "[bam_sort_core] merging from %d files...\n", n_files);
        fns = (char**)calloc(n_files, sizeof(char*));
        for (i = 0; i < n_files; ++i) {
            fns[i] = (char*)calloc(strlen(prefix) + 20, 1);
            sprintf(fns[i], "%s.%.4d.bam", prefix, i);
        }
        if (bam_merge_core2(is_by_qname, fnout, modeout, NULL, n_files, fns, MERGE_COMBINE_RG|MERGE_COMBINE_PG, NULL, n_threads) < 0) {
            // Propagate bam_merge_core2() failure; it has already emitted a
            // message explaining the failure, so no further message is needed.
            return -1;
        }
        for (i = 0; i < n_files; ++i) {
            unlink(fns[i]);
            free(fns[i]);
        }
        free(fns);
    }
    // free
    for (k = 0; k < max_k; ++k) bam_destroy1(buf[k]);
    free(buf);
    bam_hdr_destroy(header);
    sam_close(fp);
    return 0;
}

int bam_sort_core(int is_by_qname, const char *fn, const char *prefix, size_t max_mem)
{
    int ret;
    char *fnout = calloc(strlen(prefix) + 4 + 1, 1);
    sprintf(fnout, "%s.bam", prefix);
    ret = bam_sort_core_ext(is_by_qname, fn, prefix, fnout, "wb", max_mem, 0);
    free(fnout);
    return ret;
}

static int sort_usage(FILE *fp, int status)
{
    fprintf(fp,
"Usage: samtools sort [options...] [in.bam]\n"
"Options:\n"
"  -l INT     Set compression level, from 0 (uncompressed) to 9 (best)\n"
"  -m INT     Set maximum memory per thread; suffix K/M/G recognized [768M]\n"
"  -n         Sort by read name\n"
"  -o FILE    Write final output to FILE rather than standard output\n"
"  -O FORMAT  Write output as FORMAT ('sam'/'bam'/'cram')   (either -O or\n"
"  -T PREFIX  Write temporary files to PREFIX.nnnn.bam       -T is required)\n"
"  -@ INT     Set number of sorting and compression threads [1]\n"
"\n"
"Legacy usage: samtools sort [options...] <in.bam> <out.prefix>\n"
"Options:\n"
"  -f         Use <out.prefix> as full final filename rather than prefix\n"
"  -o         Write final output to stdout rather than <out.prefix>.bam\n"
"  -l,m,n,@   Similar to corresponding options above\n");
    return status;
}

int bam_sort(int argc, char *argv[])
{
    size_t max_mem = 768<<20; // 512MB
    int c, i, modern, nargs, is_by_qname = 0, is_stdout = 0, ret = EXIT_SUCCESS, n_threads = 0, level = -1, full_path = 0;
    char *fnout = "-", *fmtout = NULL, modeout[12], *tmpprefix = NULL;
    kstring_t fnout_buffer = { 0, 0, NULL };

    modern = 0;
    for (i = 1; i < argc; ++i)
        if (argv[i][0] == '-' && strpbrk(argv[i], "OT")) { modern = 1; break; }

    while ((c = getopt(argc, argv, modern? "l:m:no:O:T:@:" : "fnom:@:l:")) >= 0) {
        switch (c) {
        case 'f': full_path = 1; break;
        case 'o': if (modern) fnout = optarg; else is_stdout = 1; break;
        case 'n': is_by_qname = 1; break;
        case 'm': {
                char *q;
                max_mem = strtol(optarg, &q, 0);
                if (*q == 'k' || *q == 'K') max_mem <<= 10;
                else if (*q == 'm' || *q == 'M') max_mem <<= 20;
                else if (*q == 'g' || *q == 'G') max_mem <<= 30;
                break;
            }
        case 'O': fmtout = optarg; break;
        case 'T': tmpprefix = optarg; break;
        case '@': n_threads = atoi(optarg); break;
        case 'l': level = atoi(optarg); break;
        default: return sort_usage(stderr, EXIT_FAILURE);
        }
    }

    nargs = argc - optind;
    if (argc == 1)
        return sort_usage(stdout, EXIT_SUCCESS);
    else if (modern? (nargs > 1) : (nargs != 2))
        return sort_usage(stderr, EXIT_FAILURE);

    if (!modern) {
        fmtout = "bam";
        if (is_stdout) fnout = "-";
        else if (full_path) fnout = argv[optind+1];
        else {
            ksprintf(&fnout_buffer, "%s.%s", argv[optind+1], fmtout);
            fnout = fnout_buffer.s;
        }
        tmpprefix = argv[optind+1];
    }

    strcpy(modeout, "w");
    if (sam_open_mode(&modeout[1], fnout, fmtout) < 0) {
        if (fmtout) fprintf(stderr, "[bam_sort] can't parse output format \"%s\"\n", fmtout);
        else fprintf(stderr, "[bam_sort] can't determine output format\n");
        ret = EXIT_FAILURE;
        goto sort_end;
    }
    if (level >= 0) sprintf(strchr(modeout, '\0'), "%d", level < 9? level : 9);

    if (tmpprefix == NULL) {
        fprintf(stderr, "[bam_sort] no prefix specified for temporary files (use -T option)\n");
        ret = EXIT_FAILURE;
        goto sort_end;
    }

    if (bam_sort_core_ext(is_by_qname, (nargs > 0)? argv[optind] : "-", tmpprefix, fnout, modeout, max_mem, n_threads) < 0) ret = EXIT_FAILURE;

sort_end:
    free(fnout_buffer.s);
    return ret;
}
