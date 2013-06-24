#include <stdlib.h>
#include <ctype.h>
#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <string.h>
#include <regex.h>
#include <unistd.h>
#include "bam.h"
#include "htslib/ksort.h"
#include "htslib/khash.h"
#include "htslib/klist.h"

KHASH_INIT(c2c, char*, char*, 1, kh_str_hash_func, kh_str_hash_equal)

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

#define HEAP_EMPTY 0xffffffffffffffffull

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
		if (a.b == 0 || b.b == 0) return a.b == 0? 1 : 0;
		t = strnum_cmp(bam1_qname(a.b), bam1_qname(b.b));
		return (t > 0 || (t == 0 && (a.b->core.flag&0xc0) > (b.b->core.flag&0xc0)));
	} else return __pos_cmp(a, b);
}

KSORT_INIT(heap, heap1_t, heap_lt)

static void swap_header_targets(bam_header_t *h1, bam_header_t *h2)
{
	bam_header_t t;
	t.n_targets = h1->n_targets, h1->n_targets = h2->n_targets, h2->n_targets = t.n_targets;
	t.target_name = h1->target_name, h1->target_name = h2->target_name, h2->target_name = t.target_name;
	t.target_len = h1->target_len, h1->target_len = h2->target_len, h2->target_len = t.target_len;
}

static void swap_header_text(bam_header_t *h1, bam_header_t *h2)
{
	int tempi;
	char *temps;
	tempi = h1->l_text, h1->l_text = h2->l_text, h2->l_text = tempi;
	temps = h1->text, h1->text = h2->text, h2->text = temps;
}

typedef struct trans_tbl {
	int* tid_trans;
	kh_c2c_t* rg_trans;
	kh_c2c_t* pg_trans;
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

static void trans_tbl_init(bam_header_t* out, bam_header_t* translate, trans_tbl_t* tbl)
{
	tbl->tid_trans = (int*)calloc(translate->n_targets, sizeof(int));
	tbl->rg_trans = kh_init(c2c);
	tbl->pg_trans = kh_init(c2c);
	if (!tbl->tid_trans || !tbl->rg_trans || !tbl->pg_trans) { perror("out of memory"); exit(-1); }
	
	// Naive way of doing this but meh
	int i, j;
	for (i = 0; i < translate->n_targets; ++i) {
		int tid = -1;
		// Search 'out' for entries in 'translate' and map them
		for (j = 0; j < out->n_targets; j++) {
			if (!strcmp(translate->target_name[i],out->target_name[j])) {
				tid = j;
				break;
			}
		}
		if (tid == -1) { // Append missing entries to out
			tbl->tid_trans[i] = out->n_targets++;
			out->target_name = (char**)realloc(out->target_name, sizeof(char*)*out->n_targets);
			out->target_name[out->n_targets-1] = strdup(translate->target_name[i]);
			out->target_len = (uint32_t*)realloc(out->target_len, sizeof(uint32_t)*out->n_targets);
			out->target_len[out->n_targets-1] = translate->target_len[i];
			// grep line with regex '^@SQ.*\tID:%s(\t.*$|$)', translate->target_name[i]
			// from translate->text
			regex_t sq_id;
			regmatch_t* matches = calloc(2, sizeof(regmatch_t));
			if (matches == NULL) { perror("out of memory"); exit(-1); }
			char* seq_regex = NULL;
			asprintf(&seq_regex, "^@SQ.*\tID:%s(\t.*$|$)",translate->target_name[i]);
			regcomp(&sq_id, seq_regex, REG_EXTENDED|REG_NEWLINE);
			free(seq_regex);
			if (regexec(&sq_id, translate->text, 1, matches, 0) != 0)
			{
				fprintf(stderr, "[trans_tbl_init] @SQ ID (%s) found in binary header but not text header.\n",translate->target_name[i]);
				exit(-1);
			}
			regfree(&sq_id);
			
			// Produce our output line
			char* match_line = strndup(translate->text+matches[0].rm_so, matches[0].rm_eo-matches[0].rm_so);
			// and append it to out->text
			char* newtext = NULL;
			asprintf(&newtext, "%s\n%s",out->text,match_line);
			free(out->text);
			out->text = newtext;
			free(match_line);
			free(matches);
		} else {
			tbl->tid_trans[i] = tid;
		}
	}
	
	// grep @RG id's
	regex_t rg_id;
	regmatch_t* matches = calloc(2, sizeof(regmatch_t));
	if (matches == NULL) { perror("out of memory"); exit(-1); }
	regcomp(&rg_id, "^@RG.*\tID:([A-Za-z0-9]*)(\t.*$|$)", REG_EXTENDED|REG_NEWLINE);
	char* text = translate->text;
	while(1) { //	foreach rg id in translate's header
		if (regexec(&rg_id, text, 2, matches, 0) != 0) break;
		char* match_id = strndup(text+matches[1].rm_so, matches[1].rm_eo-matches[1].rm_so);
		
		// is our matched ID in our output list already
		regex_t rg_id_search;
		char* rg_regex = NULL;
		asprintf(&rg_regex, "^@RG.*\tID:%s(\t.*$|$)",match_id);
		regcomp(&rg_id_search, rg_regex, REG_EXTENDED|REG_NEWLINE|REG_NOSUB);
		free(rg_regex);
		char* transformed_id = NULL;
		if (regexec(&rg_id_search, out->text, 0, NULL, 0) != 0) {
			// Not in there so can add it as 1-1 mapping
			transformed_id = match_id;
		} else {
			// It's in there so we need to transform it by appending random number to id
			if (asprintf(&transformed_id, "%s-%ld",match_id, lrand48()) == -1) { perror("out of memory"); exit(-1); }
		}
		regfree(&rg_id_search);

		// Insert it into our translation map
		int in_there = 0;
		khiter_t iter = kh_put(c2c, tbl->rg_trans, strdup(match_id), &in_there);
		kh_value(tbl->rg_trans,iter) = strdup(transformed_id);
		// take matched line and replace ID with transformed_id
		char* transformed_line = NULL;
		if (match_id != transformed_id) {
			char *fmt = NULL;
			asprintf(&fmt, "%%.%ds%%s%%.%ds",matches[1].rm_so-matches[0].rm_so,matches[0].rm_eo-matches[1].rm_eo);
			asprintf(&transformed_line,fmt,text+matches[0].rm_so,transformed_id,text+matches[1].rm_eo);
			free(fmt);
			free(transformed_id);
		} else {
			transformed_line = strndup(text+matches[0].rm_so, matches[0].rm_eo-matches[0].rm_so);
		}
		
		// and append it to out->text
		char* newtext = NULL;
		asprintf(&newtext, "%s\n%s",out->text,transformed_line);
		free(out->text);
		out->text = newtext;
		
		free(transformed_line);
		free(match_id); match_id = NULL;
		text += matches[0].rm_eo;
	}
	regfree(&rg_id);
	
	// Do same for PG id's
	regex_t pg_id;
	regcomp(&pg_id, "^@PG.*\tID:([A-Za-z0-9]*)(\t.*$|$)", REG_EXTENDED|REG_NEWLINE);
	text = translate->text;
	klist_t(hdrln) *pg_list = kl_init(hdrln);
	while(1) { //	foreach pg id in translate's header
		if (regexec(&pg_id, text, 2, matches, 0) != 0) break;
		char* match_id = strndup(text+matches[1].rm_so, matches[1].rm_eo-matches[1].rm_so);

		// is our matched ID in our output list already
		regex_t pg_id_search;
		char* pg_regex = NULL;
		asprintf(&pg_regex, "^@PG.*\tID:%s(\t.*$|$)",match_id);
		regcomp(&pg_id_search, pg_regex, REG_EXTENDED|REG_NEWLINE|REG_NOSUB);
		free(pg_regex);
		char* transformed_id = NULL;
		if (regexec(&pg_id_search, out->text, 0, NULL, 0) != 0) {
			// Not in there so can add it as 1-1 mapping
			transformed_id = match_id;
		} else {
			// It's in there so we need to transform it by appending random number to id
			if (asprintf(&transformed_id, "%s-%ld",match_id, lrand48()) == -1) { perror("out of memory"); exit(-1); }
		}
		regfree(&pg_id_search);
		
		// Insert it into our translation map
		int in_there = 0;
		khiter_t iter = kh_put(c2c, tbl->pg_trans, strdup(match_id), &in_there);
		kh_value(tbl->pg_trans,iter) = strdup(transformed_id);
		// take matched line and replace ID with transformed_id
		char* transformed_line = NULL;
		if (match_id != transformed_id) {
			char *fmt = NULL;
			asprintf(&fmt, "%%.%ds%%s%%.%ds",matches[1].rm_so-matches[0].rm_so,matches[0].rm_eo-matches[1].rm_eo);
			asprintf(&transformed_line,fmt,text+matches[0].rm_so,transformed_id,text+matches[1].rm_eo);
			free(fmt);
			free(transformed_id);
		} else {
			transformed_line = strndup(text+matches[0].rm_so, matches[0].rm_eo-matches[0].rm_so);
		}

		// append line to linked list for PP processing
		char** ln = kl_pushp(hdrln, pg_list);
		*ln = strdup(transformed_line);  // Give linked list it's own copy
		free(transformed_line);
		free(match_id);
		text += matches[0].rm_eo;
	}
	regfree(&pg_id);
	// need to translate PP's on the fly (in second pass because they may not be in correct order and need complete tbl->pg_trans to do this
	// for each line {
	// with ID replaced with tranformed_id and PP's transformed using the translation table
	// }
	regex_t pg_pp;
	regcomp(&pg_pp, "^@PG.*\tPP:([A-Za-z0-9]*)(\t.*$|$)", REG_EXTENDED|REG_NEWLINE);
	kliter_t(hdrln) *iter = kl_begin(pg_list);
	while (iter != kl_end(pg_list)) {
		char* data = kl_val(iter);

		char *transformed_line = NULL;
		// Find PP tag
		if (regexec(&pg_pp, data, 2, matches, 0) == 0) {
			// Lookup in hash table
			char* pp_id = strndup(data+matches[1].rm_so, matches[1].rm_eo-matches[1].rm_so);

			khiter_t k = kh_get(c2c, tbl->pg_trans, pp_id);
			free(pp_id);
			char* transformed_id = kh_value(tbl->pg_trans,k);
			// Replace
			char *fmt = NULL;
			asprintf(&fmt, "%%.%ds%%s%%.%ds",matches[1].rm_so-matches[0].rm_so,matches[0].rm_eo-matches[1].rm_eo);
			asprintf(&transformed_line,fmt,data,transformed_id,data+matches[1].rm_eo);
			free(fmt);
		} else { transformed_line = data; }
		// and append it to out->text
		char* newtext = NULL;
		asprintf(&newtext, "%s\n%s",out->text,transformed_line);
		free(out->text);
		out->text = newtext;
		if (transformed_line != data) { free(transformed_line); transformed_line = NULL; }
		free(data);
		iter = kl_next(iter);
	}
	regfree(&pg_pp);
	
	kl_destroy(hdrln,pg_list);
	
	free(matches);
	out->l_text = strlen(out->text);
}

static void bam_translate(bam1_t* b, trans_tbl_t* tbl)
{
	b->core.tid = tbl->tid_trans[b->core.tid];
	b->core.mtid = tbl->tid_trans[b->core.mtid];
	uint8_t *rg = bam_aux_get(b, "RG");
	if (rg) {
		char* decoded_rg = bam_aux2Z(rg);
		khiter_t k = kh_get(c2c, tbl->rg_trans, decoded_rg);
		if (k == kh_end(tbl->rg_trans)) abort();
		char* translate_rg = kh_value(tbl->rg_trans,k);
		bam_aux_del(b, rg);
		bam_aux_append(b, "RG", 'Z', strlen(translate_rg) + 1, (uint8_t*)translate_rg);
	}

	uint8_t *pg = bam_aux_get(b, "PG");
	if (pg) {
		char* decoded_pg = bam_aux2Z(pg);
		khiter_t k = kh_get(c2c, tbl->pg_trans, decoded_pg);
		if (k == kh_end(tbl->pg_trans)) abort();
		char* translate_pg = kh_value(tbl->pg_trans,k);
		bam_aux_del(b, pg);
		bam_aux_append(b, "PG", 'Z', strlen(translate_pg) + 1, (uint8_t*)translate_pg);
	}
}


#define MERGE_RG     1 // Attach RG tag based on filename
#define MERGE_UNCOMP 2 // Generate uncompressed BAM
#define MERGE_LEVEL1 4 // Compress the BAM at level 1 (fast) mode
#define MERGE_FORCE  8

/*
 * How merging is handled
 *
 * if a hheader is defined use that for output header
 * otherwise take the first header from the first file
 * 
 * Now go through each file and create a translation table for that file for:
 * -RG
 * -tid
 * -PG tags
 *
 * Then when reading a record from a bam translate the read before stashing it in the hash.
 */

/*!
  @abstract    Merge multiple sorted BAM.
  @param  is_by_qname whether to sort by query name
  @param  out  output BAM file name
  @param  headers  name of SAM file from which to copy '@' header lines,
                   or NULL to copy them from the first file to be merged
  @param  n    number of files to be merged
  @param  fn   names of files to be merged

  @discussion Padding information may NOT correctly maintained. This
  function is NOT thread safe.
 */
int bam_merge_core2(int by_qname, const char *out, const char *headers, int n, char * const *fn, int flag, const char *reg, int n_threads, int level)
{
	bamFile fpout, *fp;
	heap1_t *heap;
	bam_header_t *hout = NULL;
	int i, j, *RG_len = NULL;
	uint64_t idx = 0;
	char **RG = NULL, mode[8];
	bam_iter_t *iter = NULL;
	trans_tbl_t *translation_tbl = NULL;

    // Is there a specified pre-prepared header to use for output?
	if (headers) {
		tamFile fpheaders = sam_open(headers);
		if (fpheaders == 0) {
			const char *message = strerror(errno);
			fprintf(stderr, "[bam_merge_core] cannot open '%s': %s\n", headers, message);
			return -1;
		}
		hout = sam_header_read(fpheaders);
		sam_close(fpheaders);
	}

	g_is_by_qname = by_qname;
	fp = (bamFile*)calloc(n, sizeof(bamFile));
	heap = (heap1_t*)calloc(n, sizeof(heap1_t));
	iter = (bam_iter_t*)calloc(n, sizeof(bam_iter_t));
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
			RG[i] = calloc(l + 1, 1);
			RG_len[i] = l;
			strncpy(RG[i], s + j, l);
		}
	}
	// open and read the header from each file
	for (i = 0; i < n; ++i) {
		bam_header_t *hin;
		fp[i] = bam_open(fn[i], "r");
		if (fp[i] == 0) {
			int j;
			fprintf(stderr, "[bam_merge_core] fail to open file %s\n", fn[i]);
			for (j = 0; j < i; ++j) bam_close(fp[j]);
			free(fp); free(heap);
			// FIXME: possible memory leak
			return -1;
		}
		hin = bam_header_read(fp[i]);
		if (hout == NULL) hout = hin;
		trans_tbl_init(hout, hin, translation_tbl+i);
	}

    // If we're only merging a specified region move our iters to start at that point
	if (reg) {
		int tid, beg, end;
		if (bam_parse_region(hout, reg, &tid, &beg, &end) < 0) {
			fprintf(stderr, "[%s] Malformated region string or undefined reference name\n", __func__);
			return -1;
		}
		for (i = 0; i < n; ++i) {
			bam_index_t *idx;
			idx = bam_index_load(fn[i]);
			iter[i] = bam_iter_query(idx, tid, beg, end);
			bam_index_destroy(idx);
		}
	}

    // Load the first read from each file into the heap
	for (i = 0; i < n; ++i) {
		heap1_t *h = heap + i;
		h->i = i;
		h->b = (bam1_t*)calloc(1, sizeof(bam1_t));
		if (bam_iter_read(fp[i], iter[i], h->b) >= 0) {
			h->pos = ((uint64_t)h->b->core.tid<<32) | (uint32_t)((int32_t)h->b->core.pos+1)<<1 | bam1_strand(h->b);
			h->idx = idx++;
			bam_translate(h->b, translation_tbl + i);
		}
		else {
			h->pos = HEAP_EMPTY;
			free(heap->b);
			heap->b = NULL;
		}
	}
    // Open output file and write header
	if (flag & MERGE_UNCOMP) level = 0;
	else if (flag & MERGE_LEVEL1) level = 1;
	strcpy(mode, "w");
	if (level >= 0) sprintf(mode + 1, "%d", level < 9? level : 9);
	if ((fpout = strcmp(out, "-")? bam_open(out, "w") : bam_dopen(fileno(stdout), "w")) == 0) {
		fprintf(stderr, "[%s] fail to create the output file.\n", __func__);
		return -1;
	}
	bam_header_write(fpout, hout);
	bam_header_destroy(hout);
	if (!(flag & MERGE_UNCOMP)) bgzf_mt(fpout, n_threads, 256);

    // Begin the actual merge
	ks_heapmake(heap, n, heap);
	while (heap->pos != HEAP_EMPTY) {
		bam1_t *b = heap->b;
		if (flag & MERGE_RG) {
			uint8_t *rg = bam_aux_get(b, "RG");
			if (rg) bam_aux_del(b, rg);
			bam_aux_append(b, "RG", 'Z', RG_len[heap->i] + 1, (uint8_t*)RG[heap->i]);
		}
		bam_write1(fpout, b);
		if ((j = bam_iter_read(fp[heap->i], iter[heap->i], b)) >= 0) {
			heap->pos = ((uint64_t)b->core.tid<<32) | (uint32_t)((int)b->core.pos+1)<<1 | bam1_strand(b);
			heap->idx = idx++;
			bam_translate(b, translation_tbl + heap->i);
		} else if (j == -1) {
			heap->pos = HEAP_EMPTY;
			free(heap->b->data); free(heap->b);
			heap->b = NULL;
		} else fprintf(stderr, "[bam_merge_core] '%s' is truncated. Continue anyway.\n", fn[heap->i]);
		ks_heapadjust(heap, 0, n, heap);
	}

    // Clean up and close
	if (flag & MERGE_RG) {
		for (i = 0; i != n; ++i) free(RG[i]);
		free(RG); free(RG_len);
	}
	for (i = 0; i != n; ++i) {
		trans_tbl_destroy(translation_tbl + i);
		bam_iter_destroy(iter[i]);
		bam_close(fp[i]);
	}
	bam_close(fpout);
	free(fp); free(heap); free(iter);
	return 0;
}

int bam_merge_core(int by_qname, const char *out, const char *headers, int n, char * const *fn, int flag, const char *reg)
{
	return bam_merge_core2(by_qname, out, headers, n, fn, flag, reg, 0, -1);
}

int bam_merge(int argc, char *argv[])
{
	int c, is_by_qname = 0, flag = 0, ret = 0, n_threads = 0, level = -1;
	char *fn_headers = NULL, *reg = 0;

	while ((c = getopt(argc, argv, "h:nru1R:f@:l:")) >= 0) {
		switch (c) {
		case 'r': flag |= MERGE_RG; break;
		case 'f': flag |= MERGE_FORCE; break;
		case 'h': fn_headers = strdup(optarg); break;
		case 'n': is_by_qname = 1; break;
		case '1': flag |= MERGE_LEVEL1; break;
		case 'u': flag |= MERGE_UNCOMP; break;
		case 'R': reg = strdup(optarg); break;
		case 'l': level = atoi(optarg); break;
		case '@': n_threads = atoi(optarg); break;
		}
	}
	if (optind + 2 >= argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   samtools merge [-nr] [-h inh.sam] <out.bam> <in1.bam> <in2.bam> [...]\n\n");
		fprintf(stderr, "Options: -n       sort by read names\n");
		fprintf(stderr, "         -r       attach RG tag (inferred from file names)\n");
		fprintf(stderr, "         -u       uncompressed BAM output\n");
		fprintf(stderr, "         -f       overwrite the output BAM if exist\n");
		fprintf(stderr, "         -1       compress level 1\n");
		fprintf(stderr, "         -l INT   compression level, from 0 to 9 [-1]\n");
		fprintf(stderr, "         -@ INT   number of BAM compression threads [0]\n");
		fprintf(stderr, "         -R STR   merge file in the specified region STR [all]\n");
		fprintf(stderr, "         -h FILE  copy the header in FILE to <out.bam> [in1.bam]\n\n");
		fprintf(stderr, "Note: Samtools' merge does not reconstruct the @RG dictionary in the header. Users\n");
		fprintf(stderr, "      must provide the correct header with -h, or uses Picard which properly maintains\n");
		fprintf(stderr, "      the header dictionary in merging.\n\n");
		return 1;
	}
	if (!(flag & MERGE_FORCE) && strcmp(argv[optind], "-")) {
		FILE *fp = fopen(argv[optind], "rb");
		if (fp != NULL) {
			fclose(fp);
			fprintf(stderr, "[%s] File '%s' exists. Please apply '-f' to overwrite. Abort.\n", __func__, argv[optind]);
			return 1;
		}
	}
	if (bam_merge_core2(is_by_qname, argv[optind], fn_headers, argc - optind - 1, argv + optind + 1, flag, reg, n_threads, level) < 0) ret = 1;
	free(reg);
	free(fn_headers);
	return ret;
}

/***************
 * BAM sorting *
 ***************/

#include <pthread.h>

typedef bam1_t *bam1_p;

static int change_SO(bam_header_t *h, const char *so)
{
	char *p, *q, *beg = 0, *end = 0, *newtext;
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
	if (beg == 0) { // no @HD
		h->l_text += strlen(so) + 15;
		newtext = malloc(h->l_text + 1);
		sprintf(newtext, "@HD\tVN:1.3\tSO:%s\n", so);
		strcat(newtext, h->text);
	} else { // has @HD but different or no SO
		h->l_text = (beg - h->text) + (4 + strlen(so)) + (h->text + h->l_text - end);
		newtext = malloc(h->l_text + 1);
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
		int t = strnum_cmp(bam1_qname(a), bam1_qname(b));
		return (t < 0 || (t == 0 && (a->core.flag&0xc0) < (b->core.flag&0xc0)));
	} else return (((uint64_t)a->core.tid<<32|(a->core.pos+1)<<1|bam1_strand(a)) < ((uint64_t)b->core.tid<<32|(b->core.pos+1)<<1|bam1_strand(b)));
}
KSORT_INIT(sort, bam1_p, bam1_lt)

typedef struct {
	size_t buf_len;
	const char *prefix;
	bam1_p *buf;
	const bam_header_t *h;
	int index;
} worker_t;

static void write_buffer(const char *fn, const char *mode, size_t l, bam1_p *buf, const bam_header_t *h, int n_threads)
{
	size_t i;
	bamFile fp;
	fp = strcmp(fn, "-")? bam_open(fn, mode) : bam_dopen(fileno(stdout), mode);
	if (fp == 0) return;
	bam_header_write(fp, h);
	if (n_threads > 1) bgzf_mt(fp, n_threads, 256);
	for (i = 0; i < l; ++i)
		bam_write1(fp, buf[i]);
	bam_close(fp);
}

static void *worker(void *data)
{
	worker_t *w = (worker_t*)data;
	char *name;
	ks_mergesort(sort, w->buf_len, w->buf, 0);
	name = (char*)calloc(strlen(w->prefix) + 20, 1);
	sprintf(name, "%s.%.4d.bam", w->prefix, w->index);
	write_buffer(name, "w1", w->buf_len, w->buf, w->h, 0);
	free(name);
	return 0;
}

static int sort_blocks(int n_files, size_t k, bam1_p *buf, const char *prefix, const bam_header_t *h, int n_threads)
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
	w = calloc(n_threads, sizeof(worker_t));
	tid = calloc(n_threads, sizeof(pthread_t));
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
  @param  max_mem  approxiate maximum memory (very inaccurate)
  @return 0 for successful sorting, negative on errors

  @discussion It may create multiple temporary subalignment files
  and then merge them by calling bam_merge_core(). This function is
  NOT thread safe.
 */
int bam_sort_core_ext(int is_by_qname, const char *fn, const char *prefix, const char *fnout, size_t _max_mem, int n_threads, int level)
{
	int ret, i, n_files = 0;
	size_t mem, max_k, k, max_mem;
	bam_header_t *header;
	bamFile fp;
	bam1_t *b, **buf;

	if (n_threads < 2) n_threads = 1;
	g_is_by_qname = is_by_qname;
	max_k = k = 0; mem = 0;
	max_mem = _max_mem * n_threads;
	buf = 0;
	fp = strcmp(fn, "-")? bam_open(fn, "r") : bam_dopen(fileno(stdin), "r");
	if (fp == 0) {
		fprintf(stderr, "[bam_sort_core] fail to open file %s\n", fn);
		return -1;
	}
	header = bam_header_read(fp);
	if (is_by_qname) change_SO(header, "queryname");
	else change_SO(header, "coordinate");
	// write sub files
	for (;;) {
		if (k == max_k) {
			size_t old_max = max_k;
			max_k = max_k? max_k<<1 : 0x10000;
			buf = realloc(buf, max_k * sizeof(bam1_t*));
			memset(buf + old_max, 0, sizeof(void*) * (max_k - old_max));
		}
		if (buf[k] == 0) buf[k] = (bam1_t*)calloc(1, sizeof(bam1_t));
		b = buf[k];
		if ((ret = bam_read1(fp, b)) < 0) break;
		if (b->data_len < b->m_data>>2) { // shrink
			b->m_data = b->data_len;
			kroundup32(b->m_data);
			b->data = realloc(b->data, b->m_data);
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
		char mode[8];
		strcpy(mode, "w");
		if (level >= 0) sprintf(mode + 1, "%d", level < 9? level : 9);
		ks_mergesort(sort, k, buf, 0);
		write_buffer(fnout, mode, k, buf, header, n_threads);
	} else { // then merge
		char **fns;
		n_files = sort_blocks(n_files, k, buf, prefix, header, n_threads);
		fprintf(stderr, "[bam_sort_core] merging from %d files...\n", n_files);
		fns = (char**)calloc(n_files, sizeof(char*));
		for (i = 0; i < n_files; ++i) {
			fns[i] = (char*)calloc(strlen(prefix) + 20, 1);
			sprintf(fns[i], "%s.%.4d.bam", prefix, i);
		}
		if (bam_merge_core2(is_by_qname, fnout, 0, n_files, fns, 0, 0, n_threads, level) < 0) {
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
	for (k = 0; k < max_k; ++k) {
		if (!buf[k]) continue;
		free(buf[k]->data);
		free(buf[k]);
	}
	free(buf);
	bam_header_destroy(header);
	bam_close(fp);
	return 0;
}

int bam_sort_core(int is_by_qname, const char *fn, const char *prefix, size_t max_mem)
{
	int ret;
	char *fnout = calloc(strlen(prefix) + 4 + 1, 1);
	sprintf(fnout, "%s.bam", prefix);
	ret = bam_sort_core_ext(is_by_qname, fn, prefix, fnout, max_mem, 0, -1);
	free(fnout);
	return ret;
}

int bam_sort(int argc, char *argv[])
{
	size_t max_mem = 768<<20; // 512MB
	int c, is_by_qname = 0, is_stdout = 0, ret = 0, n_threads = 0, level = -1, full_path = 0;
	char *fnout;
	while ((c = getopt(argc, argv, "fnom:@:l:")) >= 0) {
		switch (c) {
		case 'f': full_path = 1; break;
		case 'o': is_stdout = 1; break;
		case 'n': is_by_qname = 1; break;
		case 'm': {
				char *q;
				max_mem = strtol(optarg, &q, 0);
				if (*q == 'k' || *q == 'K') max_mem <<= 10;
				else if (*q == 'm' || *q == 'M') max_mem <<= 20;
				else if (*q == 'g' || *q == 'G') max_mem <<= 30;
				break;
			}
		case '@': n_threads = atoi(optarg); break;
		case 'l': level = atoi(optarg); break;
		}
	}
	if (optind + 2 > argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   samtools sort [options] <in.bam> <out.prefix>\n\n");
		fprintf(stderr, "Options: -n        sort by read name\n");
		fprintf(stderr, "         -f        use <out.prefix> as full file name instead of prefix\n");
		fprintf(stderr, "         -o        final output to stdout\n");
		fprintf(stderr, "         -l INT    compression level, from 0 to 9 [-1]\n");
		fprintf(stderr, "         -@ INT    number of sorting and compression threads [1]\n");
		fprintf(stderr, "         -m INT    max memory per thread; suffix K/M/G recognized [768M]\n");
		fprintf(stderr, "\n");
		return 1;
	}

	if (is_stdout) fnout = strdup("-");
	else {
		fnout = calloc(strlen(argv[optind+1]) + 4 + 1, 1);
		sprintf(fnout, "%s%s", argv[optind+1], full_path? "" : ".bam");
	}

	if (bam_sort_core_ext(is_by_qname, argv[optind], argv[optind+1], fnout, max_mem, n_threads, level) < 0) ret = 1;
	free(fnout);
	return ret;
}
