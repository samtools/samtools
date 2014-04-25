#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <stdio.h>
#include <zlib.h>
#include <unistd.h>
#include <inttypes.h>
#include <htslib/sam.h>
#include <htslib/khash.h>
#include "hdr_idx_priv.h"

typedef bam1_t *bam1_p;

KHASH_SET_INIT_STR(name)
KHASH_MAP_INIT_INT64(pos, bam1_p)

#define BUFFER_SIZE 0x40000

typedef struct {
	uint64_t n_checked, n_removed;
	khash_t(pos) *best_hash;
} lib_aux_t;
KHASH_MAP_INIT_STR(lib, lib_aux_t)

typedef struct {
	int n, max;
	bam1_t **a;
} tmp_stack_t;

// Adds a read to the stack
static inline void stack_insert(tmp_stack_t *stack, bam1_t *b)
{
	if (stack->n == stack->max) {
		stack->max = stack->max? stack->max<<1 : 0x10000;
		stack->a = (bam1_t**)realloc(stack->a, sizeof(bam1_t*) * stack->max);
	}
	stack->a[stack->n++] = b;
}

// Writes out all reads in the stack to file out and deallocates the memory they're using
static inline void dump_best(tmp_stack_t *stack, samFile* out, const bam_hdr_t *h)
{
	int i;
	for (i = 0; i != stack->n; ++i) {
		sam_write1(out, h, stack->a[i]);
		bam_destroy1(stack->a[i]);
	}
	stack->n = 0;
}

// Deallocate memory used by inside of kh del_set
static void clear_del_set(khash_t(name) *del_set)
{
	khint_t k;
	for (k = kh_begin(del_set); k < kh_end(del_set); ++k)
		if (kh_exist(del_set, k))
			free((char*)kh_key(del_set, k));
	kh_clear(name, del_set);
}

// Given hash aux, and lib find lib in hash and return lib_aux_t that corresponds to it
static lib_aux_t *get_aux(khash_t(lib) *aux, const char *lib)
{
	khint_t k = kh_get(lib, aux, lib);
	if (k == kh_end(aux)) {
		int ret;
		char *p = strdup(lib);
		lib_aux_t *q;
		k = kh_put(lib, aux, p, &ret);
		q = &kh_val(aux, k);
		q->n_checked = q->n_removed = 0;
		q->best_hash = kh_init(pos);
		return q;
	} else return &kh_val(aux, k);
}

static void clear_best(khash_t(lib) *aux, int max)
{
	khint_t k;
	for (k = kh_begin(aux); k != kh_end(aux); ++k) {
		if (kh_exist(aux, k)) {
			lib_aux_t *q = &kh_val(aux, k);
			if (kh_size(q->best_hash) >= max)
				kh_clear(pos, q->best_hash);
		}
	}
}

// Returns total of all quality values in read
static inline int sum_qual(const bam1_t *b)
{
	int i, q;
	uint8_t *qual = bam_get_qual(b);
	for (i = q = 0; i < b->core.l_qseq; ++i) q += qual[i];
	return q;
}


void bam_rmdup_core(samFile* in, samFile* out, bam_hdr_t *h)
{
	bam1_t *b;
	int last_tid = -1, last_pos = -1;
	tmp_stack_t stack;
	khint_t k;
	khash_t(lib) *aux;
	khash_t(name) *del_set;
	library_index_t* li = bam_library_index_init(h);
	
	aux = kh_init(lib);
	del_set = kh_init(name);
	b = bam_init1();
	memset(&stack, 0, sizeof(tmp_stack_t));

	kh_resize(name, del_set, 4 * BUFFER_SIZE);
	while (sam_read1(in, h, b) >= 0) {
		bam1_core_t *c = &b->core;
		// If we are
		if (c->tid != last_tid || last_pos != c->pos) {
			dump_best(&stack, out, h); // write the result
			clear_best(aux, BUFFER_SIZE);
			if (c->tid != last_tid) {
				clear_best(aux, 0);
				if (kh_size(del_set)) { // check
					fprintf(stderr, "[bam_rmdup_core] %"PRIu32" unmatched pairs\n", kh_size(del_set));
					clear_del_set(del_set);
				}
				if ((int)c->tid == -1) { // append unmapped reads
					sam_write1(out, h, b);
					while (sam_read1(in, h, b) >= 0) sam_write1(out, h, b);
					break;
				}
				last_tid = c->tid;
				fprintf(stderr, "[bam_rmdup_core] processing reference %s...\n", h->target_name[c->tid]);
			}
		}
		if (!(c->flag&BAM_FPAIRED) || (c->flag&(BAM_FUNMAP|BAM_FMUNMAP)) || (c->mtid >= 0 && c->tid != c->mtid)) {
			sam_write1(out, h, b);
		} else if (c->isize > 0) { // paired, head
			uint64_t key = (uint64_t)c->pos<<32 | c->isize;
			const char *lib;
			lib_aux_t *q;
			int ret;
			lib = bam_search_library_index(li, b);
			q = lib? get_aux(aux, lib) : get_aux(aux, "\t");
			++q->n_checked;
			k = kh_put(pos, q->best_hash, key, &ret);
			if (ret == 0) { // found in best_hash
				bam1_t *p = kh_val(q->best_hash, k);
				++q->n_removed;
				if (sum_qual(p) < sum_qual(b)) { // the current alignment is better; this can be accelerated in principle
					kh_put(name, del_set, strdup(bam_get_qname(p)), &ret); // p will be removed
					bam_copy1(p, b); // replaced as b
				} else kh_put(name, del_set, strdup(bam_get_qname(b)), &ret); // b will be removed
				if (ret == 0)
					fprintf(stderr, "[bam_rmdup_core] inconsistent BAM file for pair '%s'. Continue anyway.\n", bam_get_qname(b));
			} else { // not found in best_hash
				kh_val(q->best_hash, k) = bam_dup1(b);
				stack_insert(&stack, kh_val(q->best_hash, k));
			}
		} else { // paired, tail
			k = kh_get(name, del_set, bam_get_qname(b));
			if (k != kh_end(del_set)) {
				free((char*)kh_key(del_set, k));
				kh_del(name, del_set, k);
			} else sam_write1(out, h, b);
		}
		last_pos = c->pos;
	}

	for (k = kh_begin(aux); k != kh_end(aux); ++k) {
		if (kh_exist(aux, k)) {
			lib_aux_t *q = &kh_val(aux, k);			
			dump_best(&stack, out, h);
			fprintf(stderr, "[bam_rmdup_core] %" PRId64 " / %" PRId64 " = %.4lf in library '%s'\n", q->n_removed,
					q->n_checked, (double)q->n_removed/q->n_checked, kh_key(aux, k));
			kh_destroy(pos, q->best_hash);
			free((char*)kh_key(aux, k));
		}
	}
	kh_destroy(lib, aux);
	bam_library_index_destroy(li);

	clear_del_set(del_set);
	kh_destroy(name, del_set);
	free(stack.a);
	bam_destroy1(b);
}

void bam_rmdupse_core(samFile* in, samFile* out, bam_hdr_t *h, bool force_se);

static void usage(bool error)
{
	FILE* out = error ? stderr : stdout;

	fprintf(out,
			"Usage:  samtools rmdup [-sS] <input.srt.bam> <output.bam>\n\n"
			"Option: -s    rmdup for SE reads\n"
			"        -S    treat PE reads as SE in rmdup (force -s)\n\n");

}

int bam_rmdup(int argc, char *argv[])
{
	int c;
	bool is_se = false, force_se = false;
	if (argc == 1) {
		usage(false);
		return 0;
	}
	while ((c = getopt(argc, argv, "sSh")) >= 0) {
		switch (c) {
			case 's':
				is_se = true;
				break;
			case 'S':
				force_se = is_se = true;
				break;
			case 'h':
				usage(false);
				return 0;
			case '?':
			default:
				usage(true);
				return 1;
		}
	}
	if (argc-optind != 2) {
		fprintf(stderr, "[bam_rmdup] Invalid number of arguments.\n\n");
		usage(true);
		return 1;
	}
	samFile* in = sam_open(argv[optind], "rb");
	samFile* out = sam_open(argv[optind+1], "wb");
	if (in == NULL || out == NULL) {
		fprintf(stderr, "[bam_rmdup] Failed to read/write input files.\n");
		return 1;
	}
	bam_hdr_t* h = sam_hdr_read(in);
	if (is_se) bam_rmdupse_core(in, out, h, force_se);
	else bam_rmdup_core(in, out, h);
	sam_close(in); sam_close(out);
	bam_hdr_destroy(h);
	return 0;
}
