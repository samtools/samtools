#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <zlib.h>
#include "sam.h"

typedef bam1_t *bam1_p;

#include "khash.h"
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

static inline void stack_insert(tmp_stack_t *stack, bam1_t *b)
{
	if (stack->n == stack->max) {
		stack->max = stack->max? stack->max<<1 : 0x10000;
		stack->a = (bam1_t**)realloc(stack->a, sizeof(bam1_t*) * stack->max);
	}
	stack->a[stack->n++] = b;
}

static inline void dump_best(tmp_stack_t *stack, bamFile out)
{
	int i;
	for (i = 0; i != stack->n; ++i) {
		bam_write1(out, stack->a[i]);
		bam_destroy1(stack->a[i]);
	}
	stack->n = 0;
}

static void clear_del_set(khash_t(name) *del_set)
{
	khint_t k;
	for (k = kh_begin(del_set); k < kh_end(del_set); ++k)
		if (kh_exist(del_set, k))
			free((char*)kh_key(del_set, k));
	kh_clear(name, del_set);
}

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

static inline int sum_qual(const bam1_t *b)
{
	int i, q;
	uint8_t *qual = bam1_qual(b);
	for (i = q = 0; i < b->core.l_qseq; ++i) q += qual[i];
	return q;
}

void bam_rmdup_core(bamFile in, bamFile out)
{
	bam_header_t *header;
	bam1_t *b;
	int last_tid = -1, last_pos = -1;
	tmp_stack_t stack;
	khint_t k;
	khash_t(lib) *aux;
	khash_t(name) *del_set;
	
	aux = kh_init(lib);
	del_set = kh_init(name);
	b = bam_init1();
	memset(&stack, 0, sizeof(tmp_stack_t));
	header = bam_header_read(in);
	sam_header_parse_rg(header);
	bam_header_write(out, header);

	kh_resize(name, del_set, 4 * BUFFER_SIZE);
	while (bam_read1(in, b) >= 0) {
		bam1_core_t *c = &b->core;
		if (c->tid != last_tid || last_pos != c->pos) {
			dump_best(&stack, out); // write the result
			clear_best(aux, BUFFER_SIZE);
			if (c->tid != last_tid) {
				clear_best(aux, 0);
				if (kh_size(del_set)) { // check
					fprintf(stderr, "[bam_rmdup_core] %llu unmatched pairs\n", (long long)kh_size(del_set));
					clear_del_set(del_set);
				}
				if ((int)c->tid == -1) { // append unmapped reads
					bam_write1(out, b);
					while (bam_read1(in, b) >= 0) bam_write1(out, b);
					break;
				}
				last_tid = c->tid;
				fprintf(stderr, "[bam_rmdup_core] processing reference %s...\n", header->target_name[c->tid]);
			}
		}
		if (!(c->flag&BAM_FPAIRED) || (c->flag&(BAM_FUNMAP|BAM_FMUNMAP)) || (c->mtid >= 0 && c->tid != c->mtid)) {
			bam_write1(out, b);
		} else if (c->isize > 0) { // paired, head
			uint64_t key = (uint64_t)c->pos<<32 | c->isize;
			const char *lib;
			const uint8_t *rg;
			lib_aux_t *q;
			int ret;
			rg = bam_aux_get(b, "RG");
			lib = (rg == 0)? 0 : bam_strmap_get(header->rg2lib, (char*)(rg + 1));
			q = lib? get_aux(aux, lib) : get_aux(aux, "\t");
			++q->n_checked;
			k = kh_put(pos, q->best_hash, key, &ret);
			if (ret == 0) { // found in best_hash
				bam1_t *p = kh_val(q->best_hash, k);
				++q->n_removed;
				if (sum_qual(p) < sum_qual(b)) { // the current alignment is better; this can be accelerated in principle
					kh_put(name, del_set, strdup(bam1_qname(p)), &ret); // p will be removed
					bam_copy1(p, b); // replaced as b
				} else kh_put(name, del_set, strdup(bam1_qname(b)), &ret); // b will be removed
				if (ret == 0)
					fprintf(stderr, "[bam_rmdup_core] inconsistent BAM file for pair '%s'. Continue anyway.\n", bam1_qname(b));
			} else { // not found in best_hash
				kh_val(q->best_hash, k) = bam_dup1(b);
				stack_insert(&stack, kh_val(q->best_hash, k));
			}
		} else { // paired, tail
			k = kh_get(name, del_set, bam1_qname(b));
			if (k != kh_end(del_set)) {
				free((char*)kh_key(del_set, k));
				kh_del(name, del_set, k);
			} else bam_write1(out, b);
		}
		last_pos = c->pos;
	}

	for (k = kh_begin(aux); k != kh_end(aux); ++k) {
		if (kh_exist(aux, k)) {
			lib_aux_t *q = &kh_val(aux, k);			
			dump_best(&stack, out);
			fprintf(stderr, "[bam_rmdup_core] %lld / %lld = %.4lf in library '%s'\n", (long long)q->n_removed,
					(long long)q->n_checked, (double)q->n_removed/q->n_checked, kh_key(aux, k));
			kh_destroy(pos, q->best_hash);
			free((char*)kh_key(aux, k));
		}
	}
	kh_destroy(lib, aux);

	bam_header_destroy(header);
	clear_del_set(del_set);
	kh_destroy(name, del_set);
	free(stack.a);
	bam_destroy1(b);
}
int bam_rmdup(int argc, char *argv[])
{
	bamFile in, out;
	if (argc < 3) {
		fprintf(stderr, "Usage: samtools rmdup <input.srt.bam> <output.bam>\n\n");
		fprintf(stderr, "Note: Picard is recommended for this task.\n");
		return 1;
	}
	in = (strcmp(argv[1], "-") == 0)? bam_dopen(fileno(stdin), "r") : bam_open(argv[1], "r");
	out = (strcmp(argv[2], "-") == 0)? bam_dopen(fileno(stdout), "w") : bam_open(argv[2], "w");
	if (in == 0 || out == 0) {
		fprintf(stderr, "[bam_rmdup] fail to read/write input files\n");
		return 1;
	}
	bam_rmdup_core(in, out);
	bam_close(in);
	bam_close(out);
	return 0;
}
