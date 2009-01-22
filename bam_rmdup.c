#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <zlib.h>
#include "bam.h"

typedef bam1_t *bam1_p;
#include "khash.h"
KHASH_SET_INIT_STR(name)
KHASH_MAP_INIT_INT64(pos, bam1_p)

#define BUFFER_SIZE 0x40000

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

static inline void dump_best(tmp_stack_t *stack, khash_t(pos) *best_hash, bamFile out)
{
	int i;
	for (i = 0; i != stack->n; ++i) {
		bam_write1(out, stack->a[i]);
		bam_destroy1(stack->a[i]);
	}
	stack->n = 0;
	if (kh_size(best_hash) > BUFFER_SIZE) kh_clear(pos, best_hash);
}

static void clear_del_set(khash_t(name) *del_set)
{
	khint_t k;
	for (k = kh_begin(del_set); k < kh_end(del_set); ++k)
		if (kh_exist(del_set, k))
			free((char*)kh_key(del_set, k));
	kh_clear(name, del_set);
}

void bam_rmdup_core(bamFile in, bamFile out)
{
	bam_header_t *header;
	bam1_t *b;
	int last_tid = -1, last_pos = -1;
	uint64_t n_checked = 0, n_removed = 0;
	tmp_stack_t stack;
	khint_t k;
	khash_t(pos) *best_hash;
	khash_t(name) *del_set;

	best_hash = kh_init(pos);
	del_set = kh_init(name);
	b = bam_init1();
	memset(&stack, 0, sizeof(tmp_stack_t));
	header = bam_header_read(in);
	bam_header_write(out, header);

	kh_resize(name, del_set, 4 * BUFFER_SIZE);
	kh_resize(pos, best_hash, 3 * BUFFER_SIZE);
	while (bam_read1(in, b) >= 0) {
		bam1_core_t *c = &b->core;
		if (c->tid != last_tid || last_pos != c->pos) {
			dump_best(&stack, best_hash, out); // write the result
			if (c->tid != last_tid) {
				kh_clear(pos, best_hash);
				if (kh_size(del_set)) { // check
					fprintf(stderr, "[bam_rmdup_core] %llu unmatched pairs\n", (long long)kh_size(del_set));
					clear_del_set(del_set);
				}
				last_tid = c->tid;
				fprintf(stderr, "[bam_rmdup_core] processing reference %s...\n", header->target_name[c->tid]);
			}
		}
		if (!(c->flag&BAM_FPAIRED) || (c->flag&(BAM_FUNMAP|BAM_FMUNMAP)) || (c->mtid >= 0 && c->tid != c->mtid)) {
			bam_write1(out, b);
		} else if (c->isize > 0) { // paired, head
			uint64_t key = (uint64_t)c->pos<<32 | c->isize;
			int ret;
			++n_checked;
			k = kh_put(pos, best_hash, key, &ret);
			if (ret == 0) { // found in best_hash
				bam1_t *p = kh_val(best_hash, k);
				++n_removed;
				if (p->core.qual < c->qual) { // the current alignment is better
					kh_put(name, del_set, strdup(bam1_qname(p)), &ret); // p will be removed
					bam_copy1(p, b); // replaced as b
				} else kh_put(name, del_set, strdup(bam1_qname(b)), &ret); // b will be removed
				if (ret == 0)
					fprintf(stderr, "[bam_rmdup_core] inconsistent BAM file for pair '%s'. Continue anyway.\n", bam1_qname(b));
			} else { // not found in best_hash
				kh_val(best_hash, k) = bam_dup1(b);
				stack_insert(&stack, kh_val(best_hash, k));
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
	dump_best(&stack, best_hash, out);

	bam_header_destroy(header);
	clear_del_set(del_set);
	kh_destroy(name, del_set);
	kh_destroy(pos, best_hash);
	free(stack.a);
	bam_destroy1(b);
	fprintf(stderr, "[bam_rmdup_core] %lld / %lld = %.4lf\n", (long long)n_removed, (long long)n_checked,
			(double)n_removed/n_checked);
}
int bam_rmdup(int argc, char *argv[])
{
	bamFile in, out;
	if (argc < 3) {
		fprintf(stderr, "Usage: samtools rmdup <input.srt.bam> <output.bam>\n");
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
