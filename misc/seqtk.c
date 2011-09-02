#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <stdint.h>
#include <zlib.h>
#include <string.h>
#include <unistd.h>
#include <limits.h>

#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

typedef struct {
	int n, m;
	uint64_t *a;
} reglist_t;

#include "khash.h"
KHASH_MAP_INIT_STR(reg, reglist_t)

typedef kh_reg_t reghash_t;

reghash_t *stk_reg_read(const char *fn)
{
	reghash_t *h = kh_init(reg);
	gzFile fp;
	kstream_t *ks;
	int dret;
	kstring_t *str;
	// read the list
	str = calloc(1, sizeof(kstring_t));
	fp = strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
	ks = ks_init(fp);
	while (ks_getuntil(ks, 0, str, &dret) >= 0) {
		int beg = -1, end = -1;
		reglist_t *p;
		khint_t k = kh_get(reg, h, str->s);
		if (k == kh_end(h)) {
			int ret;
			char *s = strdup(str->s);
			k = kh_put(reg, h, s, &ret);
			memset(&kh_val(h, k), 0, sizeof(reglist_t));
		}
		p = &kh_val(h, k);
		if (dret != '\n') {
			if (ks_getuntil(ks, 0, str, &dret) > 0 && isdigit(str->s[0])) {
				beg = atoi(str->s);
				if (dret != '\n') {
					if (ks_getuntil(ks, 0, str, &dret) > 0 && isdigit(str->s[0])) {
						end = atoi(str->s);
						if (end < 0) end = -1;
					}
				}
			}
		}
		// skip the rest of the line
		if (dret != '\n') while ((dret = ks_getc(ks)) > 0 && dret != '\n');
		if (end < 0 && beg > 0) end = beg, beg = beg - 1; // if there is only one column
		if (beg < 0) beg = 0, end = INT_MAX;
		if (p->n == p->m) {
			p->m = p->m? p->m<<1 : 4;
			p->a = realloc(p->a, p->m * 8);
		}
		p->a[p->n++] = (uint64_t)beg<<32 | end;
	}
	ks_destroy(ks);
	gzclose(fp);
	free(str->s); free(str);
	return h;
}

void stk_reg_destroy(reghash_t *h)
{
	khint_t k;
	for (k = 0; k < kh_end(h); ++k) {
		if (kh_exist(h, k)) {
			free(kh_val(h, k).a);
			free((char*)kh_key(h, k));
		}
	}
	kh_destroy(reg, h);
}

/* constant table */

unsigned char seq_nt16_table[256] = {
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15 /*'-'*/,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15, 1,14, 2, 13,15,15, 4, 11,15,15,12, 15, 3,15,15,
	15,15, 5, 6,  8,15, 7, 9,  0,10,15,15, 15,15,15,15,
	15, 1,14, 2, 13,15,15, 4, 11,15,15,12, 15, 3,15,15,
	15,15, 5, 6,  8,15, 7, 9,  0,10,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15
};

char *seq_nt16_rev_table = "XACMGRSVTWYHKDBN";
unsigned char seq_nt16to4_table[] = { 4, 0, 1, 4, 2, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4 };
int bitcnt_table[] = { 4, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4 };

/* composition */
int stk_comp(int argc, char *argv[])
{
	gzFile fp;
	kseq_t *seq;
	int l, c, upper_only = 0;
	reghash_t *h = 0;
	reglist_t dummy;
	while ((c = getopt(argc, argv, "ur:")) >= 0) {
		switch (c) {
			case 'u': upper_only = 1; break;
			case 'r': h = stk_reg_read(optarg); break;
		}
	}
	if (argc == optind) {
		fprintf(stderr, "Usage:  seqtk comp [-u] [-r in.bed] <in.fa>\n\n");
		fprintf(stderr, "Output format: chr, length, #A, #C, #G, #T, #2, #3, #4, #CpG, #tv, #ts, #CpG-ts\n");
		return 1;
	}
	fp = (strcmp(argv[optind], "-") == 0)? gzdopen(fileno(stdin), "r") : gzopen(argv[optind], "r");
	seq = kseq_init(fp);
	dummy.n= dummy.m = 1; dummy.a = calloc(1, 8);
	while ((l = kseq_read(seq)) >= 0) {
		int i, k;
		reglist_t *p = 0;
		if (h) {
			khint_t k = kh_get(reg, h, seq->name.s);
			if (k != kh_end(h)) p = &kh_val(h, k);
		} else {
			p = &dummy;
			dummy.a[0] = l;
		}
		for (k = 0; p && k < p->n; ++k) {
			int beg = p->a[k]>>32, end = p->a[k]&0xffffffff;
			int la, lb, lc, na, nb, nc, cnt[11];
			if (beg > 0) la = seq->seq.s[beg-1], lb = seq_nt16_table[la], lc = bitcnt_table[lb];
			else la = 'a', lb = -1, lc = 0;
			na = seq->seq.s[beg]; nb = seq_nt16_table[na]; nc = bitcnt_table[nb];
			memset(cnt, 0, 11 * sizeof(int));
			for (i = beg; i < end; ++i) {
				int is_CpG = 0, a, b, c;
				a = na; b = nb; c = nc;
				na = seq->seq.s[i+1]; nb = seq_nt16_table[na]; nc = bitcnt_table[nb];
				if (b == 2 || b == 10) { // C or Y
					if (nb == 4 || nb == 5) is_CpG = 1;
				} else if (b == 4 || b == 5) { // G or R
					if (lb == 2 || lb == 10) is_CpG = 1;
				}
				if (upper_only == 0 || isupper(a)) {
					if (c > 1) ++cnt[c+2];
					if (c == 1) ++cnt[seq_nt16to4_table[b]];
					if (b == 10 || b == 5) ++cnt[9];
					else if (c == 2) {
						++cnt[8];
					}
					if (is_CpG) {
						++cnt[7];
						if (b == 10 || b == 5) ++cnt[10];
					}
				}
				la = a; lb = b; lc = c;
			}
			if (h) printf("%s\t%d\t%d", seq->name.s, beg, end);
			else printf("%s\t%d", seq->name.s, l);
			for (i = 0; i < 11; ++i) printf("\t%d", cnt[i]);
			putchar('\n');
		}
		fflush(stdout);
	}
	free(dummy.a);
	kseq_destroy(seq);
	gzclose(fp);
	return 0;
}

int stk_randbase(int argc, char *argv[])
{
	gzFile fp;
	kseq_t *seq;
	int l;
	if (argc == 1) {
		fprintf(stderr, "Usage: seqtk randbase <in.fa>\n");
		return 1;
	}
	fp = (strcmp(argv[1], "-") == 0)? gzdopen(fileno(stdin), "r") : gzopen(argv[1], "r");
	seq = kseq_init(fp);
	while ((l = kseq_read(seq)) >= 0) {
		int i;
		printf(">%s", seq->name.s);
		for (i = 0; i < l; ++i) {
			int c, b, a, j, k, m;
			b = seq->seq.s[i];
			c = seq_nt16_table[b];
			a = bitcnt_table[c];
			if (a == 2) {
				m = (drand48() < 0.5);
				for (j = k = 0; j < 4; ++j) {
					if ((1<<j & c) == 0) continue;
					if (k == m) break;
					++k;
				}
				seq->seq.s[i] = islower(b)? "acgt"[j] : "ACGT"[j];
			}
			if (i%60 == 0) putchar('\n');
			putchar(seq->seq.s[i]);
		}
		putchar('\n');
	}
	kseq_destroy(seq);
	gzclose(fp);
	return 0;
}

int stk_hety(int argc, char *argv[])
{
	gzFile fp;
	kseq_t *seq;
	int l, c, win_size = 50000, n_start = 5, win_step, is_lower_mask = 0;
	char *buf;
	uint32_t cnt[3];
	if (argc == 1) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   seqtk hety [options] <in.fa>\n\n");
		fprintf(stderr, "Options: -w INT   window size [%d]\n", win_size);
		fprintf(stderr, "         -t INT   # start positions in a window [%d]\n", n_start);
		fprintf(stderr, "         -m       treat lowercases as masked\n");
		fprintf(stderr, "\n");
		return 1;
	}
	while ((c = getopt(argc, argv, "w:t:m")) >= 0) {
		switch (c) {
		case 'w': win_size = atoi(optarg); break;
		case 't': n_start = atoi(optarg); break;
		case 'm': is_lower_mask = 1; break;
		}
	}
	fp = (strcmp(argv[optind], "-") == 0)? gzdopen(fileno(stdin), "r") : gzopen(argv[optind], "r");
	seq = kseq_init(fp);
	win_step = win_size / n_start;
	buf = calloc(win_size, 1);
	while ((l = kseq_read(seq)) >= 0) {
		int x, i, y, z, next = 0;
		cnt[0] = cnt[1] = cnt[2] = 0;
		for (i = 0; i <= l; ++i) {
			if ((i >= win_size && i % win_step == 0) || i == l) {
				if (i == l && l >= win_size) {
					for (y = l - win_size; y < next; ++y) --cnt[(int)buf[y % win_size]];
				}
				if (cnt[1] + cnt[2] > 0)
					printf("%s\t%d\t%d\t%.2lf\t%d\t%d\n", seq->name.s, next, i,
						   (double)cnt[2] / (cnt[1] + cnt[2]) * win_size, cnt[1] + cnt[2], cnt[2]);
				next = i;
			}
			if (i < l) {
				y = i % win_size;
				c = seq->seq.s[i];
				if (is_lower_mask && islower(c)) c = 'N';
				c = seq_nt16_table[c];
				x = bitcnt_table[c];
				if (i >= win_size) --cnt[(int)buf[y]];
				buf[y] = z = x > 2? 0 : x == 2? 2 : 1;
				++cnt[z];
			}
		}
	}
	free(buf);
	kseq_destroy(seq);
	gzclose(fp);
	return 0;
}

/* fq2fa */
int stk_fq2fa(int argc, char *argv[])
{
	gzFile fp;
	kseq_t *seq;
	char *buf;
	int l, i, c, qual_thres = 0, linelen = 60;
	while ((c = getopt(argc, argv, "q:l:")) >= 0) {
		switch (c) {
			case 'q': qual_thres = atoi(optarg); break;
			case 'l': linelen = atoi(optarg); break;
		}
	}
	if (argc == optind) {
		fprintf(stderr, "Usage: seqtk fq2fa [-q qualThres=0] [-l lineLen=60] <in.fq>\n");
		return 1;
	}
	buf = linelen > 0? malloc(linelen + 1) : 0;
	fp = strcmp(argv[optind], "-")? gzopen(argv[optind], "r") : gzdopen(fileno(stdin), "r");
	seq = kseq_init(fp);
	while ((l = kseq_read(seq)) >= 0) {
		if (seq->qual.l && qual_thres > 0) {
			for (i = 0; i < l; ++i)
				if (seq->qual.s[i] - 33 < qual_thres)
					seq->seq.s[i] = tolower(seq->seq.s[i]);
		}
		putchar('>');
		if (seq->comment.l) {
			fputs(seq->name.s, stdout);
			putchar(' ');
			puts(seq->comment.s);
		} else puts(seq->name.s);
		if (buf) { // multi-line
			for (i = 0; i < l; i += linelen) {
				int x = i + linelen < l? linelen : l - i;
				memcpy(buf, seq->seq.s + i, x);
				buf[x] = 0;
				puts(buf);
			}
		} else puts(seq->seq.s);
	}
	free(buf);
	kseq_destroy(seq);
	gzclose(fp);
	return 0;
}

int stk_maskseq(int argc, char *argv[])
{
	khash_t(reg) *h = kh_init(reg);
	gzFile fp;
	kseq_t *seq;
	int l, i, j, c, is_complement = 0, is_lower = 0;
	khint_t k;
	while ((c = getopt(argc, argv, "cl")) >= 0) {
		switch (c) {
		case 'c': is_complement = 1; break;
		case 'l': is_lower = 1; break;
		}
	}
	if (argc - optind < 2) {
		fprintf(stderr, "Usage:   seqtk maskseq [-cl] <in.fa> <in.bed>\n\n");
		fprintf(stderr, "Options: -c     mask the complement regions\n");
		fprintf(stderr, "         -l     soft mask (to lower cases)\n");
		return 1;
	}
	h = stk_reg_read(argv[optind+1]);
	// maskseq
	fp = strcmp(argv[optind], "-")? gzopen(argv[optind], "r") : gzdopen(fileno(stdin), "r");
	seq = kseq_init(fp);
	while ((l = kseq_read(seq)) >= 0) {
		k = kh_get(reg, h, seq->name.s);
		if (k == kh_end(h)) { // not found in the hash table
			if (is_complement) {
				for (j = 0; j < l; ++j)
					seq->seq.s[j] = is_lower? tolower(seq->seq.s[j]) : 'N';
			}
		} else {
			reglist_t *p = &kh_val(h, k);
			if (!is_complement) {
				for (i = 0; i < p->n; ++i) {
					int beg = p->a[i]>>32, end = p->a[i];
					if (beg >= seq->seq.l) {
						fprintf(stderr, "[maskseq] start position >= the sequence length.\n");
						continue;
					}
					if (end >= seq->seq.l) end = seq->seq.l;
					if (is_lower) for (j = beg; j < end; ++j) seq->seq.s[j] = tolower(seq->seq.s[j]);
					else for (j = beg; j < end; ++j) seq->seq.s[j] = 'N';
				}
			} else {
				int8_t *mask = calloc(seq->seq.l, 1);
				for (i = 0; i < p->n; ++i) {
					int beg = p->a[i]>>32, end = p->a[i];
					if (end >= seq->seq.l) end = seq->seq.l;
					for (j = beg; j < end; ++j) mask[j] = 1;
				}
				for (j = 0; j < l; ++j)
					if (mask[j] == 0) seq->seq.s[j] = is_lower? tolower(seq->seq.s[j]) : 'N';
				free(mask);
			}
		}
		printf(">%s", seq->name.s);
		for (j = 0; j < seq->seq.l; ++j) {
			if (j%60 == 0) putchar('\n');
			putchar(seq->seq.s[j]);
		}
		putchar('\n');
	}
	// free
	kseq_destroy(seq);
	gzclose(fp);
	stk_reg_destroy(h);
	return 0;
}

/* subseq */

int stk_subseq(int argc, char *argv[])
{
	khash_t(reg) *h = kh_init(reg);
	gzFile fp;
	kseq_t *seq;
	int l, i, j, c, is_tab = 0;
	khint_t k;
	while ((c = getopt(argc, argv, "t")) >= 0) {
		switch (c) {
		case 't': is_tab = 1; break;
		}
	}
	if (optind + 2 > argc) {
		fprintf(stderr, "Usage: seqtk subseq [-t] <in.fa> <in.bed>\n\n");
		fprintf(stderr, "Note: Use 'samtools faidx' if only a few regions are intended.\n");
		return 1;
	}
	h = stk_reg_read(argv[optind+1]);
	// subseq
	fp = strcmp(argv[optind], "-")? gzopen(argv[optind], "r") : gzdopen(fileno(stdin), "r");
	seq = kseq_init(fp);
	while ((l = kseq_read(seq)) >= 0) {
		reglist_t *p;
		k = kh_get(reg, h, seq->name.s);
		if (k == kh_end(h)) continue;
		p = &kh_val(h, k);
		for (i = 0; i < p->n; ++i) {
			int beg = p->a[i]>>32, end = p->a[i];
			if (beg >= seq->seq.l) {
				fprintf(stderr, "[subseq] %s: %d >= %ld\n", seq->name.s, beg, seq->seq.l);
				continue;
			}
			if (end > seq->seq.l) end = seq->seq.l;
			if (is_tab == 0) {
				printf("%c%s", seq->qual.l == seq->seq.l? '@' : '>', seq->name.s);
				if (end == INT_MAX) {
					if (beg) printf(":%d", beg+1);
				} else printf(":%d-%d", beg+1, end);
			} else printf("%s\t%d\t", seq->name.s, beg + 1);
			if (end > seq->seq.l) end = seq->seq.l;
			for (j = 0; j < end - beg; ++j) {
				if (is_tab == 0 && j % 60 == 0) putchar('\n');
				putchar(seq->seq.s[j + beg]);
			}
			putchar('\n');
			if (seq->qual.l != seq->seq.l || is_tab) continue;
			printf("+");
			for (j = 0; j < end - beg; ++j) {
				if (j % 60 == 0) putchar('\n');
				putchar(seq->qual.s[j + beg]);
			}
			putchar('\n');
		}
	}
	// free
	kseq_destroy(seq);
	gzclose(fp);
	stk_reg_destroy(h);
	return 0;
}

/* mergefa */
int stk_mergefa(int argc, char *argv[])
{
	gzFile fp[2];
	kseq_t *seq[2];
	int i, l, c, is_intersect = 0, is_haploid = 0, qual = 0, is_mask = 0;
	while ((c = getopt(argc, argv, "himq:")) >= 0) {
		switch (c) {
			case 'i': is_intersect = 1; break;
			case 'h': is_haploid = 1; break;
			case 'm': is_mask = 1; break;
			case 'q': qual = atoi(optarg); break;
		}
	}
	if (is_mask && is_intersect) {
		fprintf(stderr, "[%s] `-i' and `-h' cannot be applied at the same time.\n", __func__);
		return 1;
	}
	if (optind + 2 > argc) {
		fprintf(stderr, "\nUsage: seqtk mergefa [options] <in1.fa> <in2.fa>\n\n");
		fprintf(stderr, "Options: -q INT   quality threshold [0]\n");
		fprintf(stderr, "         -i       take intersection\n");
		fprintf(stderr, "         -m       convert to lowercase when one of the input base is N.\n");
		fprintf(stderr, "         -h       suppress hets in the input\n\n");
		return 1;
	}
	for (i = 0; i < 2; ++i) {
		fp[i] = strcmp(argv[optind+i], "-")? gzopen(argv[optind+i], "r") : gzdopen(fileno(stdin), "r");
		seq[i] = kseq_init(fp[i]);
	}
	while (kseq_read(seq[0]) >= 0) {
		int min_l, c[2], is_upper;
		kseq_read(seq[1]);
		if (strcmp(seq[0]->name.s, seq[1]->name.s))
			fprintf(stderr, "[%s] Different sequence names: %s != %s\n", __func__, seq[0]->name.s, seq[1]->name.s);
		if (seq[0]->seq.l != seq[1]->seq.l)
			fprintf(stderr, "[%s] Unequal sequence length: %ld != %ld\n", __func__, seq[0]->seq.l, seq[1]->seq.l);
		min_l = seq[0]->seq.l < seq[1]->seq.l? seq[0]->seq.l : seq[1]->seq.l;
		printf(">%s", seq[0]->name.s);
		for (l = 0; l < min_l; ++l) {
			c[0] = seq[0]->seq.s[l]; c[1] = seq[1]->seq.s[l];
			if (seq[0]->qual.l && seq[0]->qual.s[l] - 33 < qual) c[0] = tolower(c[0]);
			if (seq[1]->qual.l && seq[1]->qual.s[l] - 33 < qual) c[1] = tolower(c[1]);
			if (is_intersect) is_upper = (isupper(c[0]) || isupper(c[1]))? 1 : 0;
			else if (is_mask) is_upper = (isupper(c[0]) || isupper(c[1]))? 1 : 0;
			else is_upper = (isupper(c[0]) && isupper(c[1]))? 1 : 0;
			c[0] = seq_nt16_table[c[0]]; c[1] = seq_nt16_table[c[1]];
			if (c[0] == 0) c[0] = 15;
			if (c[1] == 0) c[1] = 15;
			if (is_haploid && (bitcnt_table[c[0]] > 1 || bitcnt_table[c[1]] > 1)) is_upper = 0;
			if (is_intersect) {
				c[0] = c[0] & c[1];
				if (c[0] == 0) is_upper = 0;
			} else if (is_mask) {
				if (c[0] == 15 || c[1] == 15) is_upper = 0;
				c[0] = c[0] & c[1];
				if (c[0] == 0) is_upper = 0;
			} else c[0] = c[0] | c[1];
			c[0] = seq_nt16_rev_table[c[0]];
			if (!is_upper) c[0] = tolower(c[0]);
			if (l%60 == 0) putchar('\n');
			putchar(c[0]);
		}
		putchar('\n');
	}
	return 0;
}

int stk_famask(int argc, char *argv[])
{
	gzFile fp[2];
	kseq_t *seq[2];
	int i, l;
	if (argc < 3) {
		fprintf(stderr, "Usage: seqtk famask <src.fa> <mask.fa>\n");
		return 1;
	}
	for (i = 0; i < 2; ++i) {
		fp[i] = strcmp(argv[optind+i], "-")? gzopen(argv[optind+i], "r") : gzdopen(fileno(stdin), "r");
		seq[i] = kseq_init(fp[i]);
	}
	while (kseq_read(seq[0]) >= 0) {
		int min_l, c[2];
		kseq_read(seq[1]);
		if (strcmp(seq[0]->name.s, seq[1]->name.s))
			fprintf(stderr, "[%s] Different sequence names: %s != %s\n", __func__, seq[0]->name.s, seq[1]->name.s);
		if (seq[0]->seq.l != seq[1]->seq.l)
			fprintf(stderr, "[%s] Unequal sequence length: %ld != %ld\n", __func__, seq[0]->seq.l, seq[1]->seq.l);
		min_l = seq[0]->seq.l < seq[1]->seq.l? seq[0]->seq.l : seq[1]->seq.l;
		printf(">%s", seq[0]->name.s);
		for (l = 0; l < min_l; ++l) {
			c[0] = seq[0]->seq.s[l]; c[1] = seq[1]->seq.s[l];
			if (c[1] == 'x') c[0] = tolower(c[0]);
			else if (c[1] != 'X') c[0] = c[1];
			if (l%60 == 0) putchar('\n');
			putchar(c[0]);
		}
		putchar('\n');
	}
	return 0;
}

int stk_mutfa(int argc, char *argv[])
{
	khash_t(reg) *h = kh_init(reg);
	gzFile fp;
	kseq_t *seq;
	kstream_t *ks;
	int l, i, dret;
	kstring_t *str;
	khint_t k;
	if (argc < 3) {
		fprintf(stderr, "Usage: seqtk mutfa <in.fa> <in.snp>\n\n");
		fprintf(stderr, "Note: <in.snp> contains at least four columns per line which are:\n");
		fprintf(stderr, "      'chr  1-based-pos  any  base-changed-to'.\n");
		return 1;
	}
	// read the list
	str = calloc(1, sizeof(kstring_t));
	fp = strcmp(argv[2], "-")? gzopen(argv[2], "r") : gzdopen(fileno(stdin), "r");
	ks = ks_init(fp);
	while (ks_getuntil(ks, 0, str, &dret) >= 0) {
		char *s = strdup(str->s);
		int beg = 0, ret;
		reglist_t *p;
		k = kh_get(reg, h, s);
		if (k == kh_end(h)) {
			k = kh_put(reg, h, s, &ret);
			memset(&kh_val(h, k), 0, sizeof(reglist_t));
		}
		p = &kh_val(h, k);
		if (ks_getuntil(ks, 0, str, &dret) > 0) beg = atol(str->s) - 1; // 2nd col
		ks_getuntil(ks, 0, str, &dret); // 3rd col
		ks_getuntil(ks, 0, str, &dret); // 4th col
		// skip the rest of the line
		if (dret != '\n') while ((dret = ks_getc(ks)) > 0 && dret != '\n');
		if (isalpha(str->s[0]) && str->l == 1) {
			if (p->n == p->m) {
				p->m = p->m? p->m<<1 : 4;
				p->a = realloc(p->a, p->m * 8);
			}
			p->a[p->n++] = (uint64_t)beg<<32 | str->s[0];
		}
	}
	ks_destroy(ks);
	gzclose(fp);
	free(str->s); free(str);
	// mutfa
	fp = strcmp(argv[1], "-")? gzopen(argv[1], "r") : gzdopen(fileno(stdin), "r");
	seq = kseq_init(fp);
	while ((l = kseq_read(seq)) >= 0) {
		reglist_t *p;
		k = kh_get(reg, h, seq->name.s);
		if (k != kh_end(h)) {
			p = &kh_val(h, k);
			for (i = 0; i < p->n; ++i) {
				int beg = p->a[i]>>32;
				if (beg < seq->seq.l)
					seq->seq.s[beg] = (int)p->a[i];
			}
		}
		printf(">%s", seq->name.s);
		for (i = 0; i < l; ++i) {
			if (i%60 == 0) putchar('\n');
			putchar(seq->seq.s[i]);
		}
		putchar('\n');
	}
	// free
	kseq_destroy(seq);
	gzclose(fp);
	for (k = 0; k < kh_end(h); ++k) {
		if (kh_exist(h, k)) {
			free(kh_val(h, k).a);
			free((char*)kh_key(h, k));
		}
	}
	kh_destroy(reg, h);
	return 0;
}

int stk_listhet(int argc, char *argv[])
{
	gzFile fp;
	kseq_t *seq;
	int i, l;
	if (argc == 1) {
		fprintf(stderr, "Usage: seqtk listhet <in.fa>\n");
		return 1;
	}
	fp = (strcmp(argv[1], "-") == 0)? gzdopen(fileno(stdin), "r") : gzopen(argv[1], "r");
	seq = kseq_init(fp);
	while ((l = kseq_read(seq)) >= 0) {
		for (i = 0; i < l; ++i) {
			int b = seq->seq.s[i];
			if (bitcnt_table[seq_nt16_table[b]] == 2)
				printf("%s\t%d\t%c\n", seq->name.s, i+1, b);
		}
	}
	kseq_destroy(seq);
	gzclose(fp);
	return 0;
}

/* cutN */
static int cutN_min_N_tract = 1000;
static int cutN_nonN_penalty = 10;

static int find_next_cut(const kseq_t *ks, int k, int *begin, int *end)
{
	int i, b, e;
	while (k < ks->seq.l) {
		if (seq_nt16_table[(int)ks->seq.s[k]] == 15) {
			int score, max;
			score = 0; e = max = -1;
			for (i = k; i < ks->seq.l && score >= 0; ++i) { /* forward */
				if (seq_nt16_table[(int)ks->seq.s[i]] == 15) ++score;
				else score -= cutN_nonN_penalty;
				if (score > max) max = score, e = i;
			}
			score = 0; b = max = -1;
			for (i = e; i >= 0 && score >= 0; --i) { /* backward */
				if (seq_nt16_table[(int)ks->seq.s[i]] == 15) ++score;
				else score -= cutN_nonN_penalty;
				if (score > max) max = score, b = i;
			}
			if (e + 1 - b >= cutN_min_N_tract) {
				*begin = b;
				*end = e + 1;
				return *end;
			}
			k = e + 1;
		} else ++k;
	}
	return -1;
}
static void print_seq(FILE *fpout, const kseq_t *ks, int begin, int end)
{
	int i;
	if (begin >= end) return; // FIXME: why may this happen? Understand it!
	fprintf(fpout, ">%s:%d-%d", ks->name.s, begin+1, end);
	for (i = begin; i < end && i < ks->seq.l; ++i) {
		if ((i - begin)%60 == 0) fputc('\n', fpout);
		fputc(ks->seq.s[i], fpout);
	}
	fputc('\n', fpout);
}
int stk_cutN(int argc, char *argv[])
{
	int c, l, gap_only = 0;
	gzFile fp;
	kseq_t *ks;
	while ((c = getopt(argc, argv, "n:p:g")) >= 0) {
		switch (c) {
		case 'n': cutN_min_N_tract = atoi(optarg); break;
		case 'p': cutN_nonN_penalty = atoi(optarg); break;
		case 'g': gap_only = 1; break;
		default: return 1;
		}
	}
	if (argc == optind) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   seqtk cutN [options] <in.fa>\n\n");
		fprintf(stderr, "Options: -n INT    min size of N tract [%d]\n", cutN_min_N_tract);
		fprintf(stderr, "         -p INT    penalty for a non-N [%d]\n", cutN_nonN_penalty);
		fprintf(stderr, "         -g        print gaps only, no sequence\n\n");
		return 1;
	}
	fp = (strcmp(argv[optind], "-") == 0)? gzdopen(fileno(stdin), "r") : gzopen(argv[optind], "r");
	ks = kseq_init(fp);
	while ((l = kseq_read(ks)) >= 0) {
		int k = 0, begin = 0, end = 0;
		while (find_next_cut(ks, k, &begin, &end) >= 0) {
			if (begin != 0) {
				if (gap_only) printf("%s\t%d\t%d\n", ks->name.s, begin, end);
				else print_seq(stdout, ks, k, begin);
			}
			k = end;
		}
		if (!gap_only) print_seq(stdout, ks, k, l);
	}
	kseq_destroy(ks);
	gzclose(fp);
	return 0;
}

/* main function */
static int usage()
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage:   seqtk <command> <arguments>\n\n");
	fprintf(stderr, "Command: comp      get the nucleotide composite of FASTA/Q\n");
	fprintf(stderr, "         hety      regional heterozygosity\n");
	fprintf(stderr, "         fq2fa     convert FASTQ to FASTA\n");
	fprintf(stderr, "         subseq    extract subsequences from FASTA/Q\n");
	fprintf(stderr, "         maskseq   mask sequences\n");
	fprintf(stderr, "         mutfa     point mutate FASTA at specified positions\n");
	fprintf(stderr, "         mergefa   merge two FASTA/Q files\n");
	fprintf(stderr, "         randbase  choose a random base from hets\n");
	fprintf(stderr, "         cutN      cut sequence at long N\n");
	fprintf(stderr, "         listhet   extract the position of each het\n");
	fprintf(stderr, "\n");
	return 1;
}

int main(int argc, char *argv[])
{
	if (argc == 1) return usage();
	if (strcmp(argv[1], "comp") == 0) stk_comp(argc-1, argv+1);
	else if (strcmp(argv[1], "hety") == 0) stk_hety(argc-1, argv+1);
	else if (strcmp(argv[1], "fq2fa") == 0) stk_fq2fa(argc-1, argv+1);
	else if (strcmp(argv[1], "subseq") == 0) stk_subseq(argc-1, argv+1);
	else if (strcmp(argv[1], "maskseq") == 0) stk_maskseq(argc-1, argv+1);
	else if (strcmp(argv[1], "mutfa") == 0) stk_mutfa(argc-1, argv+1);
	else if (strcmp(argv[1], "mergefa") == 0) stk_mergefa(argc-1, argv+1);
	else if (strcmp(argv[1], "randbase") == 0) stk_randbase(argc-1, argv+1);
	else if (strcmp(argv[1], "cutN") == 0) stk_cutN(argc-1, argv+1);
	else if (strcmp(argv[1], "listhet") == 0) stk_listhet(argc-1, argv+1);
	else if (strcmp(argv[1], "famask") == 0) stk_famask(argc-1, argv+1);
	else {
		fprintf(stderr, "[main] unrecognized commad '%s'. Abort!\n", argv[1]);
		return 1;
	}
	return 0;
}
