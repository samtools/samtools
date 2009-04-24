#include <stdio.h>
#include <unistd.h>
#include <ctype.h>
#include "bam.h"
#include "faidx.h"
#include "bam_maqcns.h"
#include "khash.h"
#include "glf.h"
#include "kstring.h"

typedef int *indel_list_t;
KHASH_MAP_INIT_INT64(64, indel_list_t)

#define BAM_PLF_SIMPLE     0x01
#define BAM_PLF_CNS        0x02
#define BAM_PLF_INDEL_ONLY 0x04
#define BAM_PLF_GLF        0x08

typedef struct {
	bam_header_t *h;
	bam_maqcns_t *c;
	bam_maqindel_opt_t *ido;
	faidx_t *fai;
	khash_t(64) *hash;
	uint32_t format;
	int tid, len, last_pos;
	int mask;
	char *ref;
	glfFile fp; // for glf output only
} pu_data_t;

char **__bam_get_lines(const char *fn, int *_n);
void bam_init_header_hash(bam_header_t *header);
int32_t bam_get_tid(const bam_header_t *header, const char *seq_name);

static khash_t(64) *load_pos(const char *fn, bam_header_t *h)
{
	char **list;
	int i, j, n, *fields, max_fields;
	khash_t(64) *hash;
	bam_init_header_hash(h);
	list = __bam_get_lines(fn, &n);
	hash = kh_init(64);
	max_fields = 0; fields = 0;
	for (i = 0; i < n; ++i) {
		char *str = list[i];
		int chr, n_fields, ret;
		khint_t k;
		uint64_t x;
		n_fields = ksplit_core(str, 0, &max_fields, &fields);
		if (n_fields < 2) continue;
		chr = bam_get_tid(h, str + fields[0]);
		if (chr < 0) {
			fprintf(stderr, "[load_pos] unknown reference sequence name: %s\n", str + fields[0]);
			continue;
		}
		x = (uint64_t)chr << 32 | (atoi(str + fields[1]) - 1);
		k = kh_put(64, hash, x, &ret);
		if (ret == 0) {
			fprintf(stderr, "[load_pos] position %s:%s has been loaded.\n", str+fields[0], str+fields[1]);
			continue;
		}
		kh_val(hash, k) = 0;
		if (n_fields > 2) {
			// count
			for (j = 2; j < n_fields; ++j) {
				char *s = str + fields[j];
				if ((*s != '+' && *s != '-') || !isdigit(s[1])) break;
 			}
			if (j > 2) { // update kh_val()
				int *q, y, z;
				q = kh_val(hash, k) = (int*)calloc(j - 1, sizeof(int));
				q[0] = j - 2; z = j; y = 1;
				for (j = 2; j < z; ++j)
					q[y++] = atoi(str + fields[j]);
			}
		}
		free(str);
	}
	free(list); free(fields);
	return hash;
}

// an analogy to pileup_func() below
static int glt3_func(uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pu, void *data)
{
	pu_data_t *d = (pu_data_t*)data;
	bam_maqindel_ret_t *r = 0;
	int rb, *proposed_indels = 0;
	glf1_t *g;
	glf3_t *g3;

	if (d->fai == 0) {
		fprintf(stderr, "[glt3_func] reference sequence is required for generating GLT. Abort!\n");
		exit(1);
	}
	if (d->hash) { // only output a list of sites
		khint_t k = kh_get(64, d->hash, (uint64_t)tid<<32|pos);
		if (k == kh_end(d->hash)) return 0;
		proposed_indels = kh_val(d->hash, k);
	}
	g3 = glf3_init1();
	if (d->fai && (int)tid != d->tid) {
		if (d->ref) { // then write the end mark
			g3->rtype = GLF3_RTYPE_END;
			glf3_write1(d->fp, g3);
		}
		glf3_ref_write(d->fp, d->h->target_name[tid], d->h->target_len[tid]); // write reference
		free(d->ref);
		d->ref = fai_fetch(d->fai, d->h->target_name[tid], &d->len);
		d->tid = tid;
		d->last_pos = 0;
	}
	rb = (d->ref && (int)pos < d->len)? d->ref[pos] : 'N';
	g = bam_maqcns_glfgen(n, pu, bam_nt16_table[rb], d->c);
	memcpy(g3, g, sizeof(glf1_t));
	g3->rtype = GLF3_RTYPE_SUB;
	g3->offset = pos - d->last_pos;
	d->last_pos = pos;
	glf3_write1(d->fp, g3);
	if (proposed_indels)
		r = bam_maqindel(n, pos, d->ido, pu, d->ref, proposed_indels[0], proposed_indels+1);
	else r = bam_maqindel(n, pos, d->ido, pu, d->ref, 0, 0);
	if (r) { // then write indel line
		int het = 3 * n, min;
		min = het;
		if (min > r->gl[0]) min = r->gl[0];
		if (min > r->gl[1]) min = r->gl[1];
		g3->ref_base = 0;
		g3->rtype = GLF3_RTYPE_INDEL;
		memset(g3->lk, 0, 10);
		g3->lk[0] = r->gl[0] - min < 255? r->gl[0] - min : 255;
		g3->lk[1] = r->gl[1] - min < 255? r->gl[1] - min : 255;
		g3->lk[2] = het - min < 255? het - min : 255;
		g3->offset = 0;
		g3->indel_len[0] = r->indel1;
		g3->indel_len[1] = r->indel2;
		g3->min_lk = min < 255? min : 255;
		g3->max_len = (abs(r->indel1) > abs(r->indel2)? abs(r->indel1) : abs(r->indel2)) + 1;
		g3->indel_seq[0] = strdup(r->s[0]+1);
		g3->indel_seq[1] = strdup(r->s[1]+1);
		glf3_write1(d->fp, g3);
		bam_maqindel_ret_destroy(r);
	}
	free(g);
	glf3_destroy1(g3);
	return 0;
}

static int pileup_func(uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pu, void *data)
{
	pu_data_t *d = (pu_data_t*)data;
	bam_maqindel_ret_t *r = 0;
	int i, j, rb, max_mapq = 0, *proposed_indels = 0;
	uint32_t x;

	if (d->format & BAM_PLF_GLF) return glt3_func(tid, pos, n, pu, data);
	if (d->hash) { // only output a list of sites
		khint_t k = kh_get(64, d->hash, (uint64_t)tid<<32|pos);
		if (k == kh_end(d->hash)) return 0;
		proposed_indels = kh_val(d->hash, k);
	}
	if (d->fai && (int)tid != d->tid) { // then update d->ref
		free(d->ref);
		d->ref = fai_fetch(d->fai, d->h->target_name[tid], &d->len);
		d->tid = tid;
	}
	rb = (d->ref && (int)pos < d->len)? d->ref[pos] : 'N';
	if (d->format & BAM_PLF_INDEL_ONLY) {
		for (i = 0; i < n; ++i)
			if (pu[i].indel != 0) break;
		if (i == n) return 0;
	}
	printf("%s\t%d\t%c\t", d->h->target_name[tid], pos + 1, rb);
	if (d->format & BAM_PLF_CNS) { // consensus
		int ref_q, rb4 = bam_nt16_table[rb];
		x = bam_maqcns_call(n, pu, d->c);
		ref_q = 0;
		if (rb4 != 15 && x>>28 != 15 && x>>28 != rb4) { // a SNP
			ref_q = ((x>>24&0xf) == rb4)? x>>8&0xff : (x>>8&0xff) + (x&0xff);
			if (ref_q > 255) ref_q = 255;
		}
		printf("%c\t%d\t%d\t%d\t", bam_nt16_rev_table[x>>28], x>>8&0xff, ref_q, x>>16&0xff);
	}
	if ((d->format & (BAM_PLF_CNS|BAM_PLF_INDEL_ONLY)) && d->ref) {
		if (proposed_indels)
			r = bam_maqindel(n, pos, d->ido, pu, d->ref, proposed_indels[0], proposed_indels+1);
		else r = bam_maqindel(n, pos, d->ido, pu, d->ref, 0, 0);
	}
	// pileup strings
	printf("%d\t", n);
	for (i = 0; i < n; ++i) {
		const bam_pileup1_t *p = pu + i;
		if (max_mapq < p->b->core.qual) max_mapq = p->b->core.qual;
		if (p->is_head) printf("^%c", p->b->core.qual > 93? 126 : p->b->core.qual + 33);
		if (!p->is_del) {
			int c = bam_nt16_rev_table[bam1_seqi(bam1_seq(p->b), p->qpos)];
			if (c == '=' || toupper(c) == toupper(rb)) c = bam1_strand(p->b)? ',' : '.';
			else c = bam1_strand(p->b)? tolower(c) : toupper(c);
			putchar(c);
			if (p->indel > 0) {
				printf("+%d", p->indel);
				for (j = 1; j <= p->indel; ++j) {
					c = bam_nt16_rev_table[bam1_seqi(bam1_seq(p->b), p->qpos + j)];
					putchar(bam1_strand(p->b)? tolower(c) : toupper(c));
				}
			} else if (p->indel < 0) {
				printf("%d", p->indel);
				for (j = 1; j <= -p->indel; ++j) {
					c = (d->ref && (int)pos+j < d->len)? d->ref[pos+j] : 'N';
					putchar(bam1_strand(p->b)? tolower(c) : toupper(c));
				}
			}
		} else putchar('*');
		if (p->is_tail) putchar('$');
	}
	putchar('\t');
	for (i = 0; i < n; ++i) {
		const bam_pileup1_t *p = pu + i;
		int c = bam1_qual(p->b)[p->qpos] + 33;
		if (c > 126) c = 126;
		putchar(c);
	}
	if (d->format & BAM_PLF_SIMPLE) {
		putchar('\t');
		for (i = 0; i < n; ++i) {
			int c = pu[i].b->core.qual + 33;
			if (c > 126) c = 126;
			putchar(c);
		}
	}
	putchar('\n');
	if (r) { // then print indel line
		printf("%s\t%d\t*\t", d->h->target_name[tid], pos + 1);
		if (r->gt < 2) printf("%s/%s\t", r->s[r->gt], r->s[r->gt]);
		else printf("%s/%s\t", r->s[0], r->s[1]);
		printf("%d\t%d\t", r->q_cns, r->q_ref);
		printf("%d\t%d\t", max_mapq, n);
		printf("%s\t%s\t", r->s[0], r->s[1]);
		//printf("%d\t%d\t", r->gl[0], r->gl[1]);
		printf("%d\t%d\t%d\n", r->cnt1, r->cnt2, r->cnt_anti);
		bam_maqindel_ret_destroy(r);
	}
	return 0;
}

int bam_pileup(int argc, char *argv[])
{
	int c;
	char *fn_list = 0, *fn_fa = 0, *fn_pos = 0;
	pu_data_t *d = (pu_data_t*)calloc(1, sizeof(pu_data_t));
	d->tid = -1; d->mask = BAM_DEF_MASK;
	d->c = bam_maqcns_init();
	d->ido = bam_maqindel_opt_init();
	while ((c = getopt(argc, argv, "st:f:cT:N:r:l:im:gI:G:")) >= 0) {
		switch (c) {
		case 's': d->format |= BAM_PLF_SIMPLE; break;
		case 't': fn_list = strdup(optarg); break;
		case 'l': fn_pos = strdup(optarg); break;
		case 'f': fn_fa = strdup(optarg); break;
		case 'T': d->c->theta = atof(optarg); break;
		case 'N': d->c->n_hap = atoi(optarg); break;
		case 'r': d->c->het_rate = atof(optarg); break;
		case 'c': d->format |= BAM_PLF_CNS; break;
		case 'i': d->format |= BAM_PLF_INDEL_ONLY; break;
		case 'm': d->mask = atoi(optarg); break;
		case 'g': d->format |= BAM_PLF_GLF; break;
		case 'I': d->ido->q_indel = atoi(optarg); break;
		case 'G': d->ido->r_indel = atof(optarg); break;
		default: fprintf(stderr, "Unrecognizd option '-%c'.\n", c); return 1;
		}
	}
	if (optind == argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:  samtools pileup [options] <in.bam>|<in.sam>\n\n");
		fprintf(stderr, "Option: -s        simple (yet incomplete) pileup format\n");
		fprintf(stderr, "        -i        only show lines/consensus with indels\n");
		fprintf(stderr, "        -m INT    filtering reads with bits in INT [%d]\n", d->mask);
		fprintf(stderr, "        -t FILE   list of reference sequences (assume the input is in SAM)\n");
		fprintf(stderr, "        -l FILE   list of sites at which pileup is output\n");
		fprintf(stderr, "        -f FILE   reference sequence in the FASTA format\n\n");
		fprintf(stderr, "        -c        output the maq consensus sequence\n");
		fprintf(stderr, "        -g        output in the GLFv3 format (suppressing -c/-i/-s)\n");
		fprintf(stderr, "        -T FLOAT  theta in maq consensus calling model (for -c/-g) [%f]\n", d->c->theta);
		fprintf(stderr, "        -N INT    number of haplotypes in the sample (for -c/-g) [%d]\n", d->c->n_hap);
		fprintf(stderr, "        -r FLOAT  prior of a difference between two haplotypes (for -c/-g) [%f]\n", d->c->het_rate);
		fprintf(stderr, "        -G FLOAT  prior of an indel between two haplotypes (for -c/-g) [%f]\n", d->ido->r_indel);
		fprintf(stderr, "        -I INT    phred prob. of an indel in sequencing/prep. (for -c/-g) [%d]\n", d->ido->q_indel);
		fprintf(stderr, "\n");
		free(fn_list); free(fn_fa); free(d);
		return 1;
	}
	if (fn_fa) d->fai = fai_load(fn_fa);
	free(fn_fa);
	if (d->format & (BAM_PLF_CNS|BAM_PLF_GLF)) bam_maqcns_prepare(d->c);
	if (d->format & BAM_PLF_GLF) {
		glf3_header_t *h;
		h = glf3_header_init();
		d->fp = bgzf_fdopen(fileno(stdout), "w");
		glf3_header_write(d->fp, h);
		glf3_header_destroy(h);
	}
	if (d->fai == 0 && (d->format & (BAM_PLF_CNS|BAM_PLF_INDEL_ONLY)))
		fprintf(stderr, "[bam_pileup] indels will not be called when -f is absent.\n");
	if (fn_list) { // the input is SAM
		tamFile fp;
		bam1_t *b;
		int ret;
		bam_plbuf_t *buf = bam_plbuf_init(pileup_func, d);
		bam_plbuf_set_mask(buf, d->mask);
		d->h = sam_header_read2(fn_list);
		if (fn_pos) d->hash = load_pos(fn_pos, d->h);
		fp = sam_open(argv[optind]);
		b = (bam1_t*)calloc(1, sizeof(bam1_t));
		while ((ret = sam_read1(fp, d->h, b)) >= 0)
			bam_plbuf_push(b, buf);
		bam_plbuf_push(0, buf);
		bam_plbuf_destroy(buf);
		bam_destroy1(b);
		sam_close(fp);
	} else { // the input is BAM
		bamFile fp;
		fp = (strcmp(argv[optind], "-") == 0)? bam_dopen(fileno(stdin), "r") : bam_open(argv[optind], "r");
		d->h = bam_header_read(fp);
		if (d->h == 0) {
			fprintf(stderr, "[bam_pileup] fail to read the BAM header. Abort!\n");
			return 1;
		}
		if (fn_pos) d->hash = load_pos(fn_pos, d->h);
		bam_pileup_file(fp, d->mask, pileup_func, d);
		bam_close(fp);
	}
	if (d->format & BAM_PLF_GLF) bgzf_close(d->fp);
	if (fn_pos) { // free the hash table
		khint_t k;
		for (k = kh_begin(d->hash); k < kh_end(d->hash); ++k)
			if (kh_exist(d->hash, k)) free(kh_val(d->hash, k));
		kh_destroy(64, d->hash);
	}
	free(fn_pos); free(fn_list);
	bam_header_destroy(d->h);
	if (d->fai) fai_destroy(d->fai);
	bam_maqcns_destroy(d->c);
	free(d->ido); free(d->ref); free(d);
	return 0;
}
