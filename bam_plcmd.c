#include <stdio.h>
#include <unistd.h>
#include <ctype.h>
#include "bam.h"
#include "faidx.h"
#include "bam_maqcns.h"
#include "khash.h"
#include "glf.h"
KHASH_SET_INIT_INT64(64)

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
	int tid, len;
	int mask;
	char *ref;
	glfFile fp; // for glf output only
} pu_data_t;

char **bam_load_pos(const char *fn, int *_n);
void bam_init_header_hash(bam_header_t *header);
int32_t bam_get_tid(const bam_header_t *header, const char *seq_name);

static khash_t(64) *load_pos(const char *fn, bam_header_t *h)
{
	int n, tmp, i;
	char **list, *s;
	uint64_t x;
	khash_t(64) *hash;
	bam_init_header_hash(h);
	list = bam_load_pos(fn, &n);
	hash = kh_init(64);
	for (i = 0; i < n; ++i) {
		x = (uint64_t)bam_get_tid(h, list[i]) << 32;
		s = list[i];
		while (*s++);
		x |= *((uint32_t*)s) - 1;
		kh_put(64, hash, x, &tmp);
		free(list[i]);
	}
	free(list);
	return hash;
}

// an analogy to pileup_func() below
static int glt_func(uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pu, void *data)
{
	pu_data_t *d = (pu_data_t*)data;
	bam_maqindel_ret_t *r = 0;
	int rb;
	glf1_t *g;
	glf2_t g2;
	if (d->fai == 0) {
		fprintf(stderr, "[glt_func] reference sequence is required for generating GLT. Abort!\n");
		exit(1);
	}
	if (d->hash && kh_get(64, d->hash, (uint64_t)tid<<32|pos) == kh_end(d->hash)) return 0;
	if (d->fai && (int)tid != d->tid) {
		if (d->ref) { // then write the end mark
			memset(&g2, 0, sizeof(glf2_t));
			g2.type = GLF_TYPE_END;
			bgzf_write(d->fp, &g2, sizeof(glf2_t));
		}
		glf_ref_write(d->fp, d->h->target_name[tid]); // write reference
		free(d->ref);
		d->ref = fai_fetch(d->fai, d->h->target_name[tid], &d->len);
		d->tid = tid;
	}
	rb = (d->ref && (int)pos < d->len)? d->ref[pos] : 'N';
	g = bam_maqcns_glfgen(n, pu, bam_nt16_table[rb], d->c);
	memcpy(&g2, g, sizeof(glf1_t));
	g2.type = GLF_TYPE_NORMAL;
	g2.pos = pos;
	bgzf_write(d->fp, &g2, sizeof(glf2_t));
	r = bam_maqindel(n, pos, d->ido, pu, d->ref);
	if (r) { // then write indel line
		int g3 = 3 * n, min;
		min = g3;
		if (min > r->gl[0]) min = r->gl[0];
		if (min > r->gl[1]) min = r->gl[1];
		g3 -= min;
		g2.ref_base = 0;
		g2.type = GLF_TYPE_INDEL;
		memset(g2.lk, 0, 10);
		g2.lk[0] = r->gl[0] < 255? r->gl[0] : 255;
		g2.lk[1] = r->gl[1] < 255? r->gl[1] : 255;
		g2.lk[2] = g3 < 255? g3 : 255;
		*(int16_t*)(g2.lk + 3) = (int16_t)r->indel1;
		*(int16_t*)(g2.lk + 5) = (int16_t)r->indel2;
		g2.min_lk = min < 255? min : 255;
		bgzf_write(d->fp, &g2, sizeof(glf2_t));
		if (r->indel1) bgzf_write(d->fp, r->s[0]+1, r->indel1>0? r->indel1 : -r->indel1);
		if (r->indel2) bgzf_write(d->fp, r->s[1]+1, r->indel2>0? r->indel2 : -r->indel2);
		bam_maqindel_ret_destroy(r);
	}
	free(g);
	return 0;
}

static int pileup_func(uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pu, void *data)
{
	pu_data_t *d = (pu_data_t*)data;
	bam_maqindel_ret_t *r = 0;
	int i, j, rb;
	uint32_t x;
	if (d->hash && kh_get(64, d->hash, (uint64_t)tid<<32|pos) == kh_end(d->hash)) return 0;
	if (d->format & BAM_PLF_GLF) return glt_func(tid, pos, n, pu, data);
	if (d->fai && (int)tid != d->tid) {
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
	if (d->format & (BAM_PLF_CNS|BAM_PLF_INDEL_ONLY)) {
		if (d->ref) // indel calling
			r = bam_maqindel(n, pos, d->ido, pu, d->ref);
	}
	// pileup strings
	printf("%d\t", n);
	for (i = 0; i < n; ++i) {
		const bam_pileup1_t *p = pu + i;
		if (p->is_head) printf("^%c", p->b->core.qual > 93? 126 : p->b->core.qual + 33);
		if (!p->is_del) {
			int c = bam_nt16_rev_table[bam1_seqi(bam1_seq(p->b), p->qpos)];
			if (toupper(c) == toupper(rb)) c = bam1_strand(p->b)? ',' : '.';
			else bam1_strand(p->b)? tolower(c) : toupper(c);
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
		printf("%s\t%s\t", r->s[0], r->s[1]);
		//printf("%d\t%d\t", r->gl[0], r->gl[1]);
		printf("%d\t%d\t%d\t%d\n", r->cnt1, r->cnt2, r->cnt_ambi, r->cnt_anti);
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
		case 'r': d->c->het_rate = atoi(optarg); break;
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
		fprintf(stderr, "Usage:  bamtk pileup [options] <in.bam>|<in.sam>\n\n");
		fprintf(stderr, "Option: -s        simple (yet incomplete) pileup format\n");
		fprintf(stderr, "        -i        only show lines/consensus with indels\n");
		fprintf(stderr, "        -m INT    filtering reads with bits in INT [%d]\n", d->mask);
		fprintf(stderr, "        -t FILE   list of reference sequences (assume the input is in SAM)\n");
		fprintf(stderr, "        -l FILE   list of sites at which pileup is output\n");
		fprintf(stderr, "        -f FILE   reference sequence in the FASTA format\n\n");
		fprintf(stderr, "        -c        output the maq consensus sequence\n");
		fprintf(stderr, "        -g        output in the extended GLT format (suppressing -c/-i/-s)\n");
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
	bam_maqcns_prepare(d->c);
	if (d->format & BAM_PLF_GLF) {
		glf_header_t *h;
		h = glf_header_init();
		d->fp = bgzf_fdopen(fileno(stdout), "w");
		glf_header_write(d->fp, h);
		glf_header_destroy(h);
	}
	if (fn_list) {
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
	} else {
		bamFile fp;
		fp = (strcmp(argv[optind], "-") == 0)? bam_dopen(fileno(stdin), "r") : bam_open(argv[optind], "r");
		d->h = bam_header_read(fp);
		if (fn_pos) d->hash = load_pos(fn_pos, d->h);
		bam_pileup_file(fp, d->mask, pileup_func, d);
		bam_close(fp);
	}
	if (d->format & BAM_PLF_GLF) bgzf_close(d->fp);
	free(fn_pos); free(fn_list);
	kh_destroy(64, d->hash);
	bam_header_destroy(d->h);
	if (d->fai) fai_destroy(d->fai);
	bam_maqcns_destroy(d->c);
	free(d->ido); free(d->ref); free(d);
	return 0;
}
