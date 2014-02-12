#include <math.h>
#include <stdio.h>
#include <unistd.h>
#include <ctype.h>
#include <string.h>
#include <errno.h>
#include <sys/stat.h>
#include <getopt.h>
#include <htslib/sam.h>
#include <htslib/faidx.h>
#include <htslib/kstring.h>
#include <htslib/khash_str2int.h>
#include "sam_header.h"
#include "samtools.h"

static inline int printw(int c, FILE *fp)
{
	char buf[16];
	int l, x;
	if (c == 0) return fputc('0', fp);
	for (l = 0, x = c < 0? -c : c; x > 0; x /= 10) buf[l++] = x%10 + '0';
	if (c < 0) buf[l++] = '-';
	buf[l] = 0;
	for (x = 0; x < l/2; ++x) {
		int y = buf[x]; buf[x] = buf[l-1-x]; buf[l-1-x] = y;
	}
	fputs(buf, fp);
	return 0;
}

static inline void pileup_seq(const bam_pileup1_t *p, int pos, int ref_len, const char *ref)
{
	int j;
	if (p->is_head) {
		putchar('^');
		putchar(p->b->core.qual > 93? 126 : p->b->core.qual + 33);
	}
	if (!p->is_del) {
		int c = seq_nt16_str[bam_seqi(bam_get_seq(p->b), p->qpos)];
		if (ref) {
			int rb = pos < ref_len? ref[pos] : 'N';
			if (c == '=' || seq_nt16_table[c] == seq_nt16_table[rb]) c = bam_is_rev(p->b)? ',' : '.';
			else c = bam_is_rev(p->b)? tolower(c) : toupper(c);
		} else {
			if (c == '=') c = bam_is_rev(p->b)? ',' : '.';
			else c = bam_is_rev(p->b)? tolower(c) : toupper(c);
		}
		putchar(c);
	} else putchar(p->is_refskip? (bam_is_rev(p->b)? '<' : '>') : '*');
	if (p->indel > 0) {
		putchar('+'); printw(p->indel, stdout);
		for (j = 1; j <= p->indel; ++j) {
			int c = seq_nt16_str[bam_seqi(bam_get_seq(p->b), p->qpos + j)];
			putchar(bam_is_rev(p->b)? tolower(c) : toupper(c));
		}
	} else if (p->indel < 0) {
		printw(p->indel, stdout);
		for (j = 1; j <= -p->indel; ++j) {
			int c = (ref && (int)pos+j < ref_len)? ref[pos+j] : 'N';
			putchar(bam_is_rev(p->b)? tolower(c) : toupper(c));
		}
	}
	if (p->is_tail) putchar('$');
}

#include <assert.h>
#include "bam2bcf.h"
#include "sample.h"

#define MPLP_GLF        1
#define MPLP_VCF        (1<<1)
#define MPLP_NO_COMP    (1<<2)
#define MPLP_NO_ORPHAN  (1<<3)
#define MPLP_REALN      (1<<4)
#define MPLP_NO_INDEL   (1<<5)
#define MPLP_REDO_BAQ   (1<<6)
#define MPLP_ILLUMINA13 (1<<7)
#define MPLP_IGNORE_RG  (1<<8)
#define MPLP_PRINT_POS  (1<<9)
#define MPLP_PRINT_MAPQ (1<<10)
#define MPLP_PER_SAMPLE (1<<11)
#define MPLP_SMART_OVERLAPS (1<<12)

void *bed_read(const char *fn);
void bed_destroy(void *_h);
int bed_overlap(const void *_h, const char *chr, int beg, int end);

typedef struct {
	int max_mq, min_mq, flag, min_baseQ, capQ_thres, max_depth, max_indel_depth, fmt_flag;
    int rflag_require, rflag_filter;
	int openQ, extQ, tandemQ, min_support; // for indels
	double min_frac; // for indels
	char *reg, *pl_list, *fai_fname;
	faidx_t *fai;
	void *bed, *rghash;
    int argc;
    char **argv;
} mplp_conf_t;

typedef struct {
	samFile *fp;
	hts_itr_t *iter;
	bam_hdr_t *h;
	int ref_id;
	char *ref;
	const mplp_conf_t *conf;
} mplp_aux_t;

typedef struct {
	int n;
	int *n_plp, *m_plp;
	bam_pileup1_t **plp;
} mplp_pileup_t;

static int mplp_func(void *data, bam1_t *b)
{
	extern int bam_realn(bam1_t *b, const char *ref);
	extern int bam_prob_realn_core(bam1_t *b, const char *ref, int);
	extern int bam_cap_mapQ(bam1_t *b, char *ref, int thres);
	mplp_aux_t *ma = (mplp_aux_t*)data;
	int ret, skip = 0;
	do {
		int has_ref;
		ret = ma->iter? sam_itr_next(ma->fp, ma->iter, b) : sam_read1(ma->fp, ma->h, b);
		if (ret < 0) break;
        // The 'B' cigar operation is not part of the specification, considering as obsolete.
		//  bam_remove_B(b);
		if (b->core.tid < 0 || (b->core.flag&BAM_FUNMAP)) { // exclude unmapped reads
			skip = 1;
			continue;
		}
        if (ma->conf->rflag_require && !(ma->conf->rflag_require&b->core.flag)) { skip = 1; continue; }
        if (ma->conf->rflag_filter && ma->conf->rflag_filter&b->core.flag) { skip = 1; continue; }
		if (ma->conf->bed) { // test overlap
			skip = !bed_overlap(ma->conf->bed, ma->h->target_name[b->core.tid], b->core.pos, bam_endpos(b));
			if (skip) continue;
		}
		if (ma->conf->rghash) { // exclude read groups
			uint8_t *rg = bam_aux_get(b, "RG");
			skip = (rg && khash_str2int_get(ma->conf->rghash, (const char*)(rg+1), NULL)==0);
			if (skip) continue;
		}
		if (ma->conf->flag & MPLP_ILLUMINA13) {
			int i;
			uint8_t *qual = bam_get_qual(b);
			for (i = 0; i < b->core.l_qseq; ++i)
				qual[i] = qual[i] > 31? qual[i] - 31 : 0;
		}
		has_ref = (ma->ref && ma->ref_id == b->core.tid)? 1 : 0;
		skip = 0;
		if (has_ref && (ma->conf->flag&MPLP_REALN)) bam_prob_realn_core(b, ma->ref, (ma->conf->flag & MPLP_REDO_BAQ)? 7 : 3);
		if (has_ref && ma->conf->capQ_thres > 10) {
			int q = bam_cap_mapQ(b, ma->ref, ma->conf->capQ_thres);
			if (q < 0) skip = 1;
			else if (b->core.qual > q) b->core.qual = q;
		}
		if (b->core.qual < ma->conf->min_mq) skip = 1; 
		else if ((ma->conf->flag&MPLP_NO_ORPHAN) && (b->core.flag&1) && !(b->core.flag&2)) skip = 1;
	} while (skip);
	return ret;
}

static void group_smpl(mplp_pileup_t *m, bam_sample_t *sm, kstring_t *buf,
					   int n, char *const*fn, int *n_plp, const bam_pileup1_t **plp, int ignore_rg)
{
	int i, j;
	memset(m->n_plp, 0, m->n * sizeof(int));
	for (i = 0; i < n; ++i) {
		for (j = 0; j < n_plp[i]; ++j) {
			const bam_pileup1_t *p = plp[i] + j;
			uint8_t *q;
			int id = -1;
			q = ignore_rg? 0 : bam_aux_get(p->b, "RG");
			if (q) id = bam_smpl_rg2smid(sm, fn[i], (char*)q+1, buf);
			if (id < 0) id = bam_smpl_rg2smid(sm, fn[i], 0, buf);
			if (id < 0 || id >= m->n) {
				assert(q); // otherwise a bug
				fprintf(stderr, "[%s] Read group %s used in file %s but absent from the header or an alignment missing read group.\n", __func__, (char*)q+1, fn[i]);
				exit(1);
			}
			if (m->n_plp[id] == m->m_plp[id]) {
				m->m_plp[id] = m->m_plp[id]? m->m_plp[id]<<1 : 8;
				m->plp[id] = realloc(m->plp[id], sizeof(bam_pileup1_t) * m->m_plp[id]);
			}
			m->plp[id][m->n_plp[id]++] = *p;
		}
	}
}

/*
 * Performs pileup
 * @param conf configuration for this pileup
 * @param n number of files specified in fn
 * @param fn filenames
 */
static int mpileup(mplp_conf_t *conf, int n, char **fn)
{
	extern void *bcf_call_add_rg(void *rghash, const char *hdtext, const char *list);
	extern void bcf_call_del_rghash(void *rghash);
	mplp_aux_t **data;
	int i, tid, pos, *n_plp, tid0 = -1, beg0 = 0, end0 = 1u<<29, ref_len, ref_tid = -1, max_depth, max_indel_depth;
	const bam_pileup1_t **plp;
	bam_mplp_t iter;
	bam_hdr_t *h = NULL; /* header of first file in input list */
	char *ref;
	void *rghash = NULL;

	bcf_callaux_t *bca = NULL;
	bcf_callret1_t *bcr = NULL;
	bcf_call_t bc;
	htsFile *bcf_fp = NULL;
	bcf_hdr_t *bcf_hdr = NULL;

	bam_sample_t *sm = NULL;
	kstring_t buf;
	mplp_pileup_t gplp;

	memset(&gplp, 0, sizeof(mplp_pileup_t));
	memset(&buf, 0, sizeof(kstring_t));
	memset(&bc, 0, sizeof(bcf_call_t));
	data = calloc(n, sizeof(mplp_aux_t*));
	plp = calloc(n, sizeof(bam_pileup1_t*));
	n_plp = calloc(n, sizeof(int));
	sm = bam_smpl_init();

    if (n == 0) {
        fprintf(stderr,"[%s] no input file/data given\n", __func__);
        exit(1);
    }

	// read the header of each file in the list and initialize data
	for (i = 0; i < n; ++i) {
		bam_hdr_t *h_tmp;
		data[i] = calloc(1, sizeof(mplp_aux_t));
		data[i]->fp = sam_open(fn[i], "rb");
		hts_set_fai_filename(data[i]->fp, conf->fai_fname);
        if ( !data[i]->fp )
        {
            fprintf(stderr, "[%s] failed to open %s: %s\n", __func__, fn[i], strerror(errno));
            exit(1);
        }
		data[i]->conf = conf;
		h_tmp = sam_hdr_read(data[i]->fp);
        if ( !h_tmp ) {
            fprintf(stderr,"[%s] fail to read the header of %s\n", __func__, fn[i]);
            exit(1);
        }
		data[i]->h = i? h : h_tmp; // for i==0, "h" has not been set yet
		bam_smpl_add(sm, fn[i], (conf->flag&MPLP_IGNORE_RG)? 0 : h_tmp->text);
        // Collect read group IDs with PL (platform) listed in pl_list (note: fragile, strstr search)
		rghash = bcf_call_add_rg(rghash, h_tmp->text, conf->pl_list);
		if (conf->reg) {
			hts_idx_t *idx = sam_index_load(data[i]->fp, fn[i]);
			if (idx == 0) {
				fprintf(stderr, "[%s] fail to load index for %s\n", __func__, fn[i]);
				exit(1);
			}
            if ( (data[i]->iter=sam_itr_querys(idx, data[i]->h, conf->reg)) == 0) {
                fprintf(stderr, "[E::%s] fail to parse region '%s'\n", __func__, conf->reg);
                exit(1);
            }
			if (i == 0) tid0 = data[i]->iter->tid, beg0 = data[i]->iter->beg, end0 = data[i]->iter->end;
			hts_idx_destroy(idx);
		}
		if (i == 0) h = h_tmp; /* save the header of first file in list */
		else {
			// FIXME: to check consistency
			bam_hdr_destroy(h_tmp);
		}
	}
	// allocate data storage proportionate to number of samples being studied sm->n
	gplp.n = sm->n;
	gplp.n_plp = calloc(sm->n, sizeof(int));
	gplp.m_plp = calloc(sm->n, sizeof(int));
	gplp.plp = calloc(sm->n, sizeof(bam_pileup1_t*));

	fprintf(stderr, "[%s] %d samples in %d input files\n", __func__, sm->n, n);
	// write the VCF header
	if (conf->flag & MPLP_GLF) 
    {
        if ( conf->flag & MPLP_VCF )
            bcf_fp = (conf->flag&MPLP_NO_COMP) ? hts_open("-","wu") : hts_open("-","wz");   // uncompressed VCF or compressed VCF
        else
            bcf_fp = (conf->flag&MPLP_NO_COMP) ? hts_open("-","wub") : hts_open("-","wb");  // uncompressed BCF or compressed BCF

		bcf_hdr = bcf_hdr_init("w"); 
        kstring_t str = {0,0,0};

        ksprintf(&str, "##samtoolsVersion=%s+htslib-%s\n",samtools_version(),hts_version());
        bcf_hdr_append(bcf_hdr, str.s);

        str.l = 0;
        ksprintf(&str, "##samtoolsCommand=samtools mpileup");
        for (i=1; i<conf->argc; i++) ksprintf(&str, " %s", conf->argv[i]);
        kputc('\n', &str);
        bcf_hdr_append(bcf_hdr, str.s);

        if (conf->fai_fname) 
        {
            str.l = 0;
            ksprintf(&str, "##reference=file://%s\n", conf->fai_fname);
            bcf_hdr_append(bcf_hdr, str.s);
        }

        // todo: use/write new BAM header manipulation routines, fill also UR, M5
        for (i=0; i<h->n_targets; i++)
        {
            str.l = 0;
            ksprintf(&str, "##contig=<ID=%s,length=%d>", h->target_name[i], h->target_len[i]);
            bcf_hdr_append(bcf_hdr, str.s);
        }
        free(str.s);
        bcf_hdr_append(bcf_hdr,"##INFO=<ID=INDEL,Number=0,Type=Flag,Description=\"Indicates that the variant is an INDEL.\">");
        bcf_hdr_append(bcf_hdr,"##INFO=<ID=IDV,Number=1,Type=Integer,Description=\"Maximum number of reads supporting an indel\">");
        bcf_hdr_append(bcf_hdr,"##INFO=<ID=IMF,Number=1,Type=Float,Description=\"Maximum fraction of reads supporting an indel\">");
        bcf_hdr_append(bcf_hdr,"##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Raw read depth\">");
        bcf_hdr_append(bcf_hdr,"##INFO=<ID=VDB,Number=1,Type=Float,Description=\"Variant Distance Bias for filtering splice-site artefacts in RNA-seq data (bigger is better)\",Version=3>");
        bcf_hdr_append(bcf_hdr,"##INFO=<ID=RPB,Number=1,Type=Float,Description=\"Mann-Whitney U test of Read Position Bias (bigger is better)\">");
        bcf_hdr_append(bcf_hdr,"##INFO=<ID=MQB,Number=1,Type=Float,Description=\"Mann-Whitney U test of Mapping Quality Bias (bigger is better)\">");
        bcf_hdr_append(bcf_hdr,"##INFO=<ID=BQB,Number=1,Type=Float,Description=\"Mann-Whitney U test of Base Quality Bias (bigger is better)\">");
        bcf_hdr_append(bcf_hdr,"##INFO=<ID=MQSB,Number=1,Type=Float,Description=\"Mann-Whitney U test of Mapping Quality vs Strand Bias (bigger is better)\">");
#if CDF_MWU_TESTS
        bcf_hdr_append(bcf_hdr,"##INFO=<ID=RPB2,Number=1,Type=Float,Description=\"Mann-Whitney U test of Read Position Bias [CDF] (bigger is better)\">");
        bcf_hdr_append(bcf_hdr,"##INFO=<ID=MQB2,Number=1,Type=Float,Description=\"Mann-Whitney U test of Mapping Quality Bias [CDF] (bigger is better)\">");
        bcf_hdr_append(bcf_hdr,"##INFO=<ID=BQB2,Number=1,Type=Float,Description=\"Mann-Whitney U test of Base Quality Bias [CDF] (bigger is better)\">");
        bcf_hdr_append(bcf_hdr,"##INFO=<ID=MQSB2,Number=1,Type=Float,Description=\"Mann-Whitney U test of Mapping Quality vs Strand Bias [CDF] (bigger is better)\">");
#endif
        bcf_hdr_append(bcf_hdr,"##INFO=<ID=SGB,Number=1,Type=Float,Description=\"Segregation based metric.\">");
        bcf_hdr_append(bcf_hdr,"##INFO=<ID=MQ0F,Number=1,Type=Float,Description=\"Fraction of MQ0 reads (smaller is better)\">");
        bcf_hdr_append(bcf_hdr,"##INFO=<ID=I16,Number=16,Type=Float,Description=\"Auxiliary tag used for calling, see description of bcf_callret1_t in bam2bcf.h\">");
        bcf_hdr_append(bcf_hdr,"##INFO=<ID=QS,Number=.,Type=Float,Description=\"Auxiliary tag used for calling\">");
        bcf_hdr_append(bcf_hdr,"##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"List of Phred-scaled genotype likelihoods\">");
        bcf_hdr_append(bcf_hdr,"##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Number of high-quality bases\">");
        bcf_hdr_append(bcf_hdr,"##FORMAT=<ID=DV,Number=1,Type=Integer,Description=\"Number of high-quality non-reference bases\">");

		for (i=0; i<sm->n; i++)
            bcf_hdr_add_sample(bcf_hdr, sm->smpl[i]);
        bcf_hdr_write(bcf_fp, bcf_hdr);

		bca = bcf_call_init(-1., conf->min_baseQ);
		bcr = calloc(sm->n, sizeof(bcf_callret1_t));
		bca->rghash = rghash;
		bca->openQ = conf->openQ, bca->extQ = conf->extQ, bca->tandemQ = conf->tandemQ;
		bca->min_frac = conf->min_frac;
		bca->min_support = conf->min_support;
        bca->per_sample_flt = conf->flag & MPLP_PER_SAMPLE;

        bc.bcf_hdr = bcf_hdr;
		bc.n = sm->n;
		bc.PL = malloc(15 * sm->n * sizeof(*bc.PL));
        if (conf->fmt_flag & B2B_FMT_DP) bc.DP = malloc(sm->n * sizeof(*bc.DP));
        if (conf->fmt_flag & B2B_FMT_DV) bc.DV = malloc(sm->n * sizeof(*bc.DV));
	}
	if (tid0 >= 0 && conf->fai) { // region is set
		ref = faidx_fetch_seq(conf->fai, h->target_name[tid0], 0, 0x7fffffff, &ref_len);
		ref_tid = tid0;
		for (i = 0; i < n; ++i) data[i]->ref = ref, data[i]->ref_id = tid0;
	} else ref_tid = -1, ref = 0;

	// begin pileup
	iter = bam_mplp_init(n, mplp_func, (void**)data);
    if ( conf->flag & MPLP_SMART_OVERLAPS ) bam_mplp_init_overlaps(iter);
	max_depth = conf->max_depth;
	if (max_depth * sm->n > 1<<20)
		fprintf(stderr, "(%s) Max depth is above 1M. Potential memory hog!\n", __func__);
	if (max_depth * sm->n < 8000) {
		max_depth = 8000 / sm->n;
		fprintf(stderr, "<%s> Set max per-file depth to %d\n", __func__, max_depth);
	}
	max_indel_depth = conf->max_indel_depth * sm->n;
	bam_mplp_set_maxcnt(iter, max_depth);
    bcf1_t *bcf_rec = bcf_init1();
    int ret;
	while ( (ret=bam_mplp_auto(iter, &tid, &pos, n_plp, plp)) > 0) {
		if (conf->reg && (pos < beg0 || pos >= end0)) continue; // out of the region requested
		if (conf->bed && tid >= 0 && !bed_overlap(conf->bed, h->target_name[tid], pos, pos+1)) continue;
		if (tid != ref_tid) {
			free(ref); ref = 0;
			if (conf->fai) ref = faidx_fetch_seq(conf->fai, h->target_name[tid], 0, 0x7fffffff, &ref_len);
			for (i = 0; i < n; ++i) data[i]->ref = ref, data[i]->ref_id = tid;
			ref_tid = tid;
		}
		if (conf->flag & MPLP_GLF) {
			int total_depth, _ref0, ref16;
			for (i = total_depth = 0; i < n; ++i) total_depth += n_plp[i];
			group_smpl(&gplp, sm, &buf, n, fn, n_plp, plp, conf->flag & MPLP_IGNORE_RG);
			_ref0 = (ref && pos < ref_len)? ref[pos] : 'N';
			ref16 = seq_nt16_table[_ref0];
            bcf_callaux_clean(bca);
			for (i = 0; i < gplp.n; ++i)
				bcf_call_glfgen(gplp.n_plp[i], gplp.plp[i], ref16, bca, bcr + i);
            bc.tid = tid; bc.pos = pos;
			bcf_call_combine(gplp.n, bcr, bca, ref16, &bc);
			bcf_clear1(bcf_rec);
			bcf_call2bcf(&bc, bcf_rec, bcr, conf->fmt_flag, 0, 0);
			bcf_write1(bcf_fp, bcf_hdr, bcf_rec);
			// call indels
			if (!(conf->flag&MPLP_NO_INDEL) && total_depth < max_indel_depth && bcf_call_gap_prep(gplp.n, gplp.n_plp, gplp.plp, pos, bca, ref, rghash) >= 0) 
            {
                bcf_callaux_clean(bca);
				for (i = 0; i < gplp.n; ++i)
					bcf_call_glfgen(gplp.n_plp[i], gplp.plp[i], -1, bca, bcr + i);
				if (bcf_call_combine(gplp.n, bcr, bca, -1, &bc) >= 0) {
					bcf_clear1(bcf_rec);
					bcf_call2bcf(&bc, bcf_rec, bcr, conf->fmt_flag, bca, ref);
					bcf_write1(bcf_fp, bcf_hdr, bcf_rec);
				}
			}
		} else {
			printf("%s\t%d\t%c", h->target_name[tid], pos + 1, (ref && pos < ref_len)? ref[pos] : 'N');
			for (i = 0; i < n; ++i) {
				int j, cnt;
				for (j = cnt = 0; j < n_plp[i]; ++j) {
					const bam_pileup1_t *p = plp[i] + j;
					if (bam_get_qual(p->b)[p->qpos] >= conf->min_baseQ) ++cnt;
				}
				printf("\t%d\t", cnt);
				if (n_plp[i] == 0) {
					printf("*\t*"); // FIXME: printf() is very slow...
					if (conf->flag & MPLP_PRINT_POS) printf("\t*");
				} else {
					for (j = 0; j < n_plp[i]; ++j) {
						const bam_pileup1_t *p = plp[i] + j;
						if (bam_get_qual(p->b)[p->qpos] >= conf->min_baseQ)
							pileup_seq(plp[i] + j, pos, ref_len, ref);
					}
					putchar('\t');
					for (j = 0; j < n_plp[i]; ++j) {
						const bam_pileup1_t *p = plp[i] + j;
						int c = bam_get_qual(p->b)[p->qpos];
						if (c >= conf->min_baseQ) {
							c = c + 33 < 126? c + 33 : 126;
							putchar(c);
						}
					}
					if (conf->flag & MPLP_PRINT_MAPQ) {
						putchar('\t');
						for (j = 0; j < n_plp[i]; ++j) {
                            const bam_pileup1_t *p = plp[i] + j;
                            int c = bam_get_qual(p->b)[p->qpos];
                            if ( c < conf->min_baseQ ) continue;
							c = plp[i][j].b->core.qual + 33;
							if (c > 126) c = 126;
							putchar(c);
						}
					}
					if (conf->flag & MPLP_PRINT_POS) {
						putchar('\t');
						for (j = 0; j < n_plp[i]; ++j) {
							if (j > 0) putchar(',');
							printf("%d", plp[i][j].qpos + 1); // FIXME: printf() is very slow...
						}
					}
				}
			}
			putchar('\n');
		}
	}

	// clean up
    free(bc.tmp.s);
    bcf_destroy1(bcf_rec);
	if (bcf_fp)
    {
        hts_close(bcf_fp);
        bcf_hdr_destroy(bcf_hdr);
        bcf_call_destroy(bca);
        free(bc.PL);
        free(bc.DP);
        free(bc.DV);
        free(bcr);
    }
	bam_smpl_destroy(sm); free(buf.s);
	for (i = 0; i < gplp.n; ++i) free(gplp.plp[i]);
	free(gplp.plp); free(gplp.n_plp); free(gplp.m_plp);
	bcf_call_del_rghash(rghash);
	bam_mplp_destroy(iter);
	bam_hdr_destroy(h);
	for (i = 0; i < n; ++i) {
		sam_close(data[i]->fp);
		if (data[i]->iter) hts_itr_destroy(data[i]->iter);
		free(data[i]);
	}
	free(data); free(plp); free(ref); free(n_plp);
	return ret;
}

#define MAX_PATH_LEN 1024
int read_file_list(const char *file_list,int *n,char **argv[])
{
    char buf[MAX_PATH_LEN];
    int len, nfiles = 0;
    char **files = NULL;
    struct stat sb;

    *n = 0;
    *argv = NULL;

    FILE *fh = fopen(file_list,"r");
    if ( !fh )
    {
        fprintf(stderr,"%s: %s\n", file_list,strerror(errno));
        return 1;
    }

    files = calloc(nfiles,sizeof(char*));
    nfiles = 0;
    while ( fgets(buf,MAX_PATH_LEN,fh) ) 
    {
        // allow empty lines and trailing spaces
        len = strlen(buf);
        while ( len>0 && isspace(buf[len-1]) ) len--;
        if ( !len ) continue;

        // check sanity of the file list
        buf[len] = 0;
        if (stat(buf, &sb) != 0)
        {
            // no such file, check if it is safe to print its name
            int i, safe_to_print = 1;
            for (i=0; i<len; i++)
                if (!isprint(buf[i])) { safe_to_print = 0; break; } 
            if ( safe_to_print )
                fprintf(stderr,"The file list \"%s\" appears broken, could not locate: %s\n", file_list,buf);
            else
                fprintf(stderr,"Does the file \"%s\" really contain a list of files and do all exist?\n", file_list);
            return 1;
        }

        nfiles++;
        files = realloc(files,nfiles*sizeof(char*));
        files[nfiles-1] = strdup(buf);
    }
    fclose(fh);
    if ( !nfiles )
    {
        fprintf(stderr,"No files read from %s\n", file_list);
        return 1;
    }
    *argv = files;
    *n    = nfiles;
    return 0;
}
#undef MAX_PATH_LEN

int bam_mpileup(int argc, char *argv[])
{
	int c;
    const char *file_list = NULL;
    char **fn = NULL;
    int nfiles = 0, use_orphan = 0;
	mplp_conf_t mplp;
	memset(&mplp, 0, sizeof(mplp_conf_t));
	mplp.max_mq = 60;
	mplp.min_baseQ = 13;
	mplp.capQ_thres = 0;
	mplp.max_depth = 250; mplp.max_indel_depth = 250;
	mplp.openQ = 40; mplp.extQ = 20; mplp.tandemQ = 100;
	mplp.min_frac = 0.002; mplp.min_support = 1;
	mplp.flag = MPLP_NO_ORPHAN | MPLP_REALN | MPLP_SMART_OVERLAPS;
    mplp.argc = argc; mplp.argv = argv;
    mplp.rflag_filter = BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP;
    static struct option lopts[] = 
    {
        {"rf",1,0,1},   // require flag
        {"ff",1,0,2},   // filter flag
        {0,0,0,0}
    };
	while ((c = getopt_long(argc, argv, "Agf:r:l:M:q:Q:uaRC:BDSd:L:b:P:po:e:h:Im:F:EG:6OsV1:2:vx",lopts,NULL)) >= 0) {
		switch (c) {
        case 'x': mplp.flag &= ~MPLP_SMART_OVERLAPS; break;
        case  1 : 
            mplp.rflag_require = bam_str2flag(optarg); 
            if ( mplp.rflag_require<0 ) { fprintf(stderr,"Could not parse --rf %s\n", optarg); return 1; }
            break;
        case  2 : 
            mplp.rflag_filter = bam_str2flag(optarg); 
            if ( mplp.rflag_filter<0 ) { fprintf(stderr,"Could not parse --ff %s\n", optarg); return 1; }
            break;
		case 'f':
			mplp.fai = fai_load(optarg);
			if (mplp.fai == 0) return 1;
            mplp.fai_fname = optarg;
			break;
		case 'd': mplp.max_depth = atoi(optarg); break;
		case 'r': mplp.reg = strdup(optarg); break;
        case 'l': 
                  // In the original version the whole BAM was streamed which is inefficient
                  //  with few BED intervals and big BAMs. Todo: devise a heuristic to determine 
                  //  best strategy, that is streaming or jumping.
                  mplp.bed = bed_read(optarg); break;
		case 'P': mplp.pl_list = strdup(optarg); break;
		case 'p': mplp.flag |= MPLP_PER_SAMPLE; break;
		case 'g': mplp.flag |= MPLP_GLF; break;
		case 'v': mplp.flag |= MPLP_GLF | MPLP_VCF; break;
		case 'u': mplp.flag |= MPLP_NO_COMP | MPLP_GLF; break;
		case 'a': mplp.flag |= MPLP_NO_ORPHAN | MPLP_REALN; break;
		case 'B': mplp.flag &= ~MPLP_REALN; break;
		case 'D': mplp.fmt_flag |= B2B_FMT_DP; break;
		case 'S': mplp.fmt_flag |= B2B_FMT_SP; break;
		case 'V': mplp.fmt_flag |= B2B_FMT_DV; break;
		case 'I': mplp.flag |= MPLP_NO_INDEL; break;
		case 'E': mplp.flag |= MPLP_REDO_BAQ; break;
		case '6': mplp.flag |= MPLP_ILLUMINA13; break;
		case 'R': mplp.flag |= MPLP_IGNORE_RG; break;
		case 's': mplp.flag |= MPLP_PRINT_MAPQ; break;
		case 'O': mplp.flag |= MPLP_PRINT_POS; break;
		case 'C': mplp.capQ_thres = atoi(optarg); break;
		case 'M': mplp.max_mq = atoi(optarg); break;
		case 'q': mplp.min_mq = atoi(optarg); break;
		case 'Q': mplp.min_baseQ = atoi(optarg); break;
        case 'b': file_list = optarg; break;
		case 'o': mplp.openQ = atoi(optarg); break;
		case 'e': mplp.extQ = atoi(optarg); break;
		case 'h': mplp.tandemQ = atoi(optarg); break;
		case 'A': use_orphan = 1; break;
		case 'F': mplp.min_frac = atof(optarg); break;
		case 'm': mplp.min_support = atoi(optarg); break;
		case 'L': mplp.max_indel_depth = atoi(optarg); break;
		case 'G': {
				FILE *fp_rg;
				char buf[1024];
				mplp.rghash = khash_str2int_init();
				if ((fp_rg = fopen(optarg, "r")) == 0)
					fprintf(stderr, "(%s) Fail to open file %s. Continue anyway.\n", __func__, optarg);
				while (!feof(fp_rg) && fscanf(fp_rg, "%s", buf) > 0) // this is not a good style, but forgive me...
					khash_str2int_inc(mplp.rghash, strdup(buf));
				fclose(fp_rg);
			}
			break;
		}
	}
	if (use_orphan) mplp.flag &= ~MPLP_NO_ORPHAN;
	if (argc == 1) 
    {
        char *tmp_require = bam_flag2str(mplp.rflag_require);
        char *tmp_filter  = bam_flag2str(mplp.rflag_filter);
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage: samtools mpileup [options] in1.bam [in2.bam [...]]\n\n");
		fprintf(stderr, "Input options:\n\n");
		fprintf(stderr, "       -6              assume the quality is in the Illumina-1.3+ encoding\n");
		fprintf(stderr, "       -A              count anomalous read pairs\n");
		fprintf(stderr, "       -B              disable BAQ computation\n");
		fprintf(stderr, "       -b FILE         list of input BAM filenames, one per line [null]\n");
		fprintf(stderr, "       -C INT          parameter for adjusting mapQ; 0 to disable [0]\n");
		fprintf(stderr, "       -d INT          max per-BAM depth to avoid excessive memory usage [%d]\n", mplp.max_depth);
		fprintf(stderr, "       -E              recalculate extended BAQ on the fly thus ignoring existing BQs\n");
		fprintf(stderr, "       -f FILE         faidx indexed reference sequence file [null]\n");
		fprintf(stderr, "       -G FILE         exclude read groups listed in FILE [null]\n");
		fprintf(stderr, "       -l FILE         list of positions (chr pos) or regions (BED) [null]\n");
		fprintf(stderr, "       -M INT          cap mapping quality at INT [%d]\n", mplp.max_mq);
		fprintf(stderr, "       -r STR          region in which pileup is generated [null]\n");
		fprintf(stderr, "       -R              ignore RG tags\n");
		fprintf(stderr, "       -q INT          skip alignments with mapQ smaller than INT [%d]\n", mplp.min_mq);
		fprintf(stderr, "       -Q INT          skip bases with baseQ/BAQ smaller than INT [%d]\n", mplp.min_baseQ);
		fprintf(stderr, "       --rf STR|INT    required flags: skip reads with mask bits unset [%s]\n", tmp_require);
		fprintf(stderr, "       --ff STR|INT    filter flags: skip reads with mask bits set [%s]\n", tmp_filter);
		fprintf(stderr, "       -x              disable read-pair overlap detection\n");
		fprintf(stderr, "\nOutput options:\n\n");
		fprintf(stderr, "       -D/V            output per-sample DP/DV in BCF (requires -g/-v)\n");
		fprintf(stderr, "       -g/v            generate BCF/VCF output (genotype likelihoods)\n");
		fprintf(stderr, "       -O              output base positions on reads (disabled by -g/-u)\n");
		fprintf(stderr, "       -s              output mapping quality (disabled by -g/-u)\n");
		fprintf(stderr, "       -S              output per-sample strand bias P-value in BCF (require -g/-u)\n");
		fprintf(stderr, "       -u              generate uncompressed BCF/VCF output\n");
		fprintf(stderr, "\nSNP/INDEL genotype likelihoods options (effective with `-g' or `-u'):\n\n");
		fprintf(stderr, "       -e INT          Phred-scaled gap extension seq error probability [%d]\n", mplp.extQ);
		fprintf(stderr, "       -F FLOAT        minimum fraction of gapped reads for candidates [%g]\n", mplp.min_frac);
		fprintf(stderr, "       -h INT          coefficient for homopolymer errors [%d]\n", mplp.tandemQ);
		fprintf(stderr, "       -I              do not perform indel calling\n");
		fprintf(stderr, "       -L INT          max per-sample depth for INDEL calling [%d]\n", mplp.max_indel_depth);
		fprintf(stderr, "       -m INT          minimum gapped reads for indel candidates [%d]\n", mplp.min_support);
		fprintf(stderr, "       -o INT          Phred-scaled gap open sequencing error probability [%d]\n", mplp.openQ);
		fprintf(stderr, "       -p              apply -m and -F per-sample to increase sensitivity\n");
		fprintf(stderr, "       -P STR          comma separated list of platforms for indels [all]\n");
		fprintf(stderr, "\n");
		fprintf(stderr, "Notes: Assuming diploid individuals.\n\n");
        free(tmp_require); free(tmp_filter);
		return 1;
	}
    int ret;
    if (file_list) {
        if ( read_file_list(file_list,&nfiles,&fn) ) return 1;
        ret = mpileup(&mplp,nfiles,fn);
        for (c=0; c<nfiles; c++) free(fn[c]);
        free(fn);
    } 
    else 
        ret = mpileup(&mplp, argc - optind, argv + optind);
	if (mplp.rghash) khash_str2int_destroy_free(mplp.rghash);
	free(mplp.reg); free(mplp.pl_list);
	if (mplp.fai) fai_destroy(mplp.fai);
	if (mplp.bed) bed_destroy(mplp.bed);
	return ret;
}
