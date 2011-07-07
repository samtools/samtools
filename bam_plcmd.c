#include <math.h>
#include <stdio.h>
#include <unistd.h>
#include <ctype.h>
#include <string.h>
#include <errno.h>
#include "sam.h"
#include "faidx.h"
#include "kstring.h"

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
		int c = bam_nt16_rev_table[bam1_seqi(bam1_seq(p->b), p->qpos)];
		if (ref) {
			int rb = pos < ref_len? ref[pos] : 'N';
			if (c == '=' || bam_nt16_table[c] == bam_nt16_table[rb]) c = bam1_strand(p->b)? ',' : '.';
			else c = bam1_strand(p->b)? tolower(c) : toupper(c);
		} else {
			if (c == '=') c = bam1_strand(p->b)? ',' : '.';
			else c = bam1_strand(p->b)? tolower(c) : toupper(c);
		}
		putchar(c);
	} else putchar(p->is_refskip? (bam1_strand(p->b)? '<' : '>') : '*');
	if (p->indel > 0) {
		putchar('+'); printw(p->indel, stdout);
		for (j = 1; j <= p->indel; ++j) {
			int c = bam_nt16_rev_table[bam1_seqi(bam1_seq(p->b), p->qpos + j)];
			putchar(bam1_strand(p->b)? tolower(c) : toupper(c));
		}
	} else if (p->indel < 0) {
		printw(p->indel, stdout);
		for (j = 1; j <= -p->indel; ++j) {
			int c = (ref && (int)pos+j < ref_len)? ref[pos+j] : 'N';
			putchar(bam1_strand(p->b)? tolower(c) : toupper(c));
		}
	}
	if (p->is_tail) putchar('$');
}

#include <assert.h>
#include "bam2bcf.h"
#include "sample.h"

#define MPLP_GLF   0x10
#define MPLP_NO_COMP 0x20
#define MPLP_NO_ORPHAN 0x40
#define MPLP_REALN   0x80
#define MPLP_FMT_DP 0x100
#define MPLP_FMT_SP 0x200
#define MPLP_NO_INDEL 0x400
#define MPLP_EXT_BAQ 0x800
#define MPLP_ILLUMINA13 0x1000
#define MPLP_IGNORE_RG 0x2000
#define MPLP_PRINT_POS 0x4000
#define MPLP_PRINT_MAPQ 0x8000

void *bed_read(const char *fn);
void bed_destroy(void *_h);
int bed_overlap(const void *_h, const char *chr, int beg, int end);

typedef struct {
	int max_mq, min_mq, flag, min_baseQ, capQ_thres, max_depth, max_indel_depth;
	int openQ, extQ, tandemQ, min_support; // for indels
	double min_frac; // for indels
	char *reg, *pl_list;
	faidx_t *fai;
	void *bed, *rghash;
} mplp_conf_t;

typedef struct {
	bamFile fp;
	bam_iter_t iter;
	bam_header_t *h;
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
		ret = ma->iter? bam_iter_read(ma->fp, ma->iter, b) : bam_read1(ma->fp, b);
		if (ret < 0) break;
		if (b->core.tid < 0 || (b->core.flag&BAM_FUNMAP)) { // exclude unmapped reads
			skip = 1;
			continue;
		}
		if (ma->conf->bed) { // test overlap
			skip = !bed_overlap(ma->conf->bed, ma->h->target_name[b->core.tid], b->core.pos, bam_calend(&b->core, bam1_cigar(b)));
			if (skip) continue;
		}
		if (ma->conf->rghash) { // exclude read groups
			uint8_t *rg = bam_aux_get(b, "RG");
			skip = (rg && bcf_str2id(ma->conf->rghash, (const char*)(rg+1)) >= 0);
			if (skip) continue;
		}
		if (ma->conf->flag & MPLP_ILLUMINA13) {
			int i;
			uint8_t *qual = bam1_qual(b);
			for (i = 0; i < b->core.l_qseq; ++i)
				qual[i] = qual[i] > 31? qual[i] - 31 : 0;
		}
		has_ref = (ma->ref && ma->ref_id == b->core.tid)? 1 : 0;
		skip = 0;
		if (has_ref && (ma->conf->flag&MPLP_REALN)) bam_prob_realn_core(b, ma->ref, (ma->conf->flag & MPLP_EXT_BAQ)? 3 : 1);
		if (has_ref && ma->conf->capQ_thres > 10) {
			int q = bam_cap_mapQ(b, ma->ref, ma->conf->capQ_thres);
			if (q < 0) skip = 1;
			else if (b->core.qual > q) b->core.qual = q;
		}
		else if (b->core.qual < ma->conf->min_mq) skip = 1; 
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

static int mpileup(mplp_conf_t *conf, int n, char **fn)
{
	extern void *bcf_call_add_rg(void *rghash, const char *hdtext, const char *list);
	extern void bcf_call_del_rghash(void *rghash);
	mplp_aux_t **data;
	int i, tid, pos, *n_plp, tid0 = -1, beg0 = 0, end0 = 1u<<29, ref_len, ref_tid = -1, max_depth, max_indel_depth;
	const bam_pileup1_t **plp;
	bam_mplp_t iter;
	bam_header_t *h = 0;
	char *ref;
	void *rghash = 0;

	bcf_callaux_t *bca = 0;
	bcf_callret1_t *bcr = 0;
	bcf_call_t bc;
	bcf_t *bp = 0;
	bcf_hdr_t *bh = 0;

	bam_sample_t *sm = 0;
	kstring_t buf;
	mplp_pileup_t gplp;

	memset(&gplp, 0, sizeof(mplp_pileup_t));
	memset(&buf, 0, sizeof(kstring_t));
	memset(&bc, 0, sizeof(bcf_call_t));
	data = calloc(n, sizeof(void*));
	plp = calloc(n, sizeof(void*));
	n_plp = calloc(n, sizeof(int*));
	sm = bam_smpl_init();

	// read the header and initialize data
	for (i = 0; i < n; ++i) {
		bam_header_t *h_tmp;
		data[i] = calloc(1, sizeof(mplp_aux_t));
		data[i]->fp = strcmp(fn[i], "-") == 0? bam_dopen(fileno(stdin), "r") : bam_open(fn[i], "r");
		data[i]->conf = conf;
		h_tmp = bam_header_read(data[i]->fp);
		data[i]->h = i? h : h_tmp; // for i==0, "h" has not been set yet
		bam_smpl_add(sm, fn[i], (conf->flag&MPLP_IGNORE_RG)? 0 : h_tmp->text);
		rghash = bcf_call_add_rg(rghash, h_tmp->text, conf->pl_list);
		if (conf->reg) {
			int beg, end;
			bam_index_t *idx;
			idx = bam_index_load(fn[i]);
			if (idx == 0) {
				fprintf(stderr, "[%s] fail to load index for %d-th input.\n", __func__, i+1);
				exit(1);
			}
			if (bam_parse_region(h_tmp, conf->reg, &tid, &beg, &end) < 0) {
				fprintf(stderr, "[%s] malformatted region or wrong seqname for %d-th input.\n", __func__, i+1);
				exit(1);
			}
			if (i == 0) tid0 = tid, beg0 = beg, end0 = end;
			data[i]->iter = bam_iter_query(idx, tid, beg, end);
			bam_index_destroy(idx);
		}
		if (i == 0) h = h_tmp;
		else {
			// FIXME: to check consistency
			bam_header_destroy(h_tmp);
		}
	}
	gplp.n = sm->n;
	gplp.n_plp = calloc(sm->n, sizeof(int));
	gplp.m_plp = calloc(sm->n, sizeof(int));
	gplp.plp = calloc(sm->n, sizeof(void*));

	fprintf(stderr, "[%s] %d samples in %d input files\n", __func__, sm->n, n);
	// write the VCF header
	if (conf->flag & MPLP_GLF) {
		kstring_t s;
		bh = calloc(1, sizeof(bcf_hdr_t));
		s.l = s.m = 0; s.s = 0;
		bp = bcf_open("-", (conf->flag&MPLP_NO_COMP)? "wu" : "w");
		for (i = 0; i < h->n_targets; ++i) {
			kputs(h->target_name[i], &s);
			kputc('\0', &s);
		}
		bh->l_nm = s.l;
		bh->name = malloc(s.l);
		memcpy(bh->name, s.s, s.l);
		s.l = 0;
		for (i = 0; i < sm->n; ++i) {
			kputs(sm->smpl[i], &s); kputc('\0', &s);
		}
		bh->l_smpl = s.l;
		bh->sname = malloc(s.l);
		memcpy(bh->sname, s.s, s.l);
		bh->txt = malloc(strlen(BAM_VERSION) + 64);
		bh->l_txt = 1 + sprintf(bh->txt, "##samtoolsVersion=%s\n", BAM_VERSION);
		free(s.s);
		bcf_hdr_sync(bh);
		bcf_hdr_write(bp, bh);
		bca = bcf_call_init(-1., conf->min_baseQ);
		bcr = calloc(sm->n, sizeof(bcf_callret1_t));
		bca->rghash = rghash;
		bca->openQ = conf->openQ, bca->extQ = conf->extQ, bca->tandemQ = conf->tandemQ;
		bca->min_frac = conf->min_frac;
		bca->min_support = conf->min_support;
	}
	if (tid0 >= 0 && conf->fai) { // region is set
		ref = faidx_fetch_seq(conf->fai, h->target_name[tid0], 0, 0x7fffffff, &ref_len);
		ref_tid = tid0;
		for (i = 0; i < n; ++i) data[i]->ref = ref, data[i]->ref_id = tid0;
	} else ref_tid = -1, ref = 0;
	iter = bam_mplp_init(n, mplp_func, (void**)data);
	max_depth = conf->max_depth;
	if (max_depth * sm->n > 1<<20)
		fprintf(stderr, "(%s) Max depth is above 1M. Potential memory hog!\n", __func__);
	if (max_depth * sm->n < 8000) {
		max_depth = 8000 / sm->n;
		fprintf(stderr, "<%s> Set max per-file depth to %d\n", __func__, max_depth);
	}
	max_indel_depth = conf->max_indel_depth * sm->n;
	bam_mplp_set_maxcnt(iter, max_depth);
	while (bam_mplp_auto(iter, &tid, &pos, n_plp, plp) > 0) {
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
			bcf1_t *b = calloc(1, sizeof(bcf1_t));
			for (i = total_depth = 0; i < n; ++i) total_depth += n_plp[i];
			group_smpl(&gplp, sm, &buf, n, fn, n_plp, plp, conf->flag & MPLP_IGNORE_RG);
			_ref0 = (ref && pos < ref_len)? ref[pos] : 'N';
			ref16 = bam_nt16_table[_ref0];
			for (i = 0; i < gplp.n; ++i)
				bcf_call_glfgen(gplp.n_plp[i], gplp.plp[i], ref16, bca, bcr + i);
			bcf_call_combine(gplp.n, bcr, ref16, &bc);
			bcf_call2bcf(tid, pos, &bc, b, (conf->flag&(MPLP_FMT_DP|MPLP_FMT_SP))? bcr : 0,
						 (conf->flag&MPLP_FMT_SP), 0, 0);
			bcf_write(bp, bh, b);
			bcf_destroy(b);
			// call indels
			if (!(conf->flag&MPLP_NO_INDEL) && total_depth < max_indel_depth && bcf_call_gap_prep(gplp.n, gplp.n_plp, gplp.plp, pos, bca, ref, rghash) >= 0) {
				for (i = 0; i < gplp.n; ++i)
					bcf_call_glfgen(gplp.n_plp[i], gplp.plp[i], -1, bca, bcr + i);
				if (bcf_call_combine(gplp.n, bcr, -1, &bc) >= 0) {
					b = calloc(1, sizeof(bcf1_t));
					bcf_call2bcf(tid, pos, &bc, b, (conf->flag&(MPLP_FMT_DP|MPLP_FMT_SP))? bcr : 0,
								 (conf->flag&MPLP_FMT_SP), bca, ref);
					bcf_write(bp, bh, b);
					bcf_destroy(b);
				}
			}
		} else {
			printf("%s\t%d\t%c", h->target_name[tid], pos + 1, (ref && pos < ref_len)? ref[pos] : 'N');
			for (i = 0; i < n; ++i) {
				int j;
				printf("\t%d\t", n_plp[i]);
				if (n_plp[i] == 0) {
					printf("*\t*"); // FIXME: printf() is very slow...
					if (conf->flag & MPLP_PRINT_POS) printf("\t*");
				} else {
					for (j = 0; j < n_plp[i]; ++j)
						pileup_seq(plp[i] + j, pos, ref_len, ref);
					putchar('\t');
					for (j = 0; j < n_plp[i]; ++j) {
						const bam_pileup1_t *p = plp[i] + j;
						int c = bam1_qual(p->b)[p->qpos] + 33;
						if (c > 126) c = 126;
						putchar(c);
					}
					if (conf->flag & MPLP_PRINT_MAPQ) {
						putchar('\t');
						for (j = 0; j < n_plp[i]; ++j) {
							int c = plp[i][j].b->core.qual + 33;
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

	bcf_close(bp);
	bam_smpl_destroy(sm); free(buf.s);
	for (i = 0; i < gplp.n; ++i) free(gplp.plp[i]);
	free(gplp.plp); free(gplp.n_plp); free(gplp.m_plp);
	bcf_call_del_rghash(rghash);
	bcf_hdr_destroy(bh); bcf_call_destroy(bca); free(bc.PL); free(bcr);
	bam_mplp_destroy(iter);
	bam_header_destroy(h);
	for (i = 0; i < n; ++i) {
		bam_close(data[i]->fp);
		if (data[i]->iter) bam_iter_destroy(data[i]->iter);
		free(data[i]);
	}
	free(data); free(plp); free(ref); free(n_plp);
	return 0;
}

#define MAX_PATH_LEN 1024
static int read_file_list(const char *file_list,int *n,char **argv[])
{
    char buf[MAX_PATH_LEN];
    int len, nfiles;
    char **files;

    FILE *fh = fopen(file_list,"r");
    if ( !fh )
    {
        fprintf(stderr,"%s: %s\n", file_list,strerror(errno));
        return 1;
    }

    // Speed is not an issue here, determine the number of files by reading the file twice
    nfiles = 0;
    while ( fgets(buf,MAX_PATH_LEN,fh) ) nfiles++;

    if ( fseek(fh, 0L, SEEK_SET) )
    {
        fprintf(stderr,"%s: %s\n", file_list,strerror(errno));
        return 1;
    }

    files = calloc(nfiles,sizeof(char*));
    nfiles = 0;
    while ( fgets(buf,MAX_PATH_LEN,fh) ) 
    {
        len = strlen(buf);
        while ( len>0 && isspace(buf[len-1]) ) len--;
        if ( !len ) continue;

        files[nfiles] = malloc(sizeof(char)*(len+1)); 
        strncpy(files[nfiles],buf,len);
        files[nfiles][len] = 0;
        nfiles++;
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
	#define MPLP_PRINT_POS 0x4000
	mplp.max_mq = 60;
	mplp.min_baseQ = 13;
	mplp.capQ_thres = 0;
	mplp.max_depth = 250; mplp.max_indel_depth = 250;
	mplp.openQ = 40; mplp.extQ = 20; mplp.tandemQ = 100;
	mplp.min_frac = 0.002; mplp.min_support = 1;
	mplp.flag = MPLP_NO_ORPHAN | MPLP_REALN;
	while ((c = getopt(argc, argv, "Agf:r:l:M:q:Q:uaRC:BDSd:L:b:P:o:e:h:Im:F:EG:6Os")) >= 0) {
		switch (c) {
		case 'f':
			mplp.fai = fai_load(optarg);
			if (mplp.fai == 0) return 1;
			break;
		case 'd': mplp.max_depth = atoi(optarg); break;
		case 'r': mplp.reg = strdup(optarg); break;
		case 'l': mplp.bed = bed_read(optarg); break;
		case 'P': mplp.pl_list = strdup(optarg); break;
		case 'g': mplp.flag |= MPLP_GLF; break;
		case 'u': mplp.flag |= MPLP_NO_COMP | MPLP_GLF; break;
		case 'a': mplp.flag |= MPLP_NO_ORPHAN | MPLP_REALN; break;
		case 'B': mplp.flag &= ~MPLP_REALN; break;
		case 'D': mplp.flag |= MPLP_FMT_DP; break;
		case 'S': mplp.flag |= MPLP_FMT_SP; break;
		case 'I': mplp.flag |= MPLP_NO_INDEL; break;
		case 'E': mplp.flag |= MPLP_EXT_BAQ; break;
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
				mplp.rghash = bcf_str2id_init();
				if ((fp_rg = fopen(optarg, "r")) == 0)
					fprintf(stderr, "(%s) Fail to open file %s. Continue anyway.\n", __func__, optarg);
				while (!feof(fp_rg) && fscanf(fp_rg, "%s", buf) > 0) // this is not a good style, but forgive me...
					bcf_str2id_add(mplp.rghash, strdup(buf));
				fclose(fp_rg);
			}
			break;
		}
	}
	if (use_orphan) mplp.flag &= ~MPLP_NO_ORPHAN;
	if (argc == 1) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage: samtools mpileup [options] in1.bam [in2.bam [...]]\n\n");
		fprintf(stderr, "Input options:\n\n");
		fprintf(stderr, "       -6           assume the quality is in the Illumina-1.3+ encoding\n");
		fprintf(stderr, "       -A           count anomalous read pairs\n");
		fprintf(stderr, "       -B           disable BAQ computation\n");
		fprintf(stderr, "       -b FILE      list of input BAM files [null]\n");
		fprintf(stderr, "       -C INT       parameter for adjusting mapQ; 0 to disable [0]\n");
		fprintf(stderr, "       -d INT       max per-BAM depth to avoid excessive memory usage [%d]\n", mplp.max_depth);
		fprintf(stderr, "       -E           extended BAQ for higher sensitivity but lower specificity\n");
		fprintf(stderr, "       -f FILE      faidx indexed reference sequence file [null]\n");
		fprintf(stderr, "       -G FILE      exclude read groups listed in FILE [null]\n");
		fprintf(stderr, "       -l FILE      list of positions (chr pos) or regions (BED) [null]\n");
		fprintf(stderr, "       -M INT       cap mapping quality at INT [%d]\n", mplp.max_mq);
		fprintf(stderr, "       -r STR       region in which pileup is generated [null]\n");
		fprintf(stderr, "       -R           ignore RG tags\n");
		fprintf(stderr, "       -q INT       skip alignments with mapQ smaller than INT [%d]\n", mplp.min_mq);
		fprintf(stderr, "       -Q INT       skip bases with baseQ/BAQ smaller than INT [%d]\n", mplp.min_baseQ);
		fprintf(stderr, "\nOutput options:\n\n");
		fprintf(stderr, "       -D           output per-sample DP in BCF (require -g/-u)\n");
		fprintf(stderr, "       -g           generate BCF output (genotype likelihoods)\n");
		fprintf(stderr, "       -O           output base positions on reads (disabled by -g/-u)\n");
		fprintf(stderr, "       -s           output mapping quality (disabled by -g/-u)\n");
		fprintf(stderr, "       -S           output per-sample strand bias P-value in BCF (require -g/-u)\n");
		fprintf(stderr, "       -u           generate uncompress BCF output\n");
		fprintf(stderr, "\nSNP/INDEL genotype likelihoods options (effective with `-g' or `-u'):\n\n");
		fprintf(stderr, "       -e INT       Phred-scaled gap extension seq error probability [%d]\n", mplp.extQ);
		fprintf(stderr, "       -F FLOAT     minimum fraction of gapped reads for candidates [%g]\n", mplp.min_frac);
		fprintf(stderr, "       -h INT       coefficient for homopolymer errors [%d]\n", mplp.tandemQ);
		fprintf(stderr, "       -I           do not perform indel calling\n");
		fprintf(stderr, "       -L INT       max per-sample depth for INDEL calling [%d]\n", mplp.max_indel_depth);
		fprintf(stderr, "       -m INT       minimum gapped reads for indel candidates [%d]\n", mplp.min_support);
		fprintf(stderr, "       -o INT       Phred-scaled gap open sequencing error probability [%d]\n", mplp.openQ);
		fprintf(stderr, "       -P STR       comma separated list of platforms for indels [all]\n");
		fprintf(stderr, "\n");
		fprintf(stderr, "Notes: Assuming diploid individuals.\n\n");
		return 1;
	}
    if (file_list) {
        if ( read_file_list(file_list,&nfiles,&fn) ) return 1;
        mpileup(&mplp,nfiles,fn);
        for (c=0; c<nfiles; c++) free(fn[c]);
        free(fn);
    } else mpileup(&mplp, argc - optind, argv + optind);
	if (mplp.rghash) bcf_str2id_thorough_destroy(mplp.rghash);
	free(mplp.reg); free(mplp.pl_list);
	if (mplp.fai) fai_destroy(mplp.fai);
	if (mplp.bed) bed_destroy(mplp.bed);
	return 0;
}
