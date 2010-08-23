#include <unistd.h>
#include <stdlib.h>
#include <math.h>
#include <zlib.h>
#include <errno.h>
#include "bcf.h"
#include "prob1.h"
#include "kstring.h"

#include "khash.h"
KHASH_SET_INIT_INT64(set64)

#include "kseq.h"
KSTREAM_INIT(gzFile, gzread, 16384)

#define VC_NO_PL   1
#define VC_NO_GENO 2
#define VC_BCFOUT  4
#define VC_CALL    8
#define VC_VARONLY 16
#define VC_VCFIN   32
#define VC_UNCOMP  64

typedef struct {
	int flag, prior_type;
	char *fn_list, *prior_file;
	double theta, pref;
} viewconf_t;

khash_t(set64) *bcf_load_pos(const char *fn, bcf_hdr_t *_h)
{
	void *str2id;
	gzFile fp;
	kstream_t *ks;
	int ret, dret, lineno = 1;
	kstring_t *str;
	khash_t(set64) *hash = 0;

	hash = kh_init(set64);
	str2id = bcf_build_refhash(_h);
	str = calloc(1, sizeof(kstring_t));
	fp = strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
	ks = ks_init(fp);
	while (ks_getuntil(ks, 0, str, &dret) >= 0) {
		int tid = bcf_str2id(str2id, str->s);
		if (tid >= 0 && dret != '\n') {
			if (ks_getuntil(ks, 0, str, &dret) >= 0) {
				uint64_t x = (uint64_t)tid<<32 | (atoi(str->s) - 1);
				kh_put(set64, hash, x, &ret);
			} else break;
		} else fprintf(stderr, "[%s] %s is not a reference name (line %d).\n", __func__, str->s, lineno);
		if (dret != '\n') while ((dret = ks_getc(ks)) > 0 && dret != '\n');
		if (dret < 0) break;
		++lineno;
	}
	bcf_str2id_destroy(str2id);
	ks_destroy(ks);
	gzclose(fp);
	free(str->s); free(str);
	return hash;
}

static double test_hwe(const double g[3])
{
	extern double kf_gammaq(double p, double x);
	double fexp, chi2, f[3], n;
	int i;
	n = g[0] + g[1] + g[2];
	fexp = (2. * g[2] + g[1]) / (2. * n);
	if (fexp > 1. - 1e-10) fexp = 1. - 1e-10;
	if (fexp < 1e-10) fexp = 1e-10;
	f[0] = n * (1. - fexp) * (1. - fexp);
	f[1] = n * 2. * fexp * (1. - fexp);
	f[2] = n * fexp * fexp;
	for (i = 0, chi2 = 0.; i < 3; ++i)
		chi2 += (g[i] - f[i]) * (g[i] - f[i]) / f[i];
	return kf_gammaq(.5, chi2 / 2.);
}

extern double kt_fisher_exact(int n11, int n12, int n21, int n22, double *_left, double *_right, double *two);

static double test_fisher(bcf1_t *b, const char *key, int d[4], int is_single)
{
	double left, right, two;
	char *p;
	int i;
	if ((p = strstr(b->info, key)) == 0) return -1.;
	p += 4;
	for (i = 0; i < 4; ++i) {
		d[i] = strtol(p, &p, 10);
		if (d[i] == 0 && (errno == EINVAL || errno == ERANGE)) return -2.;
		++p;
	}
	kt_fisher_exact(d[0], d[1], d[2], d[3], &left, &right, &two);
	return is_single? right : two;
}

static void rm_info(int n_smpl, bcf1_t *b, const char *key)
{
	char *p, *q;
	if ((p = strstr(b->info, key)) == 0) return;
	for (q = p; *q && *q != ';'; ++q);
	if (p > b->info && *(p-1) == ';') --p;
	memmove(p, q, b->l_str - (q - b->str));
	b->l_str -= q - p;
	bcf_sync(n_smpl, b);
}

static int update_bcf1(int n_smpl, bcf1_t *b, const bcf_p1aux_t *pa, const bcf_p1rst_t *pr, double pref, int flag)
{
	kstring_t s;
	int d[4], is_var = (pr->p_ref < pref);
	double p_hwe, p_dp, p_ed, r = is_var? pr->p_ref : 1. - pr->p_ref;

	p_hwe = test_hwe(pr->g);
	p_ed = test_fisher(b, "ED4=", d, 1);
	p_dp = test_fisher(b, "DP4=", d, 0);
	rm_info(n_smpl, b, "ED4=");

	memset(&s, 0, sizeof(kstring_t));
	kputc('\0', &s); kputs(b->ref, &s); kputc('\0', &s);
	if (is_var) {
		kputs(b->alt, &s);
	}
	kputc('\0', &s); kputc('\0', &s);
	kputs(b->info, &s);
	if (b->info[0]) kputc(';', &s);
	ksprintf(&s, "AF1=%.3lf;AFE=%.3lf", 1.-pr->f_em, 1.-pr->f_exp);
	if (p_hwe <= .2) ksprintf(&s, ";GC=%.2lf,%.2lf,%.2lf;HWE=%.3lf", pr->g[2], pr->g[1], pr->g[0], p_hwe);
	if (p_dp >= 0. && p_dp <= .2) ksprintf(&s, ";TDP=%.3lf", p_dp);
	if (p_ed >= 0. && p_ed <= .2) ksprintf(&s, ";TED=%.3lf", p_ed);
	kputc('\0', &s);
	kputs(b->fmt, &s); kputc('\0', &s);
	free(b->str);
	b->m_str = s.m; b->l_str = s.l; b->str = s.s;
	b->qual = r < 1e-100? 99 : -3.434 * log(r);
	if (b->qual > 99) b->qual = 99;
	bcf_sync(n_smpl, b);
	return is_var;
}

int bcfview(int argc, char *argv[])
{
	bcf_t *bp, *bout = 0;
	bcf1_t *b;
	int c;
	uint64_t n_processed = 0;
	viewconf_t vc;
	bcf_p1aux_t *p1 = 0;
	bcf_hdr_t *h;
	int tid, begin, end;
	char moder[4], modew[4];
	khash_t(set64) *hash = 0;

	tid = begin = end = -1;
	memset(&vc, 0, sizeof(viewconf_t));
	vc.prior_type = -1; vc.theta = 1e-3; vc.pref = 0.9;
	while ((c = getopt(argc, argv, "l:cGvLbSuP:t:p:")) >= 0) {
		switch (c) {
		case 'l': vc.fn_list = strdup(optarg); break;
		case 'G': vc.flag |= VC_NO_GENO; break;
		case 'L': vc.flag |= VC_NO_PL; break;
		case 'b': vc.flag |= VC_BCFOUT; break;
		case 'S': vc.flag |= VC_VCFIN; break;
		case 'c': vc.flag |= VC_CALL; break;
		case 'v': vc.flag |= VC_VARONLY; break;
		case 'u': vc.flag |= VC_UNCOMP | VC_BCFOUT; break;
		case 't': vc.theta = atof(optarg); break;
		case 'p': vc.pref = atof(optarg); break;
		case 'P':
			if (strcmp(optarg, "full") == 0) vc.prior_type = MC_PTYPE_FULL;
			else if (strcmp(optarg, "cond2") == 0) vc.prior_type = MC_PTYPE_COND2;
			else if (strcmp(optarg, "flat") == 0) vc.prior_type = MC_PTYPE_FLAT;
			else vc.prior_file = strdup(optarg);
			break;
		}
	}
	if (argc == optind) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   bcftools view [options] <in.bcf> [reg]\n\n");
		fprintf(stderr, "Options: -c        SNP calling\n");
		fprintf(stderr, "         -b        output BCF instead of VCF\n");
		fprintf(stderr, "         -u        uncompressed BCF output\n");
		fprintf(stderr, "         -S        input is VCF\n");
		fprintf(stderr, "         -G        suppress all individual genotype information\n");
		fprintf(stderr, "         -L        discard the PL genotype field\n");
		fprintf(stderr, "         -v        output potential variant sites only\n");
		fprintf(stderr, "         -l FILE   list of sites to output [all sites]\n");
		fprintf(stderr, "         -t FLOAT  scaled mutation rate [%.4lg]\n", vc.theta);
		fprintf(stderr, "         -p FLOAT  variant if P(ref|D)<FLOAT [%.3lg]\n", vc.pref);
		fprintf(stderr, "         -P STR    type of prior: full, cond2, flat [full]\n");
		fprintf(stderr, "\n");
		return 1;
	}

	b = calloc(1, sizeof(bcf1_t));
	strcpy(moder, "r");
	if (!(vc.flag & VC_VCFIN)) strcat(moder, "b");
	strcpy(modew, "w");
	if (vc.flag & VC_BCFOUT) strcat(modew, "b");
	if (vc.flag & VC_UNCOMP) strcat(modew, "u");
	bp = vcf_open(argv[optind], moder);
	h = vcf_hdr_read(bp);
	bout = vcf_open("-", modew);
	vcf_hdr_write(bout, h);
	if (vc.flag & VC_CALL) {
		p1 = bcf_p1_init(h->n_smpl);
		if (vc.prior_file) {
			if (bcf_p1_read_prior(p1, vc.prior_file) < 0) {
				fprintf(stderr, "[%s] fail to read the prior AFS.\n", __func__);
				return 1;
			}
		} else bcf_p1_init_prior(p1, vc.prior_type, vc.theta);
	}
	if (vc.fn_list) hash = bcf_load_pos(vc.fn_list, h);
	if (optind + 1 < argc && !(vc.flag&VC_VCFIN)) {
		void *str2id = bcf_build_refhash(h);
		if (bcf_parse_region(str2id, argv[optind+1], &tid, &begin, &end) >= 0) {
			bcf_idx_t *idx;
			idx = bcf_idx_load(argv[optind]);
			if (idx) {
				uint64_t off;
				off = bcf_idx_query(idx, tid, begin, end);
				bgzf_seek(bp->fp, off, SEEK_SET);
				bcf_idx_destroy(idx);
			}
		}
	}
	while (vcf_read(bp, h, b) > 0) {
		if (hash) {
			uint64_t x = (uint64_t)b->tid<<32 | b->pos;
			khint_t k = kh_get(set64, hash, x);
			if (k == kh_end(hash)) continue;
		}
		if (tid >= 0) {
			int l = strlen(b->ref);
			l = b->pos + (l > 0? l : 1);
			if (b->tid != tid || b->pos >= end) break;
			if (!(l > begin && end > b->pos)) continue;
		}
		if (vc.flag & VC_CALL) {
			bcf_p1rst_t pr;
			bcf_p1_cal(b, p1, &pr);
			if ((n_processed + 1) % 50000 == 0) bcf_p1_dump_afs(p1);
			if (pr.p_ref >= vc.pref && (vc.flag & VC_VARONLY)) continue;
			update_bcf1(h->n_smpl, b, p1, &pr, vc.pref, vc.flag);
		}
		if (vc.flag & VC_NO_GENO) {
			b->n_gi = 0;
			b->fmt[0] = '\0';
		}
		vcf_write(bout, h, b);
		++n_processed;
	}
	if (vc.prior_file) free(vc.prior_file);
	if (vc.flag & VC_CALL) bcf_p1_dump_afs(p1);
	bcf_hdr_destroy(h);
	bcf_destroy(b);
	vcf_close(bp); vcf_close(bout);
	if (hash) kh_destroy(set64, hash);
	if (vc.fn_list) free(vc.fn_list);
	if (p1) bcf_p1_destroy(p1);
	return 0;
}

