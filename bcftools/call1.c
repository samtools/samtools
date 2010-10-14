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

#define VC_NO_GENO 2
#define VC_BCFOUT  4
#define VC_CALL    8
#define VC_VARONLY 16
#define VC_VCFIN   32
#define VC_UNCOMP  64
#define VC_HWE     128
#define VC_KEEPALT 256
#define VC_ACGT_ONLY 512
#define VC_QCALL   1024
#define VC_CALL_GT 2048
#define VC_ADJLD   4096

typedef struct {
	int flag, prior_type, n1;
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

typedef struct {
	double p[4];
	int mq, depth, is_tested, d[4];
} anno16_t;

static double ttest(int n1, int n2, int a[4])
{
	extern double kf_betai(double a, double b, double x);
	double t, v, u1, u2;
	if (n1 == 0 || n2 == 0 || n1 + n2 < 3) return 1.0;
	u1 = (double)a[0] / n1; u2 = (double)a[2] / n2;
	if (u1 <= u2) return 1.;
	t = (u1 - u2) / sqrt(((a[1] - n1 * u1 * u1) + (a[3] - n2 * u2 * u2)) / (n1 + n2 - 2) * (1./n1 + 1./n2));
	v = n1 + n2 - 2;
//	printf("%d,%d,%d,%d,%lf,%lf,%lf\n", a[0], a[1], a[2], a[3], t, u1, u2);
	return t < 0.? 1. : .5 * kf_betai(.5*v, .5, v/(v+t*t));
}

static int test16_core(int anno[16], anno16_t *a)
{
	extern double kt_fisher_exact(int n11, int n12, int n21, int n22, double *_left, double *_right, double *two);
	double left, right;
	int i;
	a->p[0] = a->p[1] = a->p[2] = a->p[3] = 1.;
	memcpy(a->d, anno, 4 * sizeof(int));
	a->depth = anno[0] + anno[1] + anno[2] + anno[3];
	a->is_tested = (anno[0] + anno[1] > 0 && anno[2] + anno[3] > 0);
	if (a->depth == 0) return -1;
	a->mq = (int)(sqrt((anno[9] + anno[11]) / a->depth) + .499);
	kt_fisher_exact(anno[0], anno[1], anno[2], anno[3], &left, &right, &a->p[0]);
	for (i = 1; i < 4; ++i)
		a->p[i] = ttest(anno[0] + anno[1], anno[2] + anno[3], anno+4*i);
	return 0;
}

static int test16(bcf1_t *b, anno16_t *a)
{
	char *p;
	int i, anno[16];
	a->p[0] = a->p[1] = a->p[2] = a->p[3] = 1.;
	a->d[0] = a->d[1] = a->d[2] = a->d[3] = 0.;
	a->mq = a->depth = a->is_tested = 0;
	if ((p = strstr(b->info, "I16=")) == 0) return -1;
	p += 4;
	for (i = 0; i < 16; ++i) {
		anno[i] = strtol(p, &p, 10);
		if (anno[i] == 0 && (errno == EINVAL || errno == ERANGE)) return -2;
		++p;
	}
	return test16_core(anno, a);
}

static void rm_info(bcf1_t *b, const char *key)
{
	char *p, *q;
	if ((p = strstr(b->info, key)) == 0) return;
	for (q = p; *q && *q != ';'; ++q);
	if (p > b->info && *(p-1) == ';') --p;
	memmove(p, q, b->l_str - (q - b->str));
	b->l_str -= q - p;
	bcf_sync(b);
}

static int update_bcf1(int n_smpl, bcf1_t *b, const bcf_p1aux_t *pa, const bcf_p1rst_t *pr, double pref, int flag)
{
	kstring_t s;
	int is_var = (pr->p_ref < pref);
	double p_hwe, r = is_var? pr->p_ref : 1. - pr->p_ref;
	anno16_t a;

	p_hwe = pr->g[0] >= 0.? test_hwe(pr->g) : 1.0; // only do HWE g[] is calculated
	test16(b, &a);
	rm_info(b, "I16=");

	memset(&s, 0, sizeof(kstring_t));
	kputc('\0', &s); kputs(b->ref, &s); kputc('\0', &s);
	kputs(b->alt, &s); kputc('\0', &s); kputc('\0', &s);
	kputs(b->info, &s);
	if (b->info[0]) kputc(';', &s);
	ksprintf(&s, "AF1=%.3lf;AFE=%.3lf", 1.-pr->f_em, 1.-pr->f_exp);
	ksprintf(&s, ";DP4=%d,%d,%d,%d;MQ=%d", a.d[0], a.d[1], a.d[2], a.d[3], a.mq);
	if (a.is_tested) {
		if (pr->pc[0] >= 0.) ksprintf(&s, ";PC4=%lg,%lg,%lg,%lg", pr->pc[0], pr->pc[1], pr->pc[2], pr->pc[3]);
		ksprintf(&s, ";PV4=%.2lg,%.2lg,%.2lg,%.2lg", a.p[0], a.p[1], a.p[2], a.p[3]);
	}
	if (pr->g[0] >= 0. && p_hwe <= .2)
		ksprintf(&s, ";GC=%.2lf,%.2lf,%.2lf;HWE=%.3lf", pr->g[2], pr->g[1], pr->g[0], p_hwe);
	kputc('\0', &s);
	kputs(b->fmt, &s); kputc('\0', &s);
	free(b->str);
	b->m_str = s.m; b->l_str = s.l; b->str = s.s;
	b->qual = r < 1e-100? 99 : -4.343 * log(r);
	if (b->qual > 99) b->qual = 99;
	bcf_sync(b);
	if (!is_var) bcf_shrink_alt(b, 1);
	else if (!(flag&VC_KEEPALT))
		bcf_shrink_alt(b, pr->rank0 < 2? 2 : pr->rank0+1);
	if (is_var && (flag&VC_CALL_GT)) { // call individual genotype
		int i, x, old_n_gi = b->n_gi;
		s.m = b->m_str; s.l = b->l_str - 1; s.s = b->str;
		kputs(":GT:GQ", &s); kputc('\0', &s);
		b->m_str = s.m; b->l_str = s.l; b->str = s.s;
		bcf_sync(b);
		for (i = 0; i < b->n_smpl; ++i) {
			x = bcf_p1_call_gt(pa, pr->f_em, i);
			((uint8_t*)b->gi[old_n_gi].data)[i] = (x&3) == 0? 1<<3|1 : (x&3) == 1? 1 : 0;
			((uint8_t*)b->gi[old_n_gi+1].data)[i] = x>>2;
		}
	}
	return is_var;
}

double bcf_ld_freq(const bcf1_t *b0, const bcf1_t *b1, double f[4]);

int bcfview(int argc, char *argv[])
{
	extern int bcf_2qcall(bcf_hdr_t *h, bcf1_t *b);
	bcf_t *bp, *bout = 0;
	bcf1_t *b, *blast;
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
	vc.prior_type = vc.n1 = -1; vc.theta = 1e-3; vc.pref = 0.5;
	while ((c = getopt(argc, argv, "N1:l:cHAGvbSuP:t:p:QgL")) >= 0) {
		switch (c) {
		case '1': vc.n1 = atoi(optarg); break;
		case 'l': vc.fn_list = strdup(optarg); break;
		case 'N': vc.flag |= VC_ACGT_ONLY; break;
		case 'G': vc.flag |= VC_NO_GENO; break;
		case 'A': vc.flag |= VC_KEEPALT; break;
		case 'b': vc.flag |= VC_BCFOUT; break;
		case 'S': vc.flag |= VC_VCFIN; break;
		case 'c': vc.flag |= VC_CALL; break;
		case 'v': vc.flag |= VC_VARONLY | VC_CALL; break;
		case 'u': vc.flag |= VC_UNCOMP | VC_BCFOUT; break;
		case 'H': vc.flag |= VC_HWE; break;
		case 'g': vc.flag |= VC_CALL_GT | VC_CALL; break;
		case 't': vc.theta = atof(optarg); break;
		case 'p': vc.pref = atof(optarg); break;
		case 'Q': vc.flag |= VC_QCALL; break;
		case 'L': vc.flag |= VC_ADJLD; break;
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
		fprintf(stderr, "         -v        output potential variant sites only (force -c)\n");
		fprintf(stderr, "         -g        call genotypes at variant sites (force -c)\n");
		fprintf(stderr, "         -b        output BCF instead of VCF\n");
		fprintf(stderr, "         -u        uncompressed BCF output (force -b)\n");
		fprintf(stderr, "         -S        input is VCF\n");
		fprintf(stderr, "         -A        keep all possible alternate alleles at variant sites\n");
		fprintf(stderr, "         -G        suppress all individual genotype information\n");
		fprintf(stderr, "         -H        perform Hardy-Weinberg test (slower)\n");
		fprintf(stderr, "         -N        skip sites where REF is not A/C/G/T\n");
		fprintf(stderr, "         -Q        output the QCALL likelihood format\n");
		fprintf(stderr, "         -L        calculate LD for adjacent sites\n");
		fprintf(stderr, "         -1 INT    number of group-1 samples [0]\n");
		fprintf(stderr, "         -l FILE   list of sites to output [all sites]\n");
		fprintf(stderr, "         -t FLOAT  scaled mutation rate [%.4lg]\n", vc.theta);
		fprintf(stderr, "         -p FLOAT  variant if P(ref|D)<FLOAT [%.3lg]\n", vc.pref);
		fprintf(stderr, "         -P STR    type of prior: full, cond2, flat [full]\n");
		fprintf(stderr, "\n");
		return 1;
	}

	b = calloc(1, sizeof(bcf1_t));
	blast = calloc(1, sizeof(bcf1_t));
	strcpy(moder, "r");
	if (!(vc.flag & VC_VCFIN)) strcat(moder, "b");
	strcpy(modew, "w");
	if (vc.flag & VC_BCFOUT) strcat(modew, "b");
	if (vc.flag & VC_UNCOMP) strcat(modew, "u");
	bp = vcf_open(argv[optind], moder);
	h = vcf_hdr_read(bp);
	bout = vcf_open("-", modew);
	if (!(vc.flag & VC_QCALL)) vcf_hdr_write(bout, h);
	if (vc.flag & VC_CALL) {
		p1 = bcf_p1_init(h->n_smpl);
		if (vc.prior_file) {
			if (bcf_p1_read_prior(p1, vc.prior_file) < 0) {
				fprintf(stderr, "[%s] fail to read the prior AFS.\n", __func__);
				return 1;
			}
		} else bcf_p1_init_prior(p1, vc.prior_type, vc.theta);
		if (vc.n1 > 0) {
			bcf_p1_set_n1(p1, vc.n1);
			bcf_p1_init_subprior(p1, vc.prior_type, vc.theta);
		}
	}
	if (vc.fn_list) hash = bcf_load_pos(vc.fn_list, h);
	if (optind + 1 < argc && !(vc.flag&VC_VCFIN)) {
		void *str2id = bcf_build_refhash(h);
		if (bcf_parse_region(str2id, argv[optind+1], &tid, &begin, &end) >= 0) {
			bcf_idx_t *idx;
			idx = bcf_idx_load(argv[optind]);
			if (idx) {
				uint64_t off;
				off = bcf_idx_query(idx, tid, begin);
				if (off == 0) {
					fprintf(stderr, "[%s] no records in the query region.\n", __func__);
					return 1; // FIXME: a lot of memory leaks...
				}
				bgzf_seek(bp->fp, off, SEEK_SET);
				bcf_idx_destroy(idx);
			}
		}
	}
	while (vcf_read(bp, h, b) > 0) {
		if (vc.flag & VC_ACGT_ONLY) {
			int x;
			if (b->ref[0] == 0 || b->ref[1] != 0) continue;
			x = toupper(b->ref[0]);
			if (x != 'A' && x != 'C' && x != 'G' && x != 'T') continue;
		}
		if (hash) {
			uint64_t x = (uint64_t)b->tid<<32 | b->pos;
			khint_t k = kh_get(set64, hash, x);
			if (kh_size(hash) == 0) break;
			if (k == kh_end(hash)) continue;
			kh_del(set64, hash, k);
		}
		if (tid >= 0) {
			int l = strlen(b->ref);
			l = b->pos + (l > 0? l : 1);
			if (b->tid != tid || b->pos >= end) break;
			if (!(l > begin && end > b->pos)) continue;
		}
		++n_processed;
		if (vc.flag & VC_QCALL) { // output QCALL format; STOP here
			bcf_2qcall(h, b);
			continue;
		}
		if (vc.flag & (VC_CALL|VC_ADJLD)) bcf_gl2pl(b);
		if (vc.flag & VC_CALL) { // call variants
			bcf_p1rst_t pr;
			bcf_p1_cal(b, p1, &pr); // pr.g[3] is not calculated here
			if (vc.flag&VC_HWE) bcf_p1_cal_g3(p1, pr.g);
			if (n_processed % 100000 == 0) {
				fprintf(stderr, "[%s] %ld sites processed.\n", __func__, (long)n_processed);
				bcf_p1_dump_afs(p1);
			}
			if (pr.p_ref >= vc.pref && (vc.flag & VC_VARONLY)) continue;
			update_bcf1(h->n_smpl, b, p1, &pr, vc.pref, vc.flag);
		}
		if (vc.flag & VC_ADJLD) { // compute LD
			double f[4], r2;
			if ((r2 = bcf_ld_freq(blast, b, f)) >= 0) {
				kstring_t s;
				s.m = s.l = 0; s.s = 0;
				if (*b->info) kputc(';', &s);
				ksprintf(&s, "NEIR=%.3lf;NEIF=%.3lf,%.3lf", r2, f[0]+f[2], f[0]+f[1]);
				bcf_append_info(b, s.s, s.l);
				free(s.s);
			}
			bcf_cpy(blast, b);
		}
		if (vc.flag & VC_NO_GENO) { // do not output GENO fields
			b->n_gi = 0;
			b->fmt[0] = '\0';
		}
		vcf_write(bout, h, b);
	}
	if (vc.prior_file) free(vc.prior_file);
	if (vc.flag & VC_CALL) bcf_p1_dump_afs(p1);
	bcf_hdr_destroy(h);
	bcf_destroy(b); bcf_destroy(blast);
	vcf_close(bp); vcf_close(bout);
	if (hash) kh_destroy(set64, hash);
	if (vc.fn_list) free(vc.fn_list);
	if (p1) bcf_p1_destroy(p1);
	return 0;
}
