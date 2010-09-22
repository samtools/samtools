#include <unistd.h>
#include <assert.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "faidx.h"
#include "sam.h"
#include "kstring.h"
#include "kaln.h"

void bam_fillmd1_core(bam1_t *b, char *ref, int is_equal, int max_nm)
{
	uint8_t *seq = bam1_seq(b);
	uint32_t *cigar = bam1_cigar(b);
	bam1_core_t *c = &b->core;
	int i, x, y, u = 0;
	kstring_t *str;
	uint8_t *old_md, *old_nm;
	int32_t old_nm_i = -1, nm = 0;

	str = (kstring_t*)calloc(1, sizeof(kstring_t));
	for (i = y = 0, x = c->pos; i < c->n_cigar; ++i) {
		int j, l = cigar[i]>>4, op = cigar[i]&0xf;
		if (op == BAM_CMATCH) {
			for (j = 0; j < l; ++j) {
				int z = y + j;
				int c1 = bam1_seqi(seq, z), c2 = bam_nt16_table[(int)ref[x+j]];
				if (ref[x+j] == 0) break; // out of boundary
				if ((c1 == c2 && c1 != 15 && c2 != 15) || c1 == 0) { // a match
					if (is_equal) seq[z/2] &= (z&1)? 0xf0 : 0x0f;
					++u;
				} else {
					ksprintf(str, "%d", u);
					kputc(ref[x+j], str);
					u = 0; ++nm;
				}
			}
			if (j < l) break;
			x += l; y += l;
		} else if (op == BAM_CDEL) {
			ksprintf(str, "%d", u);
			kputc('^', str);
			for (j = 0; j < l; ++j) {
				if (ref[x+j] == 0) break;
				kputc(ref[x+j], str);
			}
			u = 0;
			if (j < l) break;
			x += l; nm += l;
		} else if (op == BAM_CINS || op == BAM_CSOFT_CLIP) {
			y += l;
			if (op == BAM_CINS) nm += l;
		} else if (op == BAM_CREF_SKIP) {
			x += l;
		}
	}
	ksprintf(str, "%d", u);
	// apply max_nm
	if (max_nm > 0 && nm >= max_nm) {
		for (i = y = 0, x = c->pos; i < c->n_cigar; ++i) {
			int j, l = cigar[i]>>4, op = cigar[i]&0xf;
			if (op == BAM_CMATCH) {
				for (j = 0; j < l; ++j) {
					int z = y + j;
					int c1 = bam1_seqi(seq, z), c2 = bam_nt16_table[(int)ref[x+j]];
					if (ref[x+j] == 0) break; // out of boundary
					if ((c1 == c2 && c1 != 15 && c2 != 15) || c1 == 0) { // a match
						seq[z/2] |= (z&1)? 0x0f : 0xf0;
						bam1_qual(b)[z] = 0;
					}
				}
				if (j < l) break;
				x += l; y += l;
			} else if (op == BAM_CDEL || op == BAM_CREF_SKIP) x += l;
			else if (op == BAM_CINS || op == BAM_CSOFT_CLIP) y += l;
		}
	}
	// update NM
	old_nm = bam_aux_get(b, "NM");
	if (c->flag & BAM_FUNMAP) return;
	if (old_nm) old_nm_i = bam_aux2i(old_nm);
	if (!old_nm) bam_aux_append(b, "NM", 'i', 4, (uint8_t*)&nm);
	else if (nm != old_nm_i) {
		fprintf(stderr, "[bam_fillmd1] different NM for read '%s': %d -> %d\n", bam1_qname(b), old_nm_i, nm);
		bam_aux_del(b, old_nm);
		bam_aux_append(b, "NM", 'i', 4, (uint8_t*)&nm);
	}
	// update MD
	old_md = bam_aux_get(b, "MD");
	if (!old_md) bam_aux_append(b, "MD", 'Z', str->l + 1, (uint8_t*)str->s);
	else {
		int is_diff = 0;
		if (strlen((char*)old_md+1) == str->l) {
			for (i = 0; i < str->l; ++i)
				if (toupper(old_md[i+1]) != toupper(str->s[i]))
					break;
			if (i < str->l) is_diff = 1;
		} else is_diff = 1;
		if (is_diff) {
			fprintf(stderr, "[bam_fillmd1] different MD for read '%s': '%s' -> '%s'\n", bam1_qname(b), old_md+1, str->s);
			bam_aux_del(b, old_md);
			bam_aux_append(b, "MD", 'Z', str->l + 1, (uint8_t*)str->s);
		}
	}
	free(str->s); free(str);
}

void bam_fillmd1(bam1_t *b, char *ref, int is_equal)
{
	bam_fillmd1_core(b, ref, is_equal, 0);
}

int bam_cap_mapQ(bam1_t *b, char *ref, int thres)
{
	uint8_t *seq = bam1_seq(b), *qual = bam1_qual(b);
	uint32_t *cigar = bam1_cigar(b);
	bam1_core_t *c = &b->core;
	int i, x, y, mm, q, len, clip_l, clip_q;
	double t;
	if (thres < 0) thres = 40; // set the default
	mm = q = len = clip_l = clip_q = 0;
	for (i = y = 0, x = c->pos; i < c->n_cigar; ++i) {
		int j, l = cigar[i]>>4, op = cigar[i]&0xf;
		if (op == BAM_CMATCH) {
			for (j = 0; j < l; ++j) {
				int z = y + j;
				int c1 = bam1_seqi(seq, z), c2 = bam_nt16_table[(int)ref[x+j]];
				if (ref[x+j] == 0) break; // out of boundary
				if (c2 != 15 && c1 != 15 && qual[z] >= 13) { // not ambiguous
					++len;
					if (c1 && c1 != c2 && qual[z] >= 13) { // mismatch
						++mm;
						q += qual[z] > 33? 33 : qual[z];
					}
				}
			}
			if (j < l) break;
			x += l; y += l; len += l;
		} else if (op == BAM_CDEL) {
			for (j = 0; j < l; ++j)
				if (ref[x+j] == 0) break;
			if (j < l) break;
			x += l;
		} else if (op == BAM_CSOFT_CLIP) {
			for (j = 0; j < l; ++j) clip_q += qual[y+j];
			clip_l += l;
			y += l;
		} else if (op == BAM_CHARD_CLIP) {
			clip_q += 13 * l;
			clip_l += l;
		} else if (op == BAM_CINS) y += l;
		else if (op == BAM_CREF_SKIP) x += l;
	}
	for (i = 0, t = 1; i < mm; ++i)
		t *= (double)len / (i+1);
	t = q - 4.343 * log(t) + clip_q / 5.;
	if (t > thres) return -1;
	if (t < 0) t = 0;
	t = sqrt((thres - t) / thres) * thres;
//	fprintf(stderr, "%s %lf %d\n", bam1_qname(b), t, q);
	return (int)(t + .499);
}

// local realignment

#define MIN_REF_LEN 10
#define MIN_BAND_WIDTH 11

int bam_realn(bam1_t *b, const char *ref)
{
	int k, score, q[2], r[2], kk[2], kl[2], x, y, max, j, n_cigar, endx = -1;
	uint32_t *cigar = bam1_cigar(b);
	uint8_t *seq = bam1_seq(b);
	bam1_core_t *c = &b->core;
	q[0] = q[1] = r[0] = r[1] = kk[0] = kk[1] = kl[0] = kl[1] = -1;
	// find the right boundary
	for (k = 0, score = max = 0, x = c->pos, y = 0; k < c->n_cigar; ++k) {
		int op = cigar[k]&0xf;
		int ol = cigar[k]>>4;
		if (op == BAM_CMATCH) {
			for (j = 0; j < ol; ++j) {
				int c1, c2, z = y + j;
				c1 = bam_nt16_nt4_table[bam1_seqi(seq, z)];
				if (ref[x+j] == 0) return -1;
				c2 = bam_nt16_nt4_table[(int)bam_nt16_table[(int)ref[x+j]]];
				if (c1 < 3 && c2 < 3)
					score += c1 == c2? 5 : -4;
				if (score < 0) score = 0;
				if (score > max) max = score, q[1] = z, r[1] = x+j, kk[1] = k, kl[1] = j + 1;
			}
			x += ol; y += ol;
		} else if (op == BAM_CINS) {
			score -= 4 - ol * 3;
			y += ol;
			if (score < 0) score = 0;
		} else if (op == BAM_CDEL) {
			score -= 4 - ol * 3;
			x += ol;
			if (score < 0) score = 0;
		} else if (op == BAM_CSOFT_CLIP) y += ol;
		else if (op == BAM_CREF_SKIP) x += ol;
	}
	if (q[1] < 0) return -1; // no high scoring segments
	endx = x - 1;
	// find the left boundary
	for (k = c->n_cigar - 1, score = max = 0, x = x-1, y = y-1; k >= 0; --k) {
		int op = cigar[k]&0xf;
		int ol = cigar[k]>>4;
		if (op == BAM_CMATCH) {
			for (j = 0; j < ol; ++j) {
				int c1, c2, z = y - j;
				c1 = bam_nt16_nt4_table[bam1_seqi(seq, z)];
				if (ref[x+j] == 0) return -1;
				c2 = bam_nt16_nt4_table[(int)bam_nt16_table[(int)ref[x-j]]];
				if (c1 < 3 && c2 < 3)
					score += c1 == c2? 5 : -4;
				if (score < 0) score = 0;
				if (score > max) max = score, q[0] = z, r[0] = x-j, kk[0] = k, kl[0] = j + 1;
			}
			x -= ol; y -= ol;
		} else if (op == BAM_CINS) {
			score -= 4 - ol * 3;
			y -= ol;
			if (score < 0) score = 0;
		} else if (op == BAM_CDEL) {
			score -= 4 - ol * 3;
			x -= ol;
			if (score < 0) score = 0;
		} else if (op == BAM_CSOFT_CLIP) y -= ol;
		else if (op == BAM_CREF_SKIP) x -= ol;
	}
	if (q[0] < 0 || q[1] - q[0] < 15) return -1; // the high-scoring segment is too short
	// modify CIGAR
	n_cigar = 0;
	cigar = calloc(c->n_cigar + 4, 4);
	if (q[0] != 0) cigar[n_cigar++] = (uint32_t)q[0]<<4 | BAM_CSOFT_CLIP;
	if (r[0] != c->pos) cigar[n_cigar++] = (uint32_t)(r[0] - c->pos)<<4 | BAM_CREF_SKIP;
	if (kk[0] == kk[1]) {
		cigar[n_cigar++] = (uint32_t)(kl[0] + kl[1] - (bam1_cigar(b)[kk[0]]>>4))<<4 | BAM_CMATCH;
	} else {
		cigar[n_cigar++] = (uint32_t)kl[0]<<4 | BAM_CMATCH;
		for (k = kk[0] + 1; k < kk[1]; ++k)
			cigar[n_cigar++] = bam1_cigar(b)[k];
		cigar[n_cigar++] = (uint32_t)kl[1]<<4 | BAM_CMATCH; // FIXME: add ref_skip after this line
	}
	if (q[1] + 1 != c->l_qseq)
		cigar[n_cigar++] = (uint32_t)(c->l_qseq - q[1] - 1)<<4 | BAM_CSOFT_CLIP;
	// copy over the alignment
	if (4 * (n_cigar - (int)c->n_cigar) + b->data_len > b->m_data) { // enlarge b->data
		b->m_data = 4 * (n_cigar - (int)c->n_cigar) + b->data_len;
		kroundup32(b->m_data);
		b->data = realloc(b->data, b->m_data);
	}
	if (n_cigar != (int)c->n_cigar) { // move data
		memmove(b->data + c->l_qname + 4 * n_cigar, bam1_seq(b), b->data_len - c->l_qname - 4 * c->n_cigar);
		b->data_len += 4 * (n_cigar - (int)c->n_cigar);
	}
	memcpy(bam1_cigar(b), cigar, n_cigar * 4);
	c->n_cigar = n_cigar;
	free(cigar);
	return 0;
}

int bam_fillmd(int argc, char *argv[])
{
	int c, is_equal = 0, tid = -2, ret, len, is_bam_out, is_sam_in, is_uncompressed, max_nm = 0, is_realn, capQ = 0;
	samfile_t *fp, *fpout = 0;
	faidx_t *fai;
	char *ref = 0, mode_w[8], mode_r[8];
	bam1_t *b;

	is_bam_out = is_sam_in = is_uncompressed = is_realn = 0;
	mode_w[0] = mode_r[0] = 0;
	strcpy(mode_r, "r"); strcpy(mode_w, "w");
	while ((c = getopt(argc, argv, "reubSC:n:")) >= 0) {
		switch (c) {
		case 'r': is_realn = 1; break;
		case 'e': is_equal = 1; break;
		case 'b': is_bam_out = 1; break;
		case 'u': is_uncompressed = is_bam_out = 1; break;
		case 'S': is_sam_in = 1; break;
		case 'n': max_nm = atoi(optarg); break;
		case 'C': capQ = atoi(optarg); break;
		default: fprintf(stderr, "[bam_fillmd] unrecognized option '-%c'\n", c); return 1;
		}
	}
	if (!is_sam_in) strcat(mode_r, "b");
	if (is_bam_out) strcat(mode_w, "b");
	else strcat(mode_w, "h");
	if (is_uncompressed) strcat(mode_w, "u");
	if (optind + 1 >= argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   samtools fillmd [-eubrS] <aln.bam> <ref.fasta>\n\n");
		fprintf(stderr, "Options: -e       change identical bases to '='\n");
		fprintf(stderr, "         -u       uncompressed BAM output (for piping)\n");
		fprintf(stderr, "         -b       compressed BAM output\n");
		fprintf(stderr, "         -S       the input is SAM with header\n");
		fprintf(stderr, "         -r       read-independent local realignment\n\n");
		return 1;
	}
	fp = samopen(argv[optind], mode_r, 0);
	if (fp == 0) return 1;
	if (is_sam_in && (fp->header == 0 || fp->header->n_targets == 0)) {
		fprintf(stderr, "[bam_fillmd] input SAM does not have header. Abort!\n");
		return 1;
	}
	fpout = samopen("-", mode_w, fp->header);
	fai = fai_load(argv[optind+1]);

	b = bam_init1();
	while ((ret = samread(fp, b)) >= 0) {
		if (b->core.tid >= 0) {
			if (tid != b->core.tid) {
				free(ref);
				ref = fai_fetch(fai, fp->header->target_name[b->core.tid], &len);
				tid = b->core.tid;
				if (ref == 0)
					fprintf(stderr, "[bam_fillmd] fail to find sequence '%s' in the reference.\n",
							fp->header->target_name[tid]);
			}
			if (is_realn) bam_realn(b, ref);
			if (capQ > 10) {
				int q = bam_cap_mapQ(b, ref, capQ);
				if (b->core.qual > q) b->core.qual = q;
			}
			if (ref) bam_fillmd1_core(b, ref, is_equal, max_nm);
		}
		samwrite(fpout, b);
	}
	bam_destroy1(b);

	free(ref);
	fai_destroy(fai);
	samclose(fp); samclose(fpout);
	return 0;
}
