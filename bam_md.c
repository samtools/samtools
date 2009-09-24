#include <unistd.h>
#include <assert.h>
#include <string.h>
#include <ctype.h>
#include "faidx.h"
#include "sam.h"
#include "kstring.h"

void bam_fillmd1(bam1_t *b, char *ref, int is_equal)
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

int bam_fillmd(int argc, char *argv[])
{
	int c, is_equal = 0, tid = -2, ret, len, is_bam_out, is_sam_in, is_uncompressed;
	samfile_t *fp, *fpout = 0;
	faidx_t *fai;
	char *ref = 0, mode_w[8], mode_r[8];
	bam1_t *b;

	is_bam_out = is_sam_in = is_uncompressed = 0;
	mode_w[0] = mode_r[0] = 0;
	strcpy(mode_r, "r"); strcpy(mode_w, "w");
	while ((c = getopt(argc, argv, "eubS")) >= 0) {
		switch (c) {
		case 'e': is_equal = 1; break;
		case 'b': is_bam_out = 1; break;
		case 'u': is_uncompressed = is_bam_out = 1; break;
		case 'S': is_sam_in = 1; break;
		default: fprintf(stderr, "[bam_fillmd] unrecognized option '-%c'\n", c); return 1;
		}
	}
	if (!is_sam_in) strcat(mode_r, "b");
	if (is_bam_out) strcat(mode_w, "b");
	else strcat(mode_w, "h");
	if (is_uncompressed) strcat(mode_w, "u");
	if (optind + 1 >= argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   samtools fillmd [-eubS] <aln.bam> <ref.fasta>\n\n");
		fprintf(stderr, "Options: -e       change identical bases to '='\n");
		fprintf(stderr, "         -u       uncompressed BAM output (for piping)\n");
		fprintf(stderr, "         -b       compressed BAM output\n");
		fprintf(stderr, "         -S       the input is SAM with header\n\n");
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
			if (ref) bam_fillmd1(b, ref, is_equal);
		}
		samwrite(fpout, b);
	}
	bam_destroy1(b);

	free(ref);
	fai_destroy(fai);
	samclose(fp); samclose(fpout);
	return 0;
}
