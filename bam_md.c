#include <unistd.h>
#include "faidx.h"
#include "bam.h"
#include "kstring.h"

void bam_fillmd1(bam1_t *b, char *ref, int is_equal)
{
	uint8_t *seq = bam1_seq(b);
	uint32_t *cigar = bam1_cigar(b);
	bam1_core_t *c = &b->core;
	int i, x, y, u = 0, has_md = 0;
	kstring_t *str;

	if (bam_aux_get(b, "MD")) has_md = 1;
	if (has_md == 1 && !is_equal) return; // no need to add MD
	str = (kstring_t*)calloc(1, sizeof(kstring_t));
	if (c->flag & BAM_FUNMAP) return;
	for (i = y = 0, x = c->pos; i < c->n_cigar; ++i) {
		int j, l = cigar[i]>>4, op = cigar[i]&0xf;
		if (op == BAM_CMATCH) {
			for (j = 0; j < l; ++j) {
				int z = y + j;
				int c1 = bam1_seqi(seq, z), c2 = bam_nt16_table[(int)ref[x+j]];
				if (ref[x+j] == 0) break; // out of boundary
				if ((c1 == c2 && c1 != 15 && c2 != 15) || c1 == 0) {
					if (is_equal) seq[z/2] &= (z&1)? 0xf0 : 0x0f;
					++u;
				} else {
					ksprintf(str, "%d", u);
					kputc(ref[x+j], str);
					u = 0;
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
			if (j < l) break;
			x += l;
		} else if (op == BAM_CINS || op == BAM_CSOFT_CLIP) {
			y += l;
		} else if (op == BAM_CREF_SKIP) {
			x += l;
		}
	}
	ksprintf(str, "%d", u);
	if (!has_md) bam_aux_append(b, "MD", 'Z', str->l + 1, (uint8_t*)str->s);
	free(str->s); free(str);
}

int bam_fillmd(int argc, char *argv[])
{
	int c, is_equal = 0, tid = -1, ret, len;
	bamFile fp, fpout = 0;
	bam_header_t *header;
	faidx_t *fai;
	char *ref = 0;
	bam1_t *b;

	while ((c = getopt(argc, argv, "e")) >= 0) {
		switch (c) {
		case 'e': is_equal = 1; break;
		default: fprintf(stderr, "[bam_fillmd] unrecognized option '-%c'\n", c); return 1;
		}
	}
	if (optind + 1 >= argc) {
		fprintf(stderr, "Usage: bam fillmd [-e] <aln.bam> <ref.fasta>\n");
		return 1;
	}
	fp = strcmp(argv[optind], "-")? bam_open(argv[optind], "r") : bam_dopen(fileno(stdin), "r");
	assert(fp);
	header = bam_header_read(fp);
	fpout = bam_dopen(fileno(stdout), "w");
	bam_header_write(fpout, header);
	fai = fai_load(argv[optind+1]);

	b = bam_init1();
	while ((ret = bam_read1(fp, b)) >= 0) {
		if (tid != b->core.tid)
			ref = fai_fetch(fai, header->target_name[b->core.tid], &len);
		bam_fillmd1(b, ref, is_equal);
		bam_write1(fpout, b);
	}
	bam_destroy1(b);

	fai_destroy(fai);
	bam_header_destroy(header);
	bam_close(fp); bam_close(fpout);
	return 0;
}
