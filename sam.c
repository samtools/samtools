#include <string.h>
#include <unistd.h>
#include "htslib/faidx.h"
#include "sam.h"

int samthreads(samfile_t *fp, int n_threads, int n_sub_blks)
{
	if (!fp->file->is_bin || !fp->file->is_write) return -1;
	bgzf_mt(fp->x.bam, n_threads, n_sub_blks);
	return 0;
}

samfile_t *samopen(const char *fn, const char *mode, const void *aux)
{
	// hts_open() is really sam_open(), except for #define games
	samFile *hts_fp = hts_open(fn, mode, aux);
	if (hts_fp == NULL)  return NULL;

	samfile_t *fp = malloc(sizeof (samfile_t));
	fp->file = hts_fp;
	fp->x.bam = hts_fp->fp;
	if (strchr(mode, 'r')) {
		fp->header = sam_hdr_read(fp->file);  // samclose() will free this
		if (fp->header->n_targets == 0) { // no @SQ fields
			if (aux) { // check if aux is present
				bam_header_t *textheader = fp->header;
				fp->header = sam_header_read2((const char*)aux);
				// FIXME should merge in any non-@SQ headers from textheader
				bam_header_destroy(textheader);
			}
			if (fp->header->n_targets == 0 && bam_verbose >= 1)
				fprintf(stderr, "[samopen] no @SQ lines in the header.\n");
		}
	}
	else {
		fp->header = (bam_hdr_t *)aux;  // For writing, we won't free it
		if (fp->file->is_bin || strchr(mode, 'h')) sam_hdr_write(fp->file, fp->header);
	}

	return fp;
}

void samclose(samfile_t *fp)
{
	if (fp) {
		if (!fp->file->is_write && fp->header) bam_hdr_destroy(fp->header);
		sam_close(fp->file);
		free(fp);
	}
}

#if 0
int sampileup(samfile_t *fp, int mask, bam_pileup_f func, void *func_data)
{
	bam_plbuf_t *buf;
	int ret;
	bam1_t *b;
	b = bam_init1();
	buf = bam_plbuf_init(func, func_data);
	bam_plbuf_set_mask(buf, mask);
	while ((ret = samread(fp, b)) >= 0)
		bam_plbuf_push(b, buf);
	bam_plbuf_push(0, buf);
	bam_plbuf_destroy(buf);
	bam_destroy1(b);
	return 0;
}
#endif

char *samfaipath(const char *fn_ref)
{
	char *fn_list = 0;
	if (fn_ref == 0) return 0;
	fn_list = calloc(strlen(fn_ref) + 5, 1);
	strcat(strcpy(fn_list, fn_ref), ".fai");
	if (access(fn_list, R_OK) == -1) { // fn_list is unreadable
		if (access(fn_ref, R_OK) == -1) {
			fprintf(stderr, "[samfaipath] fail to read file %s.\n", fn_ref);
		} else {
			if (bam_verbose >= 3) fprintf(stderr, "[samfaipath] build FASTA index...\n");
			if (fai_build(fn_ref) == -1) {
				fprintf(stderr, "[samfaipath] fail to build FASTA index.\n");
				free(fn_list); fn_list = 0;
			}
		}
	}
	return fn_list;
}
