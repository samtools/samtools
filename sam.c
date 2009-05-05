#include <string.h>
#include "sam.h"

#define TYPE_BAM  1
#define TYPE_READ 2

bam_header_t *bam_header_dup(const bam_header_t *h0)
{
	bam_header_t *h;
	int i;
	h = bam_header_init();
	*h = *h0;
	h->hash = 0;
	h->text = (char*)calloc(h->l_text + 1, 1);
	memcpy(h->text, h0->text, h->l_text);
	h->target_len = (uint32_t*)calloc(h->n_targets, 4);
	h->target_name = (char**)calloc(h->n_targets, sizeof(void*));
	for (i = 0; i < h->n_targets; ++i) {
		h->target_len[i] = h0->target_len[i];
		h->target_name[i] = strdup(h0->target_name[i]);
	}
	return h;
}

bam_header_t *bam_header_parse(const char *text)
{
	bam_header_t *h;
	int i;
	char *s, *p, *q, *r;

	i = strlen(text);
	if (i < 3) return 0; // headerless
	h = bam_header_init();
	h->l_text = i;
	h->text = strdup(text);
	s = h->text;
	while ((s = strstr(s, "@SQ")) != 0) {
		++h->n_targets;
		s += 3;
	}
	h->target_len = (uint32_t*)calloc(h->n_targets, 4);
	h->target_name = (char**)calloc(h->n_targets, sizeof(void*));
	i = 0;
	s = h->text;
	while ((s = strstr(s, "@SQ")) != 0) {
		s += 3;
		r = s;
		if ((p = strstr(s, "SN:")) != 0) {
			q = p + 3;
			for (p = q; *p && *p != '\t'; ++p);
			h->target_name[i] = (char*)calloc(p - q + 1, 1);
			strncpy(h->target_name[i], q, p - q);
		} else goto header_err_ret;
		if (r < p) r = p;
		if ((p = strstr(s, "LN:")) != 0) h->target_len[i] = strtol(p + 3, 0, 10);
		else goto header_err_ret;
		if (r < p) r = p;
		s = r + 3;
		++i;
	}
	if (h->n_targets == 0) {
		bam_header_destroy(h);
		return 0;
	} else return h;

header_err_ret:
	fprintf(stderr, "[bam_header_parse] missing SN tag in a @SQ line.\n");
	bam_header_destroy(h);
	return 0;
}

samfile_t *samopen(const char *fn, const char *mode, const void *aux)
{
	samfile_t *fp;
	fp = (samfile_t*)calloc(1, sizeof(samfile_t));
	if (mode[0] == 'r') {
		const char *fn_list = (const char*)aux;
		fp->type |= TYPE_READ;
		if (mode[1] == 'b') { // binary
			fp->type |= TYPE_BAM;
			fp->x.bam = strcmp(fn, "-")? bam_open(fn, "r") : bam_dopen(fileno(stdin), "r");
			fp->header = bam_header_read(fp->x.bam);
		} else {
			fp->x.tamr = sam_open(fn);
			fp->header = sam_header_read2(fn_list);
		}
	} else if (mode[0] == 'w') {
		fp->header = bam_header_dup((const bam_header_t*)aux);
		if (mode[1] == 'b') { // binary
			fp->type |= TYPE_BAM;
			fp->x.bam = strcmp(fn, "-")? bam_open(fn, "w") : bam_dopen(fileno(stdout), "w");
			bam_header_write(fp->x.bam, fp->header);
		} else {
			int i;
			bam_header_t *alt = 0;
			alt = bam_header_parse(fp->header->text);
			fp->x.tamw = strcmp(fn, "-")? fopen(fn, "w") : stdout;
			if (alt) {
				if (alt->n_targets != fp->header->n_targets)
					fprintf(stderr, "[samopen] inconsistent number of target sequences.\n");
				bam_header_destroy(alt);
				fwrite(fp->header->text, 1, fp->header->l_text, fp->x.tamw);
			}
			if (strstr(mode, "h")) // print header
				for (i = 0; i < fp->header->n_targets; ++i)
					fprintf(fp->x.tamw, "@SQ\tSN:%s\tLN:%d\n", fp->header->target_name[i], fp->header->target_len[i]);
		}
	}
	return fp;
}

void samclose(samfile_t *fp)
{
	if (fp == 0) return;
	if (fp->header) bam_header_destroy(fp->header);
	if (fp->type & TYPE_BAM) bam_close(fp->x.bam);
	else if (fp->type & TYPE_READ) sam_close(fp->x.tamr);
	else fclose(fp->x.tamw);
	free(fp);
}

int samread(samfile_t *fp, bam1_t *b)
{
	if (fp == 0 || !(fp->type & TYPE_READ)) return -1; // not open for reading
	if (fp->type & TYPE_BAM) return bam_read1(fp->x.bam, b);
	else return sam_read1(fp->x.tamr, fp->header, b);
}

int samwrite(samfile_t *fp, const bam1_t *b)
{
	if (fp == 0 || (fp->type & TYPE_READ)) return -1; // not open for writing
	if (fp->type & TYPE_BAM) return bam_write1(fp->x.bam, b);
	else {
		char *s = bam_format1(fp->header, b);
		int l = strlen(s);
		fputs(s, fp->x.tamw); fputc('\n', fp->x.tamw);
		free(s);
		return l + 1;
	}
}

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
