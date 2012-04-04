#include <stdio.h>
#include <stdlib.h>
#include "knetfile.h"
#include "bgzf.h"
#include "bam.h"

#define BUF_SIZE 0x10000

int bam_reheader(BGZF *in, const bam_header_t *h, int fd)
{
	BGZF *fp;
	bam_header_t *old;
	int len;
	uint8_t *buf;
	if (in->is_write) return -1;
	buf = malloc(BUF_SIZE);
	old = bam_header_read(in);
	fp = bgzf_fdopen(fd, "w");
	bam_header_write(fp, h);
	if (in->block_offset < in->block_length) {
		bgzf_write(fp, in->uncompressed_block + in->block_offset, in->block_length - in->block_offset);
		bgzf_flush(fp);
	}
#ifdef _USE_KNETFILE
	while ((len = knet_read(in->fp, buf, BUF_SIZE)) > 0)
		fwrite(buf, 1, len, fp->fp);
#else
	while (!feof(in->file) && (len = fread(buf, 1, BUF_SIZE, in->file)) > 0)
		fwrite(buf, 1, len, fp->file);
#endif
	free(buf);
	fp->block_offset = in->block_offset = 0;
	bgzf_close(fp);
	return 0;
}

int main_reheader(int argc, char *argv[])
{
	bam_header_t *h;
	BGZF *in;
	if (argc != 3) {
		fprintf(stderr, "Usage: samtools reheader <in.header.sam> <in.bam>\n");
		return 1;
	}
	{ // read the header
		tamFile fph = sam_open(argv[1]);
		if (fph == 0) {
			fprintf(stderr, "[%s] fail to read the header from %s.\n", __func__, argv[1]);
			return 1;
		}
		h = sam_header_read(fph);
		sam_close(fph);
	}
	in = strcmp(argv[2], "-")? bam_open(argv[2], "r") : bam_dopen(fileno(stdin), "r");
	if (in == 0) {
		fprintf(stderr, "[%s] fail to open file %s.\n", __func__, argv[2]);
		return 1;
	}
	bam_reheader(in, h, fileno(stdout));
	bgzf_close(in);
	return 0;
}
