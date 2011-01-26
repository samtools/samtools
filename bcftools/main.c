#include <string.h>
#include <stdlib.h>
#include <sys/stat.h>
#include "bcf.h"

int bcfview(int argc, char *argv[]);
int bcf_main_index(int argc, char *argv[]);

#define BUF_SIZE 0x10000

int bcf_cat(int n, char * const *fn)
{
	int i;
	bcf_t *out;
	uint8_t *buf;
	buf = malloc(BUF_SIZE);
	out = bcf_open("-", "w");
	for (i = 0; i < n; ++i) {
		bcf_t *in;
		bcf_hdr_t *h;
		off_t end;
		struct stat s;
		in = bcf_open(fn[i], "r");
		h = bcf_hdr_read(in);
		if (i == 0) bcf_hdr_write(out, h);
		bcf_hdr_destroy(h);
#ifdef _USE_KNETFILE
		fstat(knet_fileno(in->fp->x.fpr), &s);
		end = s.st_size - 28;
		while (knet_tell(in->fp->x.fpr) < end) {
			int size = knet_tell(in->fp->x.fpr) + BUF_SIZE < end? BUF_SIZE : end - knet_tell(in->fp->x.fpr);
			knet_read(in->fp->x.fpr, buf, size);
			fwrite(buf, 1, size, out->fp->x.fpw);
		}
#else
		abort(); // FIXME: not implemented
#endif
		bcf_close(in);
	}
	bcf_close(out);
	free(buf);
	return 0;
}

int main(int argc, char *argv[])
{
	if (argc == 1) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   bcftools <command> <arguments>\n\n");
		fprintf(stderr, "Command: view      print, extract, convert and call SNPs from BCF\n");
		fprintf(stderr, "         index     index BCF\n");
		fprintf(stderr, "         cat       concatenate BCFs\n");
		fprintf(stderr, "\n");
		return 1;
	}
	if (strcmp(argv[1], "view") == 0) return bcfview(argc-1, argv+1);
	else if (strcmp(argv[1], "index") == 0) return bcf_main_index(argc-1, argv+1);
	else if (strcmp(argv[1], "cat") == 0) return bcf_cat(argc-2, argv+2);
	else {
		fprintf(stderr, "[main] Unrecognized command.\n");
		return 1;
	}
	return 0;
}
