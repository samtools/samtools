#include <stdio.h>
#include <unistd.h>
#include "bam.h"

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "0.1.3-3"
#endif

int bam_taf2baf(int argc, char *argv[]);
int bam_pileup(int argc, char *argv[]);
int bam_merge(int argc, char *argv[]);
int bam_index(int argc, char *argv[]);
int bam_sort(int argc, char *argv[]);
int bam_tview_main(int argc, char *argv[]);
int bam_mating(int argc, char *argv[]);
int bam_rmdup(int argc, char *argv[]);
int bam_flagstat(int argc, char *argv[]);
int bam_fillmd(int argc, char *argv[]);

int faidx_main(int argc, char *argv[]);
int glf3_view_main(int argc, char *argv[]);

static int view_aux(const bam1_t *b, void *data)
{
	bam_view1((bam_header_t*)data, b);
	return 0;
}
static int view_auxb(const bam1_t *b, void *data)
{
	bam_write1((bamFile)data, b);
	return 0;
}

int bam_view(int argc, char *argv[])
{
	bamFile fp, fpout = 0;
	bam_header_t *header;
	bam1_t *b;
	int ret, c, is_bam = 0, is_header = 0, is_headeronly = 0;
	while ((c = getopt(argc, argv, "bhH")) >= 0) {
		switch (c) {
		case 'b': is_bam = 1; break;
		case 'h': is_header = 1; break;
		case 'H': is_headeronly = 1; break;
		default: fprintf(stderr, "Unrecognized option: -%c\n", c); return 1;
		}
	}
	if (argc == optind) {
		fprintf(stderr, "Usage: samtools view [-bhH] <in.bam> [<region> [...]]\n");
		return 1;
	}
	fp = strcmp(argv[optind], "-")? bam_open(argv[optind], "r") : bam_dopen(fileno(stdin), "r");
	assert(fp);
	header = bam_header_read(fp);
	if (header == 0) {
		fprintf(stderr, "[bam_view] fail to read the BAM header. Abort!\n");
		return 1;
	}
	if (is_bam) {
		assert(fpout = bam_dopen(fileno(stdout), "w"));
		bam_header_write(fpout, header);
	}
	if (is_header || is_headeronly) {
		int i, c;
		c = header->text[header->l_text-1];
		header->text[header->l_text-1] = 0;
		printf("%s", header->text);
		if (c) putchar(c);
		header->text[header->l_text-1] = c;
		for (i = 0; i < header->n_targets; ++i)
			printf("@SQ\tSN:%s\tLN:%d\n", header->target_name[i], header->target_len[i]);
		if (is_headeronly) {
			bam_header_destroy(header);
			bam_close(fp);
			return 0;
		}
	}
	if (optind + 1 == argc) {
		b = (bam1_t*)calloc(1, sizeof(bam1_t));
		while ((ret = bam_read1(fp, b)) >= 0) bam_view1(header, b);
		if (ret < -1) fprintf(stderr, "[bam_view] truncated file? Continue anyway. (%d)\n", ret);
		free(b->data); free(b);
	} else {
		int i;
		bam_index_t *idx;
		idx = bam_index_load(argv[optind]);
		for (i = optind + 1; i < argc; ++i) {
			int tid, beg, end;
			bam_parse_region(header, argv[i], &tid, &beg, &end);
			if (tid < 0) {
				fprintf(stderr, "[bam_view] fail to get the reference name. Abort!\n");
				return 1;
			}
			if (is_bam) bam_fetch(fp, idx, tid, beg, end, fpout, view_auxb);
			else bam_fetch(fp, idx, tid, beg, end, header, view_aux);
		}
		bam_index_destroy(idx);
	}
	bam_header_destroy(header);
	bam_close(fp);
	if (is_bam) bam_close(fpout);
	return 0;
}

int bam_tagview(int argc, char *argv[])
{
	bamFile fp;
	bam_header_t *header;
	bam1_t *b;
	char tag[2];
	int ret;
	if (argc < 3) {
		fprintf(stderr, "Usage: samtools tagview <in.bam> <tag>\n");
		return 1;
	}
	fp = strcmp(argv[1], "-")? bam_open(argv[1], "r") : bam_dopen(fileno(stdin), "r");
	assert(fp);
	header = bam_header_read(fp);
	if (header == 0) {
		fprintf(stderr, "[bam_view] fail to read the BAM header. Abort!\n");
		return 1;
	}
	tag[0] = argv[2][0]; tag[1] = argv[2][1];
	b = (bam1_t*)calloc(1, sizeof(bam1_t));
	while ((ret = bam_read1(fp, b)) >= 0) {
		uint8_t *d = bam_aux_get(b, tag);
		if (d) {
			printf("%s\t%d\t", bam1_qname(b), b->core.flag);
			if (d[0] == 'Z' || d[0] == 'H') printf("%s\n", bam_aux2Z(d));
			else if (d[0] == 'f') printf("%f\n", bam_aux2f(d));
			else if (d[0] == 'd') printf("%lf\n", bam_aux2d(d));
			else if (d[0] == 'A') printf("%c\n", bam_aux2A(d));
			else if (d[0] == 'c' || d[0] == 's' || d[0] == 'i') printf("%d\n", bam_aux2i(d));
			else if (d[0] == 'C' || d[0] == 'S' || d[0] == 'I') printf("%u\n", bam_aux2i(d));
			else printf("\n");
		}
	}
	if (ret < -1) fprintf(stderr, "[bam_view] truncated file? Continue anyway. (%d)\n", ret);
	free(b->data); free(b);
	bam_header_destroy(header);
	bam_close(fp);
	return 0;
}

static int usage()
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Program: samtools (Tools for alignments in the SAM format)\n");
	fprintf(stderr, "Version: %s\n\n", PACKAGE_VERSION);
	fprintf(stderr, "Usage:   samtools <command> [options]\n\n");
	fprintf(stderr, "Command: import      import from the text format\n");
	fprintf(stderr, "         view        export to the text format\n");
	fprintf(stderr, "         sort        sort alignment file\n");
	fprintf(stderr, "         merge       merge multiple sorted alignment files\n");
	fprintf(stderr, "         pileup      generate pileup output\n");
	fprintf(stderr, "         faidx       index/extract FASTA\n");
#ifndef _NO_CURSES
	fprintf(stderr, "         tview       text alignment viewer\n");
#endif
	fprintf(stderr, "         index       index alignment\n");
	fprintf(stderr, "         fixmate     fix mate information\n");
	fprintf(stderr, "         rmdup       remove PCR duplicates\n");
	fprintf(stderr, "         glfview     print GLFv2 file\n");
	fprintf(stderr, "         flagstat    simple stats\n");
	fprintf(stderr, "         fillmd      fill the MD tag and change identical base to =\n");
	fprintf(stderr, "\n");
	return 1;
}

int main(int argc, char *argv[])
{
	if (argc < 2) return usage();
	if (strcmp(argv[1], "view") == 0) return bam_view(argc-1, argv+1);
	else if (strcmp(argv[1], "import") == 0) return bam_taf2baf(argc-1, argv+1);
	else if (strcmp(argv[1], "pileup") == 0) return bam_pileup(argc-1, argv+1);
	else if (strcmp(argv[1], "merge") == 0) return bam_merge(argc-1, argv+1);
	else if (strcmp(argv[1], "sort") == 0) return bam_sort(argc-1, argv+1);
	else if (strcmp(argv[1], "index") == 0) return bam_index(argc-1, argv+1);
	else if (strcmp(argv[1], "faidx") == 0) return faidx_main(argc-1, argv+1);
	else if (strcmp(argv[1], "fixmate") == 0) return bam_mating(argc-1, argv+1);
	else if (strcmp(argv[1], "rmdup") == 0) return bam_rmdup(argc-1, argv+1);
	else if (strcmp(argv[1], "glfview") == 0) return glf3_view_main(argc-1, argv+1);
	else if (strcmp(argv[1], "flagstat") == 0) return bam_flagstat(argc-1, argv+1);
	else if (strcmp(argv[1], "tagview") == 0) return bam_tagview(argc-1, argv+1);
	else if (strcmp(argv[1], "fillmd") == 0) return bam_fillmd(argc-1, argv+1);
#ifndef _NO_CURSES
	else if (strcmp(argv[1], "tview") == 0) return bam_tview_main(argc-1, argv+1);
#endif
	else {
		fprintf(stderr, "[main] unrecognized command '%s'\n", argv[1]);
		return 1;
	}
	return 0;	
}
