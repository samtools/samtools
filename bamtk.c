#include <stdio.h>
#include <unistd.h>
#include <assert.h>
#include <fcntl.h>
#include "bam.h"

#ifdef _USE_KNETFILE
#include "knetfile.h"
#endif

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "0.1.7 (r510)"
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

int main_samview(int argc, char *argv[]);
int main_import(int argc, char *argv[]);

int faidx_main(int argc, char *argv[]);
int glf3_view_main(int argc, char *argv[]);

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
	fprintf(stderr, "Command: view        SAM<->BAM conversion\n");
	fprintf(stderr, "         sort        sort alignment file\n");
	fprintf(stderr, "         pileup      generate pileup output\n");
	fprintf(stderr, "         faidx       index/extract FASTA\n");
#if _CURSES_LIB != 0
	fprintf(stderr, "         tview       text alignment viewer\n");
#endif
	fprintf(stderr, "         index       index alignment\n");
	fprintf(stderr, "         fixmate     fix mate information\n");
	fprintf(stderr, "         glfview     print GLFv3 file\n");
	fprintf(stderr, "         flagstat    simple stats\n");
	fprintf(stderr, "         calmd       recalculate MD/NM tags and '=' bases\n");
	fprintf(stderr, "         merge       merge sorted alignments\n");
	fprintf(stderr, "         rmdup       remove PCR duplicates\n");
	fprintf(stderr, "\n");
	return 1;
}

int main(int argc, char *argv[])
{
#ifdef _WIN32
	setmode(fileno(stdout), O_BINARY);
	setmode(fileno(stdin),  O_BINARY);
#ifdef _USE_KNETFILE
	knet_win32_init();
#endif
#endif
	if (argc < 2) return usage();
	if (strcmp(argv[1], "view") == 0) return main_samview(argc-1, argv+1);
	else if (strcmp(argv[1], "import") == 0) return main_import(argc-1, argv+1);
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
	else if (strcmp(argv[1], "calmd") == 0) return bam_fillmd(argc-1, argv+1);
	else if (strcmp(argv[1], "fillmd") == 0) return bam_fillmd(argc-1, argv+1);
#if _CURSES_LIB != 0
	else if (strcmp(argv[1], "tview") == 0) return bam_tview_main(argc-1, argv+1);
#endif
	else {
		fprintf(stderr, "[main] unrecognized command '%s'\n", argv[1]);
		return 1;
	}
	return 0;	
}
