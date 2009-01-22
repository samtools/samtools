#include <string.h>
#include <stdlib.h>
#include "glf.h"

#ifdef _NO_BGZF
// then alias bgzf_*() functions
#endif

static inline int bam_is_big_endian()
{
	long one= 1;
	return !(*((char *)(&one)));
}

glf_header_t *glf_header_init()
{
	return (glf_header_t*)calloc(1, sizeof(glf_header_t));
}

glf_header_t *glf_header_read(glfFile fp)
{
	glf_header_t *h;
	char magic[4];
	if (bam_is_big_endian()) {
		fprintf(stderr, "[glf_header_read] Big endian is detected. Abort.\n");
		exit(1);
	}
	h = (glf_header_t*)calloc(1, sizeof(glf_header_t));
	bgzf_read(fp, magic, 4);
	if (strncmp(magic, "GLF\2", 4)) {
		fprintf(stderr, "[glf_header_read] invalid magic. Abort.\n");
		exit(1);
	}
	bgzf_read(fp, &h->l_text, 4);
	if (h->l_text) {
		h->text = (uint8_t*)calloc(h->l_text + 1, 1);
		bgzf_read(fp, h->text, h->l_text);
	}
	return h;
}

void glf_header_write(glfFile fp, const glf_header_t *h)
{
	bgzf_write(fp, "GLF\2", 4);
	bgzf_write(fp, &h->l_text, 4);
	if (h->l_text) bgzf_write(fp, h->text, h->l_text);
}

void glf_header_destroy(glf_header_t *h)
{
	free(h->text);
	free(h);
}

char *glf_ref_read(glfFile fp)
{
	int32_t n;
	char *str;
	if (bgzf_read(fp, &n, 4) != 4) return 0;
	if (n < 0) {
		fprintf(stderr, "[glf_ref_read] invalid reference name length: %d.\n", n);
		return 0;
	}
	str = (char*)calloc(n + 1, 1); // not necesarily n+1 in fact
	bgzf_read(fp, str, n);
	return str;
}

void glf_ref_write(glfFile fp, const char *str)
{
	int32_t n = strlen(str);
	++n;
	bgzf_write(fp, &n, 4);
	bgzf_write(fp, str, n);
}

void glf_view_normal(const char *ref_name, glf2_t *g1)
{
	int j;
	printf("%s\t%d\t%c\t%d\t%d\t%d", ref_name, g1->pos + 1, "XACMGRSVTWYHKDBN"[g1->ref_base],
		   g1->depth, g1->max_mapQ, g1->min_lk);
	for (j = 0; j != 10; ++j) printf("\t%d", g1->lk[j]);
	printf("\n");
}

static char *glf_read_indel(glfFile fp, char *str, int *max, int16_t indel)
{
	int l = indel > 0? indel : -indel;
	if (l + 1 > *max) {
		*max = l + 1;
		str = (char*)realloc(str, *max);
	}
	bgzf_read(fp, str, l);
	str[l] = 0;
	return str;
}

void glf_view(glfFile fp)
{
	glf_header_t *h;
	char *name, *str;
	glf2_t g2;
	int max;
	
	h = glf_header_read(fp);
	str = 0; max = 0;
	while ((name = glf_ref_read(fp)) != 0) {
		while (bgzf_read(fp, &g2, sizeof(glf2_t))) {
			if (g2.type == GLF_TYPE_END) break;
			else if (g2.type == GLF_TYPE_NORMAL) glf_view_normal(name, &g2);
			else if (g2.type == GLF_TYPE_INDEL) {
				int16_t indel1, indel2;
				printf("%s\t%d\t*\t%d\t%d\t%d\t", name, g2.pos + 1, g2.depth, g2.max_mapQ, g2.min_lk);
				printf("%d\t%d\t%d\t", g2.lk[0], g2.lk[1], g2.lk[2]);
				indel1 = *(int16_t*)(g2.lk + 3);
				indel2 = *(int16_t*)(g2.lk + 5);
				printf("%d\t%d\t", indel1, indel2);
				if (indel1) {
					str = glf_read_indel(fp, str, &max, indel1);
					printf("%c%d%s\t", indel1>0? '+':'-', indel1>0?indel1:-indel1, str);
				} else printf("*\t");
				if (indel2) {
					str = glf_read_indel(fp, str, &max, indel2);
					printf("%c%d%s\n", indel2>0? '+':'-', indel2>0?indel2:-indel2, str);
				} else printf("*\n");
			}
		}
		free(name);
	}
	glf_header_destroy(h);
	free(str);
}

int glf_view_main(int argc, char *argv[])
{
	glfFile fp;
	if (argc == 1) {
		fprintf(stderr, "Usage: glfview <in.glf>\n");
		return 1;
	}
	fp = (strcmp(argv[1], "-") == 0)? bgzf_fdopen(fileno(stdin), "r") : bgzf_open(argv[1], "r");
	if (fp == 0) {
		fprintf(stderr, "Fail to open file '%s'\n", argv[1]);
		return 1;
	}
	glf_view(fp);
	bgzf_close(fp);
	return 0;
}

#ifdef GLFVIEW_MAIN
int main(int argc, char *argv[])
{
	return glf_view_main(argc, argv);
}
#endif
