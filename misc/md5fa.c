#include <stdio.h>
#include <zlib.h>
#include "md5.h"
#include "kseq.h"

#define HEX_STR "0123456789abcdef"

KSEQ_INIT(gzFile, gzread)

static void md5_one(const char *fn)
{
	MD5_CTX md5_one, md5_all;
	int l, i, k;
	gzFile fp;
	kseq_t *seq;
	unsigned char unordered[16], digest[16];

	for (l = 0; l < 16; ++l) unordered[l] = 0;
	fp = strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
	if (fp == 0) {
		fprintf(stderr, "md5fa: %s: No such file or directory\n", fn);
		exit(1);
	}
	
	MD5Init(&md5_all);
	seq = kseq_init(fp);
	while ((l = kseq_read(seq)) >= 0) {
		for (i = k = 0; i < seq->seq.l; ++i) {
			if (islower(seq->seq.s[i])) seq->seq.s[k++] = toupper(seq->seq.s[i]);
			else if (isupper(seq->seq.s[i])) seq->seq.s[k++] = seq->seq.s[i];
		}
		MD5Init(&md5_one);
		MD5Update(&md5_one, (unsigned char*)seq->seq.s, k);
		MD5Final(digest, &md5_one);
		for (l = 0; l < 16; ++l) {
			printf("%c%c", HEX_STR[digest[l]>>4&0xf], HEX_STR[digest[l]&0xf]);
			unordered[l] ^= digest[l];
		}
		printf("  %s  %s\n", fn, seq->name.s);
		MD5Update(&md5_all, (unsigned char*)seq->seq.s, k);
	}
	MD5Final(digest, &md5_all);
	kseq_destroy(seq);
	for (l = 0; l < 16; ++l)
		printf("%c%c", HEX_STR[digest[l]>>4&0xf], HEX_STR[digest[l]&0xf]);
	printf("  %s  >ordered\n", fn);
	for (l = 0; l < 16; ++l)
		printf("%c%c", HEX_STR[unordered[l]>>4&0xf], HEX_STR[unordered[l]&0xf]);
	printf("  %s  >unordered\n", fn);
}

int main(int argc, char *argv[])
{
	int i;
	if (argc == 1) md5_one("-");
	else for (i = 1; i < argc; ++i) md5_one(argv[i]);
	return 0;
}
