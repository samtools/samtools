/* The MIT License

   Copyright (c) 2008 Broad Institute / Massachusetts Institute of Technology

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in
   all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
   THE SOFTWARE.
*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <fcntl.h>
#include <unistd.h>
#include <errno.h>
#include <sys/select.h>
#include <sys/stat.h>
#include "bgzf.h"

static const int WINDOW_SIZE = 64 * 1024;

static int bgzip_main_usage()
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage:   bgzip [options] [file] ...\n\n");
	fprintf(stderr, "Options: -c      write on standard output, keep original files unchanged\n");
	fprintf(stderr, "         -d      decompress\n");
	fprintf(stderr, "         -f      overwrite files without asking\n");
	fprintf(stderr, "         -b INT  decompress at virtual file pointer INT\n");
	fprintf(stderr, "         -s INT  decompress INT bytes in the uncompressed file\n");
	fprintf(stderr, "         -h      give this help\n");
	fprintf(stderr, "\n");
	return 1;
}

static int write_open(const char *fn, int is_forced)
{
	int fd = -1;
	char c;
	if (!is_forced) {
		if ((fd = open(fn, O_WRONLY | O_CREAT | O_TRUNC | O_EXCL, 0666)) < 0 && errno == EEXIST) {
			fprintf(stderr, "[bgzip] %s already exists; do you wish to overwrite (y or n)? ", fn);
			scanf("%c", &c);
			if (c != 'Y' && c != 'y') {
				fprintf(stderr, "[bgzip] not overwritten\n");
				exit(1);
			}
		}
	}
	if (fd < 0) {
		if ((fd = open(fn, O_WRONLY | O_CREAT | O_TRUNC, 0666)) < 0) {
			fprintf(stderr, "[bgzip] %s: Fail to write\n", fn);
			exit(1);
		}
	}
	return fd;
}

static void fail(BGZF* fp)
{
    fprintf(stderr, "Error: %s\n", fp->error);
    exit(1);
}

int main(int argc, char **argv)
{
	int c, compress, pstdout, is_forced;
	BGZF *fp;
	void *buffer;
	long start, end, size;

	compress = 1; pstdout = 0; start = 0; size = -1; end = -1; is_forced = 0;
	while((c  = getopt(argc, argv, "cdhfb:s:")) >= 0){
		switch(c){
		case 'h': return bgzip_main_usage();
		case 'd': compress = 0; break;
		case 'c': pstdout = 1; break;
		case 'b': start = atol(optarg); break;
		case 's': size = atol(optarg); break;
		case 'f': is_forced = 1; break;
		}
	}
	if (size >= 0) end = start + size;
	if (end >= 0 && end < start) {
		fprintf(stderr, "[bgzip] Illegal region: [%ld, %ld]\n", start, end);
		return 1;
	}
	if (compress == 1) {
		struct stat sbuf;
		int f_src = fileno(stdin);
		int f_dst = fileno(stdout);

		if ( argc>optind )
		{
			if ( stat(argv[optind],&sbuf)<0 ) 
			{ 
				fprintf(stderr, "[bgzip] %s: %s\n", strerror(errno), argv[optind]);
				return 1; 
			}

			if ((f_src = open(argv[optind], O_RDONLY)) < 0) {
				fprintf(stderr, "[bgzip] %s: %s\n", strerror(errno), argv[optind]);
				return 1;
			}

			if (pstdout)
				f_dst = fileno(stdout);
			else
			{
				char *name = malloc(strlen(argv[optind]) + 5);
				strcpy(name, argv[optind]);
				strcat(name, ".gz");
				f_dst = write_open(name, is_forced);
				if (f_dst < 0) return 1;
				free(name);
			}
		}
		else if (!pstdout && isatty(fileno((FILE *)stdout)) )
			return bgzip_main_usage();

		fp = bgzf_fdopen(f_dst, "w");
		buffer = malloc(WINDOW_SIZE);
		while ((c = read(f_src, buffer, WINDOW_SIZE)) > 0)
			if (bgzf_write(fp, buffer, c) < 0) fail(fp);
		// f_dst will be closed here
		if (bgzf_close(fp) < 0) fail(fp);
		if (argc > optind && !pstdout) unlink(argv[optind]);
		free(buffer);
		close(f_src);
		return 0;
	} else {
		struct stat sbuf;
		int f_dst;

		if ( argc>optind )
		{
			if ( stat(argv[optind],&sbuf)<0 )
			{
				fprintf(stderr, "[bgzip] %s: %s\n", strerror(errno), argv[optind]);
				return 1;
			}
			char *name;
			int len = strlen(argv[optind]);
			if ( strcmp(argv[optind]+len-3,".gz") )
			{
				fprintf(stderr, "[bgzip] %s: unknown suffix -- ignored\n", argv[optind]);
				return 1;
			}
			fp = bgzf_open(argv[optind], "r");
			if (fp == NULL) {
				fprintf(stderr, "[bgzip] Could not open file: %s\n", argv[optind]);
				return 1;
			}

			if (pstdout) {
				f_dst = fileno(stdout);
			}
			else {
				name = strdup(argv[optind]);
				name[strlen(name) - 3] = '\0';
				f_dst = write_open(name, is_forced);
				free(name);
			}
		}
		else if (!pstdout && isatty(fileno((FILE *)stdin)) )
			return bgzip_main_usage();
		else
		{
			f_dst = fileno(stdout);
			fp = bgzf_fdopen(fileno(stdin), "r");
			if (fp == NULL) {
				fprintf(stderr, "[bgzip] Could not read from stdin: %s\n", strerror(errno));
				return 1;
			}
		}
		buffer = malloc(WINDOW_SIZE);
		if (bgzf_seek(fp, start, SEEK_SET) < 0) fail(fp);
		while (1) {
			if (end < 0) c = bgzf_read(fp, buffer, WINDOW_SIZE);
			else c = bgzf_read(fp, buffer, (end - start > WINDOW_SIZE)? WINDOW_SIZE:(end - start));
			if (c == 0) break;
			if (c < 0) fail(fp);
			start += c;
			write(f_dst, buffer, c);
			if (end >= 0 && start >= end) break;
		}
		free(buffer);
		if (bgzf_close(fp) < 0) fail(fp);
		if (!pstdout) unlink(argv[optind]);
		return 0;
	}
}
