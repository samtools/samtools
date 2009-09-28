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
#include "bgzf.h"

static const int WINDOW_SIZE = 64 * 1024;

static int bgzip_main_usage()
{
	printf("\n");
	printf("Usage:   bgzip [options] [file] ...\n\n");
	printf("Options: -c      write on standard output, keep original files unchanged\n");
	printf("         -d      decompress\n");
	// printf("         -l      list compressed file contents\n");
	printf("         -b INT  decompress at virtual file pointer INT\n");
	printf("         -s INT  decompress INT bytes in the uncompressed file\n");
	printf("         -h      give this help\n");
	printf("\n");
	return 0;
}

static int write_open(const char *fn, int is_forced)
{
	int fd = -1;
	char c;
	if (!is_forced) {
		if ((fd = open(fn, O_WRONLY | O_CREAT | O_TRUNC | O_EXCL, 0666)) < 0 && errno == EEXIST) {
			printf("bgzip: %s already exists; do you wish to overwrite (y or n)? ", fn);
			scanf("%c", &c);
			if (c != 'Y' && c != 'y') {
				printf("bgzip: not overwritten\n");
				exit(1);
			}
		}
	}
	if (fd < 0) {
		if ((fd = open(fn, O_WRONLY | O_CREAT | O_TRUNC, 0666)) < 0) {
			fprintf(stderr, "bgzip: %s: Fail to write\n", fn);
			exit(1);
		}
	}
	return fd;
}

static
void
fail(BGZF* fp)
{
    printf("Error: %s\n", fp->error);
    exit(1);
}

int main(int argc, char **argv)
{
	int c, compress, pstdout, is_forced;
	BGZF *rz;
	void *buffer;
	long start, end, size;

	compress = 1; pstdout = 0; start = 0; size = -1; end = -1; is_forced = 0;
	while((c  = getopt(argc, argv, "cdlhfb:s:")) >= 0){
		switch(c){
		case 'h': return bgzip_main_usage();
		case 'd': compress = 0; break;
		case 'c': pstdout = 1; break;
                // case 'l': compress = 2; break;
		case 'b': start = atol(optarg); break;
		case 's': size = atol(optarg); break;
		case 'f': is_forced = 1; break;
		}
	}
	if (size >= 0) end = start + size;
	if(end >= 0 && end < start){
		fprintf(stderr, " -- Illegal region: [%ld, %ld] --\n", start, end);
		return 1;
	}
	if(compress == 1){
		int f_src, f_dst = -1;
		if(argc > optind){
			if((f_src = open(argv[optind], O_RDONLY)) < 0){
				fprintf(stderr, " -- Cannot open file: %s --\n", argv[optind]);
				return 1;
			}
			if(pstdout){
				f_dst = fileno(stdout);
			} else {
				char *name = malloc(sizeof(strlen(argv[optind]) + 5));
				strcpy(name, argv[optind]);
				strcat(name, ".gz");
				f_dst = write_open(name, is_forced);
				if (f_dst < 0) return 1;
				free(name);
			}
		} else if(pstdout){ 
			f_src = fileno(stdin);
			f_dst = fileno(stdout);
		} else return bgzip_main_usage();
		rz = bgzf_fdopen(f_dst, "w");
		buffer = malloc(WINDOW_SIZE);
		while((c = read(f_src, buffer, WINDOW_SIZE)) > 0) {
                  if (bgzf_write(rz, buffer, c) < 0) {
                    fail(rz);
                  }
                }
                // f_dst will be closed here
		if (bgzf_close(rz) < 0) {
                  fail(rz);
                }
		if (argc > optind) unlink(argv[optind]);
		free(buffer);
		close(f_src);
		return 0;
	} else {
		if(argc <= optind) return bgzip_main_usage();
                int f_dst;
                if (argc > optind && !pstdout) {
                  char *name;
                  if (strstr(argv[optind], ".gz") - argv[optind] != strlen(argv[optind]) - 3) {
                    printf("bgzip: %s: unknown suffix -- ignored\n", argv[optind]);
                    return 1;
                  }
                  name = strdup(argv[optind]);
                  name[strlen(name) - 3] = '\0';
                  f_dst = write_open(name, is_forced);
                  free(name);
                } else f_dst = fileno(stdout);
                rz = bgzf_open(argv[optind], "r");
                if (rz == NULL) {
                  printf("Could not open file: %s\n", argv[optind]);
                  return 1;
                }
                buffer = malloc(WINDOW_SIZE);
                if (bgzf_seek(rz, start, SEEK_SET) < 0) {
                  fail(rz);
                }
                while(1){
                  if(end < 0) c = bgzf_read(rz, buffer, WINDOW_SIZE);
                  else c = bgzf_read(rz, buffer, (end - start > WINDOW_SIZE)? WINDOW_SIZE:(end - start));
                  if(c == 0) break;
                  if (c < 0) fail(rz);
                  start += c;
                  write(f_dst, buffer, c);
                  if(end >= 0 && start >= end) break;
                }
                free(buffer);
		if (bgzf_close(rz) < 0) {
                  fail(rz);
                }
                if (!pstdout) unlink(argv[optind]);
		return 0;
	}
}

