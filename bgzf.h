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

#ifndef __BGZF_H
#define __BGZF_H

#include <stdint.h>
#include <stdio.h>
#include <stdbool.h>
#include <zlib.h>
#ifdef _USE_KNETFILE
#include "knetfile.h"
#endif

//typedef int8_t bool;

typedef struct {
    int file_descriptor;
    char open_mode;  // 'r' or 'w'
    bool owned_file, is_uncompressed;
#ifdef _USE_KNETFILE
	union {
		knetFile *fpr;
		FILE *fpw;
	} x;
#else
    FILE* file;
#endif
    int uncompressed_block_size;
    int compressed_block_size;
    void* uncompressed_block;
    void* compressed_block;
    int64_t block_address;
    int block_length;
    int block_offset;
	int cache_size;
    const char* error;
	void *cache; // a pointer to a hash table
} BGZF;

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Open an existing file descriptor for reading or writing.
 * Mode must be either "r" or "w".
 * A subsequent bgzf_close will not close the file descriptor.
 * Returns null on error.
 */
BGZF* bgzf_fdopen(int fd, const char* __restrict mode);

/*
 * Open the specified file for reading or writing.
 * Mode must be either "r" or "w".
 * Returns null on error.
 */
BGZF* bgzf_open(const char* path, const char* __restrict mode);

/*
 * Close the BGZ file and free all associated resources.
 * Does not close the underlying file descriptor if created with bgzf_fdopen.
 * Returns zero on success, -1 on error.
 */
int bgzf_close(BGZF* fp);

/*
 * Read up to length bytes from the file storing into data.
 * Returns the number of bytes actually read.
 * Returns zero on end of file.
 * Returns -1 on error.
 */
int bgzf_read(BGZF* fp, void* data, int length);

/*
 * Write length bytes from data to the file.
 * Returns the number of bytes written.
 * Returns -1 on error.
 */
int bgzf_write(BGZF* fp, const void* data, int length);

/*
 * Return a virtual file pointer to the current location in the file.
 * No interpetation of the value should be made, other than a subsequent
 * call to bgzf_seek can be used to position the file at the same point.
 * Return value is non-negative on success.
 * Returns -1 on error.
 */
int64_t bgzf_tell(BGZF* fp);

/*
 * Set the file to read from the location specified by pos, which must
 * be a value previously returned by bgzf_tell for this file (but not
 * necessarily one returned by this file handle).
 * The where argument must be SEEK_SET.
 * Seeking on a file opened for write is not supported.
 * Returns zero on success, -1 on error.
 */
int64_t bgzf_seek(BGZF* fp, int64_t pos, int where);

/*
 * Set the cache size. Zero to disable. By default, caching is
 * disabled. The recommended cache size for frequent random access is
 * about 8M bytes.
 */
void bgzf_set_cache_size(BGZF *fp, int cache_size);

int bgzf_check_EOF(BGZF *fp);

#ifdef __cplusplus
}
#endif

#endif
