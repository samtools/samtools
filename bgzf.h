/*
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2008 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever.
 * Neither the Broad Institute nor MIT can be responsible for its use, misuse,
 * or functionality.
 */

#ifndef __BGZF_H
#define __BGZF_H

#include <stdint.h>
#include <stdio.h>
#include "zlib.h"
#include <stdbool.h>
//#include "zutil.h"

//typedef int8_t bool;

typedef struct {
    int file_descriptor;
    char open_mode;  // 'r' or 'w'
    bool owned_file;
    FILE* file;
    int uncompressed_block_size;
    int compressed_block_size;
    void* uncompressed_block;
    void* compressed_block;
    int64_t block_address;
    int block_length;
    int block_offset;
    const char* error;
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

#ifdef __cplusplus
}
#endif

#endif
