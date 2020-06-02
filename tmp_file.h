/*
    tmp_file.h - write to and read from a temporary binary file
    for fast storage plus added compression.

    Copyright (C) 2017, 2018 Genome Research Ltd.

    Author: Andrew Whitwham <aw7@sanger.ac.uk>

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
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE
*/

#ifndef _TMP_SAM_FILE_H_
#define _TMP_SAM_FILE_H_

#include <lz4.h>
#include "htslib/sam.h"

#ifdef __cplusplus
extern "C" {
#endif

// Group size that seems to give reasonable compression.
#define TMP_SAM_GROUP_SIZE 100

// Arbitrary initial size values but growable.
#define TMP_SAM_MAX_DATA 1024
#define TMP_SAM_RING_SIZE 1048576

// Error numbers.
#define TMP_SAM_OK 0
#define TMP_SAM_MEM_ERROR -1
#define TMP_SAM_FILE_ERROR -2
#define TMP_SAM_LZ4_ERROR -3
#define TMP_SAM_INPUT_ERROR -4

typedef struct {
    FILE *fp;
    LZ4_stream_t *stream;
    LZ4_streamDecode_t *dstream;
    size_t data_size;
    size_t max_data_size;
    size_t ring_buffer_size;
    size_t comp_buffer_size;
    size_t offset;
    uint8_t *ring_buffer;
    uint8_t *ring_index;
    char *comp_buffer;
    char *name;
    size_t group_size;
    size_t input_size;
    size_t read_size;
    size_t output_size;
    size_t entry_number;
    int verbose;
    char *dict;
    size_t groups_written;
} tmp_file_t;


/*
 * Opens the temp file and initialises memory.
 * Verbose mode prints out error messages to stderr.
 * Returns 0 on success, a negative number on failure.
 */
int tmp_file_open_write(tmp_file_t *tmp, char *tmp_name, int verbose);


/*
 * Stores an in memory bam structure for writing and if enough are gathered together writes
 * it to a file.  Multiple alignments compress better that single ones though after a certain number
 * there is a law of diminishing returns.
 * Returns 0 on success, a negative number on failure.
 */
int tmp_file_write(tmp_file_t *tmp, bam1_t *inbam);


/*
 * Marks the end of file writing.  Adds a size_t 0 to mark the end of
 * the file. Companion function to tmp_file_begin_read below.
 * Returns 0 on success, a negative number on failure.
 */
int tmp_file_end_write(tmp_file_t *tmp);

/*
 * Prepares the file for reading.
 * Companion function to tmp_file_end_write above.
 * Returns 0 on success, a negative number on failure.
 */
int tmp_file_begin_read(tmp_file_t *tmp);

/*
 * Read the next alignment, either from memory or from a file.
 * Returns size of entry on success, 0 on end of file or a negative on error.
 */
int tmp_file_read(tmp_file_t *tmp, bam1_t *inbam);


/*
 * Frees up memory, closes the file and deletes it.
 * Returns 0 on success or EOF on failure.
 */
int tmp_file_destroy(tmp_file_t *tmp);

#ifdef __cplusplus
}
#endif

#endif /* _TMP_SAM_FILE_H_ */
