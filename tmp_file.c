/*
    tmp_file.c - write to and read from a temporary binary file
    for fast storage plus added compression.

    Copyright (C) 2017 Genome Research Ltd.

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

#include <config.h>

#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdarg.h>
#include <fcntl.h>
#include <errno.h>

#ifdef _WIN32
#include <windows.h>
#endif /* _WIN32 */

#include "tmp_file.h"
#include "htslib/sam.h"


static void tmp_print_error(tmp_file_t *tmp, const char *fmt, ...) {
    va_list argp;

    if (tmp->verbose) {
        va_start(argp, fmt);
        vfprintf(stderr, fmt, argp);
        va_end(argp);
    }
}


static int tmp_file_init(tmp_file_t *tmp, int verbose) {
    tmp->stream       = LZ4_createStream();
    tmp->data_size    = 0;
    tmp->group_size   = TMP_SAM_GROUP_SIZE;
    tmp->input_size   = 0;
    tmp->read_size    = 0;
    tmp->output_size  = 0;
    tmp->entry_number = 0;
    tmp->offset = 0;
    tmp->max_data_size    = TMP_SAM_MAX_DATA + sizeof(bam1_t); // arbitrary but growable
    tmp->ring_buffer_size = TMP_SAM_RING_SIZE; // arbitrary (min 64K) but growable
    tmp->comp_buffer_size = LZ4_COMPRESSBOUND(tmp->max_data_size * tmp->group_size);
    tmp->data = NULL;
    tmp->ring_buffer = malloc(sizeof(uint8_t) * tmp->ring_buffer_size);
    tmp->ring_index  = tmp->ring_buffer;
    tmp->comp_buffer = malloc(tmp->comp_buffer_size);
    tmp->verbose = verbose;
    tmp->dict = NULL;
    tmp->groups_written = 0;

    if (!tmp->ring_buffer || !tmp->comp_buffer || !tmp->stream) {
        tmp_print_error(tmp, "[tmp_file] Error: unable to allocate compression buffers.\n");
        return TMP_SAM_MEM_ERROR;
    }

    return TMP_SAM_OK;
}


/*
 * Opens the temp file and initialises memory.
 * Verbose mode prints out error messages to stderr.
 * Returns 0 on success, a negative number on failure.
 */
int tmp_file_open_write(tmp_file_t *tmp, char *tmp_name, int verbose) {
    int ret;
    unsigned int count = 1;
    const unsigned int max_count = 100000; // more tries than this then something else is wrong
    int fd;

    if ((ret = tmp_file_init(tmp, verbose))) {
        return ret;
    }

    // make space to write extended file name
    if ((tmp->name = malloc(strlen(tmp_name) + 7)) == NULL) {
        tmp_print_error(tmp, "[tmp_file] Error: unable to allocate memory for %s.\n", tmp_name);
        return TMP_SAM_MEM_ERROR;
    }

    // make sure temp file has a unique name
    while (count < max_count) {
        sprintf(tmp->name, "%s.%d", tmp_name, count);


        #ifdef _WIN32
        if ((fd = _open(tmp->name, O_RDWR|O_CREAT|O_EXCL|O_BINARY|O_TEMPORARY, 0600)) == -1) {
        #else
        if ((fd = open(tmp->name, O_RDWR|O_CREAT|O_EXCL, 0600)) == -1) {
        #endif /* _WIN32 */

            if (errno != EEXIST) {
                tmp_print_error(tmp, "[tmp_file] Error: unable to create tmp file %s.\n", tmp->name);
                return TMP_SAM_FILE_ERROR;
            }

            count++;
            continue;
        }

        break;
    }

    if (count >= max_count) {
        tmp_print_error(tmp, "[tmp_file] Error: unable to create unique temp file.\n");
        return TMP_SAM_FILE_ERROR;
    }

    if ((tmp->fp = fdopen(fd, "w+b")) == NULL) {
        tmp_print_error(tmp, "[tmp_file] Error: unable to open write file %s.\n", tmp->name);
        return TMP_SAM_FILE_ERROR;
    }

    #ifndef _WIN32
    unlink(tmp->name); // should auto delete when closed on linux
    #endif

    return TMP_SAM_OK;
}


/*
 * The ring buffer stores precompressionn/post decompression data.  LZ4 requires that
 * previous data (64K worth) be available for efficient compression.  This function grows
 * the ring buffer when needed.
 * Returns 0 on success, a negative number on failure.
 */
static int tmp_file_grow_ring_buffer(tmp_file_t *tmp, size_t new_size) {
    // save the dictionary so lz4 can continue to function
    int dict_size = 64 * 1024; // 64K max size

    if (tmp->groups_written) {
        // if compression has been done then there is a dictionary to save

        if (tmp->dict == NULL) {

            if ((tmp->dict = malloc(sizeof(char) * dict_size)) == NULL) {
                tmp_print_error(tmp, "[tmp_file] Error: unable to allocate memory for compression dictionary.\n");
                return TMP_SAM_MEM_ERROR;
            }
        }

        if (LZ4_saveDict(tmp->stream, tmp->dict, dict_size) == 0) {
            tmp_print_error(tmp, "[tmp_file] Error: unable to save compression dictionary.\n");
            return TMP_SAM_LZ4_ERROR;
        }
    }

    if ((tmp->ring_buffer = realloc(tmp->ring_buffer, sizeof(char) * new_size)) == NULL) {
        tmp_print_error(tmp, "[tmp_file] Error: unable to reallocate ring buffer.\n");
        return TMP_SAM_MEM_ERROR;
    }

    tmp->ring_buffer_size = new_size;

    return TMP_SAM_OK;
}


/*
 * This does the actual compression and writing to disk.  On disk format consists of a
 * single size_t for the size of the compressed data followed by the data itself.
 * Returns 0 on success, a negative number on failure.
 */
static int tmp_file_write_to_file(tmp_file_t *tmp) {
    size_t comp_size;

    if (tmp->input_size > tmp->max_data_size) {
        tmp->max_data_size += tmp->input_size + sizeof(bam1_t);
        tmp->comp_buffer_size = LZ4_COMPRESSBOUND(tmp->max_data_size);

        if ((tmp->comp_buffer = realloc(tmp->comp_buffer, sizeof(char) * tmp->comp_buffer_size)) == NULL) {
            tmp_print_error(tmp, "[tmp_file] Error: unable to reallocate compression buffer.\n");
            return TMP_SAM_MEM_ERROR;
        }

        // make sure the ring buffer is big enough to accommodate the new max_data_size
        if (tmp->ring_buffer_size < tmp->max_data_size * 5) {
            int ret;
            if ((ret = tmp_file_grow_ring_buffer(tmp, tmp->max_data_size * 5))) {
                return ret;
            }
        }
    }

    tmp->ring_index = tmp->ring_buffer + tmp->offset;

    comp_size = LZ4_compress_fast_continue(tmp->stream, (const char *)tmp->ring_index,
                   tmp->comp_buffer, tmp->input_size, tmp->comp_buffer_size, 1);

    if (comp_size == 0) {
        tmp_print_error(tmp, "[tmp_file] Error: compression failed.\n");
        return TMP_SAM_LZ4_ERROR;
    }

    if (fwrite(&comp_size, sizeof(size_t), 1, tmp->fp) < 1) {
        tmp_print_error(tmp, "[tmp_file] Error: tmp file write size failed.\n");
        return TMP_SAM_FILE_ERROR;
    }

    if (fwrite(tmp->comp_buffer, sizeof(char), comp_size, tmp->fp) < comp_size) {
        tmp_print_error(tmp, "[tmp_file] Error: tmp file write data failed.\n");
        return TMP_SAM_FILE_ERROR;
    }

    tmp->offset += tmp->input_size;

    if (tmp->offset >= tmp->ring_buffer_size - tmp->max_data_size)
        tmp->offset = 0;

    tmp->input_size = 0;
    tmp->entry_number = 0;
    tmp->groups_written++;

    return TMP_SAM_OK;
}


/*
 * Stores an in memory bam structure for writing and if enough are gathered together writes
 * it to disk.  Mulitiple alignments compress better that single ones though after a certain number
 * there is a law of diminishing returns.
 * Returns 0 on success, a negative number on failure.
 */
int tmp_file_write(tmp_file_t *tmp, bam1_t *inbam) {

    if ((tmp->input_size + sizeof(bam1_t) + inbam->l_data) >= tmp->ring_buffer_size) {
        int ret;

        if ((ret = tmp_file_grow_ring_buffer(tmp, (tmp->input_size + sizeof(bam1_t) + inbam->l_data) * 5))) {
            tmp_print_error(tmp, "[tmp_file] Error: input line too big. (%ld).\n",
                (tmp->input_size + inbam->l_data));

            return ret;
        }
    }

    tmp->ring_index = tmp->ring_buffer + tmp->offset + tmp->input_size;

    // copy data into the ring buffer
    memcpy(tmp->ring_index, inbam, sizeof(bam1_t));
    memcpy(tmp->ring_index + sizeof(bam1_t) , inbam->data, inbam->l_data);
    tmp->input_size += sizeof(bam1_t) + inbam->l_data;
    tmp->entry_number++;

    if (tmp->entry_number == tmp->group_size) {
        // actually write out the data
        int ret;

        if ((ret = tmp_file_write_to_file(tmp))) {
            return ret;
        }
    }

    return TMP_SAM_OK;
}


/*
 * Closes the file after writing out any remaining alignments.  Adds a size_t 0 to
 * mark the end of the file.  Companion function to tmp_file_open_read below.
 * Returns 0 on success, a negative number on failure.
 */
int tmp_file_close_write(tmp_file_t *tmp) {
    size_t terminator = 0;

    if (tmp->entry_number) {
        int ret;

        if ((ret = tmp_file_write_to_file(tmp))) {
            return ret;
        }
    }

    if (fwrite(&terminator, sizeof(size_t), 1, tmp->fp) < 1) {
        tmp_print_error(tmp, "[tmp_file] Error: tmp file write terminator failed.\n");
        return TMP_SAM_FILE_ERROR;
    }

    if (fclose(tmp->fp)) {
        tmp_print_error(tmp, "[tmp_file] Error: closing tmp file %s failed.\n", tmp->name);
        return TMP_SAM_FILE_ERROR;
    }

    LZ4_freeStream(tmp->stream);

    return TMP_SAM_OK;
}


/*
 * Opens the file for reading.  Optionally, if given a pointer to an existing
 * bam1_t structure, it will free the data entry to prevent memory leaks.
 * Companion function to tmp_file_close_write above.
 * Returns 0 on success, a negative number on failure.
 */
int tmp_file_open_read(tmp_file_t *tmp, bam1_t *inbam) {

    if ((tmp->fp = fopen(tmp->name, "rb")) == NULL) {
        tmp_print_error(tmp, "[tmp_file] Error: unable to open read file %s.\n", tmp->name);
        return TMP_SAM_FILE_ERROR;
    }

    tmp->dstream = LZ4_createStreamDecode();
    tmp->offset  = 0;

    if (inbam) {
        free(inbam->data);
    }

    if (!tmp->dstream) {
        tmp_print_error(tmp, "[tmp_file] Error: unable to allocate compression stream.\n");
        return TMP_SAM_MEM_ERROR;
    }


    return TMP_SAM_OK;
}


/*
 * An alternative to tmp_file_close_write that does the same job without actually
 * closing the file. Companion function to tmp_file_begin_read below.
 * Returns 0 on success, a negative number on failure.
 */
int tmp_file_end_write(tmp_file_t *tmp) {
    size_t terminator = 0;

    if (tmp->entry_number) {
        int ret;

        if ((ret = tmp_file_write_to_file(tmp))) {
            return ret;
        }
    }

    if (fwrite(&terminator, sizeof(size_t), 1, tmp->fp) < 1) {
        tmp_print_error(tmp, "[tmp_file] Error: tmp file write terminator failed.\n");
        return TMP_SAM_FILE_ERROR;
    }

    fflush(tmp->fp);

    LZ4_freeStream(tmp->stream);

    return TMP_SAM_OK;
}


/*
 * An alternative to tmp_file_open_read but works on an open file.
 * Companion function to tmp_file_end_write above.
 * Returns 0 on success, a negative number on failure.
 */
int tmp_file_begin_read(tmp_file_t *tmp, bam1_t *inbam) {

    rewind(tmp->fp);

    tmp->dstream = LZ4_createStreamDecode();
    tmp->offset  = 0;
    tmp->entry_number = tmp->group_size;

    if (inbam) {
        free(inbam->data);
    }

    if (!tmp->dstream) {
        tmp_print_error(tmp, "[tmp_file] Error: unable to allocate compression stream.\n");
        return TMP_SAM_MEM_ERROR;
    }

    return TMP_SAM_OK;
}


/*
 * Read the next alignment, either from memory or from disk.
 * Returns size of entry on success, 0 on end of file or a negative on error.
 */
int tmp_file_read(tmp_file_t *tmp, bam1_t *inbam) {
    int entry_size;

    if (tmp->entry_number == tmp->group_size) {
        // read more data
        size_t comp_size;

        if (fread(&comp_size, sizeof(size_t), 1, tmp->fp) == 0 || comp_size == 0) {
            return TMP_SAM_OK;
        }

        if  (tmp->offset >= tmp->ring_buffer_size - tmp->max_data_size)
            tmp->offset = 0;

        tmp->ring_index = tmp->ring_buffer + tmp->offset;

        if (fread(tmp->comp_buffer, sizeof(char), comp_size, tmp->fp) > comp_size) {
            tmp_print_error(tmp, "[tmp_file] Error: error reading compressed data.\n");
            return TMP_SAM_FILE_ERROR;
        }

        tmp->output_size = LZ4_decompress_safe_continue(tmp->dstream, tmp->comp_buffer,
                        (char *)tmp->ring_index, comp_size, tmp->max_data_size);

        if (tmp->output_size == 0) {
            tmp_print_error(tmp, "[tmp_file] Error: decompression failed.\n");
            return TMP_SAM_LZ4_ERROR;
        }

        tmp->entry_number = 0;
        tmp->read_size    = 0;
    }

    tmp->ring_index = tmp->ring_buffer + tmp->offset;
    memcpy(inbam, tmp->ring_index, sizeof(bam1_t));

    if ((unsigned int)inbam->l_data > tmp->data_size) {
        if ((tmp->data = realloc(tmp->data, sizeof(uint8_t) * inbam->l_data)) == NULL) {
            tmp_print_error(tmp, "[tmp_file] Error: unable to allocate tmp data memory.\n");
            return TMP_SAM_MEM_ERROR;
        }

        tmp->data_size = inbam->l_data;
    }

    inbam->data = tmp->data;
    entry_size = sizeof(bam1_t);

    memcpy(inbam->data, tmp->ring_index + entry_size, inbam->l_data);
    entry_size += inbam->l_data;

    tmp->offset += entry_size;
    tmp->read_size += entry_size;
    tmp->entry_number++;

    if (tmp->read_size > tmp->output_size) {
        tmp_print_error(tmp, "[tmp_file] Error: wrong size of data returned RS:%ld OS:%ld EN:%ld GS:%ld.\n",
            tmp->read_size, tmp->output_size, tmp->entry_number, tmp->group_size);
        return TMP_SAM_LZ4_ERROR;
    }

    if (tmp->read_size == tmp->output_size && tmp->entry_number != tmp->group_size) {
        // hopefully the last entries in the read file
        tmp->entry_number = tmp->group_size;
    }

    return entry_size;
}


/*
 * Frees up memory, closes the file and optionally deletes it.  Giving this function
 * pointer to the bam1_t structure used for reading will set its data value to null,
 * preventing bam_destroy1() from trying to free already freed memory.
 * Returns 0 on success, a negative number or EOF on failure.
 */
int tmp_file_destroy(tmp_file_t *tmp, bam1_t *inbam, int delete) {
    int ret = 0;

    ret = fclose(tmp->fp);

    if (delete && ret == 0) {
        if (unlink(tmp->name)) {
            tmp_print_error(tmp, "[tmp_file] Error: unable to delete file %s.\n", tmp->name);
            ret = TMP_SAM_FILE_ERROR;
        }
    }

    LZ4_freeStreamDecode(tmp->dstream);
    free(tmp->ring_buffer);
    free(tmp->comp_buffer);
    free(tmp->name);
    free(tmp->data);
    free(tmp->dict);


    if (inbam) {
        inbam->data = NULL;
    }

    return ret;
}
