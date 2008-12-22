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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "bgzf.h"

extern off_t ftello(FILE *stream);
extern int fseeko(FILE *stream, off_t offset, int whence);

typedef int8_t byte;

static const int DEFAULT_BLOCK_SIZE = 64 * 1024;
static const int MAX_BLOCK_SIZE = 64 * 1024;

static const int BLOCK_HEADER_LENGTH = 18;
static const int BLOCK_FOOTER_LENGTH = 8;

static const int GZIP_ID1 = 31;
static const int GZIP_ID2 = 139;
static const int CM_DEFLATE = 8;
static const int FLG_FEXTRA = 4;
static const int OS_UNKNOWN = 255;
static const int BGZF_ID1 = 66; // 'B'
static const int BGZF_ID2 = 67; // 'C'
static const int BGZF_LEN = 2;
static const int BGZF_XLEN = 6; // BGZF_LEN+4

static const int GZIP_WINDOW_BITS = -15; // no zlib header
static const int Z_DEFAULT_MEM_LEVEL = 8;


inline
void
packInt16(uint8_t* buffer, uint16_t value)
{
    buffer[0] = value;
    buffer[1] = value >> 8;
}

inline
int
unpackInt16(const uint8_t* buffer)
{
    return (buffer[0] | (buffer[1] << 8));
}

inline
void
packInt32(uint8_t* buffer, uint32_t value)
{
    buffer[0] = value;
    buffer[1] = value >> 8;
    buffer[2] = value >> 16;
    buffer[3] = value >> 24;
}

inline
int
min(int x, int y)
{
    return (x < y) ? x : y;
}

static
void
report_error(BGZF* fp, const char* message) {
    fp->error = message;
}

static
BGZF*
open_read(int fd)
{
    FILE* file = fdopen(fd, "r");
    BGZF* fp = malloc(sizeof(BGZF));
    fp->file_descriptor = fd;
    fp->open_mode = 'r';
    fp->owned_file = 0;
    fp->file = file;
    fp->uncompressed_block_size = MAX_BLOCK_SIZE;
    fp->uncompressed_block = malloc(MAX_BLOCK_SIZE);
    fp->compressed_block_size = MAX_BLOCK_SIZE;
    fp->compressed_block = malloc(MAX_BLOCK_SIZE);
    fp->block_address = 0;
    fp->block_offset = 0;
    fp->block_length = 0;
    fp->error = NULL;
    return fp;
}

static
BGZF*
open_write(int fd)
{
    FILE* file = fdopen(fd, "w");
    BGZF* fp = malloc(sizeof(BGZF));
    fp->file_descriptor = fd;
    fp->open_mode = 'w';
    fp->owned_file = 0;
    fp->file = file;
    fp->uncompressed_block_size = DEFAULT_BLOCK_SIZE;
    fp->uncompressed_block = NULL;
    fp->compressed_block_size = MAX_BLOCK_SIZE;
    fp->compressed_block = malloc(MAX_BLOCK_SIZE);
    fp->block_address = 0;
    fp->block_offset = 0;
    fp->block_length = 0;
    fp->error = NULL;
    return fp;
}

BGZF*
bgzf_open(const char* __restrict path, const char* __restrict mode)
{
    BGZF* fp = NULL;
    if (strcasecmp(mode, "r") == 0) {
	int oflag = O_RDONLY;
	int fd = open(path, oflag);
        fp = open_read(fd);
    } else if (strcasecmp(mode, "w") == 0) {
	int oflag = O_WRONLY | O_CREAT | O_TRUNC;
	int fd = open(path, oflag, 0644);
        fp = open_write(fd);
    }
    if (fp != NULL) {
        fp->owned_file = 1;
    }
    return fp;
}

BGZF*
bgzf_fdopen(int fd, const char * __restrict mode)
{
    if (strcasecmp(mode, "r") == 0) {
        return open_read(fd);
    } else if (strcasecmp(mode, "w") == 0) {
        return open_write(fd);
    } else {
        return NULL;
    }
}

static
int
deflate_block(BGZF* fp, int block_length)
{
    // Deflate the block in fp->uncompressed_block into fp->compressed_block.
    // Also adds an extra field that stores the compressed block length.

    byte* buffer = fp->compressed_block;
    int buffer_size = fp->compressed_block_size;

    // Init gzip header
    buffer[0] = GZIP_ID1;
    buffer[1] = GZIP_ID2;
    buffer[2] = CM_DEFLATE;
    buffer[3] = FLG_FEXTRA;
    buffer[4] = 0; // mtime
    buffer[5] = 0;
    buffer[6] = 0;
    buffer[7] = 0;
    buffer[8] = 0;
    buffer[9] = OS_UNKNOWN;
    buffer[10] = BGZF_XLEN;
    buffer[11] = 0;
    buffer[12] = BGZF_ID1;
    buffer[13] = BGZF_ID2;
    buffer[14] = BGZF_LEN;
    buffer[15] = 0;
    buffer[16] = 0; // placeholder for block length
    buffer[17] = 0;

    // loop to retry for blocks that do not compress enough
    int input_length = block_length;
    int compressed_length = 0;
    while (1) {

        z_stream zs;
        zs.zalloc = NULL;
        zs.zfree = NULL;
        zs.next_in = fp->uncompressed_block;
        zs.avail_in = input_length;
        zs.next_out = (void*)&buffer[BLOCK_HEADER_LENGTH];
        zs.avail_out = buffer_size - BLOCK_HEADER_LENGTH - BLOCK_FOOTER_LENGTH;

        int status = deflateInit2(&zs, Z_DEFAULT_COMPRESSION, Z_DEFLATED,
                                  GZIP_WINDOW_BITS, Z_DEFAULT_MEM_LEVEL, Z_DEFAULT_STRATEGY);
        if (status != Z_OK) {
            report_error(fp, "deflate init failed");
            return -1;
        }
        status = deflate(&zs, Z_FINISH);
        if (status != Z_STREAM_END) {
            deflateEnd(&zs);
            if (status == Z_OK) {
                // Not enough space in buffer.
                // Can happen in the rare case the input doesn't compress enough.
                // Reduce the amount of input until it fits.
                input_length -= 1024;
                if (input_length <= 0) {
                    // should never happen
                    report_error(fp, "input reduction failed");
                    return -1;
                }
                continue;
            }
            report_error(fp, "deflate failed");
            return -1;
        }
        status = deflateEnd(&zs);
        if (status != Z_OK) {
            report_error(fp, "deflate end failed");
            return -1;
        }
        compressed_length = zs.total_out;
        compressed_length += BLOCK_HEADER_LENGTH + BLOCK_FOOTER_LENGTH;
        if (compressed_length > MAX_BLOCK_SIZE) {
            // should never happen
            report_error(fp, "deflate overflow");
            return -1;
        }
        break;
    }

    packInt16((uint8_t*)&buffer[16], compressed_length-1);
    uint32_t crc = crc32(0L, NULL, 0L);
    crc = crc32(crc, fp->uncompressed_block, input_length);
    packInt32((uint8_t*)&buffer[compressed_length-8], crc);
    packInt32((uint8_t*)&buffer[compressed_length-4], input_length);

    int remaining = block_length - input_length;
    if (remaining > 0) {
        if (remaining > input_length) {
            // should never happen (check so we can use memcpy)
            report_error(fp, "remainder too large");
            return -1;
        }
        memcpy(fp->uncompressed_block,
               fp->uncompressed_block + input_length,
               remaining);
    }
    fp->block_offset = remaining;
    return compressed_length;
}

static
int
inflate_block(BGZF* fp, int block_length)
{
    // Inflate the block in fp->compressed_block into fp->uncompressed_block

    z_stream zs;
    zs.zalloc = NULL;
    zs.zfree = NULL;
    zs.next_in = fp->compressed_block + 18;
    zs.avail_in = block_length - 16;
    zs.next_out = fp->uncompressed_block;
    zs.avail_out = fp->uncompressed_block_size;

    int status = inflateInit2(&zs, GZIP_WINDOW_BITS);
    if (status != Z_OK) {
        report_error(fp, "inflate init failed");
        return -1;
    }
    status = inflate(&zs, Z_FINISH);
    if (status != Z_STREAM_END) {
        inflateEnd(&zs);
        report_error(fp, "inflate failed");
        return -1;
    }
    status = inflateEnd(&zs);
    if (status != Z_OK) {
        report_error(fp, "inflate failed");
        return -1;
    }
    return zs.total_out;
}

static
int
check_header(const byte* header)
{
    return (header[0] == GZIP_ID1 &&
            header[1] == (byte) GZIP_ID2 &&
            header[2] == Z_DEFLATED &&
            (header[3] & FLG_FEXTRA) != 0 &&
            unpackInt16((uint8_t*)&header[10]) == BGZF_XLEN &&
            header[12] == BGZF_ID1 &&
            header[13] == BGZF_ID2 &&
            unpackInt16((uint8_t*)&header[14]) == BGZF_LEN);
}

static
int
read_block(BGZF* fp)
{
    byte header[BLOCK_HEADER_LENGTH];
    int64_t block_address = ftello(fp->file);
    int count = fread(header, 1, sizeof(header), fp->file);
    if (count == 0) {
        fp->block_length = 0;
        return 0;
    }
    if (count != sizeof(header)) {
        report_error(fp, "read failed");
        return -1;
    }
    if (!check_header(header)) {
        report_error(fp, "invalid block header");
        return -1;
    }
    int block_length = unpackInt16((uint8_t*)&header[16]) + 1;
    byte* compressed_block = (byte*) fp->compressed_block;
    memcpy(compressed_block, header, BLOCK_HEADER_LENGTH);
    int remaining = block_length - BLOCK_HEADER_LENGTH;
    count = fread(&compressed_block[BLOCK_HEADER_LENGTH], 1, remaining, fp->file);
    if (count != remaining) {
        report_error(fp, "read failed");
        return -1;
    }
    count = inflate_block(fp, block_length);
    if (count < 0) {
        return -1;
    }
    if (fp->block_length != 0) {
        // Do not reset offset if this read follows a seek.
        fp->block_offset = 0;
    }
    fp->block_address = block_address;
    fp->block_length = count;
    return 0;
}

int
bgzf_read(BGZF* fp, void* data, int length)
{
    if (length <= 0) {
        return 0;
    }
    if (fp->open_mode != 'r') {
        report_error(fp, "file not open for reading");
        return -1;
    }

    int bytes_read = 0;
    byte* output = data;
    while (bytes_read < length) {
        int available = fp->block_length - fp->block_offset;
        if (available <= 0) {
            if (read_block(fp) != 0) {
                return -1;
            }
            available = fp->block_length - fp->block_offset;
            if (available <= 0) {
                break;
            }
        }
        int copy_length = min(length-bytes_read, available);
        byte* buffer = fp->uncompressed_block;
        memcpy(output, buffer + fp->block_offset, copy_length);
        fp->block_offset += copy_length;
        output += copy_length;
        bytes_read += copy_length;
    }
    if (fp->block_offset == fp->block_length) {
        fp->block_address = ftello(fp->file);
        fp->block_offset = 0;
        fp->block_length = 0;
    }
    return bytes_read;
}

static
int
flush_block(BGZF* fp)
{
    while (fp->block_offset > 0) {
        int block_length = deflate_block(fp, fp->block_offset);
        if (block_length < 0) {
            return -1;
        }
        int count = fwrite(fp->compressed_block, 1, block_length, fp->file);
        if (count != block_length) {
            report_error(fp, "write failed");
            return -1;
        }
        fp->block_address += block_length;
    }
    return 0;
}

int
bgzf_write(BGZF* fp, const void* data, int length)
{
    if (fp->open_mode != 'w') {
        report_error(fp, "file not open for writing");
        return -1;
    }

    if (fp->uncompressed_block == NULL) {
        fp->uncompressed_block = malloc(fp->uncompressed_block_size);
    }

    const byte* input = data;
    int block_length = fp->uncompressed_block_size;
    int bytes_written = 0;
    while (bytes_written < length) {
        int copy_length = min(block_length - fp->block_offset, length - bytes_written);
        byte* buffer = fp->uncompressed_block;
        memcpy(buffer + fp->block_offset, input, copy_length);
        fp->block_offset += copy_length;
        input += copy_length;
        bytes_written += copy_length;
        if (fp->block_offset == block_length) {
            if (flush_block(fp) != 0) {
                break;
            }
        }
    }
    return bytes_written;
}

int
bgzf_close(BGZF* fp)
{
    if (fp->open_mode == 'w') {
        if (flush_block(fp) != 0) {
            return -1;
        }
        if (fflush(fp->file) != 0) {
            report_error(fp, "flush failed");
            return -1;
        }
    }
    if (fp->owned_file) {
        if (fclose(fp->file) != 0) {
            return -1;
        }
    }
    free(fp->uncompressed_block);
    free(fp->compressed_block);
    free(fp);
    return 0;
}

int64_t
bgzf_tell(BGZF* fp)
{
    return ((fp->block_address << 16) | (fp->block_offset & 0xFFFF));
}

int64_t
bgzf_seek(BGZF* fp, int64_t pos, int where)
{
    if (fp->open_mode != 'r') {
        report_error(fp, "file not open for read");
        return -1;
    }
    if (where != SEEK_SET) {
        report_error(fp, "unimplemented seek option");
        return -1;
    }
    int block_offset = pos & 0xFFFF;
    int64_t block_address = (pos >> 16) & 0xFFFFFFFFFFFFLL;
    if (fseeko(fp->file, block_address, SEEK_SET) != 0) {
        report_error(fp, "seek failed");
        return -1;
    }
    fp->block_length = 0;  // indicates current block is not loaded
    fp->block_address = block_address;
    fp->block_offset = block_offset;
    return 0;
}

