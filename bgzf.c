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

/*
  2009-06-29 by lh3: cache recent uncompressed blocks.
  2009-06-25 by lh3: optionally use my knetfile library to access file on a FTP.
  2009-06-12 by lh3: support a mode string like "wu" where 'u' for uncompressed output */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "bgzf.h"

#include "khash.h"
typedef struct {
	int size;
	uint8_t *block;
	int64_t end_offset;
} cache_t;
KHASH_MAP_INIT_INT64(cache, cache_t)

#if defined(_WIN32) || defined(_MSC_VER)
#define ftello(fp) ftell(fp)
#define fseeko(fp, offset, whence) fseek(fp, offset, whence)
#else
extern off_t ftello(FILE *stream);
extern int fseeko(FILE *stream, off_t offset, int whence);
#endif

typedef int8_t bgzf_byte_t;

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

static inline
int
bgzf_min(int x, int y)
{
    return (x < y) ? x : y;
}

static
void
report_error(BGZF* fp, const char* message) {
    fp->error = message;
}

static BGZF *bgzf_read_init()
{
	BGZF *fp;
	fp = calloc(1, sizeof(BGZF));
    fp->uncompressed_block_size = MAX_BLOCK_SIZE;
    fp->uncompressed_block = malloc(MAX_BLOCK_SIZE);
    fp->compressed_block_size = MAX_BLOCK_SIZE;
    fp->compressed_block = malloc(MAX_BLOCK_SIZE);
	fp->cache_size = 0;
	fp->cache = kh_init(cache);
	return fp;
}

static
BGZF*
open_read(int fd)
{
#ifdef _USE_KNETFILE
    knetFile *file = knet_dopen(fd, "r");
#else
    FILE* file = fdopen(fd, "r");
#endif
    BGZF* fp;
	if (file == 0) return 0;
	fp = bgzf_read_init();
    fp->file_descriptor = fd;
    fp->open_mode = 'r';
#ifdef _USE_KNETFILE
    fp->x.fpr = file;
#else
    fp->file = file;
#endif
    return fp;
}

static
BGZF*
open_write(int fd, bool is_uncompressed)
{
    FILE* file = fdopen(fd, "w");
    BGZF* fp;
	if (file == 0) return 0;
	fp = malloc(sizeof(BGZF));
    fp->file_descriptor = fd;
    fp->open_mode = 'w';
    fp->owned_file = 0; fp->is_uncompressed = is_uncompressed;
#ifdef _USE_KNETFILE
    fp->x.fpw = file;
#else
    fp->file = file;
#endif
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
    if (mode[0] == 'r' || mode[0] == 'R') { /* The reading mode is preferred. */
#ifdef _USE_KNETFILE
		knetFile *file = knet_open(path, mode);
		if (file == 0) return 0;
		fp = bgzf_read_init();
		fp->file_descriptor = -1;
		fp->open_mode = 'r';
		fp->x.fpr = file;
#else
		int fd, oflag = O_RDONLY;
#ifdef _WIN32
		oflag |= O_BINARY;
#endif
		fd = open(path, oflag);
		if (fd == -1) return 0;
        fp = open_read(fd);
#endif
    } else if (mode[0] == 'w' || mode[0] == 'W') {
		int fd, oflag = O_WRONLY | O_CREAT | O_TRUNC;
#ifdef _WIN32
		oflag |= O_BINARY;
#endif
		fd = open(path, oflag, 0666);
		if (fd == -1) return 0;
        fp = open_write(fd, strstr(mode, "u")? 1 : 0);
    }
    if (fp != NULL) fp->owned_file = 1;
    return fp;
}

BGZF*
bgzf_fdopen(int fd, const char * __restrict mode)
{
	if (fd == -1) return 0;
    if (mode[0] == 'r' || mode[0] == 'R') {
        return open_read(fd);
    } else if (mode[0] == 'w' || mode[0] == 'W') {
        return open_write(fd, strstr(mode, "u")? 1 : 0);
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

    bgzf_byte_t* buffer = fp->compressed_block;
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
		int compress_level = fp->is_uncompressed? 0 : Z_DEFAULT_COMPRESSION;
        z_stream zs;
        zs.zalloc = NULL;
        zs.zfree = NULL;
        zs.next_in = fp->uncompressed_block;
        zs.avail_in = input_length;
        zs.next_out = (void*)&buffer[BLOCK_HEADER_LENGTH];
        zs.avail_out = buffer_size - BLOCK_HEADER_LENGTH - BLOCK_FOOTER_LENGTH;

        int status = deflateInit2(&zs, compress_level, Z_DEFLATED,
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
check_header(const bgzf_byte_t* header)
{
    return (header[0] == GZIP_ID1 &&
            header[1] == (bgzf_byte_t) GZIP_ID2 &&
            header[2] == Z_DEFLATED &&
            (header[3] & FLG_FEXTRA) != 0 &&
            unpackInt16((uint8_t*)&header[10]) == BGZF_XLEN &&
            header[12] == BGZF_ID1 &&
            header[13] == BGZF_ID2 &&
            unpackInt16((uint8_t*)&header[14]) == BGZF_LEN);
}

static void free_cache(BGZF *fp)
{
	khint_t k;
	khash_t(cache) *h = (khash_t(cache)*)fp->cache;
	if (fp->open_mode != 'r') return;
	for (k = kh_begin(h); k < kh_end(h); ++k)
		if (kh_exist(h, k)) free(kh_val(h, k).block);
	kh_destroy(cache, h);
}

static int load_block_from_cache(BGZF *fp, int64_t block_address)
{
	khint_t k;
	cache_t *p;
	khash_t(cache) *h = (khash_t(cache)*)fp->cache;
	k = kh_get(cache, h, block_address);
	if (k == kh_end(h)) return 0;
	p = &kh_val(h, k);
	if (fp->block_length != 0) fp->block_offset = 0;
	fp->block_address = block_address;
	fp->block_length = p->size;
	memcpy(fp->uncompressed_block, p->block, MAX_BLOCK_SIZE);
#ifdef _USE_KNETFILE
	knet_seek(fp->x.fpr, p->end_offset, SEEK_SET);
#else
	fseeko(fp->file, p->end_offset, SEEK_SET);
#endif
	return p->size;
}

static void cache_block(BGZF *fp, int size)
{
	int ret;
	khint_t k;
	cache_t *p;
	khash_t(cache) *h = (khash_t(cache)*)fp->cache;
	if (MAX_BLOCK_SIZE >= fp->cache_size) return;
	if ((kh_size(h) + 1) * MAX_BLOCK_SIZE > fp->cache_size) {
		/* A better way would be to remove the oldest block in the
		 * cache, but here we remove a random one for simplicity. This
		 * should not have a big impact on performance. */
		for (k = kh_begin(h); k < kh_end(h); ++k)
			if (kh_exist(h, k)) break;
		if (k < kh_end(h)) {
			free(kh_val(h, k).block);
			kh_del(cache, h, k);
		}
	}
	k = kh_put(cache, h, fp->block_address, &ret);
	if (ret == 0) return; // if this happens, a bug!
	p = &kh_val(h, k);
	p->size = fp->block_length;
	p->end_offset = fp->block_address + size;
	p->block = malloc(MAX_BLOCK_SIZE);
	memcpy(kh_val(h, k).block, fp->uncompressed_block, MAX_BLOCK_SIZE);
}

int
bgzf_read_block(BGZF* fp)
{
    bgzf_byte_t header[BLOCK_HEADER_LENGTH];
	int count, size = 0;
#ifdef _USE_KNETFILE
    int64_t block_address = knet_tell(fp->x.fpr);
	if (load_block_from_cache(fp, block_address)) return 0;
    count = knet_read(fp->x.fpr, header, sizeof(header));
#else
    int64_t block_address = ftello(fp->file);
	if (load_block_from_cache(fp, block_address)) return 0;
    count = fread(header, 1, sizeof(header), fp->file);
#endif
    if (count == 0) {
        fp->block_length = 0;
        return 0;
    }
	size = count;
    if (count != sizeof(header)) {
        report_error(fp, "read failed");
        return -1;
    }
    if (!check_header(header)) {
        report_error(fp, "invalid block header");
        return -1;
    }
    int block_length = unpackInt16((uint8_t*)&header[16]) + 1;
    bgzf_byte_t* compressed_block = (bgzf_byte_t*) fp->compressed_block;
    memcpy(compressed_block, header, BLOCK_HEADER_LENGTH);
    int remaining = block_length - BLOCK_HEADER_LENGTH;
#ifdef _USE_KNETFILE
    count = knet_read(fp->x.fpr, &compressed_block[BLOCK_HEADER_LENGTH], remaining);
#else
    count = fread(&compressed_block[BLOCK_HEADER_LENGTH], 1, remaining, fp->file);
#endif
    if (count != remaining) {
        report_error(fp, "read failed");
        return -1;
    }
	size += count;
    count = inflate_block(fp, block_length);
    if (count < 0) return -1;
    if (fp->block_length != 0) {
        // Do not reset offset if this read follows a seek.
        fp->block_offset = 0;
    }
    fp->block_address = block_address;
    fp->block_length = count;
	cache_block(fp, size);
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
    bgzf_byte_t* output = data;
    while (bytes_read < length) {
        int available = fp->block_length - fp->block_offset;
        if (available <= 0) {
            if (bgzf_read_block(fp) != 0) {
                return -1;
            }
            available = fp->block_length - fp->block_offset;
            if (available <= 0) {
                break;
            }
        }
        int copy_length = bgzf_min(length-bytes_read, available);
        bgzf_byte_t* buffer = fp->uncompressed_block;
        memcpy(output, buffer + fp->block_offset, copy_length);
        fp->block_offset += copy_length;
        output += copy_length;
        bytes_read += copy_length;
    }
    if (fp->block_offset == fp->block_length) {
#ifdef _USE_KNETFILE
        fp->block_address = knet_tell(fp->x.fpr);
#else
        fp->block_address = ftello(fp->file);
#endif
        fp->block_offset = 0;
        fp->block_length = 0;
    }
    return bytes_read;
}

int bgzf_flush(BGZF* fp)
{
    while (fp->block_offset > 0) {
        int count, block_length;
		block_length = deflate_block(fp, fp->block_offset);
        if (block_length < 0) return -1;
#ifdef _USE_KNETFILE
        count = fwrite(fp->compressed_block, 1, block_length, fp->x.fpw);
#else
        count = fwrite(fp->compressed_block, 1, block_length, fp->file);
#endif
        if (count != block_length) {
            report_error(fp, "write failed");
            return -1;
        }
        fp->block_address += block_length;
    }
    return 0;
}

int bgzf_flush_try(BGZF *fp, int size)
{
	if (fp->block_offset + size > fp->uncompressed_block_size)
		return bgzf_flush(fp);
	return -1;
}

int bgzf_write(BGZF* fp, const void* data, int length)
{
    if (fp->open_mode != 'w') {
        report_error(fp, "file not open for writing");
        return -1;
    }

    if (fp->uncompressed_block == NULL)
        fp->uncompressed_block = malloc(fp->uncompressed_block_size);

    const bgzf_byte_t* input = data;
    int block_length = fp->uncompressed_block_size;
    int bytes_written = 0;
    while (bytes_written < length) {
        int copy_length = bgzf_min(block_length - fp->block_offset, length - bytes_written);
        bgzf_byte_t* buffer = fp->uncompressed_block;
        memcpy(buffer + fp->block_offset, input, copy_length);
        fp->block_offset += copy_length;
        input += copy_length;
        bytes_written += copy_length;
        if (fp->block_offset == block_length) {
            if (bgzf_flush(fp) != 0) {
                break;
            }
        }
    }
    return bytes_written;
}

int bgzf_close(BGZF* fp)
{
    if (fp->open_mode == 'w') {
        if (bgzf_flush(fp) != 0) return -1;
		{ // add an empty block
			int count, block_length = deflate_block(fp, 0);
#ifdef _USE_KNETFILE
			count = fwrite(fp->compressed_block, 1, block_length, fp->x.fpw);
#else
			count = fwrite(fp->compressed_block, 1, block_length, fp->file);
#endif
		}
#ifdef _USE_KNETFILE
        if (fflush(fp->x.fpw) != 0) {
#else
        if (fflush(fp->file) != 0) {
#endif
            report_error(fp, "flush failed");
            return -1;
        }
    }
    if (fp->owned_file) {
#ifdef _USE_KNETFILE
		int ret;
		if (fp->open_mode == 'w') ret = fclose(fp->x.fpw);
		else ret = knet_close(fp->x.fpr);
        if (ret != 0) return -1;
#else
        if (fclose(fp->file) != 0) return -1;
#endif
    }
    free(fp->uncompressed_block);
    free(fp->compressed_block);
	free_cache(fp);
    free(fp);
    return 0;
}

void bgzf_set_cache_size(BGZF *fp, int cache_size)
{
	if (fp) fp->cache_size = cache_size;
}

int bgzf_check_EOF(BGZF *fp)
{
	static uint8_t magic[28] = "\037\213\010\4\0\0\0\0\0\377\6\0\102\103\2\0\033\0\3\0\0\0\0\0\0\0\0\0";
	uint8_t buf[28];
	off_t offset;
#ifdef _USE_KNETFILE
	offset = knet_tell(fp->x.fpr);
	if (knet_seek(fp->x.fpr, -28, SEEK_END) != 0) return -1;
	knet_read(fp->x.fpr, buf, 28);
	knet_seek(fp->x.fpr, offset, SEEK_SET);
#else
	offset = ftello(fp->file);
	if (fseeko(fp->file, -28, SEEK_END) != 0) return -1;
	fread(buf, 1, 28, fp->file);
	fseeko(fp->file, offset, SEEK_SET);
#endif
	return (memcmp(magic, buf, 28) == 0)? 1 : 0;
}

int64_t bgzf_seek(BGZF* fp, int64_t pos, int where)
{
	int block_offset;
	int64_t block_address;

    if (fp->open_mode != 'r') {
        report_error(fp, "file not open for read");
        return -1;
    }
    if (where != SEEK_SET) {
        report_error(fp, "unimplemented seek option");
        return -1;
    }
    block_offset = pos & 0xFFFF;
    block_address = (pos >> 16) & 0xFFFFFFFFFFFFLL;
#ifdef _USE_KNETFILE
    if (knet_seek(fp->x.fpr, block_address, SEEK_SET) != 0) {
#else
    if (fseeko(fp->file, block_address, SEEK_SET) != 0) {
#endif
        report_error(fp, "seek failed");
        return -1;
    }
    fp->block_length = 0;  // indicates current block is not loaded
    fp->block_address = block_address;
    fp->block_offset = block_offset;
    return 0;
}
