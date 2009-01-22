/*
 * RAZF : Random Access compressed(Z) File
 * Version: 1.0
 * Release Date: 2008-10-27
 *
 * Copyright 2008, Jue Ruan <ruanjue@gmail.com>, Heng Li <lh3@sanger.ac.uk>
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 *
 * To compile razf.c, zlib-1.2.3(or greater) is required.
 */

#ifndef _NO_RAZF

#include <fcntl.h>
#include <stdio.h>
#include "razf.h"

static inline uint32_t byte_swap_4(uint32_t v){
	v = ((v & 0x0000FFFFU) << 16) | (v >> 16);
	return ((v & 0x00FF00FFU) << 8) | ((v & 0xFF00FF00U) >> 8);
}

static inline uint64_t byte_swap_8(uint64_t v){
	v = ((v & 0x00000000FFFFFFFFLLU) << 32) | (v >> 32);
	v = ((v & 0x0000FFFF0000FFFFLLU) << 16) | ((v & 0xFFFF0000FFFF0000LLU) >> 16);
	return ((v & 0x00FF00FF00FF00FFLLU) << 8) | ((v & 0xFF00FF00FF00FF00LLU) >> 8);
}

static inline int is_big_endian(){
	int x = 0x01;
	char *c = (char*)&x;
	return (c[0] != 0x01);
}

static void add_zindex(RAZF *rz, int64_t in, int64_t out){
	if(rz->index->size == rz->index->cap){
		rz->index->cap = rz->index->cap * 1.5 + 2;
		rz->index->cell_offsets = realloc(rz->index->cell_offsets, sizeof(int) * rz->index->cap);
		rz->index->bin_offsets  = realloc(rz->index->bin_offsets, sizeof(int64_t) * (rz->index->cap/RZ_BIN_SIZE + 1));
	}
	if(rz->index->size % RZ_BIN_SIZE == 0) rz->index->bin_offsets[rz->index->size / RZ_BIN_SIZE] = out;
	rz->index->cell_offsets[rz->index->size] = out - rz->index->bin_offsets[rz->index->size / RZ_BIN_SIZE];
	rz->index->size ++;
}

static void save_zindex(RAZF *rz, int fd){
	int32_t i, v32;
	int is_be;
	is_be = is_big_endian();
	if(is_be) write(fd, &rz->index->size, sizeof(int));
	else {
		v32 = byte_swap_4((uint32_t)rz->index->size);
		write(fd, &v32, sizeof(uint32_t));
	}
	v32 = rz->index->size / RZ_BIN_SIZE + 1;
	if(!is_be){
		for(i=0;i<v32;i++) rz->index->bin_offsets[i]  = byte_swap_8((uint64_t)rz->index->bin_offsets[i]);
		for(i=0;i<rz->index->size;i++) rz->index->cell_offsets[i] = byte_swap_4((uint32_t)rz->index->cell_offsets[i]);
	}
	write(fd, rz->index->bin_offsets, sizeof(int64_t) * v32);
	write(fd, rz->index->cell_offsets, sizeof(int32_t) * rz->index->size);
}

static void load_zindex(RAZF *rz, int fd){
	int32_t i, v32;
	int is_be;
	if(!rz->load_index) return;
	if(rz->index == NULL) rz->index = malloc(sizeof(ZBlockIndex));
	is_be = is_big_endian();
	read(fd, &rz->index->size, sizeof(int));
	if(!is_be) rz->index->size = byte_swap_4((uint32_t)rz->index->size);
	rz->index->cap = rz->index->size;
	v32 = rz->index->size / RZ_BIN_SIZE + 1;
	rz->index->bin_offsets  = malloc(sizeof(int64_t) * v32);
	read(fd, rz->index->bin_offsets, sizeof(int64_t) * v32);
	rz->index->cell_offsets = malloc(sizeof(int) * rz->index->size);
	read(fd, rz->index->cell_offsets, sizeof(int) * rz->index->size);
	if(!is_be){
		for(i=0;i<v32;i++) rz->index->bin_offsets[i] = byte_swap_8((uint64_t)rz->index->bin_offsets[i]);
		for(i=0;i<rz->index->size;i++) rz->index->cell_offsets[i] = byte_swap_4((uint32_t)rz->index->cell_offsets[i]);
	}
}

static RAZF* razf_open_w(int fd){
	RAZF *rz;
	rz = calloc(1, sizeof(RAZF));
	rz->mode = 'w';
	rz->filedes = fd;
	rz->stream = calloc(sizeof(z_stream), 1);
	rz->inbuf  = malloc(RZ_BUFFER_SIZE);
	rz->outbuf = malloc(RZ_BUFFER_SIZE);
	rz->index = calloc(sizeof(ZBlockIndex), 1);
	deflateInit2(rz->stream, RZ_COMPRESS_LEVEL, Z_DEFLATED, WINDOW_BITS + 16, DEF_MEM_LEVEL, Z_DEFAULT_STRATEGY);
	rz->stream->avail_out = RZ_BUFFER_SIZE;
	rz->stream->next_out  = rz->outbuf;
	rz->header = calloc(sizeof(gz_header), 1);
	rz->header->os    = 0x03; //Unix
	rz->header->text  = 0;
	rz->header->time  = 0;
	rz->header->extra = malloc(7);
	strncpy((char*)rz->header->extra, "RAZF", 4);
	rz->header->extra[4] = 1; // obsolete field
	// block size = RZ_BLOCK_SIZE, Big-Endian
	rz->header->extra[5] = RZ_BLOCK_SIZE >> 8;
	rz->header->extra[6] = RZ_BLOCK_SIZE & 0xFF;
	rz->header->extra_len = 7;
	rz->header->name = rz->header->comment  = 0;
	rz->header->hcrc = 0;
	deflateSetHeader(rz->stream, rz->header);
	rz->block_pos = rz->block_off = 0;
	return rz;
}

static void _razf_write(RAZF* rz, const void *data, int size){
	int tout;
	rz->stream->avail_in = size;
	rz->stream->next_in  = (void*)data;
	while(1){
		tout = rz->stream->avail_out;
		deflate(rz->stream, Z_NO_FLUSH);
		rz->out += tout - rz->stream->avail_out;
		if(rz->stream->avail_out) break;
		write(rz->filedes, rz->outbuf, RZ_BUFFER_SIZE - rz->stream->avail_out);
		rz->stream->avail_out = RZ_BUFFER_SIZE;
		rz->stream->next_out  = rz->outbuf;
		if(rz->stream->avail_in == 0) break;
	};
	rz->in += size - rz->stream->avail_in;
	rz->block_off += size - rz->stream->avail_in;
}

static void razf_flush(RAZF *rz){
	uint32_t tout;
	if(rz->buf_len){
		_razf_write(rz, rz->inbuf, rz->buf_len);
		rz->buf_off = rz->buf_len = 0;
	}
	if(rz->stream->avail_out){
		write(rz->filedes, rz->outbuf, RZ_BUFFER_SIZE - rz->stream->avail_out);
		rz->stream->avail_out = RZ_BUFFER_SIZE;
		rz->stream->next_out  = rz->outbuf;
	}
	while(1){
		tout = rz->stream->avail_out;
		deflate(rz->stream, Z_FULL_FLUSH);
		rz->out += tout - rz->stream->avail_out;
		if(rz->stream->avail_out == 0){
			write(rz->filedes, rz->outbuf, RZ_BUFFER_SIZE - rz->stream->avail_out);
			rz->stream->avail_out = RZ_BUFFER_SIZE;
			rz->stream->next_out  = rz->outbuf;
		} else break;
	}
	rz->block_pos = rz->out;
	rz->block_off = 0;
}

static void razf_end_flush(RAZF *rz){
	uint32_t tout;
	if(rz->buf_len){
		_razf_write(rz, rz->inbuf, rz->buf_len);
		rz->buf_off = rz->buf_len = 0;
	}
	while(1){
		tout = rz->stream->avail_out;
		deflate(rz->stream, Z_FINISH);
		rz->out += tout - rz->stream->avail_out;
		if(rz->stream->avail_out < RZ_BUFFER_SIZE){
			write(rz->filedes, rz->outbuf, RZ_BUFFER_SIZE - rz->stream->avail_out);
			rz->stream->avail_out = RZ_BUFFER_SIZE;
			rz->stream->next_out  = rz->outbuf;
		} else break;
	}
}

static void _razf_buffered_write(RAZF *rz, const void *data, int size){
	int i, n;
	while(1){
		if(rz->buf_len == RZ_BUFFER_SIZE){
			_razf_write(rz, rz->inbuf, rz->buf_len);
			rz->buf_len = 0;
		}
		if(size + rz->buf_len < RZ_BUFFER_SIZE){
			for(i=0;i<size;i++) ((char*)rz->inbuf + rz->buf_len)[i] = ((char*)data)[i];
			rz->buf_len += size;
			return;
		} else {
			n = RZ_BUFFER_SIZE - rz->buf_len;
			for(i=0;i<n;i++) ((char*)rz->inbuf + rz->buf_len)[i] = ((char*)data)[i];
			size -= n;
			data += n;
			rz->buf_len += n;
		}
	}
}

int razf_write(RAZF* rz, const void *data, int size){
	int ori_size, n;
	int64_t next_block;
	ori_size = size;
	next_block = ((rz->in / RZ_BLOCK_SIZE) + 1) * RZ_BLOCK_SIZE;
	while(rz->in + rz->buf_len + size >= next_block){
		n = next_block - rz->in - rz->buf_len;
		_razf_buffered_write(rz, data, n);
		data += n;
		size -= n;
		razf_flush(rz);
		add_zindex(rz, rz->in, rz->out);
		next_block = ((rz->in / RZ_BLOCK_SIZE) + 1) * RZ_BLOCK_SIZE;
	}
	_razf_buffered_write(rz, data, size);
	return ori_size;
}

/* gzip flag byte */
#define ASCII_FLAG   0x01 /* bit 0 set: file probably ascii text */
#define HEAD_CRC     0x02 /* bit 1 set: header CRC present */
#define EXTRA_FIELD  0x04 /* bit 2 set: extra field present */
#define ORIG_NAME    0x08 /* bit 3 set: original file name present */
#define COMMENT      0x10 /* bit 4 set: file comment present */
#define RESERVED     0xE0 /* bits 5..7: reserved */

static int _read_gz_header(unsigned char *data, int size, int *extra_off, int *extra_len){
	int method, flags, n, len;
	if(size < 2) return 0;
	if(data[0] != 0x1f || data[1] != 0x8b) return 0;
	if(size < 4) return 0;
	method = data[2];
	flags  = data[3];
	if(method != Z_DEFLATED || (flags & RESERVED)) return 0;
	n = 4 + 6; // Skip 6 bytes
	*extra_off = n + 2;
	*extra_len = 0;
	if(flags & EXTRA_FIELD){
		if(size < n + 2) return 0;
		len = ((int)data[n + 1] << 8) | data[n];
		n += 2;
		*extra_off = n;
		while(len){
			if(n >= size) return 0;
			n ++;
			len --;
		}
		*extra_len = n - (*extra_off);
	}
	if(flags & ORIG_NAME) while(n < size && data[n++]);
	if(flags & COMMENT) while(n < size && data[n++]);
	if(flags & HEAD_CRC){
		if(n + 2 > size) return 0;
		n += 2;
	}
	return n;
}

static RAZF* razf_open_r(int fd, int _load_index){
	RAZF *rz;
	int ext_off, ext_len;
	int n, is_be, ret;
	int64_t end;
	unsigned char c[] = "RAZF";
	rz = calloc(1, sizeof(RAZF));
	rz->mode = 'r';
	rz->filedes = fd;
	rz->stream = calloc(sizeof(z_stream), 1);
	rz->inbuf  = malloc(RZ_BUFFER_SIZE);
	rz->outbuf = malloc(RZ_BUFFER_SIZE);
	rz->end = rz->src_end = 0x7FFFFFFFFFFFFFFFLL;
	n = read(rz->filedes, rz->inbuf, RZ_BUFFER_SIZE);
	ret = _read_gz_header(rz->inbuf, n, &ext_off, &ext_len);
	if(ret == 0){
		PLAIN_FILE:
		rz->in = n;
		rz->file_type = FILE_TYPE_PLAIN;
		memcpy(rz->outbuf, rz->inbuf, n);
		rz->buf_len = n;
		free(rz->stream);
		rz->stream = NULL;
		return rz;
	}
	rz->header_size = ret;
	ret = inflateInit2(rz->stream, -WINDOW_BITS);
	if(ret != Z_OK){ inflateEnd(rz->stream); goto PLAIN_FILE;}
	rz->stream->avail_in = n - rz->header_size;
	rz->stream->next_in  = rz->inbuf + rz->header_size;
	rz->stream->avail_out = RZ_BUFFER_SIZE;
	rz->stream->next_out  = rz->outbuf;
	rz->file_type = FILE_TYPE_GZ;
	rz->in = rz->header_size;
	rz->block_pos = rz->header_size;
	rz->next_block_pos = rz->header_size;
	rz->block_off = 0;
	if(ext_len < 7 || memcmp(rz->inbuf + ext_off, c, 4) != 0) return rz;
	if(((((unsigned char*)rz->inbuf)[ext_off + 5] << 8) | ((unsigned char*)rz->inbuf)[ext_off + 6]) != RZ_BLOCK_SIZE){
		fprintf(stderr, " -- WARNING: RZ_BLOCK_SIZE is not %d, treat source as gz file.  in %s -- %s:%d --\n", RZ_BLOCK_SIZE, __FUNCTION__, __FILE__, __LINE__);
		return rz;
	}
	rz->load_index = _load_index;
	rz->file_type = FILE_TYPE_RZ;
	if(lseek(fd, -16, SEEK_END) == -1){
		UNSEEKABLE:
		rz->seekable = 0;
		rz->index = NULL;
		rz->src_end = rz->end = 0x7FFFFFFFFFFFFFFFLL;
	} else {
		is_be = is_big_endian();
		rz->seekable = 1;
		read(fd, &end, sizeof(int64_t));
		if(!is_be) rz->src_end = (int64_t)byte_swap_8((uint64_t)end);
		else rz->src_end = end;
		read(fd, &end, sizeof(int64_t));
		if(!is_be) rz->end = (int64_t)byte_swap_8((uint64_t)end);
		else rz->end = end;
		if(n > rz->end){
			rz->stream->avail_in -= n - rz->end;
			n = rz->end;
		}
		if(rz->end > rz->src_end){
			lseek(fd, rz->in, SEEK_SET);
			goto UNSEEKABLE;
		}
		if(lseek(fd, rz->end, SEEK_SET) != rz->end){
			lseek(fd, rz->in, SEEK_SET);
			goto UNSEEKABLE;
		}
		load_zindex(rz, fd);
		lseek(fd, n, SEEK_SET);
	}
	return rz;
}

RAZF* razf_dopen(int fd, const char *mode){
	if(strcasecmp(mode, "r") == 0) return razf_open_r(fd, 1);
	else if(strcasecmp(mode, "w") == 0) return razf_open_w(fd);
	else return NULL;
}

RAZF* razf_dopen2(int fd, const char *mode)
{
	if(strcasecmp(mode, "r") == 0) return razf_open_r(fd, 0);
	else if(strcasecmp(mode, "w") == 0) return razf_open_w(fd);
	else return NULL;
}

static inline RAZF* _razf_open(const char *filename, const char *mode, int _load_index){
	int fd;
	RAZF *rz;
	if(strcasecmp(mode, "r") == 0){
		fd = open(filename, O_RDONLY);
		rz = razf_open_r(fd, _load_index);
	} else if(strcasecmp(mode, "w") == 0){
		fd = open(filename, O_WRONLY | O_CREAT | O_TRUNC, 0644);
		rz = razf_open_w(fd);
	} else return NULL;
	return rz;
}

RAZF* razf_open(const char *filename, const char *mode){
	return _razf_open(filename, mode, 1);
}

RAZF* razf_open2(const char *filename, const char *mode){
	return _razf_open(filename, mode, 0);
}

int razf_get_data_size(RAZF *rz, int64_t *u_size, int64_t *c_size){
	int64_t n;
	if(rz->mode != 'r' && rz->mode != 'R') return 0;
	switch(rz->file_type){
		case FILE_TYPE_PLAIN:
			if(rz->end == 0x7fffffffffffffffLL){
				if((n = lseek(rz->filedes, 0, SEEK_CUR)) == -1) return 0;
				rz->end = lseek(rz->filedes, 0, SEEK_END);
				lseek(rz->filedes, n, SEEK_SET);
			}
			*u_size = *c_size = rz->end;
			return 1;
		case FILE_TYPE_GZ:
			return 0;
		case FILE_TYPE_RZ:
			if(rz->src_end == rz->end) return 0;
			*u_size = rz->src_end;
			*c_size = rz->end;
			return 1;
		default:
			return 0;
	}
}

static int _razf_read(RAZF* rz, void *data, int size){
	int ret, tin;
	if(rz->z_eof || rz->z_err) return 0;
	if (rz->file_type == FILE_TYPE_PLAIN) {
		ret = read(rz->filedes, data, size);
		if (ret == 0) rz->z_eof = 1;
		return ret;
	}
	rz->stream->avail_out = size;
	rz->stream->next_out  = data;
	while(rz->stream->avail_out){
		if(rz->stream->avail_in == 0){
			if(rz->in >= rz->end){ rz->z_eof = 1; break; }
			if(rz->end - rz->in < RZ_BUFFER_SIZE){
				rz->stream->avail_in = read(rz->filedes, rz->inbuf, rz->end -rz->in);
			} else {
				rz->stream->avail_in = read(rz->filedes, rz->inbuf, RZ_BUFFER_SIZE);
			}
			if(rz->stream->avail_in == 0){
				rz->z_eof = 1;
				break;
			}
			rz->stream->next_in = rz->inbuf;
		}
		tin = rz->stream->avail_in;
		ret = inflate(rz->stream, Z_BLOCK);
		rz->in += tin - rz->stream->avail_in;
		if(ret == Z_NEED_DICT || ret == Z_MEM_ERROR || ret == Z_DATA_ERROR){
			fprintf(stderr, "[_razf_read] inflate error: %d (at %s:%d)\n", ret, __FILE__, __LINE__);
			rz->z_err = 1;
			break;
		}
		if(ret == Z_STREAM_END){
			rz->z_eof = 1;
			break;
		}
		if ((rz->stream->data_type&128) && !(rz->stream->data_type&64)){
			rz->buf_flush = 1;
			rz->next_block_pos = rz->in;
			break;
		}
	}
	return size - rz->stream->avail_out;
}

int razf_read(RAZF *rz, void *data, int size){
	int ori_size, i;
	ori_size = size;
	while(size > 0){
		if(rz->buf_len){
			if(size < rz->buf_len){
				for(i=0;i<size;i++) ((char*)data)[i] = ((char*)rz->outbuf + rz->buf_off)[i];
				rz->buf_off += size;
				rz->buf_len -= size;
				data += size;
				rz->block_off += size;
				size = 0;
				break;
			} else {
				for(i=0;i<rz->buf_len;i++) ((char*)data)[i] = ((char*)rz->outbuf + rz->buf_off)[i];
				data += rz->buf_len;
				size -= rz->buf_len;
				rz->block_off += rz->buf_len;
				rz->buf_off = 0;
				rz->buf_len = 0;
				if(rz->buf_flush){
					rz->block_pos = rz->next_block_pos;
					rz->block_off = 0;
					rz->buf_flush = 0;
				}
			}
		} else if(rz->buf_flush){
			rz->block_pos = rz->next_block_pos;
			rz->block_off = 0;
			rz->buf_flush = 0;
		}
		if(rz->buf_flush) continue;
		rz->buf_len = _razf_read(rz, rz->outbuf, RZ_BUFFER_SIZE);
		if(rz->z_eof && rz->buf_len == 0) break;
	}
	rz->out += ori_size - size;
	return ori_size - size;
}

int razf_skip(RAZF* rz, int size){
	int ori_size;
	ori_size = size;
	while(size > 0){
		if(rz->buf_len){
			if(size < rz->buf_len){
				rz->buf_off += size;
				rz->buf_len -= size;
				rz->block_off += size;
				size = 0;
				break;
			} else {
				size -= rz->buf_len;
				rz->buf_off = 0;
				rz->buf_len = 0;
				rz->block_off += rz->buf_len;
				if(rz->buf_flush){
					rz->block_pos = rz->next_block_pos;
					rz->block_off = 0;
					rz->buf_flush = 0;
				}
			}
		} else if(rz->buf_flush){
			rz->block_pos = rz->next_block_pos;
			rz->block_off = 0;
			rz->buf_flush = 0;
		}
		if(rz->buf_flush) continue;
		rz->buf_len = _razf_read(rz, rz->outbuf, RZ_BUFFER_SIZE);
		if(rz->z_eof) break;
	}
	rz->out += ori_size - size;
	return ori_size - size;
}

static void _razf_reset_read(RAZF *rz, int64_t in, int64_t out){
	lseek(rz->filedes, in, SEEK_SET);
	rz->in  = in;
	rz->out = out;
	rz->block_pos = in;
	rz->next_block_pos = in;
	rz->block_off = 0;
	rz->buf_flush = 0;
	rz->z_eof = rz->z_err = 0;
	inflateReset(rz->stream);
	rz->stream->avail_in = 0;
	rz->buf_off = rz->buf_len = 0;
}

int64_t razf_jump(RAZF *rz, int64_t block_start, int block_offset){
	int64_t pos;
	rz->z_eof = 0;
	if(rz->file_type == FILE_TYPE_PLAIN){
		rz->buf_off = rz->buf_len = 0;
		pos = block_start + block_offset;
		pos = lseek(rz->filedes, pos, SEEK_SET);
		rz->out = rz->in = pos;
		return pos;
	}
	if(block_start == rz->block_pos && block_offset >= rz->block_off) {
		block_offset -= rz->block_off;
		goto SKIP; // Needn't reset inflate
	}
	if(block_start  == 0) block_start = rz->header_size; // Automaticly revist wrong block_start
	_razf_reset_read(rz, block_start, 0);
	SKIP:
	if(block_offset) razf_skip(rz, block_offset);
	return rz->block_off;
}

int64_t razf_seek(RAZF* rz, int64_t pos, int where){
	int64_t idx;
	int64_t seek_pos, new_out;
	rz->z_eof = 0;
	if (where == SEEK_CUR) pos += rz->out;
	else if (where == SEEK_END) pos += rz->src_end;
	if(rz->file_type == FILE_TYPE_PLAIN){
		seek_pos = lseek(rz->filedes, pos, SEEK_SET);
		rz->buf_off = rz->buf_len = 0;
		rz->out = rz->in = seek_pos;
		return seek_pos;
	} else if(rz->file_type == FILE_TYPE_GZ){
		if(pos >= rz->out) goto SKIP;
		return rz->out;
	}
	if(pos == rz->out) return pos;
	if(pos > rz->src_end) return rz->out;
	if(!rz->seekable || !rz->load_index){
		if(pos >= rz->out) goto SKIP;
	}
	idx = pos / RZ_BLOCK_SIZE - 1;
	seek_pos = (idx < 0)? rz->header_size:(rz->index->cell_offsets[idx] + rz->index->bin_offsets[idx / RZ_BIN_SIZE]);
	new_out  = (idx + 1) * RZ_BLOCK_SIZE;
	if(pos > rz->out && new_out <= rz->out) goto SKIP;
	_razf_reset_read(rz, seek_pos, new_out);
	SKIP:
	razf_skip(rz, (int)(pos - rz->out));
	return rz->out;
}

uint64_t razf_tell2(RAZF *rz)
{
	/*
	if (rz->load_index) {
		int64_t idx, seek_pos;
		idx = rz->out / RZ_BLOCK_SIZE - 1;
		seek_pos = (idx < 0)? rz->header_size:(rz->index->cell_offsets[idx] + rz->index->bin_offsets[idx / RZ_BIN_SIZE]);
		if (seek_pos != rz->block_pos || rz->out%RZ_BLOCK_SIZE != rz->block_off)
			fprintf(stderr, "[razf_tell2] inconsistent block offset: (%lld, %lld) != (%lld, %lld)\n",
					(long long)seek_pos, (long long)rz->out%RZ_BLOCK_SIZE, (long long)rz->block_pos, (long long) rz->block_off);
	}
	*/
	return (uint64_t)rz->block_pos<<16 | (rz->block_off&0xffff);
}

int64_t razf_seek2(RAZF *rz, uint64_t voffset, int where)
{
	if (where != SEEK_SET) return -1;
	return razf_jump(rz, voffset>>16, voffset&0xffff);
}

void razf_close(RAZF *rz){
	uint64_t v64;
	if(rz->mode == 'w'){
		razf_end_flush(rz);
		deflateEnd(rz->stream);
		save_zindex(rz, rz->filedes);
		if(is_big_endian()){
			write(rz->filedes, &rz->in, sizeof(int64_t));
			write(rz->filedes, &rz->out, sizeof(int64_t));
		} else {
			v64 = byte_swap_8((uint64_t)rz->in);
			write(rz->filedes, &v64, sizeof(int64_t));
			v64 = byte_swap_8((uint64_t)rz->out);
			write(rz->filedes, &v64, sizeof(int64_t));
		}
	} else if(rz->mode == 'r'){
		if(rz->stream) inflateEnd(rz->stream);
	}
	if(rz->inbuf) free(rz->inbuf);
	if(rz->outbuf) free(rz->outbuf);
	if(rz->header){
		free(rz->header->extra);
		free(rz->header->name);
		free(rz->header->comment);
		free(rz->header);
	}
	if(rz->index){
		free(rz->index->bin_offsets);
		free(rz->index->cell_offsets);
		free(rz->index);
	}
	free(rz->stream);
	close(rz->filedes);
	free(rz);
}

#endif
