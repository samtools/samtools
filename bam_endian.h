/*  bam_endian.h -- endianness conversion functions.

    Copyright (C) 2008 Genome Research Ltd.

    Author: Heng Li <lh3@sanger.ac.uk>

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
DEALINGS IN THE SOFTWARE.  */

#ifndef BAM_ENDIAN_H
#define BAM_ENDIAN_H

#include <stdint.h>

static inline int bam_is_big_endian()
{
    long one= 1;
    return !(*((char *)(&one)));
}
static inline uint16_t bam_swap_endian_2(uint16_t v)
{
    return (uint16_t)(((v & 0x00FF00FFU) << 8) | ((v & 0xFF00FF00U) >> 8));
}
static inline void *bam_swap_endian_2p(void *x)
{
    *(uint16_t*)x = bam_swap_endian_2(*(uint16_t*)x);
    return x;
}
static inline uint32_t bam_swap_endian_4(uint32_t v)
{
    v = ((v & 0x0000FFFFU) << 16) | (v >> 16);
    return ((v & 0x00FF00FFU) << 8) | ((v & 0xFF00FF00U) >> 8);
}
static inline void *bam_swap_endian_4p(void *x)
{
    *(uint32_t*)x = bam_swap_endian_4(*(uint32_t*)x);
    return x;
}
static inline uint64_t bam_swap_endian_8(uint64_t v)
{
    v = ((v & 0x00000000FFFFFFFFLLU) << 32) | (v >> 32);
    v = ((v & 0x0000FFFF0000FFFFLLU) << 16) | ((v & 0xFFFF0000FFFF0000LLU) >> 16);
    return ((v & 0x00FF00FF00FF00FFLLU) << 8) | ((v & 0xFF00FF00FF00FF00LLU) >> 8);
}
static inline void *bam_swap_endian_8p(void *x)
{
    *(uint64_t*)x = bam_swap_endian_8(*(uint64_t*)x);
    return x;
}

#endif
