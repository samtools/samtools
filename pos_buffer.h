/* The MIT License
 *
 * Copyright (c) 2014 Genome Research Limited.
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
 * BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
 * ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

/* Contact: Martin Pollard <mp15@sanger.ac.uk> */

#ifndef POS_BUFFER_H
#define POS_BUFFER_H

#include <inttypes.h>

struct pos_buffer;
typedef struct pos_buffer pos_buffer_t;

// Initialises a new pos buffer
pos_buffer_t* pos_buffer_init();
/**
 * Attempt to insert a pair of reads into the pos buffer
 *
 * @param buf the pos buffer
 * @param ltid leftmost target id of pair
 * @param lpos leftmost position of pair
 * @param rtid rightmost target id of pair
 * @param rpos rightmost position of pair
 * @param score score of pair
 * @param name read name of pair
 * @return name of read to kill, or NULL if it's new insert
 */
char* pos_buffer_insert(pos_buffer_t* buf, uint32_t ltid, uint32_t lpos, uint32_t rtid, uint32_t rpos, int score, const char* name);
// Frees memory associated with a pos buffer
void pos_buffer_destroy(pos_buffer_t* target);

#endif