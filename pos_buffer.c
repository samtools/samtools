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

#include "pos_buffer.h"
#include <assert.h>
#include <stdlib.h>
#include <string.h>

struct pos_list;
typedef struct pos_list pos_list_t;

struct pos_list {
	pos_list_t* next;
	uint32_t ltid;
	uint32_t lpos;
	uint32_t rpos;
	int score;
	char* name;
};

struct pos_tree;
typedef struct pos_tree pos_tree_t;

struct pos_tree {
	pos_tree_t* left;
	pos_tree_t* right;
	uint32_t ltid;
	uint32_t lpos;
	int score;
	char* name;
};

static const int POS_BUFFER_LENGTH = 10000;
struct pos_buffer {
	uint32_t tid; // rightmost tid
	uint32_t base_pos; // rightmost pos base
	pos_list_t* aux; // Overspill here
	size_t buffer_base; //
	pos_tree_t* right_most[POS_BUFFER_LENGTH]; // Ring buffer
};


static void reap_tree(pos_tree_t* tree)
{
	if (tree == NULL) return;
	reap_tree(tree->left);
	reap_tree(tree->right);
	free(tree->left);
	free(tree->right);
	free(tree->name);
}

pos_buffer_t* pos_buffer_init()
{
	pos_buffer_t* retval = (pos_buffer_t*) calloc(1, sizeof(pos_buffer_t));
	return retval;
}

static void pos_buffer_clear(pos_buffer_t* buf)
{
	pos_list_t* ptr;
	for (ptr = buf->aux; ptr != NULL; ptr = buf->aux) {
		buf->aux = ptr->next;
		free(ptr->name);
		free(ptr);
	}
	size_t i;
	for ( i = 0; i < POS_BUFFER_LENGTH; ++i ) {
		reap_tree(buf->right_most[i]);
		buf->right_most[i] = NULL;
	}
}

void pos_buffer_destroy(pos_buffer_t* target)
{
	pos_buffer_clear(target);
	free(target);
}

// Call whenever rtid changes
void pos_buffer_reset(pos_buffer_t* buf)
{
	pos_buffer_clear(buf);
	buf->base_pos = 0;
	buf->buffer_base = 0;
}

static size_t pos_buff_rpos(const pos_buffer_t* buf, uint32_t rpos) {
	assert(buf->base_pos >= rpos);
	size_t offset = rpos - buf->base_pos;
	size_t buffer_offset = (buf->buffer_base + offset) % POS_BUFFER_LENGTH;
	return buffer_offset;
}

static void pos_buffer_load_from_aux(pos_buffer_t* buf)
{
	pos_list_t* ptr;
	pos_list_t** prev = &buf->aux;
	for (ptr = buf->aux; ptr != NULL; ptr = buf->aux) {
		if (buf->base_pos <= ptr->rpos
			&& ptr->rpos < buf->base_pos + POS_BUFFER_LENGTH ) {
			pos_buffer_insert(buf, ptr->ltid, ptr->lpos, buf->tid, ptr->rpos, ptr->score, ptr->name); // TODO: implement direct insert
			*prev = ptr->next;
			free(ptr->name);
			free(ptr);
		} else {
			prev = &ptr->next;
		}
	}
	
}

void pos_buffer_advance(pos_buffer_t* buf, uint32_t lpos)
{
	if (lpos < buf->base_pos + POS_BUFFER_LENGTH) {
		// Advance main buffer
		uint32_t advance = lpos - buf->base_pos;
		buf->base_pos = lpos;
		size_t i;
		for ( i = 0; i < advance; ++i ) {
			size_t j = pos_buff_rpos(buf, i);
			reap_tree(buf->right_most[j]);
			buf->right_most[j] = NULL;
		}
		buf->buffer_base += advance;
	} else {
		pos_buffer_reset(buf);
		buf->base_pos = lpos;
	}
	pos_buffer_load_from_aux(buf);
}

static char* pos_aux_insert(pos_buffer_t* buf, uint32_t ltid, uint32_t lpos, uint32_t rpos, int score, const char* name)
{
	// search to see if entry exists first
	pos_list_t* ptr;
	for (ptr = buf->aux; ptr != NULL; ptr = buf->aux) {
		if (ptr->ltid == ltid && ptr->lpos == lpos && ptr->rpos == rpos) {
			if (score > ptr->score) {
				return strdup(name);
			} else {
				return strdup(ptr->name);
			}
		}
	}
	// not found
	pos_list_t* insert = (pos_list_t*)malloc(sizeof(pos_list_t));
	insert->ltid = ltid;
	insert->lpos = lpos;
	insert->rpos = rpos;
	insert->score = score;
	insert->name = strdup(name);
	insert->next = buf->aux;
	buf->aux = insert;
	return NULL;
}

static char* pos_tree_insert(pos_tree_t* tree, uint32_t ltid, uint32_t lpos, int score, const char* name)
{
	// search base to see if entry exists first
	if (tree->ltid == ltid && tree->lpos == lpos) {
		if ( score > tree->score ) {
			return strdup(name);
		} else {
			return strdup(tree->name);
		}
	} else { // Check left and right
		if ( ltid < tree->ltid || (ltid == tree->ltid && lpos < tree->lpos) ) {
			if (tree->left != NULL) {
				return pos_tree_insert(tree->left, ltid, lpos, score, name);
			} else {
				pos_tree_t* insert = (pos_tree_t*)malloc(sizeof(pos_list_t));
				insert->ltid = ltid;
				insert->lpos = lpos;
				insert->score = score;
				insert->name = strdup(name);
				insert->left = NULL;
				insert->right = NULL;
				tree->left = insert;
			}
		} else { // must be right
			if (tree->right != NULL) {
				return pos_tree_insert(tree->right, ltid, lpos, score, name);
			} else {
				pos_tree_t* insert = (pos_tree_t*)malloc(sizeof(pos_list_t));
				insert->ltid = ltid;
				insert->lpos = lpos;
				insert->score = score;
				insert->name = strdup(name);
				insert->left = NULL;
				insert->right = NULL;
				tree->right = insert;
			}
		}
	}
	return NULL;
}

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
char* pos_buffer_insert(pos_buffer_t* buf, uint32_t ltid, uint32_t lpos, uint32_t rtid, uint32_t rpos, int score, const char* name)
{
	if (rtid != buf->tid) {
		pos_buffer_reset(buf);
		buf->tid = rtid;
	}
	assert(buf->base_pos <= rpos);
	if (rpos > (buf->base_pos + POS_BUFFER_LENGTH)) {
		// it's outside our buffer
		return pos_aux_insert(buf, ltid, lpos, rpos, score, name);
	} else {
		pos_tree_t* tree = buf->right_most[pos_buff_rpos(buf, rpos)];
		if ( tree != NULL) {
			return pos_tree_insert(tree, ltid, lpos, score, name);
		} else {
			// create new tree
			pos_tree_t* insert = (pos_tree_t*)malloc(sizeof(pos_list_t));
			insert->ltid = ltid;
			insert->lpos = lpos;
			insert->score = score;
			insert->name = strdup(name);
			insert->left = NULL;
			insert->right = NULL;
			buf->right_most[pos_buff_rpos(buf, rpos)] = insert;
		}
	}
	return NULL;
}
