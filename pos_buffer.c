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
#include "read_vector.h"
#include <assert.h>
#include <stdlib.h>
#include <string.h>

static char* pos_buffer_insert_inner(pos_buffer_t* buf, uint32_t ltid, uint32_t lpos, uint32_t rtid, uint32_t rpos, int score, const char* name);

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

#define POS_BUFFER_LENGTH 10000
struct pos_buffer {
	uint32_t tid; // rightmost tid
	uint32_t base_pos; // rightmost pos base
	pos_list_t* aux; // Overspill reads here
	size_t buffer_base; //
	pos_tree_t* right_most[POS_BUFFER_LENGTH]; // Ring buffer
};

// Deallocate a tree and all subtrees recursively
static void reap_tree(pos_tree_t* tree)
{
	if (tree == NULL) return;
	reap_tree(tree->left);
	reap_tree(tree->right);
	free(tree->left);
	free(tree->right);
	free(tree->name);
}

// Initialise a pos_buffer
pos_buffer_t* pos_buffer_init()
{
	pos_buffer_t* retval = (pos_buffer_t*) calloc(1, sizeof(pos_buffer_t));
	return retval;
}

// Free all used memory inside a pos_buffer
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
		free(buf->right_most[i]);
		buf->right_most[i] = NULL;
	}
}

// Free memory associated with a pos_buffer
void pos_buffer_destroy(pos_buffer_t* target)
{
	pos_buffer_clear(target);
	free(target);
}

// Called whenever rtid changes
static void pos_buffer_reset(pos_buffer_t* buf)
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

/**
 * Load any entries in aux which are now within the ring buffer's bounds into the ring buffer
 */
static void pos_buffer_load_from_aux(pos_buffer_t* buf)
{
	pos_list_t* ptr;
	pos_list_t** prev = &buf->aux;
	for (ptr = buf->aux; ptr != NULL; ptr = buf->aux) {
		if (buf->base_pos <= ptr->rpos
			&& ptr->rpos < buf->base_pos + POS_BUFFER_LENGTH ) {
			pos_buffer_insert_inner(buf, ptr->ltid, ptr->lpos, buf->tid, ptr->rpos, ptr->score, ptr->name); // TODO: implement direct insert
			*prev = ptr->next;
			free(ptr->name);
			free(ptr);
		} else {
			prev = &ptr->next;
		}
	}
	
}

/**
 * Advance the pos buffer within a tid
 */
static void pos_buffer_advance(pos_buffer_t* buf, uint32_t lpos)
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

/**
 * Attempt to insert a pair of reads into the pos buffer's stash
 *
 * @param buf the pos buffer
 * @param ltid leftmost target id of pair
 * @param lpos leftmost position of pair
 * @param score score of pair
 * @param name read name of pair
 * @return name of read to kill, or NULL if it's new insert
 */
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

/**
 * Attempt to insert a pair of reads into the pos buffer's tree
 *
 * @param tree the tree for this rpos
 * @param ltid leftmost target id of pair
 * @param lpos leftmost position of pair
 * @param score score of pair
 * @param name read name of pair
 * @return name of read to kill, or NULL if it's new insert
 */
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
				pos_tree_t* insert = (pos_tree_t*)malloc(sizeof(pos_tree_t));
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
				pos_tree_t* insert = (pos_tree_t*)malloc(sizeof(pos_tree_t));
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
static char* pos_buffer_insert_inner(pos_buffer_t* buf, uint32_t ltid, uint32_t lpos, uint32_t rtid, uint32_t rpos, int score, const char* name)
{
	// If reads are in coordinate order it should not be possible for rpos to be less than base_pos
	assert(buf->base_pos <= rpos);
	
	// If read pair's right position is past the end of our buffer put it to one side
	if (rpos > (buf->base_pos + POS_BUFFER_LENGTH)) {
		// it's outside our buffer, stash it
		return pos_aux_insert(buf, ltid, lpos, rpos, score, name);
	} else {
		// it's within our buffer, see if there's a read pair at the
		pos_tree_t* tree = buf->right_most[pos_buff_rpos(buf, rpos)];
		if ( tree != NULL) {
			return pos_tree_insert(tree, ltid, lpos, score, name);
		} else {
			// create new tree
			pos_tree_t* insert = (pos_tree_t*)malloc(sizeof(pos_tree_t));
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

char* pos_buffer_insert(pos_buffer_t* buf, read_vector_t left, read_vector_t right, int score, const char* name, uint32_t curr_tid, uint32_t curr_pos)
{
	// If target id is not equal to our current one we need to start from scratch
	if (curr_tid != buf->tid) {
		pos_buffer_reset(buf);
		buf->tid = curr_tid;
	} else {
		pos_buffer_advance(buf, curr_pos);
	}
	if (read_vector_gt(left, right))
	{
		read_vector_t tmp = left;
		left = right;
		right = tmp;
	}
	return pos_buffer_insert_inner(buf, left.tid, left.pos, right.tid, right.pos, score, name);
}
