/* The MIT License

   Copyright (c) 2021 Genome Research Ltd (GRL).

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

#ifndef SPLAYSORT_H
#define SPLAYSORT_H

#define SPLAYSORT_INIT(name, type_t, __sort_lt)                               \
    typedef struct splaynode_##name {                                         \
        type_t value;                                                         \
        struct splaynode_##name *left;                                        \
        struct splaynode_##name *right;                                       \
        struct splaynode_##name *parent;                                      \
    } splaynode_##name;                                                       \
                                                                              \
    void rotate_left_##name(splaynode_##name *node);                          \
    void rotate_right_##name(splaynode_##name *node);                         \
    int splay_sort_##name(size_t n, type_t array[] );                         \
    int splay_flatten_##name(splaynode_##name *node, type_t dest[], size_t n);\
    splaynode_##name *splay_tree_##name(splaynode_##name *node);              \
    splaynode_##name *splay_insert_##name(splaynode_##name *node,             \
                                          type_t value,                       \
                                          splaynode_##name *node_ptr);        \
                                                                              \
    void rotate_left_##name(splaynode_##name *node) {                         \
        splaynode_##name *parent = node->parent;                              \
        splaynode_##name *grandparent = parent->parent;                       \
        parent->right = node->left;                                           \
        if (node->left != NULL) {                                             \
            node->left->parent = parent;                                      \
        }                                                                     \
        node->left = parent;                                                  \
        parent->parent = node;                                                \
        node->parent = grandparent;                                           \
                                                                              \
        if (grandparent != NULL) {                                            \
            if (grandparent->left == parent) {                                \
                grandparent->left = node;                                     \
            } else {                                                          \
                grandparent->right = node;                                    \
            }                                                                 \
        }                                                                     \
    }                                                                         \
                                                                              \
    void rotate_right_##name(splaynode_##name *node) {                        \
        splaynode_##name *parent = node->parent;                              \
        splaynode_##name *grandparent = parent->parent;                       \
        parent->left = node->right;                                           \
                                                                              \
        if (node->right != NULL) {                                            \
            node->right->parent = parent;                                     \
        }                                                                     \
        node->right = parent;                                                 \
        parent->parent = node;                                                \
        node->parent = grandparent;                                           \
                                                                              \
        if (grandparent != NULL) {                                            \
            if (grandparent->left == parent) {                                \
                grandparent->left = node;                                     \
            } else {                                                          \
                grandparent->right = node;                                    \
            }                                                                 \
        }                                                                     \
    }                                                                         \
    int splay_sort_##name(size_t n, type_t array[] ) {                        \
        if (n < 1) {                                                          \
            return 0;                                                         \
        }                                                                     \
        int i;                                                                \
        splaynode_##name *node_pool = malloc(sizeof(splaynode_##name) * n);   \
        if (node_pool == NULL) return -1;                                     \
        splaynode_##name *head = node_pool;                                   \
        head->value = array[0];                                               \
        head->left = NULL; head->right = NULL; head->parent = NULL;           \
        for (i = 1; i < n; i++) {                                             \
            head = splay_insert_##name(head, array[i], node_pool + i );       \
        }                                                                     \
                                                                              \
        if (splay_flatten_##name(head, array, n) == -1) {                     \
            free(node_pool);                                                  \
            return -1;                                                        \
        }                                                                     \
        free(node_pool);                                                      \
        return 0;                                                             \
    }                                                                         \
                                                                              \
    int splay_flatten_##name(splaynode_##name *head, type_t *dest, size_t n) {\
        int sp = 0, i = 0;                                                    \
        splaynode_##name *current = head;                                     \
        splaynode_##name **stack = malloc(sizeof(current)*n);                 \
        if (stack == NULL) return -1;                                         \
                                                                              \
        do {                                                                  \
            while (current != NULL && sp < n) {                               \
                stack[sp++] = current;                                        \
                current = current->left;                                      \
            }                                                                 \
            if (sp != 0) {                                                    \
                sp--;                                                         \
                dest[i++] = stack[sp]->value;                                 \
                current = stack[sp]->right;                                   \
            }                                                                 \
        } while (!(current == NULL && sp == 0));                              \
                                                                              \
        free(stack);                                                          \
        return 0;                                                             \
    }                                                                         \
    splaynode_##name *splay_insert_##name(splaynode_##name *head,             \
                                          type_t value,                       \
                                          splaynode_##name *node_ptr) {       \
        splaynode_##name *parent = NULL;                                      \
        while (head != NULL) {                                                \
            parent = head;                                                    \
            if (__sort_lt(value, head->value)) {                              \
                head = head->left;                                            \
            } else {                                                          \
                head = head->right;                                           \
            }                                                                 \
        }                                                                     \
        splaynode_##name *new_node = node_ptr;                                \
        new_node->value = value;                                              \
        new_node->left = NULL;                                                \
        new_node->right = NULL;                                               \
        new_node->parent = parent;                                            \
        if (parent) {                                                         \
            if (__sort_lt(value, parent->value)) {                            \
                parent->left = new_node;                                      \
            } else {                                                          \
                parent->right = new_node;                                     \
            }                                                                 \
        }                                                                     \
        new_node = splay_tree_##name(new_node);                               \
        return new_node;                                                      \
    }                                                                         \
                                                                              \
    splaynode_##name *splay_tree_##name(splaynode_##name *node) {             \
        splaynode_##name *parent = node->parent;                              \
                                                                              \
        if (node->parent == NULL) {                                           \
            return node;                                                      \
        }                                                                     \
        if (node == parent->left) {                                           \
            if (parent->parent == NULL) {                                     \
                /* zig */                                                     \
                rotate_right_##name(node);                                    \
            } else if (parent->parent->left == parent) {                      \
                /* left zig zig */                                            \
                rotate_right_##name(node);                                    \
                rotate_right_##name(node);                                    \
            } else {                                                          \
                /* right left zig zag */                                      \
                rotate_right_##name(node);                                    \
                rotate_left_##name(node);                                     \
            }                                                                 \
        } else {                                                              \
            if (parent->parent == NULL) {                                     \
            /* zig */                                                         \
            rotate_left_##name(node);                                         \
            } else if (parent->parent->right == parent) {                     \
                /* right zig zig */                                           \
                rotate_left_##name(node);                                     \
                rotate_left_##name(node);                                     \
            } else  {                                                         \
                /* left right zig zag */                                      \
                rotate_left_##name(node);                                     \
                rotate_right_##name(node);                                    \
            }                                                                 \
        }                                                                     \
                                                                              \
        if (node->parent != NULL) {                                           \
            return splay_tree_##name(node);                                   \
        }                                                                     \
        return node;                                                          \
    }                                                                         \


#define splaysort(name, n, array) splay_sort_##name(n, array)

#endif
