#ifndef __SAM_HEADER_H__
#define __SAM_HEADER_H__

#include "khash.h"
KHASH_MAP_INIT_STR(str,const char *)

// HeaderDict is a list_t of header lines. Each HeaderLine holds a list of tags.
struct _list_t
{
    struct _list_t *next;
    void *data;
};
typedef struct _list_t list_t;
typedef list_t HeaderDict;

typedef struct
{
    char key[2];
    char *value;
}
HeaderTag;

typedef struct
{
    char type[2];
    list_t *tags;
}
HeaderLine;


void debug(const char *format, ...);
void error(const char *format, ...);

HeaderDict *sam_header_parse(const char *headerText);
HeaderDict *sam_header_merge(int n, const HeaderDict **dicts);
void sam_header_free(HeaderDict *header);
char *sam_header_write(const HeaderDict *headerDict);   // returns a newly allocated string

khash_t(str) *sam_header_lookup_table(const HeaderDict *dict, char type[2], char key_tag[2], char value_tag[2]);

list_t *list_append(list_t *root, void *data);
void list_free(list_t *root);

//char *sam_header_get(const HeaderDict *d, char type[2], int i, char tag[2]);  
//int sam_header_ins(HeaderDict *d, char tp[2], int i, char tg[2], const char *s);  
//int sam_header_del(HeaderDict *d, char type[2], int i, char tag[2]);  

#endif
