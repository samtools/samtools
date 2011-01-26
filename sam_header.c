#include "sam_header.h"
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <stdarg.h>

#include "khash.h"
KHASH_MAP_INIT_STR(str, const char *)

struct _HeaderList
{
    struct _HeaderList *last;   // Hack: Used and maintained only by list_append_to_end. Maintained in the root node only.
    struct _HeaderList *next;
    void *data;
};
typedef struct _HeaderList list_t;
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

const char *o_hd_tags[] = {"SO","GO",NULL};
const char *r_hd_tags[] = {"VN",NULL};

const char *o_sq_tags[] = {"AS","M5","UR","SP",NULL};
const char *r_sq_tags[] = {"SN","LN",NULL};
const char *u_sq_tags[] = {"SN",NULL};

const char *o_rg_tags[] = {"LB","DS","PU","PI","CN","DT","PL",NULL};
const char *r_rg_tags[] = {"ID",NULL};
const char *u_rg_tags[] = {"ID",NULL};

const char *o_pg_tags[] = {"VN","CL",NULL};
const char *r_pg_tags[] = {"ID",NULL};

const char *types[]          = {"HD","SQ","RG","PG","CO",NULL};
const char **optional_tags[] = {o_hd_tags,o_sq_tags,o_rg_tags,o_pg_tags,NULL,NULL};
const char **required_tags[] = {r_hd_tags,r_sq_tags,r_rg_tags,r_pg_tags,NULL,NULL};
const char **unique_tags[]   = {NULL,     u_sq_tags,u_rg_tags,NULL,NULL,NULL};


static void debug(const char *format, ...)
{
    va_list ap;
    va_start(ap, format);
    vfprintf(stderr, format, ap);
    va_end(ap);
}

#if 0
// Replaced by list_append_to_end
static list_t *list_prepend(list_t *root, void *data)
{
    list_t *l = malloc(sizeof(list_t));
    l->next = root;
    l->data = data;
    return l;
}
#endif

// Relies on the root->last being correct. Do not use with the other list_*
//  routines unless they are fixed to modify root->last as well.
static list_t *list_append_to_end(list_t *root, void *data)
{
    list_t *l = malloc(sizeof(list_t));
    l->last = l;
    l->next = NULL;
    l->data = data;

    if ( !root )
        return l;

    root->last->next = l;
    root->last = l;
    return root;
}

static list_t *list_append(list_t *root, void *data)
{
    list_t *l = root;
    while (l && l->next)
        l = l->next;
    if ( l ) 
    {
        l->next = malloc(sizeof(list_t));
        l = l->next;
    }
    else
    {
        l = malloc(sizeof(list_t));
        root = l;
    }
    l->data = data;
    l->next = NULL;
    return root;
}

static void list_free(list_t *root)
{
    list_t *l = root;
    while (root)
    {
        l = root;
        root = root->next;
        free(l);
    }
}



// Look for a tag "XY" in a predefined const char *[] array.
static int tag_exists(const char *tag, const char **tags)
{
    int itag=0;
    if ( !tags ) return -1;
    while ( tags[itag] )
    {
        if ( tags[itag][0]==tag[0] && tags[itag][1]==tag[1] ) return itag; 
        itag++;
    }
    return -1;
}



// Mimics the behaviour of getline, except it returns pointer to the next chunk of the text
//  or NULL if everything has been read. The lineptr should be freed by the caller. The
//  newline character is stripped.
static const char *nextline(char **lineptr, size_t *n, const char *text)
{
    int len;
    const char *to = text;

    if ( !*to ) return NULL;

    while ( *to && *to!='\n' && *to!='\r' ) to++;
    len = to - text + 1;

    if ( *to )
    {
        // Advance the pointer for the next call
        if ( *to=='\n' ) to++;
        else if ( *to=='\r' && *(to+1)=='\n' ) to+=2;
    }
    if ( !len )
        return to;

    if ( !*lineptr ) 
    {
        *lineptr = malloc(len);
        *n = len;
    }
    else if ( *n<len ) 
    {
        *lineptr = realloc(*lineptr, len);
        *n = len;
    }
    if ( !*lineptr ) {
		debug("[nextline] Insufficient memory!\n");
		return 0;
	}

    memcpy(*lineptr,text,len);
    (*lineptr)[len-1] = 0;

    return to;
}

// name points to "XY", value_from points to the first character of the value string and
//  value_to points to the last character of the value string.
static HeaderTag *new_tag(const char *name, const char *value_from, const char *value_to)
{
    HeaderTag *tag = malloc(sizeof(HeaderTag));
    int len = value_to-value_from+1;

    tag->key[0] = name[0];
    tag->key[1] = name[1];
    tag->value = malloc(len+1);
    memcpy(tag->value,value_from,len+1);
    tag->value[len] = 0;
    return tag;
}

static HeaderTag *header_line_has_tag(HeaderLine *hline, const char *key)
{
    list_t *tags = hline->tags;
    while (tags)
    {
        HeaderTag *tag = tags->data;
        if ( tag->key[0]==key[0] && tag->key[1]==key[1] ) return tag;
        tags = tags->next;
    }
    return NULL;
}


// Return codes:
//   0 .. different types or unique tags differ or conflicting tags, cannot be merged
//   1 .. all tags identical -> no need to merge, drop one
//   2 .. the unique tags match and there are some conflicting tags (same tag, different value) -> error, cannot be merged nor duplicated
//   3 .. there are some missing complementary tags and no unique conflict -> can be merged into a single line
static int sam_header_compare_lines(HeaderLine *hline1, HeaderLine *hline2)
{
    HeaderTag *t1, *t2;

    if ( hline1->type[0]!=hline2->type[0] || hline1->type[1]!=hline2->type[1] )
        return 0;

    int itype = tag_exists(hline1->type,types);
    if ( itype==-1 ) {
		debug("[sam_header_compare_lines] Unknown type [%c%c]\n", hline1->type[0],hline1->type[1]);
		return -1; // FIXME (lh3): error; I do not know how this will be handled in Petr's code
	}

    if ( unique_tags[itype] )
    {
        t1 = header_line_has_tag(hline1,unique_tags[itype][0]);
        t2 = header_line_has_tag(hline2,unique_tags[itype][0]);
        if ( !t1 || !t2 ) // this should never happen, the unique tags are required
            return 2;

        if ( strcmp(t1->value,t2->value) )
            return 0;   // the unique tags differ, cannot be merged
    }
    if ( !required_tags[itype] && !optional_tags[itype] )
    {
        t1 = hline1->tags->data;
        t2 = hline2->tags->data;
        if ( !strcmp(t1->value,t2->value) ) return 1; // identical comments
        return 0;
    }

    int missing=0, itag=0;
    while ( required_tags[itype] && required_tags[itype][itag] )
    {
        t1 = header_line_has_tag(hline1,required_tags[itype][itag]);
        t2 = header_line_has_tag(hline2,required_tags[itype][itag]);
        if ( !t1 && !t2 )
            return 2;       // this should never happen
        else if ( !t1 || !t2 )
            missing = 1;    // there is some tag missing in one of the hlines
        else if ( strcmp(t1->value,t2->value) )
        {
            if ( unique_tags[itype] )
                return 2;   // the lines have a matching unique tag but have a conflicting tag
                    
            return 0;    // the lines contain conflicting tags, cannot be merged
        }
        itag++;
    }
    itag = 0;
    while ( optional_tags[itype] && optional_tags[itype][itag] )
    {
        t1 = header_line_has_tag(hline1,optional_tags[itype][itag]);
        t2 = header_line_has_tag(hline2,optional_tags[itype][itag]);
        if ( !t1 && !t2 )
        {
            itag++;
            continue;
        }
        if ( !t1 || !t2 )
            missing = 1;    // there is some tag missing in one of the hlines
        else if ( strcmp(t1->value,t2->value) )
        {
            if ( unique_tags[itype] )
                return 2;   // the lines have a matching unique tag but have a conflicting tag

            return 0;   // the lines contain conflicting tags, cannot be merged
        }
        itag++;
    }
    if ( missing ) return 3;    // there are some missing complementary tags with no conflicts, can be merged
    return 1;
}


static HeaderLine *sam_header_line_clone(const HeaderLine *hline)
{
    list_t *tags;
    HeaderLine *out = malloc(sizeof(HeaderLine));
    out->type[0] = hline->type[0];
    out->type[1] = hline->type[1];
    out->tags = NULL;

    tags = hline->tags;
    while (tags)
    {
        HeaderTag *old = tags->data;

        HeaderTag *new = malloc(sizeof(HeaderTag));
        new->key[0] = old->key[0];
        new->key[1] = old->key[1];
        new->value  = strdup(old->value);
        out->tags = list_append(out->tags, new);

        tags = tags->next;
    }
    return out;
}

static int sam_header_line_merge_with(HeaderLine *out_hline, const HeaderLine *tmpl_hline)
{
    list_t *tmpl_tags;

    if ( out_hline->type[0]!=tmpl_hline->type[0] || out_hline->type[1]!=tmpl_hline->type[1] )
        return 0;
    
    tmpl_tags = tmpl_hline->tags;
    while (tmpl_tags)
    {
        HeaderTag *tmpl_tag = tmpl_tags->data;
        HeaderTag *out_tag  = header_line_has_tag(out_hline, tmpl_tag->key);
        if ( !out_tag )
        {
            HeaderTag *tag = malloc(sizeof(HeaderTag));
            tag->key[0] = tmpl_tag->key[0];
            tag->key[1] = tmpl_tag->key[1];
            tag->value  = strdup(tmpl_tag->value);
            out_hline->tags = list_append(out_hline->tags,tag);
        }
        tmpl_tags = tmpl_tags->next;
    }
    return 1;
}


static HeaderLine *sam_header_line_parse(const char *headerLine)
{
    HeaderLine *hline;
    HeaderTag *tag;
    const char *from, *to;
    from = headerLine;

    if ( *from != '@' ) {
		debug("[sam_header_line_parse] expected '@', got [%s]\n", headerLine);
		return 0;
	}
    to = ++from;

    while (*to && *to!='\t') to++;
    if ( to-from != 2 ) {
		debug("[sam_header_line_parse] expected '@XY', got [%s]\nHint: The header tags must be tab-separated.\n", headerLine);
		return 0;
	}
    
    hline = malloc(sizeof(HeaderLine));
    hline->type[0] = from[0];
    hline->type[1] = from[1];
    hline->tags = NULL;

    int itype = tag_exists(hline->type, types);
    
    from = to;
    while (*to && *to=='\t') to++;
    if ( to-from != 1 ) {
        debug("[sam_header_line_parse] multiple tabs on line [%s] (%d)\n", headerLine,(int)(to-from));
		return 0;
	}
    from = to;
    while (*from)
    {
        while (*to && *to!='\t') to++;

        if ( !required_tags[itype] && !optional_tags[itype] )
        {
            // CO is a special case, it can contain anything, including tabs
            if ( *to ) { to++; continue; }
            tag = new_tag("  ",from,to-1);
        }
        else
            tag = new_tag(from,from+3,to-1);

        if ( header_line_has_tag(hline,tag->key) ) 
                debug("The tag '%c%c' present (at least) twice on line [%s]\n", tag->key[0],tag->key[1], headerLine);
        hline->tags = list_append(hline->tags, tag);

        from = to;
        while (*to && *to=='\t') to++;
        if ( *to && to-from != 1 ) {
			debug("[sam_header_line_parse] multiple tabs on line [%s] (%d)\n", headerLine,(int)(to-from));
			return 0;
		}

        from = to;
    }
    return hline;
}


// Must be of an existing type, all tags must be recognised and all required tags must be present
static int sam_header_line_validate(HeaderLine *hline)
{
    list_t *tags;
    HeaderTag *tag;
    int itype, itag;
    
    // Is the type correct?
    itype = tag_exists(hline->type, types);
    if ( itype==-1 ) 
    {
        debug("The type [%c%c] not recognised.\n", hline->type[0],hline->type[1]);
        return 0;
    }

    // Has all required tags?
    itag = 0;
    while ( required_tags[itype] && required_tags[itype][itag] )
    {
        if ( !header_line_has_tag(hline,required_tags[itype][itag]) )
        {
            debug("The tag [%c%c] required for [%c%c] not present.\n", required_tags[itype][itag][0],required_tags[itype][itag][1],
                hline->type[0],hline->type[1]);
            return 0;
        }
        itag++;
    }

    // Are all tags recognised?
    tags = hline->tags;
    while ( tags )
    {
        tag = tags->data;
        if ( !tag_exists(tag->key,required_tags[itype]) && !tag_exists(tag->key,optional_tags[itype]) )
        {
            debug("Unknown tag [%c%c] for [%c%c].\n", tag->key[0],tag->key[1], hline->type[0],hline->type[1]);
            return 0;
        }
        tags = tags->next;
    }

    return 1;
}


static void print_header_line(FILE *fp, HeaderLine *hline)
{
    list_t *tags = hline->tags;
    HeaderTag *tag;

    fprintf(fp, "@%c%c", hline->type[0],hline->type[1]);
    while (tags)
    {
        tag = tags->data;

        fprintf(fp, "\t");
        if ( tag->key[0]!=' ' || tag->key[1]!=' ' )
            fprintf(fp, "%c%c:", tag->key[0],tag->key[1]);
        fprintf(fp, "%s", tag->value);

        tags = tags->next;
    }
    fprintf(fp,"\n");
}


static void sam_header_line_free(HeaderLine *hline)
{
    list_t *tags = hline->tags;
    while (tags)
    {
        HeaderTag *tag = tags->data;
        free(tag->value);
        free(tag);
        tags = tags->next;
    }
    list_free(hline->tags);
    free(hline);
}

void sam_header_free(void *_header)
{
	HeaderDict *header = (HeaderDict*)_header;
    list_t *hlines = header;
    while (hlines)
    {
        sam_header_line_free(hlines->data);
        hlines = hlines->next;
    }
    list_free(header);
}

HeaderDict *sam_header_clone(const HeaderDict *dict)
{
    HeaderDict *out = NULL;
    while (dict)
    {
        HeaderLine *hline = dict->data;
        out = list_append(out, sam_header_line_clone(hline));
        dict = dict->next;
    }
    return out;
}

// Returns a newly allocated string
char *sam_header_write(const void *_header)
{
	const HeaderDict *header = (const HeaderDict*)_header;
    char *out = NULL;
    int len=0, nout=0;
    const list_t *hlines;

    // Calculate the length of the string to allocate
    hlines = header;
    while (hlines)
    {
        len += 4;   // @XY and \n

        HeaderLine *hline = hlines->data;
        list_t *tags = hline->tags;
        while (tags)
        {
            HeaderTag *tag = tags->data;
            len += strlen(tag->value) + 1;                  // \t
            if ( tag->key[0]!=' ' || tag->key[1]!=' ' )
                len += strlen(tag->value) + 3;              // XY:
            tags = tags->next;
        }
        hlines = hlines->next;
    }

    nout = 0;
    out  = malloc(len+1);
    hlines = header;
    while (hlines)
    {
        HeaderLine *hline = hlines->data;

        nout += sprintf(out+nout,"@%c%c",hline->type[0],hline->type[1]);

        list_t *tags = hline->tags;
        while (tags)
        {
            HeaderTag *tag = tags->data;
            nout += sprintf(out+nout,"\t");
            if ( tag->key[0]!=' ' || tag->key[1]!=' ' )
                nout += sprintf(out+nout,"%c%c:", tag->key[0],tag->key[1]);
            nout += sprintf(out+nout,"%s", tag->value);
            tags = tags->next;
        }
        hlines = hlines->next;
        nout += sprintf(out+nout,"\n");
    }
    out[len] = 0;
    return out;
}

void *sam_header_parse2(const char *headerText)
{
    list_t *hlines = NULL;
    HeaderLine *hline;
    const char *text;
    char *buf=NULL;
    size_t nbuf = 0;

    if ( !headerText )
		return 0;

    text = headerText;
    while ( (text=nextline(&buf, &nbuf, text)) )
    {
        hline = sam_header_line_parse(buf);
        if ( hline && sam_header_line_validate(hline) )
            // With too many (~250,000) reference sequences the header parsing was too slow with list_append.
            hlines = list_append_to_end(hlines, hline);
        else
        {
			if (hline) sam_header_line_free(hline);
			sam_header_free(hlines);
            if ( buf ) free(buf);
            return NULL;
        }
    }
    if ( buf ) free(buf);

    return hlines;
}

void *sam_header2tbl(const void *_dict, char type[2], char key_tag[2], char value_tag[2])
{
	const HeaderDict *dict = (const HeaderDict*)_dict;
    const list_t *l   = dict;
    khash_t(str) *tbl = kh_init(str);
    khiter_t k;
    int ret;

	if (_dict == 0) return tbl; // return an empty (not null) hash table
    while (l)
    {
        HeaderLine *hline = l->data;
        if ( hline->type[0]!=type[0] || hline->type[1]!=type[1] ) 
        {
            l = l->next;
            continue;
        }
        
        HeaderTag *key, *value;
        key   = header_line_has_tag(hline,key_tag);
        value = header_line_has_tag(hline,value_tag); 
        if ( !key || !value )
        {
            l = l->next;
            continue;
        }
        
        k = kh_get(str, tbl, key->value);
        if ( k != kh_end(tbl) )
            debug("[sam_header_lookup_table] They key %s not unique.\n", key->value);
        k = kh_put(str, tbl, key->value, &ret);
        kh_value(tbl, k) = value->value;

        l = l->next;
    }
    return tbl;
}

char **sam_header2list(const void *_dict, char type[2], char key_tag[2], int *_n)
{
	const HeaderDict *dict = (const HeaderDict*)_dict;
    const list_t *l   = dict;
    int max, n;
	char **ret;

	ret = 0; *_n = max = n = 0;
    while (l)
    {
        HeaderLine *hline = l->data;
        if ( hline->type[0]!=type[0] || hline->type[1]!=type[1] ) 
        {
            l = l->next;
            continue;
        }
        
        HeaderTag *key;
        key   = header_line_has_tag(hline,key_tag);
        if ( !key )
        {
            l = l->next;
            continue;
        }

		if (n == max) {
			max = max? max<<1 : 4;
			ret = realloc(ret, max * sizeof(void*));
		}
		ret[n++] = key->value;

        l = l->next;
    }
	*_n = n;
    return ret;
}

const char *sam_tbl_get(void *h, const char *key)
{
	khash_t(str) *tbl = (khash_t(str)*)h;
	khint_t k;
	k = kh_get(str, tbl, key);
	return k == kh_end(tbl)? 0 : kh_val(tbl, k);
}

int sam_tbl_size(void *h)
{
	khash_t(str) *tbl = (khash_t(str)*)h;
	return h? kh_size(tbl) : 0;
}

void sam_tbl_destroy(void *h)
{
	khash_t(str) *tbl = (khash_t(str)*)h;
	kh_destroy(str, tbl);
}

void *sam_header_merge(int n, const void **_dicts)
{
	const HeaderDict **dicts = (const HeaderDict**)_dicts;
    HeaderDict *out_dict;
    int idict, status;

    if ( n<2 ) return NULL;

    out_dict = sam_header_clone(dicts[0]);

    for (idict=1; idict<n; idict++)
    {
        const list_t *tmpl_hlines = dicts[idict];

        while ( tmpl_hlines )
        {
            list_t *out_hlines = out_dict;
            int inserted = 0;
            while ( out_hlines )
            {
                status = sam_header_compare_lines(tmpl_hlines->data, out_hlines->data);
                if ( status==0 )
                {
                    out_hlines = out_hlines->next;
                    continue;
                }
                
                if ( status==2 ) 
                {
                    print_header_line(stderr,tmpl_hlines->data);
                    print_header_line(stderr,out_hlines->data);
                    debug("Conflicting lines, cannot merge the headers.\n");
					return 0;
                }
                if ( status==3 )
                    sam_header_line_merge_with(out_hlines->data, tmpl_hlines->data);

                inserted = 1;
                break;
            }
            if ( !inserted )
                out_dict = list_append(out_dict, sam_header_line_clone(tmpl_hlines->data));

            tmpl_hlines = tmpl_hlines->next;
        }
    }

    return out_dict;
}


