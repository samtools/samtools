#include "sam_header.h"
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <stdarg.h>

const char *o_hd_tags[] = {"SO","GO",NULL};
const char *r_hd_tags[] = {"VN",NULL};
const char *o_sq_tags[] = {"AS","M5","UR","SP",NULL};
const char *r_sq_tags[] = {"SN","LN",NULL};
const char *o_rg_tags[] = {"LB","DS","PU","PI","CN","DT","PL",NULL};
const char *r_rg_tags[] = {"ID","SM",NULL};
const char *o_pg_tags[] = {"VN","CL",NULL};
const char *r_pg_tags[] = {"ID",NULL};
const char *types[]          = {"HD","SQ","RG","PG","CO",NULL};
const char **optional_tags[] = {o_hd_tags,o_sq_tags,o_rg_tags,o_pg_tags,NULL,NULL};
const char **required_tags[] = {r_hd_tags,r_sq_tags,r_rg_tags,r_pg_tags,NULL,NULL};

void debug(const char *format, ...)
{
    va_list ap;
    va_start(ap, format);
    vfprintf(stderr, format, ap);
    va_end(ap);
}

void error(const char *format, ...)
{
    va_list ap;
    va_start(ap, format);
    vfprintf(stderr, format, ap);
    va_end(ap);
    exit(-1);
}

list_t *list_append(list_t *root, void *data)
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

void list_free(list_t *root)
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
int tag_exists(const char *tag, const char **tags)
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
const char *nextline(char **lineptr, size_t *n, const char *text)
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
    if ( !*lineptr )
            error("FIXME\n");

    memcpy(*lineptr,text,len);
    (*lineptr)[len-1] = 0;

    return to;
}

// name points to "XY", value_from points to the first character of the value string and
//  value_to points to the last character of the value string.
HeaderTag *new_tag(const char *name, const char *value_from, const char *value_to)
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

HeaderTag *header_line_has_tag(HeaderLine *hline, const char *key)
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

#if 0
// Is there a HeaderLine with all required fields identical to those given in the hline?
HeaderLine *sam_header_has_line(HeaderDict *dict, HeaderLine *hline)
{
    HeaderLine *found=NULL;

    while (dict)
    {
        HeaderLine *dline = dict->data;

        if ( hline->type[0]!=dline->type[0] || hline->type[1]!=dline->type[1] )
        {
            dict = dict->next;
            continue;
        }

        int itype = tag_exists(hline->type,types);
        if ( itype==-1 ) error("[sam_header_has_line] Unknown type [%c%c]\n", hline->type[0],hline->type[1]);

        int ireq=0, differ=0;
        while ( required_tags[itype] && required_tags[itype][ireq] )
        {
            HeaderTag *t1, *t2;
            t1 = header_line_has_tag(hline,required_tags[itype][ireq]);
            t2 = header_line_has_tag(dline,required_tags[itype][ireq]);
            if ( !t1 || !t2 ) error("[sam_header_has_line] Missing a required tag [%c%c]\n",
                required_tags[itype][ireq][0],required_tags[itype][ireq][1]);
            if ( strcmp(t1->value,t2->value) )
            ireq++;
        }
        dict = dict->next; 
    }
    return found;
}
#endif

HeaderLine *sam_header_line_parse(const char *headerLine)
{
    HeaderLine *hline;
    HeaderTag *tag;
    const char *from, *to;
    from = headerLine;

    if ( *from != '@' ) error("[sam_header_line_parse] expected '@', got [%s]\n", headerLine);
    to = ++from;

    while (*to && *to!='\t') to++;
    if ( to-from != 2 ) error("[sam_header_line_parse] expected '@XY', got [%s]\n", headerLine);
    
    hline = malloc(sizeof(HeaderLine));
    hline->type[0] = from[0];
    hline->type[1] = from[1];
    hline->tags = NULL;

    int itype = tag_exists(hline->type, types);
    
    from = to;
    while (*to && *to=='\t') to++;
    if ( to-from != 1 ) 
        error("[sam_header_line_parse] multiple tabs on line [%s] (%d)\n", headerLine,(int)(to-from));
    from = to;
    while (*from)
    {
        while (*to && *to!='\t') to++;

        if ( !required_tags[itype] && !optional_tags[itype] )
            tag = new_tag("  ",from,to-1);
        else
            tag = new_tag(from,from+3,to-1);

        if ( header_line_has_tag(hline,tag->key) ) 
                debug("The tag '%c%c' present (at least) twice on line [%s]\n", tag->key[0],tag->key[1], headerLine);
        hline->tags = list_append(hline->tags, tag);

        from = to;
        while (*to && *to=='\t') to++;
        if ( *to && to-from != 1 ) 
                error("[sam_header_line_parse] multiple tabs on line [%s] (%d)\n", headerLine,(int)(to-from));

        from = to;
    }
    return hline;
}


// Must be of an existing type, all tags must be recognised and all required tags must be present
int sam_header_line_validate(HeaderLine *hline)
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

void print_header_line(HeaderLine *hline)
{
    list_t *tags = hline->tags;
    HeaderTag *tag;

    printf("@%c%c", hline->type[0],hline->type[1]);
    while (tags)
    {
        tag = tags->data;
        printf("\t%c%c:%s", tag->key[0],tag->key[1],tag->value);
        tags = tags->next;
    }
    printf("\n");
}


void sam_header_free(HeaderDict *header)
{
    list_t *hlines = header;
    while (hlines)
    {
        HeaderLine *hline = hlines->data;
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
        hlines = hlines->next;
    }
    list_free(header);
}

// Returns a newly allocated string
char *sam_header_write(const HeaderDict *header)
{
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

HeaderDict *sam_header_parse(const char *headerText)
{
    list_t *hlines = NULL;
    HeaderLine *hline;
    const char *text;
    char *buf=NULL;
    size_t nbuf = 0;

    if ( !headerText )
        error("FIXME");

    text = headerText;
    while ( (text=nextline(&buf, &nbuf, text)) )
    {
        hline = sam_header_line_parse(buf);
        if ( sam_header_line_validate(hline) )
            hlines = list_append(hlines, hline);
        else
        {
            sam_header_free(hlines);
            return NULL;
        }
    }
    if ( buf ) free(buf);

    return hlines;
}

khash_t(str) *sam_header_lookup_table(const HeaderDict *dict, char type[2], char key_tag[2], char value_tag[2])
{
    const list_t *l   = dict;
    khash_t(str) *tbl = kh_init(str);
    khiter_t k;
    int ret;

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


#if 0
TODO
HeaderDict *sam_header_merge(int n, const HeaderDict **dicts)
{
    HeaderDict *out=NULL;
    int idict;

    for (idict=0; idict<n; idict++)
    {
        list_t *hlines = dicts[idict];
        while (hlines)
        {
            HeaderLine *hline = sam_header_has_line(out, hlines->data);
            sam_header_line_merge(hline,hlines->data);
            hlines = hlines->next;
        }
    }
}
#endif

