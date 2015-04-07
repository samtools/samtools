/*  sam_header.h -- basic SAM/BAM header API.

    Copyright (C) 2009, 2012, 2013 Genome Research Ltd.

    Author: Petr Danecek <pd3@sanger.ac.uk>

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

#ifndef __SAM_HEADER_H__
#define __SAM_HEADER_H__

#ifdef __cplusplus
extern "C" {
#endif

    void *sam_header_parse2(const char *headerText);
    void *sam_header_merge(int n, const void **dicts);
    void sam_header_free(void *header);
    char *sam_header_write(const void *headerDict);   // returns a newly allocated string

    /*
        // Usage example
        const char *key, *val;
        void *iter = sam_header_parse2(bam->header->text);
        while ( iter = sam_header_key_val(iter, "RG","ID","SM" &key,&val) ) printf("%s\t%s\n", key,val);
    */
    void *sam_header2key_val(void *iter, const char type[2], const char key_tag[2], const char value_tag[2], const char **key, const char **value);
    char **sam_header2list(const void *_dict, char type[2], char key_tag[2], int *_n);

    /*
        // Usage example
        int i, j, n;
        const char *tags[] = {"SN","LN","UR","M5",NULL};
        void *dict = sam_header_parse2(bam->header->text);
        char **tbl = sam_header2tbl_n(h->dict, "SQ", tags, &n);
        for (i=0; i<n; i++)
        {
            for (j=0; j<4; j++)
                if ( tbl[4*i+j] ) printf("\t%s", tbl[4*i+j]);
                else printf("-");
            printf("\n");
        }
        if (tbl) free(tbl);
     */
    char **sam_header2tbl_n(const void *dict, const char type[2], const char *tags[], int *n);

    void *sam_header2tbl(const void *dict, char type[2], char key_tag[2], char value_tag[2]);
    const char *sam_tbl_get(void *h, const char *key);
    int sam_tbl_size(void *h);
    void sam_tbl_destroy(void *h);

#ifdef __cplusplus
}
#endif

#endif
