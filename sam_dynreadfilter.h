/*
    Author: Pierre Lindenbaum PhD @yokofakun
            Institut du Thorax - Nantes - France

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notices and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */
#ifndef SAMTOOLS_HOOK_H
#define SAMTOOLS_HOOK_H

#include "htslib/sam.h"

/** structure holding callbacks for dynamic read filtering using dynamic libary */
typedef struct sam_dyn_read_filter_t {
  /** initialize the filter */
  int (*initialize)(struct sam_dyn_read_filter_t* ,const bam_hdr_t *);
  /** shall we accept the read ? returns 0 if read should be discarded */
  int (*accept)(struct sam_dyn_read_filter_t* ,const bam_hdr_t *, bam1_t *);
  /** dispose internal data for this filter. returns 0 on failure */
  int (*dispose)(struct sam_dyn_read_filter_t* );
  /** return a name for this filter */
  char* (*get_name)(struct sam_dyn_read_filter_t* );
  /** next filter in chained list */
  struct sam_dyn_read_filter_t* next;
  /* user's data */
  void* data;
  }SamDynReadFilter, *SamDynReadFilterPtr;

/** load a dynamic filter by name */
SamDynReadFilterPtr dynreadfilter_load_by_name(const char* name);
/** dispose the filter and all its children in the linked list */
void  dynreadfilter_dispose_all(SamDynReadFilterPtr ptr);
/** return true if all filters in the linked list accept the read */
int  dynreadfilter_accept_all(SamDynReadFilterPtr hook,const bam_hdr_t *, bam1_t*);
/** appends a filter in the linked list, returns root, or child if root is null */
SamDynReadFilterPtr dynreadfilter_append(SamDynReadFilterPtr root,SamDynReadFilterPtr child);



#define DYNREADFILTER_DESC "  -k      <name> plugin name\n"
#endif

