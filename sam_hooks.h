#ifndef SAMTOOLS_HOOK_H
#define SAMTOOLS_HOOK_H

#include "htslib/sam.h"


typedef struct sam_hook_t {
  int (*initialize)(struct sam_hook_t* ,const bam_hdr_t *);
  int (*accept)(struct sam_hook_t* ,const bam_hdr_t *, bam1_t *);
  int (*dispose)(struct sam_hook_t* );
  char* (*get_name)(struct sam_hook_t* );
  
  /** next hoot in chained list */
  struct sam_hook_t* next;
  /* reserved for the hook */
  void* reserved;
  }SamHook, *SamHookPtr;

SamHookPtr hook_load_by_name(SamHookPtr root,const char* name);
void hook_dispose_all(SamHookPtr ptr);
int hook_accept_all(SamHookPtr hook,const bam_hdr_t *, bam1_t*);
#endif

