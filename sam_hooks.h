#ifndef SAMTOOLS_HOOK_H
#define SAMTOOLS_HOOK_H

struct sam_hook_t {
  int (*initialize)(struct sam_hook_t* ,const bam_hdr_t *);
  int (*accept)(struct sam_hook_t* ,const bam_hdr_t *, bam1_t *);
  int (*dispose)(struct sam_hook_t* );
  };

struct sam_hook_t* hook_load_by_name(const char* name);

#endif

