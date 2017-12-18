#ifndef HTS_GUILE_H
#define HTS_GUILE_H

#include "config.h"

#ifdef HAVE_GUILE
#include "htslib/sam.h"
#include <libguile.h>
#include <libguile/modules.h>

typedef struct hts_guile_context_t {
	const bam_hdr_t *header;
	bam1_t *b;
	} HtsGuileCtx,*HtsGuileCtxPtr;

void* hts_guile_init();



#endif

#endif
