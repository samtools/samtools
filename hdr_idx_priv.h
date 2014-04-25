#ifndef HDR_IDX_PRIV
#define HDR_IDX_PRIV

#include <htslib/khash.h>
#include <htslib/sam.h>

KHASH_DECLARE(c2s, kh_cstr_t, char*)

typedef khash_t(c2s) library_index_t;

#ifdef __cplusplus
extern "C" {
#endif

	library_index_t* bam_library_index_init(const bam_hdr_t* hdr);
	void bam_library_index_destroy(library_index_t* hdr);
	const char* bam_search_library_index(library_index_t* libs, const bam1_t* search);

#ifdef __cplusplus
}
#endif

#endif // ifndef HDR_IDX_PRIV
