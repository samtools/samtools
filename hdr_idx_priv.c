#include "hdr_idx_priv.h"
#include <regex.h>
#include <htslib/khash.h>

__KHASH_IMPL(c2s, , kh_cstr_t, char*, 1, kh_str_hash_func, kh_str_hash_equal)

/*
 * This function takes a header and parses it to extract all the readgroups which it
 * then stores in library_index_t
 */
library_index_t* bam_library_index_init(const bam_hdr_t* hdr)
{
	library_index_t* retval = kh_init(c2s);

	// Compile a regex to extract every RG line with a LB and ID
	regex_t rg_regex;
	regmatch_t matches[7];
	int error;
	if ((error = regcomp( &rg_regex, "^@RG.*\t(ID|LB):([A-Za-z0-9_]*)\t(.*\t)?(ID|LB):([A-Za-z0-9_]*)(\t.*)?$", REG_EXTENDED|REG_NEWLINE )) != 0) {
		kh_destroy(c2s, retval);
		retval = NULL;
		goto bam_library_index_init_end;
	}
	char* rg_pointer = hdr->text;
	while (hdr->text+hdr->l_text > rg_pointer && regexec( &rg_regex, rg_pointer, 7, &matches[0], 0) == 0) {
		char* id = NULL;
		char* lb = NULL;
		const uint16_t id_16 = *((uint16_t*)"ID");
		uint16_t* label = (uint16_t*)&rg_pointer[matches[1].rm_so];
		if (*label == id_16) {
			id = strndup(rg_pointer + matches[2].rm_so, matches[2].rm_eo-matches[2].rm_so);
		} else {
			lb = strndup(rg_pointer + matches[2].rm_so, matches[2].rm_eo-matches[2].rm_so);
		}

		label = (uint16_t*)&rg_pointer[matches[4].rm_so];
		if (*label == id_16) {
			id = strndup(rg_pointer + matches[5].rm_so, matches[5].rm_eo-matches[5].rm_so);
		} else {
			lb = strndup(rg_pointer + matches[5].rm_so, matches[5].rm_eo-matches[5].rm_so);
		}
		// This should only occur in the case where a line has two LB or ID tags
		if (!id || !lb)
		{
			free(id); free(lb);
		}
		else {
			int ret;
			khint_t handle = kh_put(c2s, retval, id, &ret);
			kh_val(retval,handle) = lb;
		}

		// Move beginning of search pointer forward
		rg_pointer += matches[0].rm_eo + 1;
	}

bam_library_index_init_end:
	regfree(&rg_regex);
	return retval;
}

void bam_library_index_destroy(library_index_t* hdr)
{
	khint_t k;
	for (k = 0; k < kh_end(hdr); ++k) {
		if (kh_exist(hdr, k)) {
			free((void*)kh_key(hdr, k));
			free((void*)kh_value(hdr, k));
		}
	}

	kh_destroy(c2s,hdr);
	return;
}

// Takes a read and tries to get which library it is from
const char* bam_search_library_index(library_index_t* libs, const bam1_t* search)
{
	// See if read has an RG tag
	uint8_t *rg_h;
	if ((rg_h = bam_aux_get(search,"RG")) == NULL) return NULL;

	// Okay it does, do we know the library for it?
	khint_t handle = kh_get(c2s, libs, bam_aux2Z(rg_h));
	if (handle == kh_end(libs)) return NULL;

	return kh_val(libs, handle);
}
