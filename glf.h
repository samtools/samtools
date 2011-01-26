#ifndef GLF_H_
#define GLF_H_

typedef struct {
	unsigned char ref_base:4, dummy:4; /** "XACMGRSVTWYHKDBN"[ref_base] gives the reference base */
	unsigned char max_mapQ; /** maximum mapping quality */
	unsigned char lk[10];   /** log likelihood ratio, capped at 255 */
	unsigned min_lk:8, depth:24; /** minimum lk capped at 255, and the number of mapped reads */
} glf1_t;

#include <stdint.h>
#include "bgzf.h"
typedef BGZF *glfFile;

#define GLF3_RTYPE_END   0
#define GLF3_RTYPE_SUB   1
#define GLF3_RTYPE_INDEL 2

typedef struct {
	uint8_t ref_base:4, rtype:4; /** "XACMGRSVTWYHKDBN"[ref_base] gives the reference base */
	uint8_t rms_mapQ; /** RMS mapping quality */
	uint8_t lk[10];   /** log likelihood ratio, capped at 255 */
	uint32_t min_lk:8, depth:24; /** minimum lk capped at 255, and the number of mapped reads */
	int32_t offset; /** the first base in a chromosome has offset zero. */
	// for indel (lkHom1, lkHom2 and lkHet are the first three elements in lk[10])
	int16_t indel_len[2];
	int32_t max_len; // maximum indel len; will be modified by glf3_read1()
	char *indel_seq[2];
} glf3_t;

typedef struct {
	int32_t l_text;
	uint8_t *text;
} glf3_header_t;

#ifdef __cplusplus
extern "C" {
#endif

#define glf3_init1() ((glf3_t*)calloc(1, sizeof(glf3_t)))
#define glf3_destroy1(g3) do { free((g3)->indel_seq[0]); free((g3)->indel_seq[1]); free(g3); } while (0)

	glf3_header_t *glf3_header_init();
	glf3_header_t *glf3_header_read(glfFile fp);
	void glf3_header_write(glfFile fp, const glf3_header_t *h);
	void glf3_header_destroy(glf3_header_t *h);
	char *glf3_ref_read(glfFile fp, int *len);
	void glf3_ref_write(glfFile fp, const char *name, int len);
	int glf3_write1(glfFile fp, const glf3_t *g3);
	int glf3_read1(glfFile fp, glf3_t *g3);

#ifdef __cplusplus
}
#endif

#endif
