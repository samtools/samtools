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

#define GLF_TYPE_NORMAL 0
#define GLF_TYPE_INDEL  1
#define GLF_TYPE_END   15

typedef struct {
	unsigned char ref_base:4, type:4; /** "XACMGRSVTWYHKDBN"[ref_base] gives the reference base */
	unsigned char max_mapQ; /** maximum mapping quality */
	unsigned char lk[10];   /** log likelihood ratio, capped at 255 */
	unsigned min_lk:8, depth:24; /** minimum lk capped at 255, and the number of mapped reads */
	unsigned pos; /** this is ***ZERO-BASED*** coordinate */
} glf2_t;

typedef struct {
	int32_t l_text;
	uint8_t *text;
} glf_header_t;

#ifdef __cplusplus
extern "C" {
#endif

	glf_header_t *glf_header_init();
	glf_header_t *glf_header_read(glfFile fp);
	void glf_header_write(glfFile fp, const glf_header_t *h);
	void glf_header_destroy(glf_header_t *h);
	char *glf_ref_read(glfFile fp);
	void glf_ref_write(glfFile fp, const char *str);

#ifdef __cplusplus
}
#endif

#endif
