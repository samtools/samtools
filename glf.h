#ifndef GLF_H_
#define GLF_H_

typedef struct {
	unsigned char ref_base:4, dummy:4; /** "XACMGRSVTWYHKDBN"[ref_base] gives the reference base */
	unsigned char max_mapQ; /** maximum mapping quality */
	unsigned char lk[10];   /** log likelihood ratio, capped at 255 */
	unsigned min_lk:8, depth:24; /** minimum lk capped at 255, and the number of mapped reads */
} glf1_t;

#endif
