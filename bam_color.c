#include <ctype.h>
#include "bam.h"

/*!
 @abstract     Get the color encoding the previous and current base
 @param b      pointer to an alignment
 @param i      The i-th position, 0-based
 @return       color

 @discussion   Returns 0 no color information is found.
 */
char bam_aux_getCSi(bam1_t *b, int i)
{
	uint8_t *c = bam_aux_get(b, "CS");
	char *cs = NULL;

	// return the base if the tag was not found
	if(0 == c) return 0;

	cs = bam_aux2Z(c);
	// adjust for strandedness and leading adaptor
	if(bam1_strand(b)) i = strlen(cs) - 1 - i;
	else i++;
	return cs[i];
}

/*!
 @abstract     Get the color quality of the color encoding the previous and current base
 @param b      pointer to an alignment
 @param i      The i-th position, 0-based
 @return       color quality

 @discussion   Returns 0 no color information is found.
 */
char bam_aux_getCQi(bam1_t *b, int i)
{
	uint8_t *c = bam_aux_get(b, "CQ");
	char *cq = NULL;
	
	// return the base if the tag was not found
	if(0 == c) return 0;

	cq = bam_aux2Z(c);
	// adjust for strandedness
	if(bam1_strand(b)) i = strlen(cq) - 1 - i;
	return cq[i];
}

char bam_aux_nt2int(char a)
{
	switch(toupper(a)) {
		case 'A':
			return 0;
			break;
		case 'C':
			return 1;
			break;
		case 'G':
			return 2;
			break;
		case 'T':
			return 3;
			break;
		default:
			return 4;
			break;
	}
}

char bam_aux_ntnt2cs(char a, char b)
{
	a = bam_aux_nt2int(a);
	b = bam_aux_nt2int(b);
	if(4 == a || 4 == b) return '4';
	return "0123"[(int)(a ^ b)];
}

/*!
 @abstract     Get the color error profile at the give position    
 @param b      pointer to an alignment
 @return       the original color if the color was an error, '-' (dash) otherwise

 @discussion   Returns 0 no color information is found.
 */
char bam_aux_getCEi(bam1_t *b, int i)
{
	int cs_i;
	uint8_t *c = bam_aux_get(b, "CS");
	char *cs = NULL;
	char prev_b, cur_b;
	char cur_color, cor_color;

	// return the base if the tag was not found
	if(0 == c) return 0;
	
	cs = bam_aux2Z(c);

	// adjust for strandedness and leading adaptor
	if(bam1_strand(b)) { //reverse strand
		cs_i = strlen(cs) - 1 - i;
		// get current color
		cur_color = cs[cs_i];
		// get previous base.  Note: must rc adaptor
		prev_b = (cs_i == 1) ? "TGCAN"[(int)bam_aux_nt2int(cs[0])] : bam_nt16_rev_table[bam1_seqi(bam1_seq(b), i+1)];
		// get current base
		cur_b = bam_nt16_rev_table[bam1_seqi(bam1_seq(b), i)]; 
	}
	else {
		cs_i=i+1;
		// get current color
		cur_color = cs[cs_i];
		// get previous base
		prev_b = (0 == i) ? cs[0] : bam_nt16_rev_table[bam1_seqi(bam1_seq(b), i-1)];
		// get current base
		cur_b = bam_nt16_rev_table[bam1_seqi(bam1_seq(b), i)];
	}

	// corrected color
	cor_color = bam_aux_ntnt2cs(prev_b, cur_b);

	if(cur_color == cor_color) { 
		return '-';
	}
	else {
		return cur_color;
	}
}
