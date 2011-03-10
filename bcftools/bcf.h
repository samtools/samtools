/* The MIT License

   Copyright (c) 2010 Broad Institute

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

/* Contact: Heng Li <lh3@live.co.uk> */

#ifndef BCF_H
#define BCF_H

#include <stdint.h>
#include <zlib.h>

#ifndef BCF_LITE
#include "bgzf.h"
typedef BGZF *bcfFile;
#else
typedef gzFile bcfFile;
#define bgzf_open(fn, mode) gzopen(fn, mode)
#define bgzf_fdopen(fd, mode) gzdopen(fd, mode)
#define bgzf_close(fp) gzclose(fp)
#define bgzf_read(fp, buf, len) gzread(fp, buf, len)
#define bgzf_write(fp, buf, len)
#define bgzf_flush(fp)
#endif

/*
  A member in the structs below is said to "primary" if its content
  cannot be inferred from other members in any of structs below; a
  member is said to be "derived" if its content can be derived from
  other members. For example, bcf1_t::str is primary as this comes from
  the input data, while bcf1_t::info is derived as it can always be
  correctly set if we know bcf1_t::str. Derived members are for quick
  access to the content and must be synchronized with the primary data.
 */

typedef struct {
	uint32_t fmt; // format of the block, set by bcf_str2int(). 
	int len; // length of data for each individual
	void *data; // concatenated data
	// derived info: fmt, len (<-bcf1_t::fmt)
} bcf_ginfo_t;

typedef struct {
	int32_t tid, pos; // refID and 0-based position
	int32_t l_str, m_str; // length and the allocated size of ->str
	float qual; // SNP quality
	char *str; // concatenated string of variable length strings in VCF (from col.2 to col.7)
	char *ref, *alt, *flt, *info, *fmt; // they all point to ->str; no memory allocation
	int n_gi, m_gi; // number and the allocated size of geno fields
	bcf_ginfo_t *gi; // array of geno fields
	int n_alleles, n_smpl; // number of alleles and samples
	// derived info: ref, alt, flt, info, fmt (<-str), n_gi (<-fmt), n_alleles (<-alt), n_smpl (<-bcf_hdr_t::n_smpl)
} bcf1_t;

typedef struct {
	int32_t n_ref, n_smpl; // number of reference sequences and samples
	int32_t l_nm; // length of concatenated sequence names; 0 padded
	int32_t l_smpl; // length of concatenated sample names; 0 padded
	int32_t l_txt; // length of header text (lines started with ##)
	char *name, *sname, *txt; // concatenated sequence names, sample names and header text
	char **ns, **sns; // array of sequence and sample names; point to name and sname, respectively
	// derived info: n_ref (<-name), n_smpl (<-sname), ns (<-name), sns (<-sname)
} bcf_hdr_t;

typedef struct {
	int is_vcf; // if the file in operation is a VCF
	void *v; // auxillary data structure for VCF
	bcfFile fp; // file handler for BCF
} bcf_t;

struct __bcf_idx_t;
typedef struct __bcf_idx_t bcf_idx_t;

#ifdef __cplusplus
extern "C" {
#endif

	// open a BCF file; for BCF file only
	bcf_t *bcf_open(const char *fn, const char *mode);
	// close file
	int bcf_close(bcf_t *b);
	// read one record from BCF; return -1 on end-of-file, and <-1 for errors
	int bcf_read(bcf_t *bp, const bcf_hdr_t *h, bcf1_t *b);
	// call this function if b->str is changed
	int bcf_sync(bcf1_t *b);
	// write a BCF record
	int bcf_write(bcf_t *bp, const bcf_hdr_t *h, const bcf1_t *b);
	// read the BCF header; BCF only
	bcf_hdr_t *bcf_hdr_read(bcf_t *b);
	// write the BCF header
	int bcf_hdr_write(bcf_t *b, const bcf_hdr_t *h);
	// set bcf_hdr_t::ns and bcf_hdr_t::sns
	int bcf_hdr_sync(bcf_hdr_t *b);
	// destroy the header
	void bcf_hdr_destroy(bcf_hdr_t *h);
	// destroy a record
	int bcf_destroy(bcf1_t *b);
	// BCF->VCF conversion
	char *bcf_fmt(const bcf_hdr_t *h, bcf1_t *b);
	// append more info
	int bcf_append_info(bcf1_t *b, const char *info, int l);
	// copy
	int bcf_cpy(bcf1_t *r, const bcf1_t *b);

	// open a VCF or BCF file if "b" is set in "mode"
	bcf_t *vcf_open(const char *fn, const char *mode);
	// close a VCF/BCF file
	int vcf_close(bcf_t *bp);
	// read the VCF/BCF header
	bcf_hdr_t *vcf_hdr_read(bcf_t *bp);
	// read the sequence dictionary from a separate file; required for VCF->BCF conversion
	int vcf_dictread(bcf_t *bp, bcf_hdr_t *h, const char *fn);
	// read a VCF/BCF record; return -1 on end-of-file and <-1 for errors
	int vcf_read(bcf_t *bp, bcf_hdr_t *h, bcf1_t *b);
	// write the VCF header
	int vcf_hdr_write(bcf_t *bp, const bcf_hdr_t *h);
	// write a VCF record
	int vcf_write(bcf_t *bp, bcf_hdr_t *h, bcf1_t *b);

	// keep the first n alleles and discard the rest
	int bcf_shrink_alt(bcf1_t *b, int n);
	// convert GL to PL
	int bcf_gl2pl(bcf1_t *b);
	// if the site is an indel
	int bcf_is_indel(const bcf1_t *b);
	bcf_hdr_t *bcf_hdr_subsam(const bcf_hdr_t *h0, int n, char *const* samples, int *list);
	int bcf_subsam(int n_smpl, int *list, bcf1_t *b);
	// move GT to the first FORMAT field
	int bcf_fix_gt(bcf1_t *b);
	// update PL generated by old samtools
	int bcf_fix_pl(bcf1_t *b);

	// string hash table
	void *bcf_build_refhash(bcf_hdr_t *h);
	void bcf_str2id_destroy(void *_hash);
	void bcf_str2id_thorough_destroy(void *_hash);
	int bcf_str2id_add(void *_hash, const char *str);
	int bcf_str2id(void *_hash, const char *str);
	void *bcf_str2id_init();

	// indexing related functions
	int bcf_idx_build(const char *fn);
	uint64_t bcf_idx_query(const bcf_idx_t *idx, int tid, int beg);
	int bcf_parse_region(void *str2id, const char *str, int *tid, int *begin, int *end);
	bcf_idx_t *bcf_idx_load(const char *fn);
	void bcf_idx_destroy(bcf_idx_t *idx);

#ifdef __cplusplus
}
#endif

static inline uint32_t bcf_str2int(const char *str, int l)
{
	int i;
	uint32_t x = 0;
	for (i = 0; i < l && i < 4; ++i) {
		if (str[i] == 0) return x;
		x = x<<8 | str[i];
	}
	return x;
}

#endif
