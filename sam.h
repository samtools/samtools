/*  sam.h -- format-neutral SAM/BAM API.

    Copyright (C) 2009, 2013-2015 Genome Research Ltd.

    Author: Heng Li <lh3@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#ifndef BAM_SAM_H
#define BAM_SAM_H

#include "htslib/sam.h"
#include "bam.h"

/*!
  @header

  This file provides higher level of I/O routines and unifies the APIs
  for SAM and BAM formats. These APIs are more convenient and
  recommended.

  @copyright Genome Research Ltd.
 */

/*! @typedef
  @abstract SAM/BAM file handler
  @field  type    type of the handler; bit 1 for BAM, 2 for reading and bit 3-4 for flag format
  @field  bam   BAM file handler; valid if (type&1) == 1
  @field  tamr  SAM file handler for reading; valid if type == 2
  @field  tamw  SAM file handler for writing; valid if type == 0
  @field  header  header struct
 */
typedef struct {
    samFile *file;
    struct { BGZF *bam; } x;  // Hack so that fp->x.bam still works
    bam_hdr_t *header;
    unsigned short is_write:1;
} samfile_t;

#ifdef __cplusplus
extern "C" {
#endif

    /*!
      @abstract     Open a SAM/BAM file

      @param fn SAM/BAM file name; "-" is recognized as stdin (for
      reading) or stdout (for writing).

      @param mode open mode /[rw](b?)(u?)(h?)([xX]?)/: 'r' for reading,
      'w' for writing, 'b' for BAM I/O, 'u' for uncompressed BAM output,
      'h' for outputing header in SAM, 'x' for HEX flag and 'X' for
      string flag. If 'b' present, it must immediately follow 'r' or
      'w'. Valid modes are "r", "w", "wh", "wx", "whx", "wX", "whX",
      "rb", "wb" and "wbu" exclusively.

      @param aux auxiliary data; if mode[0]=='w', aux points to
      bam_header_t; if strcmp(mode, "rb")!=0 and @SQ header lines in SAM
      are absent, aux points the file name of the list of the reference;
      aux is not used otherwise. If @SQ header lines are present in SAM,
      aux is not used, either.

      @return       SAM/BAM file handler
     */
    samfile_t *samopen(const char *fn, const char *mode, const void *aux);

    /*!
      @abstract     Close a SAM/BAM handler
      @param  fp    file handler to be closed
     */
    void samclose(samfile_t *fp);

    /*!
      @abstract     Read one alignment
      @param  fp    file handler
      @param  b     alignment
      @return       bytes read
     */
    static inline int samread(samfile_t *fp, bam1_t *b) { return sam_read1(fp->file, fp->header, b); }

    /*!
      @abstract     Write one alignment
      @param  fp    file handler
      @param  b     alignment
      @return       bytes written
     */
    static inline int samwrite(samfile_t *fp, const bam1_t *b) { return sam_write1(fp->file, fp->header, b); }

    /*!
      @abstract     Load BAM/CRAM index for use with samfetch()
      @param  fp    file handler
      @param  fn    name of the BAM or CRAM file (NOT the index file)
      @return       pointer to the index structure
     */
    static inline bam_index_t *samtools_sam_index_load(samfile_t *fp, const char *fn) { return sam_index_load(fp->file, fn); }
    #undef sam_index_load
    #define sam_index_load(fp,fn) (samtools_sam_index_load((fp), (fn)))

    /*!
      @abstract Retrieve the alignments overlapping the specified region.
      @discussion A user defined function will be called for each
      retrieved alignment ordered by its start position.
      @param  fp    file handler
      @param  idx   index returned by sam_index_load()
      @param  tid   chromosome ID as is defined in the header
      @param  beg   start coordinate, 0-based
      @param  end   end coordinate, 0-based
      @param  data  user provided data (will be transferred to func)
      @param  func  user defined function
     */
    int samfetch(samfile_t *fp, const bam_index_t *idx, int tid, int beg, int end, void *data, bam_fetch_f func);

    /*!
      @abstract     Get the pileup for a whole alignment file
      @param  fp    file handler
      @param  mask  mask transferred to bam_plbuf_set_mask()
      @param  func  user defined function called in the pileup process
      #param  data  user provided data for func()
     */
    int sampileup(samfile_t *fp, int mask, bam_pileup_f func, void *data);

    char *samfaipath(const char *fn_ref);
    int samthreads(samfile_t *fp, int n_threads, int n_sub_blks);

#ifdef __cplusplus
}
#endif

#endif
