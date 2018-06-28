/*  stats.c -- This is the former bamcheck integrated into samtools/htslib.

    Copyright (C) 2012-2015 Genome Research Ltd.

    Author: Petr Danecek <pd3@sanger.ac.uk>
    Author: Sam Nicholls <sam@samnicholls.net>

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

/*  Assumptions, approximations and other issues:
        - GC-depth graph does not split reads, the starting position determines which bin is incremented.
            There are small overlaps between bins (max readlen-1). However, the bins are big (20k).
        - coverage distribution ignores softclips and deletions
        - some stats require sorted BAMs
        - GC content graph can have an untidy, step-like pattern when BAM contains multiple read lengths.
        - 'bases mapped' (stats->nbases_mapped) is calculated from read lengths given by BAM (core.l_qseq)
        - With the -t option, the whole reads are used. Except for the number of mapped bases (cigar)
            counts, no splicing is done, no indels or soft clips are considered, even small overlap is
            good enough to include the read in the stats.
        - GC content of reads not calculated for "=" sequences

*/

#include <config.h>

#include <unistd.h> // for isatty()
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <getopt.h>
#include <errno.h>
#include <assert.h>
#include <zlib.h>   // for crc32
#include <htslib/faidx.h>
#include <htslib/sam.h>
#include <htslib/hts.h>
#include "sam_header.h"
#include <htslib/khash_str2int.h>
#include "samtools.h"
#include <htslib/khash.h>
#include <htslib/kstring.h>
#include "stats_isize.h"
#include "sam_opts.h"
#include "bedidx.h"

#define BWA_MIN_RDLEN 35
#define DEFAULT_CHUNK_NO 8
#define DEFAULT_PAIR_MAX 10000
// From the spec
// If 0x4 is set, no assumptions can be made about RNAME, POS, CIGAR, MAPQ, bits 0x2, 0x10, 0x100 and 0x800, and the bit 0x20 of the previous read in the template.
#define IS_PAIRED_AND_MAPPED(bam) (((bam)->core.flag&BAM_FPAIRED) && !((bam)->core.flag&BAM_FUNMAP) && !((bam)->core.flag&BAM_FMUNMAP))
#define IS_PROPERLYPAIRED(bam) (((bam)->core.flag&(BAM_FPAIRED|BAM_FPROPER_PAIR)) == (BAM_FPAIRED|BAM_FPROPER_PAIR) && !((bam)->core.flag&BAM_FUNMAP))
#define IS_UNMAPPED(bam) ((bam)->core.flag&BAM_FUNMAP)
#define IS_REVERSE(bam) ((bam)->core.flag&BAM_FREVERSE)
#define IS_MATE_REVERSE(bam) ((bam)->core.flag&BAM_FMREVERSE)
#define IS_READ1(bam) ((bam)->core.flag&BAM_FREAD1)
#define IS_READ2(bam) ((bam)->core.flag&BAM_FREAD2)
#define IS_DUP(bam) ((bam)->core.flag&BAM_FDUP)
#define IS_ORIGINAL(bam) (((bam)->core.flag&(BAM_FSECONDARY|BAM_FSUPPLEMENTARY)) == 0)

// The GC-depth graph works as follows: split the reference sequence into
// segments and calculate GC content and depth in each bin. Then sort
// these segments by their GC and plot the depth distribution by means
// of 10th, 25th, etc. depth percentiles.
typedef struct
{
    float gc;
    uint32_t depth;
}
gc_depth_t;

// For coverage distribution, a simple pileup
typedef struct
{
    int64_t pos;
    int size, start;
    int *buffer;
}
round_buffer_t;

typedef struct { uint32_t from, to; } pos_t;
typedef struct
{
    int npos,mpos,cpos;
    pos_t *pos;
}
regions_t;

typedef struct
{
    uint64_t a;
    uint64_t c;
    uint64_t g;
    uint64_t t;
    uint64_t n;
    uint64_t other;
}
acgtno_count_t;

typedef struct
{
    // Auxiliary data
    int flag_require, flag_filter;
    faidx_t *fai;                   // Reference sequence for GC-depth graph
    int argc;                       // Command line arguments to be printed on the output
    char **argv;
    int gcd_bin_size;           // The size of GC-depth bin
    int nisize;         // The maximum insert size that the allocated array can hold - 0 indicates no limit
    int trim_qual;      // bwa trim quality
    float isize_main_bulk;  // There are always some unrealistically big insert sizes, report only the main part
    int cov_min,cov_max,cov_step;   // Minimum, maximum coverage and size of the coverage bins
    samFile* sam;
    bam_hdr_t* sam_header;

    // Filters
    int filter_readlen;

    // Misc
    char *split_tag;      // Tag on which to perform stats splitting
    char *split_prefix;   // Path or string prefix for filenames created when splitting
    int remove_overlaps;
    int cov_threshold;
}
stats_info_t;

typedef struct
{
    // Dimensions of the quality histogram holder (quals_1st,quals_2nd), GC content holder (gc_1st,gc_2nd),
    //  insert size histogram holder
    int nquals;         // The number of quality bins
    int nbases;         // The maximum sequence length the allocated array can hold
    int ngc;            // The size of gc_1st and gc_2nd
    int nindels;        // The maximum indel length for indel distribution

    // Arrays for the histogram data
    uint64_t *quals_1st, *quals_2nd;
    uint64_t *gc_1st, *gc_2nd;
    acgtno_count_t *acgtno_cycles_1st;
    acgtno_count_t *acgtno_cycles_2nd;
    uint64_t *read_lengths, *read_lengths_1st, *read_lengths_2nd;
    uint64_t *insertions, *deletions;
    uint64_t *ins_cycles_1st, *ins_cycles_2nd, *del_cycles_1st, *del_cycles_2nd;
    isize_t *isize;

    // The extremes encountered
    int max_len;            // Maximum read length
    int max_len_1st;        // Maximum read length for forward reads
    int max_len_2nd;        // Maximum read length for reverse reads
    int max_qual;           // Maximum quality
    int is_sorted;

    // Summary numbers
    uint64_t total_len;
    uint64_t total_len_1st;
    uint64_t total_len_2nd;
    uint64_t total_len_dup;
    uint64_t nreads_1st;
    uint64_t nreads_2nd;
    uint64_t nreads_filtered;
    uint64_t nreads_dup;
    uint64_t nreads_unmapped;
    uint64_t nreads_single_mapped;
    uint64_t nreads_paired_and_mapped;
    uint64_t nreads_properly_paired;
    uint64_t nreads_paired_tech;
    uint64_t nreads_anomalous;
    uint64_t nreads_mq0;
    uint64_t nbases_mapped;
    uint64_t nbases_mapped_cigar;
    uint64_t nbases_trimmed;  // bwa trimmed bases
    uint64_t nmismatches;
    uint64_t nreads_QCfailed, nreads_secondary;
    struct {
        uint32_t names, reads, quals;
    } checksum;

    // GC-depth related data
    uint32_t ngcd, igcd;        // The maximum number of GC depth bins and index of the current bin
    gc_depth_t *gcd;            // The GC-depth bins holder
    int32_t tid, gcd_pos;       // Position of the current bin
    int32_t pos;                // Position of the last read

    // Coverage distribution related data
    int ncov;                       // The number of coverage bins
    uint64_t *cov;                  // The coverage frequencies
    round_buffer_t cov_rbuf;        // Pileup round buffer

    // Mismatches by read cycle
    uint8_t *rseq_buf;              // A buffer for reference sequence to check the mismatches against
    int mrseq_buf;                  // The size of the buffer
    int32_t rseq_pos;               // The coordinate of the first base in the buffer
    int32_t nrseq_buf;              // The used part of the buffer
    uint64_t *mpc_buf;              // Mismatches per cycle

    // Target regions
    int nregions, reg_from, reg_to;
    regions_t *regions;

    // Auxiliary data
    double sum_qual;                // For calculating average quality value
    void *rg_hash;                  // Read groups to include, the array is null-terminated

    // Split
    char* split_name;

    stats_info_t* info;             // Pointer to options and settings struct
    pos_t *chunks;
    uint32_t nchunks;

    uint32_t pair_count;          // Number of active pairs in the pairing hash table
    uint32_t target_count;        // Number of bases covered by the target file
    uint32_t last_pair_tid;
    uint32_t last_read_flush;
}
stats_t;
KHASH_MAP_INIT_STR(c2stats, stats_t*)

typedef struct {
    uint32_t first;     // 1 - first read, 2 - second read
    uint32_t n, m;      // number of chunks, allocated chunks
    pos_t *chunks;      // chunk array of size m
} pair_t;
KHASH_MAP_INIT_STR(qn2pair, pair_t*)


static void error(const char *format, ...);
int is_in_regions(bam1_t *bam_line, stats_t *stats);
void realloc_buffers(stats_t *stats, int seq_len);

static int regions_lt(const void *r1, const void *r2) {
    int64_t from_diff = (int64_t)((pos_t *)r1)->from - (int64_t)((pos_t *)r2)->from;
    int64_t to_diff = (int64_t)((pos_t *)r1)->to - (int64_t)((pos_t *)r2)->to;

    return from_diff > 0 ? 1 : from_diff < 0 ? -1 : to_diff > 0 ? 1 : to_diff < 0 ? -1 : 0;
}

// Coverage distribution methods
static inline int coverage_idx(int min, int max, int n, int step, int depth)
{
    if ( depth < min )
        return 0;

    if ( depth > max )
        return n-1;

    return 1 + (depth - min) / step;
}

static inline int round_buffer_lidx2ridx(int offset, int size, int64_t refpos, int64_t pos)
{
    return (offset + (pos-refpos) % size) % size;
}

void round_buffer_flush(stats_t *stats, int64_t pos)
{
    int ibuf,idp;

    if ( pos==stats->cov_rbuf.pos )
        return;

    int64_t new_pos = pos;
    if ( pos==-1 || pos - stats->cov_rbuf.pos >= stats->cov_rbuf.size )
    {
        // Flush the whole buffer, but in sequential order,
        pos = stats->cov_rbuf.pos + stats->cov_rbuf.size - 1;
    }

    if ( pos < stats->cov_rbuf.pos )
        error("Expected coordinates in ascending order, got %ld after %ld\n", pos,stats->cov_rbuf.pos);

    int ifrom = stats->cov_rbuf.start;
    int ito = round_buffer_lidx2ridx(stats->cov_rbuf.start,stats->cov_rbuf.size,stats->cov_rbuf.pos,pos-1);
    if ( ifrom>ito )
    {
        for (ibuf=ifrom; ibuf<stats->cov_rbuf.size; ibuf++)
        {
            if ( !stats->cov_rbuf.buffer[ibuf] )
                continue;
            idp = coverage_idx(stats->info->cov_min,stats->info->cov_max,stats->ncov,stats->info->cov_step,stats->cov_rbuf.buffer[ibuf]);
            stats->cov[idp]++;
            stats->cov_rbuf.buffer[ibuf] = 0;
        }
        ifrom = 0;
    }
    for (ibuf=ifrom; ibuf<=ito; ibuf++)
    {
        if ( !stats->cov_rbuf.buffer[ibuf] )
            continue;
        idp = coverage_idx(stats->info->cov_min,stats->info->cov_max,stats->ncov,stats->info->cov_step,stats->cov_rbuf.buffer[ibuf]);
        stats->cov[idp]++;
        stats->cov_rbuf.buffer[ibuf] = 0;
    }
    stats->cov_rbuf.start = (new_pos==-1) ? 0 : round_buffer_lidx2ridx(stats->cov_rbuf.start,stats->cov_rbuf.size,stats->cov_rbuf.pos,pos);
    stats->cov_rbuf.pos   = new_pos;
}

void round_buffer_insert_read(round_buffer_t *rbuf, int64_t from, int64_t to)
{
    if ( to-from >= rbuf->size )
        error("The read length too big (%d), please increase the buffer length (currently %d)\n", to-from+1,rbuf->size);
    if ( from < rbuf->pos )
        error("The reads are not sorted (%ld comes after %ld).\n", from,rbuf->pos);

    int ifrom,ito,ibuf;
    ifrom = round_buffer_lidx2ridx(rbuf->start,rbuf->size,rbuf->pos,from);
    ito   = round_buffer_lidx2ridx(rbuf->start,rbuf->size,rbuf->pos,to);
    if ( ifrom>ito )
    {
        for (ibuf=ifrom; ibuf<rbuf->size; ibuf++)
            rbuf->buffer[ibuf]++;
        ifrom = 0;
    }
    for (ibuf=ifrom; ibuf<=ito; ibuf++)
        rbuf->buffer[ibuf]++;
}

// Calculate the number of bases in the read trimmed by BWA
int bwa_trim_read(int trim_qual, uint8_t *quals, int len, int reverse)
{
    if ( len<BWA_MIN_RDLEN ) return 0;

    // Although the name implies that the read cannot be trimmed to more than BWA_MIN_RDLEN,
    //  the calculation can in fact trim it to (BWA_MIN_RDLEN-1). (bwa_trim_read in bwa/bwaseqio.c).
    int max_trimmed = len - BWA_MIN_RDLEN + 1;
    int l, sum=0, max_sum=0, max_l=0;

    for (l=0; l<max_trimmed; l++)
    {
        sum += trim_qual - quals[ reverse ? l : len-1-l ];
        if ( sum<0 ) break;
        if ( sum>max_sum )
        {
            max_sum = sum;
            // This is the correct way, but bwa clips from some reason one base less
            // max_l   = l+1;
            max_l   = l;
        }
    }
    return max_l;
}


void count_indels(stats_t *stats,bam1_t *bam_line)
{
    int is_fwd = IS_REVERSE(bam_line) ? 0 : 1;
    int is_1st = IS_READ1(bam_line) ? 1 : 0;
    int icig;
    int icycle = 0;
    int read_len = bam_line->core.l_qseq;
    for (icig=0; icig<bam_line->core.n_cigar; icig++)
    {
        int cig  = bam_cigar_op(bam_get_cigar(bam_line)[icig]);
        int ncig = bam_cigar_oplen(bam_get_cigar(bam_line)[icig]);
        if ( !ncig ) continue;  // curiously, this can happen: 0D

        if ( cig==BAM_CINS )
        {
            int idx = is_fwd ? icycle : read_len-icycle-ncig;
            if ( idx<0 )
                error("FIXME: read_len=%d vs icycle=%d\n", read_len,icycle);
            if ( idx >= stats->nbases || idx<0 ) error("FIXME: %d vs %d, %s:%d %s\n", idx,stats->nbases, stats->info->sam_header->target_name[bam_line->core.tid],bam_line->core.pos+1,bam_get_qname(bam_line));
            if ( is_1st )
                stats->ins_cycles_1st[idx]++;
            else
                stats->ins_cycles_2nd[idx]++;
            icycle += ncig;
            if ( ncig<=stats->nindels )
                stats->insertions[ncig-1]++;
            continue;
        }
        if ( cig==BAM_CDEL )
        {
            int idx = is_fwd ? icycle-1 : read_len-icycle-1;
            if ( idx<0 ) continue;  // discard meaningless deletions
            if ( idx >= stats->nbases ) error("FIXME: %d vs %d\n", idx,stats->nbases);
            if ( is_1st )
                stats->del_cycles_1st[idx]++;
            else
                stats->del_cycles_2nd[idx]++;
            if ( ncig<=stats->nindels )
                stats->deletions[ncig-1]++;
            continue;
        }
        if ( cig!=BAM_CREF_SKIP && cig!=BAM_CHARD_CLIP && cig!=BAM_CPAD )
            icycle += ncig;
    }
}

int unclipped_length(bam1_t *bam_line)
{
    int icig, read_len = bam_line->core.l_qseq;
    for (icig=0; icig<bam_line->core.n_cigar; icig++)
    {
        int cig = bam_cigar_op(bam_get_cigar(bam_line)[icig]);
        if ( cig==BAM_CHARD_CLIP )
            read_len += bam_cigar_oplen(bam_get_cigar(bam_line)[icig]);
    }
    return read_len;
}

void count_mismatches_per_cycle(stats_t *stats, bam1_t *bam_line, int read_len)
{
    int is_fwd = IS_REVERSE(bam_line) ? 0 : 1;
    int icig,iread=0,icycle=0;
    int iref = bam_line->core.pos - stats->rseq_pos;
    uint8_t *read  = bam_get_seq(bam_line);
    uint8_t *quals = bam_get_qual(bam_line);
    uint64_t *mpc_buf = stats->mpc_buf;
    for (icig=0; icig<bam_line->core.n_cigar; icig++)
    {
        int cig  = bam_cigar_op(bam_get_cigar(bam_line)[icig]);
        int ncig = bam_cigar_oplen(bam_get_cigar(bam_line)[icig]);
        if ( cig==BAM_CINS )
        {
            iread  += ncig;
            icycle += ncig;
            continue;
        }
        if ( cig==BAM_CDEL )
        {
            iref += ncig;
            continue;
        }
        if ( cig==BAM_CSOFT_CLIP )
        {
            icycle += ncig;
            // Soft-clips are present in the sequence, but the position of the read marks a start of the sequence after clipping
            //   iref += ncig;
            iread  += ncig;
            continue;
        }
        if ( cig==BAM_CHARD_CLIP )
        {
            icycle += ncig;
            continue;
        }
        // Ignore H and N CIGARs. The letter are inserted e.g. by TopHat and often require very large
        //  chunk of refseq in memory. Not very frequent and not noticable in the stats.
        if ( cig==BAM_CREF_SKIP || cig==BAM_CHARD_CLIP || cig==BAM_CPAD ) continue;
        if ( cig!=BAM_CMATCH && cig!=BAM_CEQUAL && cig!=BAM_CDIFF ) // not relying on precalculated diffs
            error("TODO: cigar %d, %s:%d %s\n", cig,stats->info->sam_header->target_name[bam_line->core.tid],bam_line->core.pos+1,bam_get_qname(bam_line));

        if ( ncig+iref > stats->nrseq_buf )
            error("FIXME: %d+%d > %d, %s, %s:%d\n",ncig,iref,stats->nrseq_buf, bam_get_qname(bam_line),stats->info->sam_header->target_name[bam_line->core.tid],bam_line->core.pos+1);

        int im;
        for (im=0; im<ncig; im++)
        {
            uint8_t cread = bam_seqi(read,iread);
            uint8_t cref  = stats->rseq_buf[iref];

            // ---------------15
            // =ACMGRSVTWYHKDBN
            if ( cread==15 )
            {
                int idx = is_fwd ? icycle : read_len-icycle-1;
                if ( idx>stats->max_len )
                    error("mpc: %d>%d\n",idx,stats->max_len);
                idx = idx*stats->nquals;
                if ( idx>=stats->nquals*stats->nbases )
                    error("FIXME: mpc_buf overflow\n");
                mpc_buf[idx]++;
            }
            else if ( cref && cread && cref!=cread )
            {
                uint8_t qual = quals[iread] + 1;
                if ( qual>=stats->nquals )
                    error("TODO: quality too high %d>=%d (%s %d %s)\n", qual,stats->nquals, stats->info->sam_header->target_name[bam_line->core.tid],bam_line->core.pos+1,bam_get_qname(bam_line));

                int idx = is_fwd ? icycle : read_len-icycle-1;
                if ( idx>stats->max_len )
                    error("mpc: %d>%d (%s %d %s)\n",idx,stats->max_len,stats->info->sam_header->target_name[bam_line->core.tid],bam_line->core.pos+1,bam_get_qname(bam_line));

                idx = idx*stats->nquals + qual;
                if ( idx>=stats->nquals*stats->nbases )
                    error("FIXME: mpc_buf overflow\n");
                mpc_buf[idx]++;
            }

            iref++;
            iread++;
            icycle++;
        }
    }
}

void read_ref_seq(stats_t *stats, int32_t tid, int32_t pos)
{
    int i, fai_ref_len;
    char *fai_ref = faidx_fetch_seq(stats->info->fai, stats->info->sam_header->target_name[tid], pos, pos+stats->mrseq_buf-1, &fai_ref_len);
    if ( fai_ref_len<0 ) error("Failed to fetch the sequence \"%s\"\n", stats->info->sam_header->target_name[tid]);

    uint8_t *ptr = stats->rseq_buf;
    for (i=0; i<fai_ref_len; i++)
    {
        // Conversion between uint8_t coding and ACGT
        //      -12-4---8-------
        //      =ACMGRSVTWYHKDBN
        switch (fai_ref[i])
        {
            case 'A':
            case 'a': *ptr = 1; break;
            case 'C':
            case 'c': *ptr = 2; break;
            case 'G':
            case 'g': *ptr = 4; break;
            case 'T':
            case 't': *ptr = 8; break;
            default:  *ptr = 0; break;
        }
        ptr++;
    }
    free(fai_ref);

    if ( fai_ref_len < stats->mrseq_buf ) memset(ptr,0, stats->mrseq_buf - fai_ref_len);
    stats->nrseq_buf = fai_ref_len;
    stats->rseq_pos  = pos;
    stats->tid       = tid;
}

float fai_gc_content(stats_t *stats, int pos, int len)
{
    uint32_t gc,count,c;
    int i = pos - stats->rseq_pos, ito = i + len;
    assert( i>=0 );

    if (  ito > stats->nrseq_buf ) ito = stats->nrseq_buf;

    // Count GC content
    gc = count = 0;
    for (; i<ito; i++)
    {
        c = stats->rseq_buf[i];
        if ( c==2 || c==4 )
        {
            gc++;
            count++;
        }
        else if ( c==1 || c==8 )
            count++;
    }
    return count ? (float)gc/count : 0;
}

void realloc_rseq_buffer(stats_t *stats)
{
    int n = stats->nbases*10;
    if ( stats->info->gcd_bin_size > n ) n = stats->info->gcd_bin_size;
    if ( stats->mrseq_buf<n )
    {
        stats->rseq_buf = realloc(stats->rseq_buf,sizeof(uint8_t)*n);
        stats->mrseq_buf = n;
    }
}

void realloc_gcd_buffer(stats_t *stats, int seq_len)
{
    hts_expand0(gc_depth_t,stats->igcd+1,stats->ngcd,stats->gcd);
    realloc_rseq_buffer(stats);
}

void realloc_buffers(stats_t *stats, int seq_len)
{
    int n = 2*(1 + seq_len - stats->nbases) + stats->nbases;

    stats->quals_1st = realloc(stats->quals_1st, n*stats->nquals*sizeof(uint64_t));
    if ( !stats->quals_1st )
        error("Could not realloc buffers, the sequence too long: %d (%ld)\n", seq_len,n*stats->nquals*sizeof(uint64_t));
    memset(stats->quals_1st + stats->nbases*stats->nquals, 0, (n-stats->nbases)*stats->nquals*sizeof(uint64_t));

    stats->quals_2nd = realloc(stats->quals_2nd, n*stats->nquals*sizeof(uint64_t));
    if ( !stats->quals_2nd )
        error("Could not realloc buffers, the sequence too long: %d (2x%ld)\n", seq_len,n*stats->nquals*sizeof(uint64_t));
    memset(stats->quals_2nd + stats->nbases*stats->nquals, 0, (n-stats->nbases)*stats->nquals*sizeof(uint64_t));

    if ( stats->mpc_buf )
    {
        stats->mpc_buf = realloc(stats->mpc_buf, n*stats->nquals*sizeof(uint64_t));
        if ( !stats->mpc_buf )
            error("Could not realloc buffers, the sequence too long: %d (%ld)\n", seq_len,n*stats->nquals*sizeof(uint64_t));
        memset(stats->mpc_buf + stats->nbases*stats->nquals, 0, (n-stats->nbases)*stats->nquals*sizeof(uint64_t));
    }

    stats->acgtno_cycles_1st = realloc(stats->acgtno_cycles_1st, n*sizeof(acgtno_count_t));
    if ( !stats->acgtno_cycles_1st )
        error("Could not realloc buffers, the sequence too long: %d (%ld)\n", seq_len, n*sizeof(acgtno_count_t));
    memset(stats->acgtno_cycles_1st + stats->nbases, 0, (n-stats->nbases)*sizeof(acgtno_count_t));

    stats->acgtno_cycles_2nd = realloc(stats->acgtno_cycles_2nd, n*sizeof(acgtno_count_t));
    if ( !stats->acgtno_cycles_2nd )
        error("Could not realloc buffers, the sequence too long: %d (%ld)\n", seq_len, n*sizeof(acgtno_count_t));
    memset(stats->acgtno_cycles_2nd + stats->nbases, 0, (n-stats->nbases)*sizeof(acgtno_count_t));

    stats->read_lengths = realloc(stats->read_lengths, n*sizeof(uint64_t));
    if ( !stats->read_lengths )
        error("Could not realloc buffers, the sequence too long: %d (%ld)\n", seq_len,n*sizeof(uint64_t));
    memset(stats->read_lengths + stats->nbases, 0, (n-stats->nbases)*sizeof(uint64_t));

    stats->read_lengths_1st = realloc(stats->read_lengths_1st, n*sizeof(uint64_t));
    if ( !stats->read_lengths_1st )
        error("Could not realloc buffers, the sequence too long: %d (%ld)\n", seq_len,n*sizeof(uint64_t));
    memset(stats->read_lengths_1st + stats->nbases, 0, (n-stats->nbases)*sizeof(uint64_t));

    stats->read_lengths_2nd = realloc(stats->read_lengths_2nd, n*sizeof(uint64_t));
    if ( !stats->read_lengths_2nd )
        error("Could not realloc buffers, the sequence too long: %d (%ld)\n", seq_len,n*sizeof(uint64_t));
    memset(stats->read_lengths_2nd + stats->nbases, 0, (n-stats->nbases)*sizeof(uint64_t));

    stats->insertions = realloc(stats->insertions, n*sizeof(uint64_t));
    if ( !stats->insertions )
        error("Could not realloc buffers, the sequence too long: %d (%ld)\n", seq_len,n*sizeof(uint64_t));
    memset(stats->insertions + stats->nbases, 0, (n-stats->nbases)*sizeof(uint64_t));

    stats->deletions = realloc(stats->deletions, n*sizeof(uint64_t));
    if ( !stats->deletions )
        error("Could not realloc buffers, the sequence too long: %d (%ld)\n", seq_len,n*sizeof(uint64_t));
    memset(stats->deletions + stats->nbases, 0, (n-stats->nbases)*sizeof(uint64_t));

    stats->ins_cycles_1st = realloc(stats->ins_cycles_1st, (n+1)*sizeof(uint64_t));
    if ( !stats->ins_cycles_1st )
        error("Could not realloc buffers, the sequence too long: %d (%ld)\n", seq_len,(n+1)*sizeof(uint64_t));
    memset(stats->ins_cycles_1st + stats->nbases + 1, 0, (n-stats->nbases)*sizeof(uint64_t));

    stats->ins_cycles_2nd = realloc(stats->ins_cycles_2nd, (n+1)*sizeof(uint64_t));
    if ( !stats->ins_cycles_2nd )
        error("Could not realloc buffers, the sequence too long: %d (%ld)\n", seq_len,(n+1)*sizeof(uint64_t));
    memset(stats->ins_cycles_2nd + stats->nbases + 1, 0, (n-stats->nbases)*sizeof(uint64_t));

    stats->del_cycles_1st = realloc(stats->del_cycles_1st, (n+1)*sizeof(uint64_t));
    if ( !stats->del_cycles_1st )
        error("Could not realloc buffers, the sequence too long: %d (%ld)\n", seq_len,(n+1)*sizeof(uint64_t));
    memset(stats->del_cycles_1st + stats->nbases + 1, 0, (n-stats->nbases)*sizeof(uint64_t));

    stats->del_cycles_2nd = realloc(stats->del_cycles_2nd, (n+1)*sizeof(uint64_t));
    if ( !stats->del_cycles_2nd )
        error("Could not realloc buffers, the sequence too long: %d (%ld)\n", seq_len,(n+1)*sizeof(uint64_t));
    memset(stats->del_cycles_2nd + stats->nbases + 1, 0, (n-stats->nbases)*sizeof(uint64_t));

    stats->nbases = n;

    // Realloc the coverage distribution buffer
    int *rbuffer = calloc(sizeof(int),seq_len*5);
    n = stats->cov_rbuf.size-stats->cov_rbuf.start;
    memcpy(rbuffer,stats->cov_rbuf.buffer+stats->cov_rbuf.start,n);
    if ( stats->cov_rbuf.start>1 )
        memcpy(rbuffer+n,stats->cov_rbuf.buffer,stats->cov_rbuf.start);
    stats->cov_rbuf.start = 0;
    free(stats->cov_rbuf.buffer);
    stats->cov_rbuf.buffer = rbuffer;
    stats->cov_rbuf.size = seq_len*5;

    realloc_rseq_buffer(stats);
}

void update_checksum(bam1_t *bam_line, stats_t *stats)
{
    uint8_t *name = (uint8_t*) bam_get_qname(bam_line);
    int len = 0;
    while ( name[len] ) len++;
    stats->checksum.names +=  crc32(0L, name, len);

    int seq_len = bam_line->core.l_qseq;
    if ( !seq_len ) return;

    uint8_t *seq = bam_get_seq(bam_line);
    stats->checksum.reads += crc32(0L, seq, (seq_len+1)/2);

    uint8_t *qual = bam_get_qual(bam_line);
    stats->checksum.quals += crc32(0L, qual, (seq_len+1)/2);
}

// These stats should only be calculated for the original reads ignoring
// supplementary artificial reads otherwise we'll accidentally double count
void collect_orig_read_stats(bam1_t *bam_line, stats_t *stats, int* gc_count_out)
{
    int seq_len = bam_line->core.l_qseq;
    stats->total_len += seq_len; // This ignores clipping so only count primary

    if ( bam_line->core.flag & BAM_FQCFAIL ) stats->nreads_QCfailed++;
    if ( bam_line->core.flag & BAM_FPAIRED ) stats->nreads_paired_tech++;

    // Count GC and ACGT per cycle. Note that cycle is approximate, clipping is ignored
    uint8_t *seq  = bam_get_seq(bam_line);
    int i, read_cycle, gc_count = 0, reverse = IS_REVERSE(bam_line), is_first = IS_READ1(bam_line);
    for (i=0; i<seq_len; i++)
    {
        // Read cycle for current index
        read_cycle = (reverse ? seq_len-i-1 : i);

        // Conversion from uint8_t coding:
        //      -12-4---8------5
        //      =ACMGRSVTWYHKDBN
        switch (bam_seqi(seq, i)) {
        case 1:
            is_first ? stats->acgtno_cycles_1st[ read_cycle ].a++ : stats->acgtno_cycles_2nd[ read_cycle ].a++;
            break;
        case 2:
            is_first ? stats->acgtno_cycles_1st[ read_cycle ].c++ : stats->acgtno_cycles_2nd[ read_cycle ].c++;
            gc_count++;
            break;
        case 4:
            is_first ? stats->acgtno_cycles_1st[ read_cycle ].g++ : stats->acgtno_cycles_2nd[ read_cycle ].g++;
            gc_count++;
            break;
        case 8:
            is_first ? stats->acgtno_cycles_1st[ read_cycle ].t++ : stats->acgtno_cycles_2nd[ read_cycle ].t++;
            break;
        case 15:
            is_first ? stats->acgtno_cycles_1st[ read_cycle ].n++ : stats->acgtno_cycles_2nd[ read_cycle ].n++;
            break;
        default:
            /*
             * count "=" sequences in "other" along
             * with MRSVWYHKDB ambiguity codes
             */
            is_first ? stats->acgtno_cycles_1st[ read_cycle ].other++ : stats->acgtno_cycles_2nd[ read_cycle ].other++;
            break;
        }
    }
    int gc_idx_min = gc_count*(stats->ngc-1)/seq_len;
    int gc_idx_max = (gc_count+1)*(stats->ngc-1)/seq_len;
    if ( gc_idx_max >= stats->ngc ) gc_idx_max = stats->ngc - 1;

    // Determine which array (1st or 2nd read) will these stats go to,
    //  trim low quality bases from end the same way BWA does,
    //  fill GC histogram
    uint64_t *quals;
    uint8_t *bam_quals = bam_get_qual(bam_line);
    if ( IS_READ2(bam_line) )
    {
        quals  = stats->quals_2nd;
        stats->nreads_2nd++;
        stats->total_len_2nd += seq_len;
        for (i=gc_idx_min; i<gc_idx_max; i++)
            stats->gc_2nd[i]++;
    }
    else
    {
        quals = stats->quals_1st;
        stats->nreads_1st++;
        stats->total_len_1st += seq_len;
        for (i=gc_idx_min; i<gc_idx_max; i++)
            stats->gc_1st[i]++;
    }
    if ( stats->info->trim_qual>0 )
        stats->nbases_trimmed += bwa_trim_read(stats->info->trim_qual, bam_quals, seq_len, reverse);

    // Quality histogram and average quality. Clipping is neglected.
    for (i=0; i<seq_len; i++)
    {
        uint8_t qual = bam_quals[ reverse ? seq_len-i-1 : i];
        if ( qual>=stats->nquals )
            error("TODO: quality too high %d>=%d (%s %d %s)\n", qual,stats->nquals,stats->info->sam_header->target_name[bam_line->core.tid],bam_line->core.pos+1,bam_get_qname(bam_line));
        if ( qual>stats->max_qual )
            stats->max_qual = qual;

        quals[ i*stats->nquals+qual ]++;
        stats->sum_qual += qual;
    }

    // Look at the flags and increment appropriate counters (mapped, paired, etc)
    if ( IS_UNMAPPED(bam_line) )
    {
        stats->nreads_unmapped++;
    }
    else
    {
        stats->nbases_mapped += seq_len; // This ignores clipping so only count primary

        if ( !bam_line->core.qual )
            stats->nreads_mq0++;
        if ( !IS_PAIRED_AND_MAPPED(bam_line) )
            stats->nreads_single_mapped++;
        else
        {
            stats->nreads_paired_and_mapped++;

            if (IS_PROPERLYPAIRED(bam_line)) stats->nreads_properly_paired++;

            if ( bam_line->core.tid!=bam_line->core.mtid )
                stats->nreads_anomalous++;
        }
    }
    *gc_count_out = gc_count;
}

static int cleanup_overlaps(khash_t(qn2pair) *read_pairs, int max) {
    if ( !read_pairs )
        return 0;

    int count = 0;
    khint_t k;
    for (k = kh_begin(read_pairs); k < kh_end(read_pairs); k++) {
        if ( kh_exist(read_pairs, k) ) {
            char *key = (char *)kh_key(read_pairs, k);
            pair_t *val = kh_val(read_pairs, k);
            if ( val && val->chunks ) {
                if ( val->chunks[val->n-1].to < max ) {
                    free(val->chunks);
                    free(val);
                    free(key);
                    kh_del(qn2pair, read_pairs, k);
                    count++;
                }
            } else {
                free(key);
                kh_del(qn2pair, read_pairs, k);
                count++;
            }
        }
    }
    if ( max == INT_MAX )
        kh_destroy(qn2pair, read_pairs);

    return count;
}

static void remove_overlaps(bam1_t *bam_line, khash_t(qn2pair) *read_pairs, stats_t *stats, int pmin, int pmax) {
    if ( !bam_line || !read_pairs || !stats )
        return;

    uint32_t first = (IS_READ1(bam_line) > 0 ? 1 : 0) + (IS_READ2(bam_line) > 0 ? 2 : 0) ;
    if ( !(bam_line->core.flag & BAM_FPAIRED) ||
         (bam_line->core.flag & BAM_FMUNMAP) ||
         (abs(bam_line->core.isize) >= 2*bam_line->core.l_qseq) ||
         (first != 1 && first != 2) ) {
        if ( pmin >= 0 )
            round_buffer_insert_read(&(stats->cov_rbuf), pmin, pmax-1);
        return;
    }

    char *qname = bam_get_qname(bam_line);
    if ( !qname ) {
        fprintf(stderr, "Error retrieving qname for line starting at pos %d\n", bam_line->core.pos);
        return;
    }

    khint_t k = kh_get(qn2pair, read_pairs, qname);
    if ( k == kh_end(read_pairs) ) { //first chunk from this template
        if ( pmin == -1 )
            return;

        int ret;
        char *s = strdup(qname);
        if ( !s ) {
            fprintf(stderr, "Error allocating memory\n");
            return;
        }

        k = kh_put(qn2pair, read_pairs, s, &ret);
        if ( -1 == ret ) {
            fprintf(stderr, "Error inserting read '%s' in pair hash table\n", qname);
            return;
        }

        pair_t *pc = calloc(1, sizeof(pair_t));
        if ( !pc ) {
            fprintf(stderr, "Error allocating memory\n");
            return;
        }

        pc->m = DEFAULT_CHUNK_NO;
        pc->chunks = calloc(pc->m, sizeof(pos_t));
        if ( !pc->chunks ) {
            fprintf(stderr, "Error allocating memory\n");
            return;
        }

        pc->chunks[0].from = pmin;
        pc->chunks[0].to = pmax;
        pc->n = 1;
        pc->first = first;

        kh_val(read_pairs, k) = pc;
        stats->pair_count++;
    } else { //template already present
        pair_t *pc = kh_val(read_pairs, k);
        if ( !pc ) {
            fprintf(stderr, "Invalid hash table entry\n");
            return;
        }

        if ( first == pc->first ) { //chunk from an existing line
            if ( pmin == -1 )
                return;

            if ( pc->n == pc->m ) {
                pos_t *tmp = realloc(pc->chunks, (pc->m<<1)*sizeof(pos_t));
                if ( !tmp ) {
                    fprintf(stderr, "Error allocating memory\n");
                    return;
                }
                pc->chunks = tmp;
                pc->m<<=1;
            }

            pc->chunks[pc->n].from = pmin;
            pc->chunks[pc->n].to = pmax;
            pc->n++;
        } else { //the other line, check for overlapping
            if ( pmin == -1 && kh_exist(read_pairs, k) ) { //job done, delete entry
                char *key = (char *)kh_key(read_pairs, k);
                pair_t *val = kh_val(read_pairs, k);
                if ( val) {
                    free(val->chunks);
                    free(val);
                }
                free(key);
                kh_del(qn2pair, read_pairs, k);
                stats->pair_count--;
                return;
            }

            int i;
            for (i=0; i<pc->n; i++) {
                if ( pmin >= pc->chunks[i].to )
                    continue;

                if ( pmax <= pc->chunks[i].from ) //no overlap
                    break;

                if ( pmin < pc->chunks[i].from ) { //overlap at the beginning
                    round_buffer_insert_read(&(stats->cov_rbuf), pmin, pc->chunks[i].from-1);
                    pmin = pc->chunks[i].from;
                }

                if ( pmax <= pc->chunks[i].to ) { //completely contained
                    stats->nbases_mapped_cigar -= (pmax - pmin);
                    return;
                } else {                           //overlap at the end
                    stats->nbases_mapped_cigar -= (pc->chunks[i].to - pmin);
                    pmin = pc->chunks[i].to;
                }
            }
        }
    }
    round_buffer_insert_read(&(stats->cov_rbuf), pmin, pmax-1);
}

void collect_stats(bam1_t *bam_line, stats_t *stats, khash_t(qn2pair) *read_pairs)
{
    if ( stats->rg_hash )
    {
        const uint8_t *rg = bam_aux_get(bam_line, "RG");
        if ( !rg ) return;  // certain read groups were requested but this record has none
        if ( !khash_str2int_has_key(stats->rg_hash, (const char*)(rg + 1)) ) return;
    }
    if ( stats->info->flag_require && (bam_line->core.flag & stats->info->flag_require)!=stats->info->flag_require )
    {
        stats->nreads_filtered++;
        return;
    }
    if ( stats->info->flag_filter && (bam_line->core.flag & stats->info->flag_filter) )
    {
        stats->nreads_filtered++;
        return;
    }
    if ( !is_in_regions(bam_line,stats) )
        return;
    if ( stats->info->filter_readlen!=-1 && bam_line->core.l_qseq!=stats->info->filter_readlen )
        return;

    update_checksum(bam_line, stats);

    // Secondary reads don't count for most stats purposes
    if ( bam_line->core.flag & BAM_FSECONDARY )
    {
        stats->nreads_secondary++;
        return;
    }

    // If line has no sequence cannot continue
    int seq_len = bam_line->core.l_qseq;
    if ( !seq_len ) return;

    if ( IS_DUP(bam_line) )
    {
        stats->total_len_dup += seq_len;
        stats->nreads_dup++;
    }

    int read_len = unclipped_length(bam_line);
    if ( read_len >= stats->nbases )
        realloc_buffers(stats,read_len);
    // Update max_len observed
    if ( stats->max_len<read_len )
        stats->max_len = read_len;
    if ( IS_READ1(bam_line) && stats->max_len_1st < read_len )
        stats->max_len_1st = read_len;
    if ( IS_READ2(bam_line) && stats->max_len_2nd < read_len )
        stats->max_len_2nd = read_len;

    int i;
    int gc_count = 0;

    // These stats should only be calculated for the original reads ignoring supplementary artificial reads
    // otherwise we'll accidentally double count
    if ( IS_ORIGINAL(bam_line) )
    {
        stats->read_lengths[read_len]++;
        if ( IS_READ1(bam_line) ) stats->read_lengths_1st[read_len]++;
        if ( IS_READ2(bam_line) ) stats->read_lengths_2nd[read_len]++;
        collect_orig_read_stats(bam_line, stats, &gc_count);
    }

    // Look at the flags and increment appropriate counters (mapped, paired, etc)
    if ( IS_UNMAPPED(bam_line) ) return;

    count_indels(stats, bam_line);

    if ( IS_PAIRED_AND_MAPPED(bam_line) )
    {
        // The insert size is tricky, because for long inserts the libraries are
        // prepared differently and the pairs point in other direction. BWA does
        // not set the paired flag for them.  Similar thing is true also for 454
        // reads. Mates mapped to different chromosomes have isize==0.
        int32_t isize = bam_line->core.isize;
        if ( isize<0 ) isize = -isize;
        if ( stats->info->nisize > 0 && isize > stats->info->nisize )
            isize = stats->info->nisize;
        if ( isize>0 || bam_line->core.tid==bam_line->core.mtid )
        {
            int pos_fst = bam_line->core.mpos - bam_line->core.pos;
            int is_fst  = IS_READ1(bam_line) ? 1 : -1;
            int is_fwd  = IS_REVERSE(bam_line) ? -1 : 1;
            int is_mfwd = IS_MATE_REVERSE(bam_line) ? -1 : 1;

            if ( is_fwd*is_mfwd>0 )
                stats->isize->inc_other(stats->isize->data, isize);
            else if ( is_fst*pos_fst>=0 )
            {
                if ( is_fst*is_fwd>0 )
                    stats->isize->inc_inward(stats->isize->data, isize);
                else
                    stats->isize->inc_outward(stats->isize->data, isize);
            }
            else if ( is_fst*pos_fst<0 )
            {
                if ( is_fst*is_fwd>0 )
                    stats->isize->inc_outward(stats->isize->data, isize);
                else
                    stats->isize->inc_inward(stats->isize->data, isize);
            }
        }
    }

    // Number of mismatches
    uint8_t *nm = bam_aux_get(bam_line,"NM");
    if (nm)
        stats->nmismatches += bam_aux2i(nm);

    // Number of mapped bases from cigar
    if ( bam_line->core.n_cigar == 0)
        error("FIXME: mapped read with no cigar?\n");
    int readlen=seq_len;
    if ( stats->regions )
    {
        // Count only on-target bases
        int iref = bam_line->core.pos + 1;
        for (i=0; i<bam_line->core.n_cigar; i++)
        {
            int cig  = bam_cigar_op(bam_get_cigar(bam_line)[i]);
            int ncig = bam_cigar_oplen(bam_get_cigar(bam_line)[i]);
            if ( !ncig ) continue;  // curiously, this can happen: 0D
            if ( cig==BAM_CDEL ) readlen += ncig;
            else if ( cig==BAM_CMATCH || cig==BAM_CEQUAL || cig==BAM_CDIFF )
            {
                if ( iref < stats->reg_from ) ncig -= stats->reg_from-iref;
                else if ( iref+ncig-1 > stats->reg_to ) ncig -= iref+ncig-1 - stats->reg_to;
                if ( ncig<0 ) ncig = 0;
                stats->nbases_mapped_cigar += ncig;
                iref += bam_cigar_oplen(bam_get_cigar(bam_line)[i]);
            }
            else if ( cig==BAM_CINS )
            {
                iref += ncig;
                if ( iref>=stats->reg_from && iref<=stats->reg_to )
                    stats->nbases_mapped_cigar += ncig;
            }
        }
    }
    else
    {
        // Count the whole read
        for (i=0; i<bam_line->core.n_cigar; i++)
        {
            int cig  = bam_cigar_op(bam_get_cigar(bam_line)[i]);
            if ( cig==BAM_CMATCH || cig==BAM_CINS || cig==BAM_CEQUAL || cig==BAM_CDIFF )
                stats->nbases_mapped_cigar += bam_cigar_oplen(bam_get_cigar(bam_line)[i]);
            if ( cig==BAM_CDEL )
                readlen += bam_cigar_oplen(bam_get_cigar(bam_line)[i]);
        }
    }

    if ( stats->tid==bam_line->core.tid && bam_line->core.pos<stats->pos )
        stats->is_sorted = 0;
    stats->pos = bam_line->core.pos;

    if ( stats->is_sorted )
    {
        if ( stats->tid==-1 || stats->tid!=bam_line->core.tid ) {
            round_buffer_flush(stats, -1);
        }

        //cleanup the pair hash table to free memory
        stats->last_read_flush++;
        if ( stats->pair_count > DEFAULT_PAIR_MAX && stats->last_read_flush > DEFAULT_PAIR_MAX) {
            stats->pair_count -= cleanup_overlaps(read_pairs, bam_line->core.pos);
            stats->last_read_flush = 0;
        }

        if ( stats->last_pair_tid != bam_line->core.tid) {
            stats->pair_count -= cleanup_overlaps(read_pairs, INT_MAX-1);
            stats->last_pair_tid = bam_line->core.tid;
            stats->last_read_flush = 0;
        }

        // Mismatches per cycle and GC-depth graph. For simplicity, reads overlapping GCD bins
        //  are not splitted which results in up to seq_len-1 overlaps. The default bin size is
        //  20kbp, so the effect is negligible.
        if ( stats->info->fai )
        {
            int inc_ref = 0, inc_gcd = 0;
            // First pass or new chromosome
            if ( stats->rseq_pos==-1 || stats->tid != bam_line->core.tid ) { inc_ref=1; inc_gcd=1; }
            // Read goes beyond the end of the rseq buffer
            else if ( stats->rseq_pos+stats->nrseq_buf < bam_line->core.pos+readlen ) { inc_ref=1; inc_gcd=1; }
            // Read overlaps the next gcd bin
            else if ( stats->gcd_pos+stats->info->gcd_bin_size < bam_line->core.pos+readlen )
            {
                inc_gcd = 1;
                if ( stats->rseq_pos+stats->nrseq_buf < bam_line->core.pos+stats->info->gcd_bin_size ) inc_ref = 1;
            }
            if ( inc_gcd )
            {
                stats->igcd++;
                if ( stats->igcd >= stats->ngcd )
                    realloc_gcd_buffer(stats, readlen);
                if ( inc_ref )
                    read_ref_seq(stats,bam_line->core.tid,bam_line->core.pos);
                stats->gcd_pos = bam_line->core.pos;
                stats->gcd[ stats->igcd ].gc = fai_gc_content(stats, stats->gcd_pos, stats->info->gcd_bin_size);
            }

            count_mismatches_per_cycle(stats,bam_line,read_len);
        }
        // No reference and first pass, new chromosome or sequence going beyond the end of the gcd bin
        else if ( stats->gcd_pos==-1 || stats->tid != bam_line->core.tid || bam_line->core.pos - stats->gcd_pos > stats->info->gcd_bin_size )
        {
            // First pass or a new chromosome
            stats->tid     = bam_line->core.tid;
            stats->gcd_pos = bam_line->core.pos;
            stats->igcd++;
            if ( stats->igcd >= stats->ngcd )
                realloc_gcd_buffer(stats, readlen);
        }
        stats->gcd[ stats->igcd ].depth++;
        // When no reference sequence is given, approximate the GC from the read (much shorter window, but otherwise OK)
        if ( !stats->info->fai )
            stats->gcd[ stats->igcd ].gc += (float) gc_count / seq_len;

        // Coverage distribution graph
        round_buffer_flush(stats,bam_line->core.pos);
        if ( stats->regions ) {
            uint32_t p = bam_line->core.pos, pnew, pmin, pmax, j;
            pmin = pmax = i = j = 0;
            while ( j < bam_line->core.n_cigar && i < stats->nchunks ) {
                int op = bam_cigar_op(bam_get_cigar(bam_line)[j]);
                int oplen = bam_cigar_oplen(bam_get_cigar(bam_line)[j]);
                switch(op) {
                case BAM_CMATCH:
                case BAM_CEQUAL:
                case BAM_CDIFF:
                    pmin = MAX(p, stats->chunks[i].from-1);
                    pmax = MIN(p+oplen, stats->chunks[i].to);
                    if ( pmax >= pmin ) {
                        if ( stats->info->remove_overlaps )
                            remove_overlaps(bam_line, read_pairs, stats, pmin, pmax);
                        else
                            round_buffer_insert_read(&(stats->cov_rbuf), pmin, pmax-1);
                    }
                    break;
                case BAM_CDEL:
                    break;
                }
                pnew = p + (bam_cigar_type(op)&2 ? oplen : 0); // consumes reference

                if ( pnew >= stats->chunks[i].to ) {
                    // go to the next chunk
                    i++;
                } else {
                    // go to the next CIGAR op
                    j++;
                    p = pnew;
                }
            }
        } else {
            uint32_t p = bam_line->core.pos, j;
            for (j = 0; j < bam_line->core.n_cigar; j++) {
                int op = bam_cigar_op(bam_get_cigar(bam_line)[j]);
                int oplen = bam_cigar_oplen(bam_get_cigar(bam_line)[j]);
                switch(op) {
                case BAM_CMATCH:
                case BAM_CEQUAL:
                case BAM_CDIFF:
                    if ( stats->info->remove_overlaps )
                        remove_overlaps(bam_line, read_pairs, stats, p, p+oplen);
                    else
                        round_buffer_insert_read(&(stats->cov_rbuf), p, p+oplen-1);
                    break;
                case BAM_CDEL:
                    break;
                }
                p += bam_cigar_type(op)&2 ? oplen : 0; // consumes reference
            }
        }
        if ( stats->info->remove_overlaps )
           remove_overlaps(bam_line, read_pairs, stats, -1, -1); //remove the line from the hash table
    }
}

// Sort by GC and depth
#define GCD_t(x) ((gc_depth_t *)x)
static int gcd_cmp(const void *a, const void *b)
{
    if ( GCD_t(a)->gc < GCD_t(b)->gc ) return -1;
    if ( GCD_t(a)->gc > GCD_t(b)->gc ) return 1;
    if ( GCD_t(a)->depth < GCD_t(b)->depth ) return -1;
    if ( GCD_t(a)->depth > GCD_t(b)->depth ) return 1;
    return 0;
}
#undef GCD_t

float gcd_percentile(gc_depth_t *gcd, int N, int p)
{
    float n,d;
    int k;

    n = p*(N+1)/100;
    k = n;
    if ( k<=0 )
        return gcd[0].depth;
    if ( k>=N )
        return gcd[N-1].depth;

    d = n - k;
    return gcd[k-1].depth + d*(gcd[k].depth - gcd[k-1].depth);
}

void output_stats(FILE *to, stats_t *stats, int sparse)
{
    // Calculate average insert size and standard deviation (from the main bulk data only)
    int isize, ibulk=0, icov;
    uint64_t nisize=0, nisize_inward=0, nisize_outward=0, nisize_other=0, cov_sum=0;
    double bulk=0, avg_isize=0, sd_isize=0;
    for (isize=0; isize<stats->isize->nitems(stats->isize->data); isize++)
    {
        // Each pair was counted twice
        stats->isize->set_inward(stats->isize->data, isize, stats->isize->inward(stats->isize->data, isize) * 0.5);
        stats->isize->set_outward(stats->isize->data, isize, stats->isize->outward(stats->isize->data, isize) * 0.5);
        stats->isize->set_other(stats->isize->data, isize, stats->isize->other(stats->isize->data, isize) * 0.5);

        nisize_inward += stats->isize->inward(stats->isize->data, isize);
        nisize_outward += stats->isize->outward(stats->isize->data, isize);
        nisize_other += stats->isize->other(stats->isize->data, isize);
        nisize += stats->isize->inward(stats->isize->data, isize) + stats->isize->outward(stats->isize->data, isize) + stats->isize->other(stats->isize->data, isize);
    }

    for (isize=0; isize<stats->isize->nitems(stats->isize->data); isize++)
    {
        uint64_t num = stats->isize->inward(stats->isize->data, isize) +  stats->isize->outward(stats->isize->data, isize) + stats->isize->other(stats->isize->data, isize);
        if (num > 0) ibulk = isize + 1;
        bulk += num;
        avg_isize += isize * (stats->isize->inward(stats->isize->data, isize) +  stats->isize->outward(stats->isize->data, isize) + stats->isize->other(stats->isize->data, isize));

        if ( bulk/nisize > stats->info->isize_main_bulk )
        {
            ibulk  = isize+1;
            nisize = bulk;
            break;
        }
    }
    avg_isize /= nisize ? nisize : 1;
    for (isize=1; isize<ibulk; isize++)
        sd_isize += (stats->isize->inward(stats->isize->data, isize) + stats->isize->outward(stats->isize->data, isize) +stats->isize->other(stats->isize->data, isize)) * (isize-avg_isize)*(isize-avg_isize) / (nisize ? nisize : 1);
    sd_isize = sqrt(sd_isize);

    fprintf(to, "# This file was produced by samtools stats (%s+htslib-%s) and can be plotted using plot-bamstats\n", samtools_version(), hts_version());
    if( stats->split_name != NULL ){
        fprintf(to, "# This file contains statistics only for reads with tag: %s=%s\n", stats->info->split_tag, stats->split_name);
    }
    else{
        fprintf(to, "# This file contains statistics for all reads.\n");
    }
    fprintf(to, "# The command line was:  %s",stats->info->argv[0]);
    int i;
    for (i=1; i<stats->info->argc; i++)
        fprintf(to, " %s", stats->info->argv[i]);
    fprintf(to, "\n");
    fprintf(to, "# CHK, Checksum\t[2]Read Names\t[3]Sequences\t[4]Qualities\n");
    fprintf(to, "# CHK, CRC32 of reads which passed filtering followed by addition (32bit overflow)\n");
    fprintf(to, "CHK\t%08x\t%08x\t%08x\n", stats->checksum.names,stats->checksum.reads,stats->checksum.quals);
    fprintf(to, "# Summary Numbers. Use `grep ^SN | cut -f 2-` to extract this part.\n");
    fprintf(to, "SN\traw total sequences:\t%ld\n", (long)(stats->nreads_filtered+stats->nreads_1st+stats->nreads_2nd));  // not counting excluded seqs (and none of the below)
    fprintf(to, "SN\tfiltered sequences:\t%ld\n", (long)stats->nreads_filtered);
    fprintf(to, "SN\tsequences:\t%ld\n", (long)(stats->nreads_1st+stats->nreads_2nd));
    fprintf(to, "SN\tis sorted:\t%d\n", stats->is_sorted ? 1 : 0);
    fprintf(to, "SN\t1st fragments:\t%ld\n", (long)stats->nreads_1st);
    fprintf(to, "SN\tlast fragments:\t%ld\n", (long)stats->nreads_2nd);
    fprintf(to, "SN\treads mapped:\t%ld\n", (long)(stats->nreads_paired_and_mapped+stats->nreads_single_mapped));
    fprintf(to, "SN\treads mapped and paired:\t%ld\t# paired-end technology bit set + both mates mapped\n", (long)stats->nreads_paired_and_mapped);
    fprintf(to, "SN\treads unmapped:\t%ld\n", (long)stats->nreads_unmapped);
    fprintf(to, "SN\treads properly paired:\t%ld\t# proper-pair bit set\n", (long)stats->nreads_properly_paired);
    fprintf(to, "SN\treads paired:\t%ld\t# paired-end technology bit set\n", (long)stats->nreads_paired_tech);
    fprintf(to, "SN\treads duplicated:\t%ld\t# PCR or optical duplicate bit set\n", (long)stats->nreads_dup);
    fprintf(to, "SN\treads MQ0:\t%ld\t# mapped and MQ=0\n", (long)stats->nreads_mq0);
    fprintf(to, "SN\treads QC failed:\t%ld\n", (long)stats->nreads_QCfailed);
    fprintf(to, "SN\tnon-primary alignments:\t%ld\n", (long)stats->nreads_secondary);
    fprintf(to, "SN\ttotal length:\t%ld\t# ignores clipping\n", (long)stats->total_len);
    fprintf(to, "SN\ttotal first fragment length:\t%ld\t# ignores clipping\n", (long)stats->total_len_1st);
    fprintf(to, "SN\ttotal last fragment length:\t%ld\t# ignores clipping\n", (long)stats->total_len_2nd);
    fprintf(to, "SN\tbases mapped:\t%ld\t# ignores clipping\n", (long)stats->nbases_mapped);                 // the length of the whole read goes here, including soft-clips etc.
    fprintf(to, "SN\tbases mapped (cigar):\t%ld\t# more accurate\n", (long)stats->nbases_mapped_cigar);   // only matched and inserted bases are counted here
    fprintf(to, "SN\tbases trimmed:\t%ld\n", (long)stats->nbases_trimmed);
    fprintf(to, "SN\tbases duplicated:\t%ld\n", (long)stats->total_len_dup);
    fprintf(to, "SN\tmismatches:\t%ld\t# from NM fields\n", (long)stats->nmismatches);
    fprintf(to, "SN\terror rate:\t%e\t# mismatches / bases mapped (cigar)\n", stats->nbases_mapped_cigar ? (float)stats->nmismatches/stats->nbases_mapped_cigar : 0);
    float avg_read_length = (stats->nreads_1st+stats->nreads_2nd)?stats->total_len/(stats->nreads_1st+stats->nreads_2nd):0;
    fprintf(to, "SN\taverage length:\t%.0f\n", avg_read_length);
    fprintf(to, "SN\taverage first fragment length:\t%.0f\n", stats->nreads_1st? (float)stats->total_len_1st/stats->nreads_1st:0);
    fprintf(to, "SN\taverage last fragment length:\t%.0f\n", stats->nreads_2nd? (float)stats->total_len_2nd/stats->nreads_2nd:0);
    fprintf(to, "SN\tmaximum length:\t%d\n", stats->max_len);
    fprintf(to, "SN\tmaximum first fragment length:\t%d\n", stats->max_len_1st);
    fprintf(to, "SN\tmaximum last fragment length:\t%d\n", stats->max_len_2nd);
    fprintf(to, "SN\taverage quality:\t%.1f\n", stats->total_len?stats->sum_qual/stats->total_len:0);
    fprintf(to, "SN\tinsert size average:\t%.1f\n", avg_isize);
    fprintf(to, "SN\tinsert size standard deviation:\t%.1f\n", sd_isize);
    fprintf(to, "SN\tinward oriented pairs:\t%ld\n", (long)nisize_inward);
    fprintf(to, "SN\toutward oriented pairs:\t%ld\n", (long)nisize_outward);
    fprintf(to, "SN\tpairs with other orientation:\t%ld\n", (long)nisize_other);
    fprintf(to, "SN\tpairs on different chromosomes:\t%ld\n", (long)stats->nreads_anomalous/2);
    fprintf(to, "SN\tpercentage of properly paired reads (%%):\t%.1f\n", (stats->nreads_1st+stats->nreads_2nd)? (float)(100*stats->nreads_properly_paired)/(stats->nreads_1st+stats->nreads_2nd):0);
    if ( stats->target_count ) {
        fprintf(to, "SN\tbases inside the target:\t%u\n", stats->target_count);
        for (icov=stats->info->cov_threshold+1; icov<stats->ncov; icov++)
            cov_sum += stats->cov[icov];
        fprintf(to, "SN\tpercentage of target genome with coverage > %d (%%):\t%.2f\n", stats->info->cov_threshold, (float)(100*cov_sum)/stats->target_count);
    }

    int ibase,iqual;
    if ( stats->max_len<stats->nbases ) stats->max_len++;
    if ( stats->max_qual+1<stats->nquals ) stats->max_qual++;
    fprintf(to, "# First Fragment Qualities. Use `grep ^FFQ | cut -f 2-` to extract this part.\n");
    fprintf(to, "# Columns correspond to qualities and rows to cycles. First column is the cycle number.\n");
    for (ibase=0; ibase<stats->max_len_1st; ibase++)
    {
        fprintf(to, "FFQ\t%d",ibase+1);
        for (iqual=0; iqual<=stats->max_qual; iqual++)
        {
            fprintf(to, "\t%ld", (long)stats->quals_1st[ibase*stats->nquals+iqual]);
        }
        fprintf(to, "\n");
    }
    fprintf(to, "# Last Fragment Qualities. Use `grep ^LFQ | cut -f 2-` to extract this part.\n");
    fprintf(to, "# Columns correspond to qualities and rows to cycles. First column is the cycle number.\n");
    for (ibase=0; ibase<stats->max_len_2nd; ibase++)
    {
        fprintf(to, "LFQ\t%d",ibase+1);
        for (iqual=0; iqual<=stats->max_qual; iqual++)
        {
            fprintf(to, "\t%ld", (long)stats->quals_2nd[ibase*stats->nquals+iqual]);
        }
        fprintf(to, "\n");
    }
    if ( stats->mpc_buf )
    {
        fprintf(to, "# Mismatches per cycle and quality. Use `grep ^MPC | cut -f 2-` to extract this part.\n");
        fprintf(to, "# Columns correspond to qualities, rows to cycles. First column is the cycle number, second\n");
        fprintf(to, "# is the number of N's and the rest is the number of mismatches\n");
        for (ibase=0; ibase<stats->max_len; ibase++)
        {
            fprintf(to, "MPC\t%d",ibase+1);
            for (iqual=0; iqual<=stats->max_qual; iqual++)
            {
                fprintf(to, "\t%ld", (long)stats->mpc_buf[ibase*stats->nquals+iqual]);
            }
            fprintf(to, "\n");
        }
    }
    fprintf(to, "# GC Content of first fragments. Use `grep ^GCF | cut -f 2-` to extract this part.\n");
    int ibase_prev = 0;
    for (ibase=0; ibase<stats->ngc; ibase++)
    {
        if ( stats->gc_1st[ibase]==stats->gc_1st[ibase_prev] ) continue;
        fprintf(to, "GCF\t%.2f\t%ld\n", (ibase+ibase_prev)*0.5*100./(stats->ngc-1), (long)stats->gc_1st[ibase_prev]);
        ibase_prev = ibase;
    }
    fprintf(to, "# GC Content of last fragments. Use `grep ^GCL | cut -f 2-` to extract this part.\n");
    ibase_prev = 0;
    for (ibase=0; ibase<stats->ngc; ibase++)
    {
        if ( stats->gc_2nd[ibase]==stats->gc_2nd[ibase_prev] ) continue;
        fprintf(to, "GCL\t%.2f\t%ld\n", (ibase+ibase_prev)*0.5*100./(stats->ngc-1), (long)stats->gc_2nd[ibase_prev]);
        ibase_prev = ibase;
    }
    fprintf(to, "# ACGT content per cycle. Use `grep ^GCC | cut -f 2-` to extract this part. The columns are: cycle; A,C,G,T base counts as a percentage of all A/C/G/T bases [%%]; and N and O counts as a percentage of all A/C/G/T bases [%%]\n");
    for (ibase=0; ibase<stats->max_len; ibase++)
    {
        acgtno_count_t *acgtno_count_1st = &(stats->acgtno_cycles_1st[ibase]);
        acgtno_count_t *acgtno_count_2nd = &(stats->acgtno_cycles_2nd[ibase]);
        uint64_t acgt_sum = acgtno_count_1st->a + acgtno_count_1st->c + acgtno_count_1st->g + acgtno_count_1st->t +
                acgtno_count_2nd->a + acgtno_count_2nd->c + acgtno_count_2nd->g + acgtno_count_2nd->t;
        if ( ! acgt_sum ) continue;
        fprintf(to, "GCC\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n", ibase+1,
                100.*(acgtno_count_1st->a + acgtno_count_2nd->a)/acgt_sum,
                100.*(acgtno_count_1st->c + acgtno_count_2nd->c)/acgt_sum,
                100.*(acgtno_count_1st->g + acgtno_count_2nd->g)/acgt_sum,
                100.*(acgtno_count_1st->t + acgtno_count_2nd->t)/acgt_sum,
                100.*(acgtno_count_1st->n + acgtno_count_2nd->n)/acgt_sum,
                100.*(acgtno_count_1st->other + acgtno_count_2nd->other)/acgt_sum);

    }
    fprintf(to, "# ACGT content per cycle for first fragments. Use `grep ^FBC | cut -f 2-` to extract this part. The columns are: cycle; A,C,G,T base counts as a percentage of all A/C/G/T bases [%%]; and N and O counts as a percentage of all A/C/G/T bases [%%]\n");
    for (ibase=0; ibase<stats->max_len; ibase++)
    {
        acgtno_count_t *acgtno_count_1st = &(stats->acgtno_cycles_1st[ibase]);
        uint64_t acgt_sum_1st = acgtno_count_1st->a + acgtno_count_1st->c + acgtno_count_1st->g + acgtno_count_1st->t;

        if ( acgt_sum_1st )
            fprintf(to, "FBC\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n", ibase+1,
                    100.*acgtno_count_1st->a/acgt_sum_1st,
                    100.*acgtno_count_1st->c/acgt_sum_1st,
                    100.*acgtno_count_1st->g/acgt_sum_1st,
                    100.*acgtno_count_1st->t/acgt_sum_1st,
                    100.*acgtno_count_1st->n/acgt_sum_1st,
                    100.*acgtno_count_1st->other/acgt_sum_1st);

    }
    fprintf(to, "# ACGT content per cycle for last fragments. Use `grep ^LBC | cut -f 2-` to extract this part. The columns are: cycle; A,C,G,T base counts as a percentage of all A/C/G/T bases [%%]; and N and O counts as a percentage of all A/C/G/T bases [%%]\n");
    for (ibase=0; ibase<stats->max_len; ibase++)
    {
        acgtno_count_t *acgtno_count_2nd = &(stats->acgtno_cycles_2nd[ibase]);
        uint64_t acgt_sum_2nd = acgtno_count_2nd->a + acgtno_count_2nd->c + acgtno_count_2nd->g + acgtno_count_2nd->t;

        if ( acgt_sum_2nd )
            fprintf(to, "LBC\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n", ibase+1,
                    100.*acgtno_count_2nd->a/acgt_sum_2nd,
                    100.*acgtno_count_2nd->c/acgt_sum_2nd,
                    100.*acgtno_count_2nd->g/acgt_sum_2nd,
                    100.*acgtno_count_2nd->t/acgt_sum_2nd,
                    100.*acgtno_count_2nd->n/acgt_sum_2nd,
                    100.*acgtno_count_2nd->other/acgt_sum_2nd);

    }
    fprintf(to, "# Insert sizes. Use `grep ^IS | cut -f 2-` to extract this part. The columns are: insert size, pairs total, inward oriented pairs, outward oriented pairs, other pairs\n");
    for (isize=0; isize<ibulk; isize++) {
        long in = (long)(stats->isize->inward(stats->isize->data, isize));
        long out = (long)(stats->isize->outward(stats->isize->data, isize));
        long other = (long)(stats->isize->other(stats->isize->data, isize));
        if (!sparse || in + out + other > 0) {
            fprintf(to, "IS\t%d\t%ld\t%ld\t%ld\t%ld\n", isize,  in+out+other,
                in , out, other);
        }
    }

    fprintf(to, "# Read lengths. Use `grep ^RL | cut -f 2-` to extract this part. The columns are: read length, count\n");
    int ilen;
    for (ilen=0; ilen<stats->max_len; ilen++)
    {
        if ( stats->read_lengths[ilen+1]>0 )
            fprintf(to, "RL\t%d\t%ld\n", ilen+1, (long)stats->read_lengths[ilen+1]);
    }

    fprintf(to, "# Read lengths - first fragments. Use `grep ^FRL | cut -f 2-` to extract this part. The columns are: read length, count\n");
    for (ilen=0; ilen<stats->max_len_1st; ilen++)
    {
        if ( stats->read_lengths_1st[ilen+1]>0 )
            fprintf(to, "FRL\t%d\t%ld\n", ilen+1, (long)stats->read_lengths_1st[ilen+1]);
    }

    fprintf(to, "# Read lengths - last fragments. Use `grep ^LRL | cut -f 2-` to extract this part. The columns are: read length, count\n");
    for (ilen=0; ilen<stats->max_len_2nd; ilen++)
    {
        if ( stats->read_lengths_2nd[ilen+1]>0 )
            fprintf(to, "LRL\t%d\t%ld\n", ilen+1, (long)stats->read_lengths_2nd[ilen+1]);
    }

    fprintf(to, "# Indel distribution. Use `grep ^ID | cut -f 2-` to extract this part. The columns are: length, number of insertions, number of deletions\n");

    for (ilen=0; ilen<stats->nindels; ilen++)
    {
        if ( stats->insertions[ilen]>0 || stats->deletions[ilen]>0 )
            fprintf(to, "ID\t%d\t%ld\t%ld\n", ilen+1, (long)stats->insertions[ilen], (long)stats->deletions[ilen]);
    }

    fprintf(to, "# Indels per cycle. Use `grep ^IC | cut -f 2-` to extract this part. The columns are: cycle, number of insertions (fwd), .. (rev) , number of deletions (fwd), .. (rev)\n");
    for (ilen=0; ilen<=stats->nbases; ilen++)
    {
        // For deletions we print the index of the cycle before the deleted base (1-based) and for insertions
        //  the index of the cycle of the first inserted base (also 1-based)
        if ( stats->ins_cycles_1st[ilen]>0 || stats->ins_cycles_2nd[ilen]>0 || stats->del_cycles_1st[ilen]>0 || stats->del_cycles_2nd[ilen]>0 )
            fprintf(to, "IC\t%d\t%ld\t%ld\t%ld\t%ld\n", ilen+1, (long)stats->ins_cycles_1st[ilen], (long)stats->ins_cycles_2nd[ilen], (long)stats->del_cycles_1st[ilen], (long)stats->del_cycles_2nd[ilen]);
    }

    fprintf(to, "# Coverage distribution. Use `grep ^COV | cut -f 2-` to extract this part.\n");
    if  ( stats->cov[0] )
        fprintf(to, "COV\t[<%d]\t%d\t%ld\n",stats->info->cov_min,stats->info->cov_min-1, (long)stats->cov[0]);
    for (icov=1; icov<stats->ncov-1; icov++)
        if ( stats->cov[icov] )
            fprintf(to, "COV\t[%d-%d]\t%d\t%ld\n",stats->info->cov_min + (icov-1)*stats->info->cov_step, stats->info->cov_min + icov*stats->info->cov_step-1,stats->info->cov_min + icov*stats->info->cov_step-1, (long)stats->cov[icov]);
    if ( stats->cov[stats->ncov-1] )
        fprintf(to, "COV\t[%d<]\t%d\t%ld\n",stats->info->cov_min + (stats->ncov-2)*stats->info->cov_step-1,stats->info->cov_min + (stats->ncov-2)*stats->info->cov_step-1, (long)stats->cov[stats->ncov-1]);

    // Calculate average GC content, then sort by GC and depth
    fprintf(to, "# GC-depth. Use `grep ^GCD | cut -f 2-` to extract this part. The columns are: GC%%, unique sequence percentiles, 10th, 25th, 50th, 75th and 90th depth percentile\n");
    uint32_t igcd;
    for (igcd=0; igcd<stats->igcd; igcd++)
    {
        if ( stats->info->fai )
            stats->gcd[igcd].gc = rint(100. * stats->gcd[igcd].gc);
        else
            if ( stats->gcd[igcd].depth )
                stats->gcd[igcd].gc = rint(100. * stats->gcd[igcd].gc / stats->gcd[igcd].depth);
    }
    qsort(stats->gcd, stats->igcd+1, sizeof(gc_depth_t), gcd_cmp);
    igcd = 0;
    while ( igcd < stats->igcd )
    {
        // Calculate percentiles (10,25,50,75,90th) for the current GC content and print
        uint32_t nbins=0, itmp=igcd;
        float gc = stats->gcd[igcd].gc;
        while ( itmp<stats->igcd && fabs(stats->gcd[itmp].gc-gc)<0.1 )
        {
            nbins++;
            itmp++;
        }
        fprintf(to, "GCD\t%.1f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n", gc, (igcd+nbins+1)*100./(stats->igcd+1),
                gcd_percentile(&(stats->gcd[igcd]),nbins,10) *avg_read_length/stats->info->gcd_bin_size,
                gcd_percentile(&(stats->gcd[igcd]),nbins,25) *avg_read_length/stats->info->gcd_bin_size,
                gcd_percentile(&(stats->gcd[igcd]),nbins,50) *avg_read_length/stats->info->gcd_bin_size,
                gcd_percentile(&(stats->gcd[igcd]),nbins,75) *avg_read_length/stats->info->gcd_bin_size,
                gcd_percentile(&(stats->gcd[igcd]),nbins,90) *avg_read_length/stats->info->gcd_bin_size
              );
        igcd += nbins;
    }
}

void init_regions(stats_t *stats, const char *file)
{
    FILE *fp = fopen(file,"r");
    if ( !fp ) error("%s: %s\n",file,strerror(errno));

    kstring_t line = { 0, 0, NULL };
    int warned = 0, r, p, new_p;
    int prev_tid=-1, prev_pos=-1;
    while (line.l = 0, kgetline(&line, (kgets_func *)fgets, fp) >= 0)
    {
        if ( line.s[0] == '#' ) continue;

        int i = 0;
        while ( i<line.l && !isspace(line.s[i]) ) i++;
        if ( i>=line.l ) error("Could not parse the file: %s [%s]\n", file, line.s);
        line.s[i] = '\0';

        int tid = bam_name2id(stats->info->sam_header, line.s);
        if ( tid < 0 )
        {
            if ( !warned )
                fprintf(stderr,"Warning: Some sequences not present in the BAM, e.g. \"%s\". This message is printed only once.\n", line.s);
            warned = 1;
            continue;
        }

        if ( tid >= stats->nregions )
        {
            stats->regions = realloc(stats->regions,sizeof(regions_t)*(stats->nregions+100));
            int j;
            for (j=stats->nregions; j<stats->nregions+100; j++)
            {
                stats->regions[j].npos = stats->regions[j].mpos = stats->regions[j].cpos = 0;
                stats->regions[j].pos = NULL;
            }
            stats->nregions += 100;
        }
        int npos = stats->regions[tid].npos;
        if ( npos >= stats->regions[tid].mpos )
        {
            stats->regions[tid].mpos += 1000;
            stats->regions[tid].pos = realloc(stats->regions[tid].pos,sizeof(pos_t)*stats->regions[tid].mpos);
        }

        if ( (sscanf(&line.s[i+1],"%u %u",&stats->regions[tid].pos[npos].from,&stats->regions[tid].pos[npos].to))!=2 ) error("Could not parse the region [%s]\n", &line.s[i+1]);
        if ( prev_tid==-1 || prev_tid!=tid )
        {
            prev_tid = tid;
            prev_pos = stats->regions[tid].pos[npos].from;
        }
        if ( prev_pos>stats->regions[tid].pos[npos].from )
            error("The positions are not in chromosomal order (%s:%d comes after %d)\n", line.s,stats->regions[tid].pos[npos].from,prev_pos);
        stats->regions[tid].npos++;
        if ( stats->regions[tid].npos > stats->nchunks )
            stats->nchunks = stats->regions[tid].npos;
    }
    free(line.s);
    if ( !stats->regions ) error("Unable to map the -t sequences to the BAM sequences.\n");
    fclose(fp);

    // sort region intervals and remove duplicates
    for (r = 0; r < stats->nregions; r++) {
        regions_t *reg = &stats->regions[r];
        if ( reg->npos > 1 ) {
            qsort(reg->pos, reg->npos, sizeof(pos_t), regions_lt);
            for (new_p = 0, p = 1; p < reg->npos; p++) {
                if ( reg->pos[new_p].to < reg->pos[p].from )
                    reg->pos[++new_p] = reg->pos[p];
                else if ( reg->pos[new_p].to < reg->pos[p].to )
                    reg->pos[new_p].to = reg->pos[p].to;
            }
            reg->npos = ++new_p;
        }
        for (p = 0; p < reg->npos; p++)
            stats->target_count += (reg->pos[p].to - reg->pos[p].from + 1);
    }

    stats->chunks = calloc(stats->nchunks, sizeof(pos_t));
}

void destroy_regions(stats_t *stats)
{
    int i;
    for (i=0; i<stats->nregions; i++)
    {
        if ( !stats->regions[i].mpos ) continue;
        free(stats->regions[i].pos);
    }
    if ( stats->regions ) free(stats->regions);
    if ( stats->chunks ) free(stats->chunks);
}

void reset_regions(stats_t *stats)
{
    int i;
    for (i=0; i<stats->nregions; i++)
        stats->regions[i].cpos = 0;
}

int is_in_regions(bam1_t *bam_line, stats_t *stats)
{
    if ( !stats->regions ) return 1;

    if ( bam_line->core.tid >= stats->nregions || bam_line->core.tid<0 ) return 0;
    if ( !stats->is_sorted ) error("The BAM must be sorted in order for -t to work.\n");

    regions_t *reg = &stats->regions[bam_line->core.tid];
    if ( reg->cpos==reg->npos ) return 0;       // done for this chr

    // Find a matching interval or skip this read. No splicing of reads is done, no indels or soft clips considered,
    //  even small overlap is enough to include the read in the stats.
    int i = reg->cpos;
    while ( i<reg->npos && reg->pos[i].to<=bam_line->core.pos ) i++;
    if ( i>=reg->npos ) { reg->cpos = reg->npos; return 0; }
    int64_t endpos = bam_endpos(bam_line);
    if ( endpos < reg->pos[i].from ) return 0;

    //found a read overlapping a region
    reg->cpos = i;
    stats->reg_from = reg->pos[i].from;
    stats->reg_to   = reg->pos[i].to;

    //now find all the overlapping chunks
    stats->nchunks = 0;
    while (i < reg->npos) {
        if (bam_line->core.pos < reg->pos[i].to && endpos >= reg->pos[i].from) {
            stats->chunks[stats->nchunks].from = MAX(bam_line->core.pos+1, reg->pos[i].from);
            stats->chunks[stats->nchunks].to = MIN(endpos, reg->pos[i].to);
            stats->nchunks++;
        }
        i++;
    }

    return 1;
}

int replicate_regions(stats_t *stats, hts_itr_multi_t *iter) {
    if ( !stats || !iter)
        return 1;

    int i, j, tid;
    stats->nregions = iter->n_reg;
    stats->regions = calloc(stats->nregions, sizeof(regions_t));
    stats->chunks = calloc(stats->nchunks, sizeof(pos_t));
    if ( !stats->regions || !stats->chunks )
        return 1;

    for (i = 0; i < iter->n_reg; i++) {
        tid = iter->reg_list[i].tid;
        if ( tid < 0 )
            continue;

        if ( tid >= stats->nregions ) {
            regions_t *tmp = realloc(stats->regions, (tid+10) * sizeof(regions_t));
            if ( !tmp )
                return 1;
            stats->regions = tmp;
            memset(stats->regions + stats->nregions, 0,
                   (tid+10-stats->nregions) * sizeof(regions_t));
            stats->nregions = tid+10;
        }

        stats->regions[tid].mpos = stats->regions[tid].npos = iter->reg_list[i].count;
        stats->regions[tid].pos = calloc(stats->regions[tid].mpos, sizeof(pos_t));
        if ( !stats->regions[tid].pos )
            return 1;

        for (j = 0; j < stats->regions[tid].npos; j++) {
            stats->regions[tid].pos[j].from = iter->reg_list[i].intervals[j].beg+1;
            stats->regions[tid].pos[j].to = iter->reg_list[i].intervals[j].end;

            stats->target_count += (stats->regions[tid].pos[j].to - stats->regions[tid].pos[j].from + 1);
        }
    }

    return 0;
}

void init_group_id(stats_t *stats, const char *id)
{
#if 0
    if ( !stats->sam_header->dict )
        stats->sam_header->dict = sam_header_parse2(stats->sam_header->text);
    void *iter = stats->sam_header->dict;
    const char *key, *val;
    int n = 0;
    stats->rg_hash = khash_str2int_init();
    while ( (iter = sam_header2key_val(iter, "RG","ID","SM", &key, &val)) )
    {
        if ( !strcmp(id,key) || (val && !strcmp(id,val)) )
        {
            khiter_t k = kh_get(kh_rg, stats->rg_hash, key);
            if ( k != kh_end(stats->rg_hash) )
                fprintf(stderr, "[init_group_id] The group ID not unique: \"%s\"\n", key);
            int ret;
            k = kh_put(kh_rg, stats->rg_hash, key, &ret);
            kh_value(stats->rg_hash, k) = val;
            n++;
        }
    }
    if ( !n )
        error("The sample or read group \"%s\" not present.\n", id);
#else
    fprintf(stderr, "Samtools-htslib: init_group_id() header parsing not yet implemented\n");
    abort();
#endif
}


static void error(const char *format, ...)
{
    if ( !format )
    {
        printf("About: The program collects statistics from BAM files. The output can be visualized using plot-bamstats.\n");
        printf("Usage: samtools stats [OPTIONS] file.bam\n");
        printf("       samtools stats [OPTIONS] file.bam chr:from-to\n");
        printf("Options:\n");
        printf("    -c, --coverage <int>,<int>,<int>    Coverage distribution min,max,step [1,1000,1]\n");
        printf("    -d, --remove-dups                   Exclude from statistics reads marked as duplicates\n");
        printf("    -f, --required-flag  <str|int>      Required flag, 0 for unset. See also `samtools flags` [0]\n");
        printf("    -F, --filtering-flag <str|int>      Filtering flag, 0 for unset. See also `samtools flags` [0]\n");
        printf("        --GC-depth <float>              the size of GC-depth bins (decreasing bin size increases memory requirement) [2e4]\n");
        printf("    -h, --help                          This help message\n");
        printf("    -i, --insert-size <int>             Maximum insert size [8000]\n");
        printf("    -I, --id <string>                   Include only listed read group or sample name\n");
        printf("    -l, --read-length <int>             Include in the statistics only reads with the given read length []\n");
        printf("    -m, --most-inserts <float>          Report only the main part of inserts [0.99]\n");
        printf("    -P, --split-prefix <str>            Path or string prefix for filepaths output by -S (default is input filename)\n");
        printf("    -q, --trim-quality <int>            The BWA trimming parameter [0]\n");
        printf("    -r, --ref-seq <file>                Reference sequence (required for GC-depth and mismatches-per-cycle calculation).\n");
        printf("    -s, --sam                           Ignored (input format is auto-detected).\n");
        printf("    -S, --split <tag>                   Also write statistics to separate files split by tagged field.\n");
        printf("    -t, --target-regions <file>         Do stats in these regions only. Tab-delimited file chr,from,to, 1-based, inclusive.\n");
        printf("    -x, --sparse                        Suppress outputting IS rows where there are no insertions.\n");
        printf("    -p, --remove-overlaps               Remove overlaps of paired-end reads from coverage and base count computations.\n");
        printf("    -g, --cov-threshold                 Only bases with coverage above this value will be included in the target percentage computation.\n");
        sam_global_opt_help(stdout, "-.--.@");
        printf("\n");
    }
    else
    {
        va_list ap;
        va_start(ap, format);
        vfprintf(stderr, format, ap);
        va_end(ap);
    }
    exit(1);
}

void cleanup_stats_info(stats_info_t* info){
    if (info->fai) fai_destroy(info->fai);
    sam_close(info->sam);
    free(info);
}

void cleanup_stats(stats_t* stats)
{
    free(stats->cov_rbuf.buffer); free(stats->cov);
    free(stats->quals_1st); free(stats->quals_2nd);
    free(stats->gc_1st); free(stats->gc_2nd);
    stats->isize->isize_free(stats->isize->data);
    free(stats->isize);
    free(stats->gcd);
    free(stats->rseq_buf);
    free(stats->mpc_buf);
    free(stats->acgtno_cycles_1st);
    free(stats->acgtno_cycles_2nd);
    free(stats->read_lengths);
    free(stats->read_lengths_1st);
    free(stats->read_lengths_2nd);
    free(stats->insertions);
    free(stats->deletions);
    free(stats->ins_cycles_1st);
    free(stats->ins_cycles_2nd);
    free(stats->del_cycles_1st);
    free(stats->del_cycles_2nd);
    destroy_regions(stats);
    if ( stats->rg_hash ) khash_str2int_destroy(stats->rg_hash);
    free(stats->split_name);
    free(stats);
}

void output_split_stats(khash_t(c2stats) *split_hash, char* bam_fname, int sparse)
{
    int i = 0;
    kstring_t output_filename = { 0, 0, NULL };
    stats_t *curr_stats = NULL;
    for(i = kh_begin(split_hash); i != kh_end(split_hash); ++i){
        if(!kh_exist(split_hash, i)) continue;
        curr_stats = kh_value(split_hash, i);
        round_buffer_flush(curr_stats, -1);

        output_filename.l = 0;
        if (curr_stats->info->split_prefix)
            kputs(curr_stats->info->split_prefix, &output_filename);
        else
            kputs(bam_fname, &output_filename);
        kputc('_', &output_filename);
        kputs(curr_stats->split_name, &output_filename);
        kputs(".bamstat", &output_filename);

        FILE *to = fopen(output_filename.s, "w");
        if(to == NULL){
            error("Could not open '%s' for writing.\n", output_filename.s);
        }
        output_stats(to, curr_stats, sparse);
        fclose(to);
    }

    free(output_filename.s);
}

void destroy_split_stats(khash_t(c2stats) *split_hash)
{
    int i = 0;
    stats_t *curr_stats = NULL;
    for(i = kh_begin(split_hash); i != kh_end(split_hash); ++i){
        if(!kh_exist(split_hash, i)) continue;
        curr_stats = kh_value(split_hash, i);
        cleanup_stats(curr_stats);
    }
    kh_destroy(c2stats, split_hash);
}

stats_info_t* stats_info_init(int argc, char *argv[])
{
    stats_info_t* info = calloc(1, sizeof(stats_info_t));
    info->nisize = 8000;
    info->isize_main_bulk = 0.99;   // There are always outliers at the far end
    info->gcd_bin_size = 20e3;
    info->cov_min  = 1;
    info->cov_max  = 1000;
    info->cov_step = 1;
    info->filter_readlen = -1;
    info->argc = argc;
    info->argv = argv;
    info->remove_overlaps = 0;
    info->cov_threshold = 0;

    return info;
}

int init_stat_info_fname(stats_info_t* info, const char* bam_fname, const htsFormat* in_fmt)
{
    // .. bam
    samFile* sam;
    if ((sam = sam_open_format(bam_fname, "r", in_fmt)) == 0) {
        print_error_errno("stats", "failed to open \"%s\"", bam_fname);
        return 1;
    }
    info->sam = sam;
    info->sam_header = sam_hdr_read(sam);
    if (info->sam_header == NULL) {
        print_error("stats", "failed to read header for \"%s\"", bam_fname);
        return 1;
    }
    return 0;
}

stats_t* stats_init()
{
    stats_t *stats = calloc(1,sizeof(stats_t));
    stats->ngc    = 200;
    stats->nquals = 256;
    stats->nbases = 300;
    stats->rseq_pos     = -1;
    stats->tid = stats->gcd_pos = -1;
    stats->igcd = 0;
    stats->is_sorted = 1;
    stats->nindels = stats->nbases;
    stats->split_name = NULL;
    stats->nchunks = 0;
    stats->pair_count = 0;
    stats->last_pair_tid = -2;
    stats->last_read_flush = 0;
    stats->target_count = 0;

    return stats;
}

static void init_stat_structs(stats_t* stats, stats_info_t* info, const char* group_id, const char* targets)
{
    // Give stats_t a pointer to the info struct
    // This saves us having to pass the stats_info_t to every function
    stats->info = info;

    // Init structures
    //  .. coverage bins and round buffer
    if ( info->cov_step > info->cov_max - info->cov_min + 1 )
    {
        info->cov_step = info->cov_max - info->cov_min;
        if ( info->cov_step <= 0 )
            info->cov_step = 1;
    }
    stats->ncov = 3 + (info->cov_max-info->cov_min) / info->cov_step;
    info->cov_max = info->cov_min + ((info->cov_max-info->cov_min)/info->cov_step +1)*info->cov_step - 1;
    stats->cov = calloc(sizeof(uint64_t),stats->ncov);
    stats->cov_rbuf.size = stats->nbases*5;
    stats->cov_rbuf.buffer = calloc(sizeof(int32_t),stats->cov_rbuf.size);

    if ( group_id ) init_group_id(stats, group_id);
    // .. arrays
    stats->quals_1st      = calloc(stats->nquals*stats->nbases,sizeof(uint64_t));
    stats->quals_2nd      = calloc(stats->nquals*stats->nbases,sizeof(uint64_t));
    stats->gc_1st         = calloc(stats->ngc,sizeof(uint64_t));
    stats->gc_2nd         = calloc(stats->ngc,sizeof(uint64_t));
    stats->isize          = init_isize_t(info->nisize ?info->nisize+1 :0);
    stats->gcd            = calloc(stats->ngcd,sizeof(gc_depth_t));
    stats->mpc_buf        = info->fai ? calloc(stats->nquals*stats->nbases,sizeof(uint64_t)) : NULL;
    stats->acgtno_cycles_1st  = calloc(stats->nbases,sizeof(acgtno_count_t));
    stats->acgtno_cycles_2nd  = calloc(stats->nbases,sizeof(acgtno_count_t));
    stats->read_lengths   = calloc(stats->nbases,sizeof(uint64_t));
    stats->read_lengths_1st   = calloc(stats->nbases,sizeof(uint64_t));
    stats->read_lengths_2nd   = calloc(stats->nbases,sizeof(uint64_t));
    stats->insertions     = calloc(stats->nbases,sizeof(uint64_t));
    stats->deletions      = calloc(stats->nbases,sizeof(uint64_t));
    stats->ins_cycles_1st = calloc(stats->nbases+1,sizeof(uint64_t));
    stats->ins_cycles_2nd = calloc(stats->nbases+1,sizeof(uint64_t));
    stats->del_cycles_1st = calloc(stats->nbases+1,sizeof(uint64_t));
    stats->del_cycles_2nd = calloc(stats->nbases+1,sizeof(uint64_t));
    realloc_rseq_buffer(stats);
    if ( targets )
        init_regions(stats, targets);
}

static stats_t* get_curr_split_stats(bam1_t* bam_line, khash_t(c2stats)* split_hash, stats_info_t* info, char* targets)
{
    stats_t *curr_stats = NULL;
    const uint8_t *tag_val = bam_aux_get(bam_line, info->split_tag);
    if(tag_val == 0){
        error("Tag '%s' not found in bam_line.\n", info->split_tag);
    }
    char* split_name = strdup(bam_aux2Z(tag_val));

    // New stats object, under split
    khiter_t k = kh_get(c2stats, split_hash, split_name);
    if(k == kh_end(split_hash)){
        curr_stats = stats_init(); // mallocs new instance
        init_stat_structs(curr_stats, info, NULL, targets);
        curr_stats->split_name = split_name;

        // Record index in hash
        int ret = 0;
        khiter_t iter = kh_put(c2stats, split_hash, split_name, &ret);
        if( ret < 0 ){
            error("Failed to insert key '%s' into split_hash", split_name);
        }
        kh_val(split_hash, iter) = curr_stats; // store pointer to stats
    }
    else{
        curr_stats = kh_value(split_hash, k);
        free(split_name); // don't need to hold on to this if it wasn't new
    }
    return curr_stats;
}

int main_stats(int argc, char *argv[])
{
    char *targets = NULL;
    char *bam_fname = NULL;
    char *group_id = NULL;
    int sparse = 0;
    sam_global_args ga = SAM_GLOBAL_ARGS_INIT;

    stats_info_t *info = stats_info_init(argc, argv);

    static const struct option loptions[] =
    {
        SAM_OPT_GLOBAL_OPTIONS('-', 0, '-', '-', 0, '@'),
        {"help", no_argument, NULL, 'h'},
        {"remove-dups", no_argument, NULL, 'd'},
        {"sam", no_argument, NULL, 's'},
        {"ref-seq", required_argument, NULL, 'r'},
        {"coverage", required_argument, NULL, 'c'},
        {"read-length", required_argument, NULL, 'l'},
        {"insert-size", required_argument, NULL, 'i'},
        {"most-inserts", required_argument, NULL, 'm'},
        {"trim-quality", required_argument, NULL, 'q'},
        {"target-regions", required_argument, NULL, 't'},
        {"required-flag", required_argument, NULL, 'f'},
        {"filtering-flag", required_argument, NULL, 'F'},
        {"id", required_argument, NULL, 'I'},
        {"GC-depth", required_argument, NULL, 1},
        {"sparse", no_argument, NULL, 'x'},
        {"split", required_argument, NULL, 'S'},
        {"split-prefix", required_argument, NULL, 'P'},
        {"remove-overlaps", no_argument, NULL, 'p'},
        {"cov-threshold", required_argument, NULL, 'g'},
        {NULL, 0, NULL, 0}
    };
    int opt;

    while ( (opt=getopt_long(argc,argv,"?hdsxpr:c:l:i:t:m:q:f:F:g:I:1:S:P:@:",loptions,NULL))>0 )
    {
        switch (opt)
        {
            case 'f': info->flag_require = bam_str2flag(optarg); break;
            case 'F': info->flag_filter |= bam_str2flag(optarg); break;
            case 'd': info->flag_filter |= BAM_FDUP; break;
            case 's': break;
            case 'r': info->fai = fai_load(optarg);
                      if (info->fai==NULL)
                          error("Could not load faidx: %s\n", optarg);
                      break;
            case  1 : info->gcd_bin_size = atof(optarg); break;
            case 'c': if ( sscanf(optarg,"%d,%d,%d",&info->cov_min,&info->cov_max,&info->cov_step)!= 3 )
                          error("Unable to parse -c %s\n", optarg);
                      break;
            case 'l': info->filter_readlen = atoi(optarg); break;
            case 'i': info->nisize = atoi(optarg); break;
            case 'm': info->isize_main_bulk = atof(optarg); break;
            case 'q': info->trim_qual = atoi(optarg); break;
            case 't': targets = optarg; break;
            case 'I': group_id = optarg; break;
            case 'x': sparse = 1; break;
            case 'S': info->split_tag = optarg; break;
            case 'P': info->split_prefix = optarg; break;
            case 'p': info->remove_overlaps = 1; break;
            case 'g': info->cov_threshold = atoi(optarg);
                      if ( info->cov_threshold < 0 || info->cov_threshold == INT_MAX )
                          error("Unsupported value for coverage threshold %d\n", info->cov_threshold);
                      break;
            case '?':
            case 'h': error(NULL);
            default:
                if (parse_sam_global_opt(opt, optarg, loptions, &ga) != 0)
                    error("Unknown argument: %s\n", optarg);
                break;
        }
    }
    if ( optind<argc )
        bam_fname = argv[optind++];

    if ( !bam_fname )
    {
        if ( isatty(STDIN_FILENO) )
            error(NULL);
        bam_fname = "-";
    }

    if (init_stat_info_fname(info, bam_fname, &ga.in)) {
        free(info);
        return 1;
    }
    if (ga.nthreads > 0)
        hts_set_threads(info->sam, ga.nthreads);

    stats_t *all_stats = stats_init();
    stats_t *curr_stats = NULL;
    init_stat_structs(all_stats, info, group_id, targets);
    // Init
    // .. hash
    khash_t(c2stats)* split_hash = kh_init(c2stats);

    khash_t(qn2pair)* read_pairs = kh_init(qn2pair);

    // Collect statistics
    bam1_t *bam_line = bam_init1();
    if ( optind<argc )
    {
        int filter = 1;
        // Prepare the region hash table for the multi-region iterator
        void *region_hash = bed_hash_regions(NULL, argv, optind, argc, &filter);
        if (region_hash) {

            // Collect stats in selected regions only
            hts_idx_t *bam_idx = sam_index_load(info->sam,bam_fname);
            if (bam_idx) {

                int regcount = 0;
                hts_reglist_t *reglist = bed_reglist(region_hash, ALL, &regcount);
                if (reglist) {

                    hts_itr_multi_t *iter = sam_itr_regions(bam_idx, info->sam_header, reglist, regcount);
                    if (iter) {

                        if (!targets) {
                            all_stats->nchunks = argc-optind;
                            if ( replicate_regions(all_stats, iter) )
                                fprintf(stderr, "Replications of the regions failed.");
                        }

                        if ( all_stats->nregions && all_stats->regions ) {
                            while (sam_itr_multi_next(info->sam, iter, bam_line) >= 0) {
                               if (info->split_tag) {
                                   curr_stats = get_curr_split_stats(bam_line, split_hash, info, targets);
                                   collect_stats(bam_line, curr_stats, read_pairs);
                               }
                               collect_stats(bam_line, all_stats, read_pairs);
                            }
                        }

                        hts_itr_multi_destroy(iter);
                    } else {
                       fprintf(stderr, "Creation of the region iterator failed.");
                       hts_reglist_free(reglist, regcount);
                    }
                } else {
                    fprintf(stderr, "Creation of the region list failed.");
                }

                hts_idx_destroy(bam_idx);
            } else {
                fprintf(stderr, "Random alignment retrieval only works for indexed BAM files.\n");
            }

            bed_destroy(region_hash);
        } else {
            fprintf(stderr, "Creation of the region hash table failed.\n");
        }
    }
    else
    {
        if ( info->cov_threshold > 0 && !targets ) {
            fprintf(stderr, "Coverage percentage calcuation requires a list of target regions\n");
            goto cleanup;
        }

        // Stream through the entire BAM ignoring off-target regions if -t is given
        int ret;
        while ((ret = sam_read1(info->sam, info->sam_header, bam_line)) >= 0) {
            if (info->split_tag) {
                curr_stats = get_curr_split_stats(bam_line, split_hash, info, targets);
                collect_stats(bam_line, curr_stats, read_pairs);
            }
            collect_stats(bam_line, all_stats, read_pairs);
        }

        if (ret < -1) {
            fprintf(stderr, "Failure while decoding file\n");
            return 1;
        }
    }

    round_buffer_flush(all_stats, -1);
    output_stats(stdout, all_stats, sparse);
    if (info->split_tag)
        output_split_stats(split_hash, bam_fname, sparse);

cleanup:
    bam_destroy1(bam_line);
    bam_hdr_destroy(info->sam_header);
    sam_global_args_free(&ga);

    cleanup_stats(all_stats);
    cleanup_stats_info(info);
    destroy_split_stats(split_hash);
    cleanup_overlaps(read_pairs, INT_MAX);

    return 0;
}
