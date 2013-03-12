/* 
    Author: petr.danecek@sanger
    gcc -Wall -Winline -g -O2 -I ~/git/samtools bamcheck.c -o bamcheck -lm -lz -L ~/git/samtools -lbam -lpthread

    Assumptions, approximations and other issues:
        - GC-depth graph does not split reads, the starting position determines which bin is incremented.
            There are small overlaps between bins (max readlen-1). However, the bins are big (20k).
        - coverage distribution ignores softclips and deletions
        - some stats require sorted BAMs
        - GC content graph can have an untidy, step-like pattern when BAM contains multiple read lengths.
        - 'bases mapped' (stats->nbases_mapped) is calculated from read lengths given by BAM (core.l_qseq)
        - With the -t option, the whole reads are used. Except for the number of mapped bases (cigar)
            counts, no splicing is done, no indels or soft clips are considered, even small overlap is
            good enough to include the read in the stats.

*/

#define BAMCHECK_VERSION "2012-09-04"

#define _ISOC99_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <getopt.h>
#include <errno.h>
#include <assert.h>
#include "faidx.h"
#include "khash.h"
#include "sam.h"
#include "sam_header.h"
#include "razf.h"

#define BWA_MIN_RDLEN 35
#define IS_PAIRED(bam) ((bam)->core.flag&BAM_FPAIRED && !((bam)->core.flag&BAM_FUNMAP) && !((bam)->core.flag&BAM_FMUNMAP))
#define IS_UNMAPPED(bam) ((bam)->core.flag&BAM_FUNMAP)
#define IS_REVERSE(bam) ((bam)->core.flag&BAM_FREVERSE)
#define IS_MATE_REVERSE(bam) ((bam)->core.flag&BAM_FMREVERSE)
#define IS_READ1(bam) ((bam)->core.flag&BAM_FREAD1)
#define IS_READ2(bam) ((bam)->core.flag&BAM_FREAD2)
#define IS_DUP(bam) ((bam)->core.flag&BAM_FDUP)

typedef struct 
{
    int32_t line_len, line_blen;
    int64_t len;
    uint64_t offset;
} 
faidx1_t;
KHASH_MAP_INIT_STR(kh_faidx, faidx1_t)
KHASH_MAP_INIT_STR(kh_bam_tid, int)
KHASH_MAP_INIT_STR(kh_rg, const char *)
struct __faidx_t {
    RAZF *rz;
    int n, m;
    char **name;
    khash_t(kh_faidx) *hash;
};

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
    // Parameters
    int trim_qual;      // bwa trim quality

    // Dimensions of the quality histogram holder (quals_1st,quals_2nd), GC content holder (gc_1st,gc_2nd),
    //  insert size histogram holder
    int nquals;         // The number of quality bins 
    int nbases;         // The maximum sequence length the allocated array can hold
    int nisize;         // The maximum insert size that the allocated array can hold
    int ngc;            // The size of gc_1st and gc_2nd
    int nindels;        // The maximum indel length for indel distribution

    // Arrays for the histogram data
    uint64_t *quals_1st, *quals_2nd;
    uint64_t *gc_1st, *gc_2nd;
    uint64_t *isize_inward, *isize_outward, *isize_other;
    uint64_t *acgt_cycles;
    uint64_t *read_lengths;
    uint64_t *insertions, *deletions;
    uint64_t *ins_cycles_1st, *ins_cycles_2nd, *del_cycles_1st, *del_cycles_2nd;

    // The extremes encountered
    int max_len;            // Maximum read length
    int max_qual;           // Maximum quality
    float isize_main_bulk;  // There are always some unrealistically big insert sizes, report only the main part
    int is_sorted;

    // Summary numbers
    uint64_t total_len;
    uint64_t total_len_dup;
    uint64_t nreads_1st;
    uint64_t nreads_2nd;
    uint64_t nreads_filtered;
    uint64_t nreads_dup;
    uint64_t nreads_unmapped;
    uint64_t nreads_unpaired;
    uint64_t nreads_paired;
    uint64_t nreads_anomalous;
    uint64_t nreads_mq0;
    uint64_t nbases_mapped;
    uint64_t nbases_mapped_cigar;
    uint64_t nbases_trimmed;  // bwa trimmed bases
    uint64_t nmismatches;
    uint64_t nreads_QCfailed, nreads_secondary;

    // GC-depth related data
    uint32_t ngcd, igcd;        // The maximum number of GC depth bins and index of the current bin
    gc_depth_t *gcd;            // The GC-depth bins holder
    int gcd_bin_size;           // The size of GC-depth bin
    uint32_t gcd_ref_size;      // The approximate size of the genome
    int32_t tid, gcd_pos;       // Position of the current bin
    int32_t pos;                // Position of the last read

    // Coverage distribution related data
    int ncov;                       // The number of coverage bins
    uint64_t *cov;                  // The coverage frequencies
    int cov_min,cov_max,cov_step;   // Minimum, maximum coverage and size of the coverage bins
    round_buffer_t cov_rbuf;        // Pileup round buffer

    // Mismatches by read cycle 
    uint8_t *rseq_buf;              // A buffer for reference sequence to check the mismatches against
    int mrseq_buf;                  // The size of the buffer
    int32_t rseq_pos;               // The coordinate of the first base in the buffer
    int32_t nrseq_buf;              // The used part of the buffer
    uint64_t *mpc_buf;              // Mismatches per cycle

    // Filters
    int filter_readlen;

    // Target regions
    int nregions, reg_from,reg_to;
    regions_t *regions;

    // Auxiliary data
    int flag_require, flag_filter;
    double sum_qual;                // For calculating average quality value 
    samfile_t *sam;             
    khash_t(kh_rg) *rg_hash;        // Read groups to include, the array is null-terminated
    faidx_t *fai;                   // Reference sequence for GC-depth graph
    int argc;                       // Command line arguments to be printed on the output
    char **argv;
}
stats_t;

void error(const char *format, ...);
void bam_init_header_hash(bam_header_t *header);
int is_in_regions(bam1_t *bam_line, stats_t *stats);


// Coverage distribution methods
inline int coverage_idx(int min, int max, int n, int step, int depth)
{
    if ( depth < min )
        return 0;

    if ( depth > max )
        return n-1;

    return 1 + (depth - min) / step;
}

inline int round_buffer_lidx2ridx(int offset, int size, int64_t refpos, int64_t pos)
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
            idp = coverage_idx(stats->cov_min,stats->cov_max,stats->ncov,stats->cov_step,stats->cov_rbuf.buffer[ibuf]);
            stats->cov[idp]++;
            stats->cov_rbuf.buffer[ibuf] = 0;
        }
        ifrom = 0;
    }
    for (ibuf=ifrom; ibuf<=ito; ibuf++)
    {
        if ( !stats->cov_rbuf.buffer[ibuf] )
            continue;
        idp = coverage_idx(stats->cov_min,stats->cov_max,stats->ncov,stats->cov_step,stats->cov_rbuf.buffer[ibuf]);
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
        // Conversion from uint32_t to MIDNSHP
        //  0123456
        //  MIDNSHP
        int cig  = bam1_cigar(bam_line)[icig] & BAM_CIGAR_MASK;
        int ncig = bam1_cigar(bam_line)[icig] >> BAM_CIGAR_SHIFT;

        if ( cig==1 )
        {
            int idx = is_fwd ? icycle : read_len-icycle-ncig;
            if ( idx<0 ) 
                error("FIXME: read_len=%d vs icycle=%d\n", read_len,icycle);
            if ( idx >= stats->nbases || idx<0 ) error("FIXME: %d vs %d, %s:%d %s\n", idx,stats->nbases, stats->sam->header->target_name[bam_line->core.tid],bam_line->core.pos+1,bam1_qname(bam_line));
            if ( is_1st ) 
                stats->ins_cycles_1st[idx]++;
            else
                stats->ins_cycles_2nd[idx]++;
            icycle += ncig;
            if ( ncig<=stats->nindels )
                stats->insertions[ncig-1]++;
            continue;
        }
        if ( cig==2 )
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
        if ( cig!=3 && cig!=5 )
            icycle += ncig;
    }
}

void count_mismatches_per_cycle(stats_t *stats,bam1_t *bam_line) 
{
    int is_fwd = IS_REVERSE(bam_line) ? 0 : 1;
    int icig,iread=0,icycle=0;
    int iref = bam_line->core.pos - stats->rseq_pos;
    int read_len   = bam_line->core.l_qseq;
    uint8_t *read  = bam1_seq(bam_line);
    uint8_t *quals = bam1_qual(bam_line);
    uint64_t *mpc_buf = stats->mpc_buf;
    for (icig=0; icig<bam_line->core.n_cigar; icig++) 
    {
        // Conversion from uint32_t to MIDNSHP
        //  0123456
        //  MIDNSHP
        int cig  = bam1_cigar(bam_line)[icig] & BAM_CIGAR_MASK;
        int ncig = bam1_cigar(bam_line)[icig] >> BAM_CIGAR_SHIFT;
        if ( cig==1 )
        {
            iread  += ncig;
            icycle += ncig;
            continue;
        }
        if ( cig==2 )
        {
            iref += ncig;
            continue;
        }
        if ( cig==4 )
        {
            icycle += ncig;
            // Soft-clips are present in the sequence, but the position of the read marks a start of non-clipped sequence
            //   iref += ncig;
            iread  += ncig;
            continue;
        }
        if ( cig==5 )
        {
            icycle += ncig;
            continue;
        }
        // Ignore H and N CIGARs. The letter are inserted e.g. by TopHat and often require very large
        //  chunk of refseq in memory. Not very frequent and not noticable in the stats.
        if ( cig==3 || cig==5 ) continue;
        if ( cig!=0 )
            error("TODO: cigar %d, %s:%d %s\n", cig,stats->sam->header->target_name[bam_line->core.tid],bam_line->core.pos+1,bam1_qname(bam_line));
       
        if ( ncig+iref > stats->nrseq_buf )
            error("FIXME: %d+%d > %d, %s, %s:%d\n",ncig,iref,stats->nrseq_buf, bam1_qname(bam_line),stats->sam->header->target_name[bam_line->core.tid],bam_line->core.pos+1);

        int im;
        for (im=0; im<ncig; im++)
        {
            uint8_t cread = bam1_seqi(read,iread);
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
                    error("TODO: quality too high %d>=%d (%s %d %s)\n", qual,stats->nquals, stats->sam->header->target_name[bam_line->core.tid],bam_line->core.pos+1,bam1_qname(bam_line));

                int idx = is_fwd ? icycle : read_len-icycle-1;
                if ( idx>stats->max_len )
                    error("mpc: %d>%d\n",idx,stats->max_len);

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

void read_ref_seq(stats_t *stats,int32_t tid,int32_t pos)
{
    khash_t(kh_faidx) *h;
    khiter_t iter;
    faidx1_t val;
    char *chr, c;
    faidx_t *fai = stats->fai;

    h = fai->hash;
    chr = stats->sam->header->target_name[tid];

    // ID of the sequence name
    iter = kh_get(kh_faidx, h, chr);
    if (iter == kh_end(h)) 
        error("No such reference sequence [%s]?\n", chr);
    val = kh_value(h, iter);

    // Check the boundaries
    if (pos >= val.len)
        error("Was the bam file mapped with the reference sequence supplied?"
              " A read mapped beyond the end of the chromosome (%s:%d, chromosome length %d).\n", chr,pos,val.len);
    int size = stats->mrseq_buf;
    // The buffer extends beyond the chromosome end. Later the rest will be filled with N's.
    if (size+pos > val.len) size = val.len-pos;

    // Position the razf reader
    razf_seek(fai->rz, val.offset + pos / val.line_blen * val.line_len + pos % val.line_blen, SEEK_SET);

    uint8_t *ptr = stats->rseq_buf;
    int nread = 0;
    while ( nread<size && razf_read(fai->rz,&c,1) && !fai->rz->z_err )
    {
        if ( !isgraph(c) )
            continue;

        // Conversion between uint8_t coding and ACGT
        //      -12-4---8-------
        //      =ACMGRSVTWYHKDBN
        if ( c=='A' || c=='a' )
            *ptr = 1;
        else if ( c=='C' || c=='c' )
            *ptr = 2;
        else if ( c=='G' || c=='g' )
            *ptr = 4;
        else if ( c=='T' || c=='t' )
            *ptr = 8;
        else
            *ptr = 0;
        ptr++;
        nread++;
    }
    if ( nread < stats->mrseq_buf )
    {
        memset(ptr,0, stats->mrseq_buf - nread);
        nread = stats->mrseq_buf;
    }
    stats->nrseq_buf = nread;
    stats->rseq_pos  = pos;
    stats->tid       = tid;
}

float fai_gc_content(stats_t *stats, int pos, int len)
{
    uint32_t gc,count,c;
    int i = pos - stats->rseq_pos, ito = i + len;
    assert( i>=0 && ito<=stats->nrseq_buf );

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
    if ( stats->gcd_bin_size > n ) n = stats->gcd_bin_size;
    if ( stats->mrseq_buf<n )
    {
        stats->rseq_buf = realloc(stats->rseq_buf,sizeof(uint8_t)*n);
        stats->mrseq_buf = n;
    }
}

void realloc_gcd_buffer(stats_t *stats, int seq_len)
{
    if ( seq_len >= stats->gcd_bin_size )
        error("The --GC-depth bin size (%d) is set too low for the read length %d\n", stats->gcd_bin_size, seq_len);

    int n = 1 + stats->gcd_ref_size / (stats->gcd_bin_size - seq_len);
    if ( n <= stats->igcd )
        error("The --GC-depth bin size is too small or reference genome too big; please decrease the bin size or increase the reference length\n");
        
    if ( n > stats->ngcd )
    {
        stats->gcd = realloc(stats->gcd, n*sizeof(gc_depth_t));
        if ( !stats->gcd )
            error("Could not realloc GCD buffer, too many chromosomes or the genome too long?? [%u %u]\n", stats->ngcd,n);
        memset(&(stats->gcd[stats->ngcd]),0,(n-stats->ngcd)*sizeof(gc_depth_t));
        stats->ngcd = n;
    }

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

    stats->acgt_cycles = realloc(stats->acgt_cycles, n*4*sizeof(uint64_t));
    if ( !stats->acgt_cycles )
        error("Could not realloc buffers, the sequence too long: %d (%ld)\n", seq_len,n*4*sizeof(uint64_t));
    memset(stats->acgt_cycles + stats->nbases*4, 0, (n-stats->nbases)*4*sizeof(uint64_t));

    stats->read_lengths = realloc(stats->read_lengths, n*sizeof(uint64_t));
    if ( !stats->read_lengths )
        error("Could not realloc buffers, the sequence too long: %d (%ld)\n", seq_len,n*sizeof(uint64_t));
    memset(stats->read_lengths + stats->nbases, 0, (n-stats->nbases)*sizeof(uint64_t));

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

void collect_stats(bam1_t *bam_line, stats_t *stats)
{
    if ( stats->rg_hash )
    {
        const uint8_t *rg = bam_aux_get(bam_line, "RG");
        if ( !rg ) return; 
        khiter_t k = kh_get(kh_rg, stats->rg_hash, (const char*)(rg + 1));
        if ( k == kh_end(stats->rg_hash) ) return;
    }
    if ( stats->flag_require && (bam_line->core.flag & stats->flag_require)!=stats->flag_require )
    {
        stats->nreads_filtered++;
        return;
    }
    if ( stats->flag_filter && (bam_line->core.flag & stats->flag_filter) )
    {
        stats->nreads_filtered++;
        return;
    }
    if ( !is_in_regions(bam_line,stats) )
        return;
    if ( stats->filter_readlen!=-1 && bam_line->core.l_qseq!=stats->filter_readlen ) 
        return;

    if ( bam_line->core.flag & BAM_FQCFAIL ) stats->nreads_QCfailed++;
    if ( bam_line->core.flag & BAM_FSECONDARY ) stats->nreads_secondary++;

    int seq_len = bam_line->core.l_qseq;
    if ( !seq_len ) return;

    if ( seq_len >= stats->nbases )
        realloc_buffers(stats,seq_len);
    if ( stats->max_len<seq_len )
        stats->max_len = seq_len;

    stats->read_lengths[seq_len]++;

    // Count GC and ACGT per cycle
    uint8_t base, *seq  = bam1_seq(bam_line);
    int gc_count  = 0;
    int i;
    int reverse = IS_REVERSE(bam_line);
    for (i=0; i<seq_len; i++)
    {
        // Conversion from uint8_t coding to ACGT
        //      -12-4---8-------
        //      =ACMGRSVTWYHKDBN
        //       01 2   3
        base = bam1_seqi(seq,i);
        base /= 2;
        if ( base==1 || base==2 ) gc_count++;
        else if ( base>2 ) base=3;
        if ( 4*(reverse ? seq_len-i-1 : i) + base >= stats->nbases*4 ) 
            error("FIXME: acgt_cycles\n");
        stats->acgt_cycles[ 4*(reverse ? seq_len-i-1 : i) + base ]++;
    }
    int gc_idx_min = gc_count*(stats->ngc-1)/seq_len;
    int gc_idx_max = (gc_count+1)*(stats->ngc-1)/seq_len;
    if ( gc_idx_max >= stats->ngc ) gc_idx_max = stats->ngc - 1;

    // Determine which array (1st or 2nd read) will these stats go to,
    //  trim low quality bases from end the same way BWA does, 
    //  fill GC histogram
    uint64_t *quals;
    uint8_t *bam_quals = bam1_qual(bam_line);
    if ( bam_line->core.flag&BAM_FREAD2 )
    {
        quals  = stats->quals_2nd;
        stats->nreads_2nd++;
        for (i=gc_idx_min; i<gc_idx_max; i++)
            stats->gc_2nd[i]++;
    }
    else
    {
        quals = stats->quals_1st;
        stats->nreads_1st++;
        for (i=gc_idx_min; i<gc_idx_max; i++)
            stats->gc_1st[i]++;
    }
    if ( stats->trim_qual>0 ) 
        stats->nbases_trimmed += bwa_trim_read(stats->trim_qual, bam_quals, seq_len, reverse);

    // Quality histogram and average quality
    for (i=0; i<seq_len; i++)
    {
        uint8_t qual = bam_quals[ reverse ? seq_len-i-1 : i];
        if ( qual>=stats->nquals )
            error("TODO: quality too high %d>=%d (%s %d %s)\n", qual,stats->nquals,stats->sam->header->target_name[bam_line->core.tid],bam_line->core.pos+1,bam1_qname(bam_line));
        if ( qual>stats->max_qual )
            stats->max_qual = qual;

        quals[ i*stats->nquals+qual ]++;
        stats->sum_qual += qual;
    }

    // Look at the flags and increment appropriate counters (mapped, paired, etc)
    if ( IS_UNMAPPED(bam_line) )
        stats->nreads_unmapped++;
    else
    {
        if ( !bam_line->core.qual )
            stats->nreads_mq0++;

        count_indels(stats,bam_line);

        if ( !IS_PAIRED(bam_line) )
            stats->nreads_unpaired++;
        else 
        {
            stats->nreads_paired++;

            if ( bam_line->core.tid!=bam_line->core.mtid )
                stats->nreads_anomalous++;

            // The insert size is tricky, because for long inserts the libraries are
            // prepared differently and the pairs point in other direction. BWA does
            // not set the paired flag for them.  Similar thing is true also for 454
            // reads. Mates mapped to different chromosomes have isize==0.
            int32_t isize = bam_line->core.isize;
            if ( isize<0 ) isize = -isize;
            if ( isize >= stats->nisize )
                isize = stats->nisize-1;
            if ( isize>0 || bam_line->core.tid==bam_line->core.mtid )
            {
                int pos_fst = bam_line->core.mpos - bam_line->core.pos;
                int is_fst  = IS_READ1(bam_line) ? 1 : -1;
                int is_fwd  = IS_REVERSE(bam_line) ? -1 : 1;
                int is_mfwd = IS_MATE_REVERSE(bam_line) ? -1 : 1;

                if ( is_fwd*is_mfwd>0 )
                    stats->isize_other[isize]++;
                else if ( is_fst*pos_fst>0 )
                {
                    if ( is_fst*is_fwd>0 )
                        stats->isize_inward[isize]++;
                    else
                        stats->isize_outward[isize]++;
                }
                else if ( is_fst*pos_fst<0 )
                {
                    if ( is_fst*is_fwd>0 )
                        stats->isize_outward[isize]++;
                    else
                        stats->isize_inward[isize]++;
                }
            }
        }

        // Number of mismatches 
        uint8_t *nm = bam_aux_get(bam_line,"NM");
        if (nm) 
            stats->nmismatches += bam_aux2i(nm);

        // Number of mapped bases from cigar 
        // Conversion from uint32_t to MIDNSHP
        //  012-4--
        //  MIDNSHP
        if ( bam_line->core.n_cigar == 0) 
            error("FIXME: mapped read with no cigar?\n");
        int readlen=seq_len;
        if ( stats->regions )
        {
            // Count only on-target bases
            int iref = bam_line->core.pos + 1;
            for (i=0; i<bam_line->core.n_cigar; i++)
            {
                int cig  = bam1_cigar(bam_line)[i]&BAM_CIGAR_MASK;
                int ncig = bam1_cigar(bam_line)[i]>>BAM_CIGAR_SHIFT;
                if ( cig==2 ) readlen += ncig;
                else if ( cig==0 ) 
                {
                    if ( iref < stats->reg_from ) ncig -= stats->reg_from-iref;
                    else if ( iref+ncig-1 > stats->reg_to ) ncig -= iref+ncig-1 - stats->reg_to;
                    if ( ncig<0 ) ncig = 0;
                    stats->nbases_mapped_cigar += ncig;
                    iref += bam1_cigar(bam_line)[i]>>BAM_CIGAR_SHIFT;
                }
                else if ( cig==1 )
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
                if ( (bam1_cigar(bam_line)[i]&BAM_CIGAR_MASK)==0 || (bam1_cigar(bam_line)[i]&BAM_CIGAR_MASK)==1 )
                    stats->nbases_mapped_cigar += bam1_cigar(bam_line)[i]>>BAM_CIGAR_SHIFT;
                if ( (bam1_cigar(bam_line)[i]&BAM_CIGAR_MASK)==2 )
                    readlen += bam1_cigar(bam_line)[i]>>BAM_CIGAR_SHIFT;
            }
        }
        stats->nbases_mapped += seq_len;

        if ( stats->tid==bam_line->core.tid && bam_line->core.pos<stats->pos )
            stats->is_sorted = 0;
        stats->pos = bam_line->core.pos;

        if ( stats->is_sorted )
        {
            if ( stats->tid==-1 || stats->tid!=bam_line->core.tid )
                round_buffer_flush(stats,-1);

            // Mismatches per cycle and GC-depth graph. For simplicity, reads overlapping GCD bins
            //  are not splitted which results in up to seq_len-1 overlaps. The default bin size is
            //  20kbp, so the effect is negligible.
            if ( stats->fai )
            {
                int inc_ref = 0, inc_gcd = 0;
                // First pass or new chromosome
                if ( stats->rseq_pos==-1 || stats->tid != bam_line->core.tid ) { inc_ref=1; inc_gcd=1; }
                // Read goes beyond the end of the rseq buffer
                else if ( stats->rseq_pos+stats->nrseq_buf < bam_line->core.pos+readlen ) { inc_ref=1; inc_gcd=1; }
                // Read overlaps the next gcd bin
                else if ( stats->gcd_pos+stats->gcd_bin_size < bam_line->core.pos+readlen ) 
                { 
                    inc_gcd = 1;
                    if ( stats->rseq_pos+stats->nrseq_buf < bam_line->core.pos+stats->gcd_bin_size ) inc_ref = 1;
                }
                if ( inc_gcd )
                {
                    stats->igcd++;
                    if ( stats->igcd >= stats->ngcd )
                        realloc_gcd_buffer(stats, readlen);
                    if ( inc_ref )
                        read_ref_seq(stats,bam_line->core.tid,bam_line->core.pos);
                    stats->gcd_pos = bam_line->core.pos;
                    stats->gcd[ stats->igcd ].gc = fai_gc_content(stats, stats->gcd_pos, stats->gcd_bin_size);
                }

                count_mismatches_per_cycle(stats,bam_line);
            }
            // No reference and first pass, new chromosome or sequence going beyond the end of the gcd bin
            else if ( stats->gcd_pos==-1 || stats->tid != bam_line->core.tid || bam_line->core.pos - stats->gcd_pos > stats->gcd_bin_size )
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
            if ( !stats->fai )
                stats->gcd[ stats->igcd ].gc += (float) gc_count / seq_len;

            // Coverage distribution graph
            round_buffer_flush(stats,bam_line->core.pos);
            round_buffer_insert_read(&(stats->cov_rbuf),bam_line->core.pos,bam_line->core.pos+seq_len-1);
        }
    }

    stats->total_len += seq_len;
    if ( IS_DUP(bam_line) )
    {
        stats->total_len_dup += seq_len;
        stats->nreads_dup++;
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

void output_stats(stats_t *stats)
{
    // Calculate average insert size and standard deviation (from the main bulk data only)
    int isize, ibulk=0;
    uint64_t nisize=0, nisize_inward=0, nisize_outward=0, nisize_other=0;
    for (isize=0; isize<stats->nisize; isize++)
    {
        // Each pair was counted twice
        stats->isize_inward[isize] *= 0.5;
        stats->isize_outward[isize] *= 0.5;
        stats->isize_other[isize] *= 0.5;

        nisize_inward += stats->isize_inward[isize];
        nisize_outward += stats->isize_outward[isize];
        nisize_other += stats->isize_other[isize];
        nisize += stats->isize_inward[isize] + stats->isize_outward[isize] + stats->isize_other[isize];
    }

    double bulk=0, avg_isize=0, sd_isize=0;
    for (isize=0; isize<stats->nisize; isize++)
    {
        bulk += stats->isize_inward[isize] + stats->isize_outward[isize] + stats->isize_other[isize];
        avg_isize += isize * (stats->isize_inward[isize] + stats->isize_outward[isize] + stats->isize_other[isize]);

        if ( bulk/nisize > stats->isize_main_bulk )
        {
            ibulk  = isize+1;
            nisize = bulk;
            break;
        }
    }
    avg_isize /= nisize ? nisize : 1;
    for (isize=1; isize<ibulk; isize++)
        sd_isize += (stats->isize_inward[isize] + stats->isize_outward[isize] + stats->isize_other[isize]) * (isize-avg_isize)*(isize-avg_isize) / nisize;
    sd_isize = sqrt(sd_isize);


    printf("# This file was produced by bamcheck (%s)\n",BAMCHECK_VERSION);
    printf("# The command line was:  %s",stats->argv[0]);
    int i;
    for (i=1; i<stats->argc; i++)
        printf(" %s",stats->argv[i]);
    printf("\n");
    printf("# Summary Numbers. Use `grep ^SN | cut -f 2-` to extract this part.\n");
    printf("SN\traw total sequences:\t%ld\n", (long)(stats->nreads_filtered+stats->nreads_1st+stats->nreads_2nd));
    printf("SN\tfiltered sequences:\t%ld\n", (long)stats->nreads_filtered);
    printf("SN\tsequences:\t%ld\n", (long)(stats->nreads_1st+stats->nreads_2nd));
    printf("SN\tis paired:\t%d\n", stats->nreads_1st&&stats->nreads_2nd ? 1 : 0);
    printf("SN\tis sorted:\t%d\n", stats->is_sorted ? 1 : 0);
    printf("SN\t1st fragments:\t%ld\n", (long)stats->nreads_1st);
    printf("SN\tlast fragments:\t%ld\n", (long)stats->nreads_2nd);
    printf("SN\treads mapped:\t%ld\n", (long)(stats->nreads_paired+stats->nreads_unpaired));
    printf("SN\treads unmapped:\t%ld\n", (long)stats->nreads_unmapped);
    printf("SN\treads unpaired:\t%ld\n", (long)stats->nreads_unpaired);
    printf("SN\treads paired:\t%ld\n", (long)stats->nreads_paired);
    printf("SN\treads duplicated:\t%ld\n", (long)stats->nreads_dup);
    printf("SN\treads MQ0:\t%ld\n", (long)stats->nreads_mq0);
    printf("SN\treads QC failed:\t%ld\n", (long)stats->nreads_QCfailed);
    printf("SN\tnon-primary alignments:\t%ld\n", (long)stats->nreads_secondary);
    printf("SN\ttotal length:\t%ld\n", (long)stats->total_len);
    printf("SN\tbases mapped:\t%ld\n", (long)stats->nbases_mapped);
    printf("SN\tbases mapped (cigar):\t%ld\n", (long)stats->nbases_mapped_cigar);
    printf("SN\tbases trimmed:\t%ld\n", (long)stats->nbases_trimmed);
    printf("SN\tbases duplicated:\t%ld\n", (long)stats->total_len_dup);
    printf("SN\tmismatches:\t%ld\n", (long)stats->nmismatches);
    printf("SN\terror rate:\t%e\n", (float)stats->nmismatches/stats->nbases_mapped_cigar);
    float avg_read_length = (stats->nreads_1st+stats->nreads_2nd)?stats->total_len/(stats->nreads_1st+stats->nreads_2nd):0;
    printf("SN\taverage length:\t%.0f\n", avg_read_length);
    printf("SN\tmaximum length:\t%d\n", stats->max_len);
    printf("SN\taverage quality:\t%.1f\n", stats->total_len?stats->sum_qual/stats->total_len:0);
    printf("SN\tinsert size average:\t%.1f\n", avg_isize);
    printf("SN\tinsert size standard deviation:\t%.1f\n", sd_isize);
    printf("SN\tinward oriented pairs:\t%ld\n", (long)nisize_inward);
    printf("SN\toutward oriented pairs:\t%ld\n", (long)nisize_outward);
    printf("SN\tpairs with other orientation:\t%ld\n", (long)nisize_other);
    printf("SN\tpairs on different chromosomes:\t%ld\n", (long)stats->nreads_anomalous/2);

    int ibase,iqual;
    if ( stats->max_len<stats->nbases ) stats->max_len++;
    if ( stats->max_qual+1<stats->nquals ) stats->max_qual++;
    printf("# First Fragment Qualitites. Use `grep ^FFQ | cut -f 2-` to extract this part.\n");
    printf("# Columns correspond to qualities and rows to cycles. First column is the cycle number.\n");
    for (ibase=0; ibase<stats->max_len; ibase++)
    {
        printf("FFQ\t%d",ibase+1);
        for (iqual=0; iqual<=stats->max_qual; iqual++)
        {
            printf("\t%ld", (long)stats->quals_1st[ibase*stats->nquals+iqual]);
        }
        printf("\n");
    }
    printf("# Last Fragment Qualitites. Use `grep ^LFQ | cut -f 2-` to extract this part.\n");
    printf("# Columns correspond to qualities and rows to cycles. First column is the cycle number.\n");
    for (ibase=0; ibase<stats->max_len; ibase++)
    {
        printf("LFQ\t%d",ibase+1);
        for (iqual=0; iqual<=stats->max_qual; iqual++)
        {
            printf("\t%ld", (long)stats->quals_2nd[ibase*stats->nquals+iqual]);
        }
        printf("\n");
    }
    if ( stats->mpc_buf )
    {
        printf("# Mismatches per cycle and quality. Use `grep ^MPC | cut -f 2-` to extract this part.\n");
        printf("# Columns correspond to qualities, rows to cycles. First column is the cycle number, second\n");
        printf("# is the number of N's and the rest is the number of mismatches\n");
        for (ibase=0; ibase<stats->max_len; ibase++)
        {
            printf("MPC\t%d",ibase+1);
            for (iqual=0; iqual<=stats->max_qual; iqual++)
            {
                printf("\t%ld", (long)stats->mpc_buf[ibase*stats->nquals+iqual]);
            }
            printf("\n");
        }
    }
    printf("# GC Content of first fragments. Use `grep ^GCF | cut -f 2-` to extract this part.\n");
    int ibase_prev = 0;
    for (ibase=0; ibase<stats->ngc; ibase++)
    {
        if ( stats->gc_1st[ibase]==stats->gc_1st[ibase_prev] ) continue;
        printf("GCF\t%.2f\t%ld\n", (ibase+ibase_prev)*0.5*100./(stats->ngc-1), (long)stats->gc_1st[ibase_prev]);
        ibase_prev = ibase;
    }
    printf("# GC Content of last fragments. Use `grep ^GCL | cut -f 2-` to extract this part.\n");
    ibase_prev = 0;
    for (ibase=0; ibase<stats->ngc; ibase++)
    {
        if ( stats->gc_2nd[ibase]==stats->gc_2nd[ibase_prev] ) continue;
        printf("GCL\t%.2f\t%ld\n", (ibase+ibase_prev)*0.5*100./(stats->ngc-1), (long)stats->gc_2nd[ibase_prev]);
        ibase_prev = ibase;
    }
    printf("# ACGT content per cycle. Use `grep ^GCC | cut -f 2-` to extract this part. The columns are: cycle, and A,C,G,T counts [%%]\n");
    for (ibase=0; ibase<stats->max_len; ibase++)
    {
        uint64_t *ptr = &(stats->acgt_cycles[ibase*4]);
        uint64_t  sum = ptr[0]+ptr[1]+ptr[2]+ptr[3];
        if ( ! sum ) continue;
        printf("GCC\t%d\t%.2f\t%.2f\t%.2f\t%.2f\n", ibase,100.*ptr[0]/sum,100.*ptr[1]/sum,100.*ptr[2]/sum,100.*ptr[3]/sum);
    }
    printf("# Insert sizes. Use `grep ^IS | cut -f 2-` to extract this part. The columns are: pairs total, inward oriented pairs, outward oriented pairs, other pairs\n");
    for (isize=0; isize<ibulk; isize++)
        printf("IS\t%d\t%ld\t%ld\t%ld\t%ld\n", isize, (long)(stats->isize_inward[isize]+stats->isize_outward[isize]+stats->isize_other[isize]),
            (long)stats->isize_inward[isize], (long)stats->isize_outward[isize], (long)stats->isize_other[isize]);

    printf("# Read lengths. Use `grep ^RL | cut -f 2-` to extract this part. The columns are: read length, count\n");
    int ilen;
    for (ilen=0; ilen<stats->max_len; ilen++)
    {
        if ( stats->read_lengths[ilen]>0 )
            printf("RL\t%d\t%ld\n", ilen, (long)stats->read_lengths[ilen]);
    }

    printf("# Indel distribution. Use `grep ^ID | cut -f 2-` to extract this part. The columns are: length, number of insertions, number of deletions\n");
    for (ilen=0; ilen<stats->nindels; ilen++)
    {
        if ( stats->insertions[ilen]>0 || stats->deletions[ilen]>0 )
            printf("ID\t%d\t%ld\t%ld\n", ilen+1, (long)stats->insertions[ilen], (long)stats->deletions[ilen]);
    }

    printf("# Indels per cycle. Use `grep ^IC | cut -f 2-` to extract this part. The columns are: cycle, number of insertions (fwd), .. (rev) , number of deletions (fwd), .. (rev)\n");
    for (ilen=0; ilen<=stats->nbases; ilen++)
    {
        // For deletions we print the index of the cycle before the deleted base (1-based) and for insertions
        //  the index of the cycle of the first inserted base (also 1-based)
        if ( stats->ins_cycles_1st[ilen]>0 || stats->ins_cycles_2nd[ilen]>0 || stats->del_cycles_1st[ilen]>0 || stats->del_cycles_2nd[ilen]>0 )
            printf("IC\t%d\t%ld\t%ld\t%ld\t%ld\n", ilen+1, (long)stats->ins_cycles_1st[ilen], (long)stats->ins_cycles_2nd[ilen], (long)stats->del_cycles_1st[ilen], (long)stats->del_cycles_2nd[ilen]);
    }

    printf("# Coverage distribution. Use `grep ^COV | cut -f 2-` to extract this part.\n");
    if  ( stats->cov[0] )
        printf("COV\t[<%d]\t%d\t%ld\n",stats->cov_min,stats->cov_min-1, (long)stats->cov[0]);
    int icov;
    for (icov=1; icov<stats->ncov-1; icov++)
        if ( stats->cov[icov] )
            printf("COV\t[%d-%d]\t%d\t%ld\n",stats->cov_min + (icov-1)*stats->cov_step, stats->cov_min + icov*stats->cov_step-1,stats->cov_min + icov*stats->cov_step-1, (long)stats->cov[icov]);
    if ( stats->cov[stats->ncov-1] )
        printf("COV\t[%d<]\t%d\t%ld\n",stats->cov_min + (stats->ncov-2)*stats->cov_step-1,stats->cov_min + (stats->ncov-2)*stats->cov_step-1, (long)stats->cov[stats->ncov-1]);

    // Calculate average GC content, then sort by GC and depth
    printf("# GC-depth. Use `grep ^GCD | cut -f 2-` to extract this part. The columns are: GC%%, unique sequence percentiles, 10th, 25th, 50th, 75th and 90th depth percentile\n");
    uint32_t igcd;
    for (igcd=0; igcd<stats->igcd; igcd++)
    {
        if ( stats->fai )
            stats->gcd[igcd].gc = round(100. * stats->gcd[igcd].gc);
        else
            if ( stats->gcd[igcd].depth ) 
                stats->gcd[igcd].gc = round(100. * stats->gcd[igcd].gc / stats->gcd[igcd].depth);
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
        printf("GCD\t%.1f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n", gc, (igcd+nbins+1)*100./(stats->igcd+1),
                gcd_percentile(&(stats->gcd[igcd]),nbins,10) *avg_read_length/stats->gcd_bin_size,
                gcd_percentile(&(stats->gcd[igcd]),nbins,25) *avg_read_length/stats->gcd_bin_size, 
                gcd_percentile(&(stats->gcd[igcd]),nbins,50) *avg_read_length/stats->gcd_bin_size, 
                gcd_percentile(&(stats->gcd[igcd]),nbins,75) *avg_read_length/stats->gcd_bin_size, 
                gcd_percentile(&(stats->gcd[igcd]),nbins,90) *avg_read_length/stats->gcd_bin_size 
              );
        igcd += nbins;
    }
}

size_t mygetline(char **line, size_t *n, FILE *fp)
{
    if (line == NULL || n == NULL || fp == NULL)
    {
        errno = EINVAL;
        return -1;
    }
    if (*n==0 || !*line)
    {
        *line = NULL;
        *n = 0;
    }

    size_t nread=0;
    int c;
    while ((c=getc(fp))!= EOF && c!='\n')
    {
        if ( ++nread>=*n )
        {
            *n += 255;
            *line = realloc(*line, sizeof(char)*(*n));
        }
        (*line)[nread-1] = c;
    }
    if ( nread>=*n )
    {
        *n += 255;
        *line = realloc(*line, sizeof(char)*(*n));
    }
    (*line)[nread] = 0;
    return nread>0 ? nread : -1;

}

void init_regions(stats_t *stats, char *file)
{
    khiter_t iter;
    khash_t(kh_bam_tid) *header_hash;

    bam_init_header_hash(stats->sam->header);
    header_hash = (khash_t(kh_bam_tid)*)stats->sam->header->hash;

    FILE *fp = fopen(file,"r");
    if ( !fp ) error("%s: %s\n",file,strerror(errno));

    char *line = NULL;
    size_t len = 0;
    ssize_t nread;
    int warned = 0;
    int prev_tid=-1, prev_pos=-1;
    while ((nread = mygetline(&line, &len, fp)) != -1) 
    {
        if ( line[0] == '#' ) continue;

        int i = 0;
        while ( i<nread && !isspace(line[i]) ) i++;
        if ( i>=nread ) error("Could not parse the file: %s [%s]\n", file,line);
        line[i] = 0;

        iter = kh_get(kh_bam_tid, header_hash, line);
        int tid = kh_val(header_hash, iter);
        if ( iter == kh_end(header_hash) )
        {
            if ( !warned )
                fprintf(stderr,"Warning: Some sequences not present in the BAM, e.g. \"%s\". This message is printed only once.\n", line);
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

        if ( (sscanf(line+i+1,"%d %d",&stats->regions[tid].pos[npos].from,&stats->regions[tid].pos[npos].to))!=2 ) error("Could not parse the region [%s]\n");
        if ( prev_tid==-1 || prev_tid!=tid )
        {
            prev_tid = tid;
            prev_pos = stats->regions[tid].pos[npos].from;
        }
        if ( prev_pos>stats->regions[tid].pos[npos].from )
            error("The positions are not in chromosomal order (%s:%d comes after %d)\n", line,stats->regions[tid].pos[npos].from,prev_pos);
        stats->regions[tid].npos++;
    }
    if (line) free(line);
    if ( !stats->regions ) error("Unable to map the -t sequences to the BAM sequences.\n");
    fclose(fp);
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
}

static int fetch_read(const bam1_t *bam_line, void *data)
{
    collect_stats((bam1_t*)bam_line,(stats_t*)data);
    return 1;
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
    if ( bam_line->core.pos + bam_line->core.l_qseq + 1 < reg->pos[i].from ) return 0;
    reg->cpos = i;
    stats->reg_from = reg->pos[i].from;
    stats->reg_to   = reg->pos[i].to;

    return 1;
}

void init_group_id(stats_t *stats, char *id)
{
    if ( !stats->sam->header->dict )
        stats->sam->header->dict = sam_header_parse2(stats->sam->header->text);
    void *iter = stats->sam->header->dict;
    const char *key, *val;
    int n = 0;
    stats->rg_hash = kh_init(kh_rg);
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
}


void error(const char *format, ...)
{
    if ( !format )
    {
        printf("Version: %s\n", BAMCHECK_VERSION);
        printf("About: The program collects statistics from BAM files. The output can be visualized using plot-bamcheck.\n");
        printf("Usage: bamcheck [OPTIONS] file.bam\n");
        printf("       bamcheck [OPTIONS] file.bam chr:from-to\n");
        printf("Options:\n");
        printf("    -c, --coverage <int>,<int>,<int>    Coverage distribution min,max,step [1,1000,1]\n");
        printf("    -d, --remove-dups                   Exlude from statistics reads marked as duplicates\n");
        printf("    -f, --required-flag <int>           Required flag, 0 for unset [0]\n");
        printf("    -F, --filtering-flag <int>          Filtering flag, 0 for unset [0]\n");
        printf("        --GC-depth <float,float>        Bin size for GC-depth graph and the maximum reference length [2e4,4.2e9]\n");
        printf("    -h, --help                          This help message\n");
        printf("    -i, --insert-size <int>             Maximum insert size [8000]\n");
        printf("    -I, --id <string>                   Include only listed read group or sample name\n");
        printf("    -l, --read-length <int>             Include in the statistics only reads with the given read length []\n");
        printf("    -m, --most-inserts <float>          Report only the main part of inserts [0.99]\n");
        printf("    -q, --trim-quality <int>            The BWA trimming parameter [0]\n");
        printf("    -r, --ref-seq <file>                Reference sequence (required for GC-depth calculation).\n");
        printf("    -t, --target-regions <file>         Do stats in these regions only. Tab-delimited file chr,from,to, 1-based, inclusive.\n");
        printf("    -s, --sam                           Input is SAM\n");
        printf("\n");
    }
    else
    {
        va_list ap;
        va_start(ap, format);
        vfprintf(stderr, format, ap);
        va_end(ap);
    }
    exit(-1);
}

int main(int argc, char *argv[])
{
    char *targets = NULL;
    char *bam_fname = NULL;
    char *group_id = NULL;
    samfile_t *sam = NULL;
    char in_mode[5];

    stats_t *stats = calloc(1,sizeof(stats_t));
    stats->ngc    = 200;
    stats->nquals = 256;
    stats->nbases = 300;
    stats->nisize = 8000;
    stats->max_len   = 30;
    stats->max_qual  = 40;
    stats->isize_main_bulk = 0.99;   // There are always outliers at the far end
    stats->gcd_bin_size = 20e3;
    stats->gcd_ref_size = 4.2e9;
    stats->rseq_pos     = -1;
    stats->tid = stats->gcd_pos = -1;
    stats->igcd = 0;
    stats->is_sorted = 1;
    stats->cov_min  = 1;
    stats->cov_max  = 1000;
    stats->cov_step = 1;
    stats->argc = argc;
    stats->argv = argv;
    stats->filter_readlen = -1;
    stats->nindels = stats->nbases;

    strcpy(in_mode, "rb");

    static struct option loptions[] = 
    {
        {"help",0,0,'h'},
        {"remove-dups",0,0,'d'},
        {"sam",0,0,'s'},
        {"ref-seq",1,0,'r'},
        {"coverage",1,0,'c'},
        {"read-length",1,0,'l'},
        {"insert-size",1,0,'i'},
        {"most-inserts",1,0,'m'},
        {"trim-quality",1,0,'q'},
        {"target-regions",0,0,'t'},
        {"required-flag",1,0,'f'},
        {"filtering-flag",0,0,'F'},
        {"id",1,0,'I'},
        {"GC-depth",1,0,1},
        {0,0,0,0}
    };
    int opt;
    while ( (opt=getopt_long(argc,argv,"?hdsr:c:l:i:t:m:q:f:F:I:1:",loptions,NULL))>0 )
    {
        switch (opt)
        {
            case 'f': stats->flag_require=strtol(optarg,0,0); break;
            case 'F': stats->flag_filter=strtol(optarg,0,0); break;
            case 'd': stats->flag_filter|=BAM_FDUP; break;
            case 's': strcpy(in_mode, "r"); break;
            case 'r': stats->fai = fai_load(optarg); 
                      if (stats->fai==0) 
                          error("Could not load faidx: %s\n", optarg); 
                      break;
            case  1 : {
                        float flen,fbin;
                        if ( sscanf(optarg,"%f,%f",&fbin,&flen)!= 2 ) 
                            error("Unable to parse --GC-depth %s\n", optarg); 
                        stats->gcd_bin_size = fbin;
                        stats->gcd_ref_size = flen;
                      }
                      break;
            case 'c': if ( sscanf(optarg,"%d,%d,%d",&stats->cov_min,&stats->cov_max,&stats->cov_step)!= 3 ) 
                          error("Unable to parse -c %s\n", optarg); 
                      break;
            case 'l': stats->filter_readlen = atoi(optarg); break;
            case 'i': stats->nisize = atoi(optarg); break;
            case 'm': stats->isize_main_bulk = atof(optarg); break;
            case 'q': stats->trim_qual = atoi(optarg); break;
            case 't': targets = optarg; break;
            case 'I': group_id = optarg; break;
            case '?': 
            case 'h': error(NULL);
            default: error("Unknown argument: %s\n", optarg);
        }
    }
    if ( optind<argc )
        bam_fname = argv[optind++];

    if ( !bam_fname )
    {
        if ( isatty(fileno((FILE *)stdin)) )
            error(NULL);
        bam_fname = "-";
    }

    // Init structures
    //  .. coverage bins and round buffer
    if ( stats->cov_step > stats->cov_max - stats->cov_min + 1 )
    {
        stats->cov_step = stats->cov_max - stats->cov_min;
        if ( stats->cov_step <= 0 )
            stats->cov_step = 1;
    }
    stats->ncov = 3 + (stats->cov_max-stats->cov_min) / stats->cov_step;
    stats->cov_max = stats->cov_min + ((stats->cov_max-stats->cov_min)/stats->cov_step +1)*stats->cov_step - 1;
    stats->cov = calloc(sizeof(uint64_t),stats->ncov);
    stats->cov_rbuf.size = stats->nbases*5;
    stats->cov_rbuf.buffer = calloc(sizeof(int32_t),stats->cov_rbuf.size);
    // .. bam
    if ((sam = samopen(bam_fname, in_mode, NULL)) == 0) 
        error("Failed to open: %s\n", bam_fname);
    stats->sam = sam;
    if ( group_id ) init_group_id(stats, group_id);
    bam1_t *bam_line = bam_init1();
    // .. arrays
    stats->quals_1st      = calloc(stats->nquals*stats->nbases,sizeof(uint64_t));
    stats->quals_2nd      = calloc(stats->nquals*stats->nbases,sizeof(uint64_t));
    stats->gc_1st         = calloc(stats->ngc,sizeof(uint64_t));
    stats->gc_2nd         = calloc(stats->ngc,sizeof(uint64_t));
    stats->isize_inward   = calloc(stats->nisize,sizeof(uint64_t));
    stats->isize_outward  = calloc(stats->nisize,sizeof(uint64_t));
    stats->isize_other    = calloc(stats->nisize,sizeof(uint64_t));
    stats->gcd            = calloc(stats->ngcd,sizeof(gc_depth_t));
    stats->mpc_buf        = stats->fai ? calloc(stats->nquals*stats->nbases,sizeof(uint64_t)) : NULL;
    stats->acgt_cycles    = calloc(4*stats->nbases,sizeof(uint64_t));
    stats->read_lengths   = calloc(stats->nbases,sizeof(uint64_t));
    stats->insertions     = calloc(stats->nbases,sizeof(uint64_t));
    stats->deletions      = calloc(stats->nbases,sizeof(uint64_t));
    stats->ins_cycles_1st = calloc(stats->nbases+1,sizeof(uint64_t));
    stats->ins_cycles_2nd = calloc(stats->nbases+1,sizeof(uint64_t));
    stats->del_cycles_1st = calloc(stats->nbases+1,sizeof(uint64_t));
    stats->del_cycles_2nd = calloc(stats->nbases+1,sizeof(uint64_t));
    realloc_rseq_buffer(stats);
    if ( targets )
        init_regions(stats, targets);

    // Collect statistics
    if ( optind<argc )
    {
        // Collect stats in selected regions only
        bam_index_t *bam_idx = bam_index_load(bam_fname);
        if (bam_idx == 0)
            error("Random alignment retrieval only works for indexed BAM files.\n");

        int i;
        for (i=optind; i<argc; i++) 
        {
            int tid, beg, end;
            bam_parse_region(stats->sam->header, argv[i], &tid, &beg, &end);
            if ( tid < 0 ) continue;
            reset_regions(stats);
            bam_fetch(stats->sam->x.bam, bam_idx, tid, beg, end, stats, fetch_read);
        }
        bam_index_destroy(bam_idx);
    }
    else
    {
        // Stream through the entire BAM ignoring off-target regions if -t is given
        while (samread(sam,bam_line) >= 0) 
            collect_stats(bam_line,stats);
    }
    round_buffer_flush(stats,-1);

    output_stats(stats);

    bam_destroy1(bam_line);
    samclose(stats->sam);
    if (stats->fai) fai_destroy(stats->fai);
    free(stats->cov_rbuf.buffer); free(stats->cov);
    free(stats->quals_1st); free(stats->quals_2nd); 
    free(stats->gc_1st); free(stats->gc_2nd);
    free(stats->isize_inward); free(stats->isize_outward); free(stats->isize_other);
    free(stats->gcd);
    free(stats->rseq_buf);
    free(stats->mpc_buf);
    free(stats->acgt_cycles);
    free(stats->read_lengths);
    free(stats->insertions);
    free(stats->deletions);
    free(stats->ins_cycles_1st);
    free(stats->ins_cycles_2nd);
    free(stats->del_cycles_1st);
    free(stats->del_cycles_2nd);
    destroy_regions(stats);
    free(stats);
    if ( stats->rg_hash ) kh_destroy(kh_rg, stats->rg_hash);

    return 0;
}



