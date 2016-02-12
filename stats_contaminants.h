#include <htslib/khash.h>
#include <htslib/kstring.h>
#include <htslib/sam.h>

#include <stdint.h>

/* Fixed-sized integer array that is used in the contaminants hash map
 * for mapping k-mers to sequences.
 */
typedef struct
{
    /* The number of entries */
    uint32_t len;
    /* The list of sequence ids */
    uint32_t * vals;
} contaminant_map_entry_t;

KHASH_MAP_INIT_STR(contaminants, contaminant_map_entry_t);

/* Data structure for storing k-mer information about contaminants
 */
typedef struct
{
    /* The k-mer length */
    uint32_t k;
    /* Hash map from k-mer to contaminant sequence id */
    khash_t(contaminants) * map;

    /* Sequence information from FASTA file */
    kstring_t * names;
    kstring_t * seqs;
    /* Number of shared kmers for each seq */
    uint64_t * shared_kmers;
    /* Number of sequences */
    uint32_t n_seqs;
    /* Number of reads counted in total */
    uint64_t total_reads;

    /* Whether or not contaminant was hit for current read. */
    uint32_t * _was_hit;
    /* Helper buffer */
    char * _buf, * _kmer;
}
contaminant_map_t;

/* Load contaminant_map_t from FASTA file and return pointer to it.
 *
 * Returns NULL on problems with reading the FASTA file.
 *
 * Use free_contaminant_map for freeing up memory again. */
contaminant_map_t * load_contaminant_fasta(const char * path, uint32_t k);

/* Free contaminant_map_t again */
void free_contaminant_map(contaminant_map_t * map);

/* Collect information about contaminants */
void collect_contaminant_info(contaminant_map_t * map, const bam1_t * bam_line);
