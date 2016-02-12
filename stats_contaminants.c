#include "stats_contaminants.h"

#include <zlib.h>
#include <stdio.h>

#include <config.h>
#include <htslib/kseq.h>

KSEQ_INIT(gzFile, gzread)

/* Return contaminant map entry for the given k-mer.
 *
 * Allocate new entry if necessary.
 */
khiter_t contaminant_map_get_bucket(
        contaminant_map_t * map, char * kmer)
{
    khiter_t k = kh_get(contaminants, map->map, kmer);
    if (k == kh_end(map->map)) {
        // Register new occurences entry with hash table
        int ret = 0;
        k = kh_put(contaminants, map->map, strdup(kmer), &ret);
        if (ret == -1)
            exit(1);
        contaminant_map_entry_t * ptr = &kh_value(map->map, k);
        ptr->len = 0;
        ptr->vals = NULL;
    }
    return k;
}

/* Register sequence in bucket */
void contaminant_map_entry_register_in_bucket(
        contaminant_map_entry_t * entry, uint32_t seq_id)
{
    // printf("registering entry %p for seq %u\n", entry, seq_id);
    uint32_t i;
    for (i = 0; i < entry->len; ++i)
        if (entry->vals[i] == seq_id)
            return;  // already there

    // Not in bucket yet, register
    // printf("vals size is %d\n", entry->len);
    entry->vals = realloc(entry->vals, (entry->len + 1) * sizeof(uint32_t));
    entry->vals[entry->len] = seq_id;
    entry->len += 1;
}

/** Reverse-complement str, in-place */
void revcomp(char * str, uint32_t k)
{
    uint32_t i;
    char * ptr;

    // complement first
    for (i = 0, ptr = str; i < k; ++i) {
        switch (*ptr) {
            case 'A': *ptr = 'T'; break;
            case 'C': *ptr = 'G'; break;
            case 'G': *ptr = 'C'; break;
            case 'T': *ptr = 'A'; break;
        }
    }
    // then reverse
    char tmp;
    for (i = 0; i < k / 2; ++i) {
        tmp = str[i];
        str[i] = str[k - i - 1];
        str[k - i - 1] = tmp;
    }
}

contaminant_map_t * load_contaminant_fasta(const char * path, uint32_t k)
{
    // fprintf(stderr, "Loading contaminant FASTA...\n");
    // Initial memory allocation
    contaminant_map_t * result = calloc(1, sizeof(contaminant_map_t));
    result->map = kh_init(contaminants);
    // Initialize other values
    result->k = k;
    result->names = NULL;
    result->seqs = NULL;
    result->shared_kmers = NULL;
    result->total_reads = 0;
    result->n_seqs = 0;
    result->_was_hit = NULL;
    result->_buf = NULL;
    result->_kmer = NULL;

    // Open file
    gzFile fp;
    kseq_t *seq;
    fp = strcmp(path, "-") ? gzopen(path, "r") : gzdopen(fileno(stdin), "r");
    if (fp == 0) {
        kh_destroy(contaminants, result->map);
        free(result);
        fprintf(stderr, "Could not open FASTA file %s\n", path);
        return NULL;
    }

    seq = kseq_init(fp);

    int contaminant_id = 0;
    int l;
    while ((l = kseq_read(seq)) >= 0) {
        // Allocate space for new name/sequence pair
        result->n_seqs += 1;
        result->names = realloc(result->names, result->n_seqs * sizeof(kstring_t));
        result->seqs = realloc(result->seqs, result->n_seqs * sizeof(kstring_t));
        result->shared_kmers = realloc(result->shared_kmers, result->n_seqs * sizeof(uint64_t));
        result->shared_kmers[result->n_seqs - 1] = 0;
        result->_was_hit = realloc(result->_was_hit, result->n_seqs* sizeof(uint32_t));
        result->_was_hit[result->n_seqs - 1] = 0;
        kstring_t * curr_name = &(result->names[result->n_seqs - 1]);
        kstring_t * curr_seq = &(result->seqs[result->n_seqs - 1]);
        curr_name->l = 0;
        curr_name->m = 0;
        curr_name->s = NULL;
        curr_seq->l = 0;
        curr_seq->m = 0;
        curr_seq->s = NULL;

        // Copy over name
        kputs(seq->name.s, curr_name);
        kputs(" ", curr_name);
        kputs(seq->comment.s, curr_name);
        // Copy over sequence and convert to upper case
        kputs(seq->seq.s, curr_seq);
        char * ptr;
        for (ptr = curr_seq->s; ptr != curr_seq->s + curr_seq->l; ++ptr)
            *ptr = toupper(*ptr);

        // printf("NAME: %s\tSEQ: %s\n", curr_name->s, curr_seq->s);

        // Register contaminant for the shared k-mers
        uint32_t i, j;
        char * kmer = calloc(k + 1, sizeof(char));
        kmer[k] = '\0';
        for (i = 0; i < curr_seq->l - k + 1; ++i) {
            int has_n = 0;
            for (j = 0; j < k; ++j) {
                kmer[j] = curr_seq->s[i + j];
                if (kmer[j] == 'N')
                    has_n = 1;
            }
            if (has_n)
                continue;  // Skip k-mers containing N

            // printf("\tregistering kmer %s\n", kmer);
            khiter_t iter = contaminant_map_get_bucket(result, kmer);
            contaminant_map_entry_register_in_bucket(
                    &kh_value(result->map, iter), contaminant_id);
            revcomp(kmer, k);
            // printf("\tregistering kmer (revcomp) %s\n", kmer);
            iter = contaminant_map_get_bucket(result, kmer);
            contaminant_map_entry_register_in_bucket(
                    &kh_value(result->map, iter), contaminant_id);
        }
        free(kmer);

        contaminant_id += 1;
    }
    kseq_destroy(seq);
    // fprintf(stderr, "Done loading contaminant FASTA\n");

    return result;
}

void free_contaminant_map(contaminant_map_t * map)
{
    contaminant_map_entry_t entry;
    kh_foreach_value(map->map, entry, free(entry.vals));
    kh_destroy(contaminants, map->map);

    uint32_t i;
    for (i = 0; i < map->n_seqs; ++i) {
        ks_release(&(map->names[i]));
        ks_release(&(map->seqs[i]));
    }
    free(map->names);
    free(map->seqs);
    free(map->_was_hit);
    free(map->_buf);
    free(map->_kmer);

    free(map);
}

/* Return C char for BAM char (N if not ACGT) */
char bam_char_to_c_char(int8_t c)
{
    switch (c)
    {
        case 1: return 'A';
        case 2: return 'C';
        case 4: return 'G';
        case 8: return 'T';
        default: return 'N';
    }
}

void collect_contaminant_info(contaminant_map_t * map, const bam1_t * bam_line)
{
    // Get char sequence from BAM line
    uint8_t * seq = bam_get_seq(bam_line);
    uint32_t len = bam_line->core.l_qseq;
    uint32_t i, j;
    map->_buf = realloc(map->_buf, (len + 1) * sizeof(char));
    map->_buf[len] = '\0';
    for (i = 0; i < len; ++i)
        map->_buf[i] = bam_char_to_c_char(bam_seqi(seq, i));

    map->total_reads += 1;

    // Reset _was_hit map
    for (i = 0; i < map->n_seqs; ++i)
        map->_was_hit[i] = 0;

    // Enumerate all k-mers and check for sharing k-mer with contaminants.
    map->_kmer = realloc(map->_kmer, (map->k + 1) * sizeof(char));
    for (i = 0; i < len - map->k + 1; ++i) {
        for (j = 0; j < map->k; ++j)
            map->_kmer[j] = map->_buf[i + j];
        map->_kmer[map->k] = '\0';

        khiter_t k = kh_get(contaminants, map->map, map->_kmer);
        if (k != kh_end(map->map)) {
            contaminant_map_entry_t * bucket = &kh_value(map->map, k);
            int k;
            for (k = 0; k < bucket->len; ++k) {
                map->_was_hit[bucket->vals[k]] = 1;
                // printf("Read %s hit adapter %s with kmer %s\n",
                //        map->_buf, map->names[bucket->vals[k]].s, map->_kmer);
            }
        }
    }

    // Increment counts for hit contaminants
    for (i = 0; i < map->n_seqs; ++i)
        map->shared_kmers[i] += map->_was_hit[i];
}
