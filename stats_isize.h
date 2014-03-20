#include <htslib/khash.h>
#include <stdint.h>

typedef struct
{
    int total;
    uint64_t *isize_inward, *isize_outward, *isize_other;
}
isize_dense_data_t;

typedef struct 
{ 
    uint64_t isize_inward, isize_outward, isize_other; 
}
isize_sparse_record_t;

KHASH_MAP_INIT_INT(m32, isize_sparse_record_t *)

typedef struct
{
    int max;
    khash_t(m32) *array;
}
isize_sparse_data_t;

typedef union {
    isize_sparse_data_t *sparse;
    isize_dense_data_t *dense;
} isize_data_t;

// Insert size structure
typedef struct
{
    isize_data_t data;
    
    // Maximum
    int (*nitems)(isize_data_t);

    // Fetch the number of inserts of a given size
    uint64_t (*inward)(isize_data_t, int);
    uint64_t (*outward)(isize_data_t, int);
    uint64_t (*other)(isize_data_t, int);
    uint64_t (*total)(isize_data_t, int);

    // Set the number of inserts of a given size
    void (*set_inward)(isize_data_t, int, uint64_t);
    void (*set_outward)(isize_data_t, int, uint64_t);
    void (*set_other)(isize_data_t, int, uint64_t);

    // Increment the number of inserts of a given size
    void (*inc_inward)(isize_data_t, int);
    void (*inc_outward)(isize_data_t, int);
    void (*inc_other)(isize_data_t, int);

    // Free this structure
    void (*isize_free)(isize_data_t);
}
isize_t;

isize_t *init_isize_t(int bound);
