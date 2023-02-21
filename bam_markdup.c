/*  bam_markdup.c -- Mark duplicates from a coord sorted file that has gone
                     through fixmates with the mate scoring option on.

    Copyright (C) 2017-2023 Genome Research Ltd.

    Author: Andrew Whitwham <aw7@sanger.ac.uk>

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
DEALINGS IN THE SOFTWARE

Estimate library size derived from Picard DuplicationMetrics.java
Copyright (c) 2009,2018 The Broad Institute.  MIT license.
*/

#include <config.h>

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <ctype.h>
#include <time.h>
#include <sys/stat.h>
#include <math.h>
#include <regex.h>
#include "htslib/thread_pool.h"
#include "htslib/sam.h"
#include "sam_opts.h"
#include "samtools.h"
#include "htslib/khash.h"
#include "htslib/klist.h"
#include "htslib/kstring.h"
#include "tmp_file.h"
#include "bam.h"


typedef struct {
    samFile *in;
    samFile *out;
    char *prefix;
    int remove_dups;
    int32_t max_length;
    int do_stats;
    int supp;
    int tag;
    int opt_dist;
    int no_pg;
    int clear;
    int mode;
    int write_index;
    int include_fails;
    int check_chain;
    char *stats_file;
    char *arg_list;
    char *out_fn;
    regex_t *rgx;
    int rgx_x;
    int rgx_y;
    int rgx_t;
    char *barcode;
    regex_t *bc_rgx;
    int read_groups;
    int json;
} md_param_t;

typedef struct {
    hts_pos_t this_coord;
    hts_pos_t other_coord;
    int32_t this_ref;
    int32_t other_ref;
    int32_t barcode;
    int32_t read_group;
    int8_t single;
    int8_t leftmost;
    int8_t orientation;
} key_data_t;

typedef struct read_queue_s {
    key_data_t pair_key;
    key_data_t single_key;
    bam1_t *b;
    struct read_queue_s *duplicate;
    struct read_queue_s *original;
    hts_pos_t pos;
    int dup_checked;
    int read_group;
} read_queue_t;

typedef struct {
    read_queue_t *p;
} in_hash_t;

typedef struct {
    char *name;
    char type;
    int read_group;
} dup_map_t;

typedef struct {
    bam1_t *b;
    int64_t score;
    int64_t mate_score;
    long x;
    long y;
    int opt;
    int beg;
    int end;
} check_t;

typedef struct {
    check_t *c;
    size_t size;
    size_t length;
} check_list_t;

typedef struct {
    long reading;
    long writing;
    long excluded;
    long duplicate;
    long single;
    long pair;
    long single_dup;
    long examined;
    long optical;
    long single_optical;
    long np_duplicate;
    long np_opt_duplicate;
} stats_block_t;

static khint32_t do_hash(unsigned char *key, khint32_t len);

static khint_t hash_key(key_data_t key) {
    int i = 0;
    khint_t hash;

    if (key.single) {
        unsigned char sig[21];

        memcpy(sig + i, &key.this_ref, 4);      i += 4;
        memcpy(sig + i, &key.this_coord, 8);    i += 8;
        memcpy(sig + i, &key.orientation, 1);   i += 1;
        memcpy(sig + i, &key.barcode, 4);       i += 4;
        memcpy(sig + i, &key.read_group, 4);    i += 4;

        hash = do_hash(sig, i);
    } else {
        unsigned char sig[34];

        memcpy(sig + i, &key.this_ref, 4);      i += 4;
        memcpy(sig + i, &key.this_coord, 8);    i += 8;
        memcpy(sig + i, &key.other_ref, 4);     i += 4;
        memcpy(sig + i, &key.other_coord, 8);   i += 8;
        memcpy(sig + i, &key.leftmost, 1);      i += 1;
        memcpy(sig + i, &key.orientation, 1);   i += 1;
        memcpy(sig + i, &key.barcode, 4);       i += 4;
        memcpy(sig + i, &key.read_group, 4);    i += 4;

        hash = do_hash(sig, i);
    }

    return hash;
}


static int key_equal(key_data_t a, key_data_t b) {
    int match = 1;

    if (a.this_coord != b.this_coord)
        match = 0;
    else if (a.orientation != b.orientation)
        match = 0;
    else if (a.this_ref != b.this_ref)
        match = 0;
    else if (a.single != b.single)
        match = 0;
    else if (a.barcode != b.barcode)
        match = 0;
    else if (a.read_group != b.read_group)
        match = 0;

    if (!a.single) {
        if (a.other_coord != b.other_coord)
            match = 0;
        else if (a.leftmost != b.leftmost)
            match = 0;
        else if (a.other_ref != b.other_ref)
            match = 0;
    }

    return match;
}


#define __free_queue_element(p)

// Orientations (prime numbers to feed to hashing algorithm)
#define O_FF 2
#define O_RR 3
#define O_FR 5
#define O_RF 7

// Left or rightmost
#define R_LE 11
#define R_RI 13

#define BMD_WARNING_MAX 10

#define MD_MIN_QUALITY 15

// Duplicate finding mode
#define MD_MODE_TEMPLATE 0
#define MD_MODE_SEQUENCE 1

KHASH_INIT(reads, key_data_t, in_hash_t, 1, hash_key, key_equal) // read map hash
KLIST_INIT(read_queue, read_queue_t, __free_queue_element) // the reads buffer
KHASH_MAP_INIT_STR(duplicates, dup_map_t) // map of duplicates for supplementary dup id
KHASH_MAP_INIT_STR(read_groups, int) // read group lookup

/* The Bob Jenkins one_at_a_time hash to reduce the key to a 32 bit value. */

static khint32_t do_hash(unsigned char *key, khint32_t len) {
    khint32_t   hash, i;

    for (hash = 0, i = 0; i < len; ++i) {
        hash += key[i];
        hash += (hash << 10);
        hash ^= (hash >> 6);
    }

    hash += (hash << 3);
    hash ^= (hash >> 11);
    hash += (hash << 15);

    return hash;
}


/* Get mate score from tag. */

static int64_t get_mate_score(bam1_t *b) {
    uint8_t *data;
    int64_t score;

    if ((data = bam_aux_get(b, "ms"))) {
        score = bam_aux2i(data);
    } else {
        print_error("markdup", "error, no ms score tag. Please run samtools fixmate on file first.\n");
        return -1;
    }

    return score;
}


/* Calc current score from quality. */

static int64_t calc_score(bam1_t *b)
{
    int64_t score = 0;
    uint8_t  *qual = bam_get_qual(b);
    int i;

    for (i = 0; i < b->core.l_qseq; i++) {
        if (qual[i] >= MD_MIN_QUALITY) score += qual[i];
    }

    return score;
}


/* Create a signature hash of the current read and its pair.
   Uses the unclipped start (or end depending on orientation),
   the reference id, orientation and whether the current
   read is leftmost of the pair. */


static int make_pair_key(md_param_t *param, key_data_t *key, bam1_t *bam, int rg_num, long *warnings) {
     hts_pos_t this_coord, this_end, other_coord, other_end, leftmost;
    int32_t this_ref, other_ref, barcode = 0;
    int8_t orientation, left_read;
    uint8_t *data;
    char *cig, *bar;
    long incoming_warnings = *warnings;

    this_ref    = bam->core.tid + 1; // avoid a 0 being put into the hash
    other_ref   = bam->core.mtid + 1;

    this_coord = unclipped_start(bam);
    this_end   = unclipped_end(bam);

    if ((data = bam_aux_get(bam, "MC"))) {
        if (!(cig = bam_aux2Z(data))) {
            print_error("markdup", "error, MC tag wrong type. Please use the MC tag provided by samtools fixmate.\n");
            return 1;
        }

        other_end   = unclipped_other_end(bam->core.mpos, cig);
        other_coord = unclipped_other_start(bam->core.mpos, cig);
    } else {
        print_error("markdup", "error, no MC tag. Please run samtools fixmate on file first.\n");
        return 1;
    }

    // work out orientations
    if (param->mode == MD_MODE_TEMPLATE) {

        if (this_ref != other_ref) {
            leftmost = this_ref < other_ref;
        } else {
            if (bam_is_rev(bam) == bam_is_mrev(bam)) {
                if (!bam_is_rev(bam)) {
                    leftmost = this_coord <= other_coord;
                } else {
                    leftmost = this_end <= other_end;
                }
            } else {
                if (bam_is_rev(bam)) {
                    leftmost = this_end <= other_coord;
                } else {
                    leftmost = this_coord <= other_end;
                }
            }
        }

        // pair orientation
        if (leftmost) {
            if (bam_is_rev(bam) == bam_is_mrev(bam)) {
                other_coord = other_end;

                if (!bam_is_rev(bam)) {
                    if (bam->core.flag & BAM_FREAD1) {
                        orientation = O_FF;
                    } else {
                        orientation = O_RR;
                    }
                } else {
                    if (bam->core.flag & BAM_FREAD1) {
                        orientation = O_RR;
                    } else {
                        orientation = O_FF;
                    }
                }
            } else {
                if (!bam_is_rev(bam)) {
                    orientation = O_FR;
                    other_coord = other_end;
                } else {
                    orientation = O_RF;
                    this_coord = this_end;
                }
            }
        } else {
            if (bam_is_rev(bam) == bam_is_mrev(bam)) {
                this_coord = this_end;

                if (!bam_is_rev(bam)) {
                    if (bam->core.flag & BAM_FREAD1) {
                        orientation = O_RR;
                    } else {
                        orientation = O_FF;
                    }
                } else {
                    if (bam->core.flag & BAM_FREAD1) {
                        orientation = O_FF;
                    } else {
                        orientation = O_RR;
                    }
                }
            } else {
                if (!bam_is_rev(bam)) {
                    orientation = O_RF;
                    other_coord = other_end;
                } else {
                    orientation = O_FR;
                    this_coord = this_end;
                }
            }
        }
    } else { // MD_MODE_SEQUENCE

        if (this_ref != other_ref) {
            leftmost = this_ref - other_ref;
        } else {
            if (bam_is_rev(bam) == bam_is_mrev(bam)) {
                if (!bam_is_rev(bam)) {
                    leftmost = this_coord - other_coord;
                } else {
                    leftmost = this_end - other_end;
                }
            } else {
                if (bam_is_rev(bam)) {
                    leftmost = this_end - other_coord;
                } else {
                    leftmost = this_coord - other_end;
                }
            }
        }

        if (leftmost < 0) {
            leftmost = 1;
        } else if (leftmost > 0) {
            leftmost = 0;
        } else {
            // tie breaks

            if (bam->core.pos == bam->core.mpos) {
                if (bam->core.flag & BAM_FREAD1) {
                    leftmost = 1;
                } else {
                    leftmost = 0;
                }
            } else if (bam->core.pos < bam->core.mpos) {
                leftmost = 1;
            } else {
                leftmost = 0;
            }
        }

        // pair orientation
        if (leftmost) {
            if (bam_is_rev(bam) == bam_is_mrev(bam)) {

                if (!bam_is_rev(bam)) {
                        orientation = O_FF;
                } else {
                        orientation = O_RR;
                }
            } else {
                if (!bam_is_rev(bam)) {
                    orientation = O_FR;
                } else {
                    orientation = O_RF;
                }
            }
        } else {
            if (bam_is_rev(bam) == bam_is_mrev(bam)) {

                if (!bam_is_rev(bam)) {
                        orientation = O_RR;
                } else {
                        orientation = O_FF;
                }
            } else {
                if (!bam_is_rev(bam)) {
                    orientation = O_RF;
                } else {
                    orientation = O_FR;
                }
            }
        }

        if (!bam_is_rev(bam)) {
            this_coord = unclipped_start(bam);
        } else {
            this_coord = unclipped_end(bam);
        }

        if (!bam_is_mrev(bam)) {
            other_coord = unclipped_other_start(bam->core.mpos, cig);
        } else {
            other_coord = unclipped_other_end(bam->core.mpos, cig);
        }
    }

    if (!leftmost)
        left_read = R_RI;
    else
        left_read = R_LE;

    if (param->barcode) {
        if ((data = bam_aux_get(bam, param->barcode))) {
            if (!(bar = bam_aux2Z(data))) {
                (*warnings)++;

                if (*warnings <= BMD_WARNING_MAX) {
                    print_error("markdup", "warning, %s tag wrong type. Aux tag needs to be a string type.\n", param->barcode);
                }
            } else {
                barcode = do_hash((unsigned char *)bar, strlen(bar));
            }
        }
    } else if (param->bc_rgx) {
        int result;
        regmatch_t matches[3];
        size_t max_matches = 2;
        char *qname = bam_get_qname(bam);

        if ((result = regexec(param->bc_rgx, qname, max_matches, matches, 0)) == 0) {
            int bc_start, bc_end;

            bc_start = matches[1].rm_so;
            bc_end   = matches[1].rm_eo;

            if (bc_start != -1) {
                barcode = do_hash((unsigned char *)qname + bc_start, bc_end - bc_start);
            } else {
                (*warnings)++;

                if (*warnings <= BMD_WARNING_MAX) {
                    print_error("markdup", "warning, barcode regex unable to match substring on %s.\n", qname);
                }
            }
        } else {
            (*warnings)++;

            if (*warnings <= BMD_WARNING_MAX) {
                char warn_msg[256];

                regerror(result, param->bc_rgx, warn_msg, 256);
                print_error("markdup", "warning, barcode regex match error \"%s\" on %s.\n", warn_msg, qname);
            }
        }
    }

    if ((*warnings == BMD_WARNING_MAX) && (incoming_warnings != *warnings)) {
        print_error("markdup", "warning, %ld barcode read warnings.  New warnings will not be reported.\n",
                        *warnings);
    }

    key->single        = 0;
    key->this_ref      = this_ref;
    key->this_coord    = this_coord;
    key->other_ref     = other_ref;
    key->other_coord   = other_coord;
    key->leftmost      = left_read;
    key->orientation   = orientation;
    key->barcode       = barcode;
    key->read_group    = rg_num;

    return 0;
}


/* Create a signature hash of single read (or read with an unmatched pair).
   Uses unclipped start (or end depending on orientation), reference id,
   and orientation. */

static void make_single_key(md_param_t *param, key_data_t *key, bam1_t *bam, int rg_num, long *warnings) {
    hts_pos_t this_coord;
    int32_t this_ref, barcode = 0;
    int8_t orientation;
    uint8_t *data;
    char *bar;
    long incoming_warnings = *warnings;

    this_ref = bam->core.tid + 1; // avoid a 0 being put into the hash

    if (bam_is_rev(bam)) {
        this_coord = unclipped_end(bam);
        orientation = O_RR;
    } else {
        this_coord = unclipped_start(bam);
        orientation = O_FF;
    }

    if (param->barcode) {
        if ((data = bam_aux_get(bam, param->barcode))) {
            if (!(bar = bam_aux2Z(data))) {
                (*warnings)++;

                if (*warnings <= BMD_WARNING_MAX) {
                    print_error("markdup", "warning, %s tag wrong type. Aux tag needs to be a string type.\n", param->barcode);
                }
            } else {
                barcode = do_hash((unsigned char *)bar, strlen(bar));
            }
        }
    } else if (param->bc_rgx) {
        int result;
        regmatch_t matches[3];
        size_t max_matches = 2;
        char *qname = bam_get_qname(bam);

        if ((result = regexec(param->bc_rgx, qname, max_matches, matches, 0)) == 0) {
            int bc_start, bc_end;

            bc_start = matches[1].rm_so;
            bc_end   = matches[1].rm_eo;

            if (bc_start != -1) {
                barcode = do_hash((unsigned char *)qname + bc_start, bc_end - bc_start);
            } else {
                (*warnings)++;

                if (*warnings <= BMD_WARNING_MAX) {
                    print_error("markdup", "warning, barcode regex unable to match substring on %s.\n", qname);
                }
            }
        } else {
            (*warnings)++;

            if (*warnings <= BMD_WARNING_MAX) {
                char warn_msg[256];

                regerror(result, param->bc_rgx, warn_msg, 256);
                print_error("markdup", "warning, barcode regex match error \"%s\" on %s.\n", warn_msg, qname);
            }
        }
    }

    if ((*warnings == BMD_WARNING_MAX) && (incoming_warnings != *warnings)) {
        print_error("markdup", "warning, %ld barcode read warnings.  New warnings will not be reported.\n",
                        *warnings);
    }


    key->single        = 1;
    key->this_ref      = this_ref;
    key->this_coord    = this_coord;
    key->orientation   = orientation;
    key->barcode       = barcode;
    key->read_group    = rg_num;
}


/* Add the duplicate name to a hash if it does not exist. */

static int add_duplicate(khash_t(duplicates) *d_hash, bam1_t *dupe, char *orig_name, char type, int group) {
    khiter_t d;
    int ret;

    d = kh_get(duplicates, d_hash, bam_get_qname(dupe));

    if (d == kh_end(d_hash)) {
        char *name = strdup(bam_get_qname(dupe));
        if (name) {
            d = kh_put(duplicates, d_hash, name, &ret);
        } else {
            ret = -1;
        }

        if (ret >= 0) {
            if (orig_name) {
                if (ret == 0) {
                    // replace old name
                    free(kh_value(d_hash, d).name);
                    free(name);
                }

                kh_value(d_hash, d).name = strdup(orig_name);

                if (kh_value(d_hash, d).name == NULL) {
                    print_error("markdup", "error, unable to allocate memory for duplicate original name.\n");
                    return 1;
                }
            } else {
                kh_value(d_hash, d).name = NULL;
            }

            kh_value(d_hash, d).type = type;
            kh_value(d_hash, d).read_group = group;
        } else {
            print_error("markdup", "error, unable to store supplementary duplicates.\n");
            free(name);
            return 1;
        }
    }

    return 0;
}


/* Get coordinates from the standard Illumina style read names.
   Returned values are of the x and y coordinates and a section of
   the read name to test (t) for string equality e.g. lane and tile part. */

static int get_coordinates_colons(md_param_t *param, const char *qname, int *t_beg, int *t_end, long *x_coord, long *y_coord, long *warnings) {
    int sep = 0;
    int pos = 0;
    int xpos = 0, ypos = 0;
    char *end;

    while (qname[pos]) {
        if (qname[pos] == ':') {
            sep++;

            if (sep == 2) {
                xpos = pos + 1;
            } else if (sep == 3) {
                ypos = pos + 1;
            } else if (sep == 4) { // HiSeq style names
                xpos = ypos;
                ypos = pos + 1;
            } else if (sep == 5) { // Newer Illumina format
                xpos = pos + 1;
            } else if (sep == 6) {
                ypos = pos + 1;
            }
        }

        pos++;
    }

    /* The most current Illumina read format at time of writing is:
       @machine:run:flowcell:lane:tile:x:y:UMI or
       @machine:run:flowcell:lane:tile:x:y

       Counting the separating colons gives us a quick format check.
       Older name formats have fewer elements.
    */

    if (!(sep == 3 || sep == 4 || sep == 6 || sep == 7)) {
        (*warnings)++;

        if (*warnings <= BMD_WARNING_MAX) {
            print_error("markdup", "warning, cannot decipher read name %s for optical duplicate marking.\n", qname);
        }

        return 1;
    } else {
        *x_coord = strtol(qname + xpos, &end, 10);

        if ((qname + xpos) == end) {
            (*warnings)++;

            if (*warnings <= BMD_WARNING_MAX) {
                print_error("markdup", "warning, cannot decipher x coordinate in %s .\n", qname);
            }

            return 1;
        }

        *y_coord = strtol(qname + ypos, &end, 10);

        if ((qname + ypos) == end) {
            (*warnings)++;

            if (*warnings <= BMD_WARNING_MAX) {
                print_error("markdup", "warning, cannot decipher y coordinate in %s .\n", qname);
            }

            return 1;
        }

        *t_beg = 0;
        *t_end = xpos;
    }

    return 0;
}

/* Get the coordinates from the read name.
   Returned values are of the x and y coordinates and an optional section of
   the read name to test (t) for string equality e.g. lane and tile part. */

static inline int get_coordinates_regex(md_param_t *param, const char *qname, int *t_beg, int *t_end, long *x_coord, long *y_coord, long *warnings) {
    regmatch_t matches[5];
    size_t max_matches = 5;
    int xpos, ypos, xend, yend, xlen, ylen;
    char coord[255];
    char *end;

    if (!param->rgx_t)
        max_matches = 4;

    if (regexec(param->rgx, qname, max_matches, matches, 0))
        return -1;

    xpos = matches[param->rgx_x].rm_so;
    ypos = matches[param->rgx_y].rm_so;

    if (param->rgx_t) {
        *t_beg = matches[param->rgx_t].rm_so;
        *t_end = matches[param->rgx_t].rm_eo;
    } else {
        *t_beg = *t_end = 0;
    }

    if (xpos == -1 || ypos == -1 || *t_beg == -1)
        return -1;

    xend = matches[param->rgx_x].rm_eo;
    yend = matches[param->rgx_y].rm_eo;

    if ((xlen = xend - xpos) > 254) {
        (*warnings)++;

        if (*warnings <= BMD_WARNING_MAX) {
            print_error("markdup", "warning, x coordinate string longer than allowed qname length in %s (%d long).\n", qname, xlen);
        }

        return 1;
    }

    strncpy(coord, qname + xpos, xlen);
    coord[xlen] = '\0';
    *x_coord = strtol(coord, &end, 10);

    if (coord == end) {
        (*warnings)++;

        if (*warnings <= BMD_WARNING_MAX) {
            print_error("markdup", "warning, cannot decipher x coordinate in %s (%s).\n", qname, coord);
        }

        return 1;
    }

    if ((ylen = yend - ypos) > 254) {
        (*warnings)++;

        if (*warnings <= BMD_WARNING_MAX) {
            print_error("markdup", "warning, y coordinate string longer than allowed qname length in %s (%d long).\n", qname, ylen);
        }

        return 1;
    }

    strncpy(coord, qname + ypos, ylen);
    coord[ylen] = '\0';
    *y_coord = strtol(coord, &end, 10);

    if (coord == end) {
        (*warnings)++;

        if (*warnings <= BMD_WARNING_MAX) {
            print_error("markdup", "warning, cannot decipher y coordinate in %s (%s).\n", qname, coord);
        }

        return 1;
    }

    return 0;
}


static int get_coordinates(md_param_t *param, const char *name, int *t_beg, int *t_end, long *x_coord, long *y_coord, long *warnings) {
    int ret = 1;

    if (param->rgx == NULL) {
        ret = get_coordinates_colons(param, name, t_beg, t_end, x_coord, y_coord, warnings);
    } else {
        ret = get_coordinates_regex(param, name, t_beg, t_end, x_coord, y_coord, warnings);
    }

    return ret;
}


/* Using the coordinates from the read name, see whether the duplicated read is
   close enough (set by max_dist) to the original to be counted as optical.*/

static int is_optical_duplicate(md_param_t *param, bam1_t *ori, bam1_t *dup, long max_dist, long *warnings) {
    int ret = 0;
    char *original, *duplicate;
    long ox, oy, dx, dy;
    int o_beg = 0, o_end = 0, d_beg = 0, d_end = 0;

    original  = bam_get_qname(ori);
    duplicate = bam_get_qname(dup);

    if (get_coordinates(param, original, &o_beg, &o_end, &ox, &oy, warnings)) {
        return ret;
    }

    if (get_coordinates(param, duplicate, &d_beg, &d_end, &dx, &dy, warnings)) {
        return ret;
    }

    if (strncmp(original + o_beg, duplicate + d_beg, o_end - o_beg) == 0) {
        long xdiff, ydiff;

        if (ox > dx) {
            xdiff = ox - dx;
        } else {
            xdiff = dx - ox;
        }

        if (xdiff <= max_dist) {
            // still might be optical

            if (oy > dy) {
                ydiff = oy - dy;
            } else {
                ydiff = dy - oy;
            }

            if (ydiff <= max_dist) ret = 1;
        }
    }

    return ret;
}


/* Using the coordinates from the Illumina read name, see whether the duplicated read is
   close enough (set by max_dist) to the original to be counted as optical.

   This function needs the values from the first read to be already calculated. */

static int optical_duplicate_partial(md_param_t *param, const char *name, const int o_beg, const int o_end, const long ox, const long oy, bam1_t *dup, check_t *c, long max_dist, long *warnings) {
    int ret = 0;
    char *duplicate;
    int d_beg = 0, d_end = 0;
    long dx, dy;

    duplicate = bam_get_qname(dup);

    if (get_coordinates(param, duplicate, &d_beg, &d_end, &dx, &dy, warnings)) {
        return ret;
    }

    if (strncmp(name + o_beg, duplicate + d_beg, o_end - o_beg) == 0) {
        // the initial parts match, look at the numbers
        long xdiff, ydiff;

        if (ox > dx) {
            xdiff = ox - dx;
        } else {
            xdiff = dx - ox;
        }

        if (xdiff <= max_dist) {
            // still might be optical

            if (oy > dy) {
                ydiff = oy - dy;
            } else {
                ydiff = dy - oy;
            }

            if (ydiff <= max_dist) ret = 1;
        }
    }

    c->x = dx;
    c->y = dy;
    c->beg = d_beg;
    c->end = d_end;

    return ret;
}


/* Mark the read as a duplicate and update the duplicate hash (if needed) */
static int mark_duplicates(md_param_t *param, khash_t(duplicates) *dup_hash, bam1_t *ori, bam1_t *dup,
                           int read_group, long *optical, long *warn) {
    char dup_type = 0;
    long incoming_warnings = *warn;

    dup->core.flag |= BAM_FDUP;

    if (param->tag) {
        if (bam_aux_update_str(dup, "do", strlen(bam_get_qname(ori)) + 1, bam_get_qname(ori))) {
            print_error("markdup", "error, unable to append 'do' tag.\n");
            return -1;
        }
    }

    if (param->opt_dist) { // mark optical duplicates
        if (is_optical_duplicate(param, ori, dup, param->opt_dist, warn)) {
            bam_aux_update_str(dup, "dt", 3, "SQ");
            dup_type = 'O';
            (*optical)++;
        } else {
            // not an optical duplicate
            bam_aux_update_str(dup, "dt", 3, "LB");
        }
    }

    if ((*warn == BMD_WARNING_MAX) && (incoming_warnings != *warn)) {
        print_error("markdup", "warning, %ld decipher read name warnings.  New warnings will not be reported.\n",
                        *warn);
    }

    if (param->supp) {
        if (bam_aux_get(dup, "SA") || (dup->core.flag & BAM_FMUNMAP) || bam_aux_get(dup, "XA")) {
            char *original = NULL;

            if (param->tag) {
                original = bam_get_qname(ori);
            }

            if (add_duplicate(dup_hash, dup, original, dup_type, read_group))
                return -1;
        }
    }

    return 0;
}


/* If the duplicate type has changed to optical then retag and duplicate hash. */
static inline int optical_retag(md_param_t *param, khash_t(duplicates) *dup_hash, bam1_t *b, int paired, stats_block_t *stats) {
    int ret = 0;

    if (bam_aux_update_str(b, "dt", 3, "SQ")) {
        print_error("markdup", "error, unable to update 'dt' tag.\n");
        ret = -1;
    }

    if (paired) {
        stats->optical++;
    } else {
        stats->single_optical++;
    }

    if (param->supp) {
        // Change the duplicate type

        if (bam_aux_get(b, "SA") || (b->core.flag & BAM_FMUNMAP)
            || bam_aux_get(b, "XA")) {
            khiter_t d;

            d = kh_get(duplicates, dup_hash, bam_get_qname(b));

            if (d == kh_end(dup_hash)) {
                // error, name should already be in dup hash
                print_error("markdup", "error, duplicate name %s not found in hash.\n",
                    bam_get_qname(b));
                ret = -1;
            } else {
                kh_value(dup_hash, d).type = 'O';
            }
        }
    }

    return ret;
}


/* Check all duplicates of the highest quality read (the "original") for consistancy.  Also
   pre-calculate any values for use in check_duplicate_chain later.
   Returns 0 on success, >0 on coordinate reading error (program can continue) or
   <0 on an error (program should not continue. */
static int check_chain_against_original(md_param_t *param, khash_t(duplicates) *dup_hash, read_queue_t *ori,
             check_list_t *list, long *warn, stats_block_t *stats) {

    int ret = 0, coord_fail = 0;
    char *ori_name = bam_get_qname(ori->b);
    read_queue_t *current = ori->duplicate;
    int t_beg = 0, t_end = 0;
    long x, y;

    if (param->opt_dist) {
        coord_fail = get_coordinates(param, ori_name, &t_beg, &t_end, &x, &y, warn);
    }

    list->length = 0;

    while (current) {
        check_t *c;

        if (list->length >= list->size) {
            check_t *tmp;

            list->size *= 2;

            if (!(tmp = realloc(list->c, list->size * sizeof(check_t)))) {
                print_error("markdup", "error, Unable to expand optical check list.\n");
                return -1;
            }

            list->c = tmp;
        }

        c = &list->c[list->length];

        c->b = current->b;
        c->x = -1;
        c->y = -1;
        c->opt = 0;
        c->score = 0;
        c->mate_score = 0;
        current->dup_checked = 1;

        if (param->tag) {
            uint8_t *data;

            // at this stage all duplicates should have a do tag
            if ((data = bam_aux_get(current->b, "do")) != NULL) {
                // see if we need to change the tag
                char *old_name = bam_aux2Z(data);

                if (old_name) {
                    if (strcmp(old_name, ori_name) != 0) {
                        if (bam_aux_update_str(current->b, "do", strlen(ori_name) + 1, (const char *)ori_name)) {
                            print_error("markdup", "error, unable to update 'do' tag.\n");
                            ret =  -1;
                            break;
                        }
                    }
                } else {
                    print_error("markdup", "error, 'do' tag has wrong type for read %s.\n", bam_get_qname(current->b));
                    ret = -1;
                    break;
                }
            }
        }

        if (param->opt_dist && !coord_fail) {
            uint8_t *data;
            char *dup_type;
            int is_opt = 0;
            int current_paired = (current->b->core.flag & BAM_FPAIRED) && !(current->b->core.flag & BAM_FMUNMAP);

            if ((data = bam_aux_get(current->b, "dt"))) {
                if ((dup_type = bam_aux2Z(data))) {
                    if (strcmp(dup_type, "SQ") == 0) {
                        c->opt = 1;
                    }
                }
            }

            // need to run this to get the duplicates x and y scores
            is_opt = optical_duplicate_partial(param, ori_name, t_beg, t_end, x, y, current->b, c, param->opt_dist, warn);

            if (!c->opt && is_opt) {
                if (optical_retag(param, dup_hash, current->b, current_paired, stats)) {
                    ret = -1;
                    break;
                }

                c->opt = 1;
            }

            c->score = calc_score(current->b);

            if (current_paired) {
                if ((c->mate_score = get_mate_score(current->b)) == -1) {
                     print_error("markdup", "error, no ms score tag. Please run samtools fixmate on file first.\n");
                     ret = -1;
                     break;
                }
            }
        }

        current = current->duplicate;
        list->length++;
    }

    if (!ret && coord_fail)
        ret = coord_fail;

    ori->dup_checked = 1;

    return ret;
}


static int xcoord_sort(const void *a, const void *b) {
    check_t *ac = (check_t *) a;
    check_t *bc = (check_t *) b;

    return (ac->x - bc->x);
}


/* Check all the duplicates against each other to see if they are optical duplicates. */
static int check_duplicate_chain(md_param_t *param, khash_t(duplicates) *dup_hash, check_list_t *list,
             long *warn, stats_block_t *stats) {
    int ret = 0;
    size_t curr = 0;

    qsort(list->c, list->length, sizeof(list->c[0]), xcoord_sort);

    while (curr < list->length - 1) {
        check_t *current = &list->c[curr];
        size_t count = curr;
        char *cur_name = bam_get_qname(current->b);
        int current_paired = (current->b->core.flag & BAM_FPAIRED) && !(current->b->core.flag & BAM_FMUNMAP);

        while (++count < list->length && (list->c[count].x - current->x <= param->opt_dist)) {
            // while close enough along the x coordinate
            check_t *chk = &list->c[count];

            if (current->opt && chk->opt)
                continue;

            // if both are already optical duplicates there is no need to check again, otherwise...

            long ydiff;

            if (current->y > chk->y) {
                ydiff = current->y - chk->y;
            } else {
                ydiff = chk->y - current->y;
            }

            if (ydiff > param->opt_dist)
                continue;

            // the number are right, check the names
            if (strncmp(cur_name + current->beg, bam_get_qname(chk->b) + chk->beg, current->end - current->beg) != 0)
                continue;

            // optical duplicates
            int chk_dup = 0;
            int chk_paired = (chk->b->core.flag & BAM_FPAIRED) && !(chk->b->core.flag & BAM_FMUNMAP);

            if (current_paired != chk_paired) {
                if (!chk_paired) {
                    // chk is single vs pair, this is a dup.
                    chk_dup = 1;
                }
            } else {
                // do it by scores
                int64_t cur_score, chk_score;

                if ((current->b->core.flag & BAM_FQCFAIL) != (chk->b->core.flag & BAM_FQCFAIL)) {
                    if (current->b->core.flag & BAM_FQCFAIL) {
                        cur_score = 0;
                        chk_score = 1;
                    } else {
                        cur_score = 1;
                        chk_score = 0;
                    }
                } else {
                    cur_score = current->score;
                    chk_score = chk->score;

                    if (current_paired) {
                        // they are pairs so add mate scores.
                        chk_score += chk->mate_score;
                        cur_score += current->mate_score;
                    }
                }

                if (cur_score == chk_score) {
                    if (strcmp(bam_get_qname(chk->b), cur_name) < 0) {
                        chk_score++;
                    } else {
                        chk_score--;
                    }
                }

                if (cur_score > chk_score) {
                    chk_dup = 1;
                }
            }

            if (chk_dup) {
                // the duplicate is the optical duplicate
                if (!chk->opt) { // only change if not already an optical duplicate
                    if (optical_retag(param, dup_hash, chk->b, chk_paired, stats)) {
                        ret = -1;
                        goto fail;
                    }

                    chk->opt = 1;
                }
            } else {
                if (!current->opt) {
                    if (optical_retag(param, dup_hash, current->b, current_paired, stats)) {
                        ret = -1;
                        goto fail;
                    }

                    current->opt = 1;
                }
            }
        }

        curr++;
    }

 fail:
    return ret;
}


/* Where there is more than one duplicate go down the list and check for optical duplicates and change
   do tags (where used) to point to original (non-duplicate) read. */
static int find_duplicate_chains(md_param_t *param, read_queue_t *in_read , khash_t(duplicates) *dup_hash, check_list_t *dup_list,
                                long *warn, stats_block_t *stats) {
    int ret = 0;

    while (in_read->original) in_read = in_read->original;

    // check against the original for tagging and optical duplication
    if ((ret = check_chain_against_original(param, dup_hash, in_read, dup_list, warn, stats + in_read->read_group))) {
        if (ret < 0) { // real error
            ret = -1;
        } else { // coordinate decoding error
            ret = 0;
        }
    } else {
        // check the rest of the duplicates against each other for optical duplication
        if (param->opt_dist && check_duplicate_chain(param, dup_hash, dup_list, warn, stats + in_read->read_group)) {
            ret = -1;
        }
    }

    return ret;
}


/*
  Function to use when estimating library size.

  This is based on an approximate formula for the coverage of a set
  obtained after sampling it a given number of times with replacement.

  x = number of items in the set (the number of unique fragments in the library)

  c = number of unique items (unique read pairs observed)

  n = number of items samples (total number of read pairs)

  c and n are known; x is unknown.

  As n -> infinity, the coverage (c/x) can be given as:

  c / x = 1 - exp(-n / x)  (see https://math.stackexchange.com/questions/32800)

  This needs to be solved for x, so it is rearranged to put both terms on the
  left side and estimate_library_size() finds a value of x which gives a
  result of zero (or as close as it can get).
 */
static inline double coverage_equation(double x, double c, double n) {
    return c / x - 1 + exp(-n / x);
}


/* estimate the library size, based on the Picard code in DuplicationMetrics.java*/
static unsigned long estimate_library_size(unsigned long paired_reads, unsigned long paired_duplicate_reads, unsigned long optical) {
    unsigned long estimated_size = 0;
    unsigned long non_optical_pairs = (paired_reads - optical) / 2;
    unsigned long unique_pairs = (paired_reads - paired_duplicate_reads) / 2;
    unsigned long duplicate_pairs = (paired_duplicate_reads - optical) / 2;

    if ((non_optical_pairs && duplicate_pairs && unique_pairs) && (non_optical_pairs > duplicate_pairs)) {
        double m = 1;
        double M = 100;
        int i;

        if (coverage_equation(m * (double)unique_pairs, (double)unique_pairs, (double)non_optical_pairs) < 0) {
            print_error("markdup", "warning, unable to calculate estimated library size.\n");
            return  estimated_size;
        }

        while (coverage_equation(M * (double)unique_pairs, (double)unique_pairs, (double)non_optical_pairs) > 0) {
            M *= 10;
        }

        for (i = 0; i < 40; i++) {
            double r = (m + M) / 2;
            double u = coverage_equation(r * (double)unique_pairs, (double)unique_pairs, (double)non_optical_pairs);

            if (u > 0) {
                m = r;
            } else if (u < 0) {
                M = r;
            } else {
                break;
            }
        }

        estimated_size = (unsigned long)(unique_pairs * (m + M) / 2);
    } else {
        print_error("markdup", "warning, unable to calculate estimated library size."
                        " Read pairs %ld should be greater than duplicate pairs %ld,"
                        " which should both be non zero.\n",
                        non_optical_pairs, duplicate_pairs);
    }

    return estimated_size;
}


static void write_stats(FILE *fp, const char *title,  const char *title_con, stats_block_t *stats) {
    unsigned long els;

    els = estimate_library_size(stats->pair, stats->duplicate, stats->optical);

    if (title) {
        fprintf(fp, "%s%s\n", title, title_con);
    }

    fprintf(fp,
            "READ: %ld\n"
            "WRITTEN: %ld\n"
            "EXCLUDED: %ld\n"
            "EXAMINED: %ld\n"
            "PAIRED: %ld\n"
            "SINGLE: %ld\n"
            "DUPLICATE PAIR: %ld\n"
            "DUPLICATE SINGLE: %ld\n"
            "DUPLICATE PAIR OPTICAL: %ld\n"
            "DUPLICATE SINGLE OPTICAL: %ld\n"
            "DUPLICATE NON PRIMARY: %ld\n"
            "DUPLICATE NON PRIMARY OPTICAL: %ld\n"
            "DUPLICATE PRIMARY TOTAL: %ld\n"
            "DUPLICATE TOTAL: %ld\n"
            "ESTIMATED_LIBRARY_SIZE: %ld\n", stats->reading, stats->writing, stats->excluded, stats->examined, stats->pair, stats->single,
                            stats->duplicate, stats->single_dup, stats->optical, stats->single_optical, stats->np_duplicate, stats->np_opt_duplicate,
                            stats->single_dup + stats->duplicate, stats->single_dup + stats->duplicate + stats->np_duplicate, els);
}


static void write_json_stats(FILE *fp, const char *offset, const char *group_name, stats_block_t *stats, const char *end) {
    unsigned long  els;

    els = estimate_library_size(stats->pair, stats->duplicate, stats->optical);

    if (group_name) {
        fprintf(fp, "%s\"READ GROUP\": \"%s\",\n", offset, group_name);
    }

    fprintf(fp, "%s\"READ\": %ld,\n", offset, stats->reading);
    fprintf(fp, "%s\"WRITTEN\": %ld,\n", offset, stats->writing);
    fprintf(fp, "%s\"EXCLUDED\": %ld,\n", offset, stats->excluded);
    fprintf(fp, "%s\"EXAMINED\": %ld,\n", offset, stats->examined);
    fprintf(fp, "%s\"PAIRED\": %ld,\n", offset, stats->pair);
    fprintf(fp, "%s\"SINGLE\": %ld,\n", offset, stats->single);
    fprintf(fp, "%s\"DUPLICATE PAIR\": %ld,\n", offset, stats->duplicate);
    fprintf(fp, "%s\"DUPLICATE SINGLE\": %ld,\n", offset, stats->single_dup);
    fprintf(fp, "%s\"DUPLICATE PAIR OPTICAL\": %ld,\n", offset, stats->optical);
    fprintf(fp, "%s\"DUPLICATE SINGLE OPTICAL\": %ld,\n", offset, stats->single_optical);
    fprintf(fp, "%s\"DUPLICATE NON PRIMARY\": %ld,\n", offset, stats->np_duplicate);
    fprintf(fp, "%s\"DUPLICATE NON PRIMARY OPTICAL\": %ld,\n", offset, stats->np_opt_duplicate);
    fprintf(fp, "%s\"DUPLICATE PRIMARY TOTAL\": %ld,\n", offset, stats->single_dup + stats->duplicate);
    fprintf(fp, "%s\"DUPLICATE TOTAL\": %ld,\n", offset, stats->single_dup + stats->duplicate + stats->np_duplicate);
    fprintf(fp, "%s\"ESTIMATED_LIBRARY_SIZE\": %ld", offset, els);

    if (end) {
        fprintf(fp, "%s", end);
    }
}


/* Compare the reads near each other (coordinate sorted) and try to spot the duplicates.
   Generally the highest quality scoring is chosen as the original and all others the duplicates.
   The score is based on the sum of the quality values (<= 15) of the read and its mate (if any).
   While single reads are compared to only one read of a pair, the pair will chosen as the original.
   The comparison is done on position and orientation, see above for details.

   Marking the supplementary reads of a duplicate as also duplicates takes an extra file read/write
   step.  This is because the duplicate can occur before the primary read.*/

static int bam_mark_duplicates(md_param_t *param) {
    bam_hdr_t *header = NULL;
    khiter_t k;
    khash_t(reads) *pair_hash        = kh_init(reads);
    khash_t(reads) *single_hash      = kh_init(reads);
    klist_t(read_queue) *read_buffer = kl_init(read_queue);
    kliter_t(read_queue) *rq;
    khash_t(duplicates) *dup_hash    = kh_init(duplicates);
    khash_t(read_groups) *rg_hash    = kh_init(read_groups);
    int32_t prev_tid;
    hts_pos_t prev_coord;
    read_queue_t *in_read;
    int ret;
    stats_block_t *stats, *stat_array = NULL;
    int num_groups = 0;
    long opt_warnings = 0, bc_warnings = 0;
    tmp_file_t temp;
    char *idx_fn = NULL;
    int exclude = 0;
    check_list_t dup_list = {NULL, 0, 0};

    if (!pair_hash || !single_hash || !read_buffer || !dup_hash || !rg_hash) {
        print_error("markdup", "error, unable to allocate memory to initialise structures.\n");
        goto fail;
    }

    if ((header = sam_hdr_read(param->in)) == NULL) {
        print_error("markdup", "error reading header\n");
        goto fail;
    }

    // accept unknown, unsorted or coordinate sort order, but error on queryname sorted.
    // only really works on coordinate sorted files.
    kstring_t str = KS_INITIALIZE;
    if (!sam_hdr_find_tag_hd(header, "SO", &str) && str.s && !strcmp(str.s, "queryname")) {
        print_error("markdup", "error, queryname sorted, must be sorted by coordinate.\n");
        ks_free(&str);
        goto fail;
    }
    ks_free(&str);

    if (!param->no_pg && sam_hdr_add_pg(header, "samtools", "VN", samtools_version(),
                        param->arg_list ? "CL" : NULL,
                        param->arg_list ? param->arg_list : NULL,
                        NULL) != 0) {
        print_error("markdup", "warning, unable to add @PG line to header.\n");
    }

    if (sam_hdr_write(param->out, header) < 0) {
        print_error("markdup", "error writing header.\n");
        goto fail;
    }
    if (param->write_index) {
        if (!(idx_fn = auto_index(param->out, param->out_fn, header)))
            goto fail;
    }

    if (param->read_groups) {
        num_groups = sam_hdr_count_lines(header, "RG");
        int g_ret = 0;

        if (num_groups > 0) {
            int i;

            for (i = 0; i < num_groups; i++) {
                const char *rg_key;
                khiter_t rg;

                rg_key = sam_hdr_line_name(header, "RG", i);

                if (rg_key) {
                    rg = kh_get(read_groups, rg_hash, rg_key);

                    if (rg == kh_end(rg_hash)) { // new entry
                        rg = kh_put(read_groups, rg_hash, rg_key, &g_ret);

                        if (g_ret > 0) {
                            kh_value(rg_hash, rg) = i + 1;
                        } else {
                            print_error("markdup", "error, unable to populate read group ids.  "
                                     "Read groups will not be used\n");
                            g_ret = -1;
                            break;
                        }
                    } else {
                        print_error("markdup", "error, duplicate read group ids %s."
                                  "Read groups will not be used\n", rg_key);
                        g_ret = -1;
                        break;
                    }
                } else {
                    print_error("markdup", "error, Unable to retrieve read group at position %d."
                              "Read groups will not be used\n", i);
                    g_ret = -1;
                    break;
                }
            }
        } else {
            print_error("markdup", "error, no read groups found.\n");
            g_ret = -1;
        }

        if (g_ret < 0) {
            print_error("markdup", "error, read groups will not be used.\n");
            param->read_groups = 0;
            num_groups = 0;
        }
    }

    // stat_array[0] will be for ungrouped reads
    stat_array = calloc(num_groups + 1, sizeof(stats_block_t));

    if (stat_array == NULL) {
        print_error("markdup", "error, unable to allocate memory for stats.\n");
        goto fail;
    }

    // used for coordinate order checks
    prev_tid = prev_coord = 0;

    // get the buffer going
    in_read = kl_pushp(read_queue, read_buffer);
    if (!in_read) {
        print_error("markdup", "error, unable to allocate memory to hold reads.\n");
        goto fail;
    }

    // handling supplementary reads needs a temporary file
    if (param->supp) {
        if (tmp_file_open_write(&temp, param->prefix, 1)) {
            print_error("markdup", "error, unable to open tmp file %s.\n", param->prefix);
            goto fail;
        }
    }

    if ((in_read->b = bam_init1()) == NULL) {
        print_error("markdup", "error, unable to allocate memory for alignment.\n");
        goto fail;
    }

    if (param->check_chain && !(param->tag || param->opt_dist))
        param->check_chain = 0;

    if (param->check_chain) {
        dup_list.size = 128;
        dup_list.c = NULL;

        if ((dup_list.c = malloc(dup_list.size * sizeof(check_t))) == NULL) {
            print_error("markdup", "error, unable to allocate memory for dup_list.\n");
            goto fail;
        }
    }

    while ((ret = sam_read1(param->in, header, in_read->b)) >= 0) {

        // do some basic coordinate order checks
        if (in_read->b->core.tid >= 0) { // -1 for unmapped reads
            if (in_read->b->core.tid < prev_tid ||
               ((in_read->b->core.tid == prev_tid) && (in_read->b->core.pos < prev_coord))) {
                print_error("markdup", "error, not in coordinate sorted order.\n");
                goto fail;
            }
        }

        prev_coord = in_read->pos = in_read->b->core.pos;
        prev_tid   =  in_read->b->core.tid;
        in_read->pair_key.single   = 1;
        in_read->single_key.single = 0;
        in_read->duplicate = NULL;
        in_read->original = NULL;
        in_read->dup_checked = 0;
        in_read->read_group = 0;

        if (param->read_groups) {
            uint8_t *data;
            char *rg;

            if ((data = bam_aux_get(in_read->b, "RG"))) {
                if ((rg = bam_aux2Z(data))) {
                    khiter_t r;

                    r = kh_get(read_groups, rg_hash, rg);

                    if (r != kh_end(rg_hash)) {
                        in_read->read_group = kh_value(rg_hash, r);
                    }
                }
            }
        }

        stats = stat_array + in_read->read_group;

        stats->reading++;

        if (param->clear && (in_read->b->core.flag & BAM_FDUP)) {
            uint8_t *data;

            in_read->b->core.flag ^= BAM_FDUP;

            if ((data = bam_aux_get(in_read->b, "dt")) != NULL) {
                bam_aux_del(in_read->b, data);
            }

            if ((data = bam_aux_get(in_read->b, "do")) != NULL) {
                bam_aux_del(in_read->b, data);
            }
        }

        if (param->include_fails) {
            exclude |= (BAM_FSECONDARY | BAM_FSUPPLEMENTARY | BAM_FUNMAP);
        } else {
            exclude |= (BAM_FSECONDARY | BAM_FSUPPLEMENTARY | BAM_FUNMAP | BAM_FQCFAIL);
        }

        // read must not be secondary, supplementary, unmapped or (possibly) failed QC
        if (!(in_read->b->core.flag & exclude)) {
            stats->examined++;


            // look at the pairs first
            if ((in_read->b->core.flag & BAM_FPAIRED) && !(in_read->b->core.flag & BAM_FMUNMAP)) {
                int ret, mate_tmp;
                key_data_t pair_key;
                key_data_t single_key;
                in_hash_t *bp;

                if (make_pair_key(param, &pair_key, in_read->b, in_read->read_group, &bc_warnings)) {
                    print_error("markdup", "error, unable to assign pair hash key.\n");
                    goto fail;
                }

                make_single_key(param, &single_key, in_read->b, in_read->read_group, &bc_warnings);

                stats->pair++;
                in_read->pos = single_key.this_coord; // cigar/orientation modified pos

                // put in singles hash for checking against non paired reads
                k = kh_put(reads, single_hash, single_key, &ret);

                if (ret > 0) { // new
                    // add to single duplicate hash
                    bp = &kh_val(single_hash, k);
                    bp->p = in_read;
                    in_read->single_key = single_key;
                } else if (ret == 0) { // exists
                    // look at singles only for duplication marking
                    bp = &kh_val(single_hash, k);

                    if (!(bp->p->b->core.flag & BAM_FPAIRED) || (bp->p->b->core.flag & BAM_FMUNMAP)) {
                       // singleton will always be marked duplicate even if
                       // scores more than one read of the pair
                        bam1_t *dup = bp->p->b;

                        if (param->check_chain) {
                            in_read->duplicate = bp->p;
                            bp->p->original = in_read;
                        }

                        bp->p = in_read;

                        if (mark_duplicates(param, dup_hash, bp->p->b, dup, in_read->read_group, &stats->single_optical, &opt_warnings))
                            goto fail;

                        stats->single_dup++;
                    }
                } else {
                    print_error("markdup", "error, single hashing failure for paired read.\n");
                    goto fail;
                }

                // now do the pair
                k = kh_put(reads, pair_hash, pair_key, &ret);

                if (ret > 0) { // new
                    // add to the pair hash
                    bp = &kh_val(pair_hash, k);
                    bp->p = in_read;
                    in_read->pair_key = pair_key;
                } else if (ret == 0) {
                    int64_t old_score, new_score, tie_add = 0;
                    bam1_t *dup = NULL;

                    bp = &kh_val(pair_hash, k);

                    if ((bp->p->b->core.flag & BAM_FQCFAIL) != (in_read->b->core.flag & BAM_FQCFAIL)) {
                        if (bp->p->b->core.flag & BAM_FQCFAIL) {
                            old_score = 0;
                            new_score = 1;
                        } else {
                            old_score = 1;
                            new_score = 0;
                        }
                    } else {
                        if ((mate_tmp = get_mate_score(bp->p->b)) == -1) {
                            print_error("markdup", "error, no ms score tag. Please run samtools fixmate on file first.\n");
                            goto fail;
                        } else {
                            old_score = calc_score(bp->p->b) + mate_tmp;
                        }

                        if ((mate_tmp = get_mate_score(in_read->b)) == -1) {
                            print_error("markdup", "error, no ms score tag. Please run samtools fixmate on file first.\n");
                            goto fail;
                        } else {
                            new_score = calc_score(in_read->b) + mate_tmp;
                        }
                    }

                    // choose the highest score as the original
                    // and add it to the pair hash, mark the other as duplicate

                    if (new_score == old_score) {
                        if (strcmp(bam_get_qname(in_read->b), bam_get_qname(bp->p->b)) < 0) {
                            tie_add = 1;
                        } else {
                            tie_add = -1;
                        }
                    }

                    if (new_score + tie_add > old_score) { // swap reads
                        dup = bp->p->b;

                        if (param->check_chain) {

                            if (in_read->duplicate) {
                                read_queue_t *current = in_read->duplicate;

                                while (current->duplicate) {
                                    current = current->duplicate;
                                }

                                current->duplicate = bp->p;
                            } else {
                                in_read->duplicate = bp->p;
                            }

                            bp->p->original = in_read;
                        }

                        bp->p = in_read;
                    } else {
                        if (param->check_chain) {
                            if (bp->p->duplicate) {
                                if (in_read->duplicate) {
                                    read_queue_t *current = bp->p->duplicate;

                                    while (current->duplicate) {
                                        current = current->duplicate;
                                    }

                                    current->duplicate = in_read->duplicate;
                                }

                                in_read->duplicate = bp->p->duplicate;
                            }

                            bp->p->duplicate = in_read;
                            in_read->original = bp->p;
                        }

                        dup = in_read->b;
                    }

                    if (mark_duplicates(param, dup_hash, bp->p->b, dup, in_read->read_group, &stats->optical, &opt_warnings))
                        goto fail;

                    stats->duplicate++;
                } else {
                    print_error("markdup", "error, pair hashing failure.\n");
                    goto fail;
                }
            } else { // do the single (or effectively single) reads
                int ret;
                key_data_t single_key;
                in_hash_t *bp;

                make_single_key(param, &single_key, in_read->b, in_read->read_group, &bc_warnings);

                stats->single++;
                in_read->pos = single_key.this_coord; // cigar/orientation modified pos

                k = kh_put(reads, single_hash, single_key, &ret);

                if (ret > 0) { // new
                    bp = &kh_val(single_hash, k);
                    bp->p = in_read;
                    in_read->single_key = single_key;
                } else if (ret == 0) { // exists
                    bp = &kh_val(single_hash, k);

                    if ((bp->p->b->core.flag & BAM_FPAIRED) && !(bp->p->b->core.flag & BAM_FMUNMAP)) {
                        // if matched against one of a pair just mark as duplicate

                        if (param->check_chain) {
                            if (bp->p->duplicate) {
                                in_read->duplicate = bp->p->duplicate;
                            }

                            bp->p->duplicate = in_read;
                            in_read->original = bp->p;
                        }

                        if (mark_duplicates(param, dup_hash, bp->p->b, in_read->b, in_read->read_group, &stats->single_optical, &opt_warnings))
                            goto fail;

                    } else {
                        int64_t old_score, new_score;
                        bam1_t *dup = NULL;

                        old_score = calc_score(bp->p->b);
                        new_score = calc_score(in_read->b);

                        // choose the highest score as the original, add it
                        // to the single hash and mark the other as duplicate
                        if (new_score > old_score) { // swap reads
                            dup = bp->p->b;

                            if (param->check_chain) {
                                in_read->duplicate = bp->p;
                                bp->p->original = in_read;
                            }

                            bp->p = in_read;
                        } else {
                            if (param->check_chain) {
                                if (bp->p->duplicate) {
                                    in_read->duplicate = bp->p->duplicate;
                                }

                                bp->p->duplicate = in_read;
                                in_read->original = bp->p;
                            }

                            dup = in_read->b;
                        }

                        if (mark_duplicates(param, dup_hash, bp->p->b, dup, in_read->read_group, &stats->single_optical, &opt_warnings))
                            goto fail;
                    }

                    stats->single_dup++;
                } else {
                    print_error("markdup", "error, single hashing failure for single read.\n");
                    goto fail;
                }
            }
        } else {
            stats->excluded++;
        }

        // loop through the stored reads and write out those we
        // no longer need
        rq = kl_begin(read_buffer);
        while (rq != kl_end(read_buffer)) {
            in_read = &kl_val(rq);

            /* keep a moving window of reads based on coordinates and max read length.  Any unaligned reads
               should just be written as they cannot be matched as duplicates. */
            if (in_read->pos + param->max_length > prev_coord && in_read->b->core.tid == prev_tid && (prev_tid != -1 || prev_coord != -1)) {
                break;
            }

            if (param->check_chain && !in_read->dup_checked && (in_read->original || in_read->duplicate)) {
                if (find_duplicate_chains(param, in_read, dup_hash, &dup_list, &opt_warnings, stat_array)) {
                    print_error("markdup", "error, duplicate checking failed.\n");
                    goto fail;
                }
            }

            if (!param->remove_dups || !(in_read->b->core.flag & BAM_FDUP)) {
                if (param->supp) {
                    if (tmp_file_write(&temp, in_read->b)) {
                        print_error("markdup", "error, writing temp output failed.\n");
                        goto fail;
                    }
                } else {
                    if (sam_write1(param->out, header, in_read->b) < 0) {
                        print_error("markdup", "error, writing output failed.\n");
                        goto fail;
                    }
                }

                stat_array[in_read->read_group].writing++;
            }

            // remove from hash
            if (in_read->pair_key.single == 0) {
                k = kh_get(reads, pair_hash, in_read->pair_key);
                kh_del(reads, pair_hash, k);
            }

            if (in_read->single_key.single == 1) {
                k = kh_get(reads, single_hash, in_read->single_key);
                kh_del(reads, single_hash, k);
            }

            kl_shift(read_queue, read_buffer, NULL);
            bam_destroy1(in_read->b);
            rq = kl_begin(read_buffer);
        }

        // set the next one up for reading
        in_read = kl_pushp(read_queue, read_buffer);
        if (!in_read) {
            print_error("markdup", "error, unable to allocate memory for read in queue.\n");
            goto fail;
        }

        if ((in_read->b = bam_init1()) == NULL) {
            print_error("markdup", "error, unable to allocate memory for alignment.\n");
            goto fail;
        }
    }

    if (ret < -1) {
        print_error("markdup", "error, truncated input file.\n");
        goto fail;
    }

    // write out the end of the list
    rq = kl_begin(read_buffer);
    while (rq != kl_end(read_buffer)) {
        in_read = &kl_val(rq);

        if (bam_get_qname(in_read->b)) { // last entry will be blank
            if (param->check_chain && !in_read->dup_checked && (in_read->original || in_read->duplicate)) {
                if (find_duplicate_chains(param, in_read, dup_hash, &dup_list, &opt_warnings, stat_array)) {
                    print_error("markdup", "error, duplicate checking failed.\n");
                    goto fail;
                }
            }

            if (!param->remove_dups || !(in_read->b->core.flag & BAM_FDUP)) {
                if (param->supp) {
                    if (tmp_file_write(&temp, in_read->b)) {
                        print_error("markdup", "error, writing temp output failed on final write.\n");
                        goto fail;
                    }
                } else {
                    if (sam_write1(param->out, header, in_read->b) < 0) {
                        print_error("markdup", "error, writing output failed on final write.\n");
                        goto fail;
                    }
                }

                stat_array[in_read->read_group].writing++;
            }
        }

        kl_shift(read_queue, read_buffer, NULL);
        bam_destroy1(in_read->b);
        rq = kl_begin(read_buffer);
    }

    if (param->supp) {
        bam1_t *b;

        if (tmp_file_end_write(&temp)) {
            print_error("markdup", "error, unable to end tmp writing.\n");
            goto fail;
        }

        // read data from temp file and mark duplicate supplementary alignments

        if (tmp_file_begin_read(&temp)) {
            goto fail;
        }

        b = bam_init1();

        while ((ret = tmp_file_read(&temp, b)) > 0) {

            if ((b->core.flag & BAM_FSUPPLEMENTARY) || (b->core.flag & BAM_FUNMAP) || (b->core.flag & BAM_FSECONDARY)) {

                k = kh_get(duplicates, dup_hash, bam_get_qname(b));

                if (k != kh_end(dup_hash)) {

                    b->core.flag  |= BAM_FDUP;
                    stat_array[kh_val(dup_hash, k).read_group].np_duplicate++;

                    if (param->tag && kh_val(dup_hash, k).name) {
                        if (bam_aux_update_str(b, "do", strlen(kh_val(dup_hash, k).name) + 1, (char*)kh_val(dup_hash, k).name)) {
                            print_error("markdup", "error, unable to append supplementary 'do' tag.\n");
                            goto fail;
                        }
                    }

                    if (param->opt_dist) {
                        if (kh_val(dup_hash, k).type) {
                            bam_aux_update_str(b, "dt", 3, "SQ");
                            stat_array[kh_val(dup_hash, k).read_group].np_opt_duplicate++;
                        } else {
                            bam_aux_update_str(b, "dt", 3, "LB");
                        }
                    }
                }
            }

            if (!param->remove_dups || !(b->core.flag & BAM_FDUP)) {
                if (sam_write1(param->out, header, b) < 0) {
                    print_error("markdup", "error, writing final output failed.\n");
                    goto fail;
                }
            }
        }

        if (ret == -1) {
            print_error("markdup", "error, failed to read tmp file.\n");
            goto fail;
        }

        for (k = kh_begin(dup_hash); k != kh_end(dup_hash); ++k) {
            if (kh_exist(dup_hash, k)) {
                free(kh_val(dup_hash, k).name);
                free((char *)kh_key(dup_hash, k));
                kh_key(dup_hash, k) = NULL;
            }
        }

        tmp_file_destroy(&temp);
        bam_destroy1(b);
    }

    if (opt_warnings) {
        print_error("markdup", "warning, number of failed attempts to get coordinates from read names = %ld\n",
                        opt_warnings);
    }

    if (bc_warnings) {
        print_error("markdup", "warning, number of failed attempts to get barcodes = %ld\n", bc_warnings);
    }

    if (param->do_stats) {
        FILE *fp;
        int file_open = 0;
        stats_block_t total;
        int i;

        if (param->stats_file) {
            if (NULL == (fp = fopen(param->stats_file, "w"))) {
                print_error("markdup", "warning, cannot write stats to %s.\n", param->stats_file);
                fp = stderr;
            } else {
                file_open = 1;
            }
        } else {
            fp = stderr;
        }

        total = stat_array[0];

        if (param->read_groups) {
            for (i = 1; i <= num_groups; i++) {
                total.reading += stat_array[i].reading;
                total.writing += stat_array[i].writing;
                total.excluded += stat_array[i].excluded;
                total.duplicate += stat_array[i].duplicate;
                total.single += stat_array[i].single;
                total.pair += stat_array[i].pair;
                total.single_dup += stat_array[i].single_dup;
                total.examined += stat_array[i].examined;
                total.optical += stat_array[i].optical;
                total.single_optical += stat_array[i].single_optical;
                total.np_duplicate += stat_array[i].np_duplicate;
                total.np_opt_duplicate += stat_array[i].np_opt_duplicate;
            }
        }

        if (!param->json) {
            write_stats(fp, "COMMAND: ", param->arg_list, &total);
            fprintf(fp, "\n");

            if (param->read_groups) {
                if (stat_array[0].reading) {
                    write_stats(fp, "READ GROUP: ", "ungrouped", stat_array);
                    fprintf(fp, "\n");
                }

                for (i = 0; i < num_groups; i++) {
                    write_stats(fp, "READ GROUP: ", sam_hdr_line_name(header, "RG", i), stat_array + i + 1);
                    fprintf(fp, "\n");
                }
            }
        } else {
            char space4[]  = "    ";
            char space8[]  = "        ";
            char space12[] = "            ";

            fprintf(fp, "{\n");
            fprintf(fp, "%s\"COMMAND\": \"%s\",\n", space4, param->arg_list);
            write_json_stats(fp, space4, NULL, &total, param->read_groups ? ",\n" : "\n");

            if (param->read_groups) {
                fprintf(fp, "%s\"READ GROUPS\": [\n", space4);

                if (stat_array[0].reading) {
                    fprintf(fp, "%s{\n", space8);
                    write_json_stats(fp, space12, "ungrouped", stat_array, "\n");
                    fprintf(fp, "%s},\n", space8);
                }

                for (i = 0; i < num_groups; i++) {
                    fprintf(fp, "%s{\n", space8);

                    write_json_stats(fp, space12,  sam_hdr_line_name(header, "RG", i), stat_array + i + 1, "\n");

                    if (i < num_groups -1 ) {
                        fprintf(fp, "%s},\n", space8);
                    } else {
                        fprintf(fp, "%s}\n", space8);
                    }
                }

                fprintf(fp, "%s]\n", space4);
            }

            fprintf(fp, "}\n");
        }

        if (file_open) {
            fclose(fp);
        }
    }

    if (param->write_index) {
        if (sam_idx_save(param->out) < 0) {
            print_error_errno("markdup", "error, writing index failed");
            goto fail;
        }
    }

    if (param->check_chain && (param->tag || param->opt_dist))
        free(dup_list.c);

    free(stat_array);
    kh_destroy(reads, pair_hash);
    kh_destroy(reads, single_hash);
    kl_destroy(read_queue, read_buffer);
    kh_destroy(duplicates, dup_hash);
    kh_destroy(read_groups, rg_hash);
    sam_hdr_destroy(header);

    return 0;

 fail:
    for (rq = kl_begin(read_buffer); rq != kl_end(read_buffer); rq = kl_next(rq))
        bam_destroy1(kl_val(rq).b);
    kl_destroy(read_queue, read_buffer);

    for (k = kh_begin(dup_hash); k != kh_end(dup_hash); ++k) {
        if (kh_exist(dup_hash, k)) {
            free((char *)kh_key(dup_hash, k));
        }
    }
    kh_destroy(duplicates, dup_hash);
    kh_destroy(read_groups, rg_hash);

    if (param->check_chain && (param->tag || param->opt_dist))
        free(dup_list.c);

    free(stat_array);
    kh_destroy(reads, pair_hash);
    kh_destroy(reads, single_hash);
    sam_hdr_destroy(header);
    return 1;
}


static int markdup_usage(void) {
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage:  samtools markdup <input.bam> <output.bam>\n\n");
    fprintf(stderr, "Option: \n");
    fprintf(stderr, "  -r                 Remove duplicate reads\n");
    fprintf(stderr, "  -l INT             Max read length (default 300 bases)\n");
    fprintf(stderr, "  -S                 Mark supplementary alignments of duplicates as duplicates (slower).\n");
    fprintf(stderr, "  -s                 Report stats.\n");
    fprintf(stderr, "  -f NAME            Write stats to named file.  Implies -s.\n");
    fprintf(stderr, "  --json             Output stats in JSON.  Also implies -s\n");
    fprintf(stderr, "  -T PREFIX          Write temporary files to PREFIX.samtools.nnnn.nnnn.tmp.\n");
    fprintf(stderr, "  -d INT             Optical distance (if set, marks with dt tag)\n");
    fprintf(stderr, "  -c                 Clear previous duplicate settings and tags.\n");
    fprintf(stderr, "  -m --mode TYPE     Duplicate decision method for paired reads.\n"
                    "                     TYPE = t measure positions based on template start/end (default).\n"
                    "                            s measure positions based on sequence start.\n");
    fprintf(stderr, "  -u                 Output uncompressed data\n");
    fprintf(stderr, "  --include-fails    Include quality check failed reads.\n");
    fprintf(stderr, "  --no-PG            Do not add a PG line\n");
    fprintf(stderr, "  --no-multi-dup     Reduced duplicates of duplicates checking.\n");
    fprintf(stderr, "  --read-coords STR  Regex for coords from read name.\n");
    fprintf(stderr, "  --coords-order STR Order of regex elements. txy (default).  With t being a part of\n"
                    "                     the read names that must be equal and x/y being coordinates.\n");
    fprintf(stderr, "  --barcode-tag STR  Use barcode a tag that duplicates much match.\n");
    fprintf(stderr, "  --barcode-name     Use the UMI/barcode in the read name (eigth colon delimited part).\n");
    fprintf(stderr, "  --barcode-rgx STR  Regex for barcode in the readname (alternative to --barcode-name).\n");
    fprintf(stderr, "  --use-read-groups  Use the read group tags in duplicate matching.\n");
    fprintf(stderr, "  -t                 Mark primary duplicates with the name of the original in a \'do\' tag."
                                        " Mainly for information and debugging.\n");

    sam_global_opt_help(stderr, "-.O..@..");

    fprintf(stderr, "\nThe input file must be coordinate sorted and must have gone"
                     " through fixmates with the mate scoring option on.\n");

    return 1;
}


int bam_markdup(int argc, char **argv) {
    int c, ret, bc_name = 0;
    char wmode[4] = {'w', 0, 0, 0};
    sam_global_args ga = SAM_GLOBAL_ARGS_INIT;
    htsThreadPool p = {NULL, 0};
    kstring_t tmpprefix = {0, 0, NULL};
    struct stat st;
    unsigned int t;
    char *regex = NULL, *bc_regex = NULL;
    char *regex_order = "txy";
    md_param_t param = {NULL, NULL, NULL, 0, 300, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                        1, NULL, NULL, NULL, NULL, 0, 0, 0, NULL, NULL, 0, 0};

    static const struct option lopts[] = {
        SAM_OPT_GLOBAL_OPTIONS('-', 0, 'O', 0, 0, '@'),
        {"include-fails", no_argument, NULL, 1001},
        {"no-PG", no_argument, NULL, 1002},
        {"mode", required_argument, NULL, 'm'},
        {"no-multi-dup", no_argument, NULL, 1003},
        {"read-coords", required_argument, NULL, 1004},
        {"coords-order", required_argument, NULL, 1005},
        {"barcode-tag", required_argument, NULL, 1006},
        {"barcode-name", no_argument, NULL, 1007},
        {"barcode-rgx", required_argument, NULL, 1008},
        {"use-read-groups", no_argument, NULL, 1009},
        {"json", no_argument, NULL, 1010},
        {NULL, 0, NULL, 0}
    };

    while ((c = getopt_long(argc, argv, "rsl:StT:O:@:f:d:cm:u", lopts, NULL)) >= 0) {
        switch (c) {
            case 'r': param.remove_dups = 1; break;
            case 'l': param.max_length = atoi(optarg); break;
            case 's': param.do_stats = 1; break;
            case 'T': kputs(optarg, &tmpprefix); break;
            case 'S': param.supp = 1; break;
            case 't': param.tag = 1; break;
            case 'f': param.stats_file = optarg; param.do_stats = 1; break;
            case 'd': param.opt_dist = atoi(optarg); break;
            case 'c': param.clear = 1; break;
            case 'm':
                if (strcmp(optarg, "t") == 0) {
                    param.mode = MD_MODE_TEMPLATE;
                } else if (strcmp(optarg, "s") == 0) {
                    param.mode = MD_MODE_SEQUENCE;
                } else {
                    print_error("markdup", "error, unknown mode '%s'.\n", optarg);
                    return markdup_usage();
                }

                break;
            case 'u': wmode[1] = '0'; break;
            case 1001: param.include_fails = 1; break;
            case 1002: param.no_pg = 1; break;
            case 1003: param.check_chain = 0; break;
            case 1004: regex = optarg; break;
            case 1005: regex_order = optarg; break;
            case 1006: param.barcode = optarg; break;
            case 1007: bc_name = 1; break;
            case 1008: bc_name = 1, bc_regex = optarg; break;
            case 1009: param.read_groups = 1; break;
            case 1010: param.json = 1; param.do_stats = 1; break;
            default: if (parse_sam_global_opt(c, optarg, lopts, &ga) == 0) break;
            /* else fall-through */
            case '?': return markdup_usage();
        }
    }

    if (optind + 2 > argc)
        return markdup_usage();

    if (param.barcode && bc_name) {
        print_error("markdup", "error, cannot specify --barcode-tag and "
                        "--barcode-name (or --barcode-rgx) at same time.\n");
        return 1;
    }

    if (param.opt_dist < 0) param.opt_dist = 0;
    if (param.max_length < 0) param.max_length = 300;

    if (regex) {
        int result;

        // set the order the elements of the regex are assigned to.
        // x and y being coordinates, t being any other important part of the read
        // e.g. tile and lane
        // x and y order does not matter as long as it is consistent

        if ((strncmp(regex_order, "txy", 3) == 0) || (strncmp(regex_order, "tyx", 3) == 0)) {
            param.rgx_t = 1;
            param.rgx_x = 2;
            param.rgx_y = 3;
        } else if ((strncmp(regex_order, "xyt", 3) == 0) || (strncmp(regex_order, "yxt", 3) == 0)) {
            param.rgx_x = 1;
            param.rgx_y = 2;
            param.rgx_t = 3;
        } else if ((strncmp(regex_order, "xty", 3) == 0) || (strncmp(regex_order, "ytx", 3) == 0)) {
            param.rgx_x = 1;
            param.rgx_t = 2;
            param.rgx_y = 3;
        } else if ((strncmp(regex_order, "xy", 2) == 0) || (strncmp(regex_order, "yx", 2) == 0)) {
            param.rgx_x = 1;
            param.rgx_y = 2;
            param.rgx_t = 0;
        } else {
            print_error("markdup", "error,  could not recognise regex coordinate order \"%s\".\n", regex_order);
            return 1;
        }

        if ((param.rgx = malloc(sizeof(regex_t))) == NULL) {
            print_error("markdup", "error,  could not allocate memory for regex.\n");
            return 1;
        }

        if ((result = regcomp(param.rgx, regex, REG_EXTENDED))) {
            char err_msg[256];

            regerror(result, param.rgx, err_msg, 256);
            print_error("markdup", "error, regex fail \"%s\"\n", err_msg);
            free(param.rgx);
            return 1;
        }
    }

    if (bc_name) {
        int result;

        /* From Illumina UMI documentation: "The UMI sequence is located in the
           eighth colon-delimited field of the read name (QNAME)". */
        char *rgx = "[0-9A-Za-z]+:[0-9A-Za-z]+:[0-9A-Za-z]+:[0-9A-Za-z]+:[0-9A-Za-z]+:[0-9A-Za-z]+:[0-9A-Za-z]+:([!-?A-~]+)";

        if ((param.bc_rgx = malloc(sizeof(regex_t))) == NULL) {
            print_error("markdup", "error,  could not allocate memory for barcode regex.\n");
            return 1;
        }

        if (bc_regex) {
            rgx = bc_regex;
        }

        if ((result = regcomp(param.bc_rgx, rgx, REG_EXTENDED))) {
            char err_msg[256];

            regerror(result, param.bc_rgx, err_msg, 256);
            print_error("markdup", "error, barcode regex fail \"%s\"\n", err_msg);
            free(param.bc_rgx);
            return 1;
        }
    }

    param.in = sam_open_format(argv[optind], "r", &ga.in);

    if (!param.in) {
        print_error_errno("markdup", "error, failed to open \"%s\" for input", argv[optind]);
        return 1;
    }

    strcat(wmode, "b"); // default if unknown suffix
    sam_open_mode(wmode + strlen(wmode)-1, argv[optind + 1], NULL);
    param.out = sam_open_format(argv[optind + 1], wmode, &ga.out);

    if (!param.out) {
        print_error_errno("markdup", "error, failed to open \"%s\" for output", argv[optind + 1]);
        return 1;
    }

    if (ga.nthreads > 0)  {
        if (!(p.pool = hts_tpool_init(ga.nthreads))) {
            print_error("markdup", "error creating thread pool.\n");
            return 1;
        }

        hts_set_opt(param.in,  HTS_OPT_THREAD_POOL, &p);
        hts_set_opt(param.out, HTS_OPT_THREAD_POOL, &p);
    }

    // actual stuff happens here

    // we need temp files so fix up the name here
    if (tmpprefix.l == 0) {

        if (strcmp(argv[optind + 1], "-") != 0)
            ksprintf(&tmpprefix, "%s.", argv[optind + 1]);
        else
            kputc('.', &tmpprefix);
    }

    if (stat(tmpprefix.s, &st) == 0 && S_ISDIR(st.st_mode)) {
        if (tmpprefix.s[tmpprefix.l-1] != '/') kputc('/', &tmpprefix);
    }

    t = ((unsigned) time(NULL)) ^ ((unsigned) clock());
    ksprintf(&tmpprefix, "samtools.%d.%u.tmp", (int) getpid(), t % 10000);
    param.prefix = tmpprefix.s;

    param.arg_list = stringify_argv(argc + 1, argv - 1);
    param.write_index = ga.write_index;
    param.out_fn = argv[optind + 1];

    ret = bam_mark_duplicates(&param);

    sam_close(param.in);

    if (sam_close(param.out) < 0) {
        print_error("markdup", "error closing output file.\n");
        ret = 1;
    }

    if (p.pool) hts_tpool_destroy(p.pool);

    if (param.rgx) {
        regfree(param.rgx);
        free(param.rgx);
    }

    if (param.bc_rgx) {
        regfree(param.bc_rgx);
        free(param.bc_rgx);
    }

    free(param.arg_list);
    free(tmpprefix.s);
    sam_global_args_free(&ga);

    return ret;
}
