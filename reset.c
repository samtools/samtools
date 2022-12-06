/*  reset.c --  removes aligner updates and reference data from input sam /
                bam / cram file and makes read data raw for new processing

    Copyright (C) 2022, 2023 Genome Research Ltd.

    Author: Vasudeva Sarma <vasudeva.sarma@sanger.ac.uk>

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

*/

#include "samtools.h"
#include "htslib/sam.h"
#include "sam_opts.h"
#include "htslib/thread_pool.h"
#include "htslib/khash.h"
#include "sam_utils.h"
#include <unistd.h>

#define TAGNUM(X) (((X)[0] << 8) | (X)[1])  //to create key for aux tags, like type key in htslib
#define LONG_OPT(X) (128 + (X))             //to handle long and short options with same char

typedef struct conf_data
{
    int keepRGs;                    //RG line handling
    int noPGentry;                  //PG line for reset op or not
    auxhash_t aux_keep;             //SET that holds the aux tags to be retained
    auxhash_t aux_remove;           //SET that holds the aux tags to be removed
    char *pgid;                     //PG id onwards which to be removed
} conf_data;

/// usage - print the usage
/** @param fp pointer to the file / terminal to which usage to be dumped
returns nothing
*/
static void usage(FILE *fp)
{
    fprintf(fp, "Usage: samtools reset [options]\n\
  -o FILE      Output file\n\
  -x, --remove-tag STR\n\
               Aux tags to be removed\n\
      --keep-tag STR\n\
               Aux tags to be retained. Equivalent to -x ^STR\n\
      --reject-PG ID\n\
               Removes PG line with ID matching to input and succeeding PG lines\n\
      --no-RG  To have RG lines or not\n\
      --no-PG  To have PG entry or not for reset operation\n");

    sam_global_opt_help(fp, "--O--@--");
    return;
}

/// removeauxtags - remove aux tags in bam data which are not present in acceptable tag set
/** @param bamdata - pointer to the bamdata from which needs the filtering
 *  @param config - pointer to conf_data
returns nothing
*/
void removeauxtags(bam1_t *bamdata, conf_data *config)
{
    uint8_t *auxdata = NULL;
    const char *tag = NULL, rg[] = "RG";
    khint_t iter = 0;
    int ret = 0;

    if (!bamdata || !config || (!config->aux_keep && !config->aux_remove && config->keepRGs))
        return;

    //remove RG tags from bamdata if keepRG is false
    if (!config->keepRGs) {
        if (!config->aux_keep && !config->aux_remove) {
            //none of aux tag filter in use, create remove filter
            config->aux_remove = kh_init(aux_exists);
        }

        if (config->aux_keep) {
            //keep set in use, remove RG if present
            iter = kh_get(aux_exists, config->aux_keep, TAGNUM(rg));
            if (iter != kh_end(config->aux_keep)) {
                kh_del(aux_exists, config->aux_keep, iter);
            }
        }
        if (config->aux_remove) {
            //remove set in use, add RG if not present
            iter = kh_get(aux_exists, config->aux_remove, TAGNUM(rg));
            if (iter == kh_end(config->aux_remove)) {
                kh_put(aux_exists, config->aux_remove, TAGNUM(rg), &ret);
            }
        }
    }

    for (auxdata = bam_aux_first(bamdata); auxdata; ) {
        tag = bam_aux_tag(auxdata);
        if (config->aux_keep) {                         //keep option or remove option with ^ in use
            iter = kh_get(aux_exists, config->aux_keep, TAGNUM(tag));
            if (iter == kh_end(config->aux_keep)) {     //not present in keep, remove
                auxdata = bam_aux_remove(bamdata, auxdata);
            }
            else {                                      //present, keep
                auxdata = bam_aux_next(bamdata, auxdata);
            }
        }
        else if (config->aux_remove) {                  //remove option in use
            iter = kh_get(aux_exists, config->aux_remove, TAGNUM(tag));
            if (iter != kh_end(config->aux_remove)) {   //present in remove, remove
                auxdata = bam_aux_remove(bamdata, auxdata);
            }
            else {                                      //not present, keep
                auxdata = bam_aux_next(bamdata, auxdata);
            }
        }
        //else impossible
    }
}

/// getRGlines - add RG lines from input header to output header
/** @param in_samhdr - pointer to input sam header data
 *  @param out_samhdr - pointer to output sam header data
returns 1 on failure 0 on success
*/
int getRGlines(sam_hdr_t *in_samhdr, sam_hdr_t *out_samhdr)
{
    kstring_t line = KS_INITIALIZE;
    int i = 0, ret = 0, count = 0;
    const char rg[] = "RG";

    if (!in_samhdr || !out_samhdr) {
        fprintf(stderr, "Invalid parameters in getRGlines!\n");
        return 1;
    }

    if (-1 == (count = sam_hdr_count_lines(in_samhdr, rg))) {
        fprintf(stderr, "Failed to get RG count!\n");
        return 1;
    }

    for (i = 0; i < count; ++i)
    {
        ks_clear(&line);
        if (sam_hdr_find_line_pos(in_samhdr, rg, i, &line)) {
            fprintf(stderr, "Failed to get RG data!\n");
            ret = 1;
            break;
        }
        if (sam_hdr_add_lines(out_samhdr, line.s, line.l)) {
            fprintf(stderr, "Failed to add RG data!\n");
            ret = 1;
            break;
        }
    }
    ks_free(&line);

    return ret;
}

/// getPGlines - add PG lines from input header to output header based on user option
/** @param in_samhdr - pointer to input sam header data
 *  @param out_samhdr - pointer to output sam header data
 *  @param config - pointer to internal configuration data
 *  @param argdump - string containing dump of command line invocation
returns 1 on failure 0 on success
*/
int getPGlines(sam_hdr_t *in_samhdr, sam_hdr_t *out_samhdr, conf_data *config, const char *argdump)
{
    kstring_t line = KS_INITIALIZE, id = KS_INITIALIZE;
    int i = 0, ret = 0, count = 0;
    const char pg[] = "PG";

    if (!in_samhdr || !out_samhdr || !config) {
        fprintf(stderr, "Invalid parameters in getPGlines!\n");
        return 1;
    }

    if (-1 == (count = sam_hdr_count_lines(in_samhdr, pg))) {
        fprintf(stderr, "Failed to get PG count!\n");
        return 1;
    }

    if (config->pgid && config->pgid[0]) {            //when reject-PG is given, and is not empty, remove given pg onwards
        for (i = 0; i < count; ++i) {
            if (sam_hdr_find_tag_pos(in_samhdr, pg, i, "ID", &id)) {
                fprintf(stderr, "Failed to get PG entry fields for line %d!\n", i + 1);
                break;
            }

            if (!strcmp(id.s, config->pgid))
                break;

            //either current PG is prior to rejected one or all PGs are in, get PG line and add
            ks_clear(&line);
            if (sam_hdr_find_line_pos(in_samhdr, "PG", i, &line)) {
                fprintf(stderr, "Failed to get PG data at %d!\n", i + 1);
                ret = 1;
                break;
            }

            //add to output
            if (sam_hdr_add_lines(out_samhdr, line.s, line.l)) {
                fprintf(stderr, "Failed to add PG data!\n");
                ret = 1;
                break;
            }
        }
    }
    else {        //keep all
        for (i = 0; i < count; ++i) {
            if (sam_hdr_find_line_pos(in_samhdr, "PG", i, &line)) {
                fprintf(stderr, "Failed to get PG data at %d!\n", i + 1);
                ret = 1;
                break;
            }
            //line has the required PG data
            if (sam_hdr_add_lines(out_samhdr, line.s, line.l)) {
                fprintf(stderr, "Failed to add PG data!\n");
                ret = 1;
                break;
            }
        }
    }

    if (!ret && !config->noPGentry) {
        //add PG entry with reset command
        if (-1 == (ret = sam_hdr_add_pg(out_samhdr, "samtools", "CL", argdump, NULL))) {
            fprintf(stderr, "Failed to set PG entry!\n");
        }
    }
    ks_free(&line);
    ks_free(&id);

    return ret;
}

/// reset - do the reset of data and create output; create output header with required rg/pg data, add bamdata with flags set to unmapped, pair info and orientation reset,
// reerse and complement alignment if required
/** @param infile - input samfile pointer
 *  @param outfile - output sam file pointer
 *  @param config - pointer to internal configuration data
 *  @param args - string containing dump of command line invocation
returns 1 on failure 0 on success
*/
int reset(samFile *infile, samFile *outfile, conf_data *config, char *args)
{
    sam_hdr_t *in_samhdr = NULL, *out_samhdr = NULL;
    int ret = EXIT_FAILURE, ret_r = 0, ret_w = 0, i = 0;
    bam1_t *bamdata = NULL, *outdata = NULL;
    kstring_t querydata = KS_INITIALIZE, qualdata = KS_INITIALIZE;
    char *sp = NULL, *qp = NULL;
    uint8_t *bamquery = NULL, *bamqual = NULL;

    if (!infile || !outfile) {
        fprintf(stderr, "Invalid parameters in reset!\n");
        goto error;
    }

    //read input header
    in_samhdr = sam_hdr_read(infile);
    if (!in_samhdr)
    {
        fprintf(stderr, "Failed to read header from file!\n");
        goto error;
    }
    //create output header
    if (!(out_samhdr = sam_hdr_init()))
    {
        fprintf(stderr, "Failed to create output header!\n");
        goto error;
    }

    //add version to output header
    if  (-1 == sam_hdr_add_line(out_samhdr,"HD", "VN", SAM_FORMAT_VERSION, NULL)) {
        fprintf(stderr, "Failed to set header data!\n");
        goto error;
    }
    //add RG / PG lines if configured
    if ((config->keepRGs && getRGlines(in_samhdr, out_samhdr)) ||
            getPGlines(in_samhdr, out_samhdr, config, args)) {
        goto error;
    }

    //write output header
    if (sam_hdr_write(outfile, out_samhdr)) {
        print_error_errno("reset", "Output header write failed (%d)!\n", errno);
        goto error;
    }

    bamdata = bam_init1();      //input bam
    outdata = bam_init1();      //output bam
    if (!bamdata || !outdata)
    {
        fprintf(stderr, "Failed to allocate data memory!\n");
        goto error;
    }

    errno = 0; i = 0;
    sp = NULL; qp = NULL;
    bamquery = NULL; bamqual = NULL;

    //get bam data, make updates and dump to output
    while (0 <= (ret_r = sam_read1(infile, in_samhdr, bamdata)))
    {
        sp = NULL; qp = NULL;
        bamquery = NULL; bamqual = NULL;

        // read data
        if (bamdata->core.flag & BAM_FSECONDARY || bamdata->core.flag & BAM_FSUPPLEMENTARY) {
            continue;
        }

        //update flags
        uint16_t flags = bamdata->core.flag & ~BAM_FPROPER_PAIR;    //reset pair info
        flags |= BAM_FUNMAP;                                        //mark as unmapped
        if (bamdata->core.flag & BAM_FPAIRED) {
            flags |= BAM_FMUNMAP;                                   //mark mate as unmapped, if it was a pair
        }
        flags &= ~BAM_FMREVERSE;                                    //reset mate orientation

        if (0 > ks_resize(&querydata, bamdata->core.l_qseq) ||
            0 > ks_resize(&qualdata, bamdata->core.l_qseq)) {
            fprintf(stderr, "Failed to get allocate memory!\n");
            ret_r = -4;
            break;
        }
        ks_clear(&querydata);
        ks_clear(&qualdata);

        sp = ks_str(&querydata);
        qp = ks_str(&qualdata);
        bamquery = bam_get_seq(bamdata);
        bamqual = bam_get_qual(bamdata);
        if (bamdata->core.flag & BAM_FREVERSE) {
            //sequence data ordered as reverse complemented, reorder/complement sequence and quality data as read and clear the flag
            for (i = bamdata->core.l_qseq - 1; i >= 0; --i) {
                *sp++ = "=TGKCYSBAWRDMHVN"[bam_seqi(bamquery, i)];
                *qp++ = bamqual[i];
            }
            flags &= ~BAM_FREVERSE;                                 //reset flag as well
        }
        else {
            //data in read order itself
            for (i = 0; i < bamdata->core.l_qseq ; ++i) {
                *sp++ = seq_nt16_str[bam_seqi(bamquery, i)];
            }
            memcpy(qp, bam_get_qual(bamdata), bamdata->core.l_qseq);
        }

        removeauxtags(bamdata, config);
        if (0 > (ret_w = bam_set1(outdata, bamdata->core.l_qname - bamdata->core.l_extranul - 1, bam_get_qname(bamdata), flags, -1, -1, 0, 0, NULL, -1, -1, 0, bamdata->core.l_qseq, querydata.s, qualdata.s, bam_get_l_aux(bamdata)))) {
            print_error_errno("reset", "Failed to set output data (%d)!\n", errno);
            break;
        }

        memcpy(bam_get_aux(outdata), bam_get_aux(bamdata), bam_get_l_aux(bamdata));
        outdata->l_data += bam_get_l_aux(bamdata);

        errno = 0;
        //write bam data to output
        if (0 > (ret_w = sam_write1(outfile, out_samhdr, outdata)))
        {
            print_error_errno("reset", "Failed to write output data (%d)!\n", errno);
            break;
        }
        // wrote the data, continue read/write cycle
        errno = 0;
    }

    if (-1 > ret_r || 0 > ret_w) {
        //some error
        fprintf(stderr, "Error during %s!\n", (-1 > ret_r)? "read" : "write");
    }
    else {
        // no error!
        ret = EXIT_SUCCESS;
    }

error:
    // clean up and return result
    if (in_samhdr)
        sam_hdr_destroy(in_samhdr);
    if (out_samhdr)
        sam_hdr_destroy(out_samhdr);

    if (bamdata)
        bam_destroy1(bamdata);
    if (outdata)
        bam_destroy1(outdata);

    if (qualdata.s)
        ks_free(&qualdata);
    if (querydata.s)
        ks_free(&querydata);
    return ret;
}

/// cleanup - free up allocations made
/** @param config - pointer to internal configuration data
returns nothing
*/
void cleanup(conf_data *config)
{
    if (config->aux_keep) {
        kh_destroy(aux_exists, config->aux_keep);
        config->aux_keep = NULL;
    }
    if (config->aux_remove) {
        kh_destroy(aux_exists, config->aux_remove);
        config->aux_remove = NULL;
    }
}

/// main_reset - starts the reset of data
/** @param argc - count of arguments
 *  @param argv - pointer to array of arguments
returns 1 on failure 0 on success
*/
int main_reset(int argc, char *argv[])
{
    static const struct option lopts[] = {
        SAM_OPT_GLOBAL_OPTIONS('-', '-', 'O', '-', '-', '@'),       //let output format and thread count be given by user - long options
        {"keep-tag", required_argument, NULL, LONG_OPT('x')},       //aux tags to be retained, supports ^ STR
        {"remove-tag", required_argument, NULL, 'x'},               //aux tags to be removed
        {"no-RG", no_argument, NULL, 1},                            //no RG lines in output, default is to keep them
        //reject PG lines from input, default is to keep them (i.e. option not given); without optional filename, all PGs removed and those given in file are filtered when optional filename is given
        {"reject-PG", required_argument, NULL, 'p'},                //reject entries from this PG onwards
        {"no-PG", no_argument, NULL, 2},                            //do not add PG entry for reset operation, default is to add it
        {NULL, 0, NULL, 0}
    };
    samFile *infile = NULL, *outfile = NULL;
    sam_global_args ga = SAM_GLOBAL_ARGS_INIT;
    htsThreadPool tpool = {NULL, 0};
    const char *inname = NULL, *outname = NULL;
    int c = 0, ret = EXIT_FAILURE;
    char outmode[4] = "w", *args = NULL;
    conf_data resetconf = {1, 0, NULL, NULL, NULL};                //keep RGs and PGs by default


    //samtools reset -o outfile -x/--remove-tag ... --keep-tag ... --threads=n --output-fmt=fmt --no-RG --reject-PG pgid --no-PG [<infile>]
    while ((c = getopt_long(argc, argv, "o:@:x:O:", lopts, NULL)) >= 0)
    {
        switch (c)
        {
        case 1:                             //--no-RG
            if (!resetconf.keepRGs) {
                usage(stderr);              //already given!
                goto exit;
            }
            resetconf.keepRGs = 0;
            break;
        case 2:                             //--no-PG
            if (resetconf.noPGentry) {
                usage(stderr);              //already given!
                goto exit;
            }
            resetconf.noPGentry = 1;
            break;
        case 'p':                           //--reject-PG=<id>
            if (resetconf.pgid) {
                usage(stderr);              //already given!
                goto exit;
            }
            resetconf.pgid = optarg;
            break;
        case 'o':                           //output file name
            if (outname) {                  //already given!
                usage(stderr);
                goto exit;

            }
            outname = optarg;
            break;
        case 'x':                           //remove aux tag
            if (*optarg == '^') {           //remove all except given ones!
                if (parse_aux_list(&resetconf.aux_keep, optarg+1, "main_reset")) {
                    usage(stderr);
                    goto exit;
                }
            }
            else {                          //remove given ones
                if (parse_aux_list(&resetconf.aux_remove, optarg, "main_reset")) {
                    usage(stderr);
                    goto exit;
                }
            }
            break;
        case LONG_OPT('x'):                  //keep aux tags
            if (parse_aux_list(&resetconf.aux_keep, optarg, "main_reset")) {
                usage(stderr);
                goto exit;
            }
            break;
        // handle standard samtool options like thread count, verbosity...
        default:
            if (parse_sam_global_opt(c, optarg, lopts, &ga) == 0)
                break;
            // else fall-through
            // couldn't parse or unknown options, show usage!
        case '?':                           //unknown options found!
            usage(stderr);
            goto exit;
            break;
        }
    }

    if (argc == 1 && isatty(STDIN_FILENO)) {
        //no args and input is stdin -- it is the usage check
        usage(stdout);
        ret = EXIT_SUCCESS;
        goto exit;
    }
    //else have other args or input from redirection/pipe/other device -- validate and work

    if (!outname)
        outname = "-";

    //check and fail if unnecessary parameters are given
    c = argc - optind;
    if (c > 1) {
        usage(stderr);
        goto exit;
    }

    if (c == 1) {
        inname = argv[optind];
    }
    else {
        inname = "-";
    }

    //set output file format based on name
    sam_open_mode(outmode + 1, outname, NULL);

    //open input and output files
    infile = sam_open(inname, "r");
    outfile = sam_open_format(outname, outmode, &ga.out);
    if (!infile || !outfile) {
        fprintf(stderr, "Could not open %s%s%s\n", !infile ? inname : "", (!infile && !outfile)? ", " : "", !outfile ? outname : "");
        goto exit;
    }

    // set the thread count if given as argument
    if (ga.nthreads > 0)
    {
        if (!(tpool.pool = hts_tpool_init(ga.nthreads)))
        {
            fprintf(stderr, "\nFailed to setup thread pool\n");
            goto exit;
        }

        hts_set_opt(infile, HTS_OPT_THREAD_POOL, &tpool);
        hts_set_opt(outfile, HTS_OPT_THREAD_POOL, &tpool);
    }

    args = stringify_argv(argc + 1, argv - 1);              //to dump invocation in PG line

    //do the reset!
    ret = reset(infile, outfile, &resetconf, args);

exit:
    if (args)
        free(args);
    if (infile)
        sam_close(infile);
    if (outfile)
        sam_close(outfile);
    if (tpool.pool)
        hts_tpool_destroy(tpool.pool);
    cleanup(&resetconf);
    sam_global_args_free(&ga);

    return ret;
}
