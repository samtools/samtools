/* bam_addrprg.c -- samtools command to add or replace readgroups.
 
 Copyright (c) 2013, 2015 Genome Research Limited.
 
 Author: Martin O. Pollard <mp15@sanger.ac.uk>
 
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

#include <htslib/sam.h>
#include <string.h>
#include <strings.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <limits.h>
#include <assert.h>
#include <unistd.h>

struct parsed_opts {
    char* input_name;
    char* output_name;
    char* output_mode;
};

struct state {
    samFile* input_file;
    bam_hdr_t* input_header;
    samFile* output_file;
    bam_hdr_t* output_header;
};

typedef struct parsed_opts parsed_opts_t;
typedef struct state state_t;

static void usage(FILE *fp)
{
    fprintf(fp,
            "Usage: samtools addreplacerg [options] <input.bam> <output.bam>\n"
            "Options:\n"
            "  -l INT     Set compression level, from 0 (uncompressed) to 9 (best)\n"
            "  -O FORMAT  Write output as FORMAT ('sam'/'bam'/'cram')\n");
}


static bool parse_args(int argc, char** argv, parsed_opts_t** opts)
{
    int n;
    char modeout[12], *fmtout;
    strcpy(modeout, "w");
    int level = -1;
    while ((n = getopt(argc, argv, "O:l:")) >= 0) {
        switch (n) {
            case 'O': fmtout = optarg; break;
            case 'l': level = atoi(optarg); break;
            case '?':
                usage(stdout);
                return true;
            default:
                usage(stderr);
                return false;
        }
    }

    if (argc+optind < 2) {
        usage(stdout);
        return true;
    }

    parsed_opts_t* retval = malloc(sizeof(parsed_opts_t));
    if (! retval ) {
        fprintf(stderr, "Out of memory\n");
        return false;
    }
    retval->input_name = strdup(argv[optind+0]);
    retval->output_name = strdup(argv[optind+1]);

    if (sam_open_mode(&modeout[1], retval->output_name, fmtout) < 0) {
        if (fmtout) fprintf(stderr, "[bam_addrprg] can't parse output format \"%s\"\n", fmtout);
        else fprintf(stderr, "[bam_addrprg] can't determine output format\n");
        return false;
    }

    if (level >= 0) sprintf(strchr(modeout, '\0'), "%d", level < 9? level : 9);

    *opts = retval;
    return true;
}

static bool init(parsed_opts_t* opts, state_t** state_out) {
    state_t* retval = (state_t*) malloc(sizeof(state_t));
    if (retval == NULL) {
        fprintf(stderr, "Out of memory\n");
        return false;
    }
    *state_out = retval;

    // Open files
    retval->input_file = sam_open(opts->input_name, "r");
    if (retval->input_file == NULL) {
        fprintf(stderr, "Could not open input file: %s\n", opts->input_name);
        return false;
    }
    retval->input_header = sam_hdr_read(retval->input_file);

    retval->output_header = bam_hdr_dup(retval->input_header);
    retval->output_file = sam_open(opts->output_name, opts->output_mode);
    
    if (retval->output_file == NULL) {
        fprintf(stderr, "Could not open output file: %s\n", opts->output_name);
        return false;
    }

    return true;
}

static char* getRGLine(const char* text)
{
    char* rg = strstr(text,"\n@RG");
    if (rg == NULL) {
        return NULL;
    }
    rg++;//skip initial \n
    char* end = strchr(rg, '\n');
    char* line;
    if (end) {
        line = strndup(rg,(end-rg));
    } else {
        line = strdup(rg);
    }

    return line;
}

static char* getRGID(const char* text)
{
    char *line, *next;
    line = getRGLine(text);

    assert(line!=NULL);

    next = line;
    char* token = strsep(&next, "\t");
    token = strsep(&next,"\t"); // skip first token it should always be "@RG"
    while (next != NULL) {
        char* key = strsep(&token,":");
        if (!strcmp(key,"ID")) {
            char* retval = strdup(token);
            free(line);
            return retval;
        }
        token = strsep(&next,"\t");
    }
    free(line);
    return NULL;
}

static bool readgroupise(state_t* state)
{
    char* id = getRGID(state->output_header->text);

    if ( !id ) {
        fprintf(stderr, "No RG specified in header\n");
        return false;
    }

    if (sam_hdr_write(state->output_file, state->output_header) != 0) {
        free(id);
        return false;
    }

    bam1_t* file_read = bam_init1();
    if (sam_read1(state->input_file, state->input_header, file_read) < 0) {
        bam_destroy1(file_read);
        file_read = NULL;
    }
    while (file_read != NULL) {
        uint8_t* data = (uint8_t*)strdup(id);
        int len = strlen(id)+1;
        // If the old exists delete it
        uint8_t* old = bam_aux_get(file_read, "RG");
        if (old != NULL) {
            bam_aux_del(file_read, old);
        }
        bam_aux_append(file_read, "RG",'Z',len,data);
        if (sam_write1(state->output_file, state->output_header, file_read) < 0) {
            fprintf(stderr, "Could not write read to output file.\n");
            return false;
        }
        if (sam_read1(state->input_file, state->input_header, file_read) < 0) {
            bam_destroy1(file_read);
            file_read = NULL;
        }
    }
    free(id);

    return true;
}

static void cleanup_opts(parsed_opts_t* opts)
{
    if (!opts) return;
    free(opts->output_name);
    free(opts->input_name);
}

static void cleanup_state(state_t* state)
{
    if (!state) return;
    if (state->output_file) sam_close(state->output_file);
    bam_hdr_destroy(state->output_header);
    if (state->input_file) sam_close(state->input_file);
    bam_hdr_destroy(state->input_header);
}

int main_addreplacerg(int argc, char** argv)
{
    parsed_opts_t* opts = NULL;
    state_t* state = NULL;

    if (!parse_args(argc, argv, &opts)) goto error;
    if (!opts) return 0;
    if (!opts || !init(opts, &state)) goto error;
    
    if (!readgroupise(state)) goto error;
    
    cleanup_opts(opts);
    cleanup_state(state);
    
    return 0;
error:
    cleanup_opts(opts);
    cleanup_state(state);

    return 1;
}
