/*  bam_fastaref.c -- fastaref subcommand.

    Copyright (C) 2017 Genome Research Ltd.

    Author: Thomas Hickman <th10@sanger.ac.uk>

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

#include "samtools.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>
#include <unistd.h>

#include "htslib/sam.h"
#include "htslib/ref.h"
#include "htslib/kstring.h"
#include "htslib/khash.h"
#include "htslib/hfile.h"
#include "htslib/cram.h"
#include "sam_header.h"

KHASH_INIT(s2s, char*, char*, 1, kh_str_hash_func, kh_str_hash_equal)

typedef struct{
    int maxLineLength;
    char **outputKeys;
    int numOutputKeys;
    char* samFileName;
    char* outFileName; // NULL if this is to output to stdout
} FastarefOptions;

void destroy_malloced_s2s_khash(khash_t(s2s)* khash){
    char* value;
    char* key;

    kh_foreach(khash, value, key, {
        free(value);
        free(key);
    });

    kh_destroy(s2s, khash);
}

/**
 * Given a string, pointing to the start of a SQ line, parses
 * it into a hash map, which is populated into out_khash
 * 
 * @returns The length of the string parsed if successful
 *          -1 if not successful
 */
int parseSQLine(const char* sqLine, khash_t(s2s)* out_khash){
    const char* fieldStart = sqLine + 4;// 4 == len("@SQ\t")
    const char* fieldEnd;

    while(1){
        fieldEnd = strpbrk(fieldStart, "\t\n\0");

        int valueLength = fieldEnd - fieldStart - 3;
        if(valueLength < 1){
            print_error("fastaref", "failed parsing SQ line in header");
            return -1;
        }

        char* key = malloc(3);
        memcpy(key, fieldStart, 2);
        key[2] = '\0';

        char* value = malloc(valueLength + 1);
        memcpy(value, fieldStart + 3, valueLength);
        value[valueLength] = '\0';

        int ret;
        int keyPointer = kh_put(s2s, out_khash, key, &ret);
        if(ret == -1){
            print_error("fastaref", "failed parsing SQ line in header");
            free(value);
            free(key);
            return -1;
        }
        else if(ret == 0){
            print_error("fastaref", "duplicate field name in SQ line");
            free(value);
            free(key);
            return -1;
        }
        kh_value(out_khash, keyPointer) = value;

        if(*fieldEnd != '\t'){
            // reached the end of the file or line
            break;
        }
        fieldStart = fieldEnd + 1; // move fieldStart to the end of the tab
    }

    return fieldEnd - fieldStart;
}

/**
 * Generates a fasta file according to command line options specified in options
 * 
 * @returns 0 if successful
 */
int generateFastaFile(FastarefOptions *options){
    const char* readPosition = NULL;
    int writeBufferLength = options->maxLineLength;
    char* writeBuffer = malloc(writeBufferLength);
    FILE* outFile = NULL;
    int returnCode = 0;
    bam_hdr_t* header = NULL;
    hFILE* ref = NULL;
    khash_t(s2s)* khash = kh_init(s2s);

    samFile* inFile = sam_open(options->samFileName, "r");
    if (inFile == NULL) {
        print_error_errno("fastaref", "failed to open \"%s\"", options->samFileName);
        returnCode = 1;
        goto cleanup;
    }

    header = sam_hdr_read(inFile);
    if (header == NULL) {
        print_error("fastaref", "failed to read header for \"%s\"", options->samFileName);
        returnCode = 1;
        goto cleanup;
    }

    readPosition = header->text;

    outFile = options->outFileName?fopen(options->outFileName, "w") : stdout;

    if(outFile == NULL){
        print_error_errno("fastaref", "failed to write file \"%s\"", options->outFileName);
        returnCode = 1;
        goto cleanup;
    }

    int atStartOfFile = 1;

    while(1){
        for(; !(*readPosition == '\0'
            || (atStartOfFile && strncmp(readPosition, "@SQ", 3) == 0)
            || (!atStartOfFile && strncmp(readPosition, "@SQ", 3) == 0 && *(readPosition - 1) == '\n')
        );readPosition++){
            if(atStartOfFile) atStartOfFile = 0;
        }
        if(*readPosition == '\0') break;

        int parsedLength;
        if((parsedLength = parseSQLine(readPosition, khash)) == -1){
            goto cleanup;
        }
        readPosition += parsedLength;

        khint_t M5_pointer = kh_get(s2s, khash, "M5");
        khint_t SN_pointer = kh_get(s2s, khash, "SN");

        if(SN_pointer == kh_end(khash)){
            print_error("fastaref", "error: SN field not found in SQ line in \"%s\"", options->samFileName);
            goto cleanup;
        }
        
        const char* SN = kh_value(khash, SN_pointer);
        if(M5_pointer == kh_end(khash)){
            print_error("fastaref", "warning: no M5 string found for sequence \"%s\" in \"%s\"", SN, options->samFileName);
        }
        else{
            const char* M5 = kh_value(khash, M5_pointer);

            if(fputc('>', outFile) < 0 
            || fputs(SN, outFile) < 0){
                print_error("fastaref", "failed to write to output file");
                returnCode = 1;
                goto cleanup;
            }
            
            for(int i = 0;i < options->numOutputKeys;i++){
                char* fieldName = options->outputKeys[i];
                khint_t fieldPointer = kh_get(s2s, khash, fieldName);
                if(fieldPointer != kh_end(khash)){
                    if(fputc('\t', outFile) < 0
                    || fputs(fieldName, outFile) < 0
                    || fputc(':', outFile) < 0
                    || fputs(kh_value(khash, fieldPointer), outFile) < 0){
                        print_error("fastaref", "failed to write to output file");
                        returnCode = 1;
                        goto cleanup;
                    }
                }
            }
            if(fputc('\n', outFile) < 0){
                print_error("fastaref", "failed to write to output file");
                returnCode = 1;
                goto cleanup;
            }

            if(m5_to_ref(M5, &ref) != 0){
                print_error("fastaref", "could not fetch the reference with md5 of %s", M5);
                returnCode = 1;
                goto cleanup;
            }

            while(1){
                int lengthRead = hread(ref, writeBuffer, writeBufferLength);

                if(lengthRead < 0){
                    print_error_errno("fastaref", "failed to read the reference with md5 of \"%s\"", M5);
                    returnCode = 1;
                    goto cleanup;
                }

                if(fwrite(writeBuffer, sizeof(char), lengthRead, outFile) != lengthRead
                || fputc('\n', outFile) < 0){
                    print_error_errno("fastaref", "failed to write file \"%s\"", options->outFileName);
                    returnCode = 1;
                    goto cleanup;
                }

                if(lengthRead < writeBufferLength) break;
            }

            destroy_malloced_s2s_khash(khash); khash = kh_init(s2s);;
            (void)hclose(ref); ref = NULL;
        }
    }

cleanup:
    destroy_malloced_s2s_khash(khash);
    if(ref) (void)hclose(ref);
    if(writeBuffer) free(writeBuffer);
    if(outFile) fclose(outFile);
    if(header) bam_hdr_destroy(header);
    if(inFile) sam_close(inFile);
    return returnCode;
}

/**
 * Parses keysStr to obtain a list of two letter keys
 * 
 * @returns 0 if successful
 */
char** parse_SQ_keys_list(const char* keysStr, int* numKeys){
    int keysStrLength = strlen(keysStr);

    if((keysStrLength + 1) / 3 * 3 != keysStrLength + 1){ // keysStrLength + 1 is not divisible by 3
        fprintf(stderr, "%s is not a valid set of keys\n", keysStr);
        return NULL;
    }

    *numKeys = (keysStrLength + 1) / 3;
    char** outputKeys = malloc(*numKeys * sizeof(char*));

    for(int i = 0;i < *numKeys;i++){
        outputKeys[i] = malloc(3);
        memcpy(outputKeys[i], keysStr + i*3, 2);
        outputKeys[i][2] = 0;

        if (*numKeys - 1 != i && keysStr[i*3 + 2] != ','){
            fprintf(stderr, "%s is not a valid set of keys\n", keysStr);
            return NULL;
        }
    }
    
    return outputKeys;
}

void print_usage(){
    fprintf(stderr, "Usage: samtools fastaref [options] <in.bam>\n\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  -o FILE  Output file name [stdout]\n");
    fprintf(stderr, "  -k LIST  Output the specified keys into the fasta file (separated by commas)  [LN,AH,AN,AS,M5,SP]\n");
    fprintf(stderr, "  -l INT   Maximum length of outputted lines [60]\n");
}

int main_fastaref(int argc, char *argv[])
{
    FastarefOptions options;
    options.numOutputKeys = -1;
    options.maxLineLength = 60;
    options.outFileName = NULL; // (point to stdout)

    int opt;
    while((opt = getopt(argc, argv, "k:l:o:")) >= 0){
        switch (opt) {
            case 'k':
                if((options.outputKeys = parse_SQ_keys_list(optarg, &(options.numOutputKeys))) == NULL){
                    return 1;
                }
            break;
            case 'l':
                options.maxLineLength = atoi(optarg);
                if(options.maxLineLength <= 0){
                    fprintf(stderr, "invalid maximum line length\n");
                    return 1;
                }
            break;
            case 'o':
                options.outFileName = optarg;
            break;
        }
    }

    if(options.numOutputKeys == -1){ // i.e. -k is not set
        if((options.outputKeys = parse_SQ_keys_list("LN,AH,AN,AS,M5,SP", &(options.numOutputKeys))) == NULL){
            return 1;
        }
    }

    int numOtherArgs = argc - optind;
    char** otherArgs = argv + optind;

    if(numOtherArgs == 0 && !isatty(STDIN_FILENO)){
        options.samFileName = "-";
    }
    else if(numOtherArgs == 1){
        options.samFileName = otherArgs[0];
    }
    else{
        print_usage();
        return 1;
    }

    int returnValue = generateFastaFile(&options);
    for(int i = 0;i < options.numOutputKeys;i++){
        free(options.outputKeys[i]);
    }
    free(options.outputKeys);

    return returnValue;
}

