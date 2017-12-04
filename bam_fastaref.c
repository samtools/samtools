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

#include <htslib/sam.h>
#include <htslib/ref.h>
#include <htslib/kstring.h>
#include <htslib/hfile.h>

/**
 * Given a pointer to the start of a sam field, populates samField with a
 * the malloced sam field.
 * 
 * Updates samFieldPosition to be the position after the sam field.
 * 
 * @returns 1 if successful
*/
int get_sam_field(char** samFieldPosition, char** samField){
    char* samFieldStart = *samFieldPosition;
    for(;**samFieldPosition != '\t' && **samFieldPosition != '\n' && **samFieldPosition != 0;(*samFieldPosition)++){};
    int samFieldLength = *samFieldPosition - samFieldStart;

    if((*samField = malloc(samFieldLength + 1)) == NULL) return 1;

    strncpy(*samField, samFieldStart, samFieldLength);
    (*samField)[samFieldLength] = '\0';

    return 0;
}

typedef struct{
    int maxLineLength;
    char outputKeys[8][2];
    int numOutputKeys;
    char* samFileName;
    char* outFileName; // NULL if this is to output to stdout
} FastarefOptions;

void init_kstring(kstring_t* str){
    str->l = str->m = 0; str->s = NULL;
}

/**
 * Given text of an SQ line beginning at sequenceText, populates values for the SN and M5
 * fields and returns selected portions of the header according to options
 * 
 * This updates sequenceText after reading it.
 */
int getInfoForSQLine(char** sequenceText, char** SN, char** M5, kstring_t* restHeader, FastarefOptions *options){
    init_kstring(restHeader);
    char* tmpReadPos;

    while(1){
        if(strncmp(*sequenceText, "M5:", 3) == 0){
            tmpReadPos = *sequenceText + 3;
            get_sam_field(&tmpReadPos, M5);
        }
        else if(strncmp(*sequenceText, "SN:", 3) == 0){
            tmpReadPos = *sequenceText + 3;
            get_sam_field(&tmpReadPos, SN);
        }

        int shouldCaptureKey = 0;

        for(int i = 0;i < options->numOutputKeys;i++){
            if(strncmp(*sequenceText, options->outputKeys[i], 2) == 0){
                shouldCaptureKey = 1;
                break;
            }
        }

        char* samField;
        
        char* fieldStart = *sequenceText;
        *sequenceText += 3;
        get_sam_field(sequenceText, &samField);
        
        if(shouldCaptureKey){
            kputc('\t', restHeader);
            kputsn(fieldStart, 3, restHeader);
            kputs(samField, restHeader);
        }

        free(samField);

        if(**sequenceText == '\n' || **sequenceText == '\0'){
            break;
        }
        else if(**sequenceText != '\t'){
            print_error("fastaref", "invalid header in %s unexpected '%s'", options->samFileName, *sequenceText);
            return 1;
        }

        (*sequenceText)++;
    }

    if(*SN == NULL){
        print_error("fastaref", "header in %s is missing a SN field in the SQ tag", options->samFileName);
        return 1;
    }

    return 0;
}
/**
 * Generates a fasta file according to command line options specified in options
 */
int generateFastaFile(FastarefOptions *options){
    char* readPosition = NULL;
    // the fasta header after the name
    kstring_t fqHeaderLineRest = {0, 0, NULL};
    char* SN = NULL;
    char* M5 = NULL;
    int writeBufferLength = options->maxLineLength;
    char* writeBuffer = malloc(writeBufferLength);
    FILE* outFile = NULL;
    int returnCode = 0;
    bam_hdr_t* header = NULL;
    hFILE* ref = NULL;

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

    kstring_t fqHeader = {0, 0, NULL};

    outFile = options->outFileName?fopen(options->outFileName, "w") : stdout;

    if(outFile == NULL){
        print_error_errno("fastaref", "failed to write file \"%s\"", options->outFileName);
        returnCode = 1;
        goto cleanup;
    }

    int atStartOfFile = 1;

    while(1){
        SN = NULL;
        M5 = NULL;

        for(; !(*readPosition == '\0'
            || (atStartOfFile && strncmp(readPosition, "@SQ", 3) == 0)
            || (!atStartOfFile && strncmp(readPosition, "@SQ", 3) == 0 && *(readPosition - 1) == '\n')
        );readPosition++){
            if(atStartOfFile) atStartOfFile = 0;
        }
        if(*readPosition == '\0') break;

        readPosition += 4;
        if(getInfoForSQLine(&readPosition, &SN, &M5, &fqHeaderLineRest, options) != 0){
            goto cleanup;
        }

        if(M5 == NULL){
            print_error("fastaref", "warning: no M5 string found for sequence \"%s\" in \"%s\"", SN, options->samFileName);
        }
        else{
            init_kstring(&fqHeader);
            kputc('>', &fqHeader);
            kputs(SN, &fqHeader);
            kputc('\t', &fqHeader);
            kputsn(ks_str(&fqHeaderLineRest), ks_len(&fqHeaderLineRest), &fqHeader);
            free(ks_str(&fqHeaderLineRest));

            if(m5_to_ref(M5, &ref) != 0){
                print_error("fastaref", "could not fetch the reference with md5 of %s", M5);
                returnCode = 1;
                goto cleanup;
            }

            if(fputs(ks_str(&fqHeader), outFile) < 0
            || fputc('\n', outFile) < 0){
                print_error_errno("fastaref", "failed to write file fasta header to \"%s\"", options->outFileName);
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

            free(ks_release(&fqHeader));
            free(M5); M5 = NULL;
            free(SN); SN = NULL;
            (void)hclose(ref); ref = NULL;
        }
    }

cleanup:
    if(ks_str(&fqHeader)) free(ks_str(&fqHeader));
    if(M5) free(M5);
    if(SN) free(SN);
    if(ref) (void)hclose(ref);
    if(writeBuffer) free(writeBuffer);
    if(outFile) fclose(outFile);
    if(header) bam_hdr_destroy(header);
    if(inFile) sam_close(inFile);
    return returnCode;
}

/**
 * Parses keysStr to obtain a list of two letter keys
 */
int parse_SQ_keys_list(char* keysStr, char outputKeys[][2], int* numKeys){
    int keysStrLength = strlen(keysStr);

    if((keysStrLength + 1) / 3 * 3 != keysStrLength + 1){ // keysStrLength + 1 is not divisible by 3
        fprintf(stderr, "%s is not a valid set of keys\n", keysStr);
        return 1;
    }

    if(*keysStr == 0){
        *numKeys = 0;
        return 0;
    }

    for(int i = 0;i < 8;i++){
        memcpy(outputKeys[i], keysStr + i*3, 2);

        if(keysStr[i*3 + 2] == '\0'){
            *numKeys = i + 1;
            return 0;
        }
        else if (keysStr[i*3 + 2] != ','){
            fprintf(stderr, "%s is not a valid set of keys\n", keysStr);
            return 1;
        }
    }

    fprintf(stderr, "cannot specify more than 8 keys (the SQ line should have a maximum of 8 fields)\n");
    return 1;
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
                if(parse_SQ_keys_list(optarg, options.outputKeys, &(options.numOutputKeys)) != 0){
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

    if(options.numOutputKeys == -1){
        if(parse_SQ_keys_list("LN,AH,AN,AS,M5,SP", options.outputKeys, &(options.numOutputKeys)) != 0){
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

    return generateFastaFile(&options);
}

