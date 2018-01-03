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

#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "htslib/cram.h"
#include "htslib/hfile.h"
#include "htslib/khash.h"
#include "htslib/kstring.h"
#include "htslib/ref.h"
#include "htslib/sam.h"
#include "sam_header.h"

KHASH_MAP_INIT_STR(s2s, char*)

typedef struct {
    int maxLineLength;
    char (*outputKeys)[3];
    int numOutputKeys;
    char* samFileName;
    char* outFileName; // NULL if this is to output to stdout
} fastaref_options_t;

/**
 * Replacement in C99 for strsep
 */
static char* _strsep(char** stringp, const char* delim) {
    if (*stringp == NULL) return NULL;
    char* foundChar = strpbrk(*stringp, delim);

    if (foundChar == NULL) {
        char* token = *stringp;
        *stringp = NULL;
        return token;
    }

    *foundChar = 0;

    char* token = *stringp;
    *stringp = foundChar + 1;
    return token;
}

/**
 * Given an SQ line, returns a hash map, containing it's fields.
 * Note: this modifies the text in the sqLine parameter
 */
kh_s2s_t* parseSQLine(char* sqLine) {
    char* parsePosition = sqLine + 4; // 4 = len("@SQ\t")
    char* field;
    kh_s2s_t* sqLineDict = kh_init(s2s);

    while ((field = _strsep(&parsePosition, "\t")) != NULL) {
        if (strlen(field) < 4) {
            // cannot have a field of length < 4, we need to have "FD:" present
            // where FD is the field name
            print_error("fastaref", "failed parsing SQ line in header");
            return NULL;
        }

        char* key = field;
        key[2] = '\0';

        char* value = field + 3;

        int ret;
        int valuePointer = kh_put(s2s, sqLineDict, key, &ret);
        if (ret == -1) {
            print_error("fastaref", "failed parsing SQ line in header");
            return NULL;
        }
        else if (ret == 0) {
            print_error("fastaref", "duplicate field name in SQ line");
            return NULL;
        }
        kh_value(sqLineDict, valuePointer) = value;
    }

    return sqLineDict;
}

/**
 * Generates a fasta file according to command line options specified in options
 *
 * @returns 0 if successful
 */
int generateFastaFile(fastaref_options_t* options) {
    int writeBufferLength = options->maxLineLength;
    char* writeBuffer = malloc(writeBufferLength);
    if (writeBuffer == NULL) return -1;
    FILE* outFile = NULL;
    int returnCode = 0;
    bam_hdr_t* header = NULL;
    hFILE* ref = NULL;
    kh_s2s_t* sqLineDict = NULL;
    char* headerTextCopy = NULL;

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

    headerTextCopy = strdup(header->text);
    char* readPosition = headerTextCopy;

    outFile = options->outFileName ? fopen(options->outFileName, "w") : stdout;
    if (outFile == NULL) {
        print_error_errno("fastaref", "failed to write file \"%s\"", options->outFileName);
        returnCode = 1;
        goto cleanup;
    }

    char* line;
    while ((line = _strsep(&readPosition, "\n")) != NULL) {
        if (strncmp(line, "@SQ", 3) == 0) {
            if ((sqLineDict = parseSQLine(line)) == NULL) {
                returnCode = 1;
                goto cleanup;
            }

            khint_t M5_pointer = kh_get(s2s, sqLineDict, "M5");
            khint_t SN_pointer = kh_get(s2s, sqLineDict, "SN");

            if (SN_pointer == kh_end(sqLineDict)) {
                print_error("fastaref", "error: SN field not found in SQ line in \"%s\"",
                            options->samFileName);
                goto cleanup;
            }

            const char* SN = kh_value(sqLineDict, SN_pointer);
            if (M5_pointer == kh_end(sqLineDict)) {
                print_error("fastaref", "warning: no M5 string found for sequence \"%s\" in \"%s\"", SN,
                            options->samFileName);
            }
            else {
                const char* M5 = kh_value(sqLineDict, M5_pointer);

                if (fputc('>', outFile) < 0 || fputs(SN, outFile) < 0) {
                    print_error("fastaref", "failed to write to output file");
                    returnCode = 1;
                    goto cleanup;
                }

                int i;
                for (i = 0; i < options->numOutputKeys; i++) {
                    char* fieldName = options->outputKeys[i];
                    khint_t fieldPointer = kh_get(s2s, sqLineDict, fieldName);
                    if (fieldPointer != kh_end(sqLineDict)) {
                        if (fputc('\t', outFile) < 0 || fputs(fieldName, outFile) < 0 || fputc(':', outFile) < 0 || fputs(kh_value(sqLineDict, fieldPointer), outFile) < 0) {
                            print_error("fastaref", "failed to write to output "
                                                    "file");
                            returnCode = 1;
                            goto cleanup;
                        }
                    }
                }
                if (fputc('\n', outFile) < 0) {
                    print_error("fastaref", "failed to write to output file");
                    returnCode = 1;
                    goto cleanup;
                }

                if (!(ref = m5_to_ref(M5))) {
                    print_error("fastaref", "could not fetch the reference with MD5 of %s", M5);
                    returnCode = 1;
                    goto cleanup;
                }

                while (1) {
                    int lengthRead = hread(ref, writeBuffer, writeBufferLength);

                    if (lengthRead < 0) {
                        print_error_errno("fastaref", "failed to read the reference with MD5 of \"%s\"", M5);
                        returnCode = 1;
                        goto cleanup;
                    }

                    if (fwrite(writeBuffer, sizeof(char), lengthRead, outFile) != lengthRead || fputc('\n', outFile) < 0) {
                        print_error_errno("fastaref", "failed to write file \"%s\"",
                                          options->outFileName);
                        returnCode = 1;
                        goto cleanup;
                    }

                    if (lengthRead < writeBufferLength) break;
                }
                kh_destroy(s2s, sqLineDict);
                sqLineDict = NULL;
                (void)hclose(ref);
                ref = NULL;
            }
        }
    }

cleanup:
    if (headerTextCopy) free(headerTextCopy);
    if (sqLineDict) kh_destroy(s2s, sqLineDict);
    if (ref) (void)hclose(ref);
    if (writeBuffer) free(writeBuffer);
    if (outFile) fclose(outFile);
    if (header) bam_hdr_destroy(header);
    if (inFile) sam_close(inFile);
    return returnCode;
}

void print_usage() {
    fprintf(stderr, "Usage: samtools fastaref [options] <in.bam>\n\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  -o FILE  Output file name [stdout]\n");
    fprintf(stderr, "  -k LIST  Output the specified keys into the fasta file "
                    "(separated by commas)  [LN,AH,AN,AS,M5,SP]\n");
    fprintf(stderr, "  -l INT   Maximum length of outputted lines [60]\n");
}

int main_fastaref(int argc, char* argv[]) {
    fastaref_options_t options;
    options.numOutputKeys = 0;
    options.maxLineLength = 60;
    options.outFileName = NULL; // (point to stdout)

    int opt;
    while ((opt = getopt(argc, argv, "k:l:o:")) >= 0) {
        switch (opt) {
        case 'k': {
            const char* key;
            const char* strStart = optarg;
            while ((key = _strsep(&optarg, ",")) != NULL) {
                if (strlen(key) != 2) {
                    print_error("fastaref", "invalid filtering key '%s', keys must have length 2", key);
                }
                options.numOutputKeys++;
            }
            options.outputKeys = (char(*)[3])strStart;
        } break;
        case 'l':
            options.maxLineLength = atoi(optarg);
            if (options.maxLineLength <= 0) {
                print_error("fastaref", "invalid maximum line length");
                return 1;
            }
            break;
        case 'o':
            options.outFileName = optarg;
            break;
        }
    }

    if (options.numOutputKeys == 0) { // i.e. -k is not set
        char defaultValues[6][3] = {"LN", "AH", "AN", "AS", "M5", "SP"};

        options.outputKeys = defaultValues;
        options.numOutputKeys = 6;
    }

    int numOtherArgs = argc - optind;
    char** otherArgs = argv + optind;

    if (numOtherArgs == 0 && !isatty(STDIN_FILENO)) {
        options.samFileName = "-";
    }
    else if (numOtherArgs == 1) {
        options.samFileName = otherArgs[0];
    }
    else {
        print_usage();
        return 1;
    }

    return generateFastaFile(&options);
}
