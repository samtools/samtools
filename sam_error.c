/*  sam_error.c -- centralized error handling utilities implementation.

    Copyright (C) 2025 Akshay Dedaniya <Dedaniya08@hotmail.com>

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

#include "sam_error.h"
#include <stdarg.h>
#include <time.h>
#include <stdbool.h>

// =============================================================================
// Global State
// =============================================================================

static sam_log_level_t current_log_level = SAM_LOG_WARNING;
static bool verbose_logging = false;

// =============================================================================
// Logging Configuration Functions
// =============================================================================

void sam_set_log_level(sam_log_level_t level) {
    current_log_level = level;
}

void sam_set_verbose_logging(bool enabled) {
    verbose_logging = enabled;
}

// =============================================================================
// Core Logging Implementation
// =============================================================================

void sam_log_message(sam_log_level_t level, const char *cmd, 
                     const char *file, int line, const char *fmt, ...) {
    // Skip messages below current log level
    if (level > current_log_level) {
        return;
    }

    // Determine prefix and output stream based on severity
    const char *prefix;
    FILE *output = stderr;
    
    switch (level) {
        case SAM_LOG_ERROR:
            prefix = "ERROR";
            break;
        case SAM_LOG_WARNING:
            prefix = "WARNING";
            break;
        case SAM_LOG_INFO:
            prefix = "INFO";
            output = stdout;
            break;
        case SAM_LOG_DEBUG:
            prefix = "DEBUG";
            break;
        default:
            prefix = "UNKNOWN";
            break;
    }

    // Print timestamp for debug messages
    if (level == SAM_LOG_DEBUG && verbose_logging) {
        time_t now = time(NULL);
        struct tm *tm_info = localtime(&now);
        fprintf(output, "[%02d:%02d:%02d] ", 
                tm_info->tm_hour, tm_info->tm_min, tm_info->tm_sec);
    }

    // Print command and severity
    fprintf(output, "[%s] %s: ", cmd ? cmd : "samtools", prefix);

    // Print the actual message
    va_list args;
    va_start(args, fmt);
    vfprintf(output, fmt, args);
    va_end(args);

    // Add file/line information for verbose debug output
    if (verbose_logging && (level == SAM_LOG_ERROR || level == SAM_LOG_DEBUG)) {
        const char *filename = strrchr(file, '/');
        filename = filename ? filename + 1 : file;
        fprintf(output, " (%s:%d)", filename, line);
    }

    fprintf(output, "\n");
    fflush(output);
}
