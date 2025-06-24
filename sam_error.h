/*  sam_error.h -- centralized error handling utilities for samtools.

    Copyright (C) 2025 Akshay Dedaniya <Dedaniya08@hotmail.com>

    This module provides centralized error handling, logging, and diagnostic
    utilities to reduce code duplication and improve consistency across samtools.
    
    Features:
    - Standardized error message formatting
    - Memory allocation failure handling
    - File I/O error reporting
    - Debug and verbose logging support

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

#ifndef SAM_ERROR_H
#define SAM_ERROR_H

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif

// =============================================================================
// Error Severity Levels
// =============================================================================

typedef enum {
    SAM_LOG_ERROR,      // Fatal errors that terminate operation
    SAM_LOG_WARNING,    // Non-fatal issues that may affect results
    SAM_LOG_INFO,       // General information messages
    SAM_LOG_DEBUG       // Detailed debugging information
} sam_log_level_t;

// =============================================================================
// Error Handling Macros
// =============================================================================

/**
 * Report a warning and continue execution.
 * Use for non-fatal issues that may affect results but don't prevent operation.
 */
#define sam_warn(cmd, fmt, ...) \
    sam_log_message(SAM_LOG_WARNING, cmd, __FILE__, __LINE__, fmt, ##__VA_ARGS__)

/**
 * Report an error and continue execution.
 * Use for recoverable errors where operation can continue with fallback.
 */
#define sam_error(cmd, fmt, ...) \
    sam_log_message(SAM_LOG_ERROR, cmd, __FILE__, __LINE__, fmt, ##__VA_ARGS__)

/**
 * Report an informational message.
 * Use for general status updates and user feedback.
 */
#define sam_info(cmd, fmt, ...) \
    sam_log_message(SAM_LOG_INFO, cmd, __FILE__, __LINE__, fmt, ##__VA_ARGS__)

/**
 * Report a debug message (only if debug mode enabled).
 * Use for detailed internal state information during development.
 */
#define sam_debug(cmd, fmt, ...) \
    sam_log_message(SAM_LOG_DEBUG, cmd, __FILE__, __LINE__, fmt, ##__VA_ARGS__)

// =============================================================================
// Memory Allocation Error Helpers
// =============================================================================

/**
 * Check memory allocation result and report error if failed.
 * @param ptr Result of malloc/calloc/realloc
 * @param cmd Command name for error reporting
 * @return ptr if successful, NULL if failed (after reporting error)
 */
static inline void* sam_check_alloc(void *ptr, const char *cmd) {
    if (!ptr) {
        sam_error(cmd, "Memory allocation failed");
    }
    return ptr;
}

/**
 * Safe memory allocation with automatic error reporting.
 * @param size Bytes to allocate
 * @param cmd Command name for error reporting
 * @return allocated memory or NULL on failure (after reporting error)
 */
static inline void* sam_malloc(size_t size, const char *cmd) {
    void *ptr = malloc(size);
    return sam_check_alloc(ptr, cmd);
}

/**
 * Safe zero-initialized memory allocation with automatic error reporting.
 * @param count Number of elements
 * @param size Size of each element
 * @param cmd Command name for error reporting
 * @return allocated memory or NULL on failure (after reporting error)
 */
static inline void* sam_calloc(size_t count, size_t size, const char *cmd) {
    void *ptr = calloc(count, size);
    return sam_check_alloc(ptr, cmd);
}

// =============================================================================
// Function Declarations
// =============================================================================

/**
 * Core logging function - use macros above instead of calling directly.
 * @param level Severity level of the message
 * @param cmd Command name where error occurred
 * @param file Source file name
 * @param line Line number in source
 * @param fmt Printf-style format string
 */
void sam_log_message(sam_log_level_t level, const char *cmd, 
                     const char *file, int line, const char *fmt, ...);

/**
 * Set global logging level (controls which messages are displayed).
 * @param level Minimum level to display
 */
void sam_set_log_level(sam_log_level_t level);

/**
 * Enable or disable verbose logging output.
 * @param enabled true to show file/line information, false for clean output
 */
void sam_set_verbose_logging(bool enabled);

#ifdef __cplusplus
}
#endif

#endif // SAM_ERROR_H
