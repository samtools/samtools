#ifndef SAMTOOLS_H
#define SAMTOOLS_H

const char *samtools_version(void);

void print_error(const char *format, ...);
void print_error_errno(const char *format, ...);

#endif
