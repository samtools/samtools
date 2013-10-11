#ifndef __SAMTOOLS_H__
#define __SAMTOOLS_H__

#include <stdarg.h>

const char *samtools_version(void);
void error(const char *format, ...);

#endif
