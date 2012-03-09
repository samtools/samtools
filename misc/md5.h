/*
  This file is adapted from a program in this page:

  http://www.fourmilab.ch/md5/

  The original source code does not work on 64-bit machines due to the
  wrong typedef "uint32". I also added prototypes.

  -lh3
 */

#ifndef MD5_H
#define MD5_H

/*  The following tests optimise behaviour on little-endian
    machines, where there is no need to reverse the byte order
    of 32 bit words in the MD5 computation.  By default,
    HIGHFIRST is defined, which indicates we're running on a
    big-endian (most significant byte first) machine, on which
    the byteReverse function in md5.c must be invoked. However,
    byteReverse is coded in such a way that it is an identity
    function when run on a little-endian machine, so calling it
    on such a platform causes no harm apart from wasting time. 
    If the platform is known to be little-endian, we speed
    things up by undefining HIGHFIRST, which defines
    byteReverse as a null macro.  Doing things in this manner
    insures we work on new platforms regardless of their byte
    order.  */

#define HIGHFIRST

#if __LITTLE_ENDIAN__ != 0
#undef HIGHFIRST
#endif

#include <stdint.h>

struct MD5Context {
        uint32_t buf[4];
        uint32_t bits[2];
        unsigned char in[64];
};

void MD5Init(struct MD5Context *ctx);
void MD5Update(struct MD5Context *ctx, unsigned char *buf, unsigned len);
void MD5Final(unsigned char digest[16], struct MD5Context *ctx);

/*
 * This is needed to make RSAREF happy on some MS-DOS compilers.
 */
typedef struct MD5Context MD5_CTX;

/*  Define CHECK_HARDWARE_PROPERTIES to have main,c verify
    byte order and uint32_t settings.  */
#define CHECK_HARDWARE_PROPERTIES

#endif /* !MD5_H */
