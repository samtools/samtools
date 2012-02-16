#include <stdlib.h>

// complicated expression better fits as macro (or inline in C++)
#define ByteOf(x) (((x) >> bitsOffset) & 0xff)

static void radix (short bitsOffset, int N, int *source, int *dest)
{
  // suppressed the need for index as it is reported in count
  int count[256];
  // added temp variables to simplify writing, understanding and compiler optimization job
  // most of them will be allocated as registers
  int *sp, *cp, s, c, i;

  // faster than MemSet
  cp = count;
  for (i = 256; i > 0; --i, ++cp)
    *cp = 0;

  // count occurences of every byte value
  sp = source;
  for (i = N; i > 0; --i, ++sp) {
    cp = count + ByteOf (*sp);
    ++(*cp);
  }

  // transform count into index by summing elements and storing into same array
  s = 0;
  cp = count;
  for (i = 256; i > 0; --i, ++cp) {
    c = *cp;
    *cp = s;
    s += c;
  }

  // fill dest with the right values in the right place
  sp = source;
  for (i = N; i > 0; --i, ++sp) {
    cp = count + ByteOf (*sp);
    dest[*cp] = *sp;
    ++(*cp);
  }
}

// complicated expression better fits as macro (or inline in C++)
/*#define ByteOf(x) (((x) >> bitsOffset) & 0xff)

// replaced byte with bitsOffset to avoid *8 operation in loop
static void radix (short bitsOffset, int N, int *source, int *dest)
{
  // suppressed the need for index as it is reported in count
  int count[256];
  // added temp variables to simplify writing, understanding and compiler optimization job
  // most of them will be allocated as registers
  int *cp, *sp, s, c, i;

  // faster than MemSet
  cp = count;
  for (i = 256; i > 0; --i, ++cp)
    *cp = 0;

  // count occurences of every byte value
  sp = source;
  for (i = N; i > 0; --i, ++sp) {
    cp = count + ByteOf (*sp);
    ++(*cp);
  }

  // transform count into index by summing elements and storing into same array
  s = 0;
  cp = count;
  for (i = 256; i > 0; --i, ++cp) {
    c = *cp;
    *cp = s;
    s += c;
  }

  // fill dest with the right values in the right place
  sp = source;
  for (i = N; i > 0; --i, ++sp) {
    s = *sp;
    cp = count + ByteOf (s);
    dest[*cp] = s;
    ++(*cp);
  }
}*/

static void radix_sort (int *source, int N)
{
  // allocate heap memory to avoid the need of additional parameter
  int *temp = malloc (N * sizeof (int));

  radix (0, N, source, temp);
  radix (8, N, temp, source);
  radix (16, N, source, temp);
  radix (24, N, temp, source);

  free (temp);
}
/*
static int check_order (int *data, int N)
{
  // only signal errors if any (should not be)
  --N;
  for ( ; N > 0; --N, ++data)
    if (data[0] > data[1]) {
      printf("%d %d",data[0],data[1]);
      return N;
    }

  return -1;
}
*/
