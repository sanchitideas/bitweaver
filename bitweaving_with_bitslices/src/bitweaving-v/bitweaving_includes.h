/*******************************************************************************
*                                                                              *
* Vertical Bit Parallel storage & column scan technqiues                       *
*                                                                              *
* This program generates random integers over a uniform distribution and stores*
* them in VBP. Then we can run operations corresponding to Vertical            *
* BitWeaving.                                                                  *
*                                                                              *
* The number of bit slices are equal to the number of the bits in the code.    *
* The bit-group size is thus 1. The bit-slices are stored separately,          *
* but within each bit-group, all the bits (pertaining to some particular bit of*
* column tuples) are contiguous                                                *
*                                                                              *
* Author: Sanchit Jain <sanchit@cs.wisc.edu>                                   *
*                                                                              *
*******************************************************************************/

#ifndef BITWEAVING_INCLUDES_
#define BITWEAVING_INCLUDES_
#include <cstdint>
#include <stdlib.h>
#include <cstdio>
#include <iostream>
#include <vector>
#include <immintrin.h>
#include <random>
#include <string>
#include <math.h>
#include <unistd.h>
#include <chrono>
#include <sys/mman.h>
#include <fcntl.h>
#define WORD_SIZE 64
#define BYTES_IN_WORD 8
#define PAGE_SIZE 4096
#define ADDR (void *)(0x0UL)
#define FLAGS (MAP_PRIVATE | MAP_ANONYMOUS | MAP_HUGETLB | MAP_LOCKED)
#define PROTECTION (PROT_READ | PROT_WRITE)
#define MMAPED_DATA (4UL*1024*1024*1024)
// Comment the following line to disable cache-line prefetching
#define PREFETCH 1


static inline __attribute__((always_inline)) \
  void PREFETCH_ONLY_INTO_L2 (void *x) {
  asm volatile("prefetcht1 %0" : : "m" (*(uint64_t *)x));
}

static inline __attribute__((always_inline)) \
  void PREFETCH_INTO_L1_AND_L2 (void *x) {
  asm volatile("prefetcht0 %0" : : "m" (*(uint64_t *)x)); \
}

static inline __attribute__((always_inline)) \
  void PREFETCH_INTO_L1D_NT (void *x) {
  asm volatile("prefetchnta %0" : : "m" (*(uint64_t *)x));
}

static inline __attribute__((always_inline)) \
  void PREFETCH_INTO_L1D_FOR_WRITES (void *x) {
  asm volatile("prefetchw %0" : : "m" (*(uint64_t *)x));
}

static inline uint64_t powof(uint64_t x, uint64_t y) {
  uint64_t retval = 1;
  while (y--) {
    retval *= x;
  }
  return retval;
}

static inline uint64_t rdtscp( uint32_t & aux ) {
    uint64_t rax,rdx;
    asm volatile ( "rdtscp\n" : "=a" (rax), "=d" (rdx), "=c" (aux) : : );
    return (rdx << 32) + rax;
}

/*****************************************************************************
 * numOfBits()                                                               *
 *                                                                           *
 * This function computes the number of bits required to encode an integer   *
 ****************************************************************************/
static inline int numOfBits(int number) {
  int i = 1;
  while ((number = number / 2) != 0) {
    i++;
  }
  return i;
} // end numOfBits()


/*****************************************************************************
 * decimalToBinary()                                                         *
 *                                                                           *
 * Wrote this code to convert an unsigned WORD_SIZE-bit decimal into its     *
 * binary representation because even the Microsoft Windows calculator       *
 * cannot operate on unsigned WORD_SIZE-bit numbers in Programmer mode       *
 * Nor do any available online tools.                                        *
 * They don't even mention this limitation                                   *
 ****************************************************************************/
static inline void decimalToBinary(uint64_t decimalNum) {
  printf("\n");
  for (int i = 0; i <= 63; i++) {
    uint64_t mask = ((uint64_t)1 << (63 - i));
    printf("%llu", ((mask & decimalNum) >> (63 - i)));
  }
  printf("\n");
} // end decimalToBinary()

class BitWeaver {
  private:
    uint64_t** codeVectors;
    int codeLength;
    int cardinality;
  public:
    BitWeaver(std::vector<std::string> filenames, uint64_t cardinality);
    BitWeaver(std::string filename, uint64_t cardinality, int codeSize);
    uint64_t** generateRandomNumbers(uint64_t count, 
      uint64_t start, uint64_t end);
    void lessThan(uint64_t comparedTo, uint64_t* resultBitMap,
      uint64_t cardinality);
    void greaterThan(uint64_t comparedTo, uint64_t* resultBitMap,
      uint64_t cardinality);
    void greaterThanOrEqualTo(uint64_t comparedTo, uint64_t* resultBitMap,
      uint64_t cardinality);
    void lessThanOrEqualTo(uint64_t comparedTo, uint64_t* resultBitMap,
      uint64_t cardinality);
    void equality(uint64_t comparedTo, uint64_t* resultBitMap,
      uint64_t cardinality);
    void inequality(uint64_t comparedTo, uint64_t* resultBitMap,
      uint64_t cardinality);
    void between(uint64_t value1, uint64_t value2,
      uint64_t* resultBitMap, uint64_t cardinality);
    void complexLessThan(uint64_t comparedTo, uint64_t* resultBitMap,
      uint64_t cardinality);
    void complexGreaterThan(uint64_t comparedTo, uint64_t* resultBitMap,
      uint64_t cardinality);
    void complexGreaterThanOrEqualTo(uint64_t comparedTo, 
      uint64_t* resultBitMap, uint64_t cardinality);
    void complexLessThanOrEqualTo(uint64_t comparedTo, 
      uint64_t* resultBitMap, uint64_t cardinality);
    void complexEquality(uint64_t comparedTo, uint64_t* resultBitMap,
      uint64_t cardinality);
    void complexInequality(uint64_t comparedTo, uint64_t* resultBitMap,
      uint64_t cardinality);
    void complexBetween(uint64_t value1, uint64_t value2,
      uint64_t* resultBitMap, uint64_t cardinality);
    // deconstructor should be written if smart pointers aren't used
};

#endif 