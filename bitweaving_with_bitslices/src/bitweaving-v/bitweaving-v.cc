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

#include "bitweaving_includes.h"


BitWeaver::BitWeaver(std::string filename, uint64_t cardinality, int codeSize){
  codeLength = codeSize;
  codeVectors = new uint64_t*[codeLength];
  // number of CPU words required to store a bit slice
  uint64_t num_of_cpu_words =  (cardinality / WORD_SIZE) ? 
    ((cardinality % WORD_SIZE) ? cardinality/WORD_SIZE + 1 : 
    cardinality/WORD_SIZE) : 1;
  // we need to make it a multiple of page size, in order to mmap
  int numOfFillerBytes = ((num_of_cpu_words * BYTES_IN_WORD) % PAGE_SIZE) ? \
    (PAGE_SIZE - ((num_of_cpu_words * BYTES_IN_WORD) % PAGE_SIZE)) : 0;

  int fileFD = open(filename.c_str(), O_RDONLY);
  int i = 0;
  int j = ((codeLength - 8) > 0) ? 8 : codeLength;
  for (i = 0; i < j; i++) {
    codeVectors[i] = (uint64_t *)(aligned_alloc(PAGE_SIZE * 16, 
    num_of_cpu_words * BYTES_IN_WORD));
    codeVectors[i] = (uint64_t*) mmap(codeVectors[i], 
      num_of_cpu_words * BYTES_IN_WORD, PROT_READ,
      MAP_PRIVATE, fileFD, 
      i * (num_of_cpu_words * BYTES_IN_WORD + numOfFillerBytes));
    mlock(codeVectors[i], num_of_cpu_words * BYTES_IN_WORD);
  }
  if ((i = codeLength - i) > 0) {
    j = ((codeLength - 16) > 0) ? 16 : codeLength;
    for (i = 8; i < j; i++) {
      codeVectors[i] = (uint64_t *)(aligned_alloc(PAGE_SIZE * 16, 
      num_of_cpu_words * BYTES_IN_WORD + 1024));
      codeVectors[i] += 128;
      codeVectors[i] = (uint64_t*) mmap(codeVectors[i], 
        num_of_cpu_words * BYTES_IN_WORD, PROT_READ, MAP_PRIVATE, fileFD,
        i * (num_of_cpu_words * BYTES_IN_WORD + numOfFillerBytes));
      mlock(codeVectors[i] - 128, 
        num_of_cpu_words * BYTES_IN_WORD + 1024);
    }
    if ((i = codeLength - i) > 0) {
      j = ((codeLength - 24) > 0) ? 24 : codeLength;
      for (i = 16; i < j; i++) {
        codeVectors[i] = (uint64_t *)(aligned_alloc(PAGE_SIZE * 16,
        num_of_cpu_words * BYTES_IN_WORD + 2048));
        codeVectors[i] += 256;
        codeVectors[i] = (uint64_t*) mmap(codeVectors[i], 
          num_of_cpu_words * BYTES_IN_WORD, PROT_READ, MAP_PRIVATE, fileFD,
          i * (num_of_cpu_words * BYTES_IN_WORD + numOfFillerBytes));
        mlock(codeVectors[i] - 256, 
          num_of_cpu_words * BYTES_IN_WORD + 2048);
      }
      if ((i = codeLength - i) > 0) {
        for (i = 24; i < codeLength; i++) {
          codeVectors[i] = (uint64_t *)(aligned_alloc(PAGE_SIZE * 16,
            num_of_cpu_words * BYTES_IN_WORD + 3072));
          codeVectors[i] += 384;
          codeVectors[i] = (uint64_t*) mmap(codeVectors[i], 
            num_of_cpu_words * BYTES_IN_WORD, PROT_READ,
            MAP_PRIVATE, fileFD, 
            i * \
            (num_of_cpu_words * BYTES_IN_WORD + numOfFillerBytes));
          mlock(codeVectors[i] - 384, 
            num_of_cpu_words * BYTES_IN_WORD + 3072);
        }          
      }          
    }      
  }
} // end constructor


BitWeaver::BitWeaver(std::vector<std::string> filenames, uint64_t cardinality){
  codeLength = filenames.size();
  codeVectors = new uint64_t*[codeLength];
  // number of CPU words required to store a bit slice
  uint64_t num_of_cpu_words =  (cardinality / WORD_SIZE) ? 
    ((cardinality % WORD_SIZE) ? cardinality/WORD_SIZE + 1 : 
    cardinality/WORD_SIZE) : 1;

  int i = 0;
  int j = ((codeLength - 8) > 0) ? 8 : codeLength;
  for (i = 0; i < j; i++) {
    codeVectors[i] = (uint64_t *)(aligned_alloc(PAGE_SIZE * 16, 
    num_of_cpu_words * 8));
      int fileFD = open(filenames[i].c_str(), O_RDONLY);
    codeVectors[i] = (uint64_t*) mmap(codeVectors[i], 
      num_of_cpu_words * 8, PROT_READ,
      MAP_PRIVATE, fileFD, 0);
    mlock(codeVectors[i], num_of_cpu_words * 8);
  }
  if ((i = codeLength - i) > 0) {
    j = ((codeLength - 16) > 0) ? 16 : codeLength;
    for (i = 8; i < j; i++) {
      codeVectors[i] = (uint64_t *)(aligned_alloc(PAGE_SIZE * 16, 
      num_of_cpu_words * 8 + 1024));
      codeVectors[i] += 128;
        int fileFD = open(filenames[i].c_str(), O_RDONLY);
      codeVectors[i] = (uint64_t*) mmap(codeVectors[i], 
        num_of_cpu_words * 8, PROT_READ,
        MAP_PRIVATE, fileFD, 0);
      mlock(codeVectors[i] - 128, num_of_cpu_words * 8 + 1024);
    }
    if ((i = codeLength - i) > 0) {
      j = ((codeLength - 24) > 0) ? 24 : codeLength;
      for (i = 16; i < j; i++) {
        codeVectors[i] = (uint64_t *)(aligned_alloc(PAGE_SIZE * 16,
        num_of_cpu_words * 8 + 2048));
        codeVectors[i] += 256;
          int fileFD = open(filenames[i].c_str(), O_RDONLY);
        codeVectors[i] = (uint64_t*)mmap(codeVectors[i], 
          num_of_cpu_words * 8, PROT_READ, MAP_PRIVATE,
          fileFD, 0);
        mlock(codeVectors[i] - 256, num_of_cpu_words * 8 + 2048);
      }
      if ((i = codeLength - i) > 0) {
        for (i = 24; i < codeLength; i++) {
          codeVectors[i] = (uint64_t *)(aligned_alloc(PAGE_SIZE * 16,
            num_of_cpu_words * 8 + 3072));
          codeVectors[i] += 384;
            int fileFD = open(filenames[i].c_str(), O_RDONLY);
          codeVectors[i] = (uint64_t*) mmap(codeVectors[i], 
            num_of_cpu_words * 8, 
            PROT_READ, MAP_PRIVATE, fileFD, 0);
          mlock(codeVectors[i] - 384, num_of_cpu_words * 8 + 3072);
        }          
      }          
    }      
  }
}

/*******************************************************************************
 * This function not only generates random numbers, but returns a pointer to   *
 * the VBP storage                                                             *
 *                                                                             *
 ******************************************************************************/
uint64_t** BitWeaver::generateRandomNumbers(uint64_t count, 
  unsigned long start, uint64_t end) {
  int codeVectorsize = WORD_SIZE;
  codeLength = numOfBits(end);
  /*************************************************************************** 
   *random numbers generated over a uniform distribution                     *
   * This function reference has been adapted from C++ reference at          *
   * www.cplusplus.com/reference/random/uniform_int_distribution/operator()/ *
     **************************************************************************/
  std::default_random_engine generator;

  std::uniform_int_distribution<int> distribution(start, end);
  __attribute((aligned(WORD_SIZE))) int* numbers = 
    (int*)(aligned_alloc(WORD_SIZE, count*sizeof(int)));

  int numOfCodesInAWord = WORD_SIZE;
  int numOfcodeVectorsPerSegment = codeLength;
  int bytesInOneSegment = codeLength * WORD_SIZE / 8;
  uint64_t numberOfSegments = (count / WORD_SIZE) ? 
    ((count % WORD_SIZE) ? count/WORD_SIZE + 1 : 
    count/WORD_SIZE) : 1;
  codeVectors = new uint64_t*[codeLength];

  int i = 0;
  int j = ((codeLength - 8) > 0) ? 8 : codeLength;
  for (i = 0; i < j; i++) {
    codeVectors[i] = (uint64_t *)(aligned_alloc(PAGE_SIZE * 16, 
    (bytesInOneSegment * numberOfSegments) /codeLength));
    mlock(codeVectors[i], (bytesInOneSegment * numberOfSegments) / codeLength);
  }
  if ((i = codeLength - i) > 0) {
    j = ((codeLength - 16) > 0) ? 16 : codeLength;
    for (i = 8; i < j; i++) {
      codeVectors[i] = 
        (uint64_t *)(aligned_alloc(PAGE_SIZE * 16, 
        (bytesInOneSegment * numberOfSegments) /codeLength + 1024));
      mlock(codeVectors[i], (bytesInOneSegment * numberOfSegments) / 
        codeLength + 1024);
      codeVectors[i] = codeVectors[i] + 128;
    }
    if ((i = codeLength - i) > 0) {
      j = ((codeLength - 24) > 0) ? 24 : codeLength;
      for (i = 16; i < j; i++) {
        codeVectors[i] = 
          (uint64_t *)(aligned_alloc(PAGE_SIZE * 16, 
          (bytesInOneSegment * numberOfSegments) /codeLength + 2048));
        mlock(codeVectors[i], (bytesInOneSegment * numberOfSegments) / 
          codeLength + 2048);
        codeVectors[i] = codeVectors[i] + 256;
      }
      if ((i = codeLength - i) > 0) {
        for (i = 24; i < codeLength; i++) {
          codeVectors[i] = 
            (uint64_t *)(aligned_alloc(PAGE_SIZE * 16,
            (bytesInOneSegment * numberOfSegments) / 
            codeLength + 3072));
          mlock(codeVectors[i], (bytesInOneSegment * numberOfSegments) / 
            codeLength + 3072);
          codeVectors[i] = codeVectors[i] + 384;
        }          
      }          
    }      
  }


  for (uint64_t i = 0; i < count; i++) {
      int num = distribution(generator);
      //printf("\n\n%d\n\n", num);
      numbers[i] = num;
   }

   // group size is 1. All bits are stored separately for codeVectors

   uint64_t lastSegmentNumber = numberOfSegments - 1;
   uint64_t codeNo = 0;
   int numOfCodesToRead = WORD_SIZE;
    for (int bitNumber = 0; bitNumber < codeLength; bitNumber++) {
     for (uint64_t segNo = 0; segNo < numberOfSegments; segNo++) {
       uint64_t currentWord = 0;
       // extract all bits
       if (segNo == lastSegmentNumber) {
         // the last segment might have fewer than 64 numbers
         int numCodesInLastSegment = count % WORD_SIZE;
         if (numCodesInLastSegment) {
           numOfCodesToRead = count % WORD_SIZE;
         }
       }
       for (int i = 0; i < numOfCodesToRead; i++) {
         // extract the (codeLength - bitNumber) bit from each number
         currentWord |= 
           (((((uint64_t)1 << 
           (codeLength - 1 - bitNumber))
           & numbers[segNo * numOfCodesInAWord + i]) >> 
           (codeLength - 1 - bitNumber)) << (WORD_SIZE - i - 1));
       }
       codeVectors[bitNumber][segNo] = currentWord;
     }
   }
   //decimalToBinary(codeVectors[0][0]);
   return codeVectors;
} // end generateRandomNumbers()




void BitWeaver::equality(uint64_t comparedTo, uint64_t* resultBitMap, 
  uint64_t cardinality) {
  int numOfCodesInAWord = WORD_SIZE;
  int bytesInOneSegment = codeLength * WORD_SIZE / 8;
  uint64_t numberOfSegments = (cardinality / WORD_SIZE) ? 
    ((cardinality % WORD_SIZE) ? cardinality/WORD_SIZE + 1 : 
    cardinality/WORD_SIZE) : 1;

  //bit representation of comparator
  uint64_t comparedToMasks[codeLength];
  for (int i = 0; i < codeLength; i++) {
    uint64_t comparedToMask = 
      ((1 << (codeLength - i - 1)) & comparedTo) >> 
      (codeLength - i - 1);
    if (comparedToMask) {
      comparedToMasks[i] = -1;
    } else {
      comparedToMasks[i] = 0;
    }
  }

  int prefixDistance = (codeLength > 16) ? 15 : 7; 
  __m512i zeroBitMapVector = _mm512_set1_epi64(0);
  __mmask8 pruningCheckVector;

  for (uint64_t segCount = 0; segCount < numberOfSegments; 
    segCount += 8) {
    __m512i ltBitMapVector = _mm512_set1_epi64(0);
    __m512i eqBitMapVector = _mm512_set1_epi64((uint64_t)-1);
    for (int j = 0; j < codeLength; j++) {
      #ifdef PREFETCH
        if (j < (codeLength - prefixDistance)) {
          PREFETCH_INTO_L1_AND_L2(codeVectors[j + prefixDistance] + segCount);
        } else {
          PREFETCH_INTO_L1_AND_L2(codeVectors[j - (codeLength - \
            prefixDistance)] + segCount + 8);
        }
      #endif
      __m512i avxVector = _mm512_load_epi64(codeVectors[j] + segCount);
      uint64_t attributeMask = comparedToMasks[j];
      __m512i attributeMaskVector = _mm512_set1_epi64(attributeMask);
      __m512i temp = _mm512_xor_epi64(avxVector, attributeMaskVector);
      eqBitMapVector = _mm512_andnot_epi64(temp, eqBitMapVector);

      // check for early pruning
      pruningCheckVector = _mm512_cmp_epu64_mask(eqBitMapVector,
        zeroBitMapVector, _MM_CMPINT_EQ);

      if (pruningCheckVector == -1) {
        break;
      }

    } // end inner for loop
    _mm512_store_epi64(resultBitMap + segCount, eqBitMapVector);        
  } // end outer for loop    
} // end equality()


void BitWeaver::complexEquality(uint64_t comparedTo, uint64_t* resultBitMap,
  uint64_t cardinality) {
  int numOfCodesInAWord = WORD_SIZE;
  int bytesInOneSegment = codeLength * WORD_SIZE / 8;
  uint64_t numberOfSegments = (cardinality / WORD_SIZE) ? 
    ((cardinality % WORD_SIZE) ? cardinality/WORD_SIZE + 1 : 
    cardinality/WORD_SIZE) : 1;

  //bit representation of comparator
  uint64_t comparedToMasks[codeLength];
  for (int i = 0; i < codeLength; i++) {
    uint64_t comparedToMask = 
      ((1 << (codeLength - i - 1)) & comparedTo) >> 
      (codeLength - i - 1);
    if (comparedToMask) {
      comparedToMasks[i] = -1;
    } else {
      comparedToMasks[i] = 0;
    }
  }

  int prefixDistance = (codeLength > 16) ? 15 : 7; 
  __m512i zeroBitMapVector = _mm512_set1_epi64(0);
  __mmask8 pruningCheckVector;

  for (uint64_t segCount = 0; segCount < numberOfSegments; 
    segCount += 8) {
    __m512i ltBitMapVector = _mm512_set1_epi64(0);
    __m512i eqBitMapVector = _mm512_load_epi64(resultBitMap + segCount);
    for (int j = 0; j < codeLength; j++) {
      #ifdef PREFETCH
        if (j < (codeLength - prefixDistance)) {
          PREFETCH_INTO_L1_AND_L2(codeVectors[j + prefixDistance] + 
            segCount);
        } else {
          PREFETCH_INTO_L1_AND_L2(codeVectors[j - (codeLength - \
            prefixDistance)] + segCount + 8);
        }
      #endif
      __m512i avxVector = _mm512_load_epi64(codeVectors[j] + segCount);
      uint64_t attributeMask = comparedToMasks[j];
      __m512i attributeMaskVector = _mm512_set1_epi64(attributeMask);
      __m512i temp = _mm512_xor_epi64(avxVector, attributeMaskVector);
      eqBitMapVector = _mm512_andnot_epi64(temp, eqBitMapVector);

      // check for early pruning
      pruningCheckVector = _mm512_cmp_epu64_mask(eqBitMapVector,
        zeroBitMapVector, _MM_CMPINT_EQ);

      if (pruningCheckVector == -1) {
        break;
      }

    } // end inner for loop
    _mm512_store_epi64(resultBitMap + segCount, eqBitMapVector);        
  } // end outer for loop    
} // end complexEquality()


void BitWeaver::inequality(uint64_t comparedTo, uint64_t* resultBitMap, 
  uint64_t cardinality) {
  int numOfCodesInAWord = WORD_SIZE;
  int bytesInOneSegment = codeLength * WORD_SIZE / 8;
  uint64_t numberOfSegments = (cardinality / WORD_SIZE) ? 
    ((cardinality % WORD_SIZE) ? cardinality/WORD_SIZE + 1 : 
    cardinality/WORD_SIZE) : 1;

  //bit representation of comparator
  uint64_t comparedToMasks[codeLength];
  for (int i = 0; i < codeLength; i++) {
    uint64_t comparedToMask = 
      ((1 << (codeLength - i - 1)) & comparedTo) >> 
      (codeLength - i - 1);
    comparedToMasks[i] = (comparedToMask) ? -1 : 0;
  }

  int prefixDistance = (codeLength > 16) ? 15 : 7; 
  __m512i zeroBitMapVector = _mm512_set1_epi64(0);
  __mmask8 pruningCheckVector;

  for (uint64_t segCount = 0; segCount < numberOfSegments; 
    segCount += 8) {
    __m512i ltBitMapVector = _mm512_set1_epi64(0);
    __m512i eqBitMapVector = _mm512_set1_epi64((uint64_t)-1);
    for (int j = 0; j < codeLength; j++) {
      #ifdef PREFETCH
        if (j < (codeLength - prefixDistance)) {
          PREFETCH_INTO_L1_AND_L2(codeVectors[j + prefixDistance] + segCount);
        } else {
          PREFETCH_INTO_L1_AND_L2(codeVectors[j - (codeLength - \
            prefixDistance)] + segCount + 8);
        }
      #endif
      __m512i avxVector = _mm512_load_epi64(codeVectors[j] + segCount);
      uint64_t attributeMask = comparedToMasks[j];
      __m512i attributeMaskVector = _mm512_set1_epi64(attributeMask);
      __m512i temp = _mm512_xor_epi64(avxVector, attributeMaskVector);
      eqBitMapVector = _mm512_andnot_epi64(temp, eqBitMapVector);

      // check for early pruning
      pruningCheckVector = _mm512_cmp_epu64_mask(eqBitMapVector,
        zeroBitMapVector, _MM_CMPINT_EQ);

      if (pruningCheckVector == -1) {
        break;
      }


    } // end inner for loop
    __m512i temp = _mm512_set1_epi64((uint64_t)-1);
    eqBitMapVector = _mm512_andnot_epi64(eqBitMapVector, temp);
    _mm512_store_epi64(resultBitMap + segCount, eqBitMapVector);        
  } // end outer for loop          
} // end inequality()


void BitWeaver::complexInequality(uint64_t comparedTo, uint64_t* resultBitMap,
  uint64_t cardinality) {
  int numOfCodesInAWord = WORD_SIZE;
  int bytesInOneSegment = codeLength * WORD_SIZE / 8;
  uint64_t numberOfSegments = (cardinality / WORD_SIZE) ? 
    ((cardinality % WORD_SIZE) ? cardinality/WORD_SIZE + 1 : 
    cardinality/WORD_SIZE) : 1;

  //bit representation of comparator
  uint64_t comparedToMasks[codeLength];
  for (int i = 0; i < codeLength; i++) {
    uint64_t comparedToMask = 
      ((1 << (codeLength - i - 1)) & comparedTo) >> 
      (codeLength - i - 1);
    comparedToMasks[i] = (comparedToMask) ? -1 : 0;
  }

  int prefixDistance = (codeLength > 16) ? 15 : 7; 
  __m512i zeroBitMapVector = _mm512_set1_epi64(0);
  __mmask8 pruningCheckVector;

  for (uint64_t segCount = 0; segCount < numberOfSegments; 
    segCount += 8) {
    __m512i ltBitMapVector = _mm512_set1_epi64(0);
    __m512i eqBitMapVector = _mm512_load_epi64(resultBitMap + segCount);
    for (int j = 0; j < codeLength; j++) {
      #ifdef PREFETCH
        if (j < (codeLength - prefixDistance)) {
          PREFETCH_INTO_L1_AND_L2(codeVectors[j + prefixDistance] + segCount);
        } else {
          PREFETCH_INTO_L1_AND_L2(codeVectors[j - (codeLength - \
            prefixDistance)] + segCount + 8);
        }
      #endif
      __m512i avxVector = _mm512_load_epi64(codeVectors[j] + segCount);
      uint64_t attributeMask = comparedToMasks[j];
      __m512i attributeMaskVector = _mm512_set1_epi64(attributeMask);
      __m512i temp = _mm512_xor_epi64(avxVector, attributeMaskVector);
      eqBitMapVector = _mm512_andnot_epi64(temp, eqBitMapVector);

      // check for early pruning
      pruningCheckVector = _mm512_cmp_epu64_mask(eqBitMapVector,
        zeroBitMapVector, _MM_CMPINT_EQ);

      if (pruningCheckVector == -1) {
        break;
      }


    } // end inner for loop
    __m512i temp = _mm512_set1_epi64((uint64_t)-1);
    eqBitMapVector = _mm512_andnot_epi64(eqBitMapVector, temp);
    _mm512_store_epi64(resultBitMap + segCount, eqBitMapVector);        
  } // end outer for loop          
} // end complexInequality()


void BitWeaver::lessThan(uint64_t comparedTo, uint64_t* resultBitMap,
  uint64_t cardinality) {

  int numOfCodesInAWord = WORD_SIZE;
  int bytesInOneSegment = codeLength * WORD_SIZE / 8;
  uint64_t numberOfSegments = (cardinality / WORD_SIZE) ? 
    ((cardinality % WORD_SIZE) ? cardinality/WORD_SIZE + 1 : 
    cardinality/WORD_SIZE) : 1;

  //bit representation of comparator
  uint64_t comparedToMasks[codeLength];
  for (int i = 0; i < codeLength; i++) {
    uint64_t comparedToMask = 
      ((1 << (codeLength - i - 1)) & comparedTo) >> 
      (codeLength - i - 1);
    comparedToMasks[i] = (comparedToMask) ? -1 : 0;
  }

  int prefixDistance = (codeLength > 16) ? 15 : 7; 
  __m512i zeroBitMapVector = _mm512_set1_epi64(0);
  __mmask8 pruningCheckVector;

  for (uint64_t segCount = 0; segCount < numberOfSegments; 
    segCount += 8) {
    __m512i ltBitMapVector = _mm512_set1_epi64(0);
    __m512i eqBitMapVector = _mm512_set1_epi64((uint64_t)-1);
    for (int j = 0; j < codeLength; j++) {
      #ifdef PREFETCH
        if (j < (codeLength - prefixDistance)) {
          PREFETCH_INTO_L1_AND_L2(codeVectors[j + prefixDistance] + segCount);
        } else {
          PREFETCH_INTO_L1_AND_L2(codeVectors[j - (codeLength - \
            prefixDistance)] + segCount + 8);
        }
      #endif
      __m512i avxVector = _mm512_load_epi64(codeVectors[j] + segCount);
      uint64_t attributeMask = comparedToMasks[j];
      __m512i attributeMaskVector = _mm512_set1_epi64(attributeMask);
      __m512i temp = _mm512_andnot_epi64(avxVector, attributeMaskVector);
      temp = _mm512_and_epi64(temp, eqBitMapVector);
      ltBitMapVector = _mm512_or_epi64(temp, ltBitMapVector);
      temp = _mm512_xor_epi64(avxVector, attributeMaskVector);
      eqBitMapVector = _mm512_andnot_epi64(temp, eqBitMapVector);


      // check for early pruning
      pruningCheckVector = _mm512_cmp_epu64_mask(eqBitMapVector,
        zeroBitMapVector, _MM_CMPINT_EQ);

      if (pruningCheckVector == -1) {
        break;
      }


    }
    _mm512_store_epi64(resultBitMap + segCount, ltBitMapVector);        
  }      

} // end lessThan()

void BitWeaver::complexLessThan(uint64_t comparedTo, uint64_t* resultBitMap,
  uint64_t cardinality) {

  int numOfCodesInAWord = WORD_SIZE;
  int bytesInOneSegment = codeLength * WORD_SIZE / 8;
  uint64_t numberOfSegments = (cardinality / WORD_SIZE) ? 
    ((cardinality % WORD_SIZE) ? cardinality/WORD_SIZE + 1 : 
    cardinality/WORD_SIZE) : 1;

  //bit representation of comparator
  uint64_t comparedToMasks[codeLength];
  for (int i = 0; i < codeLength; i++) {
    uint64_t comparedToMask = 
      ((1 << (codeLength - i - 1)) & comparedTo) >> 
      (codeLength - i - 1);
    comparedToMasks[i] = (comparedToMask) ? -1 : 0;
  }

  int prefixDistance = (codeLength > 16) ? 15 : 7; 
  __m512i zeroBitMapVector = _mm512_set1_epi64(0);
  __mmask8 pruningCheckVector;

  for (uint64_t segCount = 0; segCount < numberOfSegments; 
    segCount += 8) {
    __m512i ltBitMapVector = _mm512_set1_epi64(0);
    __m512i eqBitMapVector = _mm512_load_epi64(resultBitMap + segCount);
    for (int j = 0; j < codeLength; j++) {
      #ifdef PREFETCH
        if (j < (codeLength - prefixDistance)) {
          PREFETCH_INTO_L1_AND_L2(codeVectors[j + prefixDistance] + segCount);
        } else {
          PREFETCH_INTO_L1_AND_L2(codeVectors[j - (codeLength - \
            prefixDistance)] + segCount + 8);
        }
      #endif
      __m512i avxVector = _mm512_load_epi64(codeVectors[j] + segCount);
      uint64_t attributeMask = comparedToMasks[j];
      __m512i attributeMaskVector = _mm512_set1_epi64(attributeMask);
      __m512i temp = _mm512_andnot_epi64(avxVector, attributeMaskVector);
      temp = _mm512_and_epi64(temp, eqBitMapVector);
      ltBitMapVector = _mm512_or_epi64(temp, ltBitMapVector);
      temp = _mm512_xor_epi64(avxVector, attributeMaskVector);
      eqBitMapVector = _mm512_andnot_epi64(temp, eqBitMapVector);


      // check for early pruning
      pruningCheckVector = _mm512_cmp_epu64_mask(eqBitMapVector,
        zeroBitMapVector, _MM_CMPINT_EQ);

      if (pruningCheckVector == -1) {
        break;
      }


    }
    _mm512_store_epi64(resultBitMap + segCount, ltBitMapVector);        
  }      

} // end complexLessThan()

void BitWeaver::greaterThan(uint64_t comparedTo, uint64_t* resultBitMap,
  uint64_t cardinality) {
  int numOfCodesInAWord = WORD_SIZE;
  int bytesInOneSegment = codeLength * WORD_SIZE / 8;
  uint64_t numberOfSegments = (cardinality / WORD_SIZE) ? 
    ((cardinality % WORD_SIZE) ? cardinality/WORD_SIZE + 1 : 
    cardinality/WORD_SIZE) : 1;

  //bit representation of comparator
  uint64_t comparedToMasks[codeLength];
  for (int i = 0; i < codeLength; i++) {
    uint64_t comparedToMask = 
      ((1 << (codeLength - i - 1)) & comparedTo) >> 
      (codeLength - i - 1);
    comparedToMasks[i] = (comparedToMask) ? -1 : 0;
  }

  int prefixDistance = (codeLength > 16) ? 15 : 7; 
  __m512i zeroBitMapVector = _mm512_set1_epi64(0);
  __mmask8 pruningCheckVector;

  for (uint64_t segCount = 0; segCount < numberOfSegments; 
    segCount += 8) {
    __m512i gtBitMapVector = _mm512_set1_epi64(0);
    __m512i eqBitMapVector = _mm512_set1_epi64((uint64_t)-1);
    for (int j = 0; j < codeLength; j++) {
      #ifdef PREFETCH
        if (j < (codeLength - prefixDistance)) {
          PREFETCH_INTO_L1_AND_L2(codeVectors[j + prefixDistance] + segCount);
        } else {
          PREFETCH_INTO_L1_AND_L2(codeVectors[j - (codeLength - \
            prefixDistance)] + segCount + 8);
        }
      #endif
      __m512i avxVector = _mm512_load_epi64(codeVectors[j] + segCount);
      uint64_t attributeMask = comparedToMasks[j];
      __m512i attributeMaskVector = _mm512_set1_epi64(attributeMask);
      __m512i temp = _mm512_andnot_epi64(attributeMaskVector, avxVector);
      temp = _mm512_and_epi64(temp, eqBitMapVector);
      gtBitMapVector = _mm512_or_epi64(temp, gtBitMapVector);
      temp = _mm512_xor_epi64(avxVector, attributeMaskVector);
      eqBitMapVector = _mm512_andnot_epi64(temp, eqBitMapVector);


      // check for early pruning
      pruningCheckVector = _mm512_cmp_epu64_mask(eqBitMapVector,
        zeroBitMapVector, _MM_CMPINT_EQ);

      if (pruningCheckVector == -1) {
        break;
      }

    }
    _mm512_store_epi64(resultBitMap + segCount, gtBitMapVector);        
  }    

} // end greaterThan()


void BitWeaver::complexGreaterThan(uint64_t comparedTo, 
  uint64_t* resultBitMap, uint64_t cardinality) {
  int numOfCodesInAWord = WORD_SIZE;
  int bytesInOneSegment = codeLength * WORD_SIZE / 8;
  uint64_t numberOfSegments = (cardinality / WORD_SIZE) ? 
    ((cardinality % WORD_SIZE) ? cardinality/WORD_SIZE + 1 : 
    cardinality/WORD_SIZE) : 1;

  //bit representation of comparator
  uint64_t comparedToMasks[codeLength];
  for (int i = 0; i < codeLength; i++) {
    uint64_t comparedToMask = 
      ((1 << (codeLength - i - 1)) & comparedTo) >> 
      (codeLength - i - 1);
    comparedToMasks[i] = (comparedToMask) ? -1 : 0;
  }

  int prefixDistance = (codeLength > 16) ? 15 : 7; 
  __m512i zeroBitMapVector = _mm512_set1_epi64(0);
  __mmask8 pruningCheckVector;

  for (uint64_t segCount = 0; segCount < numberOfSegments; 
    segCount += 8) {
    __m512i gtBitMapVector = _mm512_set1_epi64(0);
    __m512i eqBitMapVector = _mm512_load_epi64(resultBitMap + segCount);
    for (int j = 0; j < codeLength; j++) {
      #ifdef PREFETCH
        if (j < (codeLength - prefixDistance)) {
          PREFETCH_INTO_L1_AND_L2(codeVectors[j + prefixDistance] + segCount);
        } else {
          PREFETCH_INTO_L1_AND_L2(codeVectors[j - (codeLength - \
            prefixDistance)] + segCount + 8);
        }
    #endif
      __m512i avxVector = _mm512_load_epi64(codeVectors[j] + segCount);
      uint64_t attributeMask = comparedToMasks[j];
      __m512i attributeMaskVector = _mm512_set1_epi64(attributeMask);
      __m512i temp = _mm512_andnot_epi64(attributeMaskVector, avxVector);
      temp = _mm512_and_epi64(temp, eqBitMapVector);
      gtBitMapVector = _mm512_or_epi64(temp, gtBitMapVector);
      temp = _mm512_xor_epi64(avxVector, attributeMaskVector);
      eqBitMapVector = _mm512_andnot_epi64(temp, eqBitMapVector);


      // check for early pruning
      pruningCheckVector = _mm512_cmp_epu64_mask(eqBitMapVector,
        zeroBitMapVector, _MM_CMPINT_EQ);

      if (pruningCheckVector == -1) {
        break;
      }

    }
    _mm512_store_epi64(resultBitMap + segCount, gtBitMapVector);        
  }    

} // end complexGreaterThan()


void BitWeaver::lessThanOrEqualTo(uint64_t comparedTo, uint64_t* resultBitMap,
  uint64_t cardinality) {
  int numOfCodesInAWord = WORD_SIZE;
  int bytesInOneSegment = codeLength * WORD_SIZE / 8;
  uint64_t numberOfSegments = (cardinality / WORD_SIZE) ? 
    ((cardinality % WORD_SIZE) ? cardinality/WORD_SIZE + 1 : 
    cardinality/WORD_SIZE) : 1;

  //bit representation of comparator
  uint64_t comparedToMasks[codeLength];
  for (int i = 0; i < codeLength; i++) {
    uint64_t comparedToMask = 
      ((1 << (codeLength - i - 1)) & comparedTo) >> 
      (codeLength - i - 1);
    comparedToMasks[i] = (comparedToMask) ? -1 : 0;
  }

  int prefixDistance = (codeLength > 16) ? 15 : 7; 
  __m512i zeroBitMapVector = _mm512_set1_epi64(0);
  __mmask8 pruningCheckVector;


  for (uint64_t segCount = 0; segCount < numberOfSegments; 
    segCount += 8) {
    __m512i ltBitMapVector = _mm512_set1_epi64(0);
    __m512i eqBitMapVector = _mm512_set1_epi64((uint64_t)-1);
    for (int j = 0; j < codeLength; j++) {
      #ifdef PREFETCH
        if (j < (codeLength - prefixDistance)) {
          PREFETCH_INTO_L1_AND_L2(codeVectors[j + prefixDistance] + segCount);
        } else {
          PREFETCH_INTO_L1_AND_L2(codeVectors[j - (codeLength - \
            prefixDistance)] + segCount + 8);
      }
      #endif
      __m512i avxVector = _mm512_load_epi64(codeVectors[j] + segCount);
      uint64_t attributeMask = comparedToMasks[j];
      __m512i attributeMaskVector = _mm512_set1_epi64(attributeMask);
      __m512i temp = _mm512_andnot_epi64(avxVector, attributeMaskVector);
      temp = _mm512_and_epi64(temp, eqBitMapVector);
      ltBitMapVector = _mm512_or_epi64(temp, ltBitMapVector);
      temp = _mm512_xor_epi64(avxVector, attributeMaskVector);
      eqBitMapVector = _mm512_andnot_epi64(temp, eqBitMapVector);


      // check for early pruning
      pruningCheckVector = _mm512_cmp_epu64_mask(eqBitMapVector,
        zeroBitMapVector, _MM_CMPINT_EQ);

      if (pruningCheckVector == -1) {
        break;
      }

    }
    ltBitMapVector = _mm512_or_epi64(ltBitMapVector, eqBitMapVector);
    _mm512_store_epi64(resultBitMap + segCount, ltBitMapVector);        
  }    

} // end lessThanOrEqualTo()


void BitWeaver::complexLessThanOrEqualTo(uint64_t comparedTo,
  uint64_t* resultBitMap, uint64_t cardinality) {
  int numOfCodesInAWord = WORD_SIZE;
  int bytesInOneSegment = codeLength * WORD_SIZE / 8;
  uint64_t numberOfSegments = (cardinality / WORD_SIZE) ? 
    ((cardinality % WORD_SIZE) ? cardinality/WORD_SIZE + 1 : 
    cardinality/WORD_SIZE) : 1;

  //bit representation of comparator
  uint64_t comparedToMasks[codeLength];
  for (int i = 0; i < codeLength; i++) {
    uint64_t comparedToMask = 
      ((1 << (codeLength - i - 1)) & comparedTo) >> 
      (codeLength - i - 1);
    comparedToMasks[i] = (comparedToMask) ? -1 : 0;
  }

  int prefixDistance = (codeLength > 16) ? 15 : 7; 
  __m512i zeroBitMapVector = _mm512_set1_epi64(0);
  __mmask8 pruningCheckVector;


  for (uint64_t segCount = 0; segCount < numberOfSegments; 
    segCount += 8) {
    __m512i ltBitMapVector = _mm512_set1_epi64(0);
    __m512i eqBitMapVector = _mm512_load_epi64(resultBitMap + segCount);
    for (int j = 0; j < codeLength; j++) {
      #ifdef PREFETCH
        if (j < (codeLength - prefixDistance)) {
          PREFETCH_INTO_L1_AND_L2(codeVectors[j + prefixDistance] + segCount);
        } else {
          PREFETCH_INTO_L1_AND_L2(codeVectors[j - (codeLength - \
            prefixDistance)] + segCount + 8);
        }
      #endif
      __m512i avxVector = _mm512_load_epi64(codeVectors[j] + segCount);
      uint64_t attributeMask = comparedToMasks[j];
      __m512i attributeMaskVector = _mm512_set1_epi64(attributeMask);
      __m512i temp = _mm512_andnot_epi64(avxVector, attributeMaskVector);
      temp = _mm512_and_epi64(temp, eqBitMapVector);
      ltBitMapVector = _mm512_or_epi64(temp, ltBitMapVector);
      temp = _mm512_xor_epi64(avxVector, attributeMaskVector);
      eqBitMapVector = _mm512_andnot_epi64(temp, eqBitMapVector);


      // check for early pruning
      pruningCheckVector = _mm512_cmp_epu64_mask(eqBitMapVector,
        zeroBitMapVector, _MM_CMPINT_EQ);

      if (pruningCheckVector == -1) {
        break;
      }

    }
    ltBitMapVector = _mm512_or_epi64(ltBitMapVector, eqBitMapVector);
    _mm512_store_epi64(resultBitMap + segCount, ltBitMapVector);        
  }    

} // end complexLessThanOrEqualTo()


void BitWeaver::greaterThanOrEqualTo(uint64_t comparedTo, 
  uint64_t* resultBitMap, uint64_t cardinality) {
  int numOfCodesInAWord = WORD_SIZE;
  int bytesInOneSegment = codeLength * WORD_SIZE / 8;
  uint64_t numberOfSegments = (cardinality / WORD_SIZE) ? 
    ((cardinality % WORD_SIZE) ? cardinality/WORD_SIZE + 1 : 
    cardinality/WORD_SIZE) : 1;

  //bit representation of comparator
  uint64_t comparedToMasks[codeLength];
  for (int i = 0; i < codeLength; i++) {
    uint64_t comparedToMask = 
      ((1 << (codeLength - i - 1)) & comparedTo) >> 
      (codeLength - i - 1);
    comparedToMasks[i] = (comparedToMask) ? -1 : 0;
  }
  
  uint64_t segCount = 0;

  int prefixDistance = (codeLength > 16) ? 15 : 7; 
  __m512i zeroBitMapVector = _mm512_set1_epi64(0);
  __mmask8 pruningCheckVector;


  for (segCount = 0; segCount < numberOfSegments; 
    segCount += 8) {

    /*PREFETCH_INTO_L1D_FOR_WRITES(resultBitMap + segCount); */

    __m512i gtBitMapVector = _mm512_set1_epi64(0);
    __m512i eqBitMapVector = _mm512_set1_epi64((uint64_t)-1);


    // hardware prefetcher cannot prefetch across page boundaries

    for (int j = 0; j < codeLength; j++) {
      #ifdef PREFETCH
        if (j < (codeLength - prefixDistance)) {
          PREFETCH_INTO_L1_AND_L2(codeVectors[j + prefixDistance] + segCount);
        } else {
          PREFETCH_INTO_L1_AND_L2(codeVectors[j - (codeLength - \
            prefixDistance)] + segCount + 8);
        }
      #endif
 

      __m512i avxVector = _mm512_load_epi64(codeVectors[j] + segCount);

      uint64_t attributeMask = comparedToMasks[j];
      __m512i attributeMaskVector = _mm512_set1_epi64(attributeMask);

      __m512i temp0 = _mm512_andnot_epi64(attributeMaskVector, avxVector);
      temp0 = _mm512_and_epi64(temp0, eqBitMapVector);
      gtBitMapVector = _mm512_or_epi64(temp0, gtBitMapVector);
      __m512i temp1= _mm512_xor_epi64(avxVector, attributeMaskVector);
      eqBitMapVector = _mm512_andnot_epi64(temp1, eqBitMapVector);

      // check for early pruning
      pruningCheckVector = _mm512_cmp_epu64_mask(eqBitMapVector,
        zeroBitMapVector, _MM_CMPINT_EQ);

      if (pruningCheckVector == -1) {
        break;
      }

    } // end inner for loop
    gtBitMapVector = _mm512_or_epi64(gtBitMapVector, eqBitMapVector);
    _mm512_store_epi64(resultBitMap + segCount, gtBitMapVector);

  } // end outer for loop

} // end greaterThanOrEqualTo()


void BitWeaver::complexGreaterThanOrEqualTo(uint64_t comparedTo, 
  uint64_t* resultBitMap, uint64_t cardinality) {
  int numOfCodesInAWord = WORD_SIZE;
  int bytesInOneSegment = codeLength * WORD_SIZE / 8;
  uint64_t numberOfSegments = (cardinality / WORD_SIZE) ? 
    ((cardinality % WORD_SIZE) ? cardinality/WORD_SIZE + 1 : 
    cardinality/WORD_SIZE) : 1;

  //bit representation of comparator
  uint64_t comparedToMasks[codeLength];
  for (int i = 0; i < codeLength; i++) {
    uint64_t comparedToMask = 
      ((1 << (codeLength - i - 1)) & comparedTo) >> 
      (codeLength - i - 1);
    comparedToMasks[i] = (comparedToMask) ? -1 : 0;
  }
  
  uint64_t segCount = 0;

  int prefixDistance = (codeLength > 16) ? 15 : 7; 
  __m512i zeroBitMapVector = _mm512_set1_epi64(0);
  __mmask8 pruningCheckVector;


  for (segCount = 0; segCount < numberOfSegments; 
    segCount += 8) {

    /*PREFETCH_INTO_L1D_FOR_WRITES(resultBitMap + segCount); */

    __m512i gtBitMapVector = _mm512_set1_epi64(0);
    __m512i eqBitMapVector = _mm512_load_epi64(resultBitMap + segCount);


    // hardware prefetcher cannot prefetch across page boundaries

    for (int j = 0; j < codeLength; j++) {
      #ifdef PREFETCH
        if (j < (codeLength - prefixDistance)) {
          PREFETCH_INTO_L1_AND_L2(codeVectors[j + prefixDistance] + segCount);
        } else {
          PREFETCH_INTO_L1_AND_L2(codeVectors[j - (codeLength - \
            prefixDistance)] + segCount + 8);
        }
      #endif
 

      __m512i avxVector = _mm512_load_epi64(codeVectors[j] + segCount);

      uint64_t attributeMask = comparedToMasks[j];
      __m512i attributeMaskVector = _mm512_set1_epi64(attributeMask);

      __m512i temp0 = _mm512_andnot_epi64(attributeMaskVector, avxVector);
      temp0 = _mm512_and_epi64(temp0, eqBitMapVector);
      gtBitMapVector = _mm512_or_epi64(temp0, gtBitMapVector);
      __m512i temp1= _mm512_xor_epi64(avxVector, attributeMaskVector);
      eqBitMapVector = _mm512_andnot_epi64(temp1, eqBitMapVector);

      // check for early pruning
      pruningCheckVector = _mm512_cmp_epu64_mask(eqBitMapVector,
        zeroBitMapVector, _MM_CMPINT_EQ);

      if (pruningCheckVector == -1) {
        break;
      }

    } // end inner for loop
    gtBitMapVector = _mm512_or_epi64(gtBitMapVector, eqBitMapVector);
    _mm512_store_epi64(resultBitMap + segCount, gtBitMapVector);

  } // end outer for loop

} // end complexGreaterThanOrEqualTo()


void BitWeaver::between(uint64_t value1, uint64_t value2,
  uint64_t* resultBitMap, uint64_t cardinality) {
  int numOfCodesInAWord = WORD_SIZE;
  int bytesInOneSegment = codeLength * WORD_SIZE / 8;
  uint64_t numberOfSegments = (cardinality / WORD_SIZE) ? 
    ((cardinality % WORD_SIZE) ? cardinality/WORD_SIZE + 1 : 
    cardinality/WORD_SIZE) : 1;

  //bit representation of comparators
  uint64_t value1Map[codeLength];
  uint64_t value2Map[codeLength];
  for (int i = 0; i < codeLength; i++) {
    uint64_t bitAtValue1 = 
      ((1 << (codeLength - i - 1)) & value1) >> 
      (codeLength - i - 1);
    uint64_t bitAtValue2 = 
      ((1 << (codeLength - i - 1)) & value2) >> 
      (codeLength - i - 1);
    value1Map[i] = (bitAtValue1) ? -1 : 0;
    value2Map[i] = (bitAtValue2) ? -1 : 0; 
  }

  int prefixDistance = (codeLength > 16) ? 7 : 3; 
  __m512i zeroBitMapVector = _mm512_set1_epi64(0);
  __mmask8 pruningCheckVector1;
  __mmask8 pruningCheckVector2;

  for (uint64_t segCount = 0; segCount < numberOfSegments; 
    segCount += 8) {
    __m512i ltBitMapVector = _mm512_set1_epi64(0);
    __m512i gtBitMapVector = _mm512_set1_epi64(0);
    __m512i eq1BitMapVector = _mm512_set1_epi64((uint64_t)-1);
    __m512i eq2BitMapVector = _mm512_set1_epi64((uint64_t)-1);

    for (int j = 0; j < codeLength; j++) {
      
      // prefix lines if required
      #ifdef PREFETCH
        if (j < (codeLength - prefixDistance)) {
          PREFETCH_INTO_L1_AND_L2(codeVectors[j + prefixDistance] + segCount);
        } else {
          PREFETCH_INTO_L1_AND_L2(codeVectors[j - (codeLength - \
            prefixDistance)] + segCount + 8);
        }
      #endif
      __m512i avxVector = _mm512_load_epi64(codeVectors[j] + segCount);
      uint64_t attribute1Mask = value1Map[j];
      uint64_t attribute2Mask = value2Map[j];
      __m512i attribute1MaskVector = _mm512_set1_epi64(attribute1Mask);
      __m512i attribute2MaskVector = _mm512_set1_epi64(attribute2Mask);

      __m512i temp0 = _mm512_andnot_epi64(attribute1MaskVector, avxVector);
      temp0 = _mm512_and_epi64(temp0, eq1BitMapVector);
      gtBitMapVector = _mm512_or_epi64(temp0, gtBitMapVector);
      __m512i temp1 = _mm512_xor_epi64(avxVector, attribute1MaskVector);
      eq1BitMapVector = _mm512_andnot_epi64(temp1, eq1BitMapVector);

      __m512i temp2 = _mm512_andnot_epi64(avxVector, attribute2MaskVector);      
      temp2 = _mm512_and_epi64(temp2, eq2BitMapVector);
      ltBitMapVector = _mm512_or_epi64(temp2, ltBitMapVector);
      __m512i temp3 = _mm512_xor_epi64(avxVector, attribute2MaskVector);
      eq2BitMapVector = _mm512_andnot_epi64(temp3, eq2BitMapVector);


      // check for early pruning
      pruningCheckVector1 = _mm512_cmp_epu64_mask(eq1BitMapVector,
        zeroBitMapVector, _MM_CMPINT_EQ);
      pruningCheckVector2 = _mm512_cmp_epu64_mask(eq2BitMapVector,
        zeroBitMapVector, _MM_CMPINT_EQ);

      if ((pruningCheckVector1 == -1) && (pruningCheckVector2 == -1)) {
        break;
      }


    } // end inner for loop
    gtBitMapVector = _mm512_and_epi64(gtBitMapVector, ltBitMapVector);
    _mm512_store_epi64(resultBitMap + segCount, gtBitMapVector);        
  } // end outer for loop  

} // end between()

void BitWeaver::complexBetween(uint64_t value1, uint64_t value2,
  uint64_t* resultBitMap, uint64_t cardinality) {
  int numOfCodesInAWord = WORD_SIZE;
  int bytesInOneSegment = codeLength * WORD_SIZE / 8;
  uint64_t numberOfSegments = (cardinality / WORD_SIZE) ? 
    ((cardinality % WORD_SIZE) ? cardinality/WORD_SIZE + 1 : 
    cardinality/WORD_SIZE) : 1;

  //bit representation of comparators
  uint64_t value1Map[codeLength];
  uint64_t value2Map[codeLength];
  for (int i = 0; i < codeLength; i++) {
    uint64_t bitAtValue1 = 
      ((1 << (codeLength - i - 1)) & value1) >> 
      (codeLength - i - 1);
    uint64_t bitAtValue2 = 
      ((1 << (codeLength - i - 1)) & value2) >> 
      (codeLength - i - 1);
    value1Map[i] = (bitAtValue1) ? -1 : 0;
    value2Map[i] = (bitAtValue2) ? -1 : 0; 
  }

  int prefixDistance = (codeLength > 16) ? 7 : 3; 
  __m512i zeroBitMapVector = _mm512_set1_epi64(0);
  __mmask8 pruningCheckVector1;
  __mmask8 pruningCheckVector2;

  for (uint64_t segCount = 0; segCount < numberOfSegments; 
    segCount += 8) {
    __m512i ltBitMapVector = _mm512_set1_epi64(0);
    __m512i gtBitMapVector = _mm512_set1_epi64(0);
    __m512i eq1BitMapVector = _mm512_load_epi64(resultBitMap + segCount);
    __m512i eq2BitMapVector = _mm512_load_epi64(resultBitMap + segCount);

    for (int j = 0; j < codeLength; j++) {
      
      // prefix lines if required
      #ifdef PREFETCH
        if (j < (codeLength - prefixDistance)) {
          PREFETCH_INTO_L1_AND_L2(codeVectors[j + prefixDistance] + segCount);
        } else {
          PREFETCH_INTO_L1_AND_L2(codeVectors[j - (codeLength - \
            prefixDistance)] + segCount + 8);
        }
      #endif
      __m512i avxVector = _mm512_load_epi64(codeVectors[j] + segCount);
      uint64_t attribute1Mask = value1Map[j];
      uint64_t attribute2Mask = value2Map[j];
      __m512i attribute1MaskVector = _mm512_set1_epi64(attribute1Mask);
      __m512i attribute2MaskVector = _mm512_set1_epi64(attribute2Mask);

      __m512i temp0 = _mm512_andnot_epi64(attribute1MaskVector, avxVector);
      temp0 = _mm512_and_epi64(temp0, eq1BitMapVector);
      gtBitMapVector = _mm512_or_epi64(temp0, gtBitMapVector);
      __m512i temp1 = _mm512_xor_epi64(avxVector, attribute1MaskVector);
      eq1BitMapVector = _mm512_andnot_epi64(temp1, eq1BitMapVector);

      __m512i temp2 = _mm512_andnot_epi64(avxVector, attribute2MaskVector);      
      temp2 = _mm512_and_epi64(temp2, eq2BitMapVector);
      ltBitMapVector = _mm512_or_epi64(temp2, ltBitMapVector);
      __m512i temp3 = _mm512_xor_epi64(avxVector, attribute2MaskVector);
      eq2BitMapVector = _mm512_andnot_epi64(temp3, eq2BitMapVector);


      // check for early pruning
      pruningCheckVector1 = _mm512_cmp_epu64_mask(eq1BitMapVector,
        zeroBitMapVector, _MM_CMPINT_EQ);
      pruningCheckVector2 = _mm512_cmp_epu64_mask(eq2BitMapVector,
        zeroBitMapVector, _MM_CMPINT_EQ);

      if ((pruningCheckVector1 == -1) && (pruningCheckVector2 == -1)) {
        break;
      }


    } // end inner for loop
    gtBitMapVector = _mm512_and_epi64(gtBitMapVector, ltBitMapVector);
    _mm512_store_epi64(resultBitMap + segCount, gtBitMapVector);        
  } // end outer for loop  

} // end complexBetween()