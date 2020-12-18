/*******************************************************************************
*                                                                              *
* Vertical Bitweaving format converter - creates bit slices for a column       *
*                                                                              *
* This program generates bit slices corresponding to each bit, so that they can*
* be used for BitWeaving. It's assumed that the input column is already coded. *
* Padding is inserted to ensure that each bitslice's offset starts at 4KB, so  *
* that it can be MMAPed.                                                       *
*                                                                              *
* Usage format: create_bitw_files CODED-COLUMN-FILENAME CARDINALITY CODE-LENGTH*
*                                                                              *
* The current version of this code assumes correct input                       *
*                                                                              *
* Author: Sanchit Jain <sanchit@cs.wisc.edu>                                   *
*                                                                              *
*******************************************************************************/

#include <cstdint>
#include <stdlib.h>
#include <cstdio>
#include <iostream>
#include <vector>
#include <string>
#include <unistd.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#define WORD_SIZE 64
#define BYTES_IN_WORD 8
#define PAGE_SIZE 4096


int main(int argc, char* argv[]) {
  char* inputFileName = argv[1];
  uint64_t count = (uint64_t)std::stoi(argv[2]);
  int codeLength = std::stoi(argv[3]);
  // these are the bit-slices
  uint64_t** words = new uint64_t*[codeLength];

  // number of CPU words required to store a bit slice
  uint64_t num_of_cpu_words =  (count / WORD_SIZE) ? 
    ((count % WORD_SIZE) ? count/WORD_SIZE + 1 : 
    count/WORD_SIZE) : 1;
  // we need to make it a multiple of page size, in order to mmap for reading
  int numOfFillerBytes = ((num_of_cpu_words * BYTES_IN_WORD) % PAGE_SIZE) ? \
    (PAGE_SIZE - ((num_of_cpu_words * BYTES_IN_WORD) % PAGE_SIZE)) : 0;

  // these allocations don't do anything special in this file
  int i = 0;
  int j = ((codeLength - 8) > 0) ? 8 : codeLength;
  for (i = 0; i < j; i++) {
    words[i] = (uint64_t *)(aligned_alloc(PAGE_SIZE * 16,
    num_of_cpu_words * 8));
  }
  if ((i = codeLength - i) > 0) {
    j = ((codeLength - 16) > 0) ? 16 : codeLength;
    for (i = 8; i < j; i++) {
      words[i] = 
        (uint64_t *)(aligned_alloc(PAGE_SIZE * 16,
        num_of_cpu_words * 8 + 1024));
      words[i] = words[i] + 128;
    }
    if ((i = codeLength - i) > 0) {
      j = ((codeLength - 24) > 0) ? 24 : codeLength;
      for (i = 16; i < j; i++) {
        words[i] = 
          (uint64_t *)(aligned_alloc(PAGE_SIZE * 16,
          num_of_cpu_words * 8 + 2048));
        words[i] = words[i] + 256;
      }
      if ((i = codeLength - i) > 0) {
        for (i = 24; i < codeLength; i++) {
          words[i] = 
            (uint64_t *)(aligned_alloc(PAGE_SIZE * 16,
            num_of_cpu_words * 8 + 3072));
          words[i] = words[i] + 384;
        }          
      }          
    }      
  }

  int inputFileFD = open(inputFileName, O_RDONLY);
  int* numbers = (int*)(mmap(NULL, count * sizeof(int), PROT_READ,
    MAP_PRIVATE, inputFileFD, 0));
    std::string fileName = std::string(inputFileName) + "_bitslices.bwv";
    int outputFileFD = open(fileName.c_str(), O_CREAT|O_RDWR,
        S_IRWXU|S_IROTH);

    for (int bitNumber = 0; bitNumber < codeLength; bitNumber++) {
      // bit slice file for this bit number
     uint64_t lastSegmentNumber = num_of_cpu_words - 1;
     int numOfCodesToRead = WORD_SIZE;
     for (uint64_t segNo = 0; segNo < num_of_cpu_words; segNo++) {
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
           & numbers[segNo * WORD_SIZE + i]) >> 
           (codeLength - 1 - bitNumber)) << (WORD_SIZE - i - 1));
       }
       words[bitNumber][segNo] = currentWord;
       write(outputFileFD, &(words[bitNumber][segNo]), BYTES_IN_WORD);
     }
     // write padding bytes to fill the last page of each bit-slice
     write(outputFileFD, words[bitNumber], numOfFillerBytes);
   }
   // Save the bit slice 
   close(outputFileFD);
   munmap(numbers, count * sizeof(int));
   return 0;
} // end main()