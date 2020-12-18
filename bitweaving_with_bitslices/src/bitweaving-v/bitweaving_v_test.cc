#include "bitweaving_includes.h"
using namespace std;
using namespace std::chrono;

int main(int argc, char** argv) {
  // used wc to find number of lines excluding the header
  uint64_t cardinality = 59986052;

  uint64_t numOfResultVectors = (cardinality > WORD_SIZE) ? 
    (cardinality/WORD_SIZE + 1) : 1;

  __attribute((aligned(WORD_SIZE))) uint64_t* resultVectors = 
    (uint64_t *)(aligned_alloc(PAGE_SIZE, 
    sizeof(long long int) * numOfResultVectors));


  /* 
  // uncomment the multi-line comment measure for perf
  char proceed = 'y';
  printf("\n\nAre you ready to open perf in another terminal?\n"
    "The process ID is %d.\nEnter y\n", getpid());
  scanf("%c", &proceed); 
  */

  // create BitWeaver objects
  BitWeaver discountObject{"discount.col_bitslices.bwv", cardinality, 4};
  BitWeaver shipdateObject{"shipdate.col_bitslices.bwv", cardinality, 12};
  BitWeaver quantityObject{"quantity.col_bitslices.bwv", cardinality, 6};

  uint32_t placeHolder = 0;
  auto startTime = high_resolution_clock::now();
  uint64_t start = rdtscp(placeHolder);

  // test conditions on BitWeaver objects
  discountObject.between(4, 8, resultVectors, cardinality);
  // methods with prefix 'complex' use results from previous queries
  shipdateObject.complexBetween(729, 1095, resultVectors, cardinality);
  quantityObject.complexLessThan(24, resultVectors, cardinality);

  uint64_t end = rdtscp(placeHolder);
  auto stopTime = high_resolution_clock::now();
  uint64_t diff = (uint64_t)(end - start);
  cout << "\n\nIn AVX-512 based BitWeaving/V"
    " The number of cycles are " << diff << " and the number of cycles "
    "per code are " << diff/(cardinality * 1.0) << endl;
  auto duration = duration_cast<microseconds>(stopTime - startTime); 
    cout << "\n\nTime taken by AVX-512 BitWeaving/V is " <<
    duration.count() << " microseconds\n\n" << endl;

    /*
    // uncomment the multi-line comment to print result bitmap 
    for (int i = 0; i < numOfResultVectors; i++) {
      decimalToBinary(resultVectors[i]);
    } 
    */


  return 0;
} // end main()