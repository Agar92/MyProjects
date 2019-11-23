#ifndef GLOBALS_H
#define GLOBALS_H

//Include base types
#include "T3Types.h"

#ifndef FALSE
  #define FALSE 0
#endif
#ifndef TRUE
  #define TRUE 1
#endif

// /////////////////////////////////////////////////////// //
//Here the number of bins for d-d approximation and for d-d//
//T3ElasticIonIonCSImpl.h and for T3ElasticIonIonFSImpl.h  //
//and for writing and reading files (E,CS) and             //
//(|t|, partial sums) is defined                           //
// /////////////////////////////////////////////////////// //
//number of bins for (E,CS) db:
constexpr int NumberOfBinsInDDApproximation1 = 32;
//YET NumberOfBinsInDDApproximation1 = NumberOfBinsInDDApproximation2. Is is simple.
//number of bins for E in partial sums db:
constexpr int NumberOfBinsInDDApproximation2 = 128;//256;//256;
//number of bins for 1-cos(theta_cm) in partial sums db:
constexpr int NumberOfBinsInDDApproximation3 = 128;//256;//512;
// /////////////////////////////////////////////////////// //

#include <algorithm>  // Retrieve definitions of min/max

#endif /* GLOBALS_H */
