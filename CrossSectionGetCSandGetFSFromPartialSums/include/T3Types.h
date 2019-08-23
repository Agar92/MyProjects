#pragma once
#ifndef T3TYPES_H
#define T3TYPES_H

#include <complex>
#include <iostream>
#include <string>
#include "T3Defs.h"

//********************************************************************************//
//This is a file for defining T3 types (T3int, T3double) and T3cout, T3cin, T3endl//
//********************************************************************************//

namespace t3 {

// Typedefs to decouple from library classes
// Typedefs for numeric types
//
typedef std::string T3String;
typedef double T3double;
typedef float T3float;
typedef int T3int;
typedef bool T3bool;
typedef long T3long;
typedef std::complex<T3double> T3complex;

#define T3cin  std::cin
#define T3cout std::cout
#define T3endl std::endl

} // namespace t3
#endif /* T3TYPES_H */
