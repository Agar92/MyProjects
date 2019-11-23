#ifndef T3TYPES_H
#define T3TYPES_H

#include <complex>
#include <iostream>
#include <string>
#include "T3Defs.h"

//********************************************************************************//
//This is a file for defining T3 types (T3int, T3double) and T3cout, T3cin, T3endl//
//********************************************************************************//

// Typedefs to decouple from library classes
// Typedefs for numeric types
//
typedef std::string G4String;
typedef double G4double;
typedef float G4float;
typedef int G4int;
typedef bool G4bool;
typedef long G4long;
typedef std::complex<G4double> G4complex;

#define G4cin  std::cin
#define G4cout std::cout
#define G4endl std::endl

// Forward declation of void type argument for usage in direct object
// persistency to define fake default constructors
//
//class __void__;

#endif /* T3TYPES_H */
