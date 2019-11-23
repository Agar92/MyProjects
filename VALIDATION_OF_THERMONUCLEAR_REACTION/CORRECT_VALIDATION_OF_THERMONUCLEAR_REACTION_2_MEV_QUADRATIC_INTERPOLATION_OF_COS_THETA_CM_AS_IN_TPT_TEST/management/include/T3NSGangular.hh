#pragma once
#ifndef T3NSGANGULAR_HH
#define T3NSGANGULAR_HH

//\\//#include "globals.hh"

//\\//***********************************************//
#include "T3Types.h"
#include "T3RNG.h"
using namespace t3;
//\\//***********************************************//

namespace T3NSGangular
{
  //\\//G4double RandomizeCost(G4int tgZ, G4int tgA, G4int sPDG, G4int rID,
  //\\//                       G4int level, G4double Einc, G4int incZA = 1);
  T3double RandomizeCost(RNDGenerator & generator, T3int tgZ, T3int tgA, T3int sPDG, T3int rID,
                         T3int level, T3double Einc, T3int incZA = 1);
}

#endif
