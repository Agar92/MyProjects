#pragma once
#ifndef T2NSGANGULAR_HH
#define T2NSGANGULAR_HH

#include "globals.hh"

namespace T2NSGangular
{
  G4double RandomizeCost(G4int tgZ, G4int tgA, G4int sPDG, G4int rID,
                         G4int level, G4double Einc, G4int incZA = 1);
}

#endif
