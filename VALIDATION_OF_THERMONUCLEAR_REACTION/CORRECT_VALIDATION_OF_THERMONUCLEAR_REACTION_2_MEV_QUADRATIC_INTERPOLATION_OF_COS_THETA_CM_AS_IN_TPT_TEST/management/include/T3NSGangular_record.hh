#pragma once
#ifndef T3NSGANGULAR_RECORD_HH
#define T3NSGANGULAR_RECORD_HH

#include <vector>
#include "T3NSGangular_node.hh"
#include "T3NSGangular_RW.hh"

//\\//****************************************************8//
#include "T3RNG.h"
//\\//****************************************************8//

class T3NSGangular_record
{
public:
  T3NSGangular_record(const T3NSGangular_RWrecord& rwrec);
  static const size_t _num_nodes = 512; // FIXME too big: memory consumption, find better
  //\\//G4double RandomizeCost(G4double E);
//\\//*************************************//
  T3double RandomizeCost(RNDGenerator & generator, T3double E);
//\\//*************************************//  

  
private:
  T3NSGangular_node nodes[_num_nodes];
  T3double           Einc[_num_nodes];
//   T2NSGangular_node interpolate(E);
  T3NSGangular_record(const T3NSGangular_record& rec);
  T3NSGangular_record& operator=(const T3NSGangular_record& rhs);
};

std::ostream& operator<<(std::ostream& os, const T3NSGangular_record& inst);

#endif
