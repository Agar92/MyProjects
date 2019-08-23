#pragma once
#ifndef T3R_node_H
#define T3R_node_H

#include <vector>
#include <fstream>
#include <iostream>
#include "globals.h"

namespace t3 {
//THIS CLASS
//the object (E,CS)
class T3R_node
{
public:
  T3R_node(double E = 0., double XS = 0.);
  void save_binary( std::ofstream* out) const;
  void load_binary( std::ifstream* in); // @@ tmp for old format
  bool operator==(const T3R_node rval); // used once
  // Selectors
  double E() const {return e;}            // Get the projectile neutron energy
  double XS() const {return xs;}          // Get theSummedUp (h,any) reactionCrossSection
  // Modifiers
  void SetE( double En )  {e = En;}       // Set the projectile neutron energy
  void SetXS( double CS ) {xs = CS;}      // Set the (n,n'g) reaction cross-section
private:
  //body
  double e;                // neutron energy
  double xs;               // Summed cross-section of all levels + continuum
};

std::ostream& operator<<(std::ostream& os, const T3R_node& inst);
} // namespace t3
#endif // T3R_node_H
