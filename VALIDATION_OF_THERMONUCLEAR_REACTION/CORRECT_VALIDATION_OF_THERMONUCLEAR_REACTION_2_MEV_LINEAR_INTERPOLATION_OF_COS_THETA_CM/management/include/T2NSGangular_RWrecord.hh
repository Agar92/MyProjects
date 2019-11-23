#pragma once
#ifndef T2NSGANGULAR_RWRECORD_HH
#define T2NSGANGULAR_RWRECORD_HH

#include <vector>
#include "T2NSGangular_RWnode.hh"

class T2NSGangular_RWrecord: public std::vector<T2NSGangular_RWnode>
{
public:
  G4int Get_levN() const {return levN;}
  void Set_levN(G4int val) {levN = val;}
  std::ofstream& save_binary( std::ofstream& out_stream) const;
  std::ifstream& load_binary( std::ifstream& in_stream);
  T2NSGangular_RWnode interpolate(G4double E) const;
private:
  G4int levN;       /// excited level number as in T2NuclearLevels
};

std::ostream& operator<<(std::ostream& os, const T2NSGangular_RWrecord& inst);

#endif
