#pragma once
#ifndef T2NSGANGULAR_RWNODE_HH
#define T2NSGANGULAR_RWNODE_HH

//#include "G4LorentzVector.hh"
#include "globals.hh"
#include <fstream>
#include <vector>
#include <array>

#include "T3Types.h"

class T2NSGangular_RWnode
{
public:
  static const size_t _num_point = 127;
private:
  static std::array<G4double, _num_point> isotropic();
public:
  T2NSGangular_RWnode(G4double e = 0,
                      const std::array<G4double, _num_point> src = isotropic(),
                      G4double newpr = 0, G4double newsl = 1e5);
  T2NSGangular_RWnode(T2NSGangular_RWnode const &node);
  G4double Get_E() const {return E;}
  G4double Get_pr() const {return pr;}
  G4double Get_sl() const {return sl;}
  const std::array<G4double, _num_point>& Get_V() const {return V;}
  G4double Get_V(size_t point_num) const;
  T2NSGangular_RWnode& operator=(T2NSGangular_RWnode const &node);
  std::ofstream& save_binary( std::ofstream& out_stream) const;
  std::ifstream& load_binary( std::ifstream& in_stream);
  void Set_E(G4double e) {E = e;}
  void Set_pr(G4double newpr) {pr = newpr;}
  void Set_sl(G4double newsl) {sl = newsl;}
  void Set_V(size_t point_num, G4double val);
  void Set_V(const std::array<G4double, _num_point> src);
  // from T2NElasticCrossSection
  // TODO
  G4bool is_isotropic() const;
private:
  std::array<G4double, _num_point> V; /// integrated from -1 cos(theta) distribution
  G4double E;              /// incident neutron energy
  G4double sl;             /// forward exponent slope
  G4double pr;             /// forward exponent constribution
};

// TODO overload >> and both for binary streams
std::ostream& operator<<(std::ostream& os, const T2NSGangular_RWnode& inst);

#endif
