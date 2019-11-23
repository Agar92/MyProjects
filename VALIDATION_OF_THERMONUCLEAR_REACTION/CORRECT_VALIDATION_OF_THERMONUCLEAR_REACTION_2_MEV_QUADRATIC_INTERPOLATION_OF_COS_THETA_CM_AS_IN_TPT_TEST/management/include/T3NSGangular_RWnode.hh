#pragma once
#ifndef T3NSGANGULAR_RWNODE_HH
#define T3NSGANGULAR_RWNODE_HH

//\\//#include "G4LorentzVector.hh"
//\\//#include "globals.hh"
#include "T3Globals.h"

#include <fstream>
#include <vector>
#include <array>


//\\//************************************//
//Added:
#include "T3Types.h"
#include "T3LorentzVector.h"
using namespace t3;
//\\//************************************//


class T3NSGangular_RWnode
{
public:
  static const size_t _num_point = 127;
private:
  static std::array<T3double, _num_point> isotropic();
public:
  T3NSGangular_RWnode(T3double e = 0,
                      const std::array<T3double, _num_point> src = isotropic(),
                      T3double newpr = 0, T3double newsl = 1e5);
  T3NSGangular_RWnode(T3NSGangular_RWnode const &node);
  T3double Get_E() const {return E;}
  T3double Get_pr() const {return pr;}
  T3double Get_sl() const {return sl;}
  const std::array<T3double, _num_point>& Get_V() const {return V;}
  T3double Get_V(size_t point_num) const;
  T3NSGangular_RWnode& operator=(T3NSGangular_RWnode const &node);
  std::ofstream& save_binary( std::ofstream& out_stream) const;
  std::ifstream& load_binary( std::ifstream& in_stream);
  void Set_E(T3double e) {E = e;}
  void Set_pr(T3double newpr) {pr = newpr;}
  void Set_sl(T3double newsl) {sl = newsl;}
  void Set_V(size_t point_num, T3double val);
  void Set_V(const std::array<T3double, _num_point> src);
  // from T2NElasticCrossSection
  //\\//std::vector<G4LorentzVector> CalcXYZ() const;
  std::vector<LorentzVector<FloatingType>> CalcXYZ() const;
  // TODO
  T3bool is_isotropic() const;
private:
  std::array<T3double, _num_point> V; /// integrated from -1 cos(theta) distribution
  T3double E;              /// incident neutron energy
  T3double sl;             /// forward exponent slope
  T3double pr;             /// forward exponent constribution

  //\\//G4LorentzVector CalcXYZ(size_t point_num) const;
  LorentzVector<FloatingType> CalcXYZ(size_t point_num) const;
};

// TODO overload >> and both for binary streams
std::ostream& operator<<(std::ostream& os, const T3NSGangular_RWnode& inst);

#endif
