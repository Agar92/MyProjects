#pragma once
#ifndef T3NSGANGULAR_RWNODE_H
#define T3NSGANGULAR_RWNODE_H

#include "globals.h"
#include <fstream>
#include <vector>
#include <array>

namespace t3 {
class T3NSGangular_RWnode
{
public:
// ////////
//NumberOfBinsInDDApproximation is defined in T3Types.h
// ////////
  static const size_t _num_point = NumberOfBinsInDDApproximation3;
private:
  static std::array<T3double, _num_point+1> isotropic();
public:
  T3NSGangular_RWnode(T3double e = 0,
                      const std::array<T3double, _num_point+1> src = isotropic());
  
  T3NSGangular_RWnode(T3NSGangular_RWnode const &node);
  T3double Get_E() const {return E;}
  const std::array<T3double, _num_point+1>& Get_V() const {return V;}  
  T3double Get_V(size_t point_num) const;
  T3NSGangular_RWnode& operator=(T3NSGangular_RWnode const &node);
  std::ofstream& save_binary( std::ofstream& out_stream) const;
  std::ifstream& load_binary( std::ifstream& in_stream);
  void Set_E(T3double e) {E = e;}
  void Set_V(size_t point_num, T3double val);
  void Set_V(const std::array<T3double, _num_point+1> src);
  // TODO
  T3bool is_isotropic() const;
private:
  //partial sum SUM_0^k(dsigma/d|t|*d|t|):
  std::array<T3double, _num_point+1> V; /// integrated from -1 cos(theta) distribution
  T3double E;              /// incident neutron energy
};

// TODO overload >> and both for binary streams
std::ostream& operator<<(std::ostream& os, const T3NSGangular_RWnode& inst);
} // namespace t3
#endif
