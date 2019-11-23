#pragma once
#ifndef T3NSGANGULAR_NODE_HH
#define T3NSGANGULAR_NODE_HH

#include "T3NSGangular_RWnode.hh"


//\\//************************//
#include "T3RNG.h"
//for memcpy:
#include <cstring>
//\\//************************//


class T3NSGangular_node
{
public:
  T3NSGangular_node(const T3NSGangular_RWnode& rwnode = T3NSGangular_RWnode());
  static const size_t _num_point = T3NSGangular_RWnode::_num_point;
  T3NSGangular_node& operator=(T3NSGangular_node const &node);
  //\\//G4double RandomizeCost();
  T3double RandomizeCost(RNDGenerator & generator);
  T3double Get_V(size_t point_num) const {return V[point_num];}
  T3double Get_a(size_t point_num) const {return a[point_num];}
  T3double Get_b(size_t point_num) const {return b[point_num];}
  T3double Get_c(size_t point_num) const {return c[point_num];}
  T3double Get_pr() const {return pr;}
  T3double Get_sl() const {return sl;}
private:
  T3double V[_num_point];  /// integrated from -1 cos(theta) distribution
  T3double a[_num_point];  /// interpolation coefficients
  T3double b[_num_point];  /// interpolation coefficients
  T3double c[_num_point];  /// interpolation coefficients
  T3double pr;       /// forward contribution
  T3double sl;       /// forward slope
  T3NSGangular_node(T3NSGangular_node const &node);
};

std::ostream& operator<<(std::ostream& os, const T3NSGangular_node& inst);

#endif
