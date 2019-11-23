#pragma once
#ifndef T3INELASTICDD_RANDOMIZECOST_HH
#define T3INELASTICDD_RANDOMIZECOST_HH

//necessary for u_quad_t:
#include <iostream>
#include "T3RNG.h"
#include "T3NSGangular_RWnode.hh"

class T3Inelasticdd_DB
{
public:
  const size_t _num_point = 127;
  const size_t _num_nodes = 512;
  //\\//G4double RandomizeCost();
  T3double RandomizeCost_node(RNDGenerator & generator, size_t node_num) const;
  T3double RandomizeCost(RNDGenerator & generator, T3double E) const;
  void Set_V(size_t node_num, size_t point_num, T3double val) {V[node_num][point_num]=val;}
  void Set_a(size_t node_num, size_t point_num, T3double aval) {a[node_num][point_num]=aval;}
  void Set_b(size_t node_num, size_t point_num, T3double bval) {b[node_num][point_num]=bval;}
  void Set_c(size_t node_num, size_t point_num, T3double cval) {c[node_num][point_num]=cval;}
  void Set_pr(size_t node_num, T3double prval) {pr[node_num]=prval;}
  void Set_sl(size_t node_num, T3double slval) {sl[node_num]=slval;}
  void Set_Einc(size_t node_num, T3double Eincval) {Einc[node_num]=Eincval;}
  T3double Get_V(size_t node_num, size_t point_num) const {return V[node_num][point_num];}
  T3double Get_a(size_t node_num, size_t point_num) const {return a[node_num][point_num];}
  T3double Get_b(size_t node_num, size_t point_num) const {return b[node_num][point_num];}
  T3double Get_c(size_t node_num, size_t point_num) const {return c[node_num][point_num];}
  T3double Get_pr(size_t node_num) const {return pr[node_num];}
  T3double Get_sl(size_t node_num) const {return sl[node_num];}
  T3double Get_Einc(size_t node_num) const {return Einc[node_num];}
private:
//There are _num_nodes nodes:  
//begin of node  
  T3double V[512][127];///integrated from -1 cos(theta) distribution
  T3double a[512][127];///interpolation coefficients
  T3double b[512][127];///interpolation coefficients
  T3double c[512][127];///interpolation coefficients
  T3double pr[512];           ///forward contribution
  T3double sl[512];           ///forward slope
  T3double Einc[512];
//end of node.  
};

#endif
