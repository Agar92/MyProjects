#pragma once
#ifndef T3NSGANGULAR_RW_H
#define T3NSGANGULAR_RW_H
#include <array>
#include "T3NSGangular_RWrecord.h"

namespace t3 {
class T3NSGangular_RW: public std::vector<T3NSGangular_RWrecord>
{
public:
  T3int Get_secPDG() const {return secPDG;}
  void Set_secPDG(T3int val) {secPDG = val;}
  T3int Get_rID() const {return rID;}
  void Set_rID(T3int val) {rID = val;}
  T3int Get_Z() const {return Z;}
  void Set_Z(T3int val) {Z = val;}
  T3int Get_A() const {return A;}
  void Set_A(T3int val) {A = val;}

  static auto constexpr Bin2=NumberOfBinsInDDApproximation2;
  static auto constexpr Bin3=NumberOfBinsInDDApproximation3;
  void save_binary(T3int tgZ, T3int tgA, const T3String& rid, T3int pdg = 2112,
                   T3int incZA = 1) const;
  T3bool load_binary(T3int tgZ, T3int tgA, const T3String& rid, T3int pdg = 2112,
                     T3int incZA = 1);
  std::vector<T3NSGangular_RWrecord> Get_vector_of_T3NSGangular_RWrecords()
  {
      return *this;
  }
private:
  T3int Z;//Z of inc particle
  T3int A;//A of inc particle
  T3int secPDG;/// secondary particle PDGcode
  T3int rID;   /// reaction identification
//  /////////////////
  std::ofstream& save_binary( std::ofstream& out_stream) const;
  std::ifstream& load_binary( std::ifstream& in_stream);
  void save_binary( const T3String& fname ) const;// Save Binary table to the file
  T3bool load_binary( const T3String& fname );    // Read binary table from file
  T3String default_file( T3int tgZ, T3int tgA,
                         const T3String& rid, T3int pdg,
                         T3int incZA = 1) const;  // Name definition
};

std::ostream& operator<<(std::ostream& os, const T3NSGangular_RW& inst);
} // namespace t3
#endif
