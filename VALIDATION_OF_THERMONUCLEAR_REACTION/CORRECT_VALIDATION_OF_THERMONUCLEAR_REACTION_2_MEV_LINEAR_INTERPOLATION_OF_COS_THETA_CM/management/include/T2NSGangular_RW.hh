#pragma once
#ifndef T2NSGANGULAR_RW_HH
#define T2NSGANGULAR_RW_HH

#include "T2NSGangular_RWrecord.hh"
#include <array>

class T2NSGangular_RW: public std::vector<T2NSGangular_RWrecord>
{
  //There is only 1 element in  T2NSGangular_RW
  //(std::vector<T2NSGangular_RWrecord>), because all (Ei, partial sums_i)
  //data sets are the members of T2NSGangular_RWrecord class.
  //In T2NSGangular_RWrecord there is G4int levN; private field.
  //It is the level of the nuclear excitation.
  //YET IN DATABASE THERE IS ONLY 1 LEVEL OF EXCITATION (0),
  //SO, IN T2NSGangular_RW THERE IS ONLY 1 RECORD AND THE SIZE OF
  //T2NSGangular_RW is 1.
  //So, to get separately (Ei, partial sums_i) data sets, one should get
  //nodes of T2NSGangular_RWrecord in T2NSGangular_RWrecord.!!!

public:
  G4int Get_secPDG() const {return secPDG;}
  void Set_secPDG(G4int val) {secPDG = val;}
  G4int Get_rID() const {return rID;}
  void Set_rID(G4int val) {rID = val;}
  G4int Get_Z() const {return Z;}
  void Set_Z(G4int val) {Z = val;}
  G4int Get_A() const {return A;}
  void Set_A(G4int val) {A = val;}
  void save_binary(G4int tgZ, G4int tgA, const G4String& rid, G4int pdg = 2112,
                   G4int incZA = 1) const;
  G4bool load_binary(G4int tgZ, G4int tgA, const G4String& rid, G4int pdg = 2112,
                     G4int incZA = 1);  
private:
  G4int Z;
  G4int A;
  G4int secPDG;     /// secondary particle PDGcode
  G4int rID; /// reaction identification
  std::ofstream& save_binary( std::ofstream& out_stream) const;
  std::ifstream& load_binary( std::ifstream& in_stream);
  void save_binary( const G4String& fname ) const;      // Save Binary Tables to the file
  G4bool load_binary( const G4String& fname );            // Read binary table from file
  G4String default_file( G4int tgZ, G4int tgA,
                         const G4String& rid, G4int pdg,
                         G4int incZA = 1) const;     // Name definition
//std::vector<size_t> index; // map level numbers to record numbers in the vector  
};

std::ostream& operator<<(std::ostream& os, const T2NSGangular_RW& inst);

#endif
