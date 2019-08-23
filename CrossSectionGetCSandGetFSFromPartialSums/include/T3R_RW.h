#pragma once
#ifndef T3R_RW_H
#define T3R_RW_H

#include <string>
#include <ostream>
#include "T3R_node.h"
#include "T3TabulatedFunction_RW.h"
#include "T3TabulatedCS_RW.h"

namespace t3 {

class T3R_RW
{
public:
  //constructor:
  T3R_RW(const T3TabulatedCS& tcs = T3TabulatedCS(), T3int _tgZ = 0, T3int _tgA = 0);
  //copy constructor
  //without it aT3R_DDCS.GetT3R_RW(); in constructor of T3ElasticIonIonCSimpl.h does not work.
  T3R_RW(const T3R_RW & rhs)
  {
    ND.resize(rhs.size());
    for(int i=0; i<size(); ++i)
    {
      ND.at(i) = new T3R_node();
      *ND.at(i) = (*rhs.ND.at(i));
    }
    tgZ=rhs.tgZ;
    tgA=rhs.tgA;
    T3data=rhs.T3data;
  }
  //destructor:
  ~T3R_RW();
  //WRITE TO FILE
  void save_binary(T3int targZ, T3int targA, T3int pZ, T3int pA,
                   const T3String& suffix="unknown") const;

  void save_binary( T3int targZ, T3int targA, const T3String& sID,
                    const T3String& suffix = "unknown") const;

  void save_binary( const T3String& fname ) const; // @@ Save binary table to file (?)
  //READ FROM FILE
  T3bool load_binary(T3int targZ, T3int targA,
                     T3int pZ, T3int pA, const T3String& suffix="any");

  T3bool load_binary(T3int targZ, T3int targA, const T3String& rID,
                     const T3String& suffix = "any");
  //
  T3String default_file( T3int targZ, T3int targA, T3int pZ, T3int pA,
                         const T3String& suffix) const;

  T3String default_file( T3int targZ, T3int targA, const T3String& sID,
                         const T3String& suffix) const;
  //here T3TabulatedCS table (x,y) is divided in T2R_node's (xi,yi) and they are pushed
  //in ND
  void push_back(T3R_node* N) {ND.push_back(N);}    // add in the end a new R_Node
  size_t size() const       {return ND.size();}     // Get size of the NHG-nodes vector
  const T3R_node* NDI(T3int i) const {return ND.at(i);} // Get i-th NHG-node Pointer
  T3double GetZ() const  {return tgZ;}              // Extract the A-parameter for CS=A+B/p
  T3double GetA() const  {return tgA;}              // Extract the B-parameter for CS=A+B/p
  void SetZ(T3double Z)  {tgZ = Z;}                 // Put the A-parameter for CS=A+B/p
  void SetA(T3double A)  {tgA = A;}                 // Put the B-parameter for CS=A+B/p
  inline T3R_RW& operator=(const T3R_RW& rhs);
  
  T3bool load_binary( const T3String& fname );      // Read binary table from file
private:
  // Body
  std::vector<T3R_node*> ND;                      // Vector of pointers to (n,h'g) nodes
  T3double tgZ;     // FIXME currently not saved in files (it is in the name of file)
  T3double tgA;     // maybe add to the end for compability ? @@
  T3String T3data;  // path to the data directory
};

std::ostream& operator<<(std::ostream& os, const T3R_RW& inst);

inline T3R_RW& T3R_RW::operator=(const T3R_RW & rhs)
{
//there was ERROR:
//ND = rhs.ND;
  ND.resize(rhs.size());
  for(int i=0; i<size(); ++i)
  {
    ND.at(i) = new T3R_node();
    *ND.at(i) = (*rhs.ND.at(i));
  }
  tgZ=rhs.tgZ;
  tgA=rhs.tgA;
  T3data=rhs.T3data;
  return *this;
}
} // namespace t3
#endif // T3R_RW_H
