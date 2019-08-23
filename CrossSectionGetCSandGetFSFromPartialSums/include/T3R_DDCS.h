#pragma once
// ---------------------------------------------------------------------------
//  In this file there is a d-d elastic scattering approximation
//  from /users/ALPHA_ELAST/FINAL_CORRECT_dd .
// ----------------------------------------------------------------------
#ifndef T3R_DDCS_H
#define T3R_DDCS_H

//********************************************************************************************************************//
//This is a file for writing and reading |t| and partial sums for each energy for d-d elastic scattering approximation//
//********************************************************************************************************************//
#include <string>
#include <iostream>
#include <fstream>
#include <ostream>
#include <array>
#include <vector>
#include "T3Types.h"
#include "T3R_RW.h"
#include "T3TabulatedFunction_RW.h"
#include "T3TabulatedCS_RW.h"
#include "T3NSGangular_RWnode.h"
#include "T3NSGangular_RWrecord.h"
#include "T3NSGangular_RW.h"

namespace t3 {
using namespace units;
class T3R_DDCS
{
public:
  //constructor:
  T3R_DDCS();
  //destructor:
  ~T3R_DDCS(){std::cout<<"~T3R_DDCS(): "<<std::endl;}
  //this function fills all the arrays
  //fills (E  partial sum) table
  void Fill();//it does all the work
  //returns the value of d-d differential cross-section (in inner units) from approximation
  double GetdSigmadt(double Tls/*kinetic energy of the incident deuteron in MeV*/,
                   double t/*|t| - Mandelstahm variable in inner units*/) const;
  double GetdSigmadt_DD_Rutherford(double Tls/*kinetic energy of the incident deuteron in MeV*/,
                   double t/*|t| - Mandelstahm variable in inner units*/) const;
private:
  //number of bins in E of (E, partial sums):
  static constexpr int Bin2=NumberOfBinsInDDApproximation2;
  //number of bins in partial sums of (E, partial sums):
  static constexpr int Bin3=NumberOfBinsInDDApproximation3;
  //number of subbins, which we use for more accurate  integration of dsigma/d|t| over the Bin3 bins.
  static constexpr int SubBin3=4;
public:
  //if there are n bins,
  static int GetNumberOfBins2() {return Bin2;}//bins in E
  static int GetNumberOfBins3() {return Bin3;}//bins in partial sums
  //there are n+1 walls.
  int Get_size_2() const {return Bin2+1;}//size of E array
  int Get_size_3() const {return Bin3+1;}//size of partial sum array
  T3R_RW GetT3R_RW(){return ecs;}
  T3NSGangular_RW GetT3NSGangular_RW(){return trw;}
  std::array<double, Bin2+1> GetTailCS(){return TailCS;}
  std::array<double, Bin2+1> GetRutherfordTailCS(){return RutherfordTailCS;}
//WRITE TO FILE
//1. save the file (E,CS):
  void Save_to_CS(int Z, int A/*deuteron Z and A*/, T3String suffix="approx"/*approximation*/);//CS -cross section
//1. read the file (E,CS):
  bool Load_from_CS(int Z, int A/*deuteron Z and A*/, T3String suffix="approx"/*approximation*/);
//2. save the file (E, partial sums):
  void Save_to_PS(/*target Z and A*/int tgZ, int tgA, const T3String & rid,
                          T3int pdg, T3int incZA);//PS - partial sums
//2. read the file (E, partial sums):
  bool Load_from_PS(/*target Z and A*/int tgZ, int tgA, const T3String & rid,
                          T3int pdg, T3int incZA);//PS - partial sums
  //may be called after Fill() was called:
  std::array<double, Bin2+1> GetE() const {return Energy;};
  //may be called after Fill() was called:
  double GetCS(double Tls/*LS kinetic energy of the inc deuteron in inner units (MeV)*/) const;

  friend std::ostream& operator<<(std::ostream& os, const T3R_DDCS& inst);
 
private:
  const double md=1.8756*GeV;
  const double md2=md*md;
  //double Tls;     //kinetic energy of the incident deuteron
  //double t;       //|t| - Mandelstahm variable
  //we chose the range 100 eV - 20 MeV:
  double Emin=100.0 * eV;//minimum LS kinetic energy of the inc deuteron
//in our d-d approximation the maximim LS kinetic energy of the inc deuteron is 231.8 MeV,
//so we take the maximum energy =240 MeV.
  double Emax=240.0 * MeV;//maximum LS kinetic energy of the inc deuteron
  double deltalncos;//the step of logarithmic 1-cos(theta_cm) axis
  double MinLnCosGlobal;//the minimum value 1-cos(theta_cm) = tmin/(2*pcm^2(240 MeV)).
  double deltalnE;//the step of logarithmic lnE axis for E in partial sums
  std::array<double, Bin2+1> lnE;//make logarithmic E axis for partial sums db.
  std::array<double, Bin3+1> lncos;//make logarithmic 1-cos(theta_cm) for partial sums db.
////
  static constexpr double Edisplace_deuteron = 10 * eV;
////
//--------------------------------------------------------------------------//
//1st file (E,CS):
//--------------------------------------------------------------------------//   
  std::array<double, Bin2+1> Energy;//energies:
  std::array<double, Bin2+1> CS;    //d-d elastic cross sections at each energy
//--------------------------------------------------------------------------//
//2nd file (E, partial sums):
//--------------------------------------------------------------------------//
  //here the surface of each trapeze is stored:
  std::array<std::array<double, Bin3+1>, Bin2+1> Sk;//for surface of each trapeze
  std::array<std::array<double, Bin3+1>, Bin2+1> Fk;//for partial sums
  std::array<std::array<double, Bin3+1>, Bin2+1> TailFk;//for partial sums of tail
  std::array<std::array<double, Bin3+1>, Bin2+1> RutherfordTailFk;//for partial sums of Rutherford tail
  std::array<double, Bin2+1> TailCS;//for integral of approximation tail
  std::array<double, Bin2+1> RutherfordTailCS;//for integral of Ritherford tail
//--------------------------------------------------------------------------//
  T3R_RW ecs;
  T3NSGangular_RWrecord tps;
  T3NSGangular_RW trw;
};
} // namespace t3
#endif // T3R_DDCS_H
