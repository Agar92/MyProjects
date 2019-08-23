#include <cstdlib>
#include <algorithm>
#include <assert.h>
#include <math.h>

#include <chrono>
#include <iostream>
#include <vector>
#include <unistd.h>
#include <fstream>
#include <omp.h>
#include <cstring>
#include <cstdint>
#include <climits>
#include <array>

#include "T3Defs.h"
#include "T3R_DDCS.h"
//////////////////////////////////////////////////////////////////////////////
//how to compile:
//gow to /home/70-gaa/NFbuild_script/build/tproc/intel/Release
//if the BUILD_CONFIGURATION is Release and the compiler is intel.
//$cmake --build (--build instead of make)
//$make
//$ctest
//////////////////////////////////////////////////////////////////////////////
// Program main
//////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
  //check Save_to:
  std::ofstream out_stream;
//std::ios::binary - without formatting
//std::ios::trunc - if anything was in the file it is deleted
  out_stream.open("out.dat",
                  std::ios::out | std::ios::binary | std::ios::trunc);

  out_stream.close();

  t3::T3R_DDCS addcs=t3::T3R_DDCS();
  addcs.Fill();
// /////////////////////////////
  const int N2=NumberOfBinsInDDApproximation2;
  const int N3=NumberOfBinsInDDApproximation3;
  std::array<double, N2+1> Energy=addcs.GetE();
  std::array<double, N2+1> TailCS=addcs.GetTailCS();
  std::array<double, N2+1> RutherfordTailCS=addcs.GetRutherfordTailCS();
  std::cout<<"Approx vs Rutherford:"<<std::endl;
  for(int i=0; i<N2+1; ++i)
  {
    std::cout<<"E="<<Energy.at(i)<<" "<<TailCS.at(i)/t3::barn<<" "
             <<RutherfordTailCS.at(i)/t3::barn<<std::endl;
  }
  std::ofstream out("dd_rutherford_approximation_tail_cs.dat");
  for(int i=0; i<N2+1; ++i)
  {
      std::cout<<i+1<<"   "<<Energy.at(i)/t3::MeV<<"   "<<TailCS.at(i)/t3::barn
               <<"   "<<RutherfordTailCS.at(i)/t3::barn<<std::endl;
      out<<i+1<<"   "<<Energy.at(i)/t3::MeV<<" "<<TailCS.at(i)/t3::mbarn
         <<" "<<RutherfordTailCS.at(i)/t3::mbarn<<std::endl;
  }
  out.close();
  exit(0);

// ////////////////////////////
  //acs.Save_to("file1.dat");
  //acs.Load_from("file1.dat");

  //d-d approximation:

  //t3::T3R_RW new1=acs.GetT3R_RW();
  //std::cout<<" NEW1: "<<new1<<" size="<<new1.size()<<std::endl;
  //acs.Save_to_CS(1,2,"approx");
  //t3::T3R_DDCS readtrw1;
  //std::cout<<"STEP 1"<<std::endl;
  //readtrw1.Load_from_CS(1,2,"approx");
  //std::cout<<"STEP 2"<<std::endl;
  //t3::T3R_RW trw1=readtrw1.GetT3R_RW();
  //std::cout<<"HERE WE ARE: "<<trw1<<std::endl;


  addcs.Save_to_PS(1,2,"elastic",1002,1002);
  t3::T3NSGangular_RW nsgan=t3::T3NSGangular_RW();
  t3::T3R_DDCS readtrw2;

  readtrw2.Load_from_PS(1,2,"elastic",1002,1002);
  std::vector<t3::T3NSGangular_RWrecord> v1=readtrw2.GetT3NSGangular_RW();
  std::cout<<"records: "<<v1.size()<<std::endl;
  for(int i=0; i<v1.size(); ++i)
  {
    for(int j=0; j<v1.at(i).size(); ++j) std::cout<<v1.at(i).at(j)<<std::endl;
  }
  std::cout<<"readtrw2:   "<<readtrw2.GetT3NSGangular_RW()<<std::endl;


  exit(0);
  
  //std::cout<<acs.GetdSigmadt(0.001 * MeV, 0.00000001 * GeV * GeV)/mbarn*GeV*GeV<<" mbarn"<<std::endl;
  //std::cout<<acs.GetdSigmadt(0.001 * MeV, 0.000001 * GeV * GeV)/mbarn*GeV*GeV<<" mbarn"<<std::endl;
  //1.71857e+18 mbarn/GeV^2
  //std::cout<<acs.GetdSigmadt(0.0001 * MeV, 3.7512e-08 * GeV * GeV)/mbarn*GeV*GeV<<" mbarn/GeV^2"<<std::endl;
  //1.90383e+17 mbarn/GeV^2
  //std::cout<<acs.GetdSigmadt(0.0001 * MeV, 1.8756e-07 * GeV * GeV)/mbarn*GeV*GeV<<" mbarn/GeV^2"<<std::endl;
  //------------//
  //8.43571e+16 mbarn/GeV^2
  //std::cout<<acs.GetdSigmadt(0.00211474 * MeV, 3.7512e-08 * GeV * GeV)/mbarn*GeV*GeV<<" mbarn/GeV^2"<<std::endl;
  //8.33382e+13 mbarn/GeV^2
  //std::cout<<acs.GetdSigmadt(0.00211474 * MeV, 1.18623e-06 * GeV * GeV)/mbarn*GeV*GeV<<" mbarn/GeV^2"<<std::endl;
  //2.01306e+13 mbarn/GeV^2
  //std::cout<<acs.GetdSigmadt(0.00211474 * MeV, 3.96641e-06 * GeV * GeV)/mbarn*GeV*GeV<<" mbarn/GeV^2"<<std::endl;
  //------------//
  //3.99136e+15 mbarn/GeV^2
  //std::cout<<acs.GetdSigmadt(0.0447214 * MeV, 3.7512e-08 * GeV * GeV)/mbarn*GeV*GeV<<" mbarn/GeV^2"<<std::endl;
  //3.98317e+12 mbarn/GeV^2
  //std::cout<<acs.GetdSigmadt(0.0447214 * MeV, 1.18623e-06 * GeV * GeV)/mbarn*GeV*GeV<<" mbarn/GeV^2"<<std::endl;
  //3.70185e+09 mbarn/GeV^2
  //std::cout<<acs.GetdSigmadt(0.0447214 * MeV, 3.7512e-05 * GeV * GeV)/mbarn*GeV*GeV<<" mbarn/GeV^2"<<std::endl;
  //2.17858e+09 mbarn/GeV^2
  //std::cout<<acs.GetdSigmadt(0.0447214 * MeV, 8.38794e-05 * GeV * GeV)/mbarn*GeV*GeV<<" mbarn/GeV^2"<<std::endl;
  //------------//
  //1.88743e+14 mbarn/GeV^2
  //std::cout<<acs.GetdSigmadt(0.945742 * MeV, 3.7512e-08 * GeV * GeV)/mbarn*GeV*GeV<<" mbarn/GeV^2"<<std::endl;
  //1.89211e+11 mbarn/GeV^2
  //std::cout<<acs.GetdSigmadt(0.945742 * MeV, 1.18623e-06 * GeV * GeV)/mbarn*GeV*GeV<<" mbarn/GeV^2"<<std::endl;
  //2.05774e+08 mbarn/GeV^2
  //std::cout<<acs.GetdSigmadt(0.945742 * MeV, 3.7512e-05 * GeV * GeV)/mbarn*GeV*GeV<<" mbarn/GeV^2"<<std::endl;
  //1.61906e+06 mbarn/GeV^2
  //std::cout<<acs.GetdSigmadt(0.945742 * MeV, 0.00118623 * GeV * GeV)/mbarn*GeV*GeV<<" mbarn/GeV^2"<<std::endl;
  //1.45206e+06 mbarn/GeV^2
  //std::cout<<acs.GetdSigmadt(0.945742 * MeV, 0.00177383 * GeV * GeV)/mbarn*GeV*GeV<<" mbarn/GeV^2"<<std::endl;
  //------------//
  //8.92371e+12 mbarn/GeV^2
  //std::cout<<acs.GetdSigmadt(20. * MeV, 3.7512e-08 * GeV * GeV)/mbarn*GeV*GeV<<" mbarn/GeV^2"<<std::endl;
  //8.90967e+09 mbarn/GeV^2
  //std::cout<<acs.GetdSigmadt(20. * MeV, 1.18623e-06 * GeV * GeV)/mbarn*GeV*GeV<<" mbarn/GeV^2"<<std::endl;
  //8.58607e+06 mbarn/GeV^2
  //std::cout<<acs.GetdSigmadt(20. * MeV, 3.7512e-05 * GeV * GeV)/mbarn*GeV*GeV<<" mbarn/GeV^2"<<std::endl;
  //98755.9 mbarn/GeV^2
  //std::cout<<acs.GetdSigmadt(20. * MeV, 0.00118623 * GeV * GeV)/mbarn*GeV*GeV<<" mbarn/GeV^2"<<std::endl;
  //
  //std::cout<<acs.GetdSigmadt(20. * MeV, 0.037512 * GeV * GeV)/mbarn*GeV*GeV<<" mbarn/GeV^2"<<std::endl;
  
  //1.45206e+06 mbarn/GeV^2
  //std::cout<<acs.GetdSigmadt(0.945742 * MeV, 0.00177383 * GeV * GeV)/mbarn*GeV*GeV<<" mbarn/GeV^2"<<std::endl;
  
  //std::cout<<"Bin="<<T3R_DDCS::GetNumberOfBins()<<std::endl;
  
  //std::ofstream ouputT3TabulatedCS;
  //ouputT3TabulatedCS.open("ouputT3TabulatedCS.txt");
  //tcs.Save_to(ouputT3TabulatedCS);
  //ouputT3TabulatedCS.close();
  
  return 0;
}
