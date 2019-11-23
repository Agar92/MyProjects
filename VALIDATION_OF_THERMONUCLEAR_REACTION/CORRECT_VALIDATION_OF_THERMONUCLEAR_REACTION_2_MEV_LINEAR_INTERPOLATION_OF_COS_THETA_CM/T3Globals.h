#pragma once
#ifndef T3GLOBALS_H
#define T3GLOBALS_H

#include "T3Defs.h"
#include "T3RNG.h"

#ifdef OPENACC
#include <accelmath.h>
#include <openacc.h>
#ifdef CUDA
#include <cuda.h>
#include <cuda_runtime.h>
#endif
#else
#include <omp.h>
  struct float3{
    float x, y, z;
  };
  struct float4{
    float x, y, z, w;
  };
#endif

//KNL ip: 10.41.12.79.

//#define PROBLEM

constexpr bool HISTOGRAM=true;

namespace t3{

  using FloatingType=double;
  using RandGen=t3::RNDGenerator;

//BUILD OPTIONS:  
//Intel Core i7, Nvidia GeForce GTX 650 Ti:
  
//variables:
//at GL=2000000 call to cuMemAlloc returned error 2: Out of memory
//at GL=100000; fails with:
//Failing in Thread:1
//call to cuStreamSynchronize returned error 700: Illegal address during kernel execution

//at GL>2500000 out of memory,
//but if at GL <=2500000 works.
//BUT IF TO OPEN MORE TNAN 1 EMACS WINDOWS (AT LEAST 2 WINDOWS)
//the program fails at runtime with:
//Present table dump for device[1]: NVIDIA Tesla GPU 0, compute capability 3.0, threadid=
//call to cuMemAlloc returned error 2: Out of memory
//I DON'T KNOW WHY.

//Intel KNL, Nvidia Titan V:
//cmake . -DCMAKE_C_COMPILER=pgcc -DCMAKE_CXX_COMPILER=pgc++
//-DCMAKE_C_FLAGS="-acc -Minfo=acc -mcmodel=medium -ta=tesla:cc70
//-tp=haswell -Mnollvm -Minline -Mcuda=cuda10.1"
//-DCMAKE_CXX_FLAGS="-acc -Minfo=acc -mcmodel=medium -ta=tesla:cc70
//-tp=haswell -Mnollvm -Minline -Mcuda=cuda10.1"
//-DCMAKE_CXX_STANDARD=17 -DACC=ON -DCUDA=ON

  const unsigned int GL=10000000;
  constexpr bool report = true;//print output on main() of nbody.cpp.
  constexpr long int G = 27;//2000;
  constexpr int N = 99999;//9999;//9;
  constexpr int Np = N+1;
  constexpr int INJ = 1000;
  constexpr long int K = N+1;//GL;
  constexpr unsigned int max_loop = 10;
  constexpr int cuba = 16;
  constexpr int cubn = cuba + cuba;
  constexpr int cub2 = cubn * cubn;
  constexpr int cub3 = cub2 * cubn;
  constexpr FloatingType TLS=10.0*units::MeV;
  constexpr double fuse = .25;

  //Nbin=8 for i7 (4 physical core, each physical core has 2 logical cores,
  //so there are totally 4*2=8 cores).
  //Nbin=64 for KNL (64 physical cores, hyperthreading may be on or off,
  //do not know exactly).
  //Nbin=16 for a 16-core CPU.
  //Nbin=36 for a 36-core CPU (VKPP).
  constexpr unsigned int Nbin = 64;//1;
  constexpr int DN = (N+1) / Nbin +1;
  constexpr unsigned int BLt = GL/Nbin;//200000000/Nbin;//GL/Nbin;//2000 / Nbin;
  constexpr unsigned int GL1 = BLt - 1;
  constexpr double cgam = 5.0;
  constexpr auto dcgam = cgam + cgam;

  //For filling the histogram:
  constexpr int BinNumber1=100;
  constexpr int BinNumber2=200;
  //The length of a voxel rib:
  constexpr FloatingType ag = 1.0e-5 * units::cm;//the width of the cell (voxel).
  constexpr FloatingType InitParticlex0 = 0.5;
  constexpr FloatingType InitParticley0 = 0.5;
  //End of for filling the histogram.
  
  typedef int PDG_t;
  typedef int MatID_t;

}//end of namespace t3.

#endif//T3GLOBALS_H

//Beginning with PGI 16.1 pgcollect has been superseded. Its functionally is now part of PGPROF. 
//To profile an application, enter this command: 
//pgprof -o a.prof ./a.out 
//To read the profile data, enter this command: 
//pgprof -i a.prof

//pgprof --cpu-profiling-mode top-down -i test.prof

//Can not use pgprof GUI in the same GPU without sudo:
//==36714== Warning: The user does not have permission to profile on the target device. See the following link for instructions to enable permissions and get more information: https://developer.nvidia.com/NVSOLN1000 

