#include <chrono>
#include <iostream>
#include "unistd.h"
#include <typeinfo>

#include "T3DataHolder.h"


#include "T3NSGangular_node.hh"
#include "T3Inelasticdd_DB.hh"
#include "T3InelasticddCSImpl.h"
#include "T3InelasticddFSImpl.h"
#include "T3Utility.hh"


#include "T3ElasticEMIonIonCSImpl.h"
#include "T3ElasticEMIonIonFSImpl.h"


//build options on CPU i7:
//Using Intel:
//cmake . -DCMAKE_C_COMPILER=icc -DCMAKE_CXX_COMPILER=icpc
//-DCMAKE_CXX_FLAGS="-march=native -mtune=native -O3 -ipo16
//-mcmodel=large -qopt-report=5"
//-DCMAKE_CXX_STANDARD=17 -DACC=OFF -DCUDA=OFF

//Using PGI:
//cmake . -DCMAKE_C_COMPILER=pgcc -DCMAKE_CXX_COMPILER=pgc++
//-DCMAKE_C_FLAGS="-mcmodel=medium -Minline"
//-DCMAKE_CXX_FLAGS="-mcmodel=medium -Minline"
//-DCMAKE_CXX_STANDARD=17 -DACC=OFF -DCUDA=OFF

//build options on CPU GeForce GTX 650 Ti:
//DO NOT FORGET -Minline!!! WITHOUT IT THE CODE DOES NOT WORK ON GPU!!!
//cmake . -DCMAKE_C_COMPILER=pgcc -DCMAKE_CXX_COMPILER=pgc++
//-DCMAKE_C_FLAGS="-acc -Minfo=acc -mcmodel=medium -ta=tesla:cc30
//-Minline -Mcuda=cuda10.1" -DCMAKE_CXX_FLAGS="-acc -Minfo=acc
//-mcmodel=medium -ta=tesla:cc30 -Minline -Mcuda=cuda10.1"
//-DCMAKE_CXX_STANDARD=17 -DACC=ON -DCUDA=ON


//build options on KNL:
//cmake . -DCMAKE_C_COMPILER=pgcc -DCMAKE_CXX_COMPILER=pgc++
//-DCMAKE_C_FLAGS="-acc -Minfo=acc -mcmodel=medium -ta=tesla:cc70
//-tp=haswell -Mnollvm -Minline -Mcuda=cuda10.1"
//-DCMAKE_CXX_FLAGS="-acc -Minfo=acc -mcmodel=medium -ta=tesla:cc70
//-tp=haswell -Mnollvm -Minline -Mcuda=cuda10.1"
//-DCMAKE_CXX_STANDARD=17 -DACC=ON -DCUDA=ON

//////////////////////////////////////////////////////////////////////////////
// Program main
//////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
#ifdef OPENACC
  #ifdef CUDA
    std::cout<<"OPENACC IS DEFINED, CUDA IS DEFINED"<<std::endl;
  #else
    std::cout<<"OPENACC IS DEFINED, CUDA IS NOT DEFINED"<<std::endl;
  #endif
#else
  #ifdef CUDA
    std::cout<<"OPENACC IS NOT DEFINED, CUDA IS DEFINED"<<std::endl;
  #else
    std::cout<<"OPENACC IS NOT DEFINED, CUDA IS NOT DEFINED"<<std::endl;
  #endif
#endif


  std::cout<<"ENTER:"<<std::endl;
  T3R_DDCS trdcs;
  trdcs.Fill();
  //trdcs.Save_to_CS();
  std::cout<<"EXIT"<<std::endl;
  exit(0);

    

  auto begin=std::chrono::steady_clock::now();

  DataHolder<FloatingType> d;
  //register UserSteppingActionInPropagate:
  d.RegisterUserSteppingActionInPropagate();
  //register UserSteppingActionInInject:
  d.RegisterUserSteppingActionInInject();
  
  if (report)
    std::cout << "size of DataHolder is ~"
              << sizeof(DataHolder<FloatingType>) / float(1ul << 30ul)
              << "GB, size of d is ~" << sizeof(d) / float(1ul << 30ul) << "GB"
              << std::endl;
  
  int count=0;
  
  LIFE=0;
  MAX_ELEMENT=0;
#ifdef OPENACC
#pragma acc data create(ind01,ind23,arr1,arr2,arr3,outPDG1,outPDG2,outP1,outP2,csBorderDataFS,csMultipleScattering) \
  copyin(particles,aParticleTable,aMaterialTable,d,HistogramTheta,Thetai,multiplescatteringProcess,inelasticddProcess)
  //  copyin(HistogramInelasticDD)/*copy(ag,lambda,dtheta0,da,rho0,Dan)*/
  {
#endif
    
    //\\//d.InitParticle();
    ///d.Inject();
    bool InitLoop=true;
    for(unsigned int step=1; GetNumOfAliveParticles()>0 || InitLoop==true || GetNumOfInjectedParticles()<Np; ++step)
    {
      //std::cout<<"#"<<count<<" InitLoop="<<InitLoop<<std::endl;         
      d.Inject();
      d.Propagate();
      d.Compress();
      d.React();
      if(report)
      {
        std::cout << step << "   " << GetNumOfAliveParticles() << "    "
                  << GetNoNew() <<  "   " << GetNumOfInjectedParticles() << " "
                  << GetSumDGam() <<std::endl;
      }
      if(InitLoop) InitLoop=false;
      ++count;
      //std::cout<<"#"<<count<<" InitLoop="<<InitLoop<<" LIFE="<<LIFE<<std::endl;
      //std::cout<<"4 GetNumOfAliveParticles()="<<GetNumOfAliveParticles()
      //         <<" InitLoop="<<InitLoop<<std::endl;
      //sleep(1);
      
    }
    d.Histogram_theta();

    //For check of inelastic d+d it gives zeros:
    ////d.OutputHistogramThetaUserSteppingActionInPropagate(HistogramTheta);
    d.Histogram_theta_UserSteppingActionInPropagate(HistogramTheta, Thetai);

    d.Histogram_theta_UserSteppingActionInInject();

        
#ifdef OPENACC
  }
#endif

  auto end=std::chrono::steady_clock::now();
	auto elapsed_ms=std::chrono::duration_cast<std::chrono::milliseconds>(end-begin);
	std::cout<<"time="<<elapsed_ms.count()<<" ms, G="<<G<<", K="<<K<<", Ntop="
           <<Ntop<<", SumDG="<<SumDGam<<std::endl;
  std::cout<<"Nbin="<<Nbin<<" FloatingType="<<typeid(FloatingType).name()<<std::endl;
}
