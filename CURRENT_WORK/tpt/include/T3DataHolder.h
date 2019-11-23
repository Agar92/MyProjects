#pragma once
#ifndef T3DATAHOLDER_H
#define T3DATAHOLDER_H

#include <chrono>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstring>
//#include <cstdio>
#include "unistd.h"
#include <type_traits>

#include "T3Globals.h"
#include "T3RNG.h"
#include "T3Particle.h"
#include "T3ParticleTable.h"
#include "T3MaterialTable.h"
#include "T3ThreeVector.h"
#include "T3Process.h"
#include "T3ProcessImplFromCSandFS.h"
#include "T3MultipleScatteringCSImpl.h"
#include "T3MultipleScatteringFSImpl.h"


#include "T3UserSteppingActionInPropagate.h"
#include "T3UserSteppingActionInReact.h"
#include "T3UserSteppingActionInInject.h"



#include "T3InelasticddCSImpl.h"
#include "T3InelasticddFSImpl.h"
#include "T3ElasticEMIonIonCSImpl.h"
#include "T3ElasticEMIonIonFSImpl.h"
#include "T3ElasticStrongIonIonCSImpl.h"
#include "T3ElasticStrongIonIonFSImpl.h"


using namespace t3;

//e - blue p - red photon - green n - goluboi

// build with: cmake . -DCMAKE_C_COMPILER=pgcc -DCMAKE_CXX_COMPILER=pgc++
// -DCMAKE_CXX_FLAGS="-acc -mcmodel=medium -ta=tesla:cc30 -Mnollvm -Minline
// -Mcuda=cuda10.1" -DCMAKE_CXX_STANDARD=17 -DACC=ON -DCUDA=ON

#define MIN(a, b) ((a<b)?a:b)

//=======================================================================//
using Multiple_scattering_t =
  t3::Process<t3::ProcessImplFromCSandFS<t3::MultipleScatteringCS<FloatingType>,
                                         t3::MultipleScatteringFS<true, FloatingType>>>;

// typename is necessary for the compiler
// Multiple_scattering_t::Base_t::CS_t()
// and Multiple_scattering_t::Base_t::FS_t()
// are the anonymous temporary objects (r-values).
Multiple_scattering_t multiplescatteringProcess =
      Multiple_scattering_t(typename Multiple_scattering_t::Base_t::CS_t(),
                   typename Multiple_scattering_t::Base_t::FS_t());


//=======================================================================//
//For T3InelasticddCS/FSImpl.h:
using Inelasticdd_scattering_t =
  t3::Process<t3::ProcessImplFromCSandFS<t3::InelasticddCS<FloatingType>,
                                         t3::InelasticddFS<FloatingType>>>;
                                         
//*                                         
Inelasticdd_scattering_t inelasticddProcess =
      Inelasticdd_scattering_t(typename Inelasticdd_scattering_t::Base_t::CS_t(),
                               typename Inelasticdd_scattering_t::Base_t::FS_t());
//*/
//=======================================================================//

//=======================================================================//
//For T3ElasticEMIonIonCS/FSImpl.h:
using ElasticEMIonIon_scattering_t =
  t3::Process<t3::ProcessImplFromCSandFS<t3::T3ElasticEMIonIonCS<FloatingType>,
                                         t3::T3ElasticEMIonIonFS<FloatingType>>>;
                                         
//*                                         
ElasticEMIonIon_scattering_t ElasticEMIonIonProcess =
      ElasticEMIonIon_scattering_t(typename ElasticEMIonIon_scattering_t::Base_t::CS_t(),
                                   typename ElasticEMIonIon_scattering_t::Base_t::FS_t());
//*/
//=======================================================================//

//=======================================================================//
//For T3ElasticStrongIonIonCS/FSImpl.h:
using ElasticStrongIonIon_scattering_t =
  t3::Process<t3::ProcessImplFromCSandFS<t3::T3ElasticStrongIonIonCS<FloatingType>,
                                         t3::T3ElasticStrongIonIonFS<FloatingType>>>;
                                         
//*                                         
ElasticStrongIonIon_scattering_t ElasticStrongIonIonProcess =
      ElasticStrongIonIon_scattering_t(typename ElasticStrongIonIon_scattering_t::Base_t::CS_t(),
                                   typename ElasticStrongIonIon_scattering_t::Base_t::FS_t());
//*/
//=======================================================================//

//=======================================================================//

Particle<FloatingType> particles[GL] __attribute__((aligned(64)));
int ind01[Nbin][BLt] __attribute__((aligned(64)));
int ind23[Nbin][BLt] __attribute__((aligned(64)));
Particle<FloatingType> arr1[GL] __attribute__((aligned(64)));
Particle<FloatingType> arr2[GL] __attribute__((aligned(64)));
Particle<FloatingType> arr3[GL] __attribute__((aligned(64)));
long int MAX_ELEMENT;
unsigned int SHIFT;
int INJECTED_PARTICLES=0;

unsigned int POSITION3;
unsigned int POSITION2;
unsigned int POSITION1;
unsigned int POSITION0;
unsigned int POSITION23;
int LIFE=0;
unsigned int sizep=sizeof(Particle<FloatingType>);
int push=0;
int over=0;
unsigned int Ntop=0;
double SumDGam=0.;
double GAMMA=0.;
double NoNew=0.0;

decltype(INJECTED_PARTICLES) GetNumOfInjectedParticles(){return INJECTED_PARTICLES;}
decltype(LIFE) GetNumOfAliveParticles(){return LIFE;}
decltype(NoNew) GetNoNew(){return NoNew;}
decltype(SumDGam) GetSumDGam(){return SumDGam;}
decltype(Ntop) GetNtop(){return Ntop;}

ParticleTable aParticleTable;
MaterialTable aMaterialTable;
FloatingType csMultipleScattering[GL];
t3::PDG_t outPDG1[GL];
t3::T3LorentzVector<FloatingType> outP1[GL];
t3::PDG_t outPDG2[GL];
t3::T3LorentzVector<FloatingType> outP2[GL];
//array for T3MultipleScatteringFSImpl.h for storing csBorder arrays for each particle.
FloatingType csBorderDataFS[GL];//this is for GetFS() in Reactor().

//-------------------------------------------------
//For filling the histogram:
//-------------------------------------------------
int HTheta[cubn][BinNumber1];
FloatingType OxTheta[BinNumber1];
const double HThetaMin=1.0e-4, HThetaMax=5.0e-2;
double lnHThetaMin, lnHThetaMax;
double deltaTheta;
//-------------------------------------------------
int HThetax[cubn][BinNumber2];
FloatingType OxThetax[BinNumber2];
const double HThetaxMin=-0.05, HThetaxMax=0.05;
double deltaThetax;
//-------------------------------------------------
int HThetay[cubn][BinNumber2];
FloatingType OxThetay[BinNumber2];
const double HThetayMin=-0.05, HThetayMax=0.05;
double deltaThetay;
//-------------------------------------------------
//End of for filling the histogram.
//-------------------------------------------------




//*
//-------------------------------------------------
//End of for filling the histogram.
//-------------------------------------------------
//For UserSteppingAction Histogram:
  //for theta:
int HistogramTheta[cubn][HNbinPropagate]{0};
//thetai:
FloatingType Thetai[HNbinPropagate]{0.0};
//ThetaMin, ThetaMax:
const double HistogramThetaMin=1.0e-4, HistogramThetaMax=5.0e-2;
//lnThetaMin, lnThetaMax:
double lnHistogramThetaMin, lnHistogramThetaMax;
double deltaLnTheta;

//-------------------------------------------------
//End of for filling the histogram.
//-------------------------------------------------
//*/



int mini[Nbin];
int count01[Nbin];
int count23[Nbin];
int count0[Nbin];
int count1[Nbin];
int count2[Nbin];
int pos2[Nbin];
int count3[Nbin];
int ii1[Nbin];
int ii3[Nbin];
int ii23[Nbin];

int init[Nbin];
int fin[Nbin];

int pointer1[Nbin];
int pointer2[Nbin];
int pointer3[Nbin];
int dL;
int DL;
int n;
int numbin;

template <typename Floating>
class DataHolder
{
public:
  DataHolder()
  {
//Initialize arrays, necessary for filling the histogram:
//-------------------------------------------------------
    lnHThetaMin=log(HThetaMin);
    lnHThetaMax=log(HThetaMax);
    deltaTheta = (lnHThetaMax - lnHThetaMin)/BinNumber1;
    for(int m=0; m<BinNumber1; ++m)
      OxTheta[m]=exp(lnHThetaMin+deltaTheta*(m+0.5));
//-------------------------------------------------------
    deltaThetax = (HThetaxMax - HThetaxMin)/BinNumber2;
    for(int m=0; m<BinNumber2; ++m)
      OxThetax[m]=HThetaxMin+deltaThetax*(m+0.5);
    deltaThetay = (HThetayMax - HThetayMin)/BinNumber2;
    for(int m=0; m<BinNumber2; ++m)
      OxThetay[m]=HThetayMin+deltaThetay*(m+0.5);
//-------------------------------------------------------
    //make UserAction() using UserSteppingActionInPropagate:
    actionInPropagate=UserSteppingActionInPropagate(HistogramThetaMin, HistogramThetaMax,
                                                    Thetai, HistogramTheta);


    //std::cout<<"HERE: "<<std::endl;
    //exit(0);
    //std::cout<<"actionInInject.IsRegistered()="<<actionInInject.IsRegistered()<<std::endl;
    
  }
  void Propagate();
  void React();
  void InitParticle();
  void InitParticle(int lfi/*LIFE+i*/, int mei/*MAX_ELEMENT+i*/, ParticleTable & aParticleTable);
  void Compress();
  void Inject();
  //initializes the pointer to actionInPropagate with !0 pointer:
  void RegisterUserSteppingActionInPropagate()
  {
    actionInPropagate.Register();
  }
  void OutputHistogramThetaUserSteppingActionInPropagate(int histtheta[cubn][HNbinPropagate])
  {
    actionInPropagate.OutputHistogramThetaUserSteppingActionInPropagate(histtheta);
  }
  void Histogram_theta_UserSteppingActionInPropagate(int histtheta[cubn][HNbinPropagate],
                                                     double thetai[HNbinPropagate])
  {
    actionInPropagate.Histogram_theta_UserSteppingActionInPropagate(histtheta, thetai);
  }
  void RegisterUserSteppingActionInInject()
  {
    actionInInject.Register();
  }
  void Histogram_theta_UserSteppingActionInInject()
  {
    actionInInject.Histogram_theta_UserSteppingActionInInject();
  }
  void Histogram_theta();
private:
  //this if for filling the histogram and killing the particles in Propagate():
  UserSteppingActionInPropagate actionInPropagate;
  UserSteppingActionInInject    actionInInject;
};

template <typename Floating> void DataHolder<Floating>::Histogram_theta() {
#ifdef OPENACC
#pragma acc data copyout(HTheta)
{
#endif
  std::ofstream foutne_theta1;
  foutne_theta1.open("ne_theta1.dat");
  for(int m=0; m<BinNumber1; ++m)
  {
    foutne_theta1<<m<<"   ";
    for(int j=0; j<cuba; ++j)
    {
      foutne_theta1<<std::setw(3)<<HTheta[j][m]<<" ";
    }
    foutne_theta1<<std::endl;
  }
  foutne_theta1.close();
  std::ofstream foutne_theta2;
  foutne_theta2.open("ne_theta2.dat");
  for(int m=0; m<BinNumber1; ++m)
  {
    foutne_theta2<<m<<"   ";
    for(int j=cuba; j<cubn; ++j)
    {
      foutne_theta2<<std::setw(3)<<HTheta[j][m]<<" ";
    }
    foutne_theta2<<std::endl;
  }
  foutne_theta2.close();
  std::ofstream foutne_thetai;
  foutne_thetai.open("ne_thetai.dat");
  for(int m=0; m<BinNumber1; ++m)
  {
    foutne_thetai<<m<<"   ";
    foutne_thetai<<std::setw(10)<<OxTheta[m]<<" ";
    foutne_thetai<<std::endl;
  }
  foutne_thetai.close();
#ifdef OPENACC
}
#endif
}

template <typename Floating>
void DataHolder<Floating>::Propagate()
{  
//CALCULATE CROSS SECTIONS:
  auto const material=t3::MatID_t(2u);
//  
  ///\\\///multiplescatteringProcess.GetCS(particles, material, LIFE,
  ///\\\///                                csMultipleScattering, aParticleTable, aMaterialTable);
//

  inelasticddProcess.GetCS(particles, material, LIFE,
                           csMultipleScattering, aParticleTable, aMaterialTable);

  

  //IF FloatingType=double, USE 1.0e-10
  //IF FloatingType=float,  USE 1.0e-7
  Floating da;
  if(std::is_same<double, Floating>::value)     da=ag * 1.0e-10;
  else if(std::is_same<float, Floating>::value) da=ag * 1.0e-7;
  
  constexpr Floating rho0 = 1.0;
  constexpr auto lMAX = std::numeric_limits<Floating>::max();
  Floating const ntid2 = aMaterialTable.GetConcentrations(t3::MatID_t(2u));
  Floating const rho_tid2 = aMaterialTable.GetDensity(t3::MatID_t(2u));
  Floating DE=0.0;
  
#ifdef OPENACC
#pragma acc parallel loop gang vector copy(LIFE) present(particles,csMultipleScattering) copy(lMAX,ntid2,rho_tid2)
#else
#pragma omp parallel for simd
#endif
  for(int i=0; i<LIFE; ++i)
  {

    auto En = particles[i].GetEtot();
    
    if(particles[i].ir > 0)
    {
      auto const csMultipleScatteringi = csMultipleScattering[i];
      
      T3ThreeVector<FloatingType> r0(particles[i].r.x(), particles[i].r.y(), particles[i].r.z());
      //l1 l2 l3
      Floating  l1x =
        (particles[i].vx() == 0.)
            ? lMAX
            : ((particles[i].vx() > 0.)
                     ? ((particles[i].ix() + 1) * ag - particles[i].r.x())
                     : (particles[i].ix() * ag - particles[i].r.x())) /
                        particles[i].vx() +
                    da;
      Floating  l1y =
        (particles[i].vy() == 0.)
            ? lMAX
            : ((particles[i].vy() > 0.)
                     ? ((particles[i].jy() + 1) * ag - particles[i].r.y())
                     : (particles[i].jy() * ag - particles[i].r.y())) /
                        particles[i].vy() +
                     da;
      Floating  l1z =
        (particles[i].vz() == 0.)
            ? lMAX
            : ((particles[i].vz() > 0.)
                     ? ((particles[i].kz() + 1) * ag - particles[i].r.z())
                     : (particles[i].kz() * ag - particles[i].r.z())) /
                        particles[i].vz() +
                     da;

      Floating const dEdx = 2.0 * MeV * cm * cm / gr;
      Floating const dEdxfull = dEdx*rho_tid2;
      //mass:
      ///\\\///Floating const range = (particles[i].GetEtot() - aParticleTable.GetMass(particles[i].pdg))/dEdxfull;
      Floating const range = (particles[i].GetEtot() - particles[i].m)/dEdxfull;
      Floating l0 = range;
      l0=1.0*cm;
//      
      ///\\\///Floating const csMScatteringi = csMultipleScatteringi;
      Floating const csMScatteringi = std::numeric_limits<FloatingType>::max();
//      
      Floating const lambdar = 1.0/ntid2/csMScatteringi;
//IMPORTANT!!!
//THIS LEAD TO THE ERROR: R=0 l2=Not a Number
      ///\\\///Floating const R = particles[i].GenerateCanonical();
//IMPORTANT!!!
//THIS IS HOW ITRIED TO FIX IT: 
      Floating R = particles[i].GenerateCanonical();
      if(R<1.0e-10) R+=1.0e-10;

      Floating l2 = fabs(lambdar * log(R));
      Floating l1=MIN(MIN(l1x, l1y), MIN(l1y, l1z));
      int irc = (l0 < l2 && l0 < l1) ? 0 : ((l2 < l1) ? 2 : 1);
      Floating l = MIN(MIN(l1, l2), MIN(l2, l0));
      Floating dl = fabs(l);
      auto const indexz = particles[i].kz();
      if(irc == 1 /*&& fabs(l-l1z)<da*/)
      {
        auto const modparticlei = particles[i].p.R();
        auto const costheta = particles[i].p.z()/modparticlei;
        auto const theta = acos(costheta);

        //std::cout<<"PROP: theta="<<theta<<" kz="<<indexz<<std::endl;
        
        const int id = particles[i].id;
        auto const thetax=atan( particles[i].p.x()/particles[i].p.z() );//=px/pz.
        auto const thetay=atan( particles[i].p.y()/particles[i].p.z() );//=py/pz.
        auto const x=particles[i].r.x() - InitParticlex0 * ag;
        auto const y=particles[i].r.y() - InitParticley0 * ag;
        auto const r=sqrt(x*x+y*y);
        const double logtheta=log(theta);
        const int ThetaBin=(logtheta-lnHThetaMin)/deltaTheta;
 
        if(ThetaBin>=0 && ThetaBin<BinNumber1)
        {
#ifdef OPENACC
        #pragma acc atomic update
#else
        #pragma omp atomic update
#endif
          ++HTheta[indexz+cuba][ThetaBin];
        }
//-------------------------------------------
        const int ThetaxBin=(thetax-HThetaxMin)/deltaThetax;
        if(ThetaxBin>=0 && ThetaxBin<BinNumber2)
#ifdef OPENACC
        #pragma acc atomic update
#else
        #pragma omp atomic update
#endif
          ++HThetax[indexz+cuba][ThetaxBin];          
//-------------------------------------------
        const int ThetayBin=(thetay-HThetayMin)/deltaThetay;
        if(ThetayBin>=0 && ThetayBin<BinNumber2)
#ifdef OPENACC
        #pragma acc atomic update
#else
        #pragma omp atomic update
#endif
          ++HThetay[indexz+cuba][ThetayBin];                   
      }  

      particles[i].r.SetPxPyPzE(particles[i].r.x()+particles[i].vx() * (dl + da),
                                particles[i].r.y()+particles[i].vy() * (dl + da),
                                particles[i].r.z()+particles[i].vz() * (dl + da), 0.0);//t=0.0

      //std::cout<<"i="<<i<<" LIFE="<<LIFE<<std::endl;
      
      bool const out = (particles[i].ix() >= cuba || particles[i].jy() >= cuba ||
                        particles[i].kz() >= cuba || particles[i].ix() < -cuba ||
                        particles[i].jy() < -cuba || particles[i].kz() < -cuba);
      
      Floating loss = 0.;
      if (aParticleTable.IsNucleus(particles[i].pdg))
      {
        loss = dEdxfull * dl;
        //mass:
        ///\\\///if (loss >= particles[i].GetEtot() - aParticleTable.GetMass(particles[i].pdg))
        if (loss >= particles[i].GetEtot() - particles[i].m)
        {
          //mass:
          ///\\\///particles[i].de += particles[i].GetEtot() - aParticleTable.GetMass(particles[i].pdg);
          particles[i].de += particles[i].GetEtot() - particles[i].m;
          particles[i].SetEtot(aParticleTable.GetMass(particles[i].pdg), aParticleTable);       
          particles[i].ir = 0;
        }
        else
        {
          Floating const old_energy = particles[i].GetEtot();
          Floating const new_energy = old_energy - loss;
          particles[i].SetEtot(new_energy, aParticleTable);
          particles[i].de += loss;
          particles[i].ir = irc;         
        }
      }
      if(out)
      {
        particles[i].ir = 0;
      }


//UserSteppingActionInPropagate:
//Collect the histogram using UserSteppingActionInPropagate:
      if(actionInPropagate.IsRegistered())
        actionInPropagate.UserAction(particles[i], HistogramTheta);
//End of UserSteppingActionInPropagate.

      
    }//End if ir>0    
  }//End of i-particle loop
}//End of Propagator

template <typename Floating>
void DataHolder<Floating>::React()
{
  const int CSBORDERSTEP=aMaterialTable.GetcsBorderSTEP();
  const auto material=t3::MatID_t(2);

//  
  ///\\\///multiplescatteringProcess.GetFS(particles, material, POSITION23,
  ///\\\///                                outPDG1, outP1,/* outPDG2, outP2, */csBorderDataFS,
  ///\\\///                                aParticleTable, aMaterialTable, CSBORDERSTEP);
//

  inelasticddProcess.GetFS(particles, material, POSITION23,
                           outPDG1, outP1,/* outPDG2, outP2, */csBorderDataFS,
                           aParticleTable, aMaterialTable, CSBORDERSTEP);
  
#ifdef OPENACC
#pragma acc parallel loop gang vector present(particles)
#else
#pragma omp parallel for simd
#endif
  for (unsigned int i = 0; i < POSITION23; ++i) {
    particles[i].p = outP1[i];
    particles[i].ir = 2;
    particles[i].pdg = outPDG1[i];
  }  
}//End of Reactor

template <typename Floating>
void DataHolder<Floating>::Compress()
{
  
  double sdg=SumDGam;
#ifdef OPENACC
#pragma acc parallel loop gang vector reduction(+:sdg) present(particles)
#else
#pragma omp parallel for simd reduction(+:sdg)//####???? simd ????####
#endif
  for(int i=0; i<LIFE; ++i)
  {
    sdg += particles[i].de;
  }
  SumDGam=sdg;
  dL=LIFE/Nbin;
  DL=dL+1;
  
  unsigned int const n = Nbin - LIFE % Nbin;
    
  POSITION0=POSITION1=POSITION2=POSITION3=POSITION23=0;
#ifdef OPENACC
#pragma acc parallel loop gang vector copy(GL1,dL,DL,n,count01,count23,count0,count1,count2,count3,init,fin)
#else
#pragma omp parallel for simd
#pragma distribute_point
#endif
    for(int b=0; b<Nbin; ++b)    
    {
      count01[b]=GL1;
      count23[b]=0;
      count0[b]=GL1;
      count1[b]=0;
      count2[b]=GL1;
      count3[b]=0;
      if(b<n)
      {
        init[b]=b*dL;
        fin[b]=(b+1)*dL;
      }
      else if(b==n)
      {
        init[b]=n*dL;
        fin[b]=n*dL+DL;
      }
      else if(b>n)
      {
        init[b]=n*dL+DL*(b-n);
        fin[b]=n*dL+DL*(b-n+1);
      }
    }

#ifdef OPENACC
#pragma acc parallel loop gang vector copy(count01,count23,init,fin) present(particles,ind23)
#endif
  {
#ifndef OPENACC
#pragma omp parallel for
#endif
    for(int b=0; b<Nbin; ++b)
    {
#ifdef OPENACC
      //#pragma acc loop vector reduction(+:count01,count23) ///???????????? vector//
#else
      //#pragma omp simd reduction(+:count01,count23)
#endif
      for(int i=init[b]; i<fin[b]; ++i)
      {
        if(particles[i].ir<2) ind23[b][count01[b]--]=i;
        else                  ind23[b][count23[b]++]=i;
      }
    }
  }

#ifdef OPENACC
#pragma acc parallel loop gang copy(count0,count1,count2,count3,count01,count23,mini,ii23,init,fin) present(particles,ind01,ind23)
#else
#pragma omp parallel for
#pragma distribute_point
#endif
    for(int b=0; b<Nbin; ++b)    
    {
      ii23[b]=count23[b]-1;
      mini[b]=GL1-count01[b];
      if(count23[b]<mini[b]) mini[b]=count23[b];
      int js=0;
#ifdef OPENACC
#pragma acc loop vector reduction(+:js)
#else
#pragma omp simd reduction(+:js)
#endif
      for(int j=0; j<mini[b]; ++j)
        if (ind23[b][ii23[b] - j] > ind23[b][GL1 - j]) ++js;
#ifdef OPENACC
#pragma acc loop vector
#else
#pragma omp simd
#endif
      for(int j=0; j<js; ++j) std::swap(particles[ind23[b][ii23[b]-j]],particles[ind23[b][GL1-j]]);
      
      for(int i=init[b]; i<fin[b]; ++i)
      {
        if     (particles[i].ir==0) ind01[b][count0[b]--]=i;
        else if(particles[i].ir==1) ind01[b][count1[b]++]=i;
        else if(particles[i].ir==2) ind23[b][count2[b]--]=i;
        else                        ind23[b][count3[b]++]=i;
      }
    }
#ifdef OPENACC
#pragma acc parallel loop gang copy(count0,count1,count2,count3,mini,ii1,ii3,GL1) present(particles,ind01,ind23)
#else
#pragma omp parallel for
#pragma distribute_point
#endif
    for(int b=0; b<Nbin; ++b)    
    {
      ii1[b]=count1[b]-1;
      mini[b]=GL1-count0[b];
      if(count1[b]<mini[b]) mini[b]=count1[b];
      int js=0;
#ifdef OPENACC
#pragma acc loop vector reduction(+:js)
#else
#pragma omp simd reduction(+:js)
#endif
      for(int j=0; j<mini[b]; ++j)
        if (ind01[b][ii1[b] - j] > ind01[b][GL1 - j]) ++js;
#ifdef OPENACC
#pragma acc loop vector
#else
#pragma omp simd
#endif
      for(int j=0; j<js; ++j) std::swap(particles[ind01[b][ii1[b]-j]],particles[ind01[b][GL1-j]]);
      ii3[b]=count3[b]-1;
      mini[b]=GL1-count2[b];
      if(count3[b]<mini[b]) mini[b]=count3[b];
      js=0;
#ifdef OPENACC
#pragma acc loop vector reduction(+:js)
#else
#pragma omp simd reduction(+:js)
#endif
      for(int j=0; j<mini[b]; ++j)
        if (ind23[b][ii3[b] - j] > ind23[b][GL1 - j]) ++js;
#ifdef OPENACC
#pragma acc loop vector
#else
#pragma omp simd
#endif
      for(int j=0; j<js; ++j) std::swap(particles[ind23[b][ii3[b]-j]],particles[ind23[b][GL1-j]]);
    }
  
  // Reorder the pointers limits
#ifdef OPENACC
#pragma acc parallel loop gang vector reduction(+:POSITION0,POSITION1,POSITION2,POSITION3,POSITION23)
#else
#pragma omp parallel for simd reduction(+:POSITION0,POSITION1,POSITION2,POSITION3,POSITION23)
#endif
  for(int b=0; b<Nbin; ++b)
  {
    count0[b]=GL1-count0[b];
    count2[b]=GL1-count2[b];
    POSITION0+=count0[b];
    POSITION1+=count1[b];
    POSITION2+=count2[b];
    POSITION3+=count3[b];
    POSITION23+=count23[b];
  }

  auto prevLIFE = LIFE;
  SHIFT = LIFE - POSITION0;

  if(HISTOGRAM) LIFE = LIFE - POSITION0;
  else          LIFE = LIFE + POSITION23 - POSITION0;

  /*
  std::cout<<"prevLIFE="<<prevLIFE<<std::endl;
  std::cout<<"P0="<<POSITION0<<" P1="<<POSITION1
           <<" P2="<<POSITION2<<" P3="<<POSITION3
           <<" P23="<<POSITION23<<std::endl;
  std::cout<<"LIFE="<<LIFE<<std::endl;
  */
  

//здесь должно идти слияние мини ящиков в ящики для 0, 1, 2, 3, удаление 0, и перекладывание 3, 2, 1 в исходный ящик
  pointer1[0]=pointer2[0]=pointer3[0]=0;
#ifdef OPENACC
#pragma acc serial loop /*num_gangs(1) vector_length(1)*/ copy(pointer1,pointer2,pointer3)
#endif
  for(int b=0; b<Nbin-1; ++b)
  {
    pointer1[b+1]=pointer1[b]+count1[b];
    pointer2[b+1]=pointer2[b]+count2[b];
    pointer3[b+1]=pointer3[b]+count3[b];
  }

  //DO NOT parallelize or vectorize - undefined behavior
  for(int b=0; b<Nbin; ++b)
  {
#ifdef OPENACC
#pragma acc host_data use_device(particles,arr1,arr2,arr3)
    {
      cudaMemcpy(&arr1[pointer1[b]],&particles[init[b]+count3[b]+count2[b]],count1[b]*sizep,cudaMemcpyDeviceToDevice);
      cudaMemcpy(&arr2[pointer2[b]],&particles[init[b]+count3[b]],count2[b]*sizep,cudaMemcpyDeviceToDevice);
      cudaMemcpy(&arr3[pointer3[b]],&particles[init[b]],count3[b]*sizep,cudaMemcpyDeviceToDevice);
    }
#else
    memcpy(&arr1[pointer1[b]],&particles[init[b]+count3[b]+count2[b]],count1[b]*sizep);
    memcpy(&arr2[pointer2[b]],&particles[init[b]+count3[b]],count2[b]*sizep);
    memcpy(&arr3[pointer3[b]],&particles[init[b]],count3[b]*sizep);
#endif
  }
  
  // слияние ящиков для 1, 2, 3 в массив particles
#ifdef OPENACC
#pragma acc host_data use_device(particles,arr1,arr2,arr3)
  {
    cudaMemcpy(&particles[0],&arr3[0],POSITION3*sizep,cudaMemcpyDeviceToDevice);
    cudaMemcpy(&particles[POSITION3],&arr2[0],POSITION2*sizep,cudaMemcpyDeviceToDevice);
    cudaMemcpy(&particles[POSITION23],&arr1[0],POSITION1*sizep,cudaMemcpyDeviceToDevice);
  }
#else
  memcpy(&particles[0],&arr3[0],POSITION3*sizep);
  memcpy(&particles[POSITION3],&arr2[0],POSITION2*sizep);
  memcpy(&particles[POSITION23],&arr1[0],POSITION1*sizep);
#endif

  if(!HISTOGRAM)
  {
#ifdef OPENACC
#pragma acc parallel loop gang vector copy(SHIFT,MAX_ELEMENT) present(particles)
#else
#pragma omp parallel for simd
#endif
    for(int i=0; i<POSITION3; ++i)
    {
      int i2=i+i;
      particles[i+SHIFT] = particles[i];
      particles[i].rs=static_cast<unsigned int>(MAX_ELEMENT + i2);
      particles[i+SHIFT].rs=static_cast<unsigned int>(MAX_ELEMENT + i2 + 1);
      //SDI DROPPED THESE LINES:
      particles[i].id=MAX_ELEMENT+i2;
      particles[i].ir=-1;
      particles[i+SHIFT].id=MAX_ELEMENT+i2+1;
      particles[i+SHIFT].ir=-1;
      //END OF WHAT SDI DROPPED.
      if (particles[i].pdg == aParticleTable.makePDGfromZandA(1,2))
      {
        particles[i].pdg = aParticleTable.makePDGfromZandA(1,2);
        particles[i+SHIFT].pdg = aParticleTable.makePDGfromZandA(1,2);
      }
      else if (particles[i].pdg == aParticleTable.makePDGfromZandA(22,48))
      {
        particles[i].pdg = aParticleTable.makePDGfromZandA(22,48);
        particles[i+SHIFT].pdg = aParticleTable.makePDGfromZandA(22,48);
      }
      else if (particles[i].pdg == t3::PDG_t(22))
      {
        particles[i].pdg = t3::PDG_t(11);
        particles[i+SHIFT].pdg = t3::PDG_t(-11);
      } else if (particles[i].pdg == t3::PDG_t(11))
      {
        particles[i].pdg = t3::PDG_t(22);
        particles[i+SHIFT].pdg = t3::PDG_t(11);
      } else if (particles[i].pdg == t3::PDG_t(-11))
      {
        particles[i].pdg = t3::PDG_t(22);
        particles[i+SHIFT].pdg = t3::PDG_t(22);
      } else if (particles[i].pdg == t3::PDG_t(2112))
      {
        particles[i].pdg = t3::PDG_t(2112);
        particles[i+SHIFT].pdg = t3::PDG_t(2112);
      }
    }    
    MAX_ELEMENT+=POSITION3+POSITION3;
#ifdef OPENACC
#pragma acc parallel loop gang vector copy(POSITION3,POSITION23,SHIFT,MAX_ELEMENT) present(particles)
#else
#pragma omp parallel for simd
#endif
    for(int i=POSITION3; i<POSITION23; ++i)
    {
      particles[i+SHIFT] = particles[i];
      particles[i+SHIFT].rs=static_cast<unsigned int>(MAX_ELEMENT + i);
      //WHAT SDI DROPPED:
      particles[i].ir=-1;
      particles[i+SHIFT].id=MAX_ELEMENT+i;
      particles[i+SHIFT].ir=-1;
      //END OF WHAT SDI DROPPED.    
      if (particles[i].pdg == aParticleTable.makePDGfromZandA(1,2))
        particles[i+SHIFT].pdg = aParticleTable.makePDGfromZandA(1,2);
      else if (particles[i].pdg == aParticleTable.makePDGfromZandA(22,48))
        particles[i+SHIFT].pdg = aParticleTable.makePDGfromZandA(22,48);
      else if (particles[i].pdg == t3::PDG_t(22))
        particles[i+SHIFT].pdg = t3::PDG_t(11);
      else
        particles[i+SHIFT].pdg = t3::PDG_t(22);
    }
    MAX_ELEMENT+=POSITION23-POSITION3;
  }//end of (!HISTOGRAM)
}//End of Compressor

//PUSH 1 particle:
template <typename Floating> void DataHolder<Floating>::InitParticle() {

#ifdef OPENACC
#pragma acc parallel num_gangs(1) vector_length(1) present(particles[0:GL]) copy(LIFE,MAX_ELEMENT,INJECTED_PARTICLES)
{
#endif
  t3::PDG_t initPDG = /*PDG_t(22);*/aParticleTable.makePDGfromZandA(1, 2);
  auto const m = aParticleTable.GetMass(initPDG);
  auto const Tkinls = TLS;
  auto const InitEls= m+Tkinls;
  const auto pls=sqrt(InitEls*InitEls-m*m);
  particles[LIFE] = Particle<Floating>(
     t3::T3LorentzVector<Floating>(InitParticlex0*ag,
       InitParticley0*ag,-cuba*ag,0.0),
       t3::T3LorentzVector<Floating>(0.0,0.0,pls,InitEls),
     0.0, initPDG, 1., MAX_ELEMENT, 1, MAX_ELEMENT-1, -1.0);
  
#ifdef OPENACC
#pragma acc atomic update
#endif       
  ++LIFE;
#ifdef OPENACC
#pragma acc atomic update
#endif
  ++MAX_ELEMENT;
#ifdef OPENACC
#pragma acc atomic update
#endif
  ++INJECTED_PARTICLES;
#ifdef OPENACC
  
}
#endif
}//End of InitParticle.

//PUSH INJ particles:
template <typename Floating> void DataHolder<Floating>::InitParticle(int lfi/*LIFE+i*/, int mei/*MAX_ELEMENT+i*/, ParticleTable & aParticleTable)
{
  t3::PDG_t initPDG = /*PDG_t(22);*/aParticleTable.makePDGfromZandA(1, 2);
  auto const m = aParticleTable.GetMass(initPDG);
  auto const InitEls= m+TLS;
  const auto pls=sqrt(InitEls*InitEls-m*m);
  particles[lfi] = Particle<Floating>(
       t3::T3LorentzVector<Floating>(InitParticlex0*ag,
                                   InitParticley0*ag,
                                   -cuba*ag,0.0),
       t3::T3LorentzVector<Floating>(0.0,0.0,pls,InitEls),
       0.0, initPDG, 1., mei, 1, mei-1, -1.0);  
}//End of InitParticle.

template <typename Floating> void DataHolder<Floating>::Inject() {  
  double sg = 0.;
  double wt = 0.;
  double sgCompensation = 0.;
  double wtCompensation = 0.;
#ifdef OPENACC
#pragma acc parallel loop gang vector reduction(+ : sg, wt)
#else
#pragma omp parallel for simd reduction(+ : sg, wt)
#endif
  for (unsigned int j = 0; j < LIFE; ++j)
  {
    auto sgi = particles[j].GetEtot() * static_cast<double>(particles[j].wt)-sgCompensation;
    auto sgtemp = sg + sgi;
    sgCompensation = (sgtemp - sg) - sgi;
    sg = sgtemp;
    auto wti = static_cast<double>(particles[j].wt) - wtCompensation;
    auto wttemp = wt + wti;
    wtCompensation = (wttemp - wt) - wti;
    wt = wttemp;
  }
  GAMMA = sg;
  NoNew = wt;

///*
//----------------------------------------------------------------------------//
//UserSteppingActionInInject:
//Kill the deuterons using UserSteppingActionInInject:
  if(actionInInject.IsRegistered())
  {
#ifdef OPENACC
#pragma acc parallel loop gang vector present(particles)
#else
#pragma omp parallel for simd
#endif
    for(int i=0; i<LIFE; ++i)
    {
      actionInInject.UserAction(particles[i]);
    }
  }
//End of UserSteppingActionInInject.
//----------------------------------------------------------------------------//
//*/

  
  static int push = 0;
  if (LIFE > Ntop) Ntop = LIFE;
  if (push < N)
  {
    if (LIFE < K)
    {
      /*
//SERIAL INJECTION:
      for(int i=0; i<INJ; ++i)
      {
        if(push+i<N) InitParticle();      
      }
      push+=INJ;
      */
      //*
//PARALLEL INJECTION:
      int NUMBER_OF_PARTICLES_TO_PUSH=INJ;
      if(push+INJ>Np) NUMBER_OF_PARTICLES_TO_PUSH=Np-push;
      //std::cout<<"NUMBER_OF_PARTICLES_TO_PUSH="
      //         <<NUMBER_OF_PARTICLES_TO_PUSH<<" LIFE="<<LIFE<<" K="<<K
      //         <<" Np-LIFE="<<Np-LIFE<<" push="<<push<<std::endl;
      
#ifdef OPENACC
#pragma acc parallel loop gang vector /*copy(LIFE,MAX_ELEMENT,TLS)*/ present(particles,aParticleTable)
#else
#pragma omp parallel for        
#endif
      for(int i=0; i<NUMBER_OF_PARTICLES_TO_PUSH; ++i)
#ifdef PROBLEM
        InitParticle(LIFE+i, MAX_ELEMENT+i, aParticleTable);
#else
      {
        t3::PDG_t initPDG = aParticleTable.makePDGfromZandA(1, 2);
        auto const m = aParticleTable.GetMass(initPDG);
        auto const InitEls= m+TLS;
        const auto pls=sqrt(InitEls*InitEls-m*m);
        particles[LIFE+i] = Particle<Floating>(
        t3::T3LorentzVector<Floating>(InitParticlex0*ag,
                                    InitParticley0*ag,
                                    -cuba*ag,0.0),
        t3::T3LorentzVector<Floating>(0.0,0.0,pls,InitEls),
        m, 0.0, initPDG, 1., MAX_ELEMENT+i, 1, MAX_ELEMENT+i-1, -1.0);
        ////std::cout<<"###"<<std::endl;
      }
#endif
      push+=NUMBER_OF_PARTICLES_TO_PUSH;
      LIFE+=NUMBER_OF_PARTICLES_TO_PUSH;
      MAX_ELEMENT+=NUMBER_OF_PARTICLES_TO_PUSH;
      INJECTED_PARTICLES+=NUMBER_OF_PARTICLES_TO_PUSH;
    }
  }
 }//End of Injector

#endif//T3DATAHOLDER_H
