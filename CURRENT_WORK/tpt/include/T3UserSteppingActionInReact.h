#pragma once
#ifndef T3USERSTEPPINGACTIONINREACT_H
#define T3USERSTEPPINGACTIONINREACT_H

#include <iostream>
#include "T3Globals.h"
#include "T3ParticleTable.h"

using namespace t3;

//NUMBER OF BINS IN THE HISTOGRAM:
constexpr int HNbinReact=1024;

class UserSteppingActionInReact
{
public:
  UserSteppingActionInReact():deltacoscmReact(2.0/HNbinReact),
    deuteronPDG(t3::ParticleTable().makePDGfromZandA(1, 2)),
    protonPDG(t3::PDG_t(2212)),neutronPDG(t3::PDG_t(2112))
  {
    _IsRegistered=false;
    for(int m=0; m<HNbinReact; ++m) HistogramInelasticDD[m]=0;
  }
  void Register(){_IsRegistered=true;}
  bool IsRegistered(){return _IsRegistered;}
  void UserAction(Particle<FloatingType> & particlei)
  {
    if(particlei.pdg==neutronPDG || particlei.pdg==protonPDG)
      FillHistogramThetaInelasticDD(particlei);
      
  }
  void FillHistogramThetaInelasticDD(Particle<FloatingType> & particlei)
  {
    const double costcm=particlei.tr;
    const int bin=(costcm+1.0)/deltacoscmReact;
#ifdef OPENACC
#pragma acc atomic update
#else
#pragma omp atomic update
#endif
    ++HistogramInelasticDD[bin];
  }
  void Histogram_theta_UserSteppingActionInReact()
  {
#ifdef OPENACC
#pragma acc update host(HistogramInelasticDD[0:HNbinReact])
#endif
    std::ofstream foutne_theta;
    foutne_theta.open("ne_thetaineldd.dat");
    for(int m=0; m<HNbinReact; ++m)
    {
      foutne_theta<<m<<"   ";
      const double x=-1.0+deltacoscmReact*(m+0.5);
      foutne_theta<<std::setw(8)<<x<<"   "<<static_cast<double>(HistogramInelasticDD[m])/Np/deltacoscmReact<<"   ";
      foutne_theta<<std::endl;
    }
    foutne_theta.close();
  }

private:
  //histogram:
  const double deltacoscmReact;
  int HistogramInelasticDD[HNbinReact];
  bool _IsRegistered;
  const t3::PDG_t deuteronPDG;
  const t3::PDG_t protonPDG;
  const t3::PDG_t neutronPDG;
};

#endif//T3USERSTEPPINGACTIONINREACT_H
