#pragma once
#ifndef T3USERSTEPPINGACTIONINPROPAGATE_H
#define T3USERSTEPPINGACTIONINPROPAGATE_H

#include <iostream>
#include "T3Globals.h"

using namespace t3;

//NUMBER OF BINS IN THE HISTOGRAM:
constexpr int HNbinPropagate=100;

class UserSteppingActionInPropagate
{
public:
  UserSteppingActionInPropagate()
  {
    _IsRegistered=false;
  }
  //parameters for fill the histogram:
  UserSteppingActionInPropagate(double HistogramThetaMin/*min*/, double HistogramThetaMax/*max*/,
                                FloatingType Thetai[HNbinPropagate]/*theta in the mids of bins*/,
                                int histtheta[cubn][HNbinPropagate]/*histogram array*/)
  {
    lnHistogramThetaMin=log(HistogramThetaMin);
    lnHistogramThetaMax=log(HistogramThetaMax);
    deltaLnTheta = (lnHistogramThetaMax - lnHistogramThetaMin)/BinNumber1;
    for(int m=0; m<HNbinPropagate; ++m)
      Thetai[m]=exp(lnHistogramThetaMin+deltaLnTheta*(m+0.5));
    for(int i=0; i<cubn; ++i) for(int j=0; j<HNbinPropagate; ++j) histtheta[i][j]=0;
    _IsRegistered=false;
  }

  void Register(){_IsRegistered=true;}
  bool IsRegistered(){return _IsRegistered;}

  void UserAction(Particle<FloatingType> & particle, int histtheta[cubn][HNbinPropagate])
  {
    FillHistogramTheta(particle, histtheta);
  }

  void FillHistogramTheta(Particle<FloatingType> & particle, int histtheta[cubn][HNbinPropagate])
  {
    //std::cout<<"###"<<std::endl;
    if(particle.ir==1 || particle.ir==0)
    {
      if(particle.ir!=1 && particle.ir!=0)
      {
        //std::cout<<"IR!=1 IR!=0"<<std::endl;
        //sleep(3);
      }
      int BinZ=0;
      /*if(particle.kz()<0)*/ BinZ=particle.kz()-1;
      //else                BinZ=particle.kz();

      //std::cout<<"BinZ="<<BinZ<<std::endl;
      //usleep(10000);

      const double modparticlei = particle.p.R();
      const double costheta = particle.p.z()/modparticlei;
      const double theta = acos(costheta);

      /*
      std::cout<<"pdg="<<particle.pdg
             <<" ix="<<particle.ix()
             <<" jy="<<particle.jy()
             <<" kz="<<particle.kz()
             <<" x="<<particle.r.x()
             <<" y="<<particle.r.y()
             <<" z="<<particle.r.z()
             <<" px="<<particle.p.x()
             <<" py="<<particle.p.y()
             <<" pz="<<particle.p.z()
             <<" P="<<modparticlei
             <<std::endl;
      */
      
      //std::cout<<"theta="<<theta<<std::endl;
      if(theta>0.0)
      {
        const int bin=(log(theta)-lnHistogramThetaMin)/deltaLnTheta;
        if(bin>=0 && bin<HNbinPropagate)
        {
      #ifdef OPENACC
      #pragma acc atomic update
      #else
      #pragma omp atomic update
      #endif
          ++histtheta[cuba+BinZ][bin];
          //std::cout<<cuba+BinZ<<" "<<bin<<" "<<histtheta[cuba+BinZ][bin]<<std::endl;
        }
      }
    }
  }

  void OutputHistogramThetaUserSteppingActionInPropagate(int histtheta[cubn][HNbinPropagate])
  {
#ifdef OPENACC
#pragma acc update host(histtheta[0:cubn][0:HNbinPropagate])
#endif
    for(int m=0; m<HNbinPropagate; ++m)
    {
      std::cout<<m<<" ";
      for(int j=0; j<cubn; ++j) std::cout<<histtheta[j][m]<<" ";
      std::cout<<std::endl;
    }
  }

  void Histogram_theta_UserSteppingActionInPropagate(int histtheta[cubn][HNbinPropagate],
                                                     double thetai[HNbinPropagate])
  {
#ifdef OPENACC
#pragma acc update host(histtheta[0:cubn][0:HNbinPropagate],thetai[0:HNbinPropagate])
#endif
    std::ofstream foutne_theta1;
    foutne_theta1.open("ne_thetast1.dat");
    for(int m=0; m<HNbinPropagate; ++m)
    {
      foutne_theta1<<m<<"   ";
      for(int j=0; j<cuba; ++j)
      {
        foutne_theta1<<std::setw(3)<<histtheta[j][m]<<" ";
      }
      foutne_theta1<<std::endl;
    }
    foutne_theta1.close();
    std::ofstream foutne_theta2;
    foutne_theta2.open("ne_thetast2.dat");
    for(int m=0; m<HNbinPropagate; ++m)
    {
      foutne_theta2<<m<<"   ";
      for(int j=cuba; j<cubn; ++j)
      {
        foutne_theta2<<std::setw(3)<<histtheta[j][m]<<" ";
      }
      foutne_theta2<<std::endl;
    }
    foutne_theta2.close();
    std::ofstream foutne_thetai;
    foutne_thetai.open("ne_thetasti.dat");
    for(int m=0; m<HNbinPropagate; ++m)
    {
      foutne_thetai<<m<<"   ";
      foutne_thetai<<std::setw(10)<<thetai[m]<<" ";
      foutne_thetai<<std::endl;
    }
    foutne_thetai.close();
  }

  
private:
  //histogram:
  double lnHistogramThetaMin;
  double lnHistogramThetaMax;
  double deltaLnTheta;
  bool _IsRegistered;
};


#endif//T3USERSTEPPINGACTIONINPROPAGATE_H
