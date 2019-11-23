#include <chrono>
#include <iostream>
#include "unistd.h"
#include <typeinfo>
#include <cmath>
#include <fstream>
#include <iomanip>

//#include "T3Globals.h"
#include "T3Defs.h"
#include "T3RNG.h"

//
//#include "T2NSGangular_RW.hh"
//

#include <iostream>
#include "unistd.h"
#include <iterator>

#include "T3NSGangular_RW.hh"
#include "T3NSGangular.hh"



#include "T3InelasticddCSImpl.h"
#include "T3InelasticddFSImpl.h"
#include "T3ParticleTable.h"

#include "T3Inelasticdd_DB.hh"

#include "T3Particle.h"
#include "T3DataHolder.h"



using namespace t3;

//TOTAL NUMBER OF INJECTED PARTICLES:
const int NPARTICLES=10000000;//2000000000;//1 million

///Particle<FloatingType> particles[NPARTICLES];

//the indexes in the rw.at(0).at(index) database,
//which i will use for taking partial sums at necessary energies:
//10   MeV index=66
//5    MeV index=56
//2    MeV index=39
//1    MeV index=29
//0.1  MeV index=11
//0.01 MeV index=2

//!!!!!!!!!!!!!!!!!!!!!!!!!!!1//
const double kinEnergy=2.0;//0.1;
//!!!!!!!!!!!!!!!!!!!!!!!!!!!1//

//the step is 2/128 and there are 129 points in Ox (cos(theta_CM)) array:
//from -1 to 1, there are 129 points:

//Found that PAW does not allow to plot more than 1024 points
//when i try to plot 2048 bins, it plots on ly the first 1024
//=>i can see only the left zoom and the left side of the center xoom.
constexpr int Nbin=2048;//2048;//1024;//2048;//1024;
constexpr int NpointsInPS=129;

const double delta=static_cast<double>(2.0)/128;

constexpr int NpointsInCosCM=Nbin+1;
//cos(theta_CM):
const double deltacoscm=2.0/Nbin;
double coscm[129];
//differentiated partial sums (dps):
//there are  NpointsInCosCM-1 elements because can not take the derivative of x=1.0.
double dpsmain[NpointsInPS-1]{0.0};
//dpsexp:
double dpsexp[NpointsInPS-1]{0.0};

//Ox (cos(theta_cm)) in the middle of bins (128 bins=>128 elements in Ox):
double Ox[NpointsInCosCM-1]{0.0};

const int CHANNELS=2;//2 channels: exponent and main distribution.

int HistogramNiCosCM[NpointsInCosCM-1]{0};

int InHistogramRange=0;
int OutOfHistogramRange=0;

//there are totally 128 bins.
//the length from -1 to 1 is 2
//points -1 and 1 are included in Ox (cos(theta_CM)),
void MakeEquidistantBinsCosThetaCM()
{
#ifdef OPENACC
#pragma acc parallel loop gang vector present(coscm) copy(delta)
#endif
  for(int i=0; i<129; ++i) coscm[i]=-1.0+i*delta;
}

void FillDifferentiatedPartialSums(double * ps/*read from T2_DATA/ partial sums - 67 elements in ps*/,
                                   double PS_CONTRIBUTION/*contribution of main partial sums distribution=1-pr, where pr is the contribution of the exponent*/
                                  )
{
#ifdef OPENACC
#pragma acc parallel loop gang vector present(dpsmain,ps,CosCM) 
#endif
  for(int m=0; m<NpointsInPS-1; ++m)
  {
    //differentiated partial sums:
    //std::cout<<"m="<<m<<std::endl;
    dpsmain[m]=PS_CONTRIBUTION*((ps[m+1]*ps[m+1]-ps[m]*ps[m]))/delta;

    //std::cout<<"delta="<<delta<<" PS_CONTRIBUTION="<<PS_CONTRIBUTION
    //         <<std::endl;
    //for(int j=0; j<129; ++j) std::cout<<ps[j]<<" ";
    //std::cout<<std::endl;
    //sleep(1);
    
    //std::cout<<"ps.at(m+1)="<<ps.at(m+1)<<" ps.at(m)="<<ps.at(m)<<" ps.at(m+1)-ps.at(m)="<<ps.at(m+1)-ps.at(m)<<" CosCM.at(m+1)-CosCM.at(m)="<<CosCM.at(m+1)-CosCM.at(m)
    //         <<" dpsmain.at(m)="<<dpsmain.at(m)<<std::endl;
  }
}

//in the middle of bins:
void FillOxAxisCosCM() {
  //in the middle of the bin:
  //std::cout<<"Filling Ox:"<<std::endl;

#ifdef OPENACC
#pragma acc parallel loop gang present(Ox,coscm) copy(deltacoscm)
#endif
  for(int m=0; m<NpointsInCosCM; ++m)
  {
    Ox[m]=coscm[0]+(m+0.5)*deltacoscm;
  }
}

void Histogram_OxAxisCosCM() {
  std::ofstream foutne_theta;
  foutne_theta.open("ne_thetacoscm.dat");
#ifdef OPENACC
#pragma acc update host(Ox)
#endif
  for(int m=0; m<NpointsInCosCM-1; ++m)
  {
    foutne_theta<<m<<"   ";
    //in the middle of the bin:
    foutne_theta<<std::setw(8)<<Ox[m]<<"   ";
    foutne_theta<<std::endl;
  }
  foutne_theta.close();
}

void Histogram_DifferentiatedPartialSums() {
  std::ofstream foutne_theta;
  foutne_theta.open("ne_thetadps.dat");
#ifdef OPENACC
#pragma acc update host(dpsmain)
#endif
  for(int m=0; m<NpointsInPS-1; ++m)
  {
    const double deltax=2.0/128;
    const double x=-1.0+deltax*(m+0.5);
    foutne_theta<<m<<"   "<<std::setw(8)<<x<<"   "<<dpsmain[m]<<"   ";
    foutne_theta<<std::endl;
  }
  foutne_theta.close();
}

void FillExponent(double sl, double pr/*pr*/) {

#ifdef OPENACC
#pragma acc parallel loop gang vector present(Ox,dpsexp) copy(sl,pr) 
#endif
  for(int m=0; m<NpointsInPS-1; ++m)
  {
    const double deltax=2.0/128;
    const double x=-1.0+deltax*(m+0.5);
    const double norm=sl*(1.0-exp(-2.0/sl));
    /*
      std::cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<std::endl;
      std::cout<<"slope="<<sl<<std::endl;
      std::cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<std::endl;
    */
    const double fexp=pr*exp(-(1-x)/sl)/norm;
    dpsexp[m]=fexp;
    //std::cout<<"x="<<Ox.at(m)<<" norm="<<norm<<" fexp="<<fexp
    //         <<" dpsexp.at(m)="<<dpsexp.at(m)<<std::endl;
  }
}

void Histogram_Exponent() {
  std::ofstream foutne_theta;
  foutne_theta.open("ne_thetaexp.dat");
#ifdef OPENACC
#pragma acc update host(dpsexp)
#endif
  for(int m=0; m<NpointsInPS-1; ++m)
  {
    foutne_theta<<m<<"   "<<std::setw(8)<<dpsexp[m]<<"   ";
    foutne_theta<<std::endl;
  }
  foutne_theta.close();
}

void FillHistogramNi()
{
  std::ofstream foutne_theta;
  foutne_theta.open("ne_thetaniquadr.dat");
#ifdef OPENACC
#pragma acc update host(HistogramNiCosCM)
#endif
  for(int m=0; m<NpointsInCosCM-1; ++m)
  {
    foutne_theta<<m<<"   ";
    //1)number of events in each bin, normalized by the total number of events.
    //if the surface under the graphic of the function is 1, then this is equal to
    //the surface under the function graphic in each bin:
    //2)Since we take the normalized partial sums from T2_DATA/ and
    //differentiate them (divide by deltacoscm=2.0/Nbin),
    //here we also need to divide normalized Ni/Ntot by deltacoscm=2.0/Nbin,
    //to obtain a derivative. 
    //\\//const double Hm=static_cast<double>(HistogramNiCosCM.at(j).at(m))/NPARTICLES/deltacoscm;
    const double x=-1.0+deltacoscm*(m+0.5);
    foutne_theta<<std::setw(8)<<x<<"   "<<static_cast<double>(HistogramNiCosCM[m])/NPARTICLES/deltacoscm<<"   ";
    foutne_theta<<std::endl;
  }
  foutne_theta.close();
}











//////////////////////////////////////////////////////////////////////////////
// Program main
//////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
  auto begin=std::chrono::steady_clock::now();


  
  t3::RNDGenerator gen(1);
  t3::ParticleTable aParticleTable;
  t3::MaterialTable aMaterialTable;
  const PDG_t positronPDG=t3::PDG_t(-11);
  const PDG_t deuteronPDG=aParticleTable.makePDGfromZandA(1, 2);//deuteron
  t3::MatID_t matID=t3::MatID_t(2);//TiD2.
  t3::InelasticddFS<FloatingType> infs;
  const double md=aParticleTable.GetMass(deuteronPDG);
  const double Tls=2.0 * MeV;
  const double E=md + Tls;
  const double initp0=sqrt(E*E-md*md);
  LorentzVector<FloatingType> InitP=LorentzVector<FloatingType>(0.0, 0.0, initp0, E );
  FloatingType tr=0.0;
  infs.GetFS(InitP, deuteronPDG, matID,
             gen, tr, aParticleTable, aMaterialTable);

  //t3::InelasticddCS<FloatingType> incs;
  //incs.GetCS(Tls, deuteronPDG, matID,
  //           aParticleTable, aMaterialTable);
  
  //exit(0);
//\\//******************************************************************************************************************************************//  






  
//****************************************************************************************//
//Begin of read the partial sums form T2_DATA
//****************************************************************************************//  
  T3NSGangular_RW rw;  
  //target deuteron Z=1
  static constexpr auto tgZ = 1;
  //target deuteron A=2
  static constexpr auto tgA = 2;
  //reaction product=n+3He
  //this is a neutron from a reaction product:
  static constexpr auto sPDG = 2112;
  //1000*Z+A//inc deuteron
  static constexpr auto incZA = 1002;
  //MT index of reaction in ENDF
  static constexpr auto rid = "50";
  //load the data from the binary file:
  rw.load_binary(tgZ, tgA, rid, sPDG, incZA);
  //print the data records (e, partical sums):
  //std::cout << rw << std::endl;
  std::cout<<"rw size="<<rw.size()<<std::endl;
  std::cout<<"SEE:"<<std::endl;
  //std::cout<<rw<<std::endl;
  
//****************************************************************************************//
//End of read the partial sums form T2_DATA
//****************************************************************************************//  
  

  T3NSGangular_RWnode curnode=rw.at(0).interpolate(kinEnergy);


  //std::cout<<curnode<<std::endl;
  //exit(0);

  
  std::array<T3double, 127> psarray=curnode.Get_V();
  double PsArray[129];
  PsArray[0]=0.0;
  PsArray[128]=0.0;
  for(int k=0; k<127; ++k) PsArray[k+1]=psarray.at(k);

  MakeEquidistantBinsCosThetaCM();
  FillOxAxisCosCM();
  //output bin number/Ox axis points to "ne_thetacoscm.dat" file.
  Histogram_OxAxisCosCM();
  FillDifferentiatedPartialSums(PsArray, 1.0-curnode.Get_pr());
  Histogram_DifferentiatedPartialSums();
  //Fill dpsexp:
  FillExponent(curnode.Get_sl(), curnode.Get_pr());
  //output dpsexp to "ne_thetexp.dat":
  Histogram_Exponent();









/*
  
  
  //In T2VIMultyS.cc at line 335 printed on KNL:
  //Z=1 Z+N=2 mt=50 recoilLvl=0 kinEnergy=2 pZA()=1002 N=1

  //G4double const cost = (recoilLvl == 077) ? ea.second :
  //                  T2NSGangular::RandomizeCost(Z, Z + N, 2112, mt, recoilLvl, kinEnergy,
  //                                              pZA());

  RNDGenerator generator(1);
  int mt=50;
  int recoilLvl=0;
  for(int i=0; i<NPARTICLES; ++i)
  {
    const double cost=T3NSGangular::RandomizeCost(generator, tgZ, tgA, 2112, mt, recoilLvl, kinEnergy, incZA);

    const int bin=(cost+1.0)/deltacoscm;

#ifdef OPENACC
#pragma acc atomic update
#endif
    ++HistogramNiCosCM[bin];
    
    if(bin>=0 && bin<NpointsInCosCM-1)
    {
      ++InHistogramRange;
    }
    else
    {
      ++OutOfHistogramRange;
    }
  }

  std::cout<<"InHistogramRange="<<InHistogramRange
           <<" OutOfHistogramRange="<<OutOfHistogramRange<<std::endl;


*/           


//*****************************************************************************//  


  T3NSGangular_RWrecord rwrec=rw.at(0);
  //make square of partial sums:
  for(int i=0; i<rwrec.size(); ++i)
  {
    for(int j=0; j<127; ++j) rwrec.at(i).Set_V(j+1, rwrec.at(i).Get_V(j+1)*rwrec.at(i).Get_V(j+1));
  }

  T3Inelasticdd_DB db;
  const FloatingType Emin = rwrec.front().Get_E();
  const FloatingType Emax = rwrec.back().Get_E();
  const FloatingType Ediff = Emax - Emin;
  const size_t nbins = 512 - 1;
  const FloatingType dE = Ediff / nbins;
  for(size_t i = 0; i < 512; ++i)
  {
    const FloatingType Ei = Emin + i*dE;
    const T3NSGangular_RWnode rwnode = rwrec.interpolate(Ei);    
    T3NSGangular_node nodes_i = T3NSGangular_node(rwnode);
    for(int j=0; j<127; ++j)
    {
      db.Set_V(i, j, nodes_i.Get_V(j));
      db.Set_a(i, j, nodes_i.Get_a(j));
      db.Set_b(i, j, nodes_i.Get_b(j));
      db.Set_c(i, j, nodes_i.Get_c(j));
    }
    const FloatingType pri=nodes_i.Get_pr();
    const FloatingType sli=nodes_i.Get_sl();
    db.Set_pr(i, pri);
    db.Set_sl(i, sli);
    db.Set_Einc(i, Ei);
  }


  
/*
  
  //In T2VIMultyS.cc at line 335 printed on KNL:
  //Z=1 Z+N=2 mt=50 recoilLvl=0 kinEnergy=2 pZA()=1002 N=1

  //G4double const cost = (recoilLvl == 077) ? ea.second :
  //                  T2NSGangular::RandomizeCost(Z, Z + N, 2112, mt, recoilLvl, kinEnergy,
  //                                              pZA());

  RNDGenerator generator(1);
  for(int i=0; i<NPARTICLES; ++i)
  {
    const double cost=db.RandomizeCost(generator, kinEnergy);

    const int bin=(cost+1.0)/deltacoscm;

#ifdef OPENACC
#pragma acc atomic update
#endif
    ++HistogramNiCosCM[bin];
    
    if(bin>=0 && bin<NpointsInCosCM-1)
    {
      ++InHistogramRange;
    }
    else
    {
      ++OutOfHistogramRange;
    }
  }

  std::cout<<"InHistogramRange="<<InHistogramRange
           <<" OutOfHistogramRange="<<OutOfHistogramRange<<std::endl;

*/





  //In T2VIMultyS.cc at line 335 printed on KNL:
  //Z=1 Z+N=2 mt=50 recoilLvl=0 kinEnergy=2 pZA()=1002 N=1

  //G4double const cost = (recoilLvl == 077) ? ea.second :
  //                  T2NSGangular::RandomizeCost(Z, Z + N, 2112, mt, recoilLvl, kinEnergy,
  //                                              pZA());


/*
  

//
//1. Initialize particles array:
//
  for(int i=0; i<NPARTICLES; ++i)
  {
    t3::PDG_t initPDG = aParticleTable.makePDGfromZandA(1, 2);
    auto const m = aParticleTable.GetMass(initPDG);
    auto const InitEls= m+Tls;
    const auto pls=sqrt(InitEls*InitEls-m*m);
    particles[i] = Particle<FloatingType>(
                        t3::LorentzVector<FloatingType>(InitParticlex0*ag,
                        InitParticley0*ag,
                        -cuba*ag,0.0),
                        t3::LorentzVector<FloatingType>(0.0,0.0,pls,InitEls),
                        m, 0.0, initPDG, 1., i, 1, i-1, 0.0);
  }
  
  
//
//2. Fill in the histogram:
//
  RNDGenerator generator(1);
  for(int i=0; i<NPARTICLES; ++i)
  {
    if(i%1000000==0) std::cout<<"i="<<i<<std::endl;
    /////////Four<FloatingType> four=infs.GetFS(InitP, deuteronPDG, matID, generator, aParticleTable, aMaterialTable);
    ///////Four<FloatingType> four=infs.GetFS(particles[i].p, particles[i].pdg, matID, particles[i].rs, aParticleTable, aMaterialTable);
    Four<FloatingType> four=infs.GetFS(particles[i].p, particles[i].pdg, matID, gen,//particles[i].rs,
                                       particles[i].tr, aParticleTable, aMaterialTable);

    const double costcm=particles[i].tr;
    const int bin=(costcm+1.0)/deltacoscm;

#ifdef OPENACC
#pragma acc atomic update
#endif
    ++HistogramNiCosCM[bin];
    
    if(bin>=0 && bin<NpointsInCosCM-1)
    {
      ++InHistogramRange;
    }
    else
    {
      ++OutOfHistogramRange;
    }
  }

  std::cout<<"InHistogramRange="<<InHistogramRange
           <<" OutOfHistogramRange="<<OutOfHistogramRange<<std::endl;

*/

  RNDGenerator generator(1);
  for(int i=0; i<NPARTICLES; ++i)
  {
    if(i%10000000==0) std::cout<<"i="<<i<<std::endl;
    /////////Four<FloatingType> four=infs.GetFS(InitP, deuteronPDG, matID, generator, aParticleTable, aMaterialTable);
    ///////Four<FloatingType> four=infs.GetFS(particles[i].p, particles[i].pdg, matID, particles[i].rs, aParticleTable, aMaterialTable);
    FloatingType tr=0.0;
    Four<FloatingType> four=infs.GetFS(InitP, deuteronPDG, matID, gen,
                                       tr, aParticleTable, aMaterialTable);

    const double costcm=tr;
    const int bin=(costcm+1.0)/deltacoscm;

#ifdef OPENACC
#pragma acc atomic update
#endif
    ++HistogramNiCosCM[bin];
    
    if(bin>=0 && bin<NpointsInCosCM-1)
    {
      ++InHistogramRange;
    }
    else
    {
      ++OutOfHistogramRange;
    }
  }

  std::cout<<"InHistogramRange="<<InHistogramRange
           <<" OutOfHistogramRange="<<OutOfHistogramRange<<std::endl;


  
  FillHistogramNi();
  Histogram_DifferentiatedPartialSums();


  

  auto end=std::chrono::steady_clock::now();
	auto elapsed_ms=std::chrono::duration_cast<std::chrono::milliseconds>(end-begin);
	std::cout<<"time="<<elapsed_ms.count()<<" ms"<<std::endl;
  std::cout<<"Nbin="<<Nbin<<std::endl;
}
