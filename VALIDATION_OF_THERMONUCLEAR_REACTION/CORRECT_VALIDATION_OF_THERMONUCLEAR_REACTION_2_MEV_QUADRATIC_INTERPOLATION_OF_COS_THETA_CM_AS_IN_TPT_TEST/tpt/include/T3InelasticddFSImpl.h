#pragma once
#ifndef T3INELASTICDDFSIMPL_H
#define T3INELASTICDDFSIMPL_H

#include <random>
#include "T3Defs.h"
#include "T3ThreeVector.h"
#include "T3LorentzVector.h"
#include "T3ParticleTable.h"
#include "T3MaterialTable.h"
#include "T3InelasticddCSImpl.h"

#include "T3Inelasticdd_DB.hh"
#include "T3NSGangular_RW.hh"
#include "T3NSGangular_node.hh"

#include <typeinfo>

//#define DEBUG

namespace t3 {

template <typename Floating = double>
class InelasticddFS {
public:
  InelasticddFS()
  {
//********************************************************************************//
//Initialize constants:
//********************************************************************************//    
    ParticleTable aParticleTable;
    //PDG codes of particles at page 314 of big PDG book:
    neutronPDG  = PDG_t(2112);
    protonPDG   = PDG_t(2212);
    deuteronPDG = aParticleTable.makePDGfromZandA(1, 2);
    tritonPDG   = aParticleTable.makePDGfromZandA(1, 3);
    he3PDG      = aParticleTable.makePDGfromZandA(2, 3);
    mn          = aParticleTable.GetMass(neutronPDG);  //neutron mass.
    mp          = aParticleTable.GetMass(protonPDG);   //proton mass.
    md          = aParticleTable.GetMass(deuteronPDG); //deuteron mass.
    mt          = aParticleTable.GetMass(tritonPDG);   //triton mass.
    mhe3        = aParticleTable.GetMass(he3PDG);      //3He2 mass.
    md2         = md * md;
////For the formula 38.16 from the big PDG book 
    //n+he3:
    D2Plusn_he3  = (mn + mhe3)*(mn + mhe3);
    D2Minusn_he3 = (mn - mhe3)*(mn - mhe3);
    //p+t:
    D2Plusp_t    = (mp + mt)*(mp + mt);
    D2Minusp_t   = (mp - mt)*(mp - mt);
//********************************************************************************//
//Initialize partial sums:
//********************************************************************************//      
    T3NSGangular_RW rw;  
    constexpr int tgZ = 1;//target deuteron Z=1
    constexpr int tgA = 2;//target deuteron A=2
    constexpr auto sPDG = 2112;//reaction product=n+3He//this is a neutron from a reaction product:
    constexpr int incZA = 1002;//1000*Z+A//inc deuteron
    constexpr auto rid = "50";//MT index of reaction in ENDF
    rw.load_binary(tgZ, tgA, rid, sPDG, incZA);//load the data from the binary file:
    T3NSGangular_RWrecord rwrec=rw.at(0);
//********************************************************************************//
//Fill in T3Inelasticdd_DB:
//********************************************************************************//
    //make square of partial sums:
    for(int i=0; i<rwrec.size(); ++i)
    {
      for(int j=0; j<127; ++j) rwrec.at(i).Set_V(j+1, rwrec.at(i).Get_V(j+1)*rwrec.at(i).Get_V(j+1));
    }
    //T3Inelasticdd_DB db;
    const Floating Emin = rwrec.front().Get_E();
    const Floating Emax = rwrec.back().Get_E();
    const Floating Ediff = Emax - Emin;
    const size_t nbins = 512 - 1;
    const Floating dE = Ediff / nbins;
    for(size_t i = 0; i < 512; ++i)
    {
      const Floating Ei = Emin + i*dE;
      const T3NSGangular_RWnode rwnode = rwrec.interpolate(Ei);
      T3NSGangular_node nodes_i = T3NSGangular_node(rwnode);
      for(int j=0; j<127; ++j)
      {
        db.Set_V(i, j, nodes_i.Get_V(j));
        db.Set_a(i, j, nodes_i.Get_a(j));
        db.Set_b(i, j, nodes_i.Get_b(j));
        db.Set_c(i, j, nodes_i.Get_c(j));
      }
      const Floating pri=nodes_i.Get_pr();
      const Floating sli=nodes_i.Get_sl();
      db.Set_pr(i, pri);
      db.Set_sl(i, sli);
      db.Set_Einc(i, Ei);
    }
//********************************************************************************//
//End of fill in T3Inelasticdd_DB.
//********************************************************************************//
#ifdef DEBUG
    std::cout<<"Check T3Inelasticdd_DB db:"<<std::endl;
    std::cout<<"V:"<<std::endl;
    for(int i=0; i<512; ++i)
    {
      std::cout<<"node #"<<i<<":"<<std::endl;
      std::cout<<"Einc="<<db.Get_Einc(i)<<" pr="<<db.Get_pr(i)<<" sl="<<db.Get_sl(i)<<std::endl;
      /*
      for(int j=0; j<127; ++j) std::cout<<db.Get_V(i, j)<<" ";
      std::cout<<std::endl;
      */
      //for(int j=0; j<127; ++j) std::cout<<db.Get_a(i, j)<<" ";
      //std::cout<<std::endl;
      //for(int j=0; j<127; ++j) std::cout<<db.Get_b(i, j)<<" ";
      //std::cout<<std::endl;
      for(int j=0; j<127; ++j) std::cout<<db.Get_c(i, j)<<" ";
      std::cout<<std::endl;
    }

    
    /*
    std::cout<<"Einc:"<<std::endl;
    for(size_t i = 0; i < 512; ++i) std::cout<<db.Get_Einc(i)<<" ";
    std::cout<<std::endl;
    std::cout<<"pr:"<<std::endl;
    for(size_t i = 0; i < 512; ++i) std::cout<<db.Get_pr(i)<<" ";
    std::cout<<std::endl;
    std::cout<<"sl:"<<std::endl;
    for(size_t i = 0; i < 512; ++i) std::cout<<db.Get_sl(i)<<" ";
    std::cout<<std::endl;
    */

    /*
    std::cout<<"E:"<<std::endl;
    for(size_t i = 0; i < rwrec.size(); ++i) std::cout<<rwrec.at(i).Get_E()<<" ";
    std::cout<<std::endl;
    std::cout<<"pr:"<<std::endl;
    for(size_t i = 0; i < rwrec.size(); ++i) std::cout<<rwrec.at(i).Get_pr()<<" ";
    std::cout<<std::endl;
    std::cout<<"sl:"<<std::endl;
    for(size_t i = 0; i < rwrec.size(); ++i) std::cout<<rwrec.at(i).Get_sl()<<" ";
    std::cout<<std::endl;
    */

    
    /*
    std::cout<<"pr:"<<std::endl;
    for(size_t i = 0; i < 512; ++i) std::cout<<db.Get_pr(i)<<" ";
    std::cout<<std::endl;
    std::cout<<"sl:"<<std::endl;
    for(size_t i = 0; i < 512; ++i) std::cout<<db.Get_sl(i)<<" ";
    std::cout<<std::endl;
    */


    std::cout<<"Check in constructor:"<<std::endl;
    std::cout<<"neutronPDG="<<neutronPDG<<" protonPDG="<<protonPDG
             <<" deuteronPDG="<<deuteronPDG<<" tritonPDG="<<tritonPDG
             <<" he3PDG="<<he3PDG<<std::endl;

    std::cout<<"mn="<<mn<<" mp="<<mp<<" md="<<md
             <<" mt="<<mt<<" mhe3="<<mhe3
             <<" md2="<<md2<<std::endl;

    std::cout<<"D2Plusn_he3="<<D2Plusn_he3<<" D2Minusn_he3="<<D2Minusn_he3
             <<" D2Plusp_t="<<D2Plusp_t<<" D2Minusp_t="<<D2Minusp_t
             <<std::endl;

#endif
    
//Checked all these values, which are filled in the constructor.
//I checked them, and they are right.
  }
//p - LS 4-momentum-vector of the incident particle(d).
//incPDG - PDG of the incident particle.
//matID  - id of the material.
//engine - random number generator of the incident particle.
  template <typename RandomEngine>
  inline auto GetFS(LorentzVector<Floating> const &p, PDG_t incPDG, MatID_t matID, RandomEngine &engine,
                    Floating & tr, ParticleTable & aParticleTable, MaterialTable & aMaterialTable) const;
private:
//PDG:  
  PDG_t neutronPDG;
  PDG_t protonPDG;
  PDG_t deuteronPDG;
  PDG_t tritonPDG;
  PDG_t he3PDG;
//Gamow energy for D-D:  
  const Floating Eg = 0.986 * MeV;
//masses:  
  Floating mn;
  Floating mp;
  Floating md;
  Floating mt;
  Floating mhe3;
  Floating md2;
//(m1+m2)^2 and (m1-m2)^2 for the formula (36.18) on page 321 of the big PDG book:
  Floating D2Plusn_he3;
  Floating D2Minusn_he3;
  Floating D2Plusp_t;
  Floating D2Minusp_t;
  //
  T3Inelasticdd_DB db;
};
//returns PDG of inc d, 4-momentum of inc d, PDG of target d, 4-momentum of target d
template <typename Floating>
template <typename RandomEngine>
auto InelasticddFS<Floating>::GetFS(LorentzVector<Floating> const &p, PDG_t incPDG, MatID_t matID, RandomEngine &engine,
                                    Floating & tr, ParticleTable & aParticleTable, MaterialTable & aMaterialTable) const {
    auto const generateSubCanonical = [&engine]() {
      return GenerateSubCanonical<Floating>(engine);
  };
  //This method is only for d-d inelastic scattering. If inc particle is not d, exit.
  if(incPDG != deuteronPDG)
    return Four<Floating>(PDG_t(0), PDG_t(0), LorentzVector<Floating>(0.0,0.0,0.0,0.0), LorentzVector<Floating>(0.0,0.0,0.0,0.0));
  const Floating m = aParticleTable.GetMass(incPDG);  //incident deuteron mass
  const Floating E = p.E();                           //incident deuteron LS full energy
  const Floating Tls = E - m;                         //incident deuteron LS kin energy 
  PDG_t outPDG1=-1;
  PDG_t outPDG2=-1;
  const Floating Einit=E+md;

#ifdef DEBUG  
  std::cout<<"InelasticddFS::GetFS():"<<std::endl;
  std::cout<<"incPDG="<<incPDG<<" matID="<<matID
           <<" p:"<<p.x()<<" "<<p.y()<<" "<<p.z()
           <<" "<<p.E()/MeV<<std::endl;
  std::cout<<"m="<<m/MeV<<" E="<<E<<" Tls="<<Tls<<std::endl;
  std::cout<<"outPDG1="<<outPDG1<<" outPDG2="<<outPDG2<<std::endl;
  std::cout<<"matID="<<matID<<" : n="
           <<aMaterialTable.GetNumberOfIsotopes(matID)<<std::endl;
  for(int j=0; j<aMaterialTable.GetNumberOfIsotopes(matID); ++j)
    std::cout<<aMaterialTable.GetIsotopes(matID, j)<<" ";
  std::cout<<std::endl;
#endif  
  //iterate through all isotopes and try to find deuteron:
  for(int i=0; i<aMaterialTable.GetNumberOfIsotopes(matID); ++i)
  {
    const PDG_t isotope=aMaterialTable.GetIsotopes(matID, i);
#ifdef DEBUG    
    std::cout<<"isotope="<<isotope<<std::endl;
#endif    
    if(isotope == deuteronPDG)
    {
#ifdef DEBUG  
      std::cout<<"isotope="<<isotope<<"=="<<deuteronPDG<<std::endl;
#endif      
      //CM kinetic energy of the incident deuteron (for d-d = Tls/2).
      const Floating Tcm = Tls / 2;
      const Floating rat = std::sqrt(Eg / Tcm);
      //Gamov factor for d-d.
      const Floating g = std::exp(-0.5*rat);
//Formulas for integral cross section:
//x - Tcm in MeV:
//n+He3:
//funnhe3=0.21/(1.+(2.52e-2/x)**1.5)/(1.+(4.7e-3/x)**4.5)/(1.+0.23*x) * g
//p+t:
//funpt=0.176/(1.+(2.0e-2/x)**1.7)/(1.+(4.3e-3/x)**4.4)/(1.+(0.17*x)**0.85) * g
//Here g=EXP(-SQRT(Eg/Tcm)/2) is the Gamow factor.
//
      // inelastic dd channel n+he3 cs:
      const Floating nhe31 = 1.0 + pow(2.52e-2/Tcm, 1.5);
      const Floating nhe32 = 1.0 + pow(4.7e-3/Tcm, 4.5);
      //* barn, because the formula gives ths cs in barns.
      const Floating csnhe3 = 0.21 / nhe31 / nhe32 / (1.+0.23*Tcm) * g * barn;
      // inelastic dd channel p+t cs:
      const Floating pt1 = 1.0 + pow(2.0e-2/Tcm, 1.7);
      const Floating pt2 = 1.0 + pow(4.3e-3/Tcm, 4.4);
      const Floating pt3 = 1.0 + pow(0.17*Tcm, 0.85);
      //* barn, because the formula gives ths cs in barns.
      const Floating cspt = 0.176/ pt1 / pt2 / pt3 * g * barn;
      const Floating cstot = csnhe3 + cspt; //total cs=cs(n+He3)+cs(p+t)
      const Floating R = GenerateSubCanonical<Floating>(engine);

#ifdef DEBUG      
      std::cout<<"Tcm="<<Tcm<<" Eg="<<Eg/MeV<<" rat="<<rat<<" g="<<g
               <<" R="<<R<<std::endl;

      std::cout<<"csnhe3="<<csnhe3/barn<<" barn cspt="<<cspt/barn
               <<" barn cstot="<<cstot/barn<<" barn"<<std::endl;

//Checked that integral cross sections are right.

      std::cout<<"R="<<R<<std::endl;
      std::cout<<0.0<<" "<<csnhe3/cstot<<" "<<1.0<<std::endl;

#endif
      
      if(0.0<=R && R<csnhe3/cstot)       outPDG1 = PDG_t(2112);//n+He3
      else if(csnhe3/cstot<=R && R<=1.0) outPDG1 = PDG_t(2212);//p+t
      //here we store the init LS momentum of the incident particle to get the target particle LS momentum from the momentum conservation law.
      ThreeVector<Floating> InitMomentum = p.Vect();
      const Floating plsinc=InitMomentum.R();
      const Floating Elsinc=sqrt(md2+plsinc*plsinc);//incident deuteron LS energy
      const Floating s=2*md2+2*Elsinc*md;//2*md^2+2*Elsinc*md.
      const Floating coef=md/sqrt(s);

#ifdef DEBUG      
      std::cout<<"outPDG1="<<outPDG1<<std::endl;
      std::cout<<"InitMomentum: "<<InitMomentum.x()
               <<" "<<InitMomentum.y()<<" "<<InitMomentum.z()
               <<std::endl;
      std::cout<<"plsinc="<<plsinc<<" m="<<m<<" Elsinc="<<Elsinc<<" s="<<s
               <<" coef="<<coef<<std::endl;
#endif      
      
      //LS->CM:
      ThreeVector<Floating> Pcm = InitMomentum*coef;//pcm - initial incident deuteron CM momentum vector
      //components to form 2 vectors perpendicular to pcm vector.
      const Floating xy=Pcm.x()*Pcm.y(), xz=Pcm.x()*Pcm.z(), yz=Pcm.y()*Pcm.z();
      const Floating x2=Pcm.x()*Pcm.x(), y2=Pcm.y()*Pcm.y(), z2=Pcm.z()*Pcm.z();
      ThreeVector<Floating> e1, e2, e3;//create vectors for e1, e2, e3 orts:
      //choose the minimal component of the inc deuteron CM momentum vector. it is necessary because if, for example, use
      //e1={0., vz, -vy} and e2={y2+z2, -xy, -xz}; at V=(v,0,0), then we will have e1=(0,0,0), e2=(0,0,0) - error.
      if(Pcm.x() < Pcm.y())
      {
        if(Pcm.x() < Pcm.z())
        {
          e2={0., Pcm.z(), -Pcm.y()};
          e3={y2+z2, -xy, -xz};
        }
        else
        {
          e2={Pcm.y(), -Pcm.x(), 0.};
          e3={-xz, -yz, y2+x2};
        }
      }
      else
      {
        if(Pcm.y() < Pcm.z())
        {
          e2={Pcm.z(), 0., -Pcm.x()};
          e3={xy, -x2-z2, yz};
        }
        else
        {
          e2={Pcm.y(), -Pcm.x(), 0.};
          e3={-xz, -yz, y2+x2};
        }
      }
      //e1 is || to Pcm and Vcm.
      e1=Pcm;

#ifdef DEBUG      
      std::cout<<"pcminc="<<Pcm.R()<<std::endl;
      std::cout<<"Pcm: "<<Pcm.x()<<" "<<Pcm.y()<<" "<<Pcm.z()<<std::endl;
      
      std::cout<<"e1: "<<e1.x()<<" "<<e1.y()<<" "<<e1.z()<<std::endl;
      std::cout<<"e2: "<<e2.x()<<" "<<e2.y()<<" "<<e2.z()<<std::endl;
      std::cout<<"e3: "<<e3.x()<<" "<<e3.y()<<" "<<e3.z()<<std::endl;
#endif      
      
      //normalize e1, e2, e3.
      //It is necessary to normalize pcm, because in inelastic reaction
      //the value of pcm changes.
      e1.Unit();
      e2.Unit();
      e3.Unit();

#ifdef DEBUG      
      std::cout<<"e1: "<<e1.x()<<" "<<e1.y()<<" "<<e1.z()<<std::endl;
      std::cout<<"e2: "<<e2.x()<<" "<<e2.y()<<" "<<e2.z()<<std::endl;
      std::cout<<"e3: "<<e3.x()<<" "<<e3.y()<<" "<<e3.z()<<std::endl;      
#endif
      
//*************************//      
//?????????????????????????//      
      //const Floating theta=M_PI*GenerateSubCanonical<Floating>(engine);//random theta from 0 to pi (0 to pi/2 ???)
//THIS IS cos(theta_cm), NOT theta_cm:!!!      
      const Floating cosTheta=db.RandomizeCost(engine, Tls);//random theta from 0 to pi (0 to pi/2 ???)
      tr=cosTheta;
//?????????????????????????//
//*************************//      
      
      const Floating phi=2*M_PI*GenerateSubCanonical<Floating>(engine);//random phi from 0 to 2*pi
      const Floating cosPhi=cos(phi), sinPhi=sin(phi);
      const Floating sinTheta=sqrt(1.0-cosTheta*cosTheta);

#ifdef DEBUG      
      std::cout<<"phi="<<phi<<" theta="<<theta<<std::endl;
      std::cout<<"cosPhi="<<cosPhi<<" sinPhi="<<sinPhi
               <<" cosTheta="<<cosTheta<<" sinTheta="<<sinTheta
               <<std::endl;

#endif      
      
      
//*****************************************************************//
//Here the value of pcm changes, because the reaction is inelastic:      
//*****************************************************************//
//Use the formula (36.18) on page 321 from the big PDG book:
      const Floating pcminc=Pcm.R();//initial module of the momentum vector of the incident deuteron in CM
      const Floating pcminc2=pcminc*pcminc;
      const Floating pcminc4=pcminc2*pcminc2;
      Floating num=-1.0;
      Floating prodm1=-1.0;//the mass of the product particle 1 - n or p.
      Floating prodm2=-1.0;//the mass of the product particle 2 - he3 or t.
      if(outPDG1==PDG_t(2112)/*neutron*/)//n+he3:
      {
        num = ( s - D2Plusn_he3 ) * ( s - D2Minusn_he3 );
        outPDG2=he3PDG;
        prodm1=mn;
        prodm2=mhe3;
      }
      else if(outPDG1==PDG_t(2212)/*proton*/)//p+t:
      {
        num = ( s - D2Plusp_t ) * ( s - D2Minusp_t );
        outPDG2=tritonPDG;
        prodm1=mp;
        prodm2=mt;
      }
      const Floating pcmfin=sqrt(num)/2/sqrt(s);      
//End of use the formula from the big PDG book.      
//*************************************************************************//
//End of where the value of pcm changes, because the reaction is inelastic.
//*************************************************************************//

#ifdef DEBUG
      std::cout<<"pcminc="<<pcminc<<" pcminc2="<<pcminc2
               <<" pcminc4="<<pcminc4<<std::endl;
      std::cout<<"num="<<num<<" m1="<<prodm1<<" m2="<<prodm2
               <<std::endl;
      std::cout<<"outPDG1="<<outPDG1<<" outPDG2="<<outPDG2<<std::endl;
      std::cout<<"pcmfin="<<pcmfin<<std::endl;
#endif

      
      const Floating ss=pcmfin*sinTheta*sinPhi, sc=pcmfin*sinTheta*cosPhi;
      ThreeVector<Floating> Pcmfin;//pcmfin - final incident deuteron CM momentum 3-vector:
      const Floating pcmfinCosTheta=pcmfin*cosTheta;
      Pcmfin.SetX(e1.x()*pcmfinCosTheta+e2.x()*ss+e3.x()*sc);
      Pcmfin.SetY(e1.y()*pcmfinCosTheta+e2.y()*ss+e3.y()*sc);
      Pcmfin.SetZ(e1.z()*pcmfinCosTheta+e2.z()*ss+e3.z()*sc);
      //transform momentum of the inc deuteron from CM to LS:
      //center of mass velocity does not change, because there are no external forces,
      //so we use this center of mass velocity:
      ThreeVector<Floating> Vcm=InitMomentum/(Elsinc+md);
      const Floating Ecminc=sqrt(md2+pcminc*pcminc);
      const Floating E1cmfin=sqrt(prodm1*prodm1+pcmfin*pcmfin);//incident deuteron CM full energy
      const Floating E2cmfin=sqrt(prodm2*prodm2+pcmfin*pcmfin);//incident deuteron CM full energy
      const Floating AbsVcm=Vcm.R();//module of the CM velocity
      const Floating gamma=1.0/std::sqrt(1-AbsVcm*AbsVcm);//gamma-factor
      //Lorentz transformation is only for the parallel (parallel to Vcm)
      //component of the momentum vector. The perpendicular component
      //does not change.
      //So, here i should calculate \\ and perpendicular components
      //of the Pcmfin and change according to the Lorentz transformation
      //only the || to Vcm component and leave the perpendicular component
      //unchanged.

#ifdef DEBUG
      std::cout<<"ss="<<ss<<" sc="<<sc<<" Pcmfin: "
               <<Pcmfin.x()<<" "<<Pcmfin.y()<<" "<<Pcmfin.z()
               <<std::endl;
      std::cout<<"Vcm: "<<Vcm.x()<<" "<<Vcm.y()<<" "<<Vcm.z()
               <<std::endl;

      std::cout<<"Ecminc="<<Ecminc<<" E1cmfin="<<E1cmfin
               <<" E2cmfin="<<E2cmfin<<" AbsVcm="<<AbsVcm
               <<" gamma="<<gamma<<std::endl;
#endif      

      
//Calculate the parallel and perpendicular components of Pcmfin:
      //parallel:
      ThreeVector<Floating> ParallelProjection=(e1*Pcmfin)*e1;//???? (Vcm.Unit()*Pcmfin).
      //perpendicular:
      ThreeVector<Floating> PerpendicularProjection=Pcmfin-ParallelProjection;
      ThreeVector<Floating> plsfin1=gamma*(ParallelProjection+Vcm*E1cmfin);//final LS momentum vector of the 1 product particle.
      plsfin1+=PerpendicularProjection;
      ThreeVector<Floating> plsfin2=InitMomentum-plsfin1;     //final LS momentum vector of the 2 product particle.
      //calculate the energies of 1,2 deuterons
      const Floating Absplsfin1=plsfin1.R();
      const Floating Elsfin1=std::sqrt(prodm1*prodm1+Absplsfin1*Absplsfin1);//=sqrt(m1^2+p1^2)
      const Floating Absplsfin2=plsfin2.R();
      const Floating Elsfin2=std::sqrt(prodm2*prodm2+Absplsfin2*Absplsfin2);//=sqrt(m2^2+p2^2)

      const Floating Efin=Elsfin1+Elsfin2;

#ifdef DEBUG    
      std::cout<<"ParallelProjection: "
               <<ParallelProjection.x()<<" "
               <<ParallelProjection.y()<<" "
               <<ParallelProjection.z()
               <<std::endl;

      std::cout<<"PerpendicularProjection: "
               <<PerpendicularProjection.x()<<" "
               <<PerpendicularProjection.y()<<" "
               <<PerpendicularProjection.z()
               <<std::endl;

      std::cout<<"plsfin1="<<plsfin1<<" plsfin2="<<plsfin2
               <<" Absplsfin1="<<Absplsfin1<<" Absplsfin2="<<Absplsfin2
               <<" Elsfin1="<<Elsfin1<<" Elsfin2="<<Elsfin2<<std::endl;
               
      
      //std::cout.precision(5);
      //std::cout.setf(std::ios::fixed);
      std::cout<<"Einit="<<Einit<<" Efin="<<Efin
               <<" Einit-Efin="<<Einit-Efin<<std::endl;
#endif      
      
      LorentzVector<Floating> outP1 = LorentzVector<Floating>(plsfin1.x(), plsfin1.y(), plsfin1.z(), Elsfin1);
      LorentzVector<Floating> outP2 = LorentzVector<Floating>(plsfin2.x(), plsfin2.y(), plsfin2.z(), Elsfin2);

      const LorentzVector<Floating> tot = outP1 + outP2;
#ifdef DEBUG
      std::cout<<"outP1: "<<outP1.x()<<" "<<outP1.y()<<" "
               <<outP1.z()<<" "<<outP1.E()<<std::endl;
      std::cout<<"outP2: "<<outP2.x()<<" "<<outP2.y()<<" "
               <<outP2.z()<<" "<<outP2.E()<<std::endl;
      std::cout<<"tot: "<<tot.x()<<" "<<tot.y()<<" "
               <<tot.z()<<" "<<tot.E()<<std::endl;
#endif
      
      return Four(outPDG1, outPDG2, outP1, outP2);
    }
  }
}
  
}//namespace t3.
#endif // T3INELASTICDDFSIMPL_H
