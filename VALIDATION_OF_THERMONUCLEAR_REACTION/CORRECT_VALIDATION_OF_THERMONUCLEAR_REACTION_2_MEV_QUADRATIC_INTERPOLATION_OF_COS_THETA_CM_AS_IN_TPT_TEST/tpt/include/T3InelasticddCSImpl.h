#pragma once
#ifndef T3INELASTICDDCSIMPL_H
#define T3INELASTICDDCSIMPL_H

#include <cmath>
#include "T3ParticleTable.h"
#include "T3MaterialTable.h"

namespace t3 {
  
template <typename Floating>
class InelasticddCS {
public:
  ////\\\\/////InelasticddCS() = default;
  InelasticddCS()
  {
    ParticleTable aParticleTable;
    deuteronPDG=aParticleTable.makePDGfromZandA(1,2);
  }
// PDG_t - PDG code of incident particle;
// MatID_t - target material number.
  inline Floating
  GetCS(Floating e, PDG_t incPDG, MatID_t matID,
        ParticleTable & aParticleTable, MaterialTable & aMaterialTable) const;
private:
  PDG_t deuteronPDG;
  //Gamov energy for d-d = 0.986 MeV. =2*mr*c^2*(pi*alpha*Zp*Zt)^2. mr=md/2.Zp-charge of projectile. Zt-charge of target.
  const Floating Eg = 0.986 * MeV;
};

//Formulas for integral cross section:
//x - Tcm in MeV:
//n+He3:
//funnhe3=0.21/(1.+(2.52e-2/x)**1.5)/(1.+(4.7e-3/x)**4.5)/(1.+0.23*x) * g
//p+t:
//funpt=0.176/(1.+(2.0e-2/x)**1.7)/(1.+(4.3e-3/x)**4.4)/(1.+(0.17*x)**0.85) * g
//Here g=EXP(-SQRT(Eg/Tcm)/2) is the Gamow factor.
//
  
//Tls - LS kin en of incident particle.
//incPDG - PDG of the incident particle (deuteron).
//matID  - id of the material. 
template <typename Floating>
Floating InelasticddCS<Floating>::GetCS(Floating Tls, PDG_t incPDG, MatID_t matID,
                                        ParticleTable & aParticleTable,
                                        MaterialTable & aMaterialTable) const {
  ////\\\\////static const auto aParticleTable = ParticleTable();//create ParticleTable.
  ////\\\\////static const auto aMaterialTable = MaterialTable();//create MaterialTable.
  ////\\\\////const PDG_t deuteronPDG = aParticleTable.makePDGfromZandA(1,2);
  // Here we check if the incident particle is d. If it is not d, the inelastic cs is 0.0.
  if (incPDG != deuteronPDG) return 0.;
  // WE HAVE ONLY 1 MATERIAL - Tid2. HERE WE ITERATE THROUGH ITS ISOTOPES:
  // Ti48 and d. When we find a d isotope, we calculate the total cross section
  // of 2 inelastic channels at the energy e and return it.

  std::cout<<"Tls="<<Tls<<" incPDG="<<incPDG<<" matID="<<matID
           <<" deuteronPDG="<<deuteronPDG<<std::endl;

  std::cout<<"matID="<<matID<<std::endl;
  for(int i=0; i<aMaterialTable.GetNumberOfIsotopes(matID); ++i)
    std::cout<<aMaterialTable.GetIsotopes(matID, i)<<" ";
  std::cout<<std::endl;
  
  for(int i=0; i<aMaterialTable.GetNumberOfIsotopes(matID); ++i)
  {
    PDG_t isotope=aMaterialTable.GetIsotopes(matID, i);

    std::cout<<"isotope="<<isotope<<std::endl;
    
    //If D-D:
    if (isotope == deuteronPDG)
    {
      //CM kinetic energy of the incident particle (projectile) for d-d.
      const Floating Tcm = Tls / 2;
//Gamov factor for d-d:
      const Floating rat = std::sqrt(Eg / Tcm);
      //Gamov factor for d-d.
      const Floating g = std::exp(-0.5*rat);
      //integral cross section for inelastic dd channel n+he3:
      const Floating nhe31 = 1.0 + std::pow(2.52e-2/Tcm, 1.5);
      const Floating nhe32 = 1.0 + std::pow(4.7e-3/Tcm, 4.5);
      //* barn, because the formula gives ths cs in barns.
      const Floating csnhe3 = 0.21 / nhe31 / nhe32 / (1.+0.23*Tcm) * g * barn;
      //integral cross section for inelastic dd channel p+t cs:
      const Floating pt1 = 1.0 + std::pow(2.0e-2/Tcm, 1.7);
      const Floating pt2 = 1.0 + std::pow(4.3e-3/Tcm, 4.4);
      const Floating pt3 = 1.0 + std::pow(0.17*Tcm, 0.85);
      //* barn, because the formula gives ths cs in barns.
      const Floating cspt = 0.176/ pt1 / pt2 / pt3 * g * barn;
      const Floating cs = csnhe3 + cspt; // total cs=cs(n+He3)+cs(p+t)

      std::cout<<"Tcm="<<Tcm<<" rat="<<rat<<" g="<<g
               <<" csnhe3="<<csnhe3/barn<<" cspt="<<cspt/barn
               <<" cs="<<cs/barn<<std::endl;
      
      return cs;//returns cs in barns.
    }
  }
  return 0.;
}

}//namespace t3
#endif//T3INELASTICDDCSIMPL_H
