#pragma once
#ifndef T3DATAHOLDER_H
#define T3DATAHOLDER_H

#include <cmath>
#include "T3Globals.h"
#include "T3RNG.h"
#include "T3ParticleTable.h"
#include "T3LorentzVector.h"

namespace t3 {


  
//PUSH INJ particles:
template <typename Floating>
void InitParticle(int lfi/*LIFE+i*/, int mei/*MAX_ELEMENT+i*/, Floating Tls, ParticleTable & aParticleTable)
{
  t3::PDG_t initPDG = /*PDG_t(22);*/aParticleTable.makePDGfromZandA(1, 2);
  const Floating m = aParticleTable.GetMass(initPDG);
  const Floating InitEls= m+Tls;
  const Floating pls=sqrt(InitEls*InitEls-m*m);
  particles[lfi] = Particle<Floating>(
       t3::LorentzVector<Floating>(InitParticlex0*ag,
                                   InitParticley0*ag,
                                   -cuba*ag,0.0),
       t3::LorentzVector<Floating>(0.0,0.0,pls,InitEls),
       0.0, initPDG, 1., mei, 1, mei-1);  
}//End of InitParticle.



  
}//namespace t3.

#endif//T3DATAHOLDER_H
