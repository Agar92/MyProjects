#include "T3NSGangular_record.hh"
//#include "Randomize.hh"

#include <iostream>
#include "unistd.h"

//#define debug

T3NSGangular_record::T3NSGangular_record(const T3NSGangular_RWrecord& rwrec)
{
  std::cout<<"45678"<<std::endl;
  
  const size_t rwrecsize = rwrec.size();
  if(rwrecsize < 2) T3cout << __func__ << ": Warning: bad rw record" << rwrec << T3endl;
  const T3double Emin = rwrec.front().Get_E();
  const T3double Emax = rwrec.back().Get_E();
  const T3double Ediff = Emax - Emin;
  if (Ediff < 1e-7)
  {
    T3cout << __func__ << " Warning: too small energy interval " << Ediff << " MeV"
           << T3endl;
    for(size_t ind = 0; ind < rwrecsize; ++ind )
    {
      T3cout << rwrec.at(ind).Get_E() << ' ';
    }
    T3cout << T3endl;
  }
  static const size_t nbins = _num_nodes - 1;
  const T3double dE = Ediff / nbins;
  for(size_t i = 0; i < _num_nodes; ++i)
  {
    const T3double Ei = Emin + i*dE;
    Einc[i] = Ei;

    //std::cout<<"BLUE"<<std::endl;
    //sleep(1);
    
    const T3NSGangular_RWnode rwnode = rwrec.interpolate(Ei);

    //std::cout<<"RED"<<std::endl;
    
    nodes[i] = T3NSGangular_node(rwnode);

//\\//***************************************//
    //std::cout<<"Ei="<<Einc[i]<<":"<<std::endl;
    //std::cout<<nodes[i]<<std::endl;
    //usleep(1000000);
//\\//***************************************//    
    
  }
}

//\\//G4double T2NSGangular_record::RandomizeCost(G4double E)
//\\//***************************************//
T3double T3NSGangular_record::RandomizeCost(RNDGenerator & generator, T3double E)
//\\//***************************************//
{
  const T3double Emin = Einc[0];
  const T3double Emax = Einc[_num_nodes - 1];
  const T3double Ediff = Emax - Emin;
  static const size_t nbins = _num_nodes - 1;
  const T3double dE = Ediff / nbins;
  if(E <= Emin || E >= Emax)
  {
    // isotropic
    //\\//const G4double R = G4UniformRand();
//\\//*****************************************//
    const T3double R = t3::GenerateSubCanonical<FloatingType>(generator);
//\\//*****************************************//
    return -1 + 2*R;
  }
  else
  {
    size_t const binn = static_cast<size_t>((E - Emin)*(_num_nodes-1)/Ediff);
    const T3double Elow  = Einc[binn];
//     const G4double Ehigh = Einc[binn + 1];
    const T3double W = (E - Elow)/(dE);
    //\\//const G4double R = G4UniformRand();
//\\//****************************************//
    const T3double R = t3::GenerateSubCanonical<FloatingType>(generator);
//\\//****************************************//    
    if (R < W)
    {
      // randomize by high distribution
#ifdef debug
      T3cout << "Won high bin E =" << Einc[binn + 1] << T3endl;
#endif
      //\\//return nodes[binn + 1].RandomizeCost();
//\\//*****************************************//
      return nodes[binn + 1].RandomizeCost(generator);
//\\//*****************************************//            
    }
    else
    {
      // randomize by low distribution
#ifdef debug
      T3cout << "Won low bin E =" << Elow << T3endl;
#endif
      //\\//return nodes[binn].RandomizeCost();
//\\//******************************************//
      return nodes[binn].RandomizeCost(generator);
//\\//******************************************//      
    }
  }
}

