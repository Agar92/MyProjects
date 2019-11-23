#include "T3NSGangular_node.hh"
//#include "Randomize.hh"
#include "T3Utility.hh"

#include "unistd.h"

// TODO Get/set with 1-based enumeration
T3NSGangular_node::T3NSGangular_node(const T3NSGangular_RWnode& rwnode):
                                     pr(rwnode.Get_pr()), sl(rwnode.Get_sl())
{
  //\\//std::vector<G4LorentzVector> vlv = rwnode.CalcXYZ();
//\\//***************************************//
  std::vector<LorentzVector<FloatingType>> vlv = rwnode.CalcXYZ();
//\\//***************************************//

//\\//***************************************//
  //std::cout<<rwnode<<std::endl;
//\\//***************************************//    
  
  for(size_t ind =0; ind < _num_point; ++ind)
  {
    V[ind] = vlv.at(ind).t();
    a[ind] = vlv.at(ind).x();
    b[ind] = vlv.at(ind).y();
    c[ind] = vlv.at(ind).z();
  }

  
//\\//***************************************//
/*  
  std::cout<<"V:"<<std::endl;
  for(int ind=0; ind<_num_point; ++ind) std::cout<<V[ind]<<" ";
  std::cout<<std::endl;
  std::cout<<"a:"<<std::endl;
  for(int ind=0; ind<_num_point; ++ind) std::cout<<a[ind]<<" ";
  std::cout<<std::endl;
  std::cout<<"b:"<<std::endl;
  for(int ind=0; ind<_num_point; ++ind) std::cout<<b[ind]<<" ";
  std::cout<<std::endl;
  std::cout<<"c:"<<std::endl;
  for(int ind=0; ind<_num_point; ++ind) std::cout<<c[ind]<<" ";
  std::cout<<std::endl;
  //usleep(100000);
*/  
//\\//***************************************//      
  
  
}

// copypasted from T2NElasticCS
//\\//G4double T2NSGangular_node::RandomizeCost()
T3double T3NSGangular_node::RandomizeCost(RNDGenerator & generator)
{
  //\\//G4double P = G4UniformRand();
//\\//********************************//  
  T3double P = t3::GenerateSubCanonical<FloatingType>(generator);
//\\//********************************//  
  
#ifdef debug
  T3cout << "T3NSGangular_node::RandomizeCost pr=" << pr << ", sl=" << sl << T3endl;
#endif
  if(P < pr)
  {
    // use exponent
    //\\//return 1. + sl * log(1. - ( 1. - exp(-2 / sl) ) * G4UniformRand() );
    return 1. + sl * log(1. - ( 1. - exp(-2 / sl) ) *
                         t3::GenerateSubCanonical<FloatingType>(generator) );
  }
  //else use tabulated
  static const T3double dY = 2. / (_num_point + 1); // Y step
  //\\//P = sqrt(G4UniformRand());
//\\//********************************//
  P = /*sqrt(*/t3::GenerateSubCanonical<FloatingType>(generator)/*)*/;
//\\//********************************//  
  ///\\\///const size_t hct = T3Utility::bin_search<_num_point>(V, P) - V;
  const size_t hct = T3Utility::bin_search<_num_point>(V, P) - V;
  T3int const ct1 = hct ? hct-1 : hct;                   // 0-bin default for ct1 = hct = 0
  T3int const ct2 = (hct && hct < _num_point) ? hct : -1;// Doesn't have any meaning for 0
  T3double const r1 = hct ? V[ct1] : 0.;    // 0-bin default for the Left End
  T3double const r2 = (hct < _num_point) ? V[hct] : 1.;   // the Right End
//   const G4LorentzVector cf = LVP[ct1];
  if( (P < r1 && hct) || P > r2 + .000000001 || r1 >= r2)
    T3cout << "Warning " << __func__ << ":r1=" << r1 << " < R=" << P << " < r2=" << r2
            << ", h=" << hct << T3endl;
  long double cost = a[ct1] + P*( b[ct1] + P * c[ct1] ); // X+Y*R+Z*R^2 from r1
  const long double dd = dY / ( r2 - r1);                 // the 1D slope in the bin
  if(ct2 >= 0)
  {
//     const G4LorentzVector cs=LVP[ct2];                 // Parameters of f2 (r2 function)
    const T3double dr1 = r1 + r1;
    const T3double dr2 = r2 + r2;
    const long double p12 = b[ct1] + dr2 * c[ct1];        // Derivitive of ** f1 ** in r2
    const long double p21 = b[ct2] + dr1 * c[ct2];        // Derivitive of ** f2 ** in r1
    if(p12 < 0. || p21 < 0.)                        // spline in the reduced interval
    {
      long double F  =  a[ct2] + P*( b[ct2] + P * c[ct2] ); // Original f2 curve (no cor)
      const long double fs = dY * hct - 1.;               // LinFunct = fs + x * dd
      const long double d1 = P - r1;
      const long double d2 = r2 - P;
      if( p12 < 0. )                                // 1st deriv < 0 in r2 => negPart->lin
      {
        T3double z = - ( b[ct1] + r2 * 2 * c[ct1] ) / dd;
#ifdef debug
        T3cout << __func__ << ": p1=" << b[ct1] + P * 2 * c[ct1] << ", l=" << d1
                << ", L=" << fs +( P - r1 ) * dd << ", f=" << d2 << ", F=" << cost
                <<", d="<< d1+d2 << ", pl=" << dd << ", pf=" << b[ct1] + r2 * 2 * c[ct1]
                << T3endl;
#endif
        cost = ( d1 * z * ( fs + ( P - r1 ) * dd ) + d2 * cost ) / ( d1 * z + d2 );
      }
      if( p21 < 0. )                                // 2nd deriv < 0 in r1 => negPart->lin
      {
        const long double z = - ( b[ct2] + r1 * 2 * c[ct2] ) / dd;
#ifdef debug
        T3cout << __func__ << ": p2=" << b[ct2] + P * 2 * c[ct2] << ", l=" << d2
                << ", L=" << fs + (P-r1) * dd << ", f=" << d1 << ", F=" << F <<", d="
                << d1+d2 <<", sL=" << d2 * ( fs + (P-r1) * dd ) / (d1+d2) <<", sF="
                <<  d1 * F / (d1 + d2) << ", pl=" << dd << ", pf="
                << b[ct2] + r1 * 2 * c[ct2] << T3endl;
#endif
        F = (d2 * z * ( fs + ( P - r1 ) * dd ) + d1 * F ) / ( d1 + d2 * z );
      }
#ifdef debug
      T3cout << __func__ << ": "<< d2 <<" * cost="<< cost <<" + "<< d1 <<" * F="
              << F << T3endl;
#endif
      cost = ( d1 * F + d2 * cost ) / ( d1 + d2);
    }
    else                                                  //for the right end.
    {
      const long double p11 = b[ct1] + dr1 * c[ct1];      // Derivitive of ** f1 ** in r1
      const long double p22 = b[ct2] + dr2 * c[ct2];      // Derivitive of ** f2 ** in r2
      long double d1 = P - r1;
      long double d2 = r2 - P;
      long double dr = r2 - r1;
      if      (p11 > 0 && p21 > 0 && p11 > p21)     // Left edge derivitive problem
      {
        dr /= p11 / p21;
        if(d1 < dr)
        {
#ifdef debug
          T3cout << __func__ << ":L r="<< p11/p21 << T3endl;
#endif
          d2 = r1 + dr - P;
          cost = ( d1 * ( a[ct2] + P * ( b[ct2] + P * c[ct2] ) ) + d2 * cost ) / dr;
        }
        else cost = a[ct2] + P * ( b[ct2] + P * c[ct2] ); // Rirht edge approximation
      }
      else if (p12 > 0 && p21 > 0 && p22 > p12)     // Right edge derivitive problem
      {
        dr /= p22 / p12;
        if(d2 < dr)
        {
#ifdef debug
          T3cout << __func__ << ":R r="<< p22/p12 <<T3endl;
#endif
          d1 = P - r2 + dr;
          cost = ( d1 * ( a[ct2] + P * ( b[ct2] + P * c[ct2] ) ) + d2 * cost ) / dr;
        }
        // else cost = cost;                        // Left edge solution
      }
      else cost = ( d1 * ( a[ct2] + P * ( b[ct2] + P * c[ct2] ) ) + d2 * cost ) / dr;
    }
  }
  else if(!hct && b[ct1] < 0.)                      // ==> Negative derivitive in 0.
  {
    const long double z = - b[ct1] / dd;
    const long double d1 = P;
    const long double d2 = r2 - P;
    cost = (d2 * z * (P * dd - 1.) + d1 * cost ) / ( d1 + d2 * z );
    //G4cout<<"T2NECS::CC: d="<<dd<<", x="<<[ct1].x()<<", y="<<[ct1].y()<<", z="<<[ct1].z()<<G4endl;
  }
  else if(hct == _num_point && b[ct1] + 2 * c[ct1] < 0.)    // ==> Negative derivitive in 1.
  {
    const long double z = - ( b[ct1] + 2 * c[ct1] ) / dd;
    const long double d1 = 1. - P;
    const long double d2 = P - r1;
    cost = ( d2 * z * ( 1. - ( 1. - P ) * dd ) + d1 * cost ) / ( d1 + d2 * z );
    //G4cout<<"T2NECS::CCT: Cor1 = "<< [ct1].y() + 2 * [ct1].z() <<G4endl;
  }
  return T3double(cost);
}

T3NSGangular_node& T3NSGangular_node::operator=(T3NSGangular_node const &node)
{
  memcpy(this->V, node.V, _num_point * sizeof(T3double));
  memcpy(this->a, node.a, _num_point * sizeof(T3double));
  memcpy(this->b, node.b, _num_point * sizeof(T3double));
  memcpy(this->c, node.c, _num_point * sizeof(T3double));
  sl = node.Get_sl();
  pr = node.Get_pr();
  return (*this);
}


std::ostream& operator<<(std::ostream& os, const T3NSGangular_node& inst)
{
  os << "Node" << " pr = " << inst.Get_pr() << ", sl = "
     << inst.Get_sl() << ", integrated distribution:";
  for(size_t ind =1; ind <= inst._num_point; ++ind) os << " " << inst.Get_V(ind-1);
  return os;
}
