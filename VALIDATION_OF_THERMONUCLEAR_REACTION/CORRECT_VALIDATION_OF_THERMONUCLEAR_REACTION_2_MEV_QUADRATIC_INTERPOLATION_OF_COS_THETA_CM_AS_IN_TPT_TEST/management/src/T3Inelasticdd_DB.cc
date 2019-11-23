#include "T3Inelasticdd_DB.hh"
#include "T3Utility.hh"

// copypasted from T2NElasticCS
T3double T3Inelasticdd_DB::RandomizeCost_node(RNDGenerator & generator, size_t node_num) const
{
  T3double P = t3::GenerateSubCanonical<FloatingType>(generator);  
  if(P < pr[node_num])
  {
    // use exponent
    return 1. + sl[node_num] * log(1. - ( 1. - exp(-2 / sl[node_num]) ) *
                         t3::GenerateSubCanonical<FloatingType>(generator) );
  }
  //else use tabulated
  static const T3double dY = 2. / (_num_point + 1); // Y step
  //squared the square roots of partial sums when initialized T3Inelasticdd_DB db),
  //that is why do not need sqrt() here.
  P = /*sqrt(*/t3::GenerateSubCanonical<FloatingType>(generator)/*)*/;
  const size_t hct = T3Utility::bin_search<127>(V[node_num], P) - V[node_num];
  T3int const ct1 = hct ? hct-1 : hct;                   // 0-bin default for ct1 = hct = 0
  T3int const ct2 = (hct && hct < _num_point) ? hct : -1;// Doesn't have any meaning for 0
  T3double const r1 = hct ? V[node_num][ct1] : 0.;    // 0-bin default for the Left End
  T3double const r2 = (hct < _num_point) ? V[node_num][hct] : 1.;   // the Right End
  long double cost = a[node_num][ct1] + P*( b[node_num][ct1] + P * c[node_num][ct1] ); // X+Y*R+Z*R^2 from r1
  const long double dd = dY / ( r2 - r1);                 // the 1D slope in the bin
  if(ct2 >= 0)
  {
    const T3double dr1 = r1 + r1;
    const T3double dr2 = r2 + r2;
    const long double p12 = b[node_num][ct1] + dr2 * c[node_num][ct1];        // Derivitive of ** f1 ** in r2
    const long double p21 = b[node_num][ct2] + dr1 * c[node_num][ct2];        // Derivitive of ** f2 ** in r1
    if(p12 < 0. || p21 < 0.)                        // spline in the reduced interval
    {
      long double F  =  a[node_num][ct2] + P*( b[node_num][ct2] + P * c[node_num][ct2] ); // Original f2 curve (no cor)
      const long double fs = dY * hct - 1.;               // LinFunct = fs + x * dd
      const long double d1 = P - r1;
      const long double d2 = r2 - P;
      if( p12 < 0. )                                // 1st deriv < 0 in r2 => negPart->lin
      {
        T3double z = - ( b[node_num][ct1] + r2 * 2 * c[node_num][ct1] ) / dd;
        cost = ( d1 * z * ( fs + ( P - r1 ) * dd ) + d2 * cost ) / ( d1 * z + d2 );
      }
      if( p21 < 0. )                                // 2nd deriv < 0 in r1 => negPart->lin
      {
        const long double z = - ( b[node_num][ct2] + r1 * 2 * c[node_num][ct2] ) / dd;
        F = (d2 * z * ( fs + ( P - r1 ) * dd ) + d1 * F ) / ( d1 + d2 * z );
      }
      cost = ( d1 * F + d2 * cost ) / ( d1 + d2);
    }
    else                                                  //for the right end.
    {
      const long double p11 = b[node_num][ct1] + dr1 * c[node_num][ct1];      // Derivitive of ** f1 ** in r1
      const long double p22 = b[node_num][ct2] + dr2 * c[node_num][ct2];      // Derivitive of ** f2 ** in r2
      long double d1 = P - r1;
      long double d2 = r2 - P;
      long double dr = r2 - r1;
      if(p11 > 0 && p21 > 0 && p11 > p21)     // Left edge derivitive problem
      {
        dr /= p11 / p21;
        if(d1 < dr)
        {
          d2 = r1 + dr - P;
          cost = ( d1 * ( a[node_num][ct2] + P * ( b[node_num][ct2] + P * c[node_num][ct2] ) ) + d2 * cost ) / dr;
        }
        else cost = a[node_num][ct2] + P * ( b[node_num][ct2] + P * c[node_num][ct2] ); // Rirht edge approximation
      }
      else if (p12 > 0 && p21 > 0 && p22 > p12)     // Right edge derivitive problem
      {
        dr /= p22 / p12;
        if(d2 < dr)
        {
          d1 = P - r2 + dr;
          cost = ( d1 * ( a[node_num][ct2] + P * ( b[node_num][ct2] + P * c[node_num][ct2] ) ) + d2 * cost ) / dr;
        }
        // else cost = cost;                        // Left edge solution
      }
      else cost = ( d1 * ( a[node_num][ct2] + P * ( b[node_num][ct2] + P * c[node_num][ct2] ) ) + d2 * cost ) / dr;
    }
  }
  else if(!hct && b[node_num][ct1] < 0.)                      // ==> Negative derivitive in 0.
  {
    const long double z = - b[node_num][ct1] / dd;
    const long double d1 = P;
    const long double d2 = r2 - P;
    cost = (d2 * z * (P * dd - 1.) + d1 * cost ) / ( d1 + d2 * z );
  }
  else if(hct == _num_point && b[node_num][ct1] + 2 * c[node_num][ct1] < 0.)    // ==> Negative derivitive in 1.
  {
    const long double z = - ( b[node_num][ct1] + 2 * c[node_num][ct1] ) / dd;
    const long double d1 = 1. - P;
    const long double d2 = P - r1;
    cost = ( d2 * z * ( 1. - ( 1. - P ) * dd ) + d1 * cost ) / ( d1 + d2 * z );
  }
  return T3double(cost);
}


T3double T3Inelasticdd_DB::RandomizeCost(RNDGenerator & generator, T3double E) const
{
  const T3double Emin = Einc[0];
  const T3double Emax = Einc[_num_nodes - 1];
  const T3double Ediff = Emax - Emin;
  const size_t nbins = _num_nodes - 1;
  const T3double dE = Ediff / nbins;
  if(E <= Emin || E >= Emax)
  {
    // isotropic
    const T3double R = t3::GenerateSubCanonical<FloatingType>(generator);
    return -1 + 2*R;
  }
  else
  {
    size_t const binn = static_cast<size_t>((E - Emin)*(_num_nodes-1)/Ediff);
    const T3double Elow  = Einc[binn];
    const T3double W = (E - Elow)/(dE);
    const T3double R = t3::GenerateSubCanonical<FloatingType>(generator);
    if(R < W)//randomize by high distribution
      return RandomizeCost_node(generator, binn+1);
    else//randomize by low distribution
      return RandomizeCost_node(generator, binn);
  }
}






