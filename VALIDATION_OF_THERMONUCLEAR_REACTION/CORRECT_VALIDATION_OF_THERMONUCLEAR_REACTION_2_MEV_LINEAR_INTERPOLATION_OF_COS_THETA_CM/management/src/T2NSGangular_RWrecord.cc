#include "T2NSGangular_RWrecord.hh"

std::ofstream& T2NSGangular_RWrecord::save_binary(std::ofstream& out_stream) const
{
  if( out_stream.good() )
  {
    const unsigned int vector_size = size();
    out_stream.write((const char*) &levN, sizeof(G4int) );
    out_stream.write((const char*) &vector_size, sizeof(unsigned int) );
    for(T2NSGangular_RWrecord::const_iterator it = begin(); it != end(); ++it)
    {
      it->save_binary(out_stream);
    }
  }
  else G4cout<<"-Warning-T2NSGangular_RWrecord::save_binary:*Bad stream*"<<G4endl;
  return out_stream;
}

std::ifstream& T2NSGangular_RWrecord::load_binary(std::ifstream& in_stream )
{
  if( in_stream.good() )
  {
    unsigned int vector_size = 0;
    in_stream.read((char*) &levN, sizeof(G4int) );
    in_stream.read((char*) &vector_size, sizeof(unsigned int) );
    resize(vector_size);
    for(T2NSGangular_RWrecord::iterator it = begin(); it != end(); ++it)
    {
      it->load_binary(in_stream);
    }
  }
  else G4cout<<"-Warning-T2NSGangular_RWrecord::load_binary: *Bad stream*"<<G4endl;
  return in_stream;
}

T2NSGangular_RWnode T2NSGangular_RWrecord::interpolate(G4double E) const
{
  T2NSGangular_RWnode res;
  if (!size())
  {
    G4cout << "T2NSGangular_RWrecord::interpolate: empty record" << G4endl;
    return res;
  }
  if( front().Get_E() > E)
  {
    G4cout << "T2NSGangular_RWrecord::interpolate: too high E=" << E << " > "
           << front().Get_E() << G4endl;
    return res;
  }
  // profile, maybe make logarithmic (std::lower_bound with less by E)
  // currently linear
  size_t ind =0;
  for(; ind + 1 < size(); ++ind) if( at(ind+1).Get_E() > E) break;
  if (ind + 1 >= size() ) return res;
  const G4double Elow  = at(ind    ).Get_E();
  const G4double Ehigh = at(ind + 1).Get_E();
  const G4double dE = Ehigh - Elow;
  const G4double sllow = at(ind    ).Get_sl();
  const G4double slhigh= at(ind + 1).Get_sl();
  const G4double prlow = at(ind    ).Get_pr();
  const G4double prhigh= at(ind + 1).Get_pr();
  // FIXME sum of exponents is not an exponent
  G4double sl;
  G4double pr;
  if (dE < 1e-10)
  {
    pr = (prlow + prhigh)/2;
    sl = (sllow + slhigh)/2;
  }
  else
  {
    pr = prlow + (prhigh - prlow)*(E - Elow)/dE;
    sl = sllow + (slhigh - sllow)*(E - Elow)/dE;
  }
  res.Set_E(E);
  res.Set_pr(pr);
  res.Set_sl(sl);
#ifdef debug
  G4cout << "E=" << E << " pr=" << pr << " sl=" << sl << G4endl;
#endif
  for(size_t i = 1; i <= T2NSGangular_RWnode::_num_point; ++i)
  {
    const G4double vlow  = at(ind    ).Get_V(i);
    const G4double vhigh = at(ind + 1).Get_V(i);
    if (vlow < -1 || 1 < vlow)
    { // TODO if pr > 1 - 1e-bigval fill non-exponent with isotropic
      // mybe if almost isotropic use exponent only
      G4cout << "T2NSGangular_RWrecord::interpolate: vlow=" << vlow << " i=" << i
             << " ind=" << ind << " Elow=" << Elow << " E=" << E << " Ehigh=" << Ehigh
             << " pr=" << pr << " sl=" << sl <<G4endl;
      G4cout << *this << G4endl;
    }
    if (vhigh < -1 || 1 < vhigh)
    {
      G4cout << "T2NSGangular_RWrecord::interpolate: vhigh=" << vhigh << " i=" << i
             << " ind=" << ind << " Elow=" << Elow << " E=" << E << " Ehigh=" << Ehigh
             << " pr=" << pr << " sl=" << sl <<G4endl;
      G4cout << *this << G4endl;
    }
    
    if (dE < 1e-10) res.Set_V(i, (vlow + vhigh)/2);
    else
    {
      const G4double val = vlow + (vhigh - vlow)*(E - Elow)/dE;
      if (val < -1 || 1 < val)
      {
        G4cout << "T2NSGangular_RWrecord::interpolate: val=" << val << " vlow=" << vlow
               << " vhigh=" << vhigh << " Elow=" << Elow << " Ehigh=" << Ehigh 
               << " E=" << E << G4endl;
      }
      res.Set_V(i, val );
    }
  }
  return res;
}

std::ostream& operator<<(std::ostream& os, const T2NSGangular_RWrecord& inst)
{
  os << "Record levN = " << inst.Get_levN() << ", nodes:";
  for(size_t ind =0; ind < inst.size(); ++ind) os << "\n\t" << inst.at(ind);
  return os;
}







