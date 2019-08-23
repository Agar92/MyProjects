#include "T3NSGangular_RWrecord.h"

namespace t3 {
std::ofstream& T3NSGangular_RWrecord::save_binary(std::ofstream& out_stream) const
{
  if( out_stream.good() )
  {
    const unsigned int vector_size = size();
    out_stream.write((const char*) &vector_size, sizeof(unsigned int) );
    for(T3NSGangular_RWrecord::const_iterator it = begin(); it != end(); ++it)
    {
      it->save_binary(out_stream);
    }
  }
  else T3cout<<"-Warning-T3NSGangular_RWrecord::save_binary:*Bad stream*"<<T3endl;
  return out_stream;
}

std::ifstream& T3NSGangular_RWrecord::load_binary(std::ifstream& in_stream )
{
  if( in_stream.good() )
  {
    unsigned int vector_size = 0;
    in_stream.read((char*) &vector_size, sizeof(unsigned int) );

    std::cout<<"T3NSGangular_RWrecord::load_binary(std::ifstream& in_stream) 1: vector_size="
             <<vector_size<<std::endl;

    resize(vector_size);

    std::cout<<"T3NSGangular_RWrecord::load_binary(std::ifstream& in_stream) 2: vector_size="
             <<vector_size<<std::endl;

    for(T3NSGangular_RWrecord::iterator it = begin(); it != end(); ++it)
    {
      it->load_binary(in_stream);
    }
  }
  else T3cout<<"-Warning-T2NSGangular_RWrecord::load_binary: *Bad stream*"<<T3endl;
  return in_stream;
}

T3NSGangular_RWnode T3NSGangular_RWrecord::interpolate(T3double E) const
{
  T3NSGangular_RWnode res;
  if (!size())
  {
    T3cout << "T3NSGangular_RWrecord::interpolate: empty record" << T3endl;
    return res;
  }
  if( front().Get_E() > E)
  {
    T3cout << "T3NSGangular_RWrecord::interpolate: too high E=" << E << " > "
           << front().Get_E() << T3endl;
    return res;
  }
  // profile, maybe make logarithmic (std::lower_bound with less by E)
  // currently linear
  size_t ind =0;
  for(; ind + 1 < size(); ++ind) if( at(ind+1).Get_E() > E) break;
  if (ind + 1 >= size() ) return res;
  const T3double Elow  = at(ind    ).Get_E();
  const T3double Ehigh = at(ind + 1).Get_E();
  const T3double dE = Ehigh - Elow;
#ifdef debug
  T3cout << "E=" << E << T3endl;
#endif
  for(size_t i = 1; i <= T3NSGangular_RWnode::_num_point; ++i)
  {
    const T3double vlow  = at(ind    ).Get_V(i);
    const T3double vhigh = at(ind + 1).Get_V(i);
    if (vlow < 0.0 || 1.0 < vlow)
    {
      T3cout << "T3NSGangular_RWrecord::interpolate: vlow=" << vlow << " i=" << i
             << " ind=" << ind << " Elow=" << Elow << " E=" << E << " Ehigh=" << Ehigh <<T3endl;
      T3cout << *this << T3endl;
    }
    if (vhigh < 0.0 || 1.0 < vhigh)
    {
      T3cout << "T3NSGangular_RWrecord::interpolate: vhigh=" << vhigh << " i=" << i
             << " ind=" << ind << " Elow=" << Elow << " E=" << E << " Ehigh=" << Ehigh <<T3endl;
      T3cout << *this << T3endl;
    }
    
    if (dE < 1e-10) res.Set_V(i, (vlow + vhigh)/2);
    else
    {
      const T3double val = vlow + (vhigh - vlow)*(E - Elow)/dE;
      if (val < 0.0 || 1.0 < val)
      {
        T3cout << "T3NSGangular_RWrecord::interpolate: val=" << val << " vlow=" << vlow
               << " vhigh=" << vhigh << " Elow=" << Elow << " Ehigh=" << Ehigh 
               << " E=" << E << T3endl;
      }
      res.Set_V(i, val );
    }
  }
  return res;
}

std::ostream& operator<<(std::ostream& os, const T3NSGangular_RWrecord& inst)
{
  os << "Record nodes:";
  for(size_t ind =0; ind < inst.size(); ++ind) os << "\n\t" << inst.at(ind);
  return os;
}
} // namespace t3





