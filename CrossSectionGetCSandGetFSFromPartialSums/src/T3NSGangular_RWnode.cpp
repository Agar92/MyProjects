#include "T3NSGangular_RWnode.h"

namespace t3 {
T3NSGangular_RWnode::T3NSGangular_RWnode(T3double e,
                                         const std::array<T3double, _num_point+1> src) :
                                         E(e), V(src){}

T3NSGangular_RWnode::T3NSGangular_RWnode(T3NSGangular_RWnode const &node): E(node.E),
                                         V(node.V){}

T3NSGangular_RWnode& T3NSGangular_RWnode::operator=(T3NSGangular_RWnode const &node)
{
  E = node.E;
  V = node.V;
  return *this;
}

T3double T3NSGangular_RWnode::Get_V(size_t point_num) const
{
  if (point_num > _num_point + 1)
  {
    return 0; // TODO return NaN
  }
  else if (0 == point_num)
  {
    return 0;
  }
  else if (_num_point + 2 == point_num)
  {
    return 1;
  }
  else
  {
    return V.at(point_num - 1);
  }
}

void T3NSGangular_RWnode::Set_V(size_t point_num, T3double val)
{
  if ( _num_point < point_num || 0 == point_num)
  {
    T3cout << "T3NSGangular_RWnode:Set_V(" << point_num << ") point number out of range"
           << T3endl;
    return;
  }
  else V.at(point_num - 1) = val;
}

void T3NSGangular_RWnode::Set_V(const std::array<T3double, _num_point+1> src)
{
  // TODO add 0 < x < 1; x_i < x_(i+1) checks
  V = src;
}

std::ofstream& T3NSGangular_RWnode::save_binary(std::ofstream& out_stream) const
{
  if( out_stream.good() )
  {
    unsigned int vector_size = _num_point+1; // (I(0)=0,I(16)=1)
    out_stream.write((const char*) &E, sizeof(T3double) );
    out_stream.write((const char*) V.data(), vector_size * sizeof(T3double));
  }
  else
  {
    T3cout<<"-Warning-T3NSGangular_RWnode::save_binary:*Bad stream*"<<T3endl;
  }
  return out_stream;
}

std::ifstream& T3NSGangular_RWnode::load_binary(std::ifstream& in_stream )
{
  if( in_stream.good() )
  {
    unsigned int vector_size = _num_point+1; // (I(0)=0,I(16)=1)
    in_stream.read((char*) & E, sizeof(T3double) );
    in_stream.read((char*) V.data(), vector_size * sizeof( T3double) );
  }
  else T3cout<<"-Warning-T3NSGangular_RWnode::load_binary: *Bad stream*"<<T3endl;
  return in_stream;
}

std::array<T3double, T3NSGangular_RWnode::_num_point+1> T3NSGangular_RWnode::isotropic()
{
  std::array<T3double, _num_point+1> array;
  static const T3double nbins = _num_point + 1;
  for (size_t ind = 0; ind < _num_point+1; ++ ind) array.at(ind) = T3double(ind+1)/nbins;
  return array;
}

std::ostream& operator<<(std::ostream& os, const T3NSGangular_RWnode& inst)
{
  os << "Node e = " << inst.Get_E() << ", integrated distribution:";
  for(size_t ind =1; ind <= inst._num_point+1; ++ind)
  {
    os << " " << inst.Get_V(ind);
  }
  return os;
}
} // namespace t3
