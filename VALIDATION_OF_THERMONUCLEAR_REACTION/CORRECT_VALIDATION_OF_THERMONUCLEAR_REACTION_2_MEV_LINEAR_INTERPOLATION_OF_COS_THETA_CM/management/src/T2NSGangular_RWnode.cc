#include "T2NSGangular_RWnode.hh"

// #include "Randomize.hh"
// #include "T2Utility.hh"

T2NSGangular_RWnode::T2NSGangular_RWnode(G4double e,
                                         const std::array<G4double, _num_point> src,
                                         G4double newpr, G4double newsl) : E(e),
                                         sl(newsl), pr(newpr), V(src)
{}

T2NSGangular_RWnode::T2NSGangular_RWnode(T2NSGangular_RWnode const &node): E(node.E),
                                         sl(node.sl), pr(node.pr), V(node.V){}

T2NSGangular_RWnode& T2NSGangular_RWnode::operator=(T2NSGangular_RWnode const &node)
{
  E = node.E;
  sl = node.sl;
  pr = node.pr;
  V = node.V;
  return *this;
}

G4double T2NSGangular_RWnode::Get_V(size_t point_num) const
{
  if (point_num > _num_point + 1)
  {
    G4cout << "T2NSGangular_RWnode:Get_V(" << point_num << ") point number out of range"
    << G4endl;
    return 0; // TODO return NaN
  }
  else if (0 == point_num)              return 0;
  else if (_num_point + 1 == point_num) return 1;
  else                                  return V.at(point_num - 1);
}

void T2NSGangular_RWnode::Set_V(size_t point_num, G4double val)
{
  if ( _num_point < point_num || 0 == point_num)
  {
    G4cout << "T2NSGangular_RWnode:Set_V(" << point_num << ") point number out of range"
           << G4endl;
    return;
  }
  else V.at(point_num - 1) = val;
}

void T2NSGangular_RWnode::Set_V(const std::array<G4double, _num_point> src)
{
  // TODO add 0 < x < 1; x_i < x_(i+1) checks
  V = src;
}

std::ofstream& T2NSGangular_RWnode::save_binary(std::ofstream& out_stream) const
{
  if( out_stream.good() )
  {
    unsigned int vector_size = _num_point; // (I(0)=0,I(16)=1)
    out_stream.write((const char*) &E, sizeof(G4double) );
    out_stream.write((const char*) V.data(), vector_size * sizeof(G4double));
    out_stream.write((const char*) &pr, sizeof(G4double) );
    out_stream.write((const char*) &sl, sizeof(G4double) );
  }
  else
  {
    G4cout<<"-Warning-T2NSGangular_RWnode::save_binary:*Bad stream*"<<G4endl;
  }
  return out_stream;
}

std::ifstream& T2NSGangular_RWnode::load_binary(std::ifstream& in_stream )
{
  if( in_stream.good() )
  {
    unsigned int vector_size = _num_point; // (I(0)=0,I(16)=1)
    in_stream.read((char*) & E, sizeof(G4double) );
    in_stream.read((char*) V.data(), vector_size * sizeof( G4double) );
    in_stream.read((char*) & pr, sizeof(G4double) );
    in_stream.read((char*) & sl, sizeof(G4double) );
  }
  else G4cout<<"-Warning-T2NSGangular_RWnode::load_binary: *Bad stream*"<<G4endl;
  return in_stream;
}


std::array<G4double, T2NSGangular_RWnode::_num_point> T2NSGangular_RWnode::isotropic()
{
  std::array<G4double, _num_point> array;
  static const G4double nbins = _num_point + 1;
  for (size_t ind = 0; ind < _num_point; ++ ind) array.at(ind) = G4double(ind+1)/nbins;
  return array;
}


std::ostream& operator<<(std::ostream& os, const T2NSGangular_RWnode& inst)
{
  os << "Node e = " << inst.Get_E() << ", pr = " << inst.Get_pr() << ", sl = "
     << inst.Get_sl() << ", integrated distribution:";
  for(size_t ind =1; ind <= inst._num_point; ++ind) os << " " << inst.Get_V(ind);
  return os;
}

