#include "T2NSGangular_RW.hh"
//#include <T2PDGCode.hh>
#include <map>

//#define debug

std::ofstream& T2NSGangular_RW::save_binary(std::ofstream& out_stream) const
{
  if( out_stream.good() )
  {
    const unsigned int vector_size = size();
    out_stream.write((const char*) &Z, sizeof(G4int) );
    out_stream.write((const char*) &A, sizeof(G4int) );
    out_stream.write((const char*) &secPDG, sizeof(G4int) );
    out_stream.write((const char*) &rID, sizeof(G4int) );
    out_stream.write((const char*) &vector_size, sizeof(unsigned int) );
    for(T2NSGangular_RW::const_iterator it = begin(); it != end(); ++it)
    {
      it->save_binary(out_stream);
    }
  }
  else
  {
    G4cout<<"-Warning-T2NSGangular_RW::save_binary:*Bad stream*"<<G4endl;
  }
  return out_stream;
}

std::ifstream& T2NSGangular_RW::load_binary(std::ifstream& in_stream )
{
  if( in_stream.good() )
  {
    unsigned int vector_size = 0;
    in_stream.read((char*) &Z, sizeof(G4int) );
    in_stream.read((char*) &A, sizeof(G4int) );
    in_stream.read((char*) &secPDG, sizeof(G4int) );
    in_stream.read((char*) &rID, sizeof(G4int) );
    in_stream.read((char*) &vector_size, sizeof(unsigned int) );
    resize(vector_size);
    for(T2NSGangular_RW::iterator it = begin(); it != end(); ++it)
    {
      it->load_binary(in_stream);
    }
  }
  else
  {
    G4cout<<"-Warning-T2NSGangular_RW::load_binary: *Bad stream*"<<G4endl;
  }
#ifdef debug
  G4cout << *this << G4endl;
#endif
  return in_stream;
}

void T2NSGangular_RW::save_binary( const G4String& fname ) const
{
  std::ofstream out_stream;
  out_stream.open(fname.c_str(), std::ios::out | std::ios::binary | std::ios::trunc);
  if(out_stream.good())
  {
    save_binary(out_stream);
  }
  else
  {
#ifdef debug
    G4cout << "-Warning-T2NSGangular_RW::failed to open " << fname << G4endl;
#endif
  }
  out_stream.close();
}

G4bool T2NSGangular_RW::load_binary( const G4String& fname )
{
  G4bool result = false;
  std::ifstream in_stream;
  in_stream.open(fname.c_str(), std::ifstream::in | std::ifstream::binary);
  if(in_stream.good())
  {
    load_binary(in_stream);
    result = size();
  }
  else
  {
#ifdef debug
    G4cout << "-Warning-T2NSGangular_RW::failed to open " << fname << G4endl;
#endif
  }
  in_stream.close();
  return result;
}

G4String T2NSGangular_RW::default_file(G4int tgZ, G4int tgA,
                                       const G4String& rid, G4int pdg, G4int incZA) const
{
  std::map<G4int, G4String> pnames{{1, "N"},
                                   {1001, "P"},
                                   {1002, "D"},
                                   {1003, "T"},
                                   {2003, "H"},
                                   {2004, "A"}};
//  pnames[2112] = "n";
  char buffer[1023];
  G4String T2data;                // Check that the NNG data exist
  if (getenv("T2_DATA"))
  {
    T2data =  getenv("T2_DATA");
  }
  else
  {
    T2data = G4String("./data");
    G4cout << "-Warning-T2NSGangular_RW::default_file:Cant't get $T2_DATA. Use ./data"
           << G4endl;
  }
  const G4String pname = pnames[pdg];
  G4String incFolder = (incZA == 1) ? "n" : ("inc" + pnames[incZA]);
  sprintf(buffer, "%s/angular/%s/%s/%s/T2%sSGangular_10%.3d%.3d0.bin", T2data.c_str(),
          pname.c_str(), incFolder.c_str(), rid.c_str(), pnames[incZA].c_str(), tgZ, tgA);
#ifdef debug
  G4cout << "T2NSGangular_RW::default_file: " << buffer << G4endl;
#endif

//Check:
  std::cout<<"Check: "<<buffer<<std::endl;
//
  
  return G4String(buffer);
}

G4bool T2NSGangular_RW::load_binary(G4int tgZ, G4int tgA, const G4String& rid, G4int pdg,
                                    G4int incZA)
{
  return load_binary(default_file( tgZ, tgA, rid, pdg, incZA ));
}

void T2NSGangular_RW::save_binary( G4int tgZ, G4int tgA, const G4String& rid, G4int pdg,
                                   G4int incZA) const
{
  save_binary(default_file( tgZ, tgA, rid, pdg, incZA ));
}

std::ostream& operator<<(std::ostream& os, const T2NSGangular_RW& inst)
{
  os << "Z=" << inst.Get_Z() << " A=" << inst.Get_A() << " secPDG = " << inst.Get_secPDG()
     << "rID = " << inst.Get_rID() << ", records:";
  for(size_t ind =0; ind < inst.size(); ++ind) os << "\n[" << inst.at(ind) << "]";
  return os;
}
