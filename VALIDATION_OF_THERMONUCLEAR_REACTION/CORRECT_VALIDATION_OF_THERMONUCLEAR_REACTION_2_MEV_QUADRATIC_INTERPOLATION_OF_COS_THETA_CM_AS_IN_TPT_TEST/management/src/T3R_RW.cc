// #define debug
// #define pdebug

#include "T3R_RW.hh"
#include <map>

T3R_RW::T3R_RW(const T3TabulatedCS& tcs, T3int _tgZ, T3int _tgA):
  tgZ(_tgZ), tgA(_tgA)
{
  // TODO getenv once at static initialization ----------------------------------------
  const char* const t2data = getenv("T2_DATA");
  if (t2data)
  {
#ifdef debug
    T3cout<<"T2_DATA="<< t2data <<T3endl;
#endif
    T2data =  T3String(t2data);
  }
  else
  {
    T2data = T3String("./data");
    T3cout<<"-Warning-T3R_RW::default_file: Cant't get $T2_DATA. Using ./data" <<T3endl;
  }
#ifdef debug
  T3cout<<"T2_DATA =  "<< T2data << ", str="<< T2data.c_str() << T3endl;
#endif
  //------------------------------------------------------------------------------------
  for(size_t ind = 0; ind < tcs.Get_size(); ++ind)
  {
    push_back(new T3R_node(tcs.Get_E(ind), tcs.Get_xs(ind)));
  }
}

T3R_RW::~T3R_RW()
{
  T3int iD = ND.size();
  if(iD) for(T3int jd=0; jd<iD; ++jd) if(ND[jd]) delete ND[jd];
}

void T3R_RW::save_binary(const T3String& fname) const // Write R DB to file
{
  std::ofstream out_stream;
  out_stream.open(fname.c_str(), std::ios::out | std::ios::binary | std::ios::trunc);
  if(!out_stream.good())
  {
#ifdef debug
    T3cout<<"-Warning-T3R_RW::save_binary:*Bad Stream*"<<T3endl;
#endif
    return;
  }
  unsigned int vector_size = ND.size();
#ifdef debug
  T3cout<<"T3R_RW::save_binary: ND="<<vector_size<<T3endl;
#endif
  out_stream.write(  (const char*) &( vector_size ), sizeof( unsigned int ) );
  for( unsigned int i = 0; i < vector_size; ++i)
  {
#ifdef debug
    T3cout<<"T3R_RW::save_binary: node # "<< i <<", E="<< ND[i]->E() <<", XS="
          << ND[i]->XS() << T3endl;
#endif
    ND[i]->save_binary( &out_stream );
  }
}

T3bool T3R_RW::load_binary(const T3String& fname) // Read NNG DB from file
{
  std::ifstream in_stream;
  in_stream.open( fname.c_str(), std::ios::in | std::ios::binary);
  if(!in_stream.good())
  {
#ifdef debug
    T3cout<<"-Warning-T3R_RW::load_binary:*Bad Stream*"<<T3endl;
#endif
    return false;
  }
  unsigned int vector_size;
  // Now read Nodes
  in_stream.read( (char*) &vector_size, sizeof( unsigned int ) );
#ifdef debug
  T3cout<<"T3R_RW::load_binary: ND="<< vector_size <<", fname=" << fname << T3endl;
#endif
  T3int iD = ND.size();
  if(iD) for(T3int jd=0; jd<iD; ++jd) delete ND[jd];
  ND.clear();
  T3R_node* newND = 0;
  for( unsigned int i = 0; i < vector_size; ++i)
  {
    newND = new T3R_node();
#ifdef debug
    T3cout<<"T3R_RW::load_binary: Loading node # "<< i <<T3endl;
#endif
    newND->load_binary( &in_stream ); // @@ temporary old format is used
#ifdef debug
    T3cout<<"T3R_RW::load_binary: node # "<< i <<", E="<< newND->E() <<", XS="
          << newND->XS() << T3endl;
    if(newND->XS() < 0.)
    {
      T3cout << "T3R_RW::load_binary: Warning: xs("<<newND->E()<<") =" << newND->XS()
             << " < 0" << T3endl;
    }
#endif
    ND.push_back(newND);
  }
  in_stream.close();
  return true;
}

T3String T3R_RW::default_file( T3int targZ, T3int targA,
                               T3int pZ, T3int pA, const T3String& suffix ) const
{
  std::map<T3int, T3String> names;
//   names[1000*0 + 0] = "G";
//   names[1000*0 + 1] = "N";
  names[1000*1 + 1] = "P";
  names[1000*1 + 2] = "D";
  names[1000*1 + 3] = "T";
  names[1000*2 + 3] = "H";
  names[1000*2 + 4] = "A";
  T3String particle_name = names[1000 * pZ + pA];
  return default_file(targZ, targA, particle_name + "any", suffix);
}

void T3R_RW::save_binary( T3int targZ, T3int targA,
                          T3int pZ, T3int pA, const T3String& suffix ) const
{
  save_binary( default_file( targZ, targA, pZ, pA, suffix ) );
}

void T3R_RW::save_binary( T3int targZ, T3int targA,
                          const T3String& sID, const T3String& suffix ) const
{
  save_binary( default_file( targZ, targA, sID, suffix ) );
}

T3bool T3R_RW::load_binary(T3int targZ, T3int targA,
                           T3int pZ, T3int pA, const T3String& suffix)
{
  return load_binary( default_file( targZ, targA, pZ, pA, suffix ) );
}

T3bool T3R_RW::load_binary(T3int targZ, T3int targA,
                             const T3String& rID, const T3String& suffix)
{
  return load_binary( default_file( targZ, targA, rID, suffix ) );
}

T3String T3R_RW::default_file( T3int targZ, T3int targA,
                                 const T3String& sID, const T3String& suffix ) const
{
  const T3int nbases = 4;
  T3String bases[nbases] = {"endf", "tendl", "custom", "unknown"};
  T3bool CSExist = false;
  char buffer[1023];
  T3String filename;
  T3String base;// = bases[0];
  if (suffix == "any")
  {
    for( T3int index = 0; index < nbases; ++index)
    {
      base = bases[index];
      sprintf(buffer, "%s/xs/%s/T2%s_10%.3d%.3d0_%s.bin", T2data.c_str(),
              sID.c_str(), sID.c_str(), targZ, targA, base.c_str());
      filename = T3String(buffer);
      std::ifstream try_file(filename);
      if ( !try_file.good() )
      {
#ifdef debug
        T3cout<<"*WARNING*T2R_RW::CCS:absentCSFile="<<filename<<T3endl;
#endif
        CSExist = false;
        try_file.close();
      }
      else
      {
        CSExist = true;
        break;
      }
    }
    if (CSExist) return filename;
    else
    {
      // @@ TODO: find the file with any suffix
      sprintf(buffer, "%s/xs/%s/T2%s_10%.3d%.3d0_%s.bin", T2data.c_str(),
              sID.c_str(), sID.c_str(), targZ, targA, "<db>");
      filename = T3String(buffer);
#ifdef debug
      T3cout<<"*WARNING*T3R_RW::CCS:(SYNC) no CSFile, return "<<filename<<T3endl;
#endif
      return filename; // TODO return special name or emit warning
    }
  }
  else
  {
    base = suffix;
    sprintf(buffer, "%s/xs/%s/T2%s_10%.3d%.3d0_%s.bin", T2data.c_str(),
            sID.c_str(), sID.c_str(), targZ, targA, base.c_str());
    filename = T3String(buffer);
    std::ifstream try_file(filename);
#ifdef debug
    if(!try_file.good())T3cout<<"T3R_RW::CCS:absentNecessaryCSFile="<<filename<<T3endl;
#endif
    try_file.close();
    return filename;
  }
}

// void T2R_RW::save_binary(G4int targZ, G4int targA, G4int pZ, G4int pA, G4int sZ, G4int sA,
//                          const G4String& suffix) const
// {
//   save_binary( default_file( targZ, targA, pZ, pA, sZ, sA, suffix ) );
// }

// G4bool T2R_RW::load_binary(G4int targZ, G4int targA, G4int pZ, G4int pA, G4int sZ,
//                            G4int sA, const G4String& suffix)
// {
//   return load_binary( default_file( targZ, targA, pZ, pA, sZ, sA, suffix ) );
// }

std::ostream& operator<<(std::ostream& os, const T3R_RW& inst)
{
  for(size_t ind = 0; ind < inst.size(); ++ ind)
  {
    os << "node #" << ind << ' ' << *(inst.NDI(ind)) << T3endl;
  }
  return os;
}
