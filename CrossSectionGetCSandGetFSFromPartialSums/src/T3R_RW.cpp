// #define debug
// #define pdebug

#include "T3R_RW.h"
#include <map>

namespace t3 {

T3R_RW::T3R_RW(const T3TabulatedCS& tcs, T3int _tgZ, T3int _tgA):
  tgZ(_tgZ), tgA(_tgA)
{
  //TODO getenv once at static initialization ----------------------------------------
  const char* const t3data = getenv("T3_DATA");
  int channel=-1;
  if (t3data)
  {
#ifdef debug
    T3cout<<"T3_DATA="<< t3data <<T3endl;
#endif
    T3data =  T3String(t3data);
    channel=1;
  }
  else
  {
    T3data = T3String("./data");
    T3cout<<"-Warning-T3R_RW::default_file: Cant't get $T3_DATA. Using ./data" <<T3endl;
    channel=2;
  }
#ifdef debug
  T3cout<<"T3_DATA =  "<< T3data << ", str="<< T3data.c_str() << T3endl;
#endif
//
  T3cout<<"t3data="<<t3data<<" T3data="<<T3data<<" channel="<<channel<<std::endl;
//
  //------------------------------------------------------------------------------------
  for(size_t ind = 0; ind < tcs.Get_size(); ++ind)
  {
    push_back(new T3R_node(tcs.Get_E(ind), tcs.Get_xs(ind)));
  }
  //std::cout<<"~~~~SIZE="<<tcs.Get_size()<<std::endl;
}

T3R_RW::~T3R_RW()
{
  T3int iD = ND.size();
  std::cout<<"1. ~T3R_RW(): iD="<<iD<<(*this)<<std::endl;
  if(iD) for(T3int jd=0; jd<iD; ++jd) if(ND[jd]) delete ND[jd];
  std::cout<<"2. ~T3R_RW(): iD="<<iD<<std::endl;
}
//WRITE:
void T3R_RW::save_binary(const T3String& fname) const // Write R DB to file
{
  std::ofstream out_stream;
//std::ios::binary - without formatting
//std::ios::trunc - if anything was in the file it is deleted
  out_stream.open(fname.c_str(), std::ios::out | std::ios::binary | std::ios::trunc);

  std::cout<<"save_binary(): 1"<<std::endl;
  
  if(!out_stream.good())
  {
#ifdef debug
    T3cout<<"-Warning-T3R_RW::save_binary:*Bad Stream*"<<T3endl;
#endif
    return;
  }

  std::cout<<"save_binary(): 2"<<std::endl;
  
  unsigned int vector_size = ND.size();

  std::cout<<"|||Check save_binary() vector_size="<<vector_size<<std::endl;

  std::cout<<"save_binary(): 3"<<std::endl;
  
#ifdef debug
  T3cout<<"T3R_RW::save_binary: ND="<<vector_size<<T3endl;
#endif

  std::cout<<"save_binary(): 4"<<std::endl;
  
  out_stream.write(  (const char*) &( vector_size ), sizeof( unsigned int ) );

  std::cout<<"save_binary(): 5"<<std::endl;
  
  for( unsigned int i = 0; i < vector_size; ++i)
  {
#ifdef debug
    T3cout<<"T3R_RW::save_binary: node # "<< i <<", E="<< ND[i]->E() <<", XS="
          << ND[i]->XS() << T3endl;
#endif
    ND[i]->save_binary( &out_stream );
  }

  std::cout<<"save_binary(): 6"<<std::endl;


////////////////////////////////////////////////////////
//Here there was a redirection of std::cout:          //
//psbuf - stream buffer of an out_stream              //
//std::cout.rdbuf(psbuf); - redirection of the output //
//stream out_stream to std::cout.                     //
//This piece of code does not work. It throws         //
//save_binary(): 7                                    //
//Segmentation fault (core dumped)                    //
//There is no segmentation fault only if to uncomment //
//1) and 2) or 3)                                     //
//But this piece of code does the work,               //
//out_stream.write(...) does above, for the second    //
//time. It has no sense.                              //  
////////////////////////////////////////////////////////  
//  std::streambuf * psbuf=out_stream.rdbuf();        //
//  //1)://std::streambuf * coutbuf=std::cout.rdbuf();//
//  std::cout<<"save_binary(): 7"<<std::endl;         //
//  std::cout.rdbuf(psbuf);                           //
//  //2)://std::cout.rdbuf(coutbuf);                  //
//  //3)://out_stream.close();                        //
//                                                    //
////////////////////////////////////////////////////////  
  
  std::cout<<"save_binary(): 8"<<std::endl;
  
}
//READ:
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

  std::cout<<"|||Check load_binary(): vector_size="<<vector_size<<std::endl;
  
  in_stream.close();
  return true;
}

T3String T3R_RW::default_file( T3int targZ, T3int targA,
                                 T3int pZ, T3int pA,
                                 const T3String& suffix/*name of database*/ ) const
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
  return default_file(targZ, targA, particle_name/* + "any"*/, suffix);
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
                               const T3String& sID/*name of the particle(P,D,T,H,A)+any*/,
                               const T3String& suffix/*name of database*/ ) const
{
  const T3int nbases = 4;
  T3String bases[nbases] = {"endf", "tendl", "custom", "unknown"};
  T3bool CSExist = false;
  char buffer[1023];
  T3String filename;
  T3String base;// = bases[0];  
  if (suffix == "any")//FILE OF ANY DATABASE
  {
    for( T3int index = 0; index < nbases; ++index)
    {
      base = bases[index];
      //sprintf puts the string in buffer:
      //T3data/xs/sID/T2sID_10targZtargA0_base.bin
//WAS:      
//    sprintf(buffer, "%s/xs/%s/T3%s_10%.3d%.3d0_%s.bin", T3data.c_str(),
//      sID.c_str(), sID.c_str(), targZ, targA, base.c_str());
      //There were 9 digits, but in big PDG book on page 315 there is
      //the standard code of 10 digits: 10LZZZAAAI.
      sprintf(buffer, "%s/xs/%s/T3%s_100%.3d%.3d0_%s.bin", T3data.c_str(),
              sID.c_str(), sID.c_str(), targZ, targA, base.c_str());
      filename = T3String(buffer);
      std::ifstream try_file(filename);
      if ( !try_file.good() )
      {
#ifdef debug
        T3cout<<"*WARNING*T3R_RW::CCS:absentCSFile="<<filename<<T3endl;
#endif
        CSExist = false;
        try_file.close();
      }
      else
      {
        CSExist = true;
        break;         //returns the datafile of any database
      }
    }
    if (CSExist) return filename;//returns the filename
    else//if nothing found:
    {
      // @@ TODO: find the file with any suffix
//    sprintf(buffer, "%s/xs/%s/T3%s_10%.3d%.3d0_%s.bin", T3data.c_str(),
//      sID.c_str(), sID.c_str(), targZ, targA, "<db>");
      //There were 9 digits, but in big PDG book on page 315 there is
      //the standard code of 10 digits: 10LZZZAAAI.
      sprintf(buffer, "%s/xs/%s/T3%s_100%.3d%.3d0_%s.bin", T3data.c_str(),
              sID.c_str(), sID.c_str(), targZ, targA, "<db>");
      filename = T3String(buffer);
#ifdef debug
      T3cout<<"*WARNING*T3R_RW::CCS:(SYNC) no CSFile, return "<<filename<<T3endl;
#endif
      return filename; // TODO return special name or emit warning
    }
  }
  else//FILE OF DATABASE WITH THE NAME "SUFFIX"
  {
    base = suffix;
//  sprintf(buffer, "%s/xs/%s/T3%s_10%.3d%.3d0_%s.bin", T3data.c_str(),
//      sID.c_str(), sID.c_str(), targZ, targA, base.c_str());
    //There were 9 digits, but in big PDG book on page 315 there is
    //the standard code of 10 digits: 10LZZZAAAI.
    sprintf(buffer, "%s/xs/%s/T3%s_100%.3d%.3d0_%s.bin", T3data.c_str(),
            sID.c_str(), sID.c_str(), targZ, targA, base.c_str());
    filename = T3String(buffer);
    std::ifstream try_file(filename);
#ifdef debug
    if(!try_file.good()) T3cout<<"T3R_RW::CCS:absentNecessaryCSFile="<<filename<<T3endl;
#endif
    try_file.close();
    std::cout<<"HHHH="<<filename<<" try="<<try_file.good()<<std::endl;
    return filename;
  }
}

std::ostream& operator<<(std::ostream& os, const T3R_RW& inst)
{
  for(size_t ind = 0; ind < inst.size(); ++ ind)
  {
    os << "node #" << ind << ' ' << *(inst.NDI(ind)) << T3endl;
  }
  return os;
}

} // namespace t3
