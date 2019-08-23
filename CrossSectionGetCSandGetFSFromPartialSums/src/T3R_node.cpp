//#define debug
//#define pdebug

#include "T3R_node.h"
#include <cstdio>
#include <fstream>
#include <iostream>

namespace t3 {

T3R_node::T3R_node(double Ein, double CS) // Ein & XS = 0
{
  e  = Ein;
  xs = CS;
}

void T3R_node::save_binary(std::ofstream* out_stream ) const // New format!
{
  if( !out_stream->good() )
  {
    T3cout<<"-Warning-T3R_node::save_binary:*Bad stream*"<<T3endl;
    return;
  }
  out_stream->write((const char*) &e, sizeof(double) );
  out_stream->write((const char*) &xs, sizeof(double) );
}

void T3R_node::load_binary(std::ifstream* in_stream ) // @@ old format
{
  if( !in_stream->good() )
  {
    T3cout<<"-Warning-T3R_node::load_binary: *Bad stream*"<<T3endl;
    return;
  }
  in_stream->read((char*) & e, sizeof(double) );
  in_stream->read((char*) & xs, sizeof(double) );
#ifdef debug
  T3cout<<"T3R_node::load_binary: E="<< e <<", XS="<< xs << T3endl;
#endif
}
std::ostream& operator<<(std::ostream& os, const T3R_node& inst)
{
  os << "E = " << inst.E() << ", xs = " << inst.XS() << T3endl;
  return os;
}

} // namespace t3
