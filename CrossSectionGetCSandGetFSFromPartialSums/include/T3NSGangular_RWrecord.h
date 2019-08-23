#pragma once
#ifndef T3NSGANGULAR_RWRECORD_H
#define T3NSGANGULAR_RWRECORD_H

#include <vector>
#include "T3NSGangular_RWnode.h"

namespace t3 {
class T3NSGangular_RWrecord: public std::vector<T3NSGangular_RWnode>
{
public:
  std::ofstream& save_binary( std::ofstream& out_stream) const;
  std::ifstream& load_binary( std::ifstream& in_stream);
  T3NSGangular_RWnode interpolate(T3double E) const;
  std::vector<T3NSGangular_RWnode> Get_vector_of_T3NSGangular_RWnodes(){return *this;}
  inline T3NSGangular_RWrecord& operator=(const T3NSGangular_RWrecord& rhs);
private:
};

std::ostream& operator<<(std::ostream& os, const T3NSGangular_RWrecord& inst);

inline T3NSGangular_RWrecord& T3NSGangular_RWrecord::operator=(const T3NSGangular_RWrecord & rhs)
{
  for(int i=0; i<rhs.size(); ++i) this->at(i)=rhs.at(i);
  return *this;
}
} // namespace t3
#endif
