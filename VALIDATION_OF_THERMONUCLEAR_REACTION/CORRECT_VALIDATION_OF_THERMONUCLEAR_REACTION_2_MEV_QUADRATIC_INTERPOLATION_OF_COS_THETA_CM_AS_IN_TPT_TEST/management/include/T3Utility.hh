#pragma once
#ifndef T3_UTILITY_HH
#define T3_UTILITY_HH

///\\\///#include "globals.hh"
///\\\///#include "Randomize.hh"

//\\//************************************88//
#include "T3Types.h"
//\\//************************************88//

#include <type_traits>
#include <limits>

namespace T3Utility
{
  // TODO allow usage for powers of 2 only
  template<size_t size> T3double* bin_search(T3double* ptr, T3double val);
  template<size_t size> const T3double* bin_search(const T3double* ptr, T3double val);
  // TODO const versions
  T3double* bin_search_r(size_t size, T3double* ptr, T3double val);
  T3double* bin_search(size_t size, T3double* ptr, T3double val);
  T3double const* bin_search(size_t size, T3double const* ptr, T3double val);
  //

  template<size_t size> inline T3double* bin_search(T3double* ptr, T3double val)
  {
    return val > ptr[size/2] ?  bin_search<size-size/2 - 1>(ptr+size/2 + 1, val) :
                                bin_search<size/2>(ptr, val);
  }

  template<size_t size> inline const T3double* bin_search(const T3double* ptr,
                                                          T3double val)
  {
    return val > ptr[size/2] ?  bin_search<size-size/2 - 1>(ptr+size/2 + 1, val) :
    bin_search<size/2>(ptr, val);
  }

  template<> inline const T3double* bin_search<0>(const T3double* ptr, T3double)
  {
    return ptr;
  }

  template<> inline T3double* bin_search<0>(T3double* ptr, T3double) {return ptr;}

  inline T3double* bin_search_r(size_t size, T3double* ptr, T3double val)
  {
    if (!size) return ptr;
    else return val > ptr[size/2] ?  bin_search_r(size-size/2 - 1, ptr+size/2 + 1, val) :
                                     bin_search_r(size/2, ptr, val);
  }

  inline T3double* bin_search(size_t size, T3double* ptr, T3double val)
  {
    while(size)
    {
      if (val > ptr[size/2])
      {
        ptr  += size/2 + 1;
        size -= size/2 + 1;
      }
      else
      {
        size /= 2;
      }
    }
    return ptr;
  }

  inline T3double const* bin_search(size_t size, T3double const* ptr, T3double val)
  {
    while(size)
    {
      if (val > ptr[size/2])
      {
        ptr  += size/2 + 1;
        size -= size/2 + 1;
      }
      else
      {
        size /= 2;
      }
    }
    return ptr;
  }
  
}

#endif


