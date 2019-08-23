#pragma once
#ifndef T3TABULATEDFUNCTION_RW_H
#define T3TABULATEDFUNCTION_RW_H 1

#include <vector>
#include <iostream>
#include <fstream>
#include "globals.h"

namespace t3 {
//THIS CLASS:
//object with vector of x and vector of y with length Get_size().
template<class T1, class T2> class T3TabulatedFunction
{
public:
  //constructor(xvector, yvector)
  T3TabulatedFunction(std::vector<T1> arg = std::vector<T1>(),
                      std::vector<T2> val = std::vector<T2>()): _arg(arg), _val(val){}
  //returns x[n]
  T1 Get_arg(size_t n) const {return _arg.at(n);}
  //returns y[n]
  T2 Get_val(size_t n) const {return _val.at(n);}
  const std::vector<T1>& Get_arg() const {return _arg;}
  const std::vector<T2>& Get_val() const {return _val;}
  void Set_arg(size_t n, T1 arg) {_arg.at(n) = arg;}
  void Set_val(size_t n, T2 val) {_val.at(n) = val;}
  //returns the size of the xvector.
  size_t Get_size() const;
  //linear interpolation of y by any given x
  T2 Calc_val(T1 arg, T2 by_def = T2()) const;
  void Clear() {_arg.clear(); _val.clear();}
  //erase elements of xvector and yvector
  void Erase(size_t ind_min, size_t ind_after_max);
  //resize xvector and yvector
  void Resize(size_t size){_arg.resize(size); _val.resize(size);}
  //check that the xvector is in ascending order.
  bool Is_sorted() const;
  //returns min x
  T1 Get_min_arg() const {return _arg.front();} // FIXME undefined behaviour
  //returns max x
  T1 Get_max_arg() const {return _arg.back();} // FIXME undefined behaviour
  //adds (x,y) to xvector and yvector.
  void Add_point(T1 arg, T2 val, T1 delta = T1());
  void Push_back(T1 arg, T2 val) {_arg.push_back(arg); _val.push_back(val);}
  const std::vector<T2> Get_equidistant(T1 arg_min, T1 arg_max, size_t nsteps) const;
  void Save_to(std::ostream& os) const; // FIXME works only for simple T1, T2
  void Load_from(std::istream& is);    // FIXME works only for simple T1, T2
  T3TabulatedFunction<T1, T2>& operator=(const T3TabulatedFunction<T1, T2>& tab);
  const T3TabulatedFunction<T1, T2> operator*(const T1& t1) const;
  const T3TabulatedFunction<T1, T2> operator+(const T3TabulatedFunction<T1, T2>& t) const;
  bool operator==(const T3TabulatedFunction<T1, T2>& tab) const;
private:
  //Body
  std::vector<T1> _arg;
  std::vector<T2> _val;
};

template<class T1, class T2>
std::ostream& operator<<(std::ostream& os, const T3TabulatedFunction<T1, T2>& inst)
{
  for (size_t ind = 0; ind < inst.Get_size() ; ++ind)
  {
    os << '(' << inst.Get_arg(ind) << ','  << inst.Get_val(ind) << ')' << " ";
  }
  return os;
}
} // namespace t3
#endif
