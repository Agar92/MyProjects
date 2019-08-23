#pragma once
#ifndef T3TABULATEDCS_H
#define T3TABULATEDCS_H

#include <cmath>
#include "T3TabulatedFunction_RW.h"

namespace t3 {

class T3TabulatedCS
{
public:
  //initialize T3TabulatedFunction<double, double> _data:
  T3TabulatedCS(std::vector<double> arg = std::vector<double>(),
                std::vector<double> val = std::vector<double>()): _data(arg,val){}
  //returns x[n]
  double Get_E(size_t n) const {return _data.Get_arg(n);}
  //returns y[n]
  double Get_xs(size_t n) const {return _data.Get_val(n);}
  //returns the reference to xvector
  std::vector<double> Get_E() const {return _data.Get_arg();}
  //returns the reference to yvector
  std::vector<double> Get_xs() const {return _data.Get_val();}
  //returns the reference to T3TabulatedFunction<double, double>
  const T3TabulatedFunction<double, double>& Get_data() const {return _data;}
  //setx x[n]
  void Set_E(size_t n, double E) {_data.Set_arg(n, E);}
  //sets y[n]
  void Set_xs(size_t n, double xs) {_data.Set_val(n, xs);}
  //returns the size of xvector
  size_t Get_size() const {return _data.Get_size();}
  //linear interpolation of (arg,bydef) in T3TabulatedFunction<double, double>
  double Calc_xs(double arg, double by_def = 0.) const
    {return _data.Calc_val(arg, by_def);}
  void Clear() {_data.Clear();}
  void Erase(size_t ind_min, size_t ind_after_max) {_data.Erase(ind_min, ind_after_max);}
  void Resize(size_t size){_data.Resize(size);}
  bool Is_sorted() const {return _data.Is_sorted();} // Monotone Energy check
  double Get_min_E() const {return Get_size() ? _data.Get_min_arg() : 0;}
  double Get_max_E() const {return Get_size() ? _data.Get_max_arg() : 0;}
  //adds a point (x,y) to T3TabulatedFunction<double, double>
  void Add_point(double arg, double val, double delta = 1e-11)
    {return _data.Add_point(arg, val, delta);}
  //adds a point (x,y) to the end of T3TabulatedFunction<double, double>
  void Push_back(double arg, double val) {_data.Push_back(arg, val);}
  std::vector<double> Get_equidistant(double Emin, double Emax, size_t nsteps)
    const {return _data.Get_equidistant(Emin, Emax, nsteps);}
  //WRITE:
  void Save_to(std::ostream& os) const {_data.Save_to(os);}
  //READ:
  void Load_from(std::istream& is) {_data.Load_from(is);}
  inline T3TabulatedCS& operator=(const T3TabulatedCS& tab);
  inline const T3TabulatedCS operator*(const double mul) const;
  inline const T3TabulatedCS operator/(const double mul) const;
  inline const T3TabulatedCS operator+(const T3TabulatedCS& tab) const;
  inline const T3TabulatedCS operator-(const T3TabulatedCS& tab) const;
  inline bool operator==(const T3TabulatedCS& rval) const;
private:
  T3TabulatedFunction<double, double> _data;
};

inline std::ostream& operator<<(std::ostream& os, const T3TabulatedCS& inst)
{
  os << "(E, xs): " << inst.Get_data();
  return os;
}

inline T3TabulatedCS& T3TabulatedCS::operator=(const T3TabulatedCS& tab)
{
  _data = tab._data;
  return *this;
}

inline const T3TabulatedCS T3TabulatedCS::operator+(const T3TabulatedCS& tab) const
{
  T3TabulatedCS result;
  result._data = _data + tab._data;
  return result;
}

inline const T3TabulatedCS T3TabulatedCS::operator-(const T3TabulatedCS& tab) const
{
  return (*this) + (tab * (-1));
}

inline const T3TabulatedCS T3TabulatedCS::operator*(const double mul) const
{
  T3TabulatedCS result;
  result._data = _data * mul;
  return result;
}

inline const T3TabulatedCS T3TabulatedCS::operator/(const double den) const
{
  return (*this) * (1./den);
}

inline bool T3TabulatedCS::operator==(const T3TabulatedCS& rval) const
{
  bool temp = true;
  for (size_t ind = 0; ind < Get_size(); ++ind)
  {
    temp = temp && (fabs (Get_E(ind) - rval.Get_E(ind)) <= fabs(Get_E(ind)) * 1e-7 )
    && (fabs (Get_xs(ind) - rval.Get_xs(ind)) <= fabs(Get_xs(ind)) * 1e-7 );
  }
  return temp;
}
} // namespace t3
#endif
