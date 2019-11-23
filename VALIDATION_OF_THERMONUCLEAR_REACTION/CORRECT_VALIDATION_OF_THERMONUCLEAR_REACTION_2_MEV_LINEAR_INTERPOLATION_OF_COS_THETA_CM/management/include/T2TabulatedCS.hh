#pragma once
#ifndef T2TABULATEDCS_HH
#define T2TABULATEDCS_HH

#include "T2TabulatedFunction.hh"

class T2TabulatedCS
{
public:
  T2TabulatedCS(std::vector<G4double> arg = std::vector<G4double>(),
                std::vector<G4double> val = std::vector<G4double>()): _data(arg,val){}
  G4double Get_E(size_t n) const {return _data.Get_arg(n);}
  G4double Get_xs(size_t n) const {return _data.Get_val(n);}
  std::vector<G4double> Get_E() const {return _data.Get_arg();}
  std::vector<G4double> Get_xs() const {return _data.Get_val();}
  const T2TabulatedFunction<G4double, G4double>& Get_data() const {return _data;}
  void Set_E(size_t n, G4double E) {_data.Set_arg(n, E);}
  void Set_xs(size_t n, G4double xs) {_data.Set_val(n, xs);}
  size_t Get_size() const {return _data.Get_size();}
  G4double Calc_xs(G4double arg, G4double by_def = 0.) const
    {return _data.Calc_val(arg, by_def);}
  void Clear() {_data.Clear();}
  void Erase(size_t ind_min, size_t ind_after_max) {_data.Erase(ind_min, ind_after_max);}
  void Resize(size_t size){_data.Resize(size);}
  G4bool Is_sorted() const {return _data.Is_sorted();} // Monotone Energy check
  G4double Get_min_E() const {return Get_size() ? _data.Get_min_arg() : 0;}
  G4double Get_max_E() const {return Get_size() ? _data.Get_max_arg() : 0;}
  void Add_point(G4double arg, G4double val, G4double delta = 1e-11)
    {return _data.Add_point(arg, val, delta);}
  void Push_back(G4double arg, G4double val) {_data.Push_back(arg, val);}
  std::vector<G4double> Get_equidistant(G4double Emin, G4double Emax, size_t nsteps)
    const {return _data.Get_equidistant(Emin, Emax, nsteps);}
  void Save_to(std::ostream& os) const {_data.Save_to(os);}
  void Load_from(std::istream& is) {_data.Load_from(is);}
  inline T2TabulatedCS& operator=(const T2TabulatedCS& tab);
  inline const T2TabulatedCS operator*(const G4double mul) const;
  inline const T2TabulatedCS operator/(const G4double mul) const;
  inline const T2TabulatedCS operator+(const T2TabulatedCS& tab) const;
  inline const T2TabulatedCS operator-(const T2TabulatedCS& tab) const;
  inline G4bool operator==(const T2TabulatedCS& rval) const;
private:
  T2TabulatedFunction<G4double, G4double> _data;
};

inline std::ostream& operator<<(std::ostream& os, const T2TabulatedCS& inst)
{
  os << "(E, xs): " << inst.Get_data();
  return os;
}

inline T2TabulatedCS& T2TabulatedCS::operator=(const T2TabulatedCS& tab)
{
  _data = tab._data;
  return *this;
}

inline const T2TabulatedCS T2TabulatedCS::operator+(const T2TabulatedCS& tab) const
{
  T2TabulatedCS result;
  result._data = _data + tab._data;
  return result;
}

inline const T2TabulatedCS T2TabulatedCS::operator-(const T2TabulatedCS& tab) const
{
  return (*this) + (tab * (-1));
}

inline const T2TabulatedCS T2TabulatedCS::operator*(const G4double mul) const
{
  T2TabulatedCS result;
  result._data = _data * mul;
  return result;
}

inline const T2TabulatedCS T2TabulatedCS::operator/(const G4double den) const
{
  return (*this) * (1./den);
}

inline G4bool T2TabulatedCS::operator==(const T2TabulatedCS& rval) const
{
  G4bool temp = true;
  for (size_t ind = 0; ind < Get_size(); ++ind)
  {
    temp = temp && (fabs (Get_E(ind) - rval.Get_E(ind)) <= fabs(Get_E(ind)) * 1e-7 )
    && (fabs (Get_xs(ind) - rval.Get_xs(ind)) <= fabs(Get_xs(ind)) * 1e-7 );
  }
  return temp;
}

#endif
