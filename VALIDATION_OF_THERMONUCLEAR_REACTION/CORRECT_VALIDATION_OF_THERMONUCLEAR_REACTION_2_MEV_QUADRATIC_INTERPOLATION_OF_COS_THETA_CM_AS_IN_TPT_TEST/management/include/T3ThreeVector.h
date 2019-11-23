#pragma once
#ifndef T3THREEVECTOR_H
#define T3THREEVECTOR_H

#include <cmath>

namespace t3
{

//************************************************************//
//These forward declarations are necessary for defining friend//
//functions in template class ThreeVector<T>:                 //
//************************************************************//
template <typename T>
class ThreeVector;

//scalar product:
template <typename T>
T operator*(const ThreeVector<T> & lhsvector, const ThreeVector<T> & rhsvector);
  
template <typename T>
ThreeVector<T> operator*(const T a, const ThreeVector<T> & rhsside);

template <typename T>
ThreeVector<T> operator*(const ThreeVector<T> & lhsside, const T a);

template <typename T>
ThreeVector<T> cross(const ThreeVector<T> & lhsvector, const ThreeVector<T> & rhsvector);  
//***************************************************************//

template <typename T>
class ThreeVector
{
private:
  T fx, fy, fz;
public:
//contsructors:
//\\//ThreeVector()=delete;
  ThreeVector():fx(0.0), fy(0.0), fz(0.0){}
  ThreeVector(T x, T y, T z):fx(x),fy(y),fz(z){}
  ThreeVector(ThreeVector const & vector);
//return components and vector lenghth:
//getters:
  T x() const {return fx;}
  T y() const {return fy;}
  T z() const {return fz;}  
  T R() const {return sqrt(fx*fx+fy*fy+fz*fz);}
//normalizes 3-vector:
  void Unit()
  {
    const T vector_length=sqrt(fx*fx+fy*fy+fz*fz);
    this->fx/=vector_length, this->fy/=vector_length, this->fz/=vector_length;
  }
//setters:
  void SetX(T X){fx=X;}
  void SetY(T Y){fy=Y;}
  void SetZ(T Z){fz=Z;}
//operator overloading:  
  ThreeVector & operator=(const ThreeVector & rhsvector);
  ThreeVector operator+(const ThreeVector & rhsvector);
  ThreeVector operator-(const ThreeVector & rhsvector);
  //unary operator, returns a negative copy of *this:
  ThreeVector operator-();
  ThreeVector & operator+=(const ThreeVector & rhsvector);
  ThreeVector & operator-=(const ThreeVector & rhsvector);
  //scalar product:
  friend T operator*<>(const ThreeVector & lhsvector, const ThreeVector & rhsvector);
  friend ThreeVector operator*<>(const T a, const ThreeVector & rhsside);
  friend ThreeVector operator*<>(const ThreeVector & lhsside, const T a);
  friend ThreeVector cross<>(const ThreeVector & lhsvector, const ThreeVector & rhsvector);  
  ThreeVector & operator*=(const T a);
  ThreeVector operator/(const T a);
  ThreeVector & operator/=(const T a);
};

template <typename T>
ThreeVector<T>::ThreeVector(ThreeVector const & vector):fx(vector.fx), fy(vector.fy), fz(vector.fz){}

template <typename T>
ThreeVector<T> & ThreeVector<T>::operator=(const ThreeVector & rhsvector)
{
  if(this==&rhsvector)
  {
    return *this;
  }
  this->fx=rhsvector.fx, this->fy=rhsvector.fy, this->fz=rhsvector.fz;
  return *this;
}

template <typename T>
ThreeVector<T> ThreeVector<T>::operator+(const ThreeVector & rhsvector)
{
  ThreeVector<T> result(this->fx+rhsvector.fx,this->fy+rhsvector.fy,this->fz+rhsvector.fz);
  return result;
}

template <typename T>
ThreeVector<T> ThreeVector<T>::operator-(const ThreeVector & rhsvector)
{
  ThreeVector<T> result(this->fx-rhsvector.fx,this->fy-rhsvector.fy,this->fz-rhsvector.fz);
  return result;
}

//unary operator, returns a negative copy of *this:
template <typename T>
ThreeVector<T> ThreeVector<T>::operator-()
{
  ThreeVector<T> result(-this->fx,-this->fy,-this->fz);
  return result;
}

template <typename T>
ThreeVector<T> & ThreeVector<T>::operator+=(const ThreeVector & rhsvector)
{
  this->fx+=rhsvector.fx, this->fy+=rhsvector.fy, this->fz+=rhsvector.fz;
  return *this;
}

template <typename T>
ThreeVector<T> & ThreeVector<T>::operator-=(const ThreeVector & rhsvector)
{
  this->fx-=rhsvector.fx, this->fy-=rhsvector.fy, this->fz-=rhsvector.fz;
  return *this;
}

//scalar product:
template <typename T>
T operator*(ThreeVector<T> const & lhsvector, const ThreeVector<T> & rhsvector)
{
  T result=lhsvector.fx*rhsvector.fx + lhsvector.fy*rhsvector.fy + lhsvector.fz*rhsvector.fz;
  return result;
}

template <typename T>
ThreeVector<T> operator*(const T a, const ThreeVector<T> & rhsside)
{
  ThreeVector<T> result(rhsside.fx*a, rhsside.fy*a, rhsside.fz*a);
  return result;
}

template <typename T>
ThreeVector<T> operator*(const ThreeVector<T> & lhsside, const T a)
{
  ThreeVector<T> result(lhsside.fx*a, lhsside.fy*a, lhsside.fz*a);
  return result;
}

template <typename T>
ThreeVector<T> cross(const ThreeVector<T> & lhsvector, const ThreeVector<T> & rhsvector)
{
  ThreeVector<T> result;
  result.fx=lhsvector.fy*rhsvector.fz-lhsvector.fz*rhsvector.fy;
  result.fy=lhsvector.fz*rhsvector.fx-lhsvector.fx*rhsvector.fz;
  result.fz=lhsvector.fx*rhsvector.fy-lhsvector.fy*rhsvector.fx;
  return result;
}

template <typename T>
ThreeVector<T> & ThreeVector<T>::operator*=(const T a)
{
  this->fx*=a, this->fy*=a, this->fz*=a;
  return *this;
}

template <typename T>
ThreeVector<T> ThreeVector<T>::operator/(const T a)
{
  ThreeVector<T> result(this->fx/a, this->fy/a, this->fz/a);
  return result;
}

template <typename T>
ThreeVector<T> & ThreeVector<T>::operator/=(const T a)
{
  this->fx/=a, this->fy/=a, this->fz/=a;
  return *this;
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const ThreeVector<T> & vector)
{
  os<<vector.x()<<" "<<vector.y()<<" "<<vector.z()<<" ";
  return os;
}
  
}//end of namespace t3.
#endif  //T3THREEVECTOR
