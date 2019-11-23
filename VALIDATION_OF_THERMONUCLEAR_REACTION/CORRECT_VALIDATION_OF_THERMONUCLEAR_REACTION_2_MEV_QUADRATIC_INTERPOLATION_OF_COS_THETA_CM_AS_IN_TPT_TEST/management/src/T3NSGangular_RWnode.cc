#include "T3NSGangular_RWnode.hh"

// #include "Randomize.hh"
// #include "T2Utility.hh"

namespace
{
  //\\//G4LorentzVector CalcXYZ(G4double p1, G4double p2, G4double p3, size_t i2);
  LorentzVector<FloatingType> CalcXYZ(T3double p1, T3double p2, T3double p3, size_t i2);
//   G4double RandomizeCost(const std::vector< G4LorentzVector >& LVP);
}

T3NSGangular_RWnode::T3NSGangular_RWnode(T3double e,
                                         const std::array<T3double, _num_point> src,
                                         T3double newpr, T3double newsl) : E(e),
                                         sl(newsl), pr(newpr), V(src)
{}

T3NSGangular_RWnode::T3NSGangular_RWnode(T3NSGangular_RWnode const &node): E(node.E),
                                         sl(node.sl), pr(node.pr), V(node.V){}

T3NSGangular_RWnode& T3NSGangular_RWnode::operator=(T3NSGangular_RWnode const &node)
{
  E = node.E;
  sl = node.sl;
  pr = node.pr;
  V = node.V;
  return *this;
}

T3double T3NSGangular_RWnode::Get_V(size_t point_num) const
{
  if (point_num > _num_point + 1)
  {
    T3cout << "T3NSGangular_RWnode:Get_V(" << point_num << ") point number out of range"
    << T3endl;
    return 0; // TODO return NaN
  }
  else if (0 == point_num)              return 0;
  else if (_num_point + 1 == point_num) return 1;
  else                                  return V.at(point_num - 1);
}

void T3NSGangular_RWnode::Set_V(size_t point_num, T3double val)
{
  if ( _num_point < point_num || 0 == point_num)
  {
    T3cout << "T3NSGangular_RWnode:Set_V(" << point_num << ") point number out of range"
           << T3endl;
    return;
  }
  else V.at(point_num - 1) = val;
}

void T3NSGangular_RWnode::Set_V(const std::array<T3double, _num_point> src)
{
  // TODO add 0 < x < 1; x_i < x_(i+1) checks
  V = src;
}

std::ofstream& T3NSGangular_RWnode::save_binary(std::ofstream& out_stream) const
{
  if( out_stream.good() )
  {
    unsigned int vector_size = _num_point; // (I(0)=0,I(16)=1)
    out_stream.write((const char*) &E, sizeof(T3double) );
    out_stream.write((const char*) V.data(), vector_size * sizeof(T3double));
    out_stream.write((const char*) &pr, sizeof(T3double) );
    out_stream.write((const char*) &sl, sizeof(T3double) );
  }
  else
  {
    T3cout<<"-Warning-T3NSGangular_RWnode::save_binary:*Bad stream*"<<T3endl;
  }
  return out_stream;
}

std::ifstream& T3NSGangular_RWnode::load_binary(std::ifstream& in_stream )
{
  if( in_stream.good() )
  {
    unsigned int vector_size = _num_point; // (I(0)=0,I(16)=1)
    in_stream.read((char*) & E, sizeof(T3double) );
    in_stream.read((char*) V.data(), vector_size * sizeof( T3double) );
    in_stream.read((char*) & pr, sizeof(T3double) );
    in_stream.read((char*) & sl, sizeof(T3double) );
  }
  else T3cout<<"-Warning-T3NSGangular_RWnode::load_binary: *Bad stream*"<<T3endl;

//\\//*****************************//
  //std::cout<<*this<<std::endl;
//\\//*****************************//  
  
  return in_stream;
}

//\\//G4LorentzVector T2NSGangular_RWnode::CalcXYZ(size_t point_num) const
LorentzVector<FloatingType> T3NSGangular_RWnode::CalcXYZ(size_t point_num) const
{
  if(0 < point_num && point_num <= _num_point)
    return ::CalcXYZ(Get_V(point_num - 1),Get_V(point_num),Get_V(point_num + 1),point_num);
  else
  {
    T3cout << "Warning: " <<  __func__ << "(" << point_num
           << "): bad arguement, must be in 0.." << _num_point << T3endl;
    //\\//return G4LorentzVector();
    return LorentzVector<FloatingType>();
  }
}

//\\//std::vector<G4LorentzVector> T2NSGangular_RWnode::CalcXYZ() const
std::vector<LorentzVector<FloatingType>> T3NSGangular_RWnode::CalcXYZ() const
{
  //\\//std::vector<G4LorentzVector> result;
  std::vector<LorentzVector<FloatingType>> result;
  result.reserve(_num_point);
  for(size_t ind = 1; ind <= _num_point; ++ind )
  {
    result.push_back(CalcXYZ(ind));
  }
  return result;
}


std::array<T3double, T3NSGangular_RWnode::_num_point> T3NSGangular_RWnode::isotropic()
{
  std::array<T3double, _num_point> array;
  static const T3double nbins = _num_point + 1;
  for (size_t ind = 0; ind < _num_point; ++ ind) array.at(ind) = T3double(ind+1)/nbins;
  return array;
}

namespace
{
  //\\//G4LorentzVector CalcXYZ(G4double p1, G4double p2, G4double p3,  size_t i2)
  LorentzVector<FloatingType> CalcXYZ(T3double p1, T3double p2, T3double p3,  size_t i2)
  {
    static const size_t   _num_point = T3NSGangular_RWnode::_num_point;
    static const T3double dCosT= 2. / ( _num_point + 1 );
    static const T3double tCosT= 2 * dCosT;
    //static const G4double eps= 0.001;  // calculation accuracy check (smaller->more print)
    static const T3double dop= 1.E-6;  // separate calculation accuracy
    static const T3double eql= 1.E-7;  // separate calculation accuracy
    if(i2 ==0 || i2 > _num_point)
      T3cout << "-Warning-T3NSGangular_RWnode::CXYZ:Bad i2=" << i2 << T3endl;
    if(i2==1 && p1!=0.)
      T3cout << "-Warning-T3NSGangular_RWnode::CXYZ: BadL P=" << p1 << T3endl;
    if(i2==_num_point && p3!=1.)
      T3cout << "-Warning-T3NSGangular_RWnode::CXYZ: BadR P=" << p3 << T3endl;
    if(p1==0. && p2==0. && p3==0.)
      T3cout << "-Warning-T3NSGangular_RWnode::CXYZ:AllZero p1="
            << p1 << ", p2=" << p2 << ", p3=" << p3 << T3endl;
    if(p1 > p2 || p2 > p3)
      T3cout << "-Warning-T3NSGangular_RWnode::CXYZ: NotMonotone p1="
            << p1 << ", p2=" << p2 << ", p3=" << p3 << T3endl;
    const long double ct = i2 * dCosT - 1.;         // Cost value in the central point
    //\\//G4LorentzVector res(0., 0., 0., p2);      // Prototype w/ central probability
    LorentzVector<FloatingType> res(0., 0., 0., p2);      // Prototype w/ central probability

    //\\//const long double r12 = p2-p1;
    //\\//const long double r23 = p3-p2;
    //\\//const long double r13 = p3-p1;
    
#ifdef debug
    T3cout<<"T3NSGangular_RWnode::CXYZ: p1="<< p1<<", r12="<< r12<<", r23="<< r23<<T3endl;
#endif
    
    //\\//if(std::fabs(r12) < eql && std::fabs(r23) < eql ) // p1 = p2 = p3
//\\//**************************************************//    
    if(std::fabs(p2-p1) < eql && std::fabs(p3-p2) < eql ) // p1 = p2 = p3
//\\//**************************************************//      
    {
#ifdef debug
      T3cout<<"T3NSGangular_RWnode::CXYZ:Cor,p1="<<p1<<",r12="<<r12<<",r23="<<r23<<T3endl;
#endif
      
      //\\//const long double b = tCosT/ r13;             // b = 2d / (p3 - p1)
      //\\//res.SetY(G4double(b));
      //\\//res.SetX(G4double(ct - b*(p1+p2+p3)/3));// a = ct - b * (p3+p2+p1) / 3
      //\\//return res;

//\\//**************************************************//
      p2=p1+eql;
      p3=p2+eql;
//\\//**************************************************//      
    }
    //\\//else if(std::fabs(r12) < eql)             // p1 & p2 coincide
//\\//**************************************************//
    else if(std::fabs(p2-p1) < eql)             // p1 & p2 coincide
//\\//**************************************************//        
    {
      //\\//const long double b = 3 * dCosT / r13 /2;
      //\\//res.SetX(G4double(ct + dCosT - b * p3));
      //\\//res.SetY(G4double(b));
      //\\//return res;
//\\//**************************************************//              
      p2=p1+eql;
//\\//**************************************************//              
    }
    //\\//else if(std::fabs(r23) < eql)             // p2 & p3 coincide
//\\//**************************************************//
    else if(std::fabs(p3-p2) < eql)             // p2 & p3 coincide
//\\//**************************************************//            
    {
      //\\//const long double b = 3 * dCosT / r12 /2;
      //\\//res.SetX(G4double(ct - dCosT - b * p1));
      //\\//res.SetY(G4double(b));
      //\\//return res;
//\\//**************************************************//
      p3=p2+eql;
//\\//**************************************************//              
    }

//\\//****************************************//    
    const long double r12 = p2-p1;
    const long double r23 = p3-p2;
    const long double r13 = p3-p1;
//\\//****************************************//
    
    const long double s12 = p1 + p2;
    const long double s23 = p2 + p3;
    const long double s13 = p1 + p3;
    const long double d12 = p1*p2*r12;
    const long double d23 = p2*p3*r23;
    const long double d13 = p1*p3*r13;
    const long double D   = d12 + d23 - d13;                 // Determinant
    if(std::fabs(D) < 1.E-13 * d13 || D == 0.) // if Det=0. -> use the 1-3 linear fit
    {
      long double a = ct;
      long double b = tCosT;
      if(r13 > dop)
      {
#ifdef debug
        T3cout<<"T3NSGangular_RWnode::CXYZ: *Lin* r13="<< r13 <<", r12="<< r12 <<T3endl;
#endif
        a -= dCosT*s13/r13;                    // a = ct - d * (p3+p1) / (p3-p1)
        b /= r13;                              // b = 2d / (p3-p1)
      }
      else if(r23 > r12 && r23 > dop)
      {
#ifdef debug
        T3cout<<"T3NSGangular_RWnode::CXYZ: *Lin* r23="<< r23 <<", r12="<< r12 <<T3endl;
#endif
        a -= dCosT*s23/r23;                     // a = ct - d * (p3+p2) / (p3-p2)
        b /= r23;                               // b = 2d / (p3-p2)
      }
      else if (r12 > dop)
      {
#ifdef debug
        T3cout<<"T3NSGangular_RWnode::CXYZ: *Lin* r12="<< r12 <<", r23="<< r23 <<T3endl;
#endif
        a -= dCosT*s12/r12;                     // a = ct - d * (p3+p2) / (p3-p2)
        b /= r12;                               // b = 2d / (p3-p2)
      }
      else if (r13 > 0.)
      {
#ifdef debug
        T3cout<<"T3NSGangular_RWnode::CXYZ:!,p1="<<p1<<",r13="<<r13<<",r12="<<r12<<T3endl;
#endif
        b /= r13;                              // b = 2d / (p3-p1)
        a -= b*(s13+p2)/3;                     // a = ct - b * (p3+p2+p1) / 3
      }
      else // @@ FIXME choose the right action -- > Done see above p1=p2=p3
      {
        T3cout<<"-War-T3NSGangular_RWnode::CXYZ: p1="<<p1<<",p2="<<p2<<",p3="<<p3<<T3endl;
//\\//*******************************************************//
//G4Eception is undefined out of Geant4:
//That is why i commented this line:
//\\//*******************************************************//
        ///////\\\\\\\\///////G4Exception("T2NSGangular_RWnode::CXYZ:","27",FatalException,"TPTCrash");
      }
      res.SetX(T3double(a));
      res.SetY(T3double(b));
#ifdef debug
      T3cout<<"T3NSGangular_RWnode::CXYZ: a="<<a<<",b="<<b<<", 1: "<<ct-dCosT<<"="<<a+p1*b
            <<", 2: "<<ct<<"="<<a+p2*b<<", 3: "<<ct+dCosT<<"="<<a+p3*b<<T3endl;
      if(std::fabs(ct-dCosT-a-p1*b) > eps ||
        std::fabs(ct+dCosT-a-p3*b) > eps ||
        std::fabs(ct-a-p2*b) > eps) //
        T3cout<<"-Warning-T3NSGangular_RWnode::CXYZ: a="<< a <<", b="<< b <<": "<<ct-dCosT
              <<"#"<<a+p1*b <<" || "<< ct <<"#"<< a+p2*b <<" || "<< ct+dCosT <<"#"<< a+p3*b
              <<", D="<< D <<", d13="<< d13 <<", p1="<<p1<<", p2="<<p2<<", p3="<<p3<<T3endl;
#endif
      return res;
    }

#ifdef debug
    T3cout << "T3NSGangular_RWnode::CXYZ: D = " << D << T3endl;
#endif
    const long double r = dCosT / D; // @@ FIXED SIGFPE (D = 0), but works under valgrind
    // happened when D = d13 = 0
    const long double a = ct + r * (d12-d23);
    const long double b = r * (r23*s23 - r12*s12);
    const long double c = r * (r12 - r23);
    res.SetX(T3double(a));
    res.SetY(T3double(b));
    res.SetZ(T3double(c));
#ifdef debug
    T3cout<<"T3NSGangular_RWnode::CXYZ: a="<<a<<",b="<<b<<",c="<<c<<", 1: "<<ct-dCosT<<"="
          << a+p1*(b+p1*c) << ", 2: " << ct <<"="<< a+p2*(b+p2*c) << ", 3: "<< ct+dCosT <<"="
          << a+p3*(b+p3*c) << T3endl;
    if( std::fabs( ct-dCosT-a-p1*(b+p1*c) ) > eps ||
        std::fabs( ct+dCosT-a-p3*(b+p3*c) ) > eps ||
        ( std::fabs( ct-a-p2*(b+p2*c) ) > eps && c > 0. ) )
      T3cout<<"-Warning-T3NSGangular_RWnode::CXYZ: a="<<a<<", b="<<b<<", c="<<c<<": "
            <<ct-dCosT<<"#"<<a+p1*(b+p1*c)<<" || "<<ct<<"#"<<a+p2*(b+p2*c)<<" || "
            <<ct+dCosT<<"#"<<a+p3*(b+p3*c)<<T3endl;
#endif
    return res;
  }
}



std::ostream& operator<<(std::ostream& os, const T3NSGangular_RWnode& inst)
{
  os << "Node e = " << inst.Get_E() << ", pr = " << inst.Get_pr() << ", sl = "
     << inst.Get_sl() << ", integrated distribution:";
  for(size_t ind =1; ind <= inst._num_point; ++ind) os << " " << inst.Get_V(ind);
  return os;
}

