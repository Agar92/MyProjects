#include <iostream>
#include <ostream>
#include <omp.h>
#include <vector>
#include <cmath>
#include <random>
#include <algorithm>
#include <iterator>
#include <unistd.h>
#include <array>
#include <iomanip>
#include <tuple>
#include <cstring>
#include <windows.h>

//using namespace std;

const int cuba=16;
const int cubn=cuba+cuba;

const double md=1875.6;
const double mt=44900.0;
const int N=1999;

template<class Floating> struct ThreeVector
{
  Floating x,y,z;
  ThreeVector():x(0.0),y(0.0),z(0.0){}
  ThreeVector(Floating _x, Floating _y, Floating _z):x(_x),y(_y),z(_z){}
  Floating R(){ return std::sqrt(x*x+y*y+z*z);}


  ThreeVector & operator=(ThreeVector<Floating> a)
  {
    this->x=a.x, this->y=a.y, this->z=a.z;
    return *this;
  }
  ThreeVector & operator+(ThreeVector<Floating> & a)
  {
    this->x+=a.x, this->y+=a.y, this->z+=a.z;
    return *this;
  }
  ThreeVector & operator*(Floating scalar)
  {
    this->x*=scalar, this->y*=scalar, this->z*=scalar;
    return *this;
  }
  friend std::ostream & operator<<(std::ostream & o, ThreeVector<Floating> & a)
  {
    std::cout<<a.x<<" "<<a.y<<" "<<a.z<<std::endl;
    return o;
  }
  void normalize()
  {
    Floating const r=(*this).R();
    this->x=this->x/r, this->y=this->y/r, this->z=this->z/r;
  }

};

template<class Floating> struct Particle
{
  ThreeVector<Floating> r;
  ThreeVector<Floating> p;
  int ir;
  int ip;//0 - d, 1 - Ti.
  int id;
  unsigned int rseed;
  Floating en;
  Floating de;

  std::minstd_rand rndengine;

  Particle():r(0.0,0.0,0.0), p(0.0,0.0,0.0), ir(2), ip(0), id(0), rseed(1), rndengine(rseed), en(0.0), de(0.0){}
  Particle(ThreeVector<Floating> _r, ThreeVector<Floating> _p, int _ir, int _ip, int _id, unsigned int _rseed,
           Floating _en, Floating _de):
      r(_r), p(_p), ir(_ir), ip(_ip), id(_id), rseed(_rseed), rndengine(rseed), en(_en), de(_de)
  {
    //rndengine.seed(rseed);
  }

  Floating vx(){ return p.x/p.R(); }
  Floating vy(){ return p.y/p.R(); }
  Floating vz(){ return p.z/p.R(); }
  Floating ix(){ return std::floor(r.x/a); }
  Floating jy(){ return std::floor(r.y/a); }
  Floating kz(){ return std::floor(r.z/a); }

  Floating rnd()
  {
    return static_cast<Floating>(rndengine())/std::minstd_rand().max();
  }

  static constexpr Floating a=1.e-5;
};

template<class Floating> class ParticleVector
{
public:
    using iterator = typename std::vector<Particle<Floating>>::iterator;

    Particle<Floating> & at(size_t ind){ return particles.at(ind); }
    Particle<Floating> & back(){ return particles.back(); }
    iterator begin(){ return particles.begin(); }
    iterator end(){ return particles.end(); }
    void resize(int size){ particles.resize(size); }
    size_t size(){ return particles.size(); }
private:
    std::vector<Particle<Floating>> particles;
};

template<class Floating> class DataHolder
{
public:
    size_t LIFE;
    size_t POSITION2;
    size_t POSITION0;
    size_t MAX_ELEMENT;
    static constexpr int Np=N+1;
    std::array<std::array<Floating, Np>, cubn> ArrayOfTheta;
    std::array<std::array<Floating, Np>, cubn> ArrayOfdE;
    std::array<std::array<Floating, Np>, cubn> ArrayOfE;
    int indZ;
    DataHolder():LIFE(0), POSITION2(0), POSITION0(0), MAX_ELEMENT(1), indZ(-16), ArrayOfTheta{}{}
    void InitParticle();
    void Compress();
    void Propagate();
    void React();
    void Inject();
    void FillHistogram();
    void FillHistogram_dE();
    void FillHistogram_E();
    size_t GetNumOfAliveParticles(){ return LIFE; }
    ThreeVector<Floating> interact(ThreeVector<Floating> & p, int i/*index*/);
private:
    ParticleVector<Floating> particles;
};

template<class Floating> void DataHolder<Floating>::InitParticle()
{
    const Floating a=Particle<Floating>::a;
    ThreeVector<Floating> r=ThreeVector<Floating>(0.5*a, 0.5*a, -cuba*a);
    Floating const Tls=10.0;
    Floating const pls=std::sqrt(Tls*(2*md+Tls));
    ThreeVector<Floating> p=ThreeVector<Floating>(0.0, 0.0, pls);
    Particle<Floating> newparticle=
            Particle<Floating>(r, p, 1, 0, MAX_ELEMENT-1, MAX_ELEMENT, md+Tls, 0.0);
    particles.resize(++LIFE);
    particles.back() = newparticle;
    ++MAX_ELEMENT;

    //std::cout<<"kz="<<particles.at(0).kz()<<std::endl;
    //sleep(3);
}

template<class Floating> void DataHolder<Floating>::FillHistogram()
{
    for(size_t i=0; i<Np; ++i)
    {
       for(size_t j=0; j<cubn; ++j) std::cout<<std::setw(12)<<std::setprecision(4)<<ArrayOfTheta.at(j).at(i)<<" ";
       std::cout<<std::endl;
    }
    std::cout<<std::endl;

    std::array<Floating, cubn> min;
    std::array<Floating, cubn> max;
    for(int j=0; j<cubn; ++j)
    {
      std::pair<Floating*, Floating*> mm=std::minmax_element(ArrayOfTheta.at(j).begin(),ArrayOfTheta.at(j).end());
      min.at(j)=*mm.first;
      max.at(j)=*mm.second;
    }
    std::cout<<"min: ";
    for(int j=0; j<cubn; ++j) std::cout<<min.at(j)<<" ";
    std::cout<<std::endl;
    std::cout<<"max: ";
    for(int j=0; j<cubn; ++j) std::cout<<max.at(j)<<" ";
    std::cout<<std::endl;

    const int M=20;
    const double eps=1.0e-7;
    std::array<Floating, cubn> deltaj;
    std::array<std::array<Floating, M>, cubn> InitTheta;
    std::array<std::array<Floating, M>, cubn> FinTheta;
    for(int j=0; j<cubn; ++j)
    {
      deltaj.at(j)=(max.at(j)-min.at(j))/M;
      for(int m=0; m<M; ++m)
      {
        if(m==0)   InitTheta.at(j).at(m)=min.at(j)-eps;
        else       InitTheta.at(j).at(m)=min.at(j)+m*deltaj.at(j);
        if(m==M-1) FinTheta.at(j).at(m)=max.at(j)+eps;
        else       FinTheta.at(j).at(m)=min.at(j)+(m+1)*deltaj.at(j);
      }
    }
    std::array<std::array<int, M>, cubn> NE;
    for(int j=0; j<cubn; ++j)
    {
      for(int m=0; m<M; ++m)
      {
        NE.at(j).at(m)=0;
      }
    }

    //std::cout<<"ggg"<<std::endl;
    for(int j=0; j<cubn; ++j)
    {
      for(int i=0; i<Np; ++i)
      {
        for(int m=0; m<M; ++m)
        {
          if( InitTheta.at(j).at(m) <= ArrayOfTheta.at(j).at(i)  &&  ArrayOfTheta.at(j).at(i) < FinTheta.at(j).at(m) )
              ++NE.at(j).at(m);
        }
      }
    }
    //std::cout<<"fff"<<std::endl;

    //print number of events:
    for(int m=0; m<M; ++m)
    {
      for(int j=0; j<cuba; ++j)
      {
        std::cout<<std::setw(3)<<NE.at(j).at(m)<<" ";
      }
      std::cout<<std::endl;
    }
    std::cout<<std::endl<<std::endl;
    for(int m=0; m<M; ++m)
    {
      for(int j=cuba; j<cubn; ++j)
      {
        std::cout<<std::setw(3)<<NE.at(j).at(m)<<" ";
      }
      std::cout<<std::endl;
    }
}

template<class Floating> void DataHolder<Floating>::FillHistogram_dE()
{
    for(size_t i=0; i<Np; ++i)
    {
       for(size_t j=0; j<cubn; ++j) std::cout<<std::setw(12)<<std::setprecision(4)<<ArrayOfdE.at(j).at(i)<<" ";
       std::cout<<std::endl;
    }
    std::cout<<std::endl;

    std::array<Floating, cubn> min;
    std::array<Floating, cubn> max;
    for(int j=0; j<cubn; ++j)
    {
      std::pair<Floating*, Floating*> mm=std::minmax_element(ArrayOfdE.at(j).begin(),ArrayOfdE.at(j).end());
      min.at(j)=*mm.first;
      max.at(j)=*mm.second;
    }
    std::cout<<"min: ";
    for(int j=0; j<cubn; ++j) std::cout<<min.at(j)<<" ";
    std::cout<<std::endl;
    std::cout<<"max: ";
    for(int j=0; j<cubn; ++j) std::cout<<max.at(j)<<" ";
    std::cout<<std::endl;

    const int M=100;
    const double eps=1.0e-7;
    std::array<Floating, cubn> deltaj;
    std::array<std::array<Floating, M>, cubn> InitdE;
    std::array<std::array<Floating, M>, cubn> FindE;
    for(int j=0; j<cubn; ++j)
    {
      deltaj.at(j)=(max.at(j)-min.at(j))/M;
      for(int m=0; m<M; ++m)
      {
        if(m==0)   InitdE.at(j).at(m)=min.at(j)-eps;
        else       InitdE.at(j).at(m)=min.at(j)+m*deltaj.at(j);
        if(m==M-1) FindE.at(j).at(m)=max.at(j)+eps;
        else       FindE.at(j).at(m)=min.at(j)+(m+1)*deltaj.at(j);
      }
    }
    std::array<std::array<int, M>, cubn> NE;
    for(int j=0; j<cubn; ++j)
    {
      for(int m=0; m<M; ++m)
      {
        NE.at(j).at(m)=0;
      }
    }

    //std::cout<<"ggg"<<std::endl;
    for(int j=0; j<cubn; ++j)
    {
      for(int i=0; i<Np; ++i)
      {
        for(int m=0; m<M; ++m)
        {
          if( InitdE.at(j).at(m) <= ArrayOfdE.at(j).at(i)  &&  ArrayOfdE.at(j).at(i) < FindE.at(j).at(m) )
              ++NE.at(j).at(m);
        }
      }
    }
    //std::cout<<"fff"<<std::endl;

    std::cout<<"NE dE:"<<std::endl;
    //print number of events:
    for(int m=0; m<M; ++m)
    {
      for(int j=0; j<cuba; ++j)
      {
        std::cout<<std::setw(3)<<NE.at(j).at(m)<<" ";
      }
      std::cout<<std::endl;
    }
    std::cout<<std::endl<<std::endl;
    for(int m=0; m<M; ++m)
    {
      for(int j=cuba; j<cubn; ++j)
      {
        std::cout<<std::setw(3)<<NE.at(j).at(m)<<" ";
      }
      std::cout<<std::endl;
    }
}

template<class Floating> void DataHolder<Floating>::FillHistogram_E()
{
    std::array<Floating, cubn> Eav;
    std::array<Floating, cubn> E2av;
    std::array<Floating, cubn> Eav2;
    std::array<Floating, cubn> D;
    for(int j=0; j<cubn; ++j)
    {
      Floating eav=0.0;
      Floating e2av=0.0;
      for(int i=0; i<Np; ++i)
      {
        eav+=ArrayOfE.at(j).at(i);
        e2av+=ArrayOfE.at(j).at(i) * ArrayOfE.at(j).at(i);
      }
      eav/=Np;
      e2av/=Np;
      Eav.at(j)=eav;
      E2av.at(j)=e2av;
      Eav2.at(j)=eav*eav;
      D.at(j)=E2av.at(j)-Eav2.at(j);
    }
    for(int j=0; j<cubn; ++j)
    {
        std::cout.precision(20);
        std::cout<<std::setw(10)<<j+1<<" "<<Eav.at(j)<<" "<<D.at(j)<<" "
                 <<E2av.at(j)<<" "<<Eav2.at(j)<<std::endl;
    }

}

template<class Floating> void DataHolder<Floating>::Compress()
{
    std::sort(particles.begin(), particles.end(),  [](Particle<Floating> & a, Particle<Floating> & b){ return a.ir>b.ir; });
    /*
    std::cout<<"LIFE="<<LIFE<<" id: ";
    for(size_t i=0; i<LIFE; ++i) std::cout<<particles.at(i).id<<" ";
    std::cout<<std::endl;
    std::cout<<"ir: ";
    for(size_t i=0; i<LIFE; ++i) std::cout<<particles.at(i).ir<<" ";
    std::cout<<std::endl;
    std::cout<<"z="<<particles.at(0).r.z<<std::endl;
    //sleep(1);
    */

    /*
    std::cout<<"rseed: ";
    for(size_t i=0; i<LIFE; ++i) std::cout<<particles.at(i).rseed<<" ";
    std::cout<<std::endl;
    std::cout<<"r: ";
    for(size_t i=0; i<LIFE; ++i) std::cout<<particles.at(i).r<<" ";
    std::cout<<std::endl;
    std::cout<<"p: ";
    for(size_t i=0; i<LIFE; ++i) std::cout<<particles.at(i).r<<" ";
    std::cout<<std::endl;
    //sleep(1);
    */

    //std::cout<<"FAIL COMP"<<std::endl;

    POSITION0=0;
    POSITION2=0;
    for(size_t i=0; i<LIFE; ++i)
    {
      if(particles.at(i).ir==0) ++POSITION0;
      if(particles.at(i).ir==2) ++POSITION2;
    }
    if(POSITION0)
    {
      particles.resize(LIFE-POSITION0);
      LIFE-=POSITION0;
    }

}

template<class Floating> void DataHolder<Floating>::Propagate()
{
    static Floating former_z;
    //std::cout<<"FAIL PROP LIFE="<<LIFE<<" size="<<particles.size()<<std::endl;
    for(size_t i=0; i<LIFE; ++i)
    {
      Particle<Floating> & particlei = particles.at(i);
      Floating const a=Particle<Floating>::a;
      Floating const da=1.0e-4*a;
      if(particlei.vx()<1.0e-10) particlei.p.x+=1.0e-10;
      if(particlei.vy()<1.0e-10) particlei.p.y+=1.0e-10;
      Floating l1x=(particlei.vx()>0.0?((particlei.ix()+1)*a-particlei.r.x):
                   (particlei.ix()*a-particlei.r.x))/particlei.vx()+da;
      Floating l1y=(particlei.vy()>0.0?((particlei.jy()+1)*a-particlei.r.y):
                   (particlei.jy()*a-particlei.r.y))/particlei.vy()+da;
      Floating l1z=(particlei.vz()>0.0?((particlei.kz()+1)*a-particlei.r.z):
                   (particlei.kz()*a-particlei.r.z))/particlei.vz()+da;
      /*
      std::cout<<"l1z="<<l1z<<" vz="<<particlei.vz()
               <<" a="<<a<<" z="<<particlei.r.z
               <<" kz="<<particlei.kz()<<" (kz+1)*a="<<(particlei.kz()+1)*a<<" d="
               <<((particlei.kz()+1)*a-particlei.r.z)<<" "
               <<((particlei.kz()+1)*a-particlei.r.z)/particlei.vz()
               <<std::endl;
      */

      if(indZ!=particlei.kz())
      {
        Floating d=(particlei.kz()+1)*a-particlei.r.z;
        //std::cout<<" indZ="<<indZ<<" kz="<<particlei.kz()<<std::endl;
        //std::cout<<"vz="<<particlei.vz()<<" z="<<particlei.r.z
        //         <<" (particlei.kz()+1)*a-particlei.r.z="
        //         <<d<<std::endl;
        indZ=particlei.kz();
        //sleep(3);
      }

      Floating const l1=std::min(std::min(l1x,l1y), std::min(l1y,l1z));
      //Floating const l0=1.0e-8;
      //Floating const dl=(2.0e-7-1.0e-8);
      //Floating const l=l0+dl*particlei.rnd();
      Floating r=std::log(particlei.rnd());
      Floating const l2=1.0e-7 * std::abs(r);

      Floating l=(l2<(l1-da)?l2:l1);
      size_t irc=(l2<(l1-da)?2:1);
      particlei.ir=irc;

      if(particlei.kz()==1)
      {
          //std::cout<<"kz="<<particlei.kz()<<" z="<<particlei.r.z
          //         <<std::endl;
          //sleep(2);
      }

      /*
      std::cout<<"| l1x="<<l1x<<" l1y="<<l1y<<" l1z="<<l1z
              <<" l2="<<l2<<" l="<<l<<" kz="<<particlei.kz()
              <<" z="<<particlei.r.z
              <<" ((particlei.kz()+1)*a-particlei.r.z)="
              <<((particlei.kz()+1)*a-particlei.r.z)
              <<" irc="<<irc<<" ir="<<particlei.ir<<std::endl;
      */

      Floating former_z=particlei.r.z;
      int indexZ=particlei.kz();
      //std::cout<<"HELL"<<std::endl;
      if(irc==1)
      {
        //std::cout<<"TOOP"<<std::endl;
        Floating thetai=std::acos(particlei.p.z/particlei.p.R());
        //std::cout<<"1="<<particlei.id<<" 2="<<indexZ+cuba<<" indexZ="<<indexZ<<std::endl;
        ArrayOfTheta.at(indexZ+cuba).at(particlei.id)=thetai;

        Floating const dE=particlei.de;
        ArrayOfdE.at(indexZ+cuba).at(particlei.id)=dE;

        ArrayOfE.at(indexZ+cuba).at(particlei.id)=particlei.en-md;

        //std::cout<<"ArrayOfdE.at(indexZ+cuba).at(particlei.id)="
        //        <<ArrayOfdE.at(indexZ+cuba).at(particlei.id)<<" dE="<<dE
        //        <<std::endl;

        //sleep(1);

        /*std::cout<<"irc="<<irc<<" id="<<particlei.id<<" thetai="<<thetai
                 <<" indexZ="<<indexZ<<" ic="<<indexZ+cuba
                 <<" z="<<particlei.r.z<<" z/a="<<particlei.r.z/a<<" a="<<a
                 <<" former_z="<<former_z<<std::endl;
        std::cout<<"-------------------------------------------"<<std::endl;
        sleep(1);
        */
      }


      /*
      if(irc==1)
      {
        std::cout<<"irc=1"<<" id="<<particlei.id<<" ir="<<particlei.ir
                 <<" ix="<<particlei.ix()<<" jy="<<particlei.jy()
                 <<" kz="<<particlei.kz()<<std::endl;
        //sleep(3);
      }
      */

      /*
      std::cout<<"l1="<<l1<<" l2="<<l2<<" r="<<r
               <<" ir="<<irc<<" l="<<l<<" l1x="<<l1x<<" l1y="<<l1y<<" l1z="<<l1z
               <<" vx="<<particlei.vx()<<" vy="<<particlei.vy()
               <<" vz="<<particlei.vz()
               <<" ix="<<particlei.ix()<<" jy="<<particlei.jy()
               <<" kz="<<particlei.kz()<<std::endl;
      Floating former_z=particlei.r.z;
      */
      particlei.r.x+=l*particlei.vx();
      particlei.r.y+=l*particlei.vy();
      particlei.r.z+=l*particlei.vz();
      //std::cout<<"fz="<<former_z<<" nz="<<particlei.r.z<<" vz="<<particlei.vz()
      //         <<" l="<<l<<" delta="<<(particlei.r.z-former_z)<<std::endl;
      if(particlei.ix()>=cuba || particlei.jy()>=cuba || particlei.kz()>=cuba
              || particlei.ix()<-cuba || particlei.jy()<-cuba || particlei.kz()<-cuba)
          particlei.ir=0;

      Floating const dEdx=2.0;//2*MeV/(g/cm^2).
      Floating rhotid2=3.91;//g/cm^3;
      Floating const dEdxfull=rhotid2*dEdx;
      if(particlei.ir==0)
      {
        Floating const m=(particlei.ip==0)?md:mt;
        particlei.de=particlei.en-m;
        particlei.en=m;
      }
      else
      {
          Floating loss=dEdxfull*l;
          particlei.de=loss;
          particlei.en-=loss;
      }

      if(indexZ!=particlei.kz() && irc!=1)
      {
        std::cout<<"l1z="<<l1z<<" l2="<<l2<<" l="<<l<<" kz="<<particlei.kz()
                 <<" indexZ="<<indexZ<<" irc="<<irc<<" ir="<<particlei.ir
                 <<" ((particlei.kz()+1)*a-particlei.r.z)/particlei.vz()+da="
                 <<(((indexZ+1)*a-former_z)/particlei.vz()+da)
                 <<" former_z+l2="<<(former_z+l2)
                 <<" dz="<<((indexZ+1)*a-former_z)<<std::endl;
        std::cout<<"HERE"<<std::endl;
        sleep(5);

      }
    }
    former_z=particles.at(0).r.z;

    //std::cout<<"FAIL END PROP"<<std::endl;
    //sleep(10);

}

template<class Floating>
ThreeVector<Floating> DataHolder<Floating>::interact(ThreeVector<Floating> & p, int i/*index*/)
{
  //std::cout<<"FAIL interact"<<std::endl;
  Floating r=particles.at(i).rnd();
  int isotope = (r<2.0/3.0)?0:1;//0 - d, 1 - Ti.
  Floating const m = md;
  Floating const M = (isotope==0?md:mt);
  Floating const z = 1;
  Floating const Z = (isotope==0?1:22);
  Floating const CTF = 0.885;
  Floating const alpha=1.0/137;
  Floating const ec=std::sqrt(alpha);
  Floating const me=0.511;
  Floating const a0=1.0/me/ec/ec;
  Floating const Elsinc=std::sqrt(m*m+p.R()*p.R());
  Floating const s=m*m+2*m*Elsinc+M*M;
  Floating mr=m*m/std::sqrt(s);
  Floating hc=200.0;
  Floating const Tls=Elsinc-m;
  Floating const pls2=Tls*(2*m+Tls);
  Floating const rat=m/M;
  Floating const rat2=rat*rat;
  Floating pcm2=pls2*M/(1.0+rat2+2*std::sqrt(rat2+pls2/M/M));
  Floating const pcm=std::sqrt(pcm2);
  Floating const aL=CTF*a0/std::sqrt(pow(z, 2.0/3.0)+pow(Z, 2.0/3.0));
  Floating c1=1.0/2/aL/pcm;
  Floating c12=c1*c1;
  Floating c2=alpha*z*Z;
  Floating c22=c2*c2;
  Floating AS=c12*(1.13+3.76*c22);
  Floating R=particles.at(i).rnd();
  Floating const rc1=4*pcm2*AS;
  Floating const Edisplace = (isotope==0?0.010:0.025);
  Floating const tmin=2*Edisplace*M;
  Floating const tmax=((4*pcm2<=tmin)?4*pcm2:tmin);
  Floating const rc2=tmax/(tmax+rc1);
  Floating const t=rc1*R*rc2/(1.0-R*rc2);//for multiple scattering (small |t|).
  //Floating const t=4*pcm2*R*AS/(1.0+AS-R);//this is if to integrate from 0 to 4*pcm^2 (small and big |t|).
  Floating const sinthetacm2=std::sqrt(t/4/pcm2);
  Floating const theta=2*std::asin(sinthetacm2);
  //LS->CM:
  Floating coef1=M/std::sqrt(s);
  ThreeVector<Floating> Pcm=p*coef1;
  Floating const xy=Pcm.x*Pcm.y, xz=Pcm.x*Pcm.z, yz=Pcm.x*Pcm.z, x2=Pcm.x*Pcm.x, y2=Pcm.y*Pcm.y, z2=Pcm.z*Pcm.z;
  ThreeVector<Floating> e1, e2;

  if(Pcm.x<Pcm.y)
  {
    if(Pcm.x<Pcm.z)
    {
      e1=ThreeVector<Floating>(0.0, Pcm.z, -Pcm.y);
      e2=ThreeVector<Floating>(y2+z2, -xy, -xz);
    }
    else
    {
      e1=ThreeVector<Floating>(Pcm.y, -Pcm.x, 0.0);
      e2=ThreeVector<Floating>(-xz, -yz, x2+y2);
    }
  }
  else
  {
    if(Pcm.y<Pcm.z)
    {
      e1=ThreeVector<Floating>(Pcm.z, 0.0, -Pcm.x);
      e2=ThreeVector<Floating>(-xy, x2+z2, -yz);
    }
    else
    {
      e1=ThreeVector<Floating>(Pcm.y, -Pcm.x, 0.0);
      e2=ThreeVector<Floating>(-xz, -yz, x2+y2);
    }
  }
  e1.normalize();
  e2.normalize();

  Floating const costhetacm=std::cos(theta);
  Floating const sinthetacm=std::sin(theta);
  Floating const Phi=2*M_PI*particles.at(i).rnd();
  Floating const cosPhi=std::cos(Phi);
  Floating const sinPhi=std::sin(Phi);

  Floating const sc=pcm*sinthetacm*cosPhi;
  Floating const ss=pcm*sinthetacm*sinPhi;
  ThreeVector<Floating> Pcmfin;
  Pcmfin=Pcm*costhetacm+e1*sc+e2*ss;
  //CM->LS:
  Floating coef=1.0/(Elsinc+M);
  ThreeVector<Floating> Vcm=p*coef;
  Floating const gamma=1.0/std::sqrt(1.-Vcm.R()*Vcm.R());
  Floating const Ecminc=std::sqrt(m*m+pcm2);
  ThreeVector<Floating> Plsfin=(Pcmfin+Vcm*Ecminc)*gamma;
  return Plsfin;
}

template<class Floating> void DataHolder<Floating>::React()
{
  //std::cout<<"FAIL REACT"<<std::endl;
  for(size_t i=0; i<POSITION2; ++i)
  {
    particles.at(i).p = interact(particles.at(i).p, i);
  }


}

template<class Floating> void DataHolder<Floating>::Inject()
{
  //std::cout<<"FAIL INJECT"<<std::endl;
  static int push=0;
  if(push<N)
  {
    InitParticle();
    ++push;
  }
  //std::cout<<"FAIL INJECT END"<<std::endl;

}

/*
std::tuple<int, double, std::string> func()
{
  return std::make_tuple(1, 7.0, "ser");
}
*/

/*
template<class Floating>
void test_assert()
{
  const int L=3;
  static_assert(L!=2, "The width is too small!");
}
*/

int factorial(int n)
{
  if(n==0) return 1;
  return factorial(n-1)*n;
}

typedef int (*IntIntFPointer)(int);


class Histogram
{
private:
  int * hist;
public:
  static constexpr int N=100;
  static constexpr double delta=1.0/N;
  Histogram()
  {
      hist=new int[N]{0};
  }
  ~Histogram()
  {
    delete [] hist;
  }
  void FillHistogram(double R)
  {
    if(R<0.0 || R>1.0) return;
    //std::cout<<"FillHistogram: bin="<<static_cast<int>(R/delta)<<" delta="<<delta
    //         <<" N="<<N<<std::endl;
    ++hist[static_cast<int>(R/delta)];
  }
  int & GetHist(int i)
  {
    return hist[i];
  }


};


void interpolate(double (*fun)(double), int Num)
{
  double Sum=0.0;
  for(int i=0; i<Num; ++i)
  {
    Sum+=fun(i*1.0/Num)*1.0/Num;
  }
}


int ** PrintMatrix(int **matrix)
{
  for(int i=0; i<3; ++i)
  {
    for(int j=0; j<3; ++j) std::cout<<matrix[i][j]<<" ";
    std::cout<<std::endl;
  }
  return matrix;
}

void FillMatrix(int **matrix)
{
  for(int i=0; i<3; ++i)
  {
    for(int j=0; j<3; ++j) matrix[i][j]=3*i+j;

  }
}

class Sum
{
  private:
    int sum;
  public:
    Sum():sum(0){}
    int operator () (int d){ return sum+=d; }
    int GetSum(){ return sum; }

};

using Floating = double;

int BinarySearch(Floating arr[], int N, Floating key)
{
  int l=0, r=N-1;
  int i=0;
  while(l<=r)
  {
    int mid=(l+r)/2;
    std::cout<<"i="<<i<<" l="<<l<<" r="<<r<<" mid="<<mid
             <<" arr["<<mid<<"]="<<arr[mid]<<" key="<<key<<std::endl;
    if(arr[mid]==key)
    {
      std::cout<<"@@@ "<<arr[mid]<<" "<<key<<std::endl;
      return mid;
    }
    else if(arr[mid]>key)
    {
      //std::cout<<">"<<std::endl;
      r=mid-1;
    }
    else if(arr[mid]<key)
    {
      //std::cout<<">"<<std::endl;
      l=mid+1;
    }
    sleep(1);
    ++i;
  }
  return l-1;
}

int main()
{
    int arr[10]={1,2,3,4,5,6,7,8,9,10};
    Floating farr[10]={0.0,0.04,0.1,0.25,0.4,0.43,0.5,0.8,0.98,1.0};
    //int index=BinarySearch(arr, 10, 18.0);
    int index=BinarySearch(farr, 10, 0.01);
    std::cout<<"index="<<index<<std::endl;
    exit(0);

    auto func = [](){ std::cout<<"ALL IS HELL THAT ENDS WELL"<<std::endl; };

    IntIntFPointer fun=factorial;

    std::cout<<"fun="<<(*fun)(5)<<std::endl;

    std::minstd_rand generator(1);
    Histogram * hist=new Histogram();
    #include <random>
    for(int i=0; i<1000000; ++i)
    {

      double R=static_cast<double>(generator())/generator.max();
      //std::cout<<"R="<<R<<std::endl;
      hist->FillHistogram(R);
      //std::cout<<"d"<<std::endl;
    }
    for(int i=0; i<Histogram::N; ++i) std::cout<<hist->GetHist(i)<<std::endl;

    int **matrix;
    matrix=new int * [3];
    for(int i=0; i<3; ++i)
    {
      matrix[i]=new int[3];
    }

    PrintMatrix(matrix);
    FillMatrix(matrix);
    PrintMatrix(matrix);

    for(int i=0; i<3; ++i)
    {
      delete [] matrix[i];
    }
    delete matrix;

    Sum sum;

    sum(0);
    std::cout<<sum.GetSum()<<std::endl;
    sum(10);
    std::cout<<sum.GetSum()<<std::endl;

    exit(0);

    //Sleep(2000);
    //func();

    //test_assert<double>();
    /*
    int a;
    double b;
    std::string c;
    std::tie(a, b, c)=func();
    std::cout<<"a="<<a<<" b="<<b<<" c="<<c<<std::endl;
    */
    std::cout << "Hello World!" << std::endl;
    std::cout<<"factorial="<<factorial(5)<<std::endl;
    exit(0);
    DataHolder<double> * dataholder = new DataHolder<double>();
    DataHolder<double> & d = *dataholder;
    d.InitParticle();

    //Particle<double> p =Particle<double>();
    //for(int i=0; i<10; ++i)  std::cout<<p.rnd()<<" ";


    while(d.GetNumOfAliveParticles())
    {
      d.Propagate();
      d.Compress();
      d.React();
      d.Inject();
      //std::cout<<"while:"<<" LIFE="<<d.GetNumOfAliveParticles()<<std::endl;
    }
    //std::cout<<"HERE"<<std::endl;
    d.FillHistogram();
    d.FillHistogram_dE();
    d.FillHistogram_E();
    delete dataholder;
    std::cout<<"END PRESS ANY KEY"<<std::endl;
    return 0;
}
