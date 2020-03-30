#include <iostream>
#include <cmath>
#include "unistd.h"

namespace Realization
{

  double sin(double x)
  {
    const double eps = 1e-15;
    double s = 0;
    x = x - 2*M_PI*static_cast<int>(x/2/M_PI);
    double r = x;
    int n = 1;
    while (fabs(r) > eps)
    {
      s=s+r;
      n = n + 2;
      r=-r*x*x/(n*(n-1));
    }
    return s;
  }

  double cos(double x)
  {
    const double eps = 1e-15;
    double s = 1.0;
    x = x - 2*M_PI*static_cast<int>(x/2/M_PI);
    int n = 1;
    double r = 1.0; 
    while (fabs(r) > eps)
    {
      r=-r*x*x/(n*(n+1));
      n = n + 2;
      s=s+r;
    }
    return s;
  }

  long double exp(double x)
  {
    const double eps = 1e-15;
    long double s = 1.0;
    int n = 1;
    long double r = 1.0; 
    while (fabs(r) > eps)
    {
      r=r*x/n;
      ++n;
      s=s+r;
      //usleep(100000);
      std::cout<<"n="<<n<<" r="<<r<<" x="<<x<<std::endl;
    }
    std::cout<<"END3"<<std::endl;
    return s;
  }

  long double log(double x)
  {
    if(x<0.0) return -1.0;
    const double eps = 1e-12;
    long double s = 0.0;
    int n = 1;
    int Nt=0;
    if(x>2.0)
    {
      double temp=x;
      while(temp>2.0)
      {
        temp/=2;
        ++Nt;
      }
    }
    else Nt=0;
    //int Nt=static_cast<int>(x)/2;
    if(Nt>0) std::cout<<"YES "<<static_cast<int>(x)<<" "
                              <<Nt<<std::endl;

    std::cout<<"double(1<<(Nt))="<<double(1<<Nt)<<std::endl;
    x=x/double(1<<Nt);
    std::cout<<"1 x="<<x<<std::endl;
    x=x-1;
    std::cout<<"2 x="<<x<<std::endl;
    std::cout<<"x="<<x<<" Nt="<<Nt<<std::endl;
    long double r = -1.0; 
    long double delta=1.0;
    long double prev=r/n;
    if(!fabs(x-1.0)<1.0e-9)
    {
      while (fabs(r/n) > eps)
      {
        //std::cout<<"RA"<<" "<<fabs(r-1.0)<<" x="<<x<<std::endl;
        r=-r*x;
        s=s+r/n;
        ++n;
        //usleep(100000);
        delta=prev-r/n;
        //std::cout<<"n="<<n<<" r="<<r<<" r/n="<<r/n<<" s="<<s<<" delta="<<delta<<std::endl;
        prev=r/n;
      }
    }
    else
    {
      std::cout<<"s="<<s<<std::endl;
      //sleep(3);
      for(int h=1; h<(int)1.0e8; ++h)
      {
        r=-r*x;
        s=s+r/h;
        //std::cout<<"RB"<<" x="<<x<<" r="<<r<<" h="<<h<<" s="<<s<<std::endl;
      }
    }
    std::cout<<"END1"<<std::endl;
    //std::cout<<"s="<<s<<" Nt="<<Nt<<" log(2.0)="<<std::log(2.0)<<" log(x)="<<std::log(x+1)
    //         <<" 1<<Nt="<<(1<<Nt)<<std::endl;

    //std::cout<<"Nt="<<Nt<<" log(2)="<<std::log(2.0)<<" s="<<s
    //         <<" Mult="<<Nt*std::log(2.0)<<std::endl;
    if(Nt) s=s+Nt*Realization::log(2.0);
    std::cout<<"END2"<<" s="<<s<<std::endl;
    return s;
  }

  double powR(double x, double p)
  {
    std::cout<<"1R"<<std::endl;
    std::cout<<"Realization::log(x)="<<Realization::log(x)<<std::endl;
    std::cout<<"LLL"<<std::endl;
    //exit(0);
    const double res=Realization::exp(Realization::log(x)*p);
    std::cout<<"2R"<<std::endl;
    return res;
  }

  double sqrt1(double x)
  {
    if(x<0.0) return -1.0;
    const double eps = 1e-12;
    int n = 1;
    double X = 1.0;
    double Y = x/X;
    //std::cout<<"n="<<n<<" X="<<X<<" Y="<<Y<<" X-Y="<<X-Y<<std::endl;
    while(fabs(X-Y) > eps)
    {
      X=(X+Y)/2;
      ++n;
      Y=x/X;
      //usleep(100000);
      //std::cout.precision(20);
      //std::cout<<"n="<<n<<" X="<<X<<" Y="<<Y<<" X-Y="<<X-Y<<std::endl;
    }
    return X;
  }
  
  double sqrt2(double x)
  {
    if(x<0.0) return -1.0;
    const double eps = 1e-12;
    int n = 1;
    double L = 0.0;
    double H = x>1.0?x:1.0;
    double mid=0.0;
    while(fabs(H-L) > eps)
    {
      mid=(L+H)/2;
      ++n;
      std::cout.precision(20);
      //std::cout<<"n="<<n<<" mid="<<mid<<" mid^2="<<mid*mid<<" L="<<L<<" H="<<H<<" H-L="<<H-L<<std::endl;
      if(mid*mid>x)      H=mid;
      else if(mid*mid<x) L=mid;
      else               return mid;
      //usleep(100000);
    }
    return L;
  }

  double pow(double x, int n)
  {
    if(n==0)      return 1;
    else if(n==1) return x;
    else          return x*Realization::pow(x, n-1);
  }

  double rootn(double x, int N)
  {
    if(N%2==0 && x<0.0) return -1.0;
    const double eps = 1e-12;
    int n = 1;
    double X = 1.0;
    double Y = x/X;
    //std::cout<<"\n@@n="<<n<<" N="<<N<<" X="<<X<<" Y="<<Y<<" X-Y="<<X-Y<<std::endl;
    while(fabs(X-x/Realization::pow(X, N-1)) > eps)
    {
      //std::cout<<"T="<<(N-1)*X<<"   R="<<x/Realization::pow(X, N-1)
      //         <<" G="<<((N-1)*X+x/Realization::pow(X, N-1))<<std::endl;
      X=1.0/N*((N-1)*X+x/Realization::pow(X, N-1));
      //std::cout<<"HHH="<<X<<"   "<<1/N<<std::endl;
      ++n;
      //usleep(100000);
      //std::cout.precision(20);
      //std::cout<<"#n="<<n<<" X="<<X<<" N="<<N<<std::endl;
    }
    return X;
  }

}

int main()
{

  int a=-1;
  int b=a/2;
  std::cout<<"a="<<a<<" b="<<b<<std::endl; 

  const double x=1.1e-15;
  /*
  std::cout<<std::sin(x)<<"   "<<Realization::sin(x)<<std::endl;
  std::cout<<std::cos(x)<<"   "<<Realization::cos(x)<<std::endl;
  //std::cout<<std::exp(x)<<"   "<<Realization::exp(x)<<std::endl;
  std::cout<<std::log(x)<<"   "<<log(x)<<std::endl;//" er"<<std::log(x)-Realization::log(x)<<std::endl;
  std::cout<<std::sqrt(x)<<"   "<<Realization::sqrt1(x)<<"   "<<Realization::sqrt2(x)<<std::endl;
  */
  //std::cout<<pow(20, 1.0/1.0)<<"   "<<Realization::rootn(20.0, 1)<<std::endl;

  std::cout<<pow(10.0, 1.7)<<"   "<<Realization::powR(10.0, 1.7)<<std::endl;
	
  return 0;
}

