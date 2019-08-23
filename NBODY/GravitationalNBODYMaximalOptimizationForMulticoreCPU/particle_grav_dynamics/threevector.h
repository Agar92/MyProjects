#ifndef THREEVECTOR_H
#define THREEVECTOR_H

#include <random>
#include "math.h"
#include <iostream>
using namespace std;

class ThreeVector
{
    public:
    ~ThreeVector() {}
    ThreeVector(double x = 0., double y = 0., double z = 0.) : fx(x), fy(y), fz(z){}
    // private:
    double fx;
    double fy;
    double fz;


  inline static double RND()
  {
        double r = static_cast <double> (rand())/(static_cast <double> (RAND_MAX));
        return r;
  }

  inline static ThreeVector RandomDirection()
  {
    /*default_random_engine generator1;
    uniform_real_distribution <double> distribution1(0.0, 1.0);
    double number1 = distribution1(generator1);*/
        double number1 = RND();
    //std::cout<<"RandomDirection() dir_number_1 = "<<number1<<std::endl;
        double z   = 2.0*number1 - 1.;  // z = cos(theta)
        double rho = sqrt((1.+z)*(1.-z)); // rho = sqrt(1-z*z)
        /*default_random_engine generator2;
        uniform_real_distribution <double> distribution2(0.0, 1.0);
        double number2 = distribution2(generator2);*/
        double number2 = RND();
        //std::cout<<"RandomDirection() dir_number_2 = "<<number2<<std::endl;
        double phi = M_PI*2*number2;
        return ThreeVector(rho*cos((float)phi), rho*sin((float)phi), z);
  }

  inline static long double RndMaxwell() // TheResult must be multiplied byTemperature
  {

        double number1 = RND();
        //std::cout<<"RndMaxwell() number_1 = "<<number1<<std::endl;
        double number2 = RND();
        //std::cout<<"RndMaxwell() number_1 = "<<number2<<std::endl;
        double number3 = RND();
        //std::cout<<"RndMaxwell() number_1 = "<<number3<<std::endl;
        static const double halfpi = M_PI / 2;
        const double R = cos((float)(number1 * halfpi));
        return -log((float)number2) - log((float)number3) * R * R;
  }
};
#endif // THREEVECTOR_H
