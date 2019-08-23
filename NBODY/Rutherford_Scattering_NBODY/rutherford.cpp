#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <algorithm>
#include <assert.h>
#include <math.h>

#include "stdlib.h"
#include "unistd.h"
#include <iostream>
#include <chrono>
#include <omp.h>

//handling Not a number exception:
#include <fenv.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <ostream>

///*
void handler(int sig)
{
  printf("Floating Point Exception\n");
  exit(0);
}
//*/
//

//build options pgcc: (-tp=haswell - for KNL):
//cmake3 ~/source/GPU/ework_recovery_v1_v2_gpu/ -DCMAKE_C_COMPILER=pgcc
//-DCMAKE_CXX_COMPILER=pgc++ -DCMAKE_C_FLAGS="-acc -Minfo=all -mcmodel=medium
//-ta=tesla:cc70 -tp=haswell -Minline" -DCMAKE_CXX_FLAGS="-acc -Minfo=all -mcmodel=medium
//-ta=tesla:cc70 -tp=haswell -Minline"

//either with "-fast" or without

//build options icpc: cmake3 ~/source/KNL/ework_recovery_v1_v2_KMV_random_generator/
//-DCMAKE_C_COMPILER=icc -DCMAKE_CXX_COMPILER=icpc -DCMAKE_CXX_FLAGS="-march=native
//-mtune=native -O3 -ipo16 -fp-model fast=2 -qopt-report=5 -mcmodel=large"

//#define ACCEL

#ifdef ACCEL
#include <accelmath.h>
#include <openacc.h>
#endif

//#define MAKEVRML//make vrml file
//#define DEB_NOPAR // switch off all the parallel pragmas and switch on all the diagnostic print
//#define EKCOR//switch on Ekin correction.
//#define DEB_RV
//#define DEB123
//#define DEB_VNEAR  // opens check of MAX_vx, MAX_dvx, MAX_NEAR
//#define DEB_KREL   // opens check of min, max of relative kinetic energy of i j particles in center mass
//#define DEB_KUR    // opens check of kinetic energy/potential energy of i j particles

//const float VELOCITY=2.e-2; // <=1
const float deltatime=1.0f;//1.e-5f;
const float deltatau=0.0001f;
//const float alpha=200.0f/137.0f;//197.46357f; // e^2/hc    hc=200 MeV * fm
const float alpha=197.46357/137; // e^2/hc    hc=200 MeV * fm

//512, 1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288, 1048576, 2097152, 4194304
const int numBodies=1024;//1024;//1024;//1024;//256;//2097152;//1048576;//65536;//16777216;//8388608;
const int max_loop=10000;//10;  // >=2
const float md=938.272;//1875.613;  // deuteron mass
const float Ep=0.1f;   // initial energy of deuterons = 0.1 MeV
float E_init;    // initial full energy
float Kf;  //kinetic energy at the previous step  
float Pf;  //potential energy at the previous step

unsigned int SEED1[numBodies];
unsigned int SEED2[numBodies];
unsigned int SEED3[numBodies];
unsigned int SEED4[numBodies];
unsigned int SEED5[numBodies];

#ifdef DEB_NOPAR
float ER1; // full energy in the beginning of the rotation() function
float ERB2; // full energy at the end of the rotation_back() function
float ELSCM;// full energy at ls->cm in calculation...() function

float MAX_S=0.0f;
int I_S=-1;
int J_S=-1;

float MAX_R=0.0f;
int I_R=-1;
int J_R=-1;

float MAX_R1=0.0f;
int I_R1=-1;
int J_R1=-1;

float MAX_R3=0.0f;
int I_R3=-1;
int J_R3=-1;

float MAX_R4=0.0f;
int I_R4=-1;
int J_R4=-1;

float MAX_RR2=0.0f;      // check dE=E2-E1, where E1 in the beginning of rotation() and E2 at the end of rotation().
int I_RR2=-1;
int J_RR2=-1;

float MAX_RRB2=0.0f;     // check dE=E2-E1, where E1 in the beginning of rotation_back() and E2 at the end of rotation_back().
int I_RRB2=-1;
int J_RRB2=-1;

float MAX_DERRB12=0.0f;  // check dE=E2-E1, where E1 in the beginning of rotation() and E2 at the end of rotation_back(). 
int I_DERRB12=-1;
int J_DERRB12=-1;

float MAX_C1=0.0f;       // check dE=E2-ER1, where ER1 in the beginning of rotation() and E2 int the mid of calculation...(). 
int I_C1=-1;
int J_C1=-1;

float MAX_LSCM=0.0f;       // check dE=E2-ER1, where ER1 in the beginning of rotation() and E2 int the mid of calculation...(). 
int I_LSCM=-1;
int J_LSCM=-1;

float MAX_CMLS=0.0f;       // check dE=E2-E1, where E1=Eiji[i][j] in the beginning of calculation...() and E2=ERB2 is the energy at the end of rotation_back(). 
int I_CMLS=-1;
int J_CMLS=-1;

float MAX_DA1=0.0f;       // in calculation...() before dRx[i], dRy[i]     
int I_DA1=-1;
int J_DA1=-1;

float MAX_DA2=0.0f;       // in calculation...() after the main i,j cycle        
int I_DA2=-1;
int J_DA2=-1;

float MAX_DA3=0.0f;       // in calculation...() after correction dRi, dVi        
int I_DA3=-1;
int J_DA3=-1;

float MAX_DIF1=0.0f;       // the difference betweeen r1cm and r1ls: d1=r1cm-r1ls
float ERR_DIF1=0.0f;
int I_DIF1=-1;
int J_DIF1=-1;

float MAX_DIF2=0.0f;       // the difference betweeen r2ls and r2cm: d2=r2ls-r2cm
float ERR_DIF2=0.0f;
int I_DIF2=-1;
int J_DIF2=-1;

float MAX_DIF3=0.0f;       // the difference betweeen dEk21ls and dEk21cm[i][j]=dEk21ls[i][j]-dEk21cm[i][j]
float ERR_DIF3=0.0f;
int I_DIF3=-1;
int J_DIF3=-1;

float MAX_NEAR_PAIR=1000000.0f; // the distance of the maximal closest approach of the particles in calculation...().
int I_NEAR_PAIR=-1;
int J_NEAR_PAIR=-1;

float MAX_U=0.0f;
int I_U=-1;
int J_U=-1;

float MAX_P=0.0f;
int I_P=-1;
int J_P=-1;

float MAX_E=0.0f;
float MIN_d=0.0f;
int I_E=-1;
int J_E=-1;

float MAX_D=0.0f;
int I_D=-1;
int J_D=-1;

float VCM=0.0f;
float UCM=0.0f;

float PING[numBodies];
#endif

#ifdef DEB_VNEAR
float MAX_vx=0.0f; 
float MAX_vy=0.0f; 
float MAX_vz=0.0f; 
                   
float MAX_dvx=0.0f;
float MAX_dvy=0.0f;
float MAX_dvz=0.0f;

float MAX_NEAR=1000000.0f;  // the distance of the maximal closest approach of the particles.
int I_NEAR=-1;
int J_NEAR=-1;
#endif

#ifdef DEB_KREL
float MAX_K=0.0;
int I_MAX_K=-1;
int J_MAX_K=-1;

float MIN_K=1000000.0;
int I_MIN_K=-1;
int J_MIN_K=-1;

float MAX_KUR=0.0;
int I_MAX_KUR=-1;
int J_MAX_KUR=-1;

float MIN_KUR=1000000.0;
int I_MIN_KUR=-1;
int J_MIN_KUR=-1;
#endif


int STEP=1;//0;

struct float3
{
  float x, y, z;
};

float3 operator+(const float3 v0, const float3 v1)
{
  float3 result{v0.x+v1.x, v0.y+v1.y, v0.z+v1.z};
  return result;
}
float3 operator-(const float3 v0, const float3 v1)
{
  float3 result{v0.x-v1.x, v0.y-v1.y, v0.z-v1.z};
  return result;
}
// this is wrong. Operators +=, -=, *= should receive the v0 argiment as a reference or the can't modify it.
//float3 & operator+=(const float3 & v0, const float3 v1)
//{
//  float3 result{v0.x+=v1.x, v0.y+=v1.y, v0.z+=v1.z};    /// ???
//  return result;
//}
/*
float3 operator+=(float3 v0, const float3 v1)
{
  float3 result{v0.x+=v1.x, v0.y+=v1.y, v0.z+=v1.z};    /// ???
  return result;
}
float3 operator-=(float3 v0, const float3 v1)
{
  float3 result{v0.x-=v1.x, v0.y-=v1.y, v0.z-=v1.z};    /// ???
  return result;
}
*/
float3 operator*(const float3 v, const float a)
{
  float3 result{v.x*a, v.y*a, v.z*a};
  return result;
}
float3 operator/(const float3 v, const float a)
{
  float3 result{v.x/a, v.y/a, v.z/a};
  return result;
}
float3 operator*(const float a, const float3 v)
{
  float3 result{v.x*a, v.y*a, v.z*a};
  return result;
}
float operator*(const float3 v0, const float3 v1)
{
  return v0.x*v1.x+v0.y*v1.y+v0.z*v1.z;
}
float3 normalize(float3 v)
{
  float d=sqrtf(v.x*v.x+v.y*v.y+v.z*v.z);
  float3 rt;
  rt.x=v.x/d;
  rt.y=v.y/d;
  rt.z=v.z/d;
  return rt;
}

float3 cross(float3 v0, float3 v1)
{
  float3 rt;
  rt.x = v0.y*v1.z-v0.z*v1.y;
  rt.y = v0.z*v1.x-v0.x*v1.z;
  rt.z = v0.x*v1.y-v0.y*v1.x;
  return rt;
}

struct double3
{
  double x, y, z;
};

double3 operator+(const double3 v0, const double3 v1)
{
  double3 result{v0.x+v1.x, v0.y+v1.y, v0.z+v1.z};
  return result;
}
double3 operator-(const double3 v0, const double3 v1)
{
  double3 result{v0.x-v1.x, v0.y-v1.y, v0.z-v1.z};
  return result;
}
double3 operator-(const double3 v)
{
  double3 result{-v.x, -v.y, -v.z};
  return result;
}
// this is wrong. Operators +=, -=, *= should receive the v0 argiment as a reference or the can't modify it.
/*
double3 operator+=(double3 v0, const double3 v1)
{
  double3 result{v0.x+=v1.x, v0.y+=v1.y, v0.z+=v1.z};    /// ???
  return result;
}
double3 operator-=(double3 v0, const double3 v1)
{
  double3 result{v0.x-=v1.x, v0.y-=v1.y, v0.z-=v1.z};    /// ???
  return result;
}
*/
double3 operator*(const double3 v, const double a)
{
  double3 result{v.x*a, v.y*a, v.z*a};
  return result;
}
double3 operator/(const double3 v, const double a)
{
  double3 result{v.x/a, v.y/a, v.z/a};
  return result;
}
double3 operator*(const double a, const double3 v)
{
  double3 result{v.x*a, v.y*a, v.z*a};
  return result;
}
float operator*(const double3 v0, const double3 v1)
{
  return v0.x*v1.x+v0.y*v1.y+v0.z*v1.z;
}
double3 normalize(double3 vec)
{
  double d=sqrt(vec.x*vec.x+vec.y*vec.y+vec.z*vec.z);
  double3 rt;
  rt.x=vec.x/d;
  rt.y=vec.y/d;
  rt.z=vec.z/d;
  return rt;
}
double3 cross(double3 v0, double3 v1)
{
  double3 rt;
  rt.x = v0.y*v1.z-v0.z*v1.y;
  rt.y = v0.z*v1.x-v0.x*v1.z;
  rt.z = v0.x*v1.y-v0.y*v1.x;
  return rt;
}

unsigned int Rand32(unsigned int xn)
{
  u_quad_t a=0x5DEECE66D;
  u_quad_t c=0xB;
  return (unsigned int)((a*xn+c) & 0xFFFFFFFF);
}

double rndv(unsigned int xn)
{
  return (double) xn / (double) 0x100000000LL;
}

double3 P_init;
//double3 dP;

struct Particle
{
  float x;
  float y;
  float z;
  float t;
  float vx;
  float vy;
  float vz;
  float m;
};

#ifdef DEB123
float S1[numBodies] __attribute__((aligned(64)));
float S2[numBodies] __attribute__((aligned(64)));
float S3[numBodies] __attribute__((aligned(64)));

float Sx[numBodies] __attribute__((aligned(64)));
float Sy[numBodies] __attribute__((aligned(64)));
float Sz[numBodies] __attribute__((aligned(64)));

float Sxn[numBodies] __attribute__((aligned(64)));
float Syn[numBodies] __attribute__((aligned(64)));
float Szn[numBodies] __attribute__((aligned(64)));
float Smv1[numBodies] __attribute__((aligned(64)));
#endif

#ifdef DEB_RV
float R1[numBodies] __attribute__((aligned(64)));
float R2[numBodies] __attribute__((aligned(64)));
float R3[numBodies] __attribute__((aligned(64)));
float VR1[numBodies] __attribute__((aligned(64)));
float VR2[numBodies] __attribute__((aligned(64)));
float VR3[numBodies] __attribute__((aligned(64)));
#endif

float dRxx[numBodies] __attribute__((aligned(64)));
float dRyy[numBodies] __attribute__((aligned(64)));
float dRzz[numBodies] __attribute__((aligned(64)));

float dVx[numBodies] __attribute__((aligned(64)));
float dVy[numBodies] __attribute__((aligned(64)));
float dVz[numBodies] __attribute__((aligned(64)));
// Uij
float U[numBodies][numBodies] __attribute__((aligned(64)));

#ifdef DEB_NOPAR
// dVij to check dVij=-dVji
float dPijx[numBodies][numBodies] __attribute__((aligned(64)));
float dPijy[numBodies][numBodies] __attribute__((aligned(64)));
float dPijz[numBodies][numBodies] __attribute__((aligned(64)));

float Eiji[numBodies][numBodies] __attribute__((aligned(64)));
float r1iji[numBodies][numBodies] __attribute__((aligned(64)));   // init r1 in ls
float r2cmij[numBodies][numBodies] __attribute__((aligned(64)));  // fin r2 in cm

float Eklsi[numBodies][numBodies] __attribute__((aligned(64)));  // init Eki+Ekj in ls
float Ekrli[numBodies][numBodies] __attribute__((aligned(64)));  // init Eki+Ekj in ls dv
float Ekrci[numBodies][numBodies] __attribute__((aligned(64)));  // init Eki+Ekj in cm dv
float Ekls[numBodies][numBodies] __attribute__((aligned(64)));   // fin Eki+Ekj in ls
float Ekcmi[numBodies][numBodies] __attribute__((aligned(64)));   // init Ek in cm
float dEk21cm[numBodies][numBodies] __attribute__((aligned(64)));   // fin  Ek - init Ek cm
float dEk21ls[numBodies][numBodies] __attribute__((aligned(64)));   // fin  Ek - init Ek ls
#endif

Particle p[numBodies] __attribute__((aligned(64)));

//allocation of the particle array p on the GPU:
//#pragma acc declare device_resident(p,dRxx,dRyy,dRzz,dVx,dVy,dVz,U,SEED1,SEED2,SEED3,SEED4,SEED5)
//#pragma acc declare create(p,dRxx,dRyy,dRzz,dVx,dVy,dVz,U,SEED1,SEED2,SEED3,SEED4,SEED5,S1,S2,S3,Sx,Sy,Sz,Sxn,Syn,Szn,Smv1)//,R1,R2,R3,VR1,VR2,VR3)
#ifdef ACCEL
#ifdef DEB123
#pragma acc declare create(S1,S2,S3,Sx,Sy,Sz,Sxn,Syn,Szn,Smv1)
#endif
#pragma acc declare device_resident(p,dRxx,dRyy,dRzz,dVx,dVy,dVz,U,SEED1,SEED2,SEED3,SEED4,SEED5)//,R1,R2,R3,VR1,VR2,VR3)
#endif

#ifdef MAKEVRML
struct VRML
{
  float x;
  float y;
  float z;
};
// VRML structures container array
VRML vrmlpos[numBodies][max_loop];
#endif
typedef std::pair<float3,float3> float6;
typedef std::pair<double3,double3> double6;

struct fourvec
{
  float3 rj;
  float3 vj;
  float3 ri;
  float3 vi;
};

//!!!THIS METHOD AND ROTATION MATRIXES R, Ri (inversed) MUST BE DOUBLE!!!
//USING float pivot, det=1.0; and float rotation MATRIXES
//GIVES AN ERROR: R*R^-1!=E
double inverseMat(double MatB[3][3])
{
  double pivot, det=1.0;           // changed float to double !!!
  int i, j, k;
  for(k=0; k<3; ++k)
  {
    pivot=MatB[k][k];
    det=det*pivot;
    //if(fabs(pivot)<1.e-9) return 0;
    if(fabs(pivot)<1.e-9) pivot=1.e-9;
    for(i=0; i<3; ++i)
      MatB[i][k]=-MatB[i][k]/pivot;
    for(i=0; i<3; ++i)
      if(i!=k)
        for(j=0; j<3; ++j)
          if(j!=k)
            MatB[i][j]=MatB[i][j]+MatB[k][j]*MatB[i][k];
    for(j=0; j<3; ++j) 
      MatB[k][j]=MatB[k][j]/pivot;
    MatB[k][k]=1./pivot;
    if(fabs(MatB[k][k])<1.e-9) MatB[k][k]=0.0;
  }
  return det;
}

inline float RND()
{
	float r = static_cast <float> (rand())/(static_cast <float> (RAND_MAX));
	return r;
}

inline float3 RandomDirection(unsigned int s1, unsigned int s2, int i)
{
  /*default_random_engine generator1;
  uniform_real_distribution <double> distribution1(0.0, 1.0);
  double number1 = distribution1(generator1);*/
  
	float number1 = rndv(s1);
    //std::cout<<"RandomDirection() dir_number_1 = "<<number1<<std::endl;
	float z   = 2.0*number1 - 1.;  // z = cos(theta)
	float rho = sqrtf((1.+z)*(1.-z)); // rho = sqrt(1-z*z)
#ifdef DEB123
  S1[i]=rho;
#endif
	/*default_random_engine generator2;
	uniform_real_distribution <double> distribution2(0.0, 1.0);
	double number2 = distribution2(generator2);*/
	float number2 = rndv(s2);
#ifdef DEB123
  S2[i]=number2;
#endif
	//std::cout<<"RandomDirection() dir_number_2 = "<<number2<<std::endl;
	float phi = M_PI*2.0*number2;
#ifdef DEB123
  S3[i]=phi;
#endif
  float3 result={rho*cosf(phi), rho*sinf(phi), z};
	///////return float3(rho*std::cos((float)phi), rho*std::sin((float)phi), z);
  return result;
}

//float6 rotation(float3 r1, float3 v1, float (&R)[3][3], float (&Ri)[3][3])
float6 rotation(float3 r1, float3 v1, double R[3][3], double Ri[3][3], int i, int j)
{  
  float x1=r1.x;
  float y1=r1.y;
  float z1=r1.z;
  float vx1=v1.x;
  float vy1=v1.y;
  float vz1=v1.z;

#ifdef DEB_NOPAR
  float r1dr=sqrt(x1*x1+y1*y1+z1*z1);
  float r1E=p[i].m*p[j].m/(p[i].m+p[j].m)*(vx1*vx1+vy1*vy1+vz1*vz1)/2+alpha/r1dr;
  ER1=r1E;
  //std::cout<<"r1E="<<r1E<<std::endl;
#endif
  
  float rx1, ry1, rz1;
  float rvx1, rvy1, rvz1;
  float r1ls=sqrt(r1.x*r1.x+r1.y*r1.y+r1.z*r1.z);     
  // here we calculate the orts of the coordinate system where Oxy is the plane where
  // the particles tracks and center mass in the CM CS lay.
  // ex should be || pi => || dV=Vi-Vj
  float3 ex={vx1, vy1, vz1};
  float df=sqrt(ex.x*ex.x+ex.y*ex.y+ex.z*ex.z);
  //normalize dV
  //ex=normalize(ex);
  ex.x/=df;
  ex.y/=df;
  ex.z/=df;
  float3 ey, ez;
  float3 ec={x1, y1, z1};
  //   [r x v]
  // __________
  //  |[r x v]|  
  //float3 n1=cross(ec,ex);            //??????????????????????????????????????????
  float3 n1;
  n1.x = ec.y*ex.z-ec.z*ex.y;
  n1.y = ec.z*ex.x-ec.x*ex.z;
  n1.z = ec.x*ex.y-ec.y*ex.x;
  
  float dist_n1=sqrt(n1.x*n1.x+n1.y*n1.y+n1.z*n1.z)/r1ls;
  // a vector perpendicular to r=ri-rj:
  float3 ec_p={y1-z1, -(x1+z1), x1+y1};
  //   [r_p x v]
  // ___________
  //  |[r_p x v]|
  //float3 n2=cross(ec_p,ex);           //????????????????????????????????
  float3 n2;
  n2.x = ec_p.y*ex.z-ec_p.z*ex.y;
  n2.y = ec_p.z*ex.x-ec_p.x*ex.z;
  n2.z = ec_p.x*ex.y-ec_p.y*ex.x;
  
  //float dist_n2=sqrt(n2.x*n2.x+n2.y*n2.y+n2.z*n2.z);
  //std::cout<<"dist_n1="<<dist_n1<<" dist_n2="<<dist_n2<<std::endl; 
  float eps=1.e-6;  
  if(dist_n1<eps)
  {
    //std::cout<<"****** CL"<<std::endl;
    float dn2=sqrt(n2.x*n2.x+n2.y*n2.y+n2.z*n2.z);
    n2.x/=dn2;
    n2.y/=dn2;
    n2.z/=dn2;
    
    //n2=normalize(n2);         //????????????????????????????????????????????
    ez.x=n2.x;
    ez.y=n2.y;
    ez.z=n2.z;       
  }
  else
  {
    float dn1=sqrt(n1.x*n1.x+n1.y*n1.y+n1.z*n1.z);
    n1.x/=dn1;
    n1.y/=dn1;
    n1.z/=dn1;

    //n1=normalize(n1);       //????????????????????????????????????????????
    ez.x=n1.x;
    ez.y=n1.y;
    ez.z=n1.z;
  }
  // ey=ez x ex
  //ey=cross(ez,ex);
  ey.x = ez.y*ex.z-ez.z*ex.y;
  ey.y = ez.z*ex.x-ez.x*ex.z;
  ey.z = ez.x*ex.y-ez.y*ex.x;
  // ex ey ez
  // here we fill in the rotation matrix R
  // R * (1 0 0) = (ex.x ex.y ex.z)
  R[0][0]=ex.x;
  R[1][0]=ex.y;
  R[2][0]=ex.z;
  // R * (0 1 0) = (ey.x ey.y ey.z)
  R[0][1]=ey.x;
  R[1][1]=ey.y;
  R[2][1]=ey.z;
  // R * (0 0 1) = (ez.x ez.y ez.z)
  R[0][2]=ez.x;
  R[1][2]=ez.y;
  R[2][2]=ez.z;   
  // here we reverse the rotation matrix to use it further
  // copy R to Ri
  for(int i1=0; i1<3; ++i1)
  {
    for(int j=0; j<3; ++j)
    {
      Ri[i1][j]=R[i1][j];
    }
  }
  inverseMat(Ri); // inversed rotation matrix R^-1

#ifdef DEB_NOPAR
/*
  if((i==263 && j==805) || (j==263 && i==805))
  {
    std::cout<<"*0.1 i="<<i<<" j="<<j<<" dist_n1="<<dist_n1<<" r1ls="<<r1ls<<std::endl;
    std::cout<<"ex: "<<ex.x<<" "<<ex.y<<" "<<ex.z<< " ey: "<<ey.x<<" "<<ey.y<<" "<<ey.z
             <<" ez: "<<ez.x<<" "<<ez.y<<" "<<ez.z<<std::endl;
    for(int k=0; k<3; ++k) std::cout<<R[k][0]<<" "<<R[k][1]<<" "<<R[k][2]<<std::endl;
    std::cout<<std::endl;
    for(int k=0; k<3; ++k) std::cout<<Ri[k][0]<<" "<<Ri[k][1]<<" "<<Ri[k][2]<<std::endl;
    std::cout<<std::endl;

    double Eu[3][3];    
    for(int t=0; t<3; ++t)
    {
      for(int pp=0; pp<3; ++pp)
      {
        float s1=0.0;
        for(int k=0; k<3; ++k)
        {
          s1+=R[t][k]*Ri[k][pp];
        }
        Eu[t][pp]=s1;
      }
    }
    for(int k=0; k<3; ++k) std::cout<<Eu[k][0]<<" "<<Eu[k][1]<<" "<<Eu[k][2]<<std::endl;
    std::cout<<std::endl;
    
  }
*/
#endif
       
  //Transform dr and dV from the initial CS system to the CS system where Oz is perp the plane (Oxy).
  //How to do it?
  //u1=R*(1 0 0) u2=R*(0 1 0) u3=R*(0 0 1)
  //R^-1*u1=(1 0 0) R^-1*u2=(0 1 0) R^-1*u3=(0 0 1)
  //=>
  //(rx1 ry1 rz1)=R^-1*(x1 y1 z1)
  rx1=Ri[0][0]*x1+Ri[0][1]*y1+Ri[0][2]*z1;
  ry1=Ri[1][0]*x1+Ri[1][1]*y1+Ri[1][2]*z1;
  rz1=Ri[2][0]*x1+Ri[2][1]*y1+Ri[2][2]*z1;
  //(rvx1 rvy1 rvz1)=R^-1*(vx1 vy1 vz1)
  rvx1=Ri[0][0]*vx1+Ri[0][1]*vy1+Ri[0][2]*vz1;
  rvy1=Ri[1][0]*vx1+Ri[1][1]*vy1+Ri[1][2]*vz1;
  rvz1=Ri[2][0]*vx1+Ri[2][1]*vy1+Ri[2][2]*vz1;
  float6 result;
  float3 vecr1={rx1, ry1, rz1};
  float3 vecv1={rvx1, rvy1, rvz1};
#ifdef DEB_NOPAR
  /*if((i==263 && j==805) || (j==263 && i==805))
  {
    std::cout<<"*1 i="<<i<<" j="<<j<<" r1: "<<vecr1.x<<" "<<vecr1.y<<" "<<vecr1.z
             <<" v1: "<<vecv1.x<<" "<<vecv1.y<<" "<<vecv1.z<<std::endl;
  }*/

  //float r2dr=sqrt(rx1*rx1+ry1*ry1+rz1*rz1);
  //float r2E=p[i].m*p[j].m/(p[i].m+p[j].m)*(rvx1*rvx1+rvy1*rvy1+rvz1*rvz1)/2+alpha/r2dr;
  //float dr2E=r2E-r1E;
  //std::cout<<"r2E="<<r2E<<std::endl;
  //if(fabs(dr2E)>MAX_RR2)
  //{
  //  MAX_RR2=fabs(dr2E);
  //  I_RR2=i;
  //  J_RR2=j;
  //}  
  float r2ls=sqrt(vecr1.x*vecr1.x+vecr1.y*vecr1.y+vecr1.z*vecr1.z);     
  //if(i==137 && j==645)
  //  std::cout<<"i="<<i<<" j="<<j<<" r2ls="<<r1ls<<std::endl;
#endif
  result.first=vecr1;
  result.second=vecv1;
  return result;
}
//fourvec calculation_i_j_interaction(float3 r1vec, float3 v1vec, int i, int j, int VR, int VPHI)
//fourvec calculation_i_j_interaction(float3 r1vec, float3 v1vec, int i, int j, int VR)
fourvec calculation_i_j_interaction(float3 r1vec, float3 v1vec, int i, int j)
{
  const float pi=3.14159265359;
  float rx1=r1vec.x;
  float ry1=r1vec.y;
  float rvx1=v1vec.x;
  
  //float rvy1=v1vec.y;
  //rvy1=0.0
  float mi=p[i].m;
  float mj=p[j].m;
  float mt=mi+mj;                               // total mass = m1 + m2
  float mr=mi*mj/mt;
  float ACCURACY=2.e-3;
  bool  MAX_APPROACH=true;
  if(fabs(rx1/ry1)<ACCURACY) MAX_APPROACH=false;
  // to polar coordinates
  float r1=sqrt(rx1*rx1+ry1*ry1);                   // polar radius of r1
#ifdef DEB_NOPAR
  float Ekcminit=p[i].m*p[j].m/(p[i].m+p[j].m)*rvx1*rvx1/2;
  Ekcmi[i][j]=Ekcminit;
  float r1cm_DIF1=sqrt(rx1*rx1+ry1*ry1);
  float dr1_DIF1=r1cm_DIF1-r1iji[i][j];
  //std::cout<<"dr1_DIF1="<<dr1_DIF1<<std::endl;
  if(fabs(dr1_DIF1)>MAX_DIF1)
  {
    MAX_DIF1=fabs(dr1_DIF1);
    ERR_DIF1=fabs(dr1_DIF1/r1cm_DIF1);
    I_DIF1=i;
    J_DIF1=j;
  }
  float dr1=sqrt(rx1*rx1+ry1*ry1);
  float E1=p[i].m*p[j].m/(p[i].m+p[j].m)*rvx1*rvx1/2+alpha/dr1;    
  if(r1<MAX_NEAR_PAIR)
  {
    MAX_NEAR_PAIR=r1;
    I_NEAR_PAIR=i;
    J_NEAR_PAIR=j;
  }
#endif   
  float phils;                                     // polar angle  of r1
  if(rx1>0 && ry1>=0) phils=atan(ry1/rx1);
  if(rx1<0 && ry1>=0) phils=pi+atan(ry1/rx1);
  if(rx1<0 && ry1<0)  phils=pi+atan(ry1/rx1);
  if(rx1>0 && ry1<0)  phils=2*pi+atan(ry1/rx1);  
  ///if(rx1==0.0)
  ///{
  ///  if(ry1>=0.0) phils=pi/2;
  ///  else         phils=pi*3.0/2.0;
  ///  //std::cout<<"rx1="<<rx1<<" phils="<<phils<<std::endl;
  ///}
  //if(rx1==0.0) std::cout<<"***ERROR: rx1=0!"<<std::endl;
  
  // transform from this cartesian CM coordinate system to the polar coordinate system
  // vector F={Fx,Fy}
  // Fx=Fr*cos(phi)-Fphi*sin(phi)
  // Fy=Fr*sin(phi)+Fphi*cos(phi)
  // =>
  // Fr  = Fx*cos(phi)+Fy*sin(phi)
  // Fphi=-Fx*sin(phi)+Fy*cos(phi)
  float v1r  = rvx1*cos(phils);     // v1r
  float v1phi=-rvx1*sin(phils);     // v1phi
  int VR;
  if(v1r>0.0)   VR=1;                               // ? =0 ?
  else          VR=-1;
  // E=K+U K=mu*V^2/2  U(r)=-alpha/r
  float M=fabs(r1*mr*v1phi);    // module of Oz projection of initial momentum M = r x p = m1*m2/(m1+m2)*r*vphi
  float E=M*M/(2*mr)/r1/r1+alpha/r1+mr*v1r*v1r/2; // initial E=M^2/(2*mr*r^2)+alpha/r+mr*v^2/2
  float dm;   // mr*alpha/M
  float den;  // sqrt(2*mr*E+(mr*alpha/M)^2)
  float da;   // angle to the asymptote
  if(M<1.e-10)  // vphi==0 (tracks are on the same line)
  {
    da=0.0;// vphi==0 (tracks are on the same line) // rho ->0 => da=arccos(alpha/((mu*V^2inf*rho)^2+alpha^2)^0.5)->0
#ifdef DEB_NOPAR
    std::cout<<"da=0 i="<<i<<" j="<<j<<std::endl;
#endif
  }
  else
  {
    dm=mr*alpha/M;                          // mu*alpha/M
    den=sqrt(2*mr*E+dm*dm);
    // angle of rotation da>0
    // we need a positive angle of rotation. so the integration will be from r1 to inf
    float Mr1=M/r1;
    float darg1=dm/den;
    float darg2=(Mr1+dm)/den; // does not give zero !!! r1 is not rmin
    //if(darg1>1.0) std::cout<<"***ERROR darg="<<darg1<<" > "<<darg1<<std::endl;
    //correction
    if(darg1>1.)
    {
      //std::cout<<"***darg1="<<darg1-1.0<<">1"<<std::endl;
      darg1=1.0f;
    }
    if(darg2>1.)
    {
      //std::cout<<"***darg2="<<darg2-1.0<<">1"<<std::endl;
      darg2=1.0f;
    }
    if(darg1<-1.)
    {
      //std::cout<<"***darg1="<<darg1+1.<<"<-1"<<std::endl;
      darg1=-1.0f;
    }
    if(darg2<-1.)
    {
      //std::cout<<"***darg2="<<darg2+1.<<"<-1"<<std::endl;
      darg2=-1.0f;
    }
    float da1=acos(darg1);
    float da2=acos(darg2);    // does not give zero !!! r1 is not rmin   
    da=da1-da2;// da    // does not give zero !!! r1 is not rmin        
  }
  // rmin
  float rmin=r1;
  if(MAX_APPROACH)
  {
    float ae=alpha/E; 
    rmin=ae/2+sqrt(ae*ae+2*(M*M/mr/E))/2;
  } 
  //phi1
  float phi1=phils-da;
  // coordinate system rotation angle=-da;
  float r2;
  //while(STEP<max_loop)
  //{ 
  float dr=v1r*deltatime;
  r2=r1+dr;
  //dphi=Mdt/[mu*r1*(r1+vr*dt)]
  float dph=(M/mr)*deltatime/r1/r2;          // dphi
  ////if(r2<=rmin)
  ////{
  ////  dph=(M/mr)*deltatime/rmin/rmin;
  ////  r2=r1;
  ////  VR=1;
  ////}
  if(MAX_APPROACH)
  {
    ////if(r2<=rmin && VR<0)
    if(r1+dr/2<=rmin && VR<0)
    {
      //dtau
      float dtau0=fabs(2*(r1-rmin)/v1r);
      ////float dtau=sqrt((r1-rmin)/mr/(alpha+rv2*r1*mr)/2);
      //y=sqrt(2Emr)*r-alpha*mr/sqrt(2Emr)
      float y_r1=sqrt(2*E*mr)*r1/alpha-sqrt(mr/2/E);
      float y_rmin=sqrt(2*E*mr)*rmin/alpha-sqrt(mr/2/E);
      float A2=mr/2/E+M*M/alpha/alpha;
      float rA2_r1=sqrt(fabs(y_r1*y_r1-A2));
      //std::cout<<"() A2="<<A2<<",y_rmin*y_rmin="<<y_rmin*y_rmin<<",y_rmin*y_rmin-A2="
      //         <<y_rmin*y_rmin-A2<<",y_rmin="<<y_rmin<<",rmin="<<rmin<<std::endl;
      float rrA2_rmin=y_rmin*y_rmin-A2; //=0
      //correction
      if(rrA2_rmin<0.0) rrA2_rmin=0.0f;
      float rA2_rmin=sqrt(rrA2_rmin);
      float lg_r1=log(y_r1+rA2_r1);
      //float log_rmin=log(y_rmin+rA2_rmin);
      float lg_rmin=log(y_rmin);
      float dtau=alpha/2/E*(rA2_r1-rA2_rmin)+sqrt(mr/2/E)*(lg_r1-lg_rmin);
      float rt;
      if(dtau<deltatau) rt=0.0f;
      else              rt=(deltatime/dtau-1.0);
      float dr2=(r1-rmin)*rt*rt;
      ////if(dtau!=0.0) rt=deltatime/dtau;
      ////else          rt=0.0;
      ////(deltatime/dtau-1)^2
      ////if(dtau>deltatime)
      ////{
      ////  dtau=deltatime;
      ////  dr2=0.0;
      ////}
      r2=rmin+dr2;
#ifdef DEB_NOPAR
      ////R=fabs(rx1/ry1)
      //std::cout<<"i="<<i<<",j="<<j<<",dtau="<<dtau<<",dtau0="<<dtau0<<",dr="<<r1-rmin<<",r1="<<r1
      //         <<",rmin="<<rmin<<",v1r="<<v1r<<",F="<<MAX_APPROACH<<",R="<<rx1/ry1<<",dr2="<<dr2<<std::endl;
      ////std::cout<<">>dr2="<<dr2<<std::endl;
      ////float dphdt1=M/(2*mr)*(1.0/r1/r1+1.0/rmin/rmin);
      ////float dphdt2=M/(2*mr)*(1.0/r2/r2+1.0/rmin/rmin);
      ////dph=dtau*dphdt1+(deltatime-dtau)*dphdt2;
#endif
      dph=(M/mr)*deltatime/rmin/rmin;
      VR=1;
    }
  }
  else
  {
    //std::cout<<"i="<<i<<",j="<<j<<",VR="<<VR<<",R="<<rx1/ry1<<",rx1="<<rx1<<",ry1="<<ry1<<std::endl;
    VR=0;
    ////dph=deltatime*rvx1/r1;
    dph=(M/mr)*deltatime/rmin/rmin;
    r2=r1;
  }  
  // V^2phi=M^2/(mu*r)^2
  float v2phi=M/mr/r2;
  float v2r2=2/mr*(E-alpha/r2)-v2phi*v2phi;
#ifdef DEB_NOPAR
  //float M2=fabs(mr*r2*v2phi);
  float dE=mr/2*(v2r2+v2phi*v2phi)+alpha/r2-E;
  //if(i==137 && j==645)
  //  std::cout<<"i="<<i<<" j="<<j<<" dE="<<dE<<std::endl;
  /*if(fabs(dE)>MAX_R)
  {
    std::cout<<"i="<<i<<" j="<<j<<" dE="<<dE<<std::endl;
    MAX_R=fabs(dE);
    I_R=i;
    J_R=j;   
  }*/
  
  //const float E2=M*M/(2*mr)/r2/r2+alpha/r2+mr*v2r2/2;
  //std::cout<<"E="<<E<<" E2-E="<<E2-E<<std::endl;
  //std::cout<<"v2r2="<<v2r2<<std::endl;
#endif
  float v2r=0.0f;
  if(v2r2>0.0f) v2r=VR*sqrtf(v2r2);
  float phi2=phi1+dph;                        // new angle: nphi=phi+dphi 
  // l0/l=m1/(mo+m1)
  // r0cm=-m1/mt*r
  // => r0cmx=-m1/mt*r*cos(phi)=m1/mt*r*cos(phi+pi)
  // => r0cmy=-m1/mt*r*sin(phi)=m1/mt*r*sin(phi+pi)
  float cx0=-mi/mt*r2*cos(phi2);
  float cy0=-mi/mt*r2*sin(phi2);
  float cx1=mj/mt*r2*cos(phi2);
  float cy1=mj/mt*r2*sin(phi2);
  // v2phi v2r to v2x v2y
  float v2x=v2r*cos(phi2)-v2phi*sin(phi2);
  float v2y=v2r*sin(phi2)+v2phi*cos(phi2);  
  // v0x v0y v0z
  float cvx0=-mi/mt*v2x;
  float cvy0=-mi/mt*v2y;
  //float cvz0=0.0;
  // v1x v1y v1z
  float cvx1=mj/mt*v2x;
  float cvy1=mj/mt*v2y;
  //float cvz1=0.0;
#ifdef DEB_NOPAR
  float dr3=sqrt((cx1-cx0)*(cx1-cx0)+(cy1-cy0)*(cy1-cy0));
  const float dE3=p[i].m*(cvx0*cvx0+cvy0*cvy0)/2+p[j].m*(cvx1*cvx1+cvy1*cvy1)/2+alpha/dr3-E;
  //std::cout<<"E="<<E<<" E3-E="<<dE3<<std::endl;
  if(fabs(dE3)>MAX_R3)
  {
    MAX_R3=fabs(dE3);
    I_R3=i;
    J_R3=j;
  }
#endif
  // to return to coordinates in the init LS coordinate system
  // a back rotation matrix
  // cos(phi)   -sin(phi)
  // sin(phi)    cos(phi)
  // may be used or one can add rangle to
  //float cx0=m_1/mt*r2*cos(phi2+pi+rangle);
  //float cy0=m_1/mt*r2*sin(phi2+pi+rangle);
  //float cx1=m0/mt*r2*cos(phi2+rangle);
  //float cy1=m0/mt*r2*sin(phi2+rangle);
  // it is right
  float cx0ls=cx0*cos(da)-cy0*sin(da);     // 
  float cy0ls=cx0*sin(da)+cy0*cos(da);     //                                  // z=0     
  float cx1ls=cx1*cos(da)-cy1*sin(da);     // 
  float cy1ls=cx1*sin(da)+cy1*cos(da);     //        
  float cvx0ls=cvx0*cos(da)-cvy0*sin(da);     // 
  float cvy0ls=cvx0*sin(da)+cvy0*cos(da);     //       
  float cvx1ls=cvx1*cos(da)-cvy1*sin(da);     // 
  float cvy1ls=cvx1*sin(da)+cvy1*cos(da);     // 
#ifdef DEB_NOPAR
  r2cmij[i][j]=sqrt((cx1ls-cx0ls)*(cx1ls-cx0ls)+(cy1ls-cy0ls)*(cy1ls-cy0ls));
  float dr4=sqrt((cx1ls-cx0ls)*(cx1ls-cx0ls)+(cy1ls-cy0ls)*(cy1ls-cy0ls));
  const float E4=p[i].m*(cvx0ls*cvx0ls+cvy0ls*cvy0ls)/2+p[j].m*(cvx1ls*cvx1ls+cvy1ls*cvy1ls)/2+alpha/dr4;
  const float dE4=E4-E;
  //std::cout<<"E="<<E<<" E1="<<E1<<" E4="<<E4<<" dE4="<<dE4<<" "<<" E4-E1="<<E4-E1<<std::endl; 
  if(fabs(dE4)>MAX_R4)
  {
    MAX_R4=fabs(dE4);
    I_R4=i;
    J_R4=j;
  }
  if(fabs(E4-E1)>MAX_R1)
  {
    MAX_R1=fabs(E4-E1);
    I_R1=i;
    J_R1=j;
  }
#endif 
  float3 vr0r={cx0ls, cy0ls, 0.0};
  float3 vr1r={cx1ls, cy1ls, 0.0};
  float3 vv0r={cvx0ls, cvy0ls, 0.0};
  float3 vv1r={cvx1ls, cvy1ls, 0.0};
  
#ifdef DEB_NOPAR
  float3 vcmr=(vv0r+vv1r)/2;
  float3 Vij=vv1r-vv0r;
  float Ekcmf=p[i].m*p[j].m/(p[i].m+p[j].m)*(Vij.x*Vij.x+Vij.y*Vij.y+Vij.z*Vij.z)/2;
  dEk21cm[i][j]=Ekcmf-Ekcmi[i][j];
  //std::cout<<"dEk21cm="<<dEk21cm[i][j]<<std::endl;
  //std::cout<<"Ekcmi="<<Ekcmi[i][j]<<" Ekcm="<<Ekcm[i][j]<<std::endl;
  
  //if(i==137 && j==645)
  //  std::cout<<"i="<<i<<" j="<<j<<" dE1="<<dE1<<" dR="<<distr2-r2<<" r1="<<r1<<" r2="<<r2<<" dr="<<dr<<" rmin="<<rmin
  //           <<" dvx="<<vcmr.x<<" dvy="<<vcmr.y<<" dvz="<<vcmr.z<<std::endl;
  /*if(fabs(dE)>MAX_R)
  {
    std::cout<<"i="<<i<<" j="<<j<<" dE="<<dE<<std::endl;
    MAX_R=fabs(dE);
    I_R=i;
    J_R=j;
    
  }*/
#endif
  
  //if(i==0 && j==2)
  //  std::cout<<"2. i="<<i<<" j="<<j<<" 2*vcmr.x-vv0r.x-vv1.x="<<vcmr.x-vv0r.x-vv1r.x<<" 2*vcmr.y-vv0r.y-vv1r.y="<<2*vcmr.y-vv0r.y-vv1r.y
  //                    <<" 2*vcmr.z-vv0r.z-vv1r.z="<<2*vcmr.z-vv0r.z-vv1r.z<<std::endl;

  // fabs(vv0r +vv1r=0 check
/*
//////////////////      
// !!!NO OUTPUT!!!
  if((i==263 && j==805) || (j==263 && i==805))
  {
  //if(fabs(vv0r.x+vv1r.x)>1.e-13)
  {
    std::cout<<"*c* vv1rx="<<vv1r.x<<" vv0rx="<<vv0r.x<<" i="<<i<<" j="<<j<<std::endl;
    //sleep(1);
  }
  //if(fabs(vv0r.y+vv1r.y)>1.e-13)
  {
    std::cout<<"*c* vv1ry="<<vv1r.y<<" vv0ry="<<vv0r.y<<" i="<<i<<" j="<<j<<std::endl;
    //sleep(1);
  }
  }
// !!!NO OUTPUT!!!
//////////////////
*/  
  fourvec result;
  result.rj=vr0r;
  result.vj=vv0r;
  result.ri=vr1r;
  result.vi=vv1r;
  return result;
}
/*
//We created this function to work with such positions of particles, where
//uij<<kij, but it gave a strange error when we tried to turn on all the parallel for and simd pragmas in (i,j).
//It was decided that only a for loop, inside which the same operations are performed in each thread, can be vectorised.
//What we vectorized was the for loop, in different iterations of which  very-very different task loads are (very heavy
//calculation_i_j_interaction() and very light field_i_j_interaction()). If we want to vectorize the for loop, all iterations
//of it should do the same work. That's why the error occured. 
//simple interaction if uij/kij < eps
fourvec field_i_j_interaction(float3 r1vec, float3 v1vec, int i, int j)
{
  float x1=r1vec.x;
  float y1=r1vec.y;
  float vx1=v1vec.x;
  float vy1=v1vec.y;
  float rdist1sqr=x1*x1+y1*y1;
  float rdist1=sqrt(rdist1sqr);
  float mi=p[i].m;
  float mj=p[j].m;
  float mt=mi+mj;                               // total mass = m1 + m2
  float mr=mi*mj/mt;
  float F=alpha/rdist1sqr;
  float Fx=F*x1/rdist1;
  float Fy=F*y1/rdist1;
  float dPx=Fx*deltatime;
  float dPy=Fy*deltatime;
  float dVx=dPx/mr;
  float dVy=dPy/mr;
  float vx2=vx1+dVx;
  float vy2=vy1+dVy;
  float x2=x1+(vx1+vx2)/2*deltatime;
  float y2=y1+(vy1+vy2)/2*deltatime;
  float rdist2sqr=x2*x2+y2*y2;
  float rdist2=sqrt(rdist2sqr);
  //energy conservation law: K1+P1=K2+P2
  float K2=mr*(vx2*vx2+vy2*vy2)/2;
  float K1=mr*(vx1*vx1+vy1*vy1)/2;
  float P1=alpha/rdist1;
  //correct r
  //float rcor=alpha/(K1+P1-K2);
  //float rcor=rdist2;
  //correction
  //float cor=rcor/rdist2;
  //x2*=cor;
  //y2*=cor;
  //coordinates of i and j particles
  float cx0=-mi/mt*x2;
  float cy0=-mi/mt*y2;
  float cx1=mj/mt*x2;
  float cy1=mj/mt*y2;
  //velocities of i and j particles
  float cvx0=-mi/mt*vx2;
  float cvy0=-mi/mt*vy2;
  float cvx1=mj/mt*vx2;
  float cvy1=mj/mt*vy2;
  float3 vr0r={cx0, cy0, 0.0};
  float3 vr1r={cx1, cy1, 0.0};
  float3 vv0r={cvx0, cvy0, 0.0};
  float3 vv1r={cvx1, cvy1, 0.0};
  fourvec result;
  result.rj=vr0r;
  result.vj=vv0r;
  result.ri=vr1r;
  result.vi=vv1r;
  return result;
}
*/

inline fourvec field_i_j_interaction(float3 r1vec, float3 v1vec, int i, int j)
{
  float x1=r1vec.x;
  float y1=r1vec.y;
  float z1=r1vec.z;
  float vx1=v1vec.x;
  float vy1=v1vec.y;
  float vz1=v1vec.z;
  float rdist1sqr=x1*x1+y1*y1+z1*z1;
  float rdist1=sqrt(rdist1sqr);
  float mi=p[i].m;
  float mj=p[j].m;
  float mt=mi+mj;                               // total mass = m1 + m2
  float mr=mi*mj/mt;
  float F=alpha/rdist1sqr;
  float Fx=F*x1/rdist1;
  float Fy=F*y1/rdist1;
  float Fz=F*z1/rdist1;
  float dPx=Fx*deltatime;
  float dPy=Fy*deltatime;
  float dPz=Fz*deltatime;
  float dVx=dPx/mr;
  float dVy=dPy/mr;
  float dVz=dPz/mr;
  float vx2=vx1+dVx;
  float vy2=vy1+dVy;
  float vz2=vz1+dVz;
  float x2=x1+(vx1+dVx/2)*deltatime;
  float y2=y1+(vy1+dVy/2)*deltatime;
  float z2=z1+(vz1+dVz/2)*deltatime;
  float rdist2sqr=x2*x2+y2*y2+z2*z2;
  float rdist2=sqrt(rdist2sqr);
  //energy conservation law: K1+P1=K2+P2
  float K2=mr*(vx2*vx2+vy2*vy2+vz2*vz2)/2;
  float K1=mr*(vx1*vx1+vy1*vy1+vz1*vz1)/2;
  float P1=alpha/rdist1;
  //correct r
  float rcor=alpha/(K1+P1-K2);
  //float rcor=rdist2;
  //correction
  float cor=rcor/rdist2;
  x2*=cor;
  y2*=cor;
  z2*=cor;
  //coordinates of i and j particles
  float cx0=-mi/mt*x2;
  float cy0=-mi/mt*y2;
  float cz0=-mi/mt*z2;
  float cx1=mj/mt*x2;
  float cy1=mj/mt*y2;
  float cz1=mj/mt*z2;
  //velocities of i and j particles
  float cvx0=-mi/mt*vx2;
  float cvy0=-mi/mt*vy2;
  float cvz0=-mi/mt*vz2;
  float cvx1=mj/mt*vx2;
  float cvy1=mj/mt*vy2;
  float cvz1=mj/mt*vz2;
  float3 vr0r={cx0, cy0, cz0};
  float3 vr1r={cx1, cy1, cz1};
  float3 vv0r={cvx0, cvy0, cvz0};
  float3 vv1r={cvx1, cvy1, cvz1};
  fourvec result;
  result.rj=vr0r;
  result.vj=vv0r;
  result.ri=vr1r;
  result.vi=vv1r;
  return result;
}

//fourvec rotation_back(float3 r0ls, float3 v0ls, float3 r1ls, float3 v1ls, float R[3][3])
fourvec rotation_back(fourvec rv01ls, double R[3][3], int i, int j)
{
  float3 r0ls=rv01ls.rj;
  float3 v0ls=rv01ls.vj;
  float3 r1ls=rv01ls.ri;
  float3 v1ls=rv01ls.vi;
  float cx0ls=r0ls.x;
  float cy0ls=r0ls.y;
  float cz0ls=r0ls.z;
        
  float cx1ls=r1ls.x;
  float cy1ls=r1ls.y;
  float cz1ls=r1ls.z;
        
  float cvx0ls=v0ls.x;
  float cvy0ls=v0ls.y;
  float cvz0ls=v0ls.z;
        
  float cvx1ls=v1ls.x;
  float cvy1ls=v1ls.y;
  float cvz1ls=v1ls.z;
  
#ifdef DEB_NOPAR
  float rb1dr=sqrt((cx1ls-cx0ls)*(cx1ls-cx0ls)+(cy1ls-cy0ls)*(cy1ls-cy0ls)+(cz1ls-cz0ls)*(cz1ls-cz0ls));
  float rb1E=md*(cvx1ls*cvx1ls+cvy1ls*cvy1ls+cvz1ls*cvz1ls)/2+md*(cvx0ls*cvx0ls+cvy0ls*cvy0ls+cvz0ls*cvz0ls)/2+alpha/rb1dr;
  //std::cout<<"rb1E="<<rb1E<<std::endl;
#endif  
  
  float brx0, bry0, brz0, brx1, bry1, brz1;
  float brvx0, brvy0, brvz0, brvx1, brvy1, brvz1;
  brx0=R[0][0]*cx0ls+R[0][1]*cy0ls+R[0][2]*cz0ls;
  bry0=R[1][0]*cx0ls+R[1][1]*cy0ls+R[1][2]*cz0ls;
  brz0=R[2][0]*cx0ls+R[2][1]*cy0ls+R[2][2]*cz0ls;
  // (brx1 bry1 brz1)=R*(cx1ls cy1ls cz1ls)
  brx1=R[0][0]*cx1ls+R[0][1]*cy1ls+R[0][2]*cz1ls;
  bry1=R[1][0]*cx1ls+R[1][1]*cy1ls+R[1][2]*cz1ls;
  brz1=R[2][0]*cx1ls+R[2][1]*cy1ls+R[2][2]*cz1ls;
  // (brvx0 brvy0 brvz0)=R*(cvx0ls cvy0ls cvz0ls)
  brvx0=R[0][0]*cvx0ls+R[0][1]*cvy0ls+R[0][2]*cvz0ls;
  brvy0=R[1][0]*cvx0ls+R[1][1]*cvy0ls+R[1][2]*cvz0ls;
  brvz0=R[2][0]*cvx0ls+R[2][1]*cvy0ls+R[2][2]*cvz0ls;
  // (brvx1 brvy1 brvz1)=R*(cvx1ls cvy1ls cvz1ls)
  brvx1=R[0][0]*cvx1ls+R[0][1]*cvy1ls+R[0][2]*cvz1ls;
  brvy1=R[1][0]*cvx1ls+R[1][1]*cvy1ls+R[1][2]*cvz1ls;
  brvz1=R[2][0]*cvx1ls+R[2][1]*cvy1ls+R[2][2]*cvz1ls;
#ifdef DEB_NOPAR
  float rb2dr=sqrt((brx1-brx0)*(brx1-brx0)+(bry1-bry0)*(bry1-bry0)+(brz1-brz0)*(brz1-brz0));
  float rb2E=md*(brvx1*brvx1+brvy1*brvy1+brvz1*brvz1)/2+md*(brvx0*brvx0+brvy0*brvy0+brvz0*brvz0)/2+alpha/rb1dr;
  float drb2E=rb2E-rb1E;
  //std::cout<<"rb2E="<<rb2E<<std::endl;
  if(fabs(drb2E)>MAX_RRB2)
  {
    MAX_RRB2=fabs(drb2E);
    I_RRB2=i;
    J_RRB2=j;
  }
  ERB2=rb2E;
  float dE12rbr=ERB2-ER1;
  if(fabs(dE12rbr)>MAX_DERRB12)
  {
    MAX_DERRB12=dE12rbr;
    I_DERRB12=i;
    J_DERRB12=j;
  }
#endif
  float3 br0={brx0, bry0, brz0};
  float3 br1={brx1, bry1, brz1};
  float3 bv0={brvx0, brvy0, brvz0};
  float3 bv1={brvx1, brvy1, brvz1};
    
  fourvec result;
  result.rj=br0;
  result.vj=bv0;
  result.ri=br1;
  result.vi=bv1;
  return result;
}
//================================================================
void interaction()
{
  float EPS_INTERACTION=1.e-27;
  //dP={0.0, 0.0, 0.0}; 
  ///const double pi=3.14159265359; 
#ifdef DEB_VNEAR
  MAX_vx=0.0f; 
  MAX_vy=0.0f; 
  MAX_vz=0.0f; 
               
  MAX_dvx=0.0f;
  MAX_dvy=0.0f;
  MAX_dvz=0.0f;
#endif

#ifdef DEB_KREL
  MAX_K=0.0;
  MIN_K=1000000.0;

  MAX_KUR=0.0;
  MIN_KUR=1000000.0;
#endif

  float Ke=0.0;
#ifdef ACCEL
#pragma acc parallel loop copy(Ke) present(p)
#else
//#pragma omp parallel for schedule(dynamic)
#pragma simd reduction(+:Ke)
#endif
  for(int i=0; i<numBodies; ++i)
  {
    Ke+=p[i].m*(p[i].vx*p[i].vx+p[i].vy*p[i].vy+p[i].vz*p[i].vz)/2;
  }
  
/*  
  std::cout<<"@K="<<Ke<<std::endl;
#ifdef ACCEL
#pragma acc update host(p)
#endif
  //CHeck3
  for(int i=0; i<numBodies; ++i) std::cout<<"3. r: "<<p[i].x<<" "<<p[i].y<<" "<<p[i].z<<" v: "<<p[i].vx<<" "<<p[i].vy<<" "<<p[i].vz<<std::endl;
*/
#ifndef DEB_NOPAR
#ifdef ACCEL
#pragma acc parallel loop present(p)
#else
  //#pragma omp parallel for schedule(dynamic)
  //#pragma distribute_point
#endif
#endif
  for(int i=0; i<numBodies; ++i) //if(i==1)
  { 
    /*
    dRx[i]=0.0;
    dRy[i]=0.0;
    dRz[i]=0.0;
    dVx[i]=0.0;
    dVy[i]=0.0;
    dVz[i]=0.0;
    */
    //float Eu[3][3]; // unit matrix     
    //float3 ri={p[i].x,p[i].y,p[i].z};
    //float3 vi={p[i].vx,p[i].vy,p[i].vz};
    // instead of dPi write dVi !
    float dRxi=0.0,dRyi=0.0,dRzi=0.0;
    float dVxi=0.0,dVyi=0.0,dVzi=0.0;
    float mi=p[i].m;
#ifdef ACCEL
#pragma acc loop vector reduction(+:dRxi,dRyi,dRzi,dVxi,dVyi,dVzi)
#else
     //#pragma simd reduction(+:dRxi,dRyi,dRzi,dVxi,dVyi,dVzi) //private(dRxi,dRyi,dRzi,dVxi,dVyi,dVzi)//,ri,vi)
     //#pragma vector aligned
#endif
    for(int j=0; j<numBodies; ++j) if(i!=j)
    {
      //bool FLAG=true;
      //Right is to put the matrixes in j for cycle!!!
      double R[3][3];  // rotation matrix        
      double Ri[3][3]; // inverse rotation matrix 
    #ifdef DEB_NOPAR
      float Ekii=p[i].m*(p[i].vx*p[i].vx+p[i].vy*p[i].vy+p[i].vz*p[i].vz)/2;
      float Ekji=p[j].m*(p[j].vx*p[j].vx+p[j].vy*p[j].vy+p[j].vz*p[j].vz)/2;
      float distiji=sqrtf((p[i].x-p[j].x)*(p[i].x-p[j].x)+(p[i].y-p[j].y)*(p[i].y-p[j].y)+(p[i].z-p[j].z)*(p[i].z-p[j].z));
      //std::cout<<"d="<<distiji<<" i="<<i<<" j="<<j<<std::endl;
      float Epiji=alpha/distiji;
      Eiji[i][j]=Ekii+Ekji+Epiji;
      Eklsi[i][j]=Ekii+Ekji;
    #endif   
      //std::cout<<"Init: i="<<i<<" j="<<j<<" Eiji="<<Eiji<<std::endl;   
      //int VR;    // sign of the initial vr projection    // put the definitiona of VR in calculation_i_j(...).
      //int VPHI;  // sign of the initial vhi projection      
      //float dRxj=0.0,dRyj=0.0,dRzj=0.0;
      //float dVxj=0.0,dVyj=0.0,dVzj=0.0;
      float3 ri={p[i].x,p[i].y,p[i].z};
      float3 vi={p[i].vx,p[i].vy,p[i].vz};
      float3 rj={p[j].x,p[j].y,p[j].z};
      float3 vj={p[j].vx,p[j].vy,p[j].vz};
      
      ///float6 rvcm=getcmrv(i, r1, v1, j, r0, v0);
      ///float3 rcm=rvcm.first; 
      ///float3 vcm=rvcm.second;

      float3 rcm=(p[i].m*ri+p[j].m*rj)/(p[i].m+p[j].m); 
      float3 vcm=(p[i].m*vi+p[j].m*vj)/(p[i].m+p[j].m);

      //float mvi=sqrt(vi.x*vi.x+vi.y*vi.y+vi.z*vi.z);
      //float mvj=sqrt(vj.x*vj.x+vj.y*vj.y+vj.z*vj.z);
      //float mvcm=sqrt(vcm.x*vcm.x+vcm.y*vcm.y+vcm.z*vcm.z);    
      ///float6 rijcm=rLS2CM(ri, rj, rcm);
      ///float3 ricm=rijcm.first;
      ///float3 rjcm=rijcm.second;

      ///float6 vijcm=pLS2CM(vi, vj, vcm);
      ///float3 vicm=vijcm.first;
      ///float3 vjcm=vijcm.second;
      
      // dr
      float3 r1vec=ri-rj;
      // dV
      float3 v1vec=vi-vj;

      float x1=r1vec.x;
      float y1=r1vec.y;
      float z1=r1vec.z;
      float vx1=v1vec.x;
      float vy1=v1vec.y;
      float vz1=v1vec.z;
      //float rdist1sqr=x1*x1+y1*y1+z1*z1;
      //float PK1=alpha/rdist1;
      //float rdist3=rdist1*rdist1sqr;
      
      mi=p[i].m;
      float mj=p[j].m;
      float mt=mi+mj;                               // total mass = m1 + m2
      float mr=mi*mj/mt;
      //field_i_j_interaction - simple calculation.
      float kij=mr*((vi.x-vj.x)*(vi.x-vj.x)+(vi.y-vj.y)*(vi.y-vj.y)+(vi.z-vj.z)*(vi.z-vj.z))/2;
      float rdist1=sqrtf(x1*x1+y1*y1+z1*z1);
      float KUR=alpha/kij/rdist1;
      //here we define, we will use the results of calculation_i_j_interaction(...)
      //(FLAG=true) or field_i_j_interaction(...) (FLAG=false).
      //if(KUR<EPS_INTERACTION) FLAG=false;
      //float F=alpha/rdist1sqr;
      //float Fx=F*x1/rdist1;
      //float Fy=F*y1/rdist1;
      //float Fz=F*z1/rdist1;
      //std::cout<<" "<<rdist3<<std::endl;
      float ad=alpha*deltatime;
      float adr=ad/mr;
      float adr1=adr/rdist1;
      float adr2=adr1/rdist1;
      std::cout<<"r1="<<rdist1<<",a2="<<adr2<<std::endl;
      float dtmr=adr2/rdist1;
      float dVx=x1*dtmr;
      float dVy=y1*dtmr;
      float dVz=z1*dtmr;
      float vx2=vx1+dVx;
      float vy2=vy1+dVy;
      float vz2=vz1+dVz;
      float x2=x1+(vx1+dVx/2)*deltatime;
      float y2=y1+(vy1+dVy/2)*deltatime;
      float z2=z1+(vz1+dVz/2)*deltatime;
      float rdist2sqr=x2*x2+y2*y2+z2*z2;
      float rdist2=sqrt(rdist2sqr);
      //energy conservation law: K1+P1=K2+P2
      float KK2=mr*(vx2*vx2+vy2*vy2+vz2*vz2)/2;
      float KK1=mr*(vx1*vx1+vy1*vy1+vz1*vz1)/2;
      //correct r
      float rcor=alpha/(KK1+alpha/rdist1-KK2);
      //float rcor=rdist2;
      //correction
      float cor=rcor/rdist2;
      x2*=cor;
      y2*=cor;
      z2*=cor;
      //coordinates of i and j particles
      float mit=-mi/mt;
      float mjt=mj/mt;
      float3 frj={mit*x2, mit*y2, mit*z2};
      float3 fri={mjt*x2, mjt*y2, mjt*z2};
      //velocities of i and j particles
      float3 fvj={mit*vx2, mit*vy2, mit*vz2};
      float3 fvi={mjt*vx2, mjt*vy2, mjt*vz2};
      
#ifdef DEB_NOPAR
      Ekrli[i][j]=md/4*(v1vec.x*v1vec.x+v1vec.y*v1vec.y+v1vec.z*v1vec.z);     
      r1iji[i][j]=sqrt(r1vec.x*r1vec.x+r1vec.y*r1vec.y+r1vec.z*r1vec.z);
#endif
      float6 resvecr1v1=rotation(r1vec, v1vec, R, Ri, i, j);
      float3 r1resvec=resvecr1v1.first;
      float3 v1resvec=resvecr1v1.second;
      ////std::cout<<"345: "<<v1resvec.x<<" "<<v1resvec.y<<" "<<v1resvec.z<<std::endl;
      ////Ekrci[i][j]=md/4*(v1resvec.x*v1resvec.x+v1resvec.y*v1resvec.y+v1resvec.z*v1resvec.z);

      //float vi2=vi.x*vi.x+vi.y*vi.y+vi.z*vi.z;
      //float vj2=vj.x*vj.x+vj.y*vj.y+vj.z*vj.z;
      //LS->CM: Ep=const=alpha/rij, dEk=-1/2*(mi+mj)*Vcm^2.
      ////float dEk=md*(vcm.x*vcm.x+vcm.y*vcm.y+vcm.z*vcm.z);
      //should be Eiji[i][j]=ER1+dEk
      ////float dEijiC=Eiji[i][j]-ER1-dEk;
      //std::cout<<"dEijiC="<<dEijiC<<std::endl;
#ifdef DEB_NOPAR
      ELSCM=Eiji[i][j];
      //std::cout<<" Eiji[i][j]="<<Eiji[i][j]<<" ER1="<<ER1<<" dEijiC="<<dEijiC<<" dEk="<<dEk<<std::endl;
      ////if(fabs(dEijiC)>MAX_LSCM)
      ////{
      ////  MAX_LSCM=fabs(dEijiC);
      ////  I_LSCM=i;
      ////  J_LSCM=j;
      ////}

      ////if(fabs(dEijiC)>MAX_C1)
      ////{
      ////  MAX_C1=fabs(dEijiC);
      ////  I_C1=i;
      ////  J_C1=j;
      ////}
#endif
      
/*
      std::cout<<"CHECK:"<<std::endl;
      std::cout<<"WAS: (1,2,3) (0,0,1)"<<std::endl;
      float3 ch1={1, 2, 3};
      float3 ch2={0, 0, 1};
      float6 ch12i=rotation(ch1, ch2, R, Ri, Eu);
      fourvec chch;
      chch.rj=ch12i.first;
      chch.vj=ch12i.second;
      chch.ri=ch12i.first;
      chch.vi=ch12i.second;
      
      //std::cout<<" ||1="<<chch.rj.x*chch.rj.x+chch.rj.y*chch.rj.y+chch.rj.z*chch.rj.z<<std::endl;
      //std::cout<<" ||2="<<chch.vj.x*chch.vj.x+chch.vj.y*chch.vj.y+chch.vj.z*chch.vj.z<<std::endl;
      
      fourvec ch12f=rotation_back(chch, R);
      float3 t1=ch12f.rj;
      float3 t2=ch12f.vj;
      float3 t3=ch12f.ri;
      float3 t4=ch12f.vi;
      std::cout<<" t1: "<<t1.x<<" "<<t1.y<<" "<<t1.z<<std::endl;
      std::cout<<" t2: "<<t2.x<<" "<<t2.y<<" "<<t2.z<<std::endl;
      std::cout<<" t3: "<<t3.x<<" "<<t3.y<<" "<<t3.z<<std::endl;
      std::cout<<" t4: "<<t4.x<<" "<<t4.y<<" "<<t4.z<<std::endl;
*/

      //std::cout<<"x="<<r1resvec.x<<",y="<<r1resvec.y<<std::endl;
      //float rdistmr=sqrt(r1resvec.x*r1resvec.x+r1resvec.y*r1resvec.y);//+r1resvec.z*r1resvec.z);
      //std::cout<<"rdistmr="<<rdistmr<<std::endl;
      //if(rdistmr<1.e-10) std::cout<<"***ERROR rdistmr=0"<<std::endl;
      //std::cout<<"rdistmr="<<rdistmr<<std::endl;
      //if(rdistmr<0.001) rdistmr=0.001;
      //float pij=alpha/rdistmr;
      //float kij=p[i].m*p[j].m/(p[i].m+p[j].m)*(rdistmr*rdistmr)/2;//+v1resvec.z*v1resvec.z)/2;
      //float UKR=pij/kij;
      fourvec resrjvjrivi=calculation_i_j_interaction(r1resvec, v1resvec, i, j);
      //fourvec resrjvjrivi;
      //if(UKR<EPS_INTERACTION) resrjvjrivi=field_i_j_interaction(r1resvec, v1resvec, i, j);
      //else                    resrjvjrivi=calculation_i_j_interaction(r1resvec, v1resvec, i, j);     
      
#ifdef DEB_NOPAR
      float3 resv0=resrjvjrivi.vj;
      float3 resv1=resrjvjrivi.vi;
      float Ekcmf=p[i].m*(resv1.x*resv1.x+resv1.y*resv1.y+resv1.z*resv1.z)/2+p[j].m*(resv0.x*resv0.x+resv0.y*resv0.y+resv0.z*resv0.z)/2;
      float dEkcmf=Ekcmf-Ekrci[i][j]; // ................print    
      Ekrci[i][j]-md/4*(v1resvec.x*v1resvec.x+v1resvec.y*v1resvec.y+v1resvec.z*v1resvec.z);
#endif
      fourvec bbrjvjrivi=rotation_back(resrjvjrivi, R, i, j);

      float3 brj=bbrjvjrivi.rj;
      float3 bvj=bbrjvjrivi.vj;
      float3 bri=bbrjvjrivi.ri;
      float3 bvi=bbrjvjrivi.vi;
      
#ifdef DEB_NOPAR
      float Ekrlf=p[i].m*(bvi.x*bvi.x+bvi.y*bvi.y+bvi.z*bvi.z)/2+p[j].m*(bvj.x*bvj.x+bvj.y*bvj.y+bvj.z*bvj.z)/2;
      float dEkrlf=Ekrlf-Ekrli[i][j]; // ................print    
      float r2ls_DIF2=sqrt((bri.x-brj.x)*(bri.x-brj.x)+(bri.y-brj.y)*(bri.y-brj.y)+(bri.z-brj.z)*(bri.z-brj.z));
      float dr2_DIF2=r2ls_DIF2-r2cmij[i][j];
      if(fabs(dr2_DIF2)>MAX_DIF2)
      {
        MAX_DIF2=fabs(dr2_DIF2);
        ERR_DIF2=fabs(dr2_DIF2/r2ls_DIF2);
        I_DIF2=i;
        J_DIF2=j;
      }      
      float dEf=ERB2-Eiji[i][j]+md*(vcm.x*vcm.x+vcm.y*vcm.y+vcm.z*vcm.z);
      //std::cout<<"dEf="<<dEf<<std::endl;
      if(fabs(dEf)>MAX_CMLS)
      {
        MAX_CMLS=fabs(dEf);
        I_CMLS=i;
        J_CMLS=j; 
      }
#endif

      // bvi +bvj=0 check
/*
//////////////////      
// !!!NO OUTPUT!!!
      if((i==263 && j==805) || (j==263 && i==805))
      {
        float fl=bvj.x+bvi.x;
        //if(fabs(fl)>1.e-13)
        {
          std::cout<<"*l* fl="<<fl<<" bvix="<<bvi.x<<" bvjx="<<bvj.x<<" i="<<i<<" j="<<j<<std::endl;
          //sleep(1);
        }
        //if(fabs(bvj.y+bvi.y)>1.e-13)
        //{
        //  std::cout<<"*l* bviy="<<bvi.y<<" bvjy="<<bvj.y<<" i="<<i<<" j="<<std::endl;
        //  sleep(1);
        //}
        //float fl=fabs(bvj.z+bvi.z);
        ////////if(fl>1.e-13)
        //{
        //  std::cout<<"*l* fl="<<fl<<" bviz="<<bvi.z<<" bvjz="<<bvj.z<<" i="<<i<<" j="<<j<<std::endl;
        //  sleep(1);
        //}
      }
// !!!NO OUTPUT!!!
//////////////////
*/

/*      
//////////////////
//HERE
      if((i==263 && j==805) || (j==263 && i==805))
      {
        float dcx=2*vcm.x-(vi.x+vj.x);
        //float dcx=2*vcm.x-vi.x-vj.x;
        if(i==263 && j==805)
        {
          //std::cout.precision(15);
          VCM=2*vcm.x;
          UCM=vi.x+vj.x;
          //std::cout<<"vi: "<<vi.x<<" "<<vi.y<<" "<<vi.z<<std::endl;
          //std::cout<<"vj: "<<vj.x<<" "<<vj.y<<" "<<vj.z<<std::endl;
        }
        //if(fabs(2*vcm.x-vi.x-vj.x)>1.e-10)
        {
          //std::cout.precision(15);
          std::cout<<"*v* d="<<dcx<<" vcm.x="<<vcm.x<<" vij="<<(vi.x+vj.x)/2<<" i="<<i<<" j="<<j;
          if(j==263 && i==805)
          {
            std::cout<<" DCM="<<2*vcm.x-VCM<<" dCM="<<vi.x+vj.x-UCM;
            //std::cout<<"vi: "<<vi.x<<" "<<vi.y<<" "<<vi.z<<std::endl;
            //std::cout<<"vj: "<<vj.x<<" "<<vj.y<<" "<<vj.z<<std::endl;
          }
          std::cout<<std::endl;
          //sleep(1);
        }
        //if(fabs(2*vcm.y-vi.y-vj.y)>1.e-10)
        //{
        //  std::cout<<"*v* d="<<2*vcm.y-vi.y-vj.y<<" vcm.y="<<vcm.y<<" vij="<<(vi.y+vj.y)/2<<" i="<<i<<" j="<<j<<std::endl;
        //  //sleep(1);
        //}
        ////////if(fabs(2*vcm.z-vi.z-vj.z)>1.e-10)
        //{
        //  std::cout<<"*v* d="<<2*vcm.z-vi.z-vj.z<<" vcm.z="<<vcm.z<<" vij="<<(vi.z+vj.z)/2<<" i="<<i<<" j="<<j<<std::endl;
        //  //sleep(1);
        //}
        //sleep(1);
        
      }
//////////////////
////////i=0 j=2///
//////////////////
*/
      

      ///float6 bbrjri=rCM2LS(bri, brj, rcm, vcm);

      ///float6 bbvjvi=pCM2LS(bvi, bvj, vcm);

      ///float3 bbri=bbrjri.first;
      ///float3 bbrj=bbrjri.second;
      ///float3 bbvi=bbvjvi.first;
      ///float3 bbvj=bbvjvi.second;

      float3 b2ri=bri+rcm+vcm*deltatime;
      float3 b2rj=brj+rcm+vcm*deltatime;
      float3 b2vi=bvi+vcm;
      float3 b2vj=bvj+vcm;
#ifdef DEB_NOPAR
      float Eki=p[i].m*(b2vi.x*b2vi.x+b2vi.y*b2vi.y+b2vi.z*b2vi.z)/2;
      float Ekj=p[j].m*(b2vj.x*b2vj.x+b2vj.y*b2vj.y+b2vj.z*b2vj.z)/2;
      dEk21ls[i][j]=Eki+Ekj-Eklsi[i][j];
      float dEklscm=dEk21ls[i][j]-dEk21cm[i][j];     
      //std::cout<<"dEklscm="<<dEklscm<<std::endl;
      /*if(fabs(dEklscm)>MAX_DIF3)
      {
        MAX_DIF3=fabs(dEklscm);
        if(fabs(dEk21ls[i][j])<1.e-5) std::cout<<"dEk21ls[i][j]="<<dEk21ls[i][j]<<std::endl;
        ERR_DIF3=fabs(dEklscm/dEk21ls[i][j]);
        I_DIF3=i;
        J_DIF3=j;
      }*/
      //float 
      float dra1=sqrt((b2ri.x-b2rj.x)*(b2ri.x-b2rj.x)+(b2ri.y-b2rj.y)*(b2ri.y-b2rj.y)+(b2ri.z-b2rj.z)*(b2ri.z-b2rj.z));
      float Ea1=p[i].m*(b2vi.x*b2vi.x+b2vi.y*b2vi.y+b2vi.z*b2vi.z)/2+p[j].m*(b2vj.x*b2vj.x+b2vj.y*b2vj.y+b2vj.z*b2vj.z)/2+alpha/dra1;
      float dEa1=Ea1-Eiji[i][j];
      //std::cout<<"dEa1="<<dEa1<<std::endl;
      if(fabs(dEa1)>MAX_DA1)
      {
        MAX_DA1=fabs(dEa1);
        I_DA1=i;
        J_DA1=j;
      }
#endif
      // j -> 0
      // i -> 1
      float3 bbrj, bbvj, bbri, bbvi;
      if(KUR<EPS_INTERACTION) bbrj=b2rj;      
      else bbrj=frj;
      if(KUR<EPS_INTERACTION) bbvj=b2vj;      
      else bbvj=fvj;       
      if(KUR<EPS_INTERACTION)bbri=b2ri;       
      else bbri=fri;      
      if(KUR<EPS_INTERACTION) bbvi=b2vi;
      else bbvi=fvi;
      
      if(KUR<EPS_INTERACTION) std::cout<<"***KUR="<<KUR<<std::endl;
      /*bbrj=b2rj;
      bbvj=b2vj;
      bbri=b2ri;
      bbvi=b2vi;*/
      
      
      float3 dRj=bbrj-rj-vj*deltatime;
      float3 dVj=bbvj-vj;
      float3 dRi=bbri-ri-vi*deltatime;
      float3 dVi=bbvi-vi;
      
      
      //float correction dPij=-dPji
      dVi.x=(dVi.x-dVj.x)/2;
      dVi.y=(dVi.y-dVj.y)/2;
      dVi.z=(dVi.z-dVj.z)/2;     
      dVj.x=-dVi.x;
      dVj.y=-dVi.y;
      dVj.z=-dVi.z;
#ifdef DEB_NOPAR
//////FINAL CHECK/////////////////
      float3 nri=ri+vi*deltatime+dRi;
      float3 nrj=rj+vj*deltatime+dRj;
      float3 nvi=vi+dVi;
      float3 nvj=vj+dVj;
      float dra3=sqrt((nri.x-nrj.x)*(nri.x-nrj.x)+(nri.y-nrj.y)*(nri.y-nrj.y)+(nri.z-nrj.z)*(nri.z-nrj.z));
      float Ea3=p[i].m*(nvi.x*nvi.x+nvi.y*nvi.y+nvi.z*nvi.z)/2+p[j].m*(nvj.x*nvj.x+nvj.y*nvj.y+nvj.z*nvj.z)/2+alpha/dra3;
      float dEa3=Ea3-Eiji[i][j];
      //std::cout<<"dEa3="<<dEa3<<std::endl;
      if(fabs(dEa3)>MAX_DA3)
      {
        MAX_DA3=fabs(dEa3);
        I_DA3=i;
        J_DA3=j;
      }
//////FINAL CHECK/////////////////    
      //End of float correction dPij=-dPji
      float fs=dVj.x+dVi.x;
      //if(fabs(dVj.x+dVi.x)>1.e-10) std::cout<<"*s* dVix="<<dVi.x<<" dVjx="<<dVj.x<<" i="<<i<<" j="<<std::endl;
      //if(fabs(dVj.y+dVi.y)>1.e-10) std::cout<<"*s* dViy="<<dVi.y<<" dVjy="<<dVj.y<<" i="<<i<<" j="<<std::endl;
      float fsx=dVj.x+dVi.x;
      float fsy=dVj.y+dVi.y;
      float fsz=dVj.z+dVi.z;
      /*
      //if(fabs(fs)>1.e-10)
      {
        if(fabs(fs)>MAX_S)
        {
          MAX_S=fs;
          I_S=i;
          J_S=j;
        }
        if((i==263 && j==805) || (j==263 && i==805))
        {
          std::cout<<"*s* fsx="<<fsx<<" fsy="<<fsy<<" fsz="<<fsz<<" dViz="<<dVi.z<<" dVjz="<<dVj.z<<" i="<<i<<" j="<<j<<std::endl;
        }
      }
      */
#endif       
      dRxi+=dRi.x;
      dRyi+=dRi.y;
      dRzi+=dRi.z;
      dVxi+=dVi.x;
      dVyi+=dVi.y;
      dVzi+=dVi.z;

#ifdef DEB_NOPAR
      dPijx[i][j]=dVxi;
      dPijy[i][j]=dVyi;
      dPijz[i][j]=dVzi;
#endif
      //std::cout<<"i="<<i<<" j="<<j<<" dVxi="<<dVxi<<" dVyi="<<dVyi
      //         <<" dVzi="<<dVzi<<std::endl;     
    //}//End of while
    }//End of j
    dRxx[i]=dRxi;
    dRyy[i]=dRyi;
    dRzz[i]=dRzi;    
    dVx[i]=dVxi;
    dVy[i]=dVyi;
    dVz[i]=dVzi;
#ifdef DEB_RV
    R1[i]=dRxx[i];
    R2[i]=dRyy[i];
    R3[i]=dRzz[i];
    VR1[i]=dVx[i];
    VR2[i]=dVy[i];
    VR3[i]=dVz[i];
#endif
  }//End of i
/*
  // check:
  for(int i=0; i<numBodies; ++i)
  {
    for(int j=0; j<numBodies; ++j) if(i!=j)
    {
      if((i==263 && j==805) || (i==805 && j==263))
      {
        std::cout<<"C i="<<i<<" j="<<j<<" dPijx[i][j]+dPijx[j][i]="<<dPijx[i][j]+dPijx[j][i]<<std::endl;
        std::cout<<"C i="<<i<<" j="<<j<<" dPijy[i][j]+dPijy[j][i]="<<dPijy[i][j]+dPijy[j][i]<<std::endl;
        std::cout<<"C i="<<i<<" j="<<j<<" dPijz[i][j]+dPijz[j][i]="<<dPijz[i][j]+dPijz[j][i]<<std::endl;
      }
    }
  }
*/
/*
  for(int i=1; i<numBodies; ++i)
  {
    for(int j=i+1; j<numBodies; ++j)
    {
      //if(i==137 && j==645)
      {
        //if(fabs(dPijx[i][j]+dPijx[j][i])>1.e-12)
        {
          float dp=dPijx[i][j]+dPijx[j][i];
          if(fabs(dp)>MAX_D)
          {
            MAX_D=fabs(dp);
            I_D=i;
            J_D=j;
          }
          //////if(i==263 && j==805) std::cout<<"*x* i="<<i<<" j="<<j<<" "<<dPijx[i][j]<<" "<<dPijx[j][i]<<" dx="<<dp<<std::endl;
        }
        //////if(fabs(dPijy[i][j]+dPijy[j][i])>1.e-12)
        //{
        //  std::cout<<"*y* i="<<i<<" j="<<j<<" "<<dPijy[i][j]<<" "<<dPijy[j][i]<<" dy="<<dPijy[i][j]+dPijy[j][i]<<std::endl;
        //  //exit(0);
        //}
        ////if(fabs(dPijz[i][j]+dPijz[j][i])>1.e-12)
        //{
        //  std::cout<<"*z* i="<<i<<" j="<<j<<" "<<dPijz[i][j]<<" "<<dPijz[j][i]<<" dz="<<dPijz[i][j]+dPijz[j][i]<<std::endl;
        //  //exit(0);
        //}
      }
    }
  }
*/
#ifdef DEB_NOPAR
  for(int i=0; i<numBodies; ++i)
  {
    for(int j=0; j<numBodies; ++j) if(i!=j)
    {
      float rda2=sqrt((p[i].x-p[j].x)*(p[i].x-p[j].x)+(p[i].y-p[j].y)*(p[i].y-p[j].y)+(p[i].z-p[j].z)*(p[i].z-p[j].z));
      float Ea2=p[i].m*(p[i].vx*p[i].vx+p[i].vy*p[i].vy+p[i].vz*p[i].vz)/2+p[j].m*(p[j].vx*p[j].vx+p[j].vy*p[j].vy+p[j].vz*p[j].vz)/2+alpha/rda2;
      float dEa2=Ea2-Eiji[i][j];
      //std::cout<<"dEa2="<<dEa2<<std::endl;
      if(fabs(dEa2)>MAX_DA2)
      {
        MAX_DA2=fabs(dEa2);
        I_DA2=i;
        J_DA2=j;
      }
    }
  }
  /*for(int i=0; i<numBodies; ++i)
  {
    for(int j=0; j<numBodies; ++j) if(i!=j)
    {
      float dra4=sqrt((p[i].x-p[j].x)*(p[i].x-p[j].x)+(p[i].y-p[j].y)*(p[i].y-p[j].y)+(p[i].z-p[j].z)*(p[i].z-p[j].z));
      float Ea4=p[i].m*(p[i].vx*p[i].vx+p[i].vy*p[i].vy+p[i].vz*p[i].vz)/2+p[j].m*(p[j].vx*p[j].vx+p[j].vy*p[j].vy+p[j].vz*p[j].vz)/2+alpha/dra4;
      float dEa4=Ea4-Eiji[i][j];
      std::cout<<"dEa4="<<dEa4<<std::endl;
    }
  }*/
#endif
    
#ifdef ACCEL
#pragma acc parallel loop present(p,dRxx,dRyy,dRzz,dVx,dVy,dVz)
#else
#pragma omp parallel for //reduction(max:MAX_vx,MAX_vy,MAX_vz,MAX_dvx,MAX_dvy,MAX_dvz)                     //T=1744 micros
  //#pragma omp parallel for schedule(dynamic)  //T=3490 micros
  //#pragma simd                                //T=6064 microseconds
  //#pragma distribute_point                    //T=3519 micros
#endif
  for(int i=0; i<numBodies; ++i)
  {
    p[i].x+=p[i].vx*deltatime+dRxx[i];
    p[i].y+=p[i].vy*deltatime+dRyy[i];
    p[i].z+=p[i].vz*deltatime+dRzz[i];
    
    p[i].vx+=dVx[i];
    p[i].vy+=dVy[i];
    p[i].vz+=dVz[i];
    //std::cout<<"V: "<<p[i].vx<<" "<<p[i].vy<<" "<<p[i].vz<<std::endl;
    //std::cout<<"dV: "<<dVx[i]<<" "<<dVy[i]<<" "<<dVz[i]<<std::endl;
    //std::cout<<"K="<<K<<" md="<<md<<std::endl;
    //float sum=p[i].vx*p[i].vx+p[i].vy*p[i].vy+p[i].vz*p[i].vz;
    //std::cout<<"sum="<<p[i].vx*p[i].vx+p[i].vy*p[i].vy+p[i].vz*p[i].vz<<std::endl;
    //if(sum>1.)
    //{
    //  std::cout<<"sum="<<sum<<std::endl;
    //  exit(0);
    //}
#ifdef MAKEVRML    
    vrmlpos[i][STEP].x=p[i].x;
    vrmlpos[i][STEP].y=p[i].y;
    vrmlpos[i][STEP].z=p[i].z;
#endif
#ifdef DEB_VNEAR
    if(fabs(p[i].vx)>MAX_vx) MAX_vx=fabs(p[i].vx);
    if(fabs(p[i].vy)>MAX_vy) MAX_vy=fabs(p[i].vy);
    if(fabs(p[i].vz)>MAX_vz) MAX_vz=fabs(p[i].vz);
    if(fabs(dVx[i])>MAX_dvx) MAX_dvx=fabs(dVx[i]);
    if(fabs(dVy[i])>MAX_dvy) MAX_dvy=fabs(dVy[i]);
    if(fabs(dVz[i])>MAX_dvz) MAX_dvz=fabs(dVz[i]);
#endif
  }
  float K=0.0;
  float P=0.0;
  P_init={0.0, 0.0, 0.0};
  double P_x=0.0;
  double P_y=0.0;
  double P_z=0.0;
  float p_en=0.0;
#ifdef ACCEL
#pragma acc parallel loop copy(p_en,K,P_x,P_y,P_z)
#else
#pragma omp parallel for reduction(+:p_en,K,P_x,P_y,P_z) //reduction(min:MAX_NEAR)
  //#pragma omp parallel for schedule(dynamic)  
  //#pragma simd                                
#endif  
  for(int i=0; i<numBodies; ++i)
  {
    //E_p[i]=0.0;
    P_x+=p[i].vx*p[i].m;
    P_y+=p[i].vy*p[i].m;
    P_z+=p[i].vz*p[i].m;
    float PP=0.0;
#ifdef ACCEL
#pragma acc loop vector reduction(+:PP)
#else
    //#pragma omp simd reduction(+:PP)
    //#pragma simd reduction(+:PP) //reduction(min:MAX_NEAR,MIN_K,MIN_KUR) reduction(max:MAX_K,MAX_KUR)
    //#pragma simd // works if to replace PP with p_en
#endif  
    for(int j=0; j<numBodies; ++j) if(i!=j)
    {
      float rdist1=sqrtf((p[i].x-p[j].x)*(p[i].x-p[j].x)+(p[i].y-p[j].y)*(p[i].y-p[j].y)+(p[i].z-p[j].z)*(p[i].z-p[j].z));
      //if(rdist1==0.0) std::cout<<"***ERROR: rdist1=0"<<std::endl;
      PP+=alpha/rdist1;
#ifdef DEB_KREL
      float Kij=md/4*((p[i].vx-p[j].vx)*(p[i].vx-p[j].vx)+(p[i].vy-p[j].vy)*(p[i].vy-p[j].vy)+(p[i].vz-p[j].vz)*(p[i].vz-p[j].vz));
      if(Kij>MAX_K)
      {
        MAX_K=Kij;
        I_MAX_K=i;
        J_MAX_K=j;
      }
      if(Kij<MIN_K)
      {
        MIN_K=Kij;
        I_MIN_K=i;
        J_MIN_K=j;
      }
      float Uij=alpha/rdist1;
      float KUR=Uij/Kij;
      if(KUR>MAX_KUR)
      {
        MAX_KUR=KUR;
        I_MAX_KUR=i;
        J_MAX_KUR=j;
      }
      if(KUR<MIN_KUR)
      {
        MIN_KUR=KUR;
        I_MIN_KUR=i;
        J_MIN_KUR=j;
      }
#endif

#ifdef DEB_NOPAR
      float deltaU=alpha/rdist1-U[i][j];
      //if(i==137 && j==645) std::cout<<"dU["<<i<<","<<j<<"]="<<deltaU<<" Uf="<<alpha/rdist<<" Ui="<<U[i][j]<<std::endl;
      if(fabs(deltaU)>MAX_U)
      {
        MAX_U=fabs(deltaU);
        I_U=i;
        J_U=j;
      }
      //U[i][j]=alpha/rdist;
      float Efin1=p[i].m*(p[i].vx*p[i].vx+p[i].vy*p[i].vy+p[i].vz*p[i].vz)/2+p[j].m*(p[j].vx*p[j].vx+p[j].vy*p[j].vy+p[j].vz*p[j].vz)/2+alpha/rdist1;
      float deltaE1=Efin1-Eiji[i][j];
      //std::cout<<"dE1="<<deltaE1<<std::endl;      
      float Ekif=p[i].m*(p[i].vx*p[i].vx+p[i].vy*p[i].vy+p[i].vz*p[i].vz)/2;
      float Ekjf=p[j].m*(p[j].vx*p[j].vx+p[j].vy*p[j].vy+p[j].vz*p[j].vz)/2;
      float distijf=sqrtf((p[i].x-p[j].x)*(p[i].x-p[j].x)+(p[i].y-p[j].y)*(p[i].y-p[j].y)+(p[i].z-p[j].z)*(p[i].z-p[j].z));
      float Epijf=alpha/distijf;
      float Eijf=Ekif+Ekjf+Epijf;
      float dFE=Eijf-Eiji[i][j];
      //std::cout<<"Fin: i="<<i<<" j="<<j<<" dFE="<<dFE<<" rrat="<<dFE/(Eijf+Eiji[i][j])<<std::endl;     
      if(fabs(dFE)>MAX_R)
      {
        MAX_R=fabs(dFE);
        I_R=i;
        J_R=j;
      }
#endif
#ifdef DEB_VNEAR
      if(MAX_NEAR>rdist1)
      {
        MAX_NEAR=rdist1;
        I_NEAR=i;
        J_NEAR=j;
      }
#endif
    }
    p_en+=PP;
    float v2=p[i].vx*p[i].vx+p[i].vy*p[i].vy+p[i].vz*p[i].vz;
    //Check |v|>1
    /*if(v2>=1.)
    {
      //std::cout<<"***ERROR: v2["<<i<<"]="<<v2<<" > 1"<<std::endl;
      //printf("***ERROR: v2[%d]=%f > 1",i,v2);
      exit(0);
    }*/
    //***ERROR: v2[945]=1.00552 > 1
    K+=p[i].m*v2/2;
#ifdef DEB_NOPAR
    //float dEU=E_p[i]-PING[i];
    //if(fabs(dEU)>MAX_P)
    //{
    //  MAX_P=fabs(dEU);
    //  I_P=i;
    //}
#endif
    //p_en+=E_p[i];
  }
  P_init.x=P_x;
  P_init.y=P_y;
  P_init.z=P_z;
  P=p_en;
/*
  float DV=0.0;
  int DI;
  for(int i=0; i<numBodies; ++i)
  {
    
    if(fabs(dVx[i])>DV)
    {
      DV=fabs(dVx[i]);
      DI=i;
    }
  }
  std::cout<<"DV="<<DV<<" DI="<<DI<<std::endl;
  
  std::cout<<"dVx:"<<std::endl;
  for(int i=0; i<numBodies; ++i) std::cout<<dVx[i]<<" ";
  std::cout<<std::endl;
*/

/*
///////////////////
//CHECK:///////////
  float dVxs=0.0;
  float dVys=0.0;
  float dVzs=0.0;
  for(int i=0; i<numBodies; ++i)
  {
    dVxs+=dVx[i];
    dVys+=dVy[i];
    dVzs+=dVz[i];
  }
  std::cout<<"1. Sum dVx[i], dVy[i], dVz[i] : dVxs="<<dVxs<<" dVys="<<dVys<<" dVzs="<<dVzs<<std::endl;

  float SdVxij=0.0;
  float SdVyij=0.0;
  float SdVzij=0.0;
  for(int i=0; i<numBodies; ++i)
  {
    for(int j=0; j<numBodies; ++j) if(i!=j)
    {
      SdVxij+=dPijx[i][j];
      SdVyij+=dPijy[i][j];
      SdVzij+=dPijz[i][j];
    }
  }
  std::cout<<"2. Sum dPijx[i][j], dPijy[i][j], dPijz: SdVxij="<<SdVxij<<" SdVyij="<<SdVyij<<" SdVzij="<<SdVzij<<std::endl;

  float S1dVxij=0.0;
  float S1dVyij=0.0;
  float S1dVzij=0.0;
  for(int i=0; i<numBodies; ++i)
  {
    for(int j=0; j<numBodies; ++j) if(i!=j)
    {
      S1dVxij+=dPijx[i][j]+dPijx[j][i];
      S1dVyij+=dPijy[i][j]+dPijy[j][i];
      S1dVzij+=dPijz[i][j]+dPijz[j][i];
    }
  }
  std::cout<<"3. Sum dPijx[i][j]+dPijx[j][i], dPijy[i][j]+dPijy[j][i], dPijz[i][j]+dPijz[j][i]: S1dVxij="
           <<S1dVxij<<" S1dVyij="<<S1dVyij<<" S1dVzij="<<S1dVzij<<std::endl;

  std::cout<<"The problem is here:"<<std::endl;
  for(int i=0; i<numBodies; ++i)
  {
    for(int j=0; j<numBodies; ++j) if (i!=j)
    {
      if((dPijx[i][j]+dPijx[j][i])!=0)
        std::cout<<"i="<<i<<" j="<<j<<" dPijx[i][j]+dPijx[j][i]="<<dPijx[i][j]+dPijx[j][i]<<"!=0"<<std::endl;
    }
  }
  std::cout<<std::endl;

  for(int i=0; i<numBodies; ++i)
  {
    for(int j=0; j<numBodies; ++j)
    {
      float dxx=dPijx[i][j]+dPijx[j][i];
      if(dxx!=0.0)
      {
        std::cout<<"***ERROR: dPijx+dPjix="<<dxx<<" i="<<i<<" j="<<j<<"!=0"<<std::endl;
        //exit(0);
      }
      float dyy=dPijy[i][j]+dPijy[j][i];
      if(dyy!=0.0)
      {
        std::cout<<"***ERROR: dPijy+dPjiy="<<dyy<<" i="<<i<<" j="<<j<<"!=0"<<std::endl;
        //exit(0);
      }
      float dzz=dPijz[i][j]+dPijz[j][i];
      if(dzz!=0.0)
      {
        std::cout<<"***ERROR: dPijz+dPjiz="<<dzz<<" i="<<i<<" j="<<j<<"!=0"<<std::endl;
        //exit(0);
      }
    }
  }
 
  std::cout<<"CHECK diag elements:"<<std::endl;
  for(int i=0; i<numBodies; ++i)
  {
    if(dPijx[i][i]!=0.0) std::cout<<"***ERROR: dPijx["<<i<<"]["<<i<<"]="<<dPijx[i][i]<<"!=0"<<std::endl;
    if(dPijy[i][i]!=0.0) std::cout<<"***ERROR: dPijy["<<i<<"]["<<i<<"]="<<dPijy[i][i]<<"!=0"<<std::endl;
    if(dPijz[i][i]!=0.0) std::cout<<"***ERROR: dPijz["<<i<<"]["<<i<<"]="<<dPijz[i][i]<<"!=0"<<std::endl;
  }
  
  
  //dPijx[i][j]=dVxi;
  //dPijy[i][j]=dVyi;
  //dPijz[i][j]=dVzi;


//CHECK////////////
///////////////////
*/

  P/=2;
  float dK=K-Kf;
  float dP=P-Pf;
  float dE=dK+dP;
  float E=K+P;
  
  //float dE21=K+P-Kf-Pf;
  //float dE21=E_init-P;
  float r=(E_init-P)/K;
  float rc=sqrt(r);
  ///std::cout<<"1. Ekin CORRECTION:"<<std::endl;
  ///std::cout<<"Ecurfull="<<K+P<<",Effull="<<Kf+Pf<<",dE21="<<dE21<<",r="<<r<<",(r)^0.5="<<sqrt(r)<<std::endl;
  ///std::cout<<"K="<<K<<",P="<<P<<",Kf="<<Kf<<",Pf="<<Pf<<",dE21="<<dE21<<",K-dE21="<<K-dE21
  ///         <<",r="<<r<<",r*K="<<r*K<<",r*K+P-Kf-Pf="<<r*K+P-Kf-Pf<<std::endl;

#ifdef EKCOR
//Kinetic energy correction
#ifdef ACCEL
#pragma acc parallel loop present(p) copy(rc)
#else
  //#pragma omp parallel for schedule(dynamic)
#pragma simd 
#endif
  for(int i=0; i<numBodies; ++i)
  {   
    p[i].vx*=rc;
    p[i].vy*=rc;
    p[i].vz*=rc;
  }
//End of full energy correction
#endif

  K=0.0;
#ifdef ACCEL
#pragma acc parallel loop copy(K) present(p)
#else
//#pragma omp parallel for schedule(dynamic)
#pragma simd reduction(+:K)
#endif
  for(int i=0; i<numBodies; ++i)
  {
    K+=p[i].m*(p[i].vx*p[i].vx+p[i].vy*p[i].vy+p[i].vz*p[i].vz)/2;
  }

///dK=K-Kf;

#ifdef DEB_RV
#ifdef ACCEL
#pragma acc update host(R1,R2,R3,VR1,VR2,VR3)
#endif
  for(int i=0; i<numBodies; ++i)
  {
    std::cout<<"R1 VR1: "<<R1[i]<<" "<<R2[i]<<" "<<R3[i]<<" "<<VR1[i]<<" "<<VR2[i]<<" "<<VR3[i]<<std::endl;
  }
#endif

  ///std::cout<<"2. Ekin CORRECTION:"<<std::endl;
  ///std::cout<<"Ecurfull="<<K+P<<",Effull="<<Kf+Pf<<std::endl;

  //std::cout<<"MAX_S="<<MAX_S<<" I_S="<<I_S<<" J_S="<<J_S<<" dPx="<<md*dP.x<<" dPy="<<md*dP.y<<" dPz="<<md*dP.z<<std::endl;
  /*
  std::cout<<"MAX_R="<<MAX_R<<" I_R="<<I_R<<" J_R="<<J_R<<std::endl;
  std::cout<<"MAX_R1="<<MAX_R1<<" I_R1="<<I_R1<<" J_R1="<<J_R1<<std::endl;
  std::cout<<"MAX_C1="<<MAX_C1<<" I_C1="<<I_C1<<" J_C1="<<J_C1<<std::endl;
  std::cout<<"MAX_RR2="<<MAX_RR2<<" I_RR2="<<I_RR2<<" J_RR2="<<J_RR2<<std::endl;
  std::cout<<"MAX_RRB2="<<MAX_RRB2<<" I_RRB2="<<I_RRB2<<" J_RRB2="<<J_RRB2<<std::endl;
  std::cout<<"MAX_DERRB12="<<MAX_DERRB12<<" I_DERRB12="<<I_DERRB12<<" J_DERRB12="<<J_DERRB12<<std::endl;
  std::cout<<"MAX_R3="<<MAX_R3<<" I_R3="<<I_R3<<" J_R3="<<J_R3<<std::endl;
  std::cout<<"MAX_R4="<<MAX_R4<<" I_R4="<<I_R4<<" J_R4="<<J_R4<<std::endl;
  std::cout<<"MAX_LSCM="<<MAX_LSCM<<" I_LSCM="<<I_LSCM<<" J_LSCM="<<J_LSCM<<std::endl;
  std::cout<<"MAX_CMLS="<<MAX_CMLS<<" I_CMLS="<<I_CMLS<<" J_CMLS="<<J_CMLS<<std::endl;
  std::cout<<"MAX_DA1="<<MAX_DA1<<" I_DA1="<<I_DA1<<" J_DA1="<<J_DA1<<std::endl;
  std::cout<<"MAX_DA2="<<MAX_DA2<<" I_DA2="<<I_DA2<<" J_DA2="<<J_DA2<<std::endl;
  std::cout<<"MAX_DA3="<<MAX_DA3<<" I_DA3="<<I_DA3<<" J_DA3="<<J_DA3<<std::endl;
  */
  //std::cout<<"MAX_DIF1="<<MAX_DIF1<<" ERR_DIF1="<<ERR_DIF1<<" I_DIF1="<<I_DIF1<<" J_DIF1="<<J_DIF1<<std::endl;
  //std::cout<<"MAX_DIF2="<<MAX_DIF2<<" ERR_DIF2="<<ERR_DIF2<<" I_DIF2="<<I_DIF2<<" J_DIF2="<<J_DIF2<<std::endl;
  //std::cout<<"MAX_DIF3="<<MAX_DIF3<<" ERR_DIF3="<<ERR_DIF3<<" I_DIF3="<<I_DIF3<<" J_DIF3="<<J_DIF3<<std::endl;
  //std::cout<<"MAX_NEAR_PAIR="<<MAX_NEAR_PAIR<<" I_NEAR_PAIR="<<I_NEAR_PAIR<<" J_NEAR_PAIR="<<J_NEAR_PAIR<<std::endl;
  //std::cout<<"MAX_NEAR="<<MAX_NEAR<<" I_NEAR="<<I_NEAR<<" J_NEAR="<<J_NEAR<<std::endl;
  //std::cout<<"MAX_U="<<MAX_U<<" I_U="<<I_U<<" J_U="<<J_U<<std::endl;
  //std::cout<<"MAX_P="<<MAX_P<<" I_P="<<I_P<<" J_P="<<J_P<<std::endl;
  //std::cout<<"MAX_E="<<MAX_E<<" I_E="<<I_E<<" J_E="<<J_E<<" MIN_d="<<MIN_d<<std::endl;
  //std::cout<<"MAX_D="<<MAX_D<<" I_D="<<I_D<<" J_D="<<J_D<<std::endl;
  //|Vi|=0.01414
  //std::cout<<"MAX_K="<<MAX_K<<" I_MAX_K="<<I_MAX_K<<" J_MAX_K="<<J_MAX_K<<std::endl;
  //std::cout<<"MIN_K="<<MIN_K<<" I_MIN_K="<<I_MIN_K<<" J_MIN_K="<<J_MIN_K<<std::endl;
  //std::cout<<"MAX_KUR="<<MAX_KUR<<" I_MAX_KUR="<<I_MAX_KUR<<" J_MAX_KUR="<<J_MAX_KUR<<std::endl;
  //std::cout<<"MIN_KUR="<<MIN_KUR<<" I_MIN_KUR="<<I_MIN_KUR<<" J_MIN_KUR="<<J_MIN_KUR<<std::endl;
  //std::cout<<"MAX_vx="<<MAX_vx<<" MAX_vy="<<MAX_vy<<" MAX_vz="<<MAX_vz
  //         <<" MAX_dvx="<<MAX_dvx<<" MAX_dvy="<<MAX_dvy<<" MAX_dvz="<<MAX_dvz<<std::endl;
  std::cout<<"E="<<K+P<<" dE/E="<<dE/E<<" K="<<K<<" dK="<<dK<<" dk="<<dK/numBodies<<" P="<<P<<"\n dP="<<dP<<" dp="
           <<2*dP/(numBodies*(numBodies-1))<<" Px="<<P_init.x<<" Py="<<P_init.y<<" Pz="<<P_init.z<<std::endl;
  //sleep(1);
  Kf=K;
  Pf=P;
}//End of interaction()

#ifdef MAKEVRML
void makeVRML()
{
  static int vrmlinit=0;
  std::ofstream fout;
  if(vrmlinit==0)
  {
    if(!fout.is_open()) fout.open("vrmlpositions");
    //std::cout<<"OUT"<<std::endl;
    fout<<"#VRML V2.0 utf8"<<std::endl;
    //else{pause();}
    //fout<<"Group {\n"<<"children [\n"<<"Shape {\n"<<"appearance DEF White Appearance {\n"<<"material Material {}\n}\n"<<"geometry Sphere { radius 0.1}\n}, \n";
    
    //fout<<"Shape {\n appearance Appearance {\n material Material { emissiveColor 1 0 0 }\n}\n";
    //fout<<"geometry IndexedLineSet {\n coord Coordinate {\n point [ 10 0 10, -10 0 10, -10 0 -10, 10 0 -10, 0 20 0]\n}\n";
    //fout<<" coordIndex [ 0 1 2 3 0 -1, 0 4 -1, 1 4 -1, 2 4 -1, 3 4 -1 ]\n}\n}";
    
    // Ox axis red
    fout<<"Shape {\n appearance Appearance {\n material Material { emissiveColor 1 0 0 }\n}\n";
    fout<<"geometry IndexedLineSet {\n coord Coordinate {\n point [ -100 0 0, 100 0 0 ]\n}\n";
    fout<<" coordIndex [ 0 1 -1 ]\n}\n}";
    //Oy axis blue
    fout<<"Shape {\n appearance Appearance {\n material Material { emissiveColor 0 0 1 }\n}\n";
    fout<<"geometry IndexedLineSet {\n coord Coordinate {\n point [ 0 -100 0, 0 100 0 ]\n}\n";
    fout<<" coordIndex [ 0 1 -1 ]\n}\n}";
    //1 particle axis green
    fout<<"Shape {\n appearance Appearance {\n material Material { emissiveColor 0 1 0 }\n}\n";
    fout<<"geometry IndexedLineSet {\n coord Coordinate {\n point [ -100 0 1.0, 100 0 1.0 ]\n}\n";
    fout<<" coordIndex [ 0 1 -1 ]\n}\n}";
    
    fout<<"Group {\n"<<"children [\n";
    fout<<"DEF VIEW Viewpoint { position 0 0 50 orientation 0 -1 0 0 }\n";
    
    for(int i=0; i< numBodies; ++i)
    {
      if(i==0)
      {
        fout<<"DEF Atom"<<i<<" Transform {\n"<<"translation "<<vrmlpos[i][0].x<<" "<<vrmlpos[i][0].y<<" "<<vrmlpos[i][0].z<<"\n";
        fout<<"center "<<-vrmlpos[i][0].x<<" "<<-vrmlpos[i][0].y<<" "<<-vrmlpos[i][0].z<<"\n";
        //fout<<"center "<<0<<" "<<0<<" "<<0<<"\n";
        //fout<<"children Shape {\n appearance USE White\n geometry Sphere {radius 0.1} \n} \n";
        fout<<"children [ Shape {\n appearance Appearance {\n material Material { emissiveColor 1 0 0 ambientIntensity 1.0}\n}\n geometry Sphere { radius 1.0 } \n} \nDEF Atom"<<i<<"_sensor TouchSensor {}\n]";
        fout<<"},\n";
 			}
      else if(i==1)
      {
        fout<<"DEF Atom"<<i<<" Transform {\n"<<"translation "<<vrmlpos[i][0].x<<" "<<vrmlpos[i][0].y<<" "<<vrmlpos[i][0].z<<"\n";
        fout<<"center "<<-vrmlpos[i][0].x<<" "<<-vrmlpos[i][0].y<<" "<<-vrmlpos[i][0].z<<"\n";
        //fout<<"center "<<0<<" "<<0<<" "<<0<<"\n";
        //fout<<"children Shape {\n appearance USE White\n geometry Sphere {radius 0.1} \n} \n";
        fout<<"children [ Shape {\n appearance Appearance {\n material Material { emissiveColor 1 1 1 ambientIntensity 1.0}\n}\n geometry Sphere { radius 1.0 } \n} \nDEF Atom"<<i<<"_sensor TouchSensor {}\n]";
        fout<<"},\n";
      }
      else if(i==2)
      {
        fout<<"DEF Atom"<<i<<" Transform {\n"<<"translation "<<vrmlpos[i][0].x<<" "<<vrmlpos[i][0].y<<" "<<vrmlpos[i][0].z<<"\n";
        fout<<"center "<<-vrmlpos[i][0].x<<" "<<-vrmlpos[i][0].y<<" "<<-vrmlpos[i][0].z<<"\n";
        //fout<<"center "<<0<<" "<<0<<" "<<0<<"\n";
        //fout<<"children Shape {\n appearance USE White\n geometry Sphere {radius 0.1} \n} \n";
        fout<<"children [ Shape {\n appearance Appearance {\n material Material { emissiveColor 0 0 1 ambientIntensity 1.0}\n}\n geometry Sphere { radius 1.0 } \n} \nDEF Atom"<<i<<"_sensor TouchSensor {}\n]";
        fout<<"},\n";
      }
      else if(i==3)
      {
        fout<<"DEF Atom"<<i<<" Transform {\n"<<"translation "<<vrmlpos[i][0].x<<" "<<vrmlpos[i][0].y<<" "<<vrmlpos[i][0].z<<"\n";
        fout<<"center "<<-vrmlpos[i][0].x<<" "<<-vrmlpos[i][0].y<<" "<<-vrmlpos[i][0].z<<"\n";
        //fout<<"center "<<0<<" "<<0<<" "<<0<<"\n";
        //fout<<"children Shape {\n appearance USE White\n geometry Sphere {radius 0.1} \n} \n";
        fout<<"children [ Shape {\n appearance Appearance {\n material Material { emissiveColor 0 1 0 ambientIntensity 1.0}\n}\n geometry Sphere { radius 1.0 } \n} \nDEF Atom"<<i<<"_sensor TouchSensor {}\n]";
        fout<<"},\n";
      }
      else if(i==4)
      {
        fout<<"DEF Atom"<<i<<" Transform {\n"<<"translation "<<vrmlpos[i][0].x<<" "<<vrmlpos[i][0].y<<" "<<vrmlpos[i][0].z<<"\n";
        fout<<"center "<<-vrmlpos[i][0].x<<" "<<-vrmlpos[i][0].y<<" "<<-vrmlpos[i][0].z<<"\n";
        //fout<<"center "<<0<<" "<<0<<" "<<0<<"\n";
        //fout<<"children Shape {\n appearance USE White\n geometry Sphere {radius 0.1} \n} \n";
        fout<<"children [ Shape {\n appearance Appearance {\n material Material { emissiveColor 0 0.5 0 ambientIntensity 1.0}\n}\n geometry Sphere { radius 1.0 } \n} \nDEF Atom"<<i<<"_sensor TouchSensor {}\n]";
        fout<<"},\n";
      }
      else if(i==5)
      {
        fout<<"DEF Atom"<<i<<" Transform {\n"<<"translation "<<vrmlpos[i][0].x<<" "<<vrmlpos[i][0].y<<" "<<vrmlpos[i][0].z<<"\n";
        fout<<"center "<<-vrmlpos[i][0].x<<" "<<-vrmlpos[i][0].y<<" "<<-vrmlpos[i][0].z<<"\n";
        //fout<<"center "<<0<<" "<<0<<" "<<0<<"\n";
        //fout<<"children Shape {\n appearance USE White\n geometry Sphere {radius 0.1} \n} \n";
        fout<<"children [ Shape {\n appearance Appearance {\n material Material { emissiveColor 1 1 1 ambientIntensity 1.0}\n}\n geometry Sphere { radius 1.0 } \n} \nDEF Atom"<<i<<"_sensor TouchSensor {}\n]";
        fout<<"},\n";
      }
      else
      {
        fout<<"DEF Atom"<<i<<" Transform {\n"<<"translation "<<vrmlpos[i][0].x<<" "<<vrmlpos[i][0].y<<" "<<vrmlpos[i][0].z<<"\n";
        fout<<"center "<<-vrmlpos[i][0].x<<" "<<-vrmlpos[i][0].y<<" "<<-vrmlpos[i][0].z<<"\n";
        //fout<<"center "<<0<<" "<<0<<" "<<0<<"\n";
        //fout<<"children Shape {\n appearance USE White\n geometry Sphere {radius 0.1} \n} \n";
        fout<<"children [ Shape {\n appearance Appearance {\n material Material { emissiveColor 1 1 1 ambientIntensity 1.0}\n}\n geometry Sphere { radius 1.0 } \n} \nDEF Atom"<<i<<"_sensor TouchSensor {}\n]";
        fout<<"},\n";
      }
    }
    //shininess 0.5
    //ambientIntensity 0.5
    //for(int i=0; i< numBodies; ++i) fout<<"DEF Clock"<<i<<" TimeSensor {\n"<<"cycleInterval 2.0 \n"<<"loop TRUE\n}, \n";
    fout<<"DEF Clock TimeSensor {\n"<<"cycleInterval 20.0 \n"<<"loop FALSE\n}, \n";
    for(int i=0; i< numBodies; ++i)
    {
      fout<<"DEF AtomPath"<<i<<" PositionInterpolator {\n"<<"key [ ";
      for(int k=0; k<max_loop; ++k)
      {
        fout<<k*1.0/(max_loop-1)<<" ";
      }
      fout<<"]\n"<<"keyValue [ ";
      for(int j=0; j<max_loop; ++j) fout<<vrmlpos[i][j].x<<" "<<vrmlpos[i][j].y<<" "<<vrmlpos[i][j].z<<"\n";
      fout<<"\n]\n}";
      if(i<numBodies-1) fout<<", \n";
    }
    fout<<"\n]\n}\n";
    for(int i=0; i< numBodies; ++i) fout<<"ROUTE Atom"<<i<<"_sensor.touchTime   TO Clock.set_startTime\n";
    for(int i=0; i< numBodies; ++i) fout<<"ROUTE Clock.fraction_changed   TO AtomPath"<<i<<".set_fraction\n";
    for(int i=0; i< numBodies; ++i) fout<<"ROUTE AtomPath"<<i<<".value_changed   TO Atom"<<i<<".set_translation\n";
  }
  vrmlinit++;
  fout.close();
  //system("freewrl /home/70-gaa/projects/grav_CPU/build-new-Desktop-Default/vrmlpositions");
  
/*
  std::cout<<" CHECK vrmlpos:"<<std::endl;
  for(int i=0; i<max_loop; ++i)
  {
    std::cout<<vrmlpos[0][i].x<<" "<<vrmlpos[0][i].y<<" "<<vrmlpos[0][i].z<<"     "
             <<vrmlpos[1][i].x<<" "<<vrmlpos[1][i].y<<" "<<vrmlpos[1][i].z<<std::endl;
  }
*/
}
#endif

void randomizeBodies()
{
  float kin_en=0.;
  float pot_en=0.;
  
  const double pi=3.1415926536;   
  float RADIUS=pow(30.0*numBodies,1.0/3.0);
  //std::cout<<"R="<<RADIUS<<std::endl;
  //sleep(1);
  P_init={0.0, 0.0, 0.0};

  unsigned int SEED=2000000000;
  
  //Preparing initial seeds:
#ifdef ACCEL
#pragma acc parallel num_gangs(1) vector_length(1) copyin(SEED) present(SEED1,SEED2,SEED3,SEED4,SEED5)
#else
#pragma novector
#endif
{
  for(int i=0; i<numBodies; ++i)
  {
    SEED=Rand32(SEED);
    SEED1[i]=SEED;
    SEED=Rand32(SEED);
    SEED2[i]=SEED;
    SEED=Rand32(SEED);
    SEED3[i]=SEED;
    SEED=Rand32(SEED);
    SEED4[i]=SEED;
    SEED=Rand32(SEED);
    SEED5[i]=SEED;
  }
}

#ifdef DEB123
  //Check 1
#ifdef ACCEL
#pragma acc update host(SEED1,SEED2,SEED3,SEED4,SEED5)
#endif
  for(int i=0; i<numBodies; ++i)
  {
    std::cout<<SEED1[i]<<" "<<SEED2[i]<<" "<<SEED3[i]<<" "<<SEED4[i]<<" "<<SEED5[i]<<std::endl;
  }
#endif
  
  double P_x=0.0, P_y=0.0, P_z=0.0;
#ifdef ACCEL
//#pragma acc parallel loop copy(kin_en,Ep,pi,md,P_x,P_y,P_z) present(p,SEED1,SEED2,SEED3,SEED4,SEED5) //reduction(+:kin_en)
#pragma acc parallel num_gangs(1) vector_length(1) copy(kin_en,Ep,pi,md,P_x,P_y,P_z) present(p,SEED1,SEED2,SEED3,SEED4,SEED5)
#else
//#pragma omp parallel for reduction(+:kin_en,P_x,P_y,P_z)
#endif
  for(int i=0; i<numBodies; ++i)
  {
    //std::cout<<"1. i="<<i<<std::endl;
    float DISTANCE=0.0f;
    bool start=false;
    if(i>0)
    {
      while(DISTANCE<=1.0f)  // 3 fermi
      {
        if(start) //does not work in parallel on GPU.!!!
        {
          SEED1[i]=Rand32(SEED1[i]);
          SEED2[i]=Rand32(SEED2[i]);
          SEED3[i]=Rand32(SEED3[i]);
        }
        float r1=rndv(SEED1[i]);
        float r2=rndv(SEED2[i]);
        float r3=rndv(SEED3[i]);        
        //p.x[i] = XYLAG*(2.0f*static_cast <float> (rand())/(static_cast <float> (RAND_MAX))-1.0f);
        //float theta=pi*static_cast <float> (rand())/(static_cast <float> (RAND_MAX));
        //float phi=2*pi*static_cast <float> (rand())/(static_cast <float> (RAND_MAX));
        //float rr=RADIUS*static_cast <float> (rand())/(static_cast <float> (RAND_MAX));
        float theta=r1*pi;
        float phi=2*pi*r2;
        float rr=RADIUS*r3;        
        p[i].x =rr*sin(theta)*cos(phi);
        p[i].y =rr*sin(theta)*sin(phi);
        p[i].z =rr*cos(theta);
        DISTANCE=10000.0f;
#ifdef ACCEL
#pragma acc loop vector reduction(min:DISTANCE)
#else
#pragma simd reduction(min:DISTANCE)      
#endif
        for(int j=0; j<i; ++j)
        {
          float dij=sqrt((p[i].x-p[j].x)*(p[i].x-p[j].x)+(p[i].y-p[j].y)*(p[i].y-p[j].y)+(p[i].z-p[j].z)*(p[i].z-p[j].z));
          //std::cout<<"dij="<<dij<<std::endl;
          if(dij<DISTANCE) DISTANCE=dij;
          //std::cout<<"1. i="<<i<<",j="<<j<<std::endl;
        }
        start=true;
      }
    }
    else
    {
      float r1=rndv(SEED1[i]);
      float r2=rndv(SEED2[i]);
      float r3=rndv(SEED3[i]);
      float theta=r1*pi;
      float phi=2*pi*r2;
      float rr=RADIUS*r3;
      p[i].x =rr*sin(theta)*cos(phi);
      p[i].y =rr*sin(theta)*sin(phi);
      p[i].z =rr*cos(theta);
    }
    //std::cout<<"2. i="<<i<<std::endl;
    //sleep(1);
    float modv=sqrt(2.0*Ep/md);
    
    float3 res=RandomDirection(SEED4[i], SEED5[i], i);
#ifdef DEB123   
    Sx[i]=res.x;
    Sy[i]=res.y;
    Sz[i]=res.z;
#endif

    float d=sqrtf(res.x*res.x+res.y*res.y+res.z*res.z);   
    //res.x/=d;
    //res.y/=d;
    //res.z/=d;
    //float3 resv=normalize(res,d);               // ???????????????????????????????????????????????

#ifdef DEB123
    Smv1[i]=d;
    Sxn[i]=res.x;
    Syn[i]=res.y;
    Szn[i]=res.z;
#endif
    
    double3 v;
    v.x=modv*res.x;
    v.y=modv*res.y;
    v.z=modv*res.z; 
    //v=-normalize(v)*VELOCITY*10;
    ///p[i].vx =-p[i].x/dist_r*VELOCITY*10;
    ///p[i].vy =-p[i].y/dist_r*VELOCITY*10;
    ///p[i].vz =-p[i].z/dist_r*VELOCITY*10;

    p[i].vx =v.x;
    p[i].vy =v.y;
    p[i].vz =v.z;
    p[i].m=md;  //940;               // deuteron mass // proton mass   (MeV)
    kin_en+=p[i].vx*p[i].vx;
    kin_en+=p[i].vy*p[i].vy;
    kin_en+=p[i].vz*p[i].vz;
#ifdef MAKEVRML
    vrmlpos[i][0].x=p[i].x;
    vrmlpos[i][0].y=p[i].y;
    vrmlpos[i][0].z=p[i].z;
#endif
    P_x+=p[i].vx*p[i].m;
    P_y+=p[i].vy*p[i].m;
    P_z+=p[i].vz*p[i].m;   
  }
#ifdef DEB123
#ifdef ACCEL
#pragma acc update host(p,S1,S2,S3,Sx,Sy,Sz,Sxn,Syn,Szn,Smv1)
#endif 
  for(int i=0; i<numBodies; ++i)
  {
    std::cout<<"1. S: "<<S1[i]<<" "<<S2[i]<<" "<<S3[i]<<std::endl;
  }

  for(int i=0; i<numBodies; ++i)
  {
    std::cout<<"Sxyz: "<<Sx[i]<<" "<<Sy[i]<<" "<<Sz[i]<<std::endl;
  }

  for(int i=0; i<numBodies; ++i)
  {
    std::cout<<"Sxyzn: "<<Sxn[i]<<" "<<Syn[i]<<" "<<Szn[i]<<std::endl;
  }

  for(int i=0; i<numBodies; ++i)
  {
    std::cout<<"Smv1: "<<Smv1[i]<<std::endl;
  }

  std::cout<<"kin_en="<<kin_en<<std::endl;
#endif
  P_x/=numBodies;
  P_y/=numBodies;
  P_z/=numBodies;
  //P_init.x=P_init.x/numBodies;
  //P_init.y=P_init.y/numBodies;
  //P_init.z=P_init.z/numBodies;
#ifdef ACCEL
#pragma acc parallel loop copy(P_x,P_y,P_z)
#else
  //#pragma omp parallel for
#pragma simd
#endif 
  for(int i=0; i<numBodies; ++i)
  {
    p[i].vx-=P_x/p[i].m;
    p[i].vy-=P_y/p[i].m;
    p[i].vz-=P_z/p[i].m;
  }
  
  //P_init={0.0, 0.0, 0.0};
  P_x=0.0, P_y=0.0, P_z=0.0;
  //double P_xin=0.0;
  //double P_yin=0.0;
  //double P_zin=0.0;
  float K_init=0.0;
#ifdef ACCEL
#pragma acc parallel loop copy(P_x,P_y,P_z,K_init,pot_en) 
#else
//#pragma omp parallel for reduction(+:P_x,P_y,P_z,K_init,pot_en)
//#pragma simd reduction(+:P_x,P_y,P_z,K_init,pot_en)
#endif
  for(int i=0; i<numBodies; ++i)
  {
    float Pp=0.0f;
    P_x+=p[i].vx*p[i].m;
    P_y+=p[i].vy*p[i].m;
    P_z+=p[i].vz*p[i].m;
    K_init+=p[i].vx*p[i].vx+p[i].vy*p[i].vy+p[i].vz*p[i].vz;
#ifdef ACCEL
#pragma acc loop vector reduction(+:Pp)
#else
#pragma simd reduction(+:Pp)
#endif    
    for(int j=0; j<numBodies; ++j) if(i!=j)
    {
      float rd=sqrtf((p[i].x-p[j].x)*(p[i].x-p[j].x)+(p[i].y-p[j].y)*(p[i].y-p[j].y)+(p[i].z-p[j].z)*(p[i].z-p[j].z));
      //std::cout<<"rd="<<rd<<std::endl;
      U[i][j]=alpha/rd;
      //if(i==137 && j==645) std::cout<<"@rd="<<rd<<std::endl;
      Pp+=U[i][j];
#ifdef DEB_NOPAR
      if(U[i][j]>MAX_E)
      {
        MAX_E=U[i][j];
        I_E=i;
        J_E=j;
        MIN_d=rd;
      }
#endif
    }
    pot_en+=Pp;
  }
  P_init.x=P_x;
  P_init.y=P_y;
  P_init.z=P_z;
  K_init*=md/2;
  pot_en/=2;
#ifdef DEB123
//Check 2
#ifdef ACCEL
#pragma acc update host(p,S1,S2,S3,Sx,Sy,Sz,Sxn,Syn,Szn,Smv1)
#endif
  for(int i=0; i<numBodies; ++i)
  {
    std::cout<<"2. i="<<i<<" r: "<<p[i].x<<" "<<p[i].y<<" "<<p[i].z<<" v: "<<p[i].vx<<" "<<p[i].vy<<" "<<p[i].vz<<std::endl;
  }
#endif  
  std::cout<<"IN Px="<<P_init.x<<" Py="<<P_init.y<<" Pz="<<P_init.z<<" K="<<K_init<<" P="<<pot_en<<" E="<<K_init+pot_en<<std::endl;
  //std::cout<<"I pi: "<<p[137].x<<" "<<p[137].y<<" "<<p[137].z<<"pj: "<<p[645].x<<" "<<p[645].y<<" "<<p[645].z<<std::endl;

  /*for(int i=0; i<numBodies; ++i)
  {
    for(int j=0; j<numBodies; ++j)
    {
      float r=sqrt((p[i].x-p[j].x)*(p[i].x-p[j].x)+(p[i].y-p[j].y)*(p[i].y-p[j].y)+(p[i].z-p[j].z)*(p[i].z-p[j].z));
      //std::cout<<"r"<<i<<j<<"="<<r<<std::endl;
    }
  }*/
  E_init=K_init+pot_en;
  Kf=K_init;
  Pf=pot_en;
  
/*
  std::cout<<"randomizebodies(): CHECK initial v:"<<std::endl;
  for(int i=0; i<numBodies; ++i) std::cout<<"p["<<i<<"].vx="<<p[i].vx<<" p["<<i<<"].vy="<<p[i].vy
                                          <<" p["<<i<<"].vz="<<p[i].vz<<" |v["<<i<<"]|="
                                          <<sqrt(p[i].vx*p[i].vx+p[i].vy*p[i].vy+p[i].vz*p[i].vz)
                                          <<" sqrt(2*E/md)="<<sqrt(2*E/md)<<std::endl;
  std::cout<<std::endl;
*/
}


/*
#include <array>
std::array<int, 3> arr={1, 2, 3};

void printar(std::array<int, 3> & ar)
{
  int size=ar.size();
  std::cout<<"size="<<size<<std::endl;
  for(int i=0; i<size; ++i) std::cout<<ar[i]<<" ";
  std::cout<<std::endl;
}
*/

void hist()
{
  int Nn=50;
  float Ee=Nn*2.01;
  int BIN[Nn];
  float inite[Nn];
  float fine[Nn];
  float dE=Ee/Nn;
  for(int j=0; j<Nn; ++j) BIN[j]=0;
  for(int j=0; j<Nn; ++j)
  {
    inite[j]=j*dE;
    fine[j]=(j+1)*dE;
  }
  //for(int i=0; i<Nn; ++i)
  //{
  //  std::cout<<"BIN["<<i<<"] inite="<<inite[i]<<" fine="<<fine[i]<<std::endl;
  //}
  //#pragma novector
  for(int i=0; i<numBodies; ++i)
  {
    float ep=p[i].vx*p[i].vx+p[i].vy*p[i].vy+p[i].vz*p[i].vz;
    ep*=md/2;
    //std::cout.precision(10);
    //std::cout<<"ep="<<ep<<std::endl;
    //printf("ep=%.20f", ep);
    //if(ep>2.0) printf("\nERROR: ep=%.20\n",ep);
    //#pragma novector
    for(int k=0; k<Nn; ++k)
    {
      if(ep>inite[k] && ep<=fine[k])
      {
        //#pragma omp atomic
        ++BIN[k];
      }
    }   
  }
  std::cout<<"BIN:"<<std::endl;
  for(int i=0; i<Nn; ++i) std::cout<<BIN[i]<<" ";
  std::cout<<std::endl;
}

void dispersion()
{
  float Ekav=0.0;
  float Ek2av=0.0;
  float DEk=0.0;   // kinetic enrgy dispersion
  float RDE=0.0;
  //float SQRTDE=0.0;
#ifdef ACCEL
#pragma acc parallel loop present(p) copy(Ekav,Ek2av)
#else
#pragma omp parallel for simd reduction(+:Ekav,Ek2av)
#endif
  for(int i=0; i<numBodies; ++i)
  {
    float eki=p[i].m*(p[i].vx*p[i].vx+p[i].vy*p[i].vy+p[i].vz*p[i].vz)/2;
    Ekav+=eki;
    Ek2av+=eki*eki;
  }
  Ekav/=numBodies;
  Ek2av/=numBodies;
  double Ekav2=Ekav*Ekav;
  DEk=Ek2av-Ekav2;
  //std::cout<<"DEk="<<DEk<<std::endl;
  RDE=DEk/Ekav2;
  std::cout<<"RDE="<<RDE<<std::endl;
  //SQRTDE=sqrtf(DEk)/Ekav;
  //std::cout<<"SQRTDE="<<SQRTDE<<std::endl;
}

//////////////////////////////////////////////////////////////////////////////
// Program main
//////////////////////////////////////////////////////////////////////////////
int
main(int argc, char **argv)
{
//handling Not a number exception:
  feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);
  signal(SIGFPE, handler);
//
  //std::cout<<"WE WORK WITH THIS FILE"<<std::endl;
  randomizeBodies();
  static int count=0;
  auto begin=std::chrono::steady_clock::now();
  while(STEP<max_loop)
  {
    //hist();
    dispersion();
    std::cout<<"STEP="<<STEP<<std::endl;
    interaction();
    ++STEP;
    ++count;
  }
  auto end=std::chrono::steady_clock::now();
  auto elapsed_ms=std::chrono::duration_cast<std::chrono::milliseconds>(end-begin);
#ifdef MAKEVRML
  makeVRML();
#endif
  std::cout<<"chrono exe time="<<elapsed_ms.count()<<" ms"<<std::endl;
  //std::cout<<"V="<<VELOCITY<<std::endl;
  //std::cout<<"N="<<numBodies<<" max_loop="<<max_loop<<std::endl;
}

//cpu nonparallel:
//STEP=999
//MAX_NEAR_PAIR=0.552322 I_NEAR_PAIR=212 J_NEAR_PAIR=695
//MAX_NEAR=0.552322 I_NEAR=212 J_NEAR=695
//MAX_vx=0.411507 MAX_vy=0.422237 MAX_vz=0.445483
//MAX_dvx=0.000158254 MAX_dvy=0.000178501 MAX_dvz=0.000323982
//E=40550.8 dE/E=-1.01448e-05 K=37803.1 dK=2.62891 dk=0.00256729 P=2747.75
// dP=-3.04028 dp=-5.80455e-06 Px=0.00919482 Py=-0.0067539 Pz=0.00119138
//chrono exe time=1611747 ms

//gpu:
//STEP=9999
//E=40550.8 dE/E=-7.5578e-06 K=40306.3 dK=-0.28125 dk=-0.000274658 P=244.549
// dP=-0.0252228 dp=-4.81557e-08 Px=0.638783 Py=0.377057 Pz=1.83122
//chrono exe time=227921 ms

//cpu:
//STEP=9999
//E=40550.8 dE/E=-7.75159e-06 K=40286.7 dK=-0.289062 dk=-0.000282288 P=264.069
// dP=-0.0252686 dp=-4.82431e-08 Px=0.263413 Py=0.964223 Pz=0.264539
//chrono exe time=603141 ms

//cpu knl:
//STEP=9999
//MAX_NEAR=0.607688 I_NEAR=755 J_NEAR=1000
//MAX_K=227.651 I_MAX_K=104 J_MAX_K=600
//MIN_K=0.124741 I_MIN_K=89 J_MIN_K=907
//MAX_KUR=1.16104e+06 I_MAX_KUR=104 J_MAX_KUR=600
//MIN_KUR=14.8925 I_MIN_KUR=89 J_MIN_KUR=907
//MAX_vx=0.39368 MAX_vy=0.418785 MAX_vz=0.406448 MAX_dvx=2.27101e-05 MAX_dvy=2.30018e-05 MAX_dvz=1.90463e-05
//E=40550.8 dE/E=-0.000100289 K=40227.4 Kf=40227.4 K-Kf=0.0117188 n=1024 dK=0.0117188 K-Kf-dK=0 dk=1.14441e-05 P=323.422
// dP=-0.0243225 dp=-4.64369e-08 Px=-0.421302 Py=-0.218671 Pz=2.27614 r=1.0001
//chrono exe time=864143 ms
