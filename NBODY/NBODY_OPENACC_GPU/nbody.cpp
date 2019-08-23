#include <stdlib.h> 
#include <stdio.h> 
#include <sys/time.h> 
#include <math.h>
#include <unistd.h>
#include <random>
#include <memory.h>
#include <cstdio>
#include <cstdlib>
#include <chrono>
#include <iostream>
#include <curand.h>

#define ACCEL

#ifdef ACCEL
/* OpenACC math library */ 
#include <accelmath.h>
#include <openacc.h>
#else
#include <omp.h>
#endif

//cmake ../ -DCMAKE_CXX_COMPILER=pgc++ -DCMAKE_C_COMPILER=pgcc -DCMAKE_C_FLAGS="-acc -Minline -Minfo -ta=nvidia,kepler,keepptx,keepgpu"
//pgcc -acc -Minline -Minfo -ta=nvidia,kepler,keepptx,keepgpu -o nbody nbody.c

//launch:
//cmake ../ -DCMAKE_CXX_COMPILER=pgc++ -DCMAKE_CXX_FLAGS=" -acc -Minline -Minfo=accel -ta=tesla:cc30"

const float CLUSTERSCALE=100.0f;
const float b=CLUSTERSCALE/2;
const float MASS = 1.0f;
const float DELTATIME = 0.016f;
const float HDELTA=DELTATIME/2;
const float SOFTENING = 0.1f;
const float Temperature =0.1f;
const int numIterations = 10;
const int n = 16384;
int counter=0;

typedef struct
{ 
  double x; 
  double y; 
  double z; 
} Point;

typedef struct
{ 
  double x; 
  double y; 
  double z;
  double w;
} Point4;

Point4  positions   [n];
Point4  velocities  [n];
Point4  acceleration[n];
float   E_pot       [n];
#pragma acc declare device_resident(positions,velocities,acceleration,E_pot)

#ifdef ACCEL
#pragma acc routine seq
#endif
unsigned int Rand32(unsigned int xn)
{
  u_quad_t a=0x5DEECE66D;
  u_quad_t c=0xB;
  return (unsigned int)((a*xn+c) & 0xFFFFFFFF);
}

#ifdef ACCEL
#pragma acc routine seq
#endif
double rndv(unsigned int xn)
{
  return (double) xn / (double) 0x100000000LL;
}

#ifdef ACCEL
#pragma acc routine seq
#endif
Point RandomDirection(float r1, float r2)
{
	double z  =2.0*r1-1.;  // z = cos(theta)
	double rho=sqrt((1.+z)*(1.-z)); // rho = sqrt(1-z*z)
  double phi = M_PI*2*r2;
  Point point;
  point.x=rho*cosf((float)phi);
  point.y=rho*sinf((float)phi);
  point.z=z;
	return point;
}
#ifdef ACCEL
#pragma acc routine seq
#endif
double RndMaxwell(float r1, float r2, float r3) // TheResult must be multiplied byTemperature
{
	static const double halfpi = M_PI / 2;
	const double R =cos((float)(r1 * halfpi));
	return -log((float)r2) - log((float)r3) * R * R;
}
#ifdef ACCEL
#pragma acc routine seq
#endif
Point scalevec(Point vector, double scalar)
{
	vector.x *= scalar;
	vector.y *= scalar;
	vector.z *= scalar;
  return vector;
}

void randomizeBodies()
{
  float kin_energy=0;
#ifdef ACCEL
#pragma acc parallel loop reduction(+:kin_energy)
#else
#pragma simd reduction(+:kin_energy)
#endif
  for(int i=0;i<n; ++i)
  {
    float x, y, z, r1, r2, r3, r4, r5;
    //x = b*(1.0-2*(static_cast <float> (rand())/(static_cast <float> (RAND_MAX))));
    //y = b*(1.0-2*(static_cast <float> (rand())/(static_cast <float> (RAND_MAX))));
    //z = b*(1.0-2*(static_cast <float> (rand())/(static_cast <float> (RAND_MAX))));
    int indexi=5*i;
    unsigned int a1=Rand32(indexi);
    x=b*(1.0-2*(rndv(a1)));
    unsigned int a2=Rand32(indexi+1);
    y=b*(1.0-2*(rndv(a2)));
    unsigned int a3=Rand32(indexi+2);
    z=b*(1.0-2*(rndv(a3)));
    unsigned int a4=Rand32(indexi+3);
    r1=rndv(a4);
    unsigned int a5=Rand32(indexi+4);
    r2=rndv(a5);
    unsigned int a6=Rand32(indexi+5);
    r3=rndv(a6);
    unsigned int a7=Rand32(indexi+6);
    r4=rndv(a7);
    unsigned int a8=Rand32(indexi+7);
    r5=rndv(a8);
    
    positions[i].x=x;
    positions[i].y=y;
    positions[i].z=z;
    positions[i].w=MASS;

    double modvelocity = (double)sqrtf((2*(RndMaxwell(r1,r2,r3))*(Temperature/MASS)));
    Point velocity=scalevec((RandomDirection(r4,r5)),modvelocity);
    kin_energy+=velocity.x*velocity.x;
    kin_energy+=velocity.y*velocity.y;
    kin_energy+=velocity.z*velocity.z;
    velocities[i].x = velocity.x ;
    velocities[i].y = velocity.y ;
    velocities[i].z = velocity.z ;
    //printf("vx=%f",velocities[i].x);
  }
  kin_energy*=MASS;
  std::cout<<"\ninit_K="<<kin_energy<<std::endl;
}

float calcAccels()
{
#ifdef ACCEL
#pragma acc parallel loop //independent
#else
#pragma omp parallel for
#endif
  for(int i = 0; i < n; ++i)
    {
      float Fx=0.0, Fy=0.0, Fz=0.0;
      float E_p=0.0;

// Uncommenting this pragma exposes more parallelism, but in my test it actually slows things down.  
//#ifdef ACCEL
//#pragma acc kernels loop
//#endif
//#pragma acc kernels loop gang vector                    \
      //device_type(acc_device_nvidia) vector_length(128) \
      //device_type(acc_device_nvidia) vector_length(128)
#ifdef ACCEL
#pragma acc loop vector
#else
#pragma  simd reduction(+:Fx,Fy,Fz,E_p)
#endif
      for (int j=0; j<n; ++j)
      {
        float dx=positions[j].x - positions[i].x;
        float dy=positions[j].y - positions[i].y;
        float dz=positions[j].z - positions[i].z;
        const float drSquared = dx*dx+dy*dy+dz*dz+SOFTENING;
        const float rdistance=1.0f/sqrtf(drSquared);
        if(i!=j) E_p+=-rdistance;
        const float drPower32 = rdistance/drSquared;
        Fx += dx*drPower32;
        Fy += dy*drPower32;
        Fz += dz*drPower32;
      }
      acceleration[i].x = Fx;
      acceleration[i].y = Fy; 
      acceleration[i].z = Fz;
      E_pot[i] = E_p;
    }
    
    float potential = 0.0;    
#ifdef ACCEL  
#pragma acc parallel loop reduction(+:potential)
#else
#pragma omp simd reduction(+:potential)
#endif
    for(int k=0; k<n; ++k)
    {
      potential += E_pot[k];
    }
    return potential;
 } 

void update()
{
   float P = 0.0;
   float K = 0.0;
#ifdef ACCEL
#pragma acc parallel loop present(positions,velocities,acceleration)
#else
#pragma omp simd
#endif
   for(int j=0; j<n; ++j)
   {
     Point4 position = positions    [j];
     Point4 velocity = velocities   [j];
     velocity.x+=acceleration[j].x*HDELTA;
     velocity.y+=acceleration[j].y*HDELTA;
     velocity.z+=acceleration[j].z*HDELTA;
     velocities[j] = velocity;

     position.x += velocity.x*DELTATIME;
     position.y += velocity.y*DELTATIME;
     position.z += velocity.z*DELTATIME;
     positions[j] = position;
   }
   
   P = calcAccels();
   
#ifdef ACCEL
#pragma acc parallel loop //independent
#else
#pragma omp simd
#endif
   for(int j=0; j<n; ++j)
   {
     Point4 velocity = velocities  [j];
     velocity.x     += acceleration[j].x*HDELTA;
     velocity.y     += acceleration[j].y*HDELTA;
     velocity.z     += acceleration[j].z*HDELTA;
     velocities[j]   = velocity;
   }
#ifdef ACCEL
#pragma acc parallel loop reduction(+:K)
#else
#pragma omp simd reduction(+:K)
#endif
   for(int i=0; i<n; ++i)
   {
     K += velocities[i].x*velocities[i].x+velocities[i].y*velocities[i].y+velocities[i].z*velocities[i].z;
   }
   K*=MASS;
   float E = K + P;
   printf("count=%i, n=%i, E=%f, K=%f, P=%f\n", counter, n, E, K, P);
   ++counter;
}
   
int main(int argc, char **argv)
{
  printf("nbody_acc\n");
  auto start=std::chrono::steady_clock::now();
  randomizeBodies();
  for(int i=0; i<numIterations; ++i) update();
  auto final=std::chrono::steady_clock::now();
  auto elapsed_ms=std::chrono::duration_cast<std::chrono::milliseconds>(final-start);
  std::cout<<"chrono exe time="<<elapsed_ms.count()<<" ms"<<std::endl;
  printf("Good Bye\n");
  return 0;
}
