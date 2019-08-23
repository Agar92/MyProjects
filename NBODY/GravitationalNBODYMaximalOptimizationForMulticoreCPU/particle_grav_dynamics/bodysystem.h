#ifndef BODYSYSTEM_H
#define BODYSYSTEM_H

#include <algorithm>
#include<unistd.h>
#include "threevector.h"

extern float CLUSTERSCALE;
extern float HALF;
extern float MASS;
extern int numBodies;

enum NBodyConfig
{
    NBODY_CONFIG_RANDOM,
    NBODY_CONFIG_SHELL,
    NBODY_CONFIG_EXPAND,
    NBODY_NUM_CONFIGS
};

enum BodyArray
{
    BODYSYSTEM_POSITION,
    BODYSYSTEM_VELOCITY,
};

struct float3
{
    float x, y, z;
};
struct double3
{
    double x, y, z;
};
struct ParticleType
{
    float *x;
    float *y;
    float *z;
    float *vx;
    float *vy;
    float *vz;
};




class string;

// BodySystem abstract base class
template <typename T>
class BodySystem
{
    public: // methods
        BodySystem(int numBodies) {}
        virtual ~BodySystem() {}

        virtual void update(T deltaTime) = 0;

        virtual void setSoftening(T softening) = 0;
        virtual void setDamping(T damping) = 0;

        virtual T *getArray(BodyArray array) = 0;
        virtual void   setArray(BodyArray array, const T *data) = 0;

        virtual unsigned int getCurrentReadBuffer() const = 0;

        virtual unsigned int getNumBodies() const = 0;

        virtual void   synchronizeThreads() const {}

    protected: // methods
        BodySystem() {} // default constructor

        virtual void _initialize(int numBodies) = 0;
        virtual void _finalize() = 0;
};

inline float3
scalevec(float3 &vector, float scalar)
{
    float3 rt = vector;
    rt.x *= scalar;
    rt.y *= scalar;
    rt.z *= scalar;
    return rt;
}

inline ThreeVector
_scalevec(ThreeVector vector, double scalar)
{
    vector.fx *= scalar;
    vector.fy *= scalar;
    vector.fz *= scalar;
    return vector;
}
// utility function
template <typename T>
//void randomizeBodies(T *pos, T *vel, float *color, float clusterScale,
//                     float Temperature, int numBodies, bool vec4vel)
void randomizeBodies(ParticleType *p, float clusterScale,
                     float Temperature, int numBodies, bool vec4vel)
{
                float scale = 1;
                float vscale = 1;
                float inner = 2.5f * scale;
                float outer = 4.0f * scale;

                int px = 0, py = 0, pz = 0, vx=0, vy=0, vz=0;
                int i = 0;
                T kin_en=0;

                while (i < numBodies)//for(int i=0; i < numBodies; i++)
                {
                    float x, y, z;
                    x = HALF*(1.0-2*(static_cast <float> (rand())/(static_cast <float> (RAND_MAX))));
                    y = HALF*(1.0-2*(static_cast <float> (rand())/(static_cast <float> (RAND_MAX))));
                    z = HALF*(1.0-2*(static_cast <float> (rand())/(static_cast <float> (RAND_MAX))));
                    float3 point = {x, y, z};

                    p->x[i] =  point.x ;
                    p->y[i] =  point.y ;
                    p->z[i] =  point.z ;         
                    double modvelocity = (double)sqrt((2*(ThreeVector::RndMaxwell())*(Temperature/MASS)));
                    ThreeVector container =_scalevec((ThreeVector::RandomDirection()), modvelocity);
                    double3 velocity;
                    velocity.x = container.fx;
                    kin_en+=velocity.x*velocity.x;
                    velocity.y = container.fy;
                    kin_en+=velocity.y*velocity.y;
                    velocity.z = container.fz;
                    kin_en+=velocity.z*velocity.z;
		    
                    p->vx[i] = velocity.x ;
                    p->vy[i] = velocity.y ;
                    p->vz[i] = velocity.z ;  
		    /*std::cout<<"x="<<p[i].x<<",y="<<p[i].y<<",z="<<p[i].z<<"vx="<<p[i].vx<<",vy="<<p[i].vy<<",vz="<<p[i].vz<<std::endl;*/                        
                    i++;
                }
                kin_en*=MASS;
                std::cout<<"K="<<kin_en<<std::endl;
}

#endif // BODYSYSTEM_H
