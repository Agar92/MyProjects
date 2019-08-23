/*
 * Copyright 1993-2015 NVIDIA Corporation.  All rights reserved.
 *
 * Please refer to the NVIDIA end user license agreement (EULA) associated
 * with this source code for terms and conditions that govern your use of
 * this software. Any use, reproduction, disclosure, or distribution of
 * this software and related documentation outside the terms of the EULA
 * is strictly prohibited.
 *
 */

#ifndef __BODYSYSTEM_H__
#define __BODYSYSTEM_H__

#include <algorithm>
#include "ThreeVector.h"

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

template <typename T> struct vec3
{
    typedef float   Type;
}; // dummy
template <>           struct vec3<float>
{
    typedef float3  Type;
};
template <>           struct vec3<double>
{
    typedef double3 Type;
};

template <typename T> struct vec4
{
    typedef float   Type;
}; // dummy
template <>           struct vec4<float>
{
    typedef float4  Type;
};
template <>           struct vec4<double>
{
    typedef double4 Type;
};

extern int numBodies;

class string;

// BodySystem abstract base class
template <typename T>
class BodySystem
{
    public: // methods
        BodySystem(int numBodies) {}
        virtual ~BodySystem() {}

        virtual void loadTipsyFile(const std::string &filename) = 0;

        virtual void update(T deltaTime) = 0;

        virtual void setSoftening(T softening) = 0;
        virtual void setDamping(T damping) = 0;

        virtual T *getArray(BodyArray array) = 0;
        virtual void   setArray(BodyArray array, const T *data) = 0;

        virtual unsigned int getCurrentReadBuffer() const = 0;

        virtual unsigned int getNumBodies() const = 0;

        virtual void   synchronizeThreads() const {};

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

inline float
normalize(float3 &vector)
{
    float dist = sqrtf(vector.x*vector.x + vector.y*vector.y + vector.z*vector.z);

    if (dist > 1e-6)
    {
        vector.x /= dist;
        vector.y /= dist;
        vector.z /= dist;
    }

    return dist;
}

inline float
dot(float3 v0, float3 v1)
{
    return v0.x*v1.x+v0.y*v1.y+v0.z*v1.z;
}

inline float3
cross(float3 v0, float3 v1)
{
    float3 rt;
    rt.x = v0.y*v1.z-v0.z*v1.y;
    rt.y = v0.z*v1.x-v0.x*v1.z;
    rt.z = v0.x*v1.y-v0.y*v1.x;
    return rt;
}

inline static double RND()
{
  double r = static_cast <double> (rand())/(static_cast <double> (RAND_MAX));
  return r;
}


// utility function
template <typename T>
void randomizeBodies(NBodyConfig config, T *pos, T *vel, float *color, float clusterScale,
                     float velocityScale, int numBodies, float damping,  bool vec4vel)
{
    switch (config)
    {
        default:
        case NBODY_CONFIG_RANDOM:
            {
                float scale = clusterScale * std::max<float>(1.0f, numBodies / (1024.0f));
                float vscale = velocityScale * scale;

                int p = 0, v = 0;
                int i = 0;

                while (i < numBodies)
                {
                    float3 point;
                    //const int scale = 16;
                    point.x = rand() / (float) RAND_MAX * 2 - 1;
                    point.y = rand() / (float) RAND_MAX * 2 - 1;
                    point.z = rand() / (float) RAND_MAX * 2 - 1;
                    float lenSqr = dot(point, point);

                    if (lenSqr > 1)
                        continue;

                    float3 velocity;
                    velocity.x = rand() / (float) RAND_MAX * 2 - 1;
                    velocity.y = rand() / (float) RAND_MAX * 2 - 1;
                    velocity.z = rand() / (float) RAND_MAX * 2 - 1;
                    lenSqr = dot(velocity, velocity);

                    if (lenSqr > 1)
                        continue;

                    pos[p++] = point.x * scale; // pos.x
                    pos[p++] = point.y * scale; // pos.y
                    pos[p++] = point.z * scale; // pos.z
                    pos[p++] = 1.0f; // mass

                    vel[v++] = velocity.x * vscale; // pos.x
                    vel[v++] = velocity.y * vscale; // pos.x
                    vel[v++] = velocity.z * vscale; // pos.x

                    if (vec4vel) vel[v++] = 1.0f; // inverse mass

                    i++;
                }
            }
            break;

        case NBODY_CONFIG_SHELL:
            {
            	float mass=0.01f;
            	float MASS=10000.0f*numBodies*mass;
                /*float scale = clusterScale;
                float vscale = scale * velocityScale;
                float inner = 2.5f * scale;
                float outer = 4.0f * scale;*/

                int p = 0, v=0;
                float3 axis1 = {0, 0, 1.};
                float3 axis2 = {0, 1., 0.};
                int nb2 = numBodies/2;
                //float s = 1.0;
                float d = 10.;
                float b = 3.;
                float u = 500.;

			    pos[p++] =  d ;
			    pos[p++] = -b;
			    pos[p++] =  0;
			    //std::cout<<"first x="<<d<<", y="<<-b<<", vz="<<-u<<std::endl;
			    //mass
			    pos[p++] = MASS;

			    vel[v++] = -u;
			    vel[v++] = 0;
			    vel[v++] = 0;
			    vel[v++] = MASS;
                for(int i=1; i < nb2; ++i)//for(int i=0; i < numBodies; i++)
                {
                	float number1 = RND();
				    float r = clusterScale*sqrt(number1);
                	//float r = clusterScale*number1;
                	float number2 = RND();
                	float phi = 2*M_PI*number2;
                	float x = r*sin(phi);
                	float y = r*cos(phi);
                	//std::cout<<"x="<<x<<", y="<<y<<std::endl;
                  pos[p++] =  d + x;
				  pos[p++] =  -b + y;
				  pos[p++] =  0;
				  //mass
                  pos[p++] = mass;

			      float3 vv = {x, y, 0.};
			      //float k = std::sqrt(numBodies/(r))/clusterScale;
			      //float k = s*std::sqrt(nb2*damping/clusterScale)/r;
			      float k = std::sqrt((MASS+r*r*numBodies*mass/(clusterScale*clusterScale))/r)/r;
				  float3 velocity = cross(vv, axis1);
				  vel[v++] = velocity.x*k-u;
				  vel[v++] = velocity.y*k;
				  vel[v++] = velocity.z*k;
				  vel[v++] = mass;
				  //std::cout<<"f="<<i<<", x="<<d+x<<", y="<<-b+y<<", vx="<<velocity.x*k-u<<", vy="<<velocity.y*k<<", vz="<<velocity.z*k<<std::endl;
                }

			    pos[p++] =  -d ;
			    pos[p++] = 0.;
			    pos[p++] =  b;
			    //mass
			    pos[p++] = MASS;
			    //std::cout<<"second x="<<-d<<", z="<<b<<", vx="<<+u<<std::endl;

			    vel[v++] = u;
			    vel[v++] = 0;
			    vel[v++] = 0;
			    vel[v++] = MASS;
                for(int i=nb2+1; i < numBodies; ++i)//for(int i=0; i < numBodies; i++)
				{
                	float number1 = RND();
				    float r = clusterScale*sqrt(number1);
				    //float r = clusterScale*number1;
				    float number2 = RND();
				    float phi = 2*M_PI*number2;
				    float x = r*sin(phi);
				    float z = r*cos(phi);
				    pos[p++] =  -d + x;
				    pos[p++] =   0.;
				    pos[p++] =  b+z;
				    //mass
				    pos[p++] = mass;

				    float3 vv = {x, 0., z};
				    //float k = std::sqrt(numBodies/(r))/clusterScale;
				    //float k = s*std::sqrt(nb2*damping/clusterScale)/r;
				    float k = std::sqrt((MASS+r*r*numBodies*mass/(clusterScale*clusterScale))/r)/r;
				    float3 velocity = cross(vv, axis2);
				    vel[v++] = velocity.x*k+u;
				    vel[v++] = velocity.y*k;
				    vel[v++] = velocity.z*k;
				    vel[v++] = mass;
				    //std::cout<<"s="<<i<<", x="<<-d+x<<", z="<<b+z<<", vx="<<velocity.x*k+u<<", vy="<<velocity.y*k<<", vz="<<velocity.z*k<<std::endl;
				}
            }
            break;

        case NBODY_CONFIG_EXPAND:
            {
                float scale = clusterScale * numBodies / (1024.f);

                if (scale < 1.0f)
                    scale = clusterScale;

                float vscale = scale * velocityScale;

                int p = 0, v = 0;

                for (int i=0; i < numBodies;)
                {
                    float3 point;

                    point.x = rand() / (float) RAND_MAX * 2 - 1;
                    point.y = rand() / (float) RAND_MAX * 2 - 1;
                    point.z = rand() / (float) RAND_MAX * 2 - 1;

                    float lenSqr = dot(point, point);

                    if (lenSqr > 1)
                        continue;

                    pos[p++] = point.x * scale; // pos.x
                    pos[p++] = point.y * scale; // pos.y
                    pos[p++] = point.z * scale; // pos.z
                    pos[p++] = 1.0f; // mass
                    vel[v++] = point.x * vscale; // pos.x
                    vel[v++] = point.y * vscale; // pos.x
                    vel[v++] = point.z * vscale; // pos.x

                    if (vec4vel) vel[v++] = 1.0f; // inverse mass

                    i++;
                }
            }
            break;
    }

    if (color)
    {
        int v = 0;

        for (int i=0; i < numBodies; i++)
        {
            //const int scale = 16;
            color[v++] = rand() / (float) RAND_MAX;
            color[v++] = rand() / (float) RAND_MAX;
            color[v++] = rand() / (float) RAND_MAX;
            color[v++] = 1.0f;
        }
    }

}



#endif // __BODYSYSTEM_H__
