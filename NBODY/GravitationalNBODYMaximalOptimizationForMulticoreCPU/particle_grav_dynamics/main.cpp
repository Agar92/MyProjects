#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <algorithm>
#include <assert.h>
#include <math.h>
#include "bodysystemcpu.h"
#include "stdlib.h"
#include "unistd.h"
#include <chrono>

float CLUSTERSCALE=100.0f;
float HALF=CLUSTERSCALE/2;
float SOFTENING=0.1f;
float MASS=1.0f;

float SIGMA=3.1f;
float EPSILON=0.03f;

bool useCpu = true;
bool bSupportDouble = true;
int flopsPerInteraction = 20;

int numBodies = 8192;

////////////////////////////////////////
// Demo Parameters
////////////////////////////////////////
struct NBodyParams
{
    float m_timestep;
    float m_clusterScale;
    float m_Temperature;
    float m_softening;
    float m_damping;
    float m_pointSize;
    float m_x, m_y, m_z;

    void print()
    {
        printf("{ %f, %f, %f, %f, %f, %f, %f, %f, %f },\n",
               m_timestep, m_clusterScale, m_Temperature,
               m_softening, m_damping, m_pointSize, m_x, m_y, m_z);
    }
};

NBodyParams demoParams[] =
{
    { 0.016f, 100.f, 0.1f, 0.1f, 1.0f, 1.0f, 0, -2, -100}/*,
    { 0.016f, 0.68f, 20.0f, 0.1f, 1.0f, 0.8f, 0, -2, -30},
    { 0.0006f, 0.16f, 1000.0f, 1.0f, 1.0f, 0.07f, 0, 0, -1.5f},
    { 0.0006f, 0.16f, 1000.0f, 1.0f, 1.0f, 0.07f, 0, 0, -1.5f},
    { 0.0019f, 0.32f, 276.0f, 1.0f, 1.0f, 0.07f, 0, 0, -5},
    { 0.0016f, 0.32f, 272.0f, 0.145f, 1.0f, 0.08f, 0, 0, -5},
    { 0.016000f, 6.040000f, 0.000000f, 1.000000f, 1.000000f, 0.760000f, 0, 0, -50},*/
};

int numDemos = sizeof(demoParams) / sizeof(NBodyParams);
int activeDemo = 0;
float demoTime = 10000.0f; // ms

// run multiple iterations to compute an average sort time

NBodyParams activeParams = demoParams[activeDemo];
template <typename T>
class NBodyDemo
{
    public:
        static void Create()
        {
            m_singleton = new NBodyDemo;
            std::cout<<m_singleton<<std::endl;
        }
        static void Destroy()
        {
            delete m_singleton;
        }

        static void init(int numBodies)
        {
            m_singleton->_init(numBodies);
        }

        static void reset(int numBodies)
        {
            m_singleton->_reset(numBodies);
        }

        static void selectDemo(int index)
        {
            m_singleton->_selectDemo(index);
        }

        static void updateSimulation()
        {
	    std::printf("FELL");
 	    //m_singleton->m_nbody->update(activeParams.m_timestep);
            m_singleton->m_nbodyCpu->_update(m_singleton->m_particles, activeParams.m_timestep);
        }

        static void display()
        {
	    static int count=0;
            //m_singleton->m_renderer->setSpriteSize(activeParams.m_pointSize);
            //m_singleton->m_renderer->setPositions(
	    std::cout<<"count="<<count<<std::endl;
            //for(int i=0;i<numBodies;i++){
                //int index=4*i;
                /*std::cout<<"x"<<i<<"="<<m_singleton->m_particles[i].x<<",y"<<i<<"="<<m_singleton->m_particles[i].y<<",z"<<i<<"="<<m_singleton->m_particles[i].z<<",vx="<<m_singleton->m_particles[i].vx<<",vy="<<m_singleton->m_particles[i].vy<<",vz="<<m_singleton->m_particles[i].vz<<std::endl;*/
	   //}
	    count++;
        }

        static void getArrays(T *pos, T *vel)
        {
            T *_pos = m_singleton->m_nbody->getArray(BODYSYSTEM_POSITION);
            T *_vel = m_singleton->m_nbody->getArray(BODYSYSTEM_VELOCITY);
            memcpy(pos, _pos, m_singleton->m_nbody->getNumBodies() * 4 * sizeof(T));
            memcpy(vel, _vel, m_singleton->m_nbody->getNumBodies() * 4 * sizeof(T));
        }

        static void setArrays(const T *pos, const T *vel)
        {
            if (pos != m_singleton->m_hPos)
            {
                memcpy(m_singleton->m_hPos, pos, numBodies * 4 * sizeof(T));
            }

            if (vel != m_singleton->m_hVel)
            {
                memcpy(m_singleton->m_hVel, vel, numBodies * 4 * sizeof(T));
            }
            m_singleton->m_nbody->setArray(BODYSYSTEM_POSITION, m_singleton->m_hPos);
            m_singleton->m_nbody->setArray(BODYSYSTEM_VELOCITY, m_singleton->m_hVel);
        }

    private:
        static NBodyDemo *m_singleton;
	ParticleType *m_particles;
        BodySystem<T>     *m_nbody;
        BodySystemCPU<T>  *m_nbodyCpu;

        T *m_hPos;
        T *m_hVel;
        float *m_hColor;

    private:
        NBodyDemo()
            : m_nbody(0),
              m_nbodyCpu(0),
              m_hPos(0),
              m_hVel(0),
              m_hColor(0)
        {

        }

        ~NBodyDemo()
        {
            if (m_nbodyCpu)
            {
                delete m_nbodyCpu;
            }

            if (m_hPos)
            {
                delete [] m_hPos;
            }

            if (m_hVel)
            {
                delete [] m_hVel;
            }

            if (m_hColor)
            {
                delete [] m_hColor;
            }
        }

        void _init(int numBodies)
        {
            m_nbodyCpu = new BodySystemCPU<T>(numBodies);
            m_nbody = m_nbodyCpu;
            // allocate host memory
            m_hPos = new T[numBodies*4];
            m_hVel = new T[numBodies*4];
	    std::cout<<"FAIL"<<std::endl;
	    m_particles=(ParticleType *)malloc(sizeof(ParticleType));
	    m_particles->x=new float[numBodies];
	    m_particles->y=new float[numBodies];
	    m_particles->z=new float[numBodies];
	    m_particles->vx=new float[numBodies];
	    m_particles->vy=new float[numBodies];
	    m_particles->vz=new float[numBodies];
	    //m_particles.x = new float[numBodies];
	    //m_particles.y = new float[numBodies];
	    //m_particles.z = new float[numBodies];
	    //m_particles.vx = new float[numBodies];
	    //m_particles->vy = new float[numBodies];
	    //m_particles->vz = new float[numBodies];
            m_hColor = new float[numBodies*4];
            m_nbody->setSoftening(activeParams.m_softening);
            m_nbody->setDamping(activeParams.m_damping);
        }

        void _reset(int numBodies)
        {
                randomizeBodies<float>(m_singleton->m_particles,
                                activeParams.m_clusterScale,
                                activeParams.m_Temperature,
                                numBodies, true);
                //setArrays(m_hPos, m_hVel);
        }

        void _selectDemo(int index)
        {
            assert(index < numDemos);
            activeParams = demoParams[index];
            reset(numBodies);
        }
};

void finalize()
{
    NBodyDemo<float>::Destroy();
    //if (bSupportDouble) NBodyDemo<double>::Destroy();
}

//template <> NBodyDemo<double> *NBodyDemo<double>::m_singleton = 0;
template <> NBodyDemo<float> *NBodyDemo<float>::m_singleton = 0;

void selectDemo(int activeDemo)
{
        NBodyDemo<float>::selectDemo(activeDemo);
}

void updateSimulation()
{
        NBodyDemo<float>::updateSimulation();
}

void displayNBodySystem()
{
        NBodyDemo<float>::display();
}
void display()
{
    updateSimulation();
    displayNBodySystem();
}

//////////////////////////////////////////////////////////////////////////////
// Program main
//////////////////////////////////////////////////////////////////////////////
int
main(int argc, char **argv)
{
    clock_t start, stop;
    start=clock();
    auto begin=std::chrono::steady_clock::now();
    std::cout<<"HALLO"<<std::endl;
    printf("Run \"nbody -benchmark [-numbodies=<numBodies>]\" to measure performance.\n");
    printf("NOTE: The CUDA Samples are not meant for performance measurements. Results may vary when GPU Boost is enabled.\n\n");
    NBodyDemo<float>::Create();

    NBodyDemo<float>::init(numBodies);
    NBodyDemo<float>::reset(numBodies);
    std::cout<<"ERTO"<<std::endl;
    //sleep(1);

    printf("WE ARE HERE !");
    /*NBodyDemo<double>::Create();
    NBodyDemo<double>::init(numBodies);
    NBodyDemo<double>::reset(numBodies);*/
    static int count=0;
    while(count<10)
    {
	//std::cout<<"ERTO"<<std::endl;
	//sleep(1);
	display();
	//std::cout<<count<<std::endl;
	count++;
    }
    finalize();
    std::cout<<"GOOD BYE"<<std::endl;
    double time=1000*(clock()-start)/(double) CLOCKS_PER_SEC;
    std::cout<<"exe time="<<time<<" ms"<<std::endl;
    auto end=std::chrono::steady_clock::now();
    auto elapsed_ms=std::chrono::duration_cast<std::chrono::milliseconds>(end-begin);
    std::cout<<"chrono exe time="<<elapsed_ms.count()<<" ms"<<std::endl;
    //std::cout<<"FAIL HERE2"<<std::endl;

}
