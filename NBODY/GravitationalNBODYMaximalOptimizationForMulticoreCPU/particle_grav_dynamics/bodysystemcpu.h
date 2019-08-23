#ifndef BODYSYSTEMCPU_H
#define BODYSYSTEMCPU_H

#include "bodysystem.h"

// CPU Body System
template <typename T>
class BodySystemCPU : public BodySystem<T>
{
    public:
        BodySystemCPU(int numBodies);
        virtual ~BodySystemCPU();

        virtual void update(T deltaTime);
	virtual void _update(ParticleType *p, T deltaTime);

        virtual void setSoftening(T softening)
        {
            m_softeningSquared = softening * softening;
        }
        virtual void setDamping(T damping)
        {
            m_damping = damping;
        }

        virtual T *getArray(BodyArray array);
        virtual void   setArray(BodyArray array, const T *data);

        virtual unsigned int getCurrentReadBuffer() const
        {
            return 0;
        }

        virtual unsigned int getNumBodies() const
        {
            return m_numBodies;
        }

    protected: // methods
        BodySystemCPU() {} // default constructor

        virtual void _initialize(int numBodies);
        virtual void _finalize();

        void all_particle_vrml_import();
        void _computeNBodyGravitation();
        void _integrateNBodySystem(T deltaTime);
	//void _MoveParticles(int numBodies, ParticleType *p,float deltatime);
	double _CalcAccel(ParticleType *p, float deltatime);
        void _integrateNParticles(ParticleType *p, float deltatime);

    protected: // data
        int m_numBodies;
        bool m_bInitialized;

        T *m_pos;
        T *m_vel;
	ParticleType *m_aos;
        T *m_force;
        T *m_acc_container;
//V*(t)=V(t)+a(t)*dt/2
        T *m_vel_help;

        T m_softeningSquared;
        T m_damping;
};

#include "bodysystemcpu_impl.h"


#endif // BODYSYSTEMCPU_H
