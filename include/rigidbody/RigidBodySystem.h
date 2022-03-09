#pragma once

#include <Eigen/Dense>

#include <memory>
#include <vector>

class Contact;
class CollisionDetect;
class Solver;
class RigidBody;

typedef std::function<void(std::vector<RigidBody*>&)> PreStepFunc;
typedef std::function<void()> ResetFunc;

enum eSolverType { kPGS = 0, kBPP };
enum eSamplingType { kCorners = 0, kGrid, kSampling };
enum eFrictionDistribution { kUniform = 0, kGaussian, kRamp, kStep };

class RigidBodySystem
{
public:

    RigidBodySystem();

    virtual ~RigidBodySystem();

    // Advance the simulation.
    void step(float _dt);

    // Remove all rigid bodies and cleanup the memory.
    void clear();

    // Add rigid body @a _b to the system. The rigid body is owned and managed by this RigidBodySystem.
    void addBody(RigidBody* _b);

    // Accessors for the body array.
    const std::vector<RigidBody*>& getBodies() const { return m_bodies; }
    std::vector<RigidBody*>& getBodies() { return m_bodies; }

    // Accessors for the contact array.
    const std::vector<Contact*>& getContacts() const;
    std::vector<Contact*>& getContacts();    
    const std::vector<Contact*>& getSubContacts() const;
    std::vector<Contact*>& getSubContacts();

    // Frame counter.
    int getFrameCount() const { return m_frame; }

    // Contact solver params
    void setContactStiffness(float _stiffness) { m_contactStiffness = _stiffness; }
    void setContactDamping(float _damping) { m_contactDamping = _damping; }
    void setFrictionCoefficient(float _mu) { m_mu = _mu; }
    float getFrictionCoefficient() { return m_mu; }
    void setSolverIterations(int _solverIter) { m_solverIter = _solverIter; }

    void setSolverType(eSolverType _solverType) { m_solverType = _solverType; }
    void setSamplingType(eSamplingType _samplingType) { m_samplingType = _samplingType; }
    eSamplingType getSamplingType() { return m_samplingType; }
    void setFrictionDistribution(eFrictionDistribution _frictionDistribution) { m_frictionDistribution = _frictionDistribution; }
    eFrictionDistribution getFrictionDistribution() { return m_frictionDistribution; }

    // Callbacks.
    void setPreStepFunc(PreStepFunc _func) { m_preStepFunc = _func; }
    void setResetFunc(ResetFunc _func) { m_resetFunc = _func; }

private:

    std::vector<RigidBody*> m_bodies;
    std::unique_ptr<CollisionDetect> m_collisionDetect;
    float m_contactStiffness;
    float m_contactDamping;
    float m_mu;
    int m_solverIter;
    eSolverType m_solverType;
    int m_frame;
    eSamplingType m_samplingType;
    eFrictionDistribution m_frictionDistribution;

    // Compute the world-space inverse inertia matrices for all bodies.
    //  This function also updates the rotation matrices using the quaternions.
    void computeInertias();

    // Compute the constraint forces.
    void calcConstraintForces(float dt);

    // Compute the friction coefficient based on the baseline coefficient m_mu
    // the friction distribution m_frictionDistribuation, and the contact position pos
    float computeFrictionCoefficient(const Eigen::Vector3f& pos);

    PreStepFunc m_preStepFunc;
    ResetFunc m_resetFunc;
};


template<typename ScalarType, int Options>
Eigen::Quaternion<ScalarType, Options> operator*(ScalarType a, const Eigen::Quaternion<ScalarType, Options>& q)
{
    return Eigen::Quaternion<ScalarType, Options>(a*q.w(),a*q.x(),a*q.y(),a*q.z());
}

template<typename ScalarType, int Options>
Eigen::Quaternion<ScalarType, Options> operator+(const Eigen::Quaternion<ScalarType, Options>& q1, const Eigen::Quaternion<ScalarType, Options>& q2)
{
    return Eigen::Quaternion<ScalarType, Options>(q1.w()+q2.w(),q1.x()+q2.x(),q1.y()+q2.y(),q1.z()+q2.z());
}
