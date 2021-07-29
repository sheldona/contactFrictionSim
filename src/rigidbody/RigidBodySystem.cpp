#include "rigidbody/RigidBodySystem.h"

#include "rigidbody/RigidBody.h"

#include "collision/CollisionDetect.h"
#include "contact/Contact.h"
#include "util/MeshAssets.h"

#include "solvers/SolverBoxPGS.h"
#include "solvers/SolverBoxBPP.h"

namespace Eigen
{
    typedef Matrix<float, 6, 1, ColMajor> Vector6f;
}


RigidBodySystem::RigidBodySystem() :
    m_frame(0),
    m_contactStiffness(1e5f), m_contactDamping(1e4f),
    m_mu(0.4f), m_pgsIter(10), m_preStepFunc(nullptr), m_resetFunc(nullptr)
{
    m_collisionDetect = std::make_unique<CollisionDetect>(this);
    m_solver = new SolverBoxPGS(this);
    //m_solver = new SolverBoxBPP(this);
}

RigidBodySystem::~RigidBodySystem()
{
    clear();
    delete m_solver;
}

void RigidBodySystem::addBody(RigidBody *_b)
{
    _b->R = _b->q.toRotationMatrix();
    m_bodies.push_back(_b);
}

void RigidBodySystem::step(float dt)
{
    // Initialize the system.
    // Apply gravitional forces and reset angular forces.
    // Cleanup contacts from the previous time step.
    //
    for(auto b : m_bodies)
    {
        b->f = b->mass * Eigen::Vector3f(0.f, -9.81f, 0.f);
        b->tau.setZero();
        b->fc.setZero();
        b->tauc.setZero();
        b->contacts.clear();
    }

    // Standard simulation pipeline.
    computeInertias();

    if( m_preStepFunc )
    {
        m_preStepFunc(m_bodies);
    }

    m_collisionDetect->clear();
    m_collisionDetect->detectCollisions();
    m_collisionDetect->computeContactJacobians();

    // Update contact stiffness and damping
    //
    auto contacts = m_collisionDetect->getContacts();
    for(auto c : contacts)
    {
        c->k = m_contactStiffness;
        c->b = m_contactDamping;
        c->mu = m_mu;
    }

    for(auto b : m_bodies)
    {
        b->fc.setZero();
        b->tauc.setZero();
    }

    calcConstraintForces(dt);

    // Perform numerical integration of each rigid body in @a m_bodies
    // using the equations provided in the course notes and
    // the Baraff 2001 SIGGRAPH course.
    //
    for(RigidBody* b : m_bodies)
    {
        if( !b->fixed )
        {
            b->xdot += dt * (1.0f/b->mass) * (b->f + b->fc);
            b->omega += dt * b->Iinv * (b->tau + b->tauc);
        }
        else
        {
            b->xdot.setZero();
            b->omega.setZero();
        }

        b->x += dt * b->xdot;
        b->q = b->q + 0.5f * dt * Eigen::Quaternionf(0, b->omega[0], b->omega[1], b->omega[2])*b->q;
        b->q.normalize();
    }

    ++m_frame;
}

void RigidBodySystem::clear()
{
    if( m_resetFunc )
    {
        m_resetFunc();
    }

    m_collisionDetect->clear();

    for(auto b : m_bodies)
    {
        delete b;
    }
    m_bodies.clear();

    m_frame = 0;
}



void RigidBodySystem::computeInertias()
{
    for(RigidBody* b : m_bodies)
    {
        b->R = b->q.toRotationMatrix();
        b->updateInertiaMatrix();
    }
}

// Accessors for the contact constraint array.
const std::vector<Contact*>& RigidBodySystem::getContacts() const
{
    return m_collisionDetect->getContacts();
}

std::vector<Contact*>& RigidBodySystem::getContacts()
{
    return m_collisionDetect->getContacts();
}

void RigidBodySystem::calcConstraintForces(float dt)
{
    // Solve for the constraint forces lambda
    //
    //m_solver->setMaxIter(m_pgsIter);
    m_solver->solve(dt);

    // Apply the constraint forces as forces and torques acting on each body.
    // Essentially, we compute contact forces by computing :
    //
    //       f_contact = J^T * lambda
    //
    // for each contact constraint.
    //
    auto contacts = m_collisionDetect->getContacts();
    for(const auto c : contacts)
    {
        const Eigen::Vector6f f0 = c->J0.transpose() * c->lambda / dt;
        const Eigen::Vector6f f1 = c->J1.transpose() * c->lambda / dt;
        c->body0->fc += f0.head<3>();
        c->body0->tauc += f0.tail<3>();
        c->body1->fc += f1.head<3>();
        c->body1->tauc += f1.tail<3>();
    }
}
