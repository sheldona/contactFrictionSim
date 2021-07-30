#include "solvers/SolverBoxPGS.h"

#include "contact/Contact.h"
#include "rigidbody/RigidBody.h"
#include "rigidbody/RigidBodySystem.h"

#include <Eigen/Dense>


namespace
{
    static inline void multAndSub(const JBlock& G, const Eigen::Vector3f& x, const Eigen::Vector3f& y, const float a, Eigen::VectorXf& b)
    {
            b -= a * G.col(0) * x(0);
            b -= a * G.col(1) * x(1);
            b -= a * G.col(2) * x(2);
            b -= a * G.col(3) * y(0);
            b -= a * G.col(4) * y(1);
            b -= a * G.col(5) * y(2);
    }

    static inline void buildRHS(Contact* c, float h, Eigen::VectorXf& b)
    {
        const float hinv = 1.0f / h;
        const float gamma = h * c->k / (h * c->k + c->b); // error reduction parameter
        b = hinv * gamma * c->phi;

        multAndSub(c->J0, c->body0->xdot, c->body0->omega, 1.0f, b);
        multAndSub(c->J1, c->body1->xdot, c->body1->omega, 1.0f, b);

        if( !c->body0->fixed )
        {
            multAndSub(c->J0Minv, c->body0->f, c->body0->tau, h, b);
        }
        if( !c->body1->fixed )
        {
            multAndSub(c->J1Minv, c->body1->f, c->body1->tau, h, b);
        }
    }

    // Loop over all other contacts involving _body
    //      and compute : x -= (J0*Minv0*Jother^T) * lambda_other
    //
    static inline void accumulateCoupledContacts(Contact* c, const JBlock& JMinv, RigidBody* body, Eigen::VectorXf& x)
    {
        if( body->fixed )
            return;

        for(Contact* cc : body->contacts)
        {
            if( cc != c )
            {
                if( body == cc->body0 )
                {
                    x -= JMinv * (cc->J0.transpose() * cc->lambda);
                }
                else
                {
                    x -= JMinv * (cc->J1.transpose() * cc->lambda);
                }
            }
        }
    }

    static inline void solveContact(const Eigen::Matrix3f& A, const Eigen::VectorXf& b, Eigen::VectorXf& x, const float mu)
    {
        x(0) = std::max(0.0f, (b(0) - A(0,1) * x(1) - A(0,2) * x(2) ) / A(0,0) );
        x(1) = std::max(-mu*x(0), std::min(mu*x(0), ( b(1) - A(1,0) * x(0) - A(1,2) * x(2) ) / A(1,1) ));
        x(2) = std::max(-mu*x(0), std::min(mu*x(0), ( b(2) - A(2,0) * x(0) - A(2,1) * x(1) ) / A(2,2) ));
    }
}

SolverBoxPGS::SolverBoxPGS(RigidBodySystem* _rigidBodySystem) : Solver(_rigidBodySystem)
{

}

void SolverBoxPGS::solve(float h)
{
    std::vector<Contact*>& contacts = m_rigidBodySystem->getContacts();
    const int numContacts = contacts.size();

    // Build diagonal matrices of contacts
    std::vector<Eigen::MatrixXf> Acontactii;
    if( numContacts > 0 )
    {
        // Build diagonal matrices
        Acontactii.resize(numContacts);
        for(int i = 0; i < numContacts; ++i)
        {
            Contact* c = contacts[i];
            const float eps = 1.0f / (h * h * c->k + h * c->b);    // constraint force mixing

            // Compute the diagonal term : Aii = J0*Minv0*J0^T + J1*Minv1*J1^T
            //
            Acontactii[i].setZero(3,3);
            Acontactii[i](0,0) += eps;

            if( !c->body0->fixed )
            {
                Acontactii[i] += c->J0Minv * c->J0.transpose();
            }
            if( !c->body1->fixed )
            {
                Acontactii[i] += c->J1Minv * c->J1.transpose();
            }
        }

        std::vector<Eigen::VectorXf> b;
        b.resize(numContacts);

        // Compute the right-hand side vector : b = phi - J*vel - dt*JMinvJT*force
        //
        for(int i = 0; i < numContacts; ++i)
        {
            Contact* c = contacts[i];
            buildRHS(c, h, b[i]);
            c->lambda.setZero();
        }

        // PGS main loop.
        // There is no convergence test here.
        // Stop after @a maxIter iterations.
        //
        for(int iter = 0; iter < m_maxIter; ++iter)
        {
            // For each contact, compute an updated value of contacts[i]->lambda
            //      using matrix-free pseudo-code provided in the course notes.
            //
            for(int i = 0; i < numContacts; ++i)
            {
                Contact* c = contacts[i];

                // Initialize current solution as x = b[i]
                Eigen::VectorXf x = b[i];

                accumulateCoupledContacts(c, c->J0Minv, c->body0, x);
                accumulateCoupledContacts(c, c->J1Minv, c->body1, x);
                solveContact(Acontactii[i], x, c->lambda, c->mu);
            }
        }
    }
}

