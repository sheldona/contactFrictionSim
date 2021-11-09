#include "solvers/SolverPROX.h"

#include "contact/Contact.h"
#include "rigidbody/RigidBody.h"
#include "rigidbody/RigidBodySystem.h"

namespace  {

static inline void multAndAdd(const Eigen::MatrixXf& G, const Eigen::Vector3f& x, const Eigen::Vector3f& y, const float a, Eigen::Vector3f& b)
{
    b += a * G.col(0) * x(0);
    b += a * G.col(1) * x(1);
    b += a * G.col(2) * x(2);
    b += a * G.col(3) * y(0);
    b += a * G.col(4) * y(1);
    b += a * G.col(5) * y(2);
}

// Computes the right-hand side vector of the Schur complement system: 
static inline void buildRHS(Contact* c, float h, Eigen::Vector3f& b)
{
    const float gamma = h * c->k / (h * c->k + c->b);       // error reduction parameter
    b = gamma * c->phi / h;

    multAndAdd(c->J0, c->body0->xdot, c->body0->omega, 1.0f, b);
    multAndAdd(c->J1, c->body1->xdot, c->body1->omega, 1.0f, b);

    if (!c->body0->fixed)
    {
        multAndAdd(c->J0Minv, c->body0->f, c->body0->tau, h, b);
    }
    if (!c->body1->fixed)
    {
        multAndAdd(c->J1Minv, c->body1->f, c->body1->tau, h, b);
    }
}

static inline float infiniteNorm(const std::vector<Eigen::Vector3f>& v) {
    float max = -1;
    for (size_t i=0; i<v.size(); i++) {
        for (int j=0; j<v[i].size(); j++) {
            if (std::fabs(v[i][j]) > max) {
                max = std::fabs(v[i][j]);
            }
        }
    }
    return max;
}

// z_k = x_k - R_kk * ( wbar + b_k )
static inline void computeZk(Eigen::Vector3f& x_k, Eigen::VectorXf& w0, Eigen::VectorXf& w1, const Eigen::MatrixXf& R_k, 
    const Eigen::Vector3f& b, Eigen::Vector3f& z_k, JBlock& J0, JBlock& J1) {

    z_k = Eigen::Vector3f::Zero();

    z_k = J0 * w0 + J1 * w1;

    // Add b_k
    z_k += b;

    Eigen::Vector3f temp = Eigen::Vector3f::Zero();
    temp = R_k * z_k;

    // Substract x_k - temp
    z_k = x_k - temp;

}

static inline void initializeR(const std::vector<Eigen::Matrix3f>& A, std::vector<Eigen::MatrixXf>& R, std::vector<Eigen::MatrixXf>& nu) {
    for (size_t i=0; i<R.size(); i++) {
        // Set to zero
        R[i] = Eigen::Matrix3f::Zero();
        nu[i] = Eigen::Matrix3f::Zero();

        for (int u=0; u<R[i].rows(); u++) {
            R[i](u,u) = 1.0f / A[i](u,u);
            nu[i](u,u) = 0.9f;
        }
    }
}

// lambda_n = prox_{R^+} (lambda_n - r (A * lambda_n + b))
static inline void normalSolver(Eigen::Vector3f& z_k, Eigen::Vector3f& x_k) {
    // Max strategy
    x_k[0] = std::max(0.0f, z_k[0]);
    x_k[1] = z_k[1];
    x_k[2] = z_k[2];
}

// lambda_f = prox_C (lambda_f - r (A * lambda_f + b))
static inline void frictionSolver(Eigen::Vector3f& z_k, Eigen::Vector3f& x_k, Eigen::Vector4f& coeffs) {
    // 0: n
    // 1: s
    // 2: t
    // Set all directions to the same mu for now

    const float z_s = z_k[1];
    const float z_t = z_k[2];

    // Set contraints to zero
    x_k[1] = 0; x_k[2] = 0;

    if (x_k[0] <= 0) return;

    // Compute a,b,c constants
    const float a = coeffs[0] * x_k[0];
    const float b = coeffs[1] * x_k[0];

    if (a <= 0 || b <= 0) return;

    const float tol = 10e-10f;
    const size_t max_k = 100;

    // Scale problem
    const float scale = std::max(1.0f, std::max( a, std::max( b, std::max(std::fabs( z_s ), std::fabs( z_t )) ) ) );

    const float sx = z_s / scale;
    const float sy = z_t / scale;

    const float sa = a / scale;
    const float sb = b / scale;

    const float aa = sa*sa;
    const float bb = sb*sb;

    const float xx = sx*sx;
    const float yy = sy*sy;

    const float f0 = (xx/aa) + (yy/bb) - 1;

    if (f0 < tol) {
        x_k[1] = z_s;
        x_k[2] = z_t;
        return;
    }

    const float aaxx = aa*xx;
    const float bbyy = bb*yy;

    float t0 = 0.0f;
    float t1 = std::max(sa, sb) * std::sqrt(xx + yy);
    float g0 = (aaxx) / ((aa+t0)*(aa+t0)) + (bbyy) / ((bb+t0)*(bb+t0)) - 1.0f;
    float g1 = (aaxx) / ((aa+t1)*(aa+t1)) + (bbyy) / ((bb+t1)*(bb+t1)) - 1.0f;

    // Scalar expansion coefficient of the interval
    const float expansion = 1.5f;
    while(g1 > 0) {
        t1 *= expansion;
        g1 = (aaxx) / ((aa+t1)*(aa+t1)) + (bbyy) / ((bb+t1)*(bb+t1)) - 1.0f;
    }

    // Perform binary search
    float t_k = (t0+t1) * 0.5f;
    for (size_t k=0; k<max_k; k++) {
        if (std::fabs(t1-t0) < tol) break;

        const float aat = aa + t_k;
        const float bbt = bb + t_k;
        const float g_k = aaxx/(aat*aat) + bbyy/(bbt*bbt) - 1.0f;
        if (std::fabs(g_k) < tol) break;
        if (g_k > 0) {
            t0 = t_k; g0 = g_k;
        } else {
            t1 = t_k; g1 = g_k;
        }
        t_k = (t0+t1) * 0.5;
    }

    // Convert root to unscaled problem
    t_k *= scale*scale;

    // Compute closest point
    x_k[1] = (a*a*z_s) / (a*a + t_k);
    x_k[2] = (b*b*z_t) / (b*b + t_k);
}

static inline void updateW(Eigen::VectorXf& w, const Eigen::Matrix3f& MinvJT, const Eigen::Matrix3f& IinvJT, const Eigen::Vector3f& lambda, bool fixedFlag)
{
    if (!fixedFlag) {
        w.head(3) += MinvJT * lambda;
        w.tail(3) += IinvJT * lambda;
    }
}
}


SolverPROX::SolverPROX(RigidBodySystem* _rigidBodySystem) : Solver(_rigidBodySystem)
{

}

void SolverPROX::solve(float h)
{
    std::vector<Contact*>& contacts = m_rigidBodySystem->getContacts();
    std::vector<RigidBody*>& bodies = m_rigidBodySystem->getBodies();
    const int numBodies = bodies.size();
    const int numContacts = contacts.size();

    // Give an identifier to each body to be able to index the w matrix
    for (int i = 0; i < numBodies; i++) {
        bodies[i]->index = i;
    }

    // Build diagonal matrices of contacts
    std::vector<Eigen::Matrix3f> Acontacts;
    std::vector<Eigen::Vector3f> lambdasContacts;
    std::vector<Eigen::Vector4f> coeffs;
    std::vector<Eigen::Vector3f> bContacts;
    std::vector<Eigen::Matrix3f> MinvJ0T;
    std::vector<Eigen::Matrix3f> MinvJ1T;
    std::vector<Eigen::Matrix3f> IinvJ0T;
    std::vector<Eigen::Matrix3f> IinvJ1T;

    if( numContacts > 0)
    {
        // Build diagonal matrices
        Acontacts.resize(numContacts);
        lambdasContacts.resize(numContacts);
        coeffs.resize(numContacts);
        MinvJ0T.resize(numContacts);
        MinvJ1T.resize(numContacts);
        IinvJ0T.resize(numContacts);
        IinvJ1T.resize(numContacts);

        // Build components for Minv * JT
        for (int i = 0; i < numContacts; i++) {
            MinvJ0T[i] = (1.0f / contacts[i]->body0->mass) * contacts[i]->J0.transpose().block(0, 0, 3, 3);
            IinvJ0T[i] = contacts[i]->body0->Iinv * contacts[i]->J0.transpose().block(3, 0, 3, 3);

            MinvJ1T[i] = (1.0f / contacts[i]->body1->mass) * contacts[i]->J1.transpose().block(0, 0, 3, 3);
            IinvJ1T[i] = contacts[i]->body1->Iinv * contacts[i]->J1.transpose().block(3, 0, 3, 3);
        }

        for(int i = 0; i < numContacts; ++i)
        {
            Contact* c = contacts[i];
            c->tag = i;
            contacts[i]->w = Eigen::Vector3f::Zero();

            c->w = Eigen::VectorXf::Zero(6);

            // Compute the diagonal term : Aii = J0*Minv0*J0^T + J1*Minv1*J1^T
            //
            Acontacts[i] = Eigen::Matrix3f::Zero();
            Acontacts[i](0,0) += 1.0f / (h * h * c->k + h * c->b); // constraint force mixing

            if( !c->body0->fixed )
            {
                JBlockTranspose JT = c->J0.transpose();
                Acontacts[i] += c->J0Minv * JT;
            }
            if( !c->body1->fixed )
            {
                JBlockTranspose JT = c->J1.transpose();
                Acontacts[i] += c->J1Minv * JT;
            }

            lambdasContacts[i] = c->lambda;

            // Isotropic friction
            coeffs[i][0] = c->mu;
            coeffs[i][1] = c->mu;
            coeffs[i][2] = c->mu;
            coeffs[i][3] = c->mu;
        }

        // Assemble b
        bContacts.resize(numContacts);
        for (int i=0; i<numContacts; ++i) {
            Contact* c = contacts[i];
            bContacts[i] = Eigen::Vector3f::Zero();
            buildRHS(c, h, bContacts[i]);
            c->lambda.setZero();
        }
    }

    // PROX solver
    if (numContacts > 0) {
        float residual_norm = 0;
        float last_residual_norm = FLT_MAX;
        Eigen::Vector3f z_k = Eigen::Vector3f::Zero();

        float relative_threshold = 0.0001;
        float absolute_threshold = 0.0001;

        std::vector<Eigen::Vector3f> lambdaCandidate;
        std::vector<Eigen::MatrixXf> R;
        std::vector<Eigen::MatrixXf> nu;
        std::vector<Eigen::Vector3f> residual;

        lambdaCandidate.resize(Acontacts.size());
        residual.resize(Acontacts.size());
        R.resize(Acontacts.size());
        nu.resize(Acontacts.size());
        w.resize(numBodies);

        initializeR(Acontacts, R, nu);

        for (int iter = 0; iter < m_maxIter; ++iter) {
            for (int i = 0; i < numBodies; i++) {
                w[i] = Eigen::VectorXf::Zero(6);
            }

            // Initialize w = Minv * JT * lambda_k
            for (int i = 0; i < numContacts; ++i) {
                lambdaCandidate[i] = Eigen::Vector3f::Zero();
                updateW(w[contacts[i]->body0->index], MinvJ0T[i], IinvJ0T[i], lambdasContacts[i], contacts[i]->body0->fixed);
                updateW(w[contacts[i]->body1->index], MinvJ1T[i], IinvJ1T[i], lambdasContacts[i], contacts[i]->body1->fixed);
            }

            // Solve for each contact
            for (int i = 0; i < numContacts; ++i) {
                computeZk(lambdasContacts[i], w[contacts[i]->body0->index], w[contacts[i]->body1->index], R[i], bContacts[i], z_k, 
                    contacts[i]->J0, contacts[i]->J1);
                
                normalSolver(z_k, lambdaCandidate[i]);

                frictionSolver(z_k, lambdaCandidate[i], coeffs[i]);

                residual[i] = lambdaCandidate[i] - lambdasContacts[i];
                
                // Update w = Minv * JT * (lambda(k+1) - lambda(k))
                updateW(w[contacts[i]->body0->index], MinvJ0T[i], IinvJ0T[i], lambdaCandidate[i] - lambdasContacts[i], contacts[i]->body0->fixed);
                updateW(w[contacts[i]->body1->index], MinvJ1T[i], IinvJ1T[i], lambdaCandidate[i] - lambdasContacts[i], contacts[i]->body1->fixed);
            }

            residual_norm = infiniteNorm(residual);
            
            // Check for convergence for early termination
            if (residual_norm < absolute_threshold)
                break;
            if (std::fabs(residual_norm - last_residual_norm) < relative_threshold * last_residual_norm)
                break;

            // Only update lambda if solver is converging
            if (residual_norm > last_residual_norm) {
                for (size_t j = 0; j < R.size(); j++) {
                    R[j] = nu[j] * R[j];
                }
            }
            else {
                last_residual_norm = residual_norm;
                for (int j = 0; j < numContacts; j++) {
                    lambdasContacts[j] = lambdaCandidate[j];
                }
            }
        }

        // Set the new lambdas to the old contact lambdas
        for (size_t i=0; i<lambdasContacts.size(); i++) {
            contacts[i]->lambda = lambdasContacts[i];
        }
    }
}