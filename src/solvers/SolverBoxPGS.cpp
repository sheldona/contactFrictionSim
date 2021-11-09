#include "solvers/SolverBoxPGS.h"

#include "contact/Contact.h"
#include "rigidbody/RigidBody.h"
#include "rigidbody/RigidBodySystem.h"

#include <Eigen/Dense>

std::ofstream pgs_merit_file{ "pgs_merit.csv" };
std::ofstream pgs_conv_file{ "convergence_pgs.csv" };
std::ofstream pgs_lcp_energy_error_file{ "pgs_lcp_energy_error.csv" };
std::set<int> pgs_timestep_samples{ 1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 100, 120, 140, 160, 180 };
bool PGS_ALL_SAMPLES = true;


namespace
{
    static inline float infiniteNorm(const std::vector<Eigen::VectorXf>& v) {
        float max = -1;
        for (size_t i = 0; i < v.size(); i++) {
            for (int j = 0; j < v[i].size(); j++) {
                if (std::fabs(v[i][j]) > max) {
                    max = std::fabs(v[i][j]);
                }
            }
        }
        return max;
    }

    // LCP errors that just checks and make sure the error of the friction components matches its bounds and error in the free component is close to 0
    float get_complementarity_error(Contact* c, Eigen::Matrix3f Aii, Eigen::VectorXf b)
    {
        static const float inf_tol = FLT_MAX;

        float err = 0.0;
        float d;
        for (int i = 0; i < Aii.rows(); ++i)
        {
            float lower = -c->mu * c->lambda(i);
            float upper = c->mu * c->lambda(i);
            float w = Aii(i, i) * c->lambda(i) + b(i);

            if (lower <= -inf_tol)
            {
                if (upper >= inf_tol)       // free
                {
                    d = w;
                }
                else                        // bounded above
                {
                    d = std::min(Aii(i,i) * (upper - c->lambda(i)), -w);
                }
            }
            else
            {
                if (upper >= inf_tol)       // bounded below
                {
                    d = std::min(Aii(i,i) * (c->lambda(i) - lower), w);
                }
                else                        // box bounds
                {
                    const float wplus = std::max(w, 0.0f);
                    const float wminus = -std::min(w, 0.0f);
                    d = std::fabs(std::min(Aii(i,i) * (c->lambda(i) - lower), wplus)) + std::fabs(std::min(Aii(i,i) * (upper - c->lambda(i)), wminus));
                }
            }

            err += d * d;
        }
        return std::sqrt(err);
    }

    // LCP error that is in units of Joules
    float unit_consistent_energy_lcp_error(Contact* c, Eigen::Matrix3f Aii, Eigen::VectorXf b)
    {
        float err = 0.0f;
        for (int i = 0; i < c->lambda.size(); i++) {
            float lower = -c->mu * c->lambda(i);
            float upper = c->mu * c->lambda(i);
            float w = Aii(i, i) * c->lambda(i) + b(i);

            // Compute violation of the upperand lower bound
            float delta_x_u = std::max(c->lambda(i) - upper, 0.0f);
            float delta_x_l = std::max(lower - c->lambda(i), 0.0f);

            // Compute feasible component x0
            float x0 = std::max(lower, std::min(c->lambda(i), upper));

            // Compute upper impulse energy error
            float delta_e_x_u = 0.5f * Aii(i, i) * (delta_x_u * delta_x_u);

            // Compute lower impulse energy error
            float delta_e_x_l = 0.5f * Aii(i, i) * (delta_x_l * delta_x_l);

            // Compute w_+ and w_-
            float w_plus = std::max(w, 0.0f);
            float w_minus = std::max(w, 0.0f);

            // Compute the saturation of the upperand lower bound
            float sigma_l = (x0 + delta_x_u) - lower;
            float sigma_u = upper - (x0 - delta_x_l);

            // Compute the positive velocity energy error
            float delta_e_w_plus = std::min(((1.0f / (2.0f * (Aii(i,i) + 1e-10f))) * (w_plus * w_plus)), (0.5f * Aii(i,i) * (sigma_l * sigma_l)));
        
            // Compute the negative velocity energy error
            float delta_e_w_minus = std::min(((1.0f / (2.0f * (Aii(i,i) + 1e-10f))) * (w_minus * w_minus)), (0.5f * Aii(i,i) * (sigma_u * sigma_u)));

            // Compute the energy error per constraint
            float delta_e = std::max(std::max(std::abs(delta_e_x_u), std::abs(delta_e_x_l)), std::max(std::abs(delta_e_w_plus), std::abs(delta_e_w_minus)));

            err += std::abs(delta_e);
        }

        return err;
    }


    static inline void multAndSub(const JBlock& G, const Eigen::Vector3f& x, const Eigen::Vector3f& y, const float a, Eigen::VectorXf& b)
    {
            b -= a * G.col(0) * x(0);
            b -= a * G.col(1) * x(1);
            b -= a * G.col(2) * x(2);
            b -= a * G.col(3) * y(0);
            b -= a * G.col(4) * y(1);
            b -= a * G.col(5) * y(2);
    }


	// Computes the right-hand side vector of the Schur complement system: 
    //      b = gamma*phi/h - J*vel - dt*JMinvJT*force
    //
    static inline void buildRHS(Contact* c, float h, Eigen::VectorXf& b)
    {
        const float gamma = h * c->k / (h * c->k + c->b);       // error reduction parameter
        b = -gamma * c->phi / h;

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

    // Loop over all other contacts for a body and compute modifications to the rhs vector b: 
    //           x -= (JMinv*Jother^T) * lambda_other
    //
    static inline void accumulateCoupledContacts(Contact* c, const JBlock& JMinv, RigidBody* body, Eigen::VectorXf& b)
    {
        if( body->fixed )
            return;

        for(Contact* cc : body->contacts)
        {
            if( cc != c )
            {
                if( body == cc->body0 )
                    b -= JMinv * (cc->J0.transpose() * cc->lambda);
                else
                    b -= JMinv * (cc->J1.transpose() * cc->lambda);
            }
        }
    }

    // Solve the Boxed LCP problem for a single contact and isotropic Coulomb friction.
    // The solution vector, @a x, contains the impulse the non-interpenetration constraint in x(0), and
    // the friction constraints in x(1) and x(2)
    // 
    // The solution is projected to the lower and upper bounds imposed by the box model.
    // 
    // Inputs: 
    //    x - contains three impulse variables (non-interpenetration + two friction)
    //    b - rhs vector
    //    mu - the friction coefficient
    static inline void solveContact(const Eigen::Matrix3f& A, const Eigen::VectorXf& b, Eigen::VectorXf& x, const float mu)
    {
        // Normal impulse is projected to [0, inf]
        //
        x(0) = std::max(0.0f, (b(0) - A(0,1) * x(1) - A(0,2) * x(2) ) / A(0,0) );

        // Next, friction impulses are projected to [-mu * x(0), mu * x(1)]
        //
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

    // Build array of 3x3 diagonal matrices, one for each contact.
    // 
    std::vector<Eigen::Matrix3f> Acontactii;
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

        // Compute the right-hand side vector : 
        //      b = -gamma*phi/h - J*vel - dt*JMinvJT*force
        //
        for(int i = 0; i < numContacts; ++i)
        {
            Contact* c = contacts[i];
            buildRHS(c, h, b[i]);
            c->lambda.setZero();
        }

        // Setup residual calculations
        std::vector<Eigen::VectorXf> x_pre;
        x_pre.resize(numContacts);
        for (size_t i = 0; i < x_pre.size(); i++) {
            x_pre[i] = contacts[i]->lambda;
        }
        
        std::vector<Eigen::VectorXf> residual;
        residual.resize(numContacts);
        for (size_t i = 0; i < residual.size(); i++) {
            residual[i] = Eigen::Vector3f::Zero();
        }

        float residual_norm = 0;
        float relative_threshold = 0.0001;
        float absolute_threshold = 0.0001;

        int abs_conv_in_iteration = -1;
        int rel_conv_in_iteration = -1;
        size_t count_divergence = 0;

        float last_residual_norm = FLT_MAX;

        const bool is_in = pgs_timestep_samples.find(m_rigidBodySystem->getFrameCount()) != pgs_timestep_samples.end();
        if (is_in || PGS_ALL_SAMPLES)
            pgs_lcp_energy_error_file << "#" << m_rigidBodySystem->getFrameCount() << "\n";
            pgs_merit_file << "#" << m_rigidBodySystem->getFrameCount() << "\n";

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

                residual[i] = x_pre[i] - c->lambda;
                if (is_in || PGS_ALL_SAMPLES)
                    pgs_lcp_energy_error_file << unit_consistent_energy_lcp_error(contacts[i], Acontactii[i], b[i]) << " ";
            }

            residual_norm = infiniteNorm(residual);
            if (is_in || PGS_ALL_SAMPLES)
                pgs_merit_file << residual_norm << "\n";

            // Check if residual is under abs threshold
            if (residual_norm < absolute_threshold) {
                abs_conv_in_iteration = iter;
                break;
            }

            // Check if residual norm diff is under relative threshold
            if (std::fabs(residual_norm - last_residual_norm) < relative_threshold * last_residual_norm) {
                rel_conv_in_iteration = iter;
                break;
            }

            // Check for divergence
            if (residual_norm > last_residual_norm) {
                ++count_divergence;
            }
            last_residual_norm = residual_norm;
            for (int k = 0; k < numContacts; k++) {
                x_pre[k] = contacts[k]->lambda;
            }
            if (is_in || PGS_ALL_SAMPLES)
                pgs_lcp_energy_error_file << '\n';
        }

        pgs_conv_file << abs_conv_in_iteration << "," << rel_conv_in_iteration << "," << count_divergence << "," << residual_norm << '\n';

        if (is_in || PGS_ALL_SAMPLES)
            pgs_lcp_energy_error_file << '\n';
    }
}

