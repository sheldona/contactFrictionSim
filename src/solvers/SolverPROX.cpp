#include "solvers/SolverPROX.h"

#include "contact/Contact.h"
#include "rigidbody/RigidBody.h"
#include "rigidbody/RigidBodySystem.h"


namespace 
{
    static inline bool checkInVector(std::vector<size_t> &v, size_t tag)
    {
        for (size_t i = 0; i < v.size(); i++)
        {
            if (v[i] == tag)
            {
                return true;
            }
        }
        return false;
    }
    
    static inline float infiniteNorm(const std::vector<Eigen::Vector3f> &v)
    {
        float max = -1;
        for (size_t i = 0; i < v.size(); i++)
        {
            for (int j = 0; j < v[i].size(); j++)
            {
                if (std::fabs(v[i][j]) > max)
                {
                    max = std::fabs(v[i][j]);
                }
            }
        }
        return max;
    }

    // lambda_f = prox_C (lambda_f - r (A * lambda_f + b))
    static inline void friction_solver(Eigen::Vector3f &z_k, Eigen::Vector3f &x_k, Eigen::Vector4f &coeffs)
    {
        // 0: n
        // 1: s
        // 2: t
        // Set all directions to the same mu for now

        using std::fabs;
        using std::max;
        using std::sqrt;

        const float z_s = z_k[1];
        const float z_t = z_k[2];

        // Set contraints to zero
        x_k[1] = 0;
        x_k[2] = 0;

        if (x_k[0] <= 0)
            return;

        // Compute a,b,c constants
        const float a = coeffs[0] * x_k[0];
        const float b = coeffs[1] * x_k[0];

        if (a <= 0 || b <= 0)
            return;

        const float tol = 10e-10f;
        const size_t max_k = 100;

        // Scale problem
        const float scale = max(1.0f, max(a, max(b, max(fabs(z_s), fabs(z_t)))));

        const float sx = z_s / scale;
        const float sy = z_t / scale;

        const float sa = a / scale;
        const float sb = b / scale;

        const float aa = sa * sa;
        const float bb = sb * sb;

        const float xx = sx * sx;
        const float yy = sy * sy;

        const float f0 = (xx / aa) + (yy / bb) - 1;

        if (f0 < tol)
        {
            x_k[1] = z_s;
            x_k[2] = z_t;
            return;
        }

        const float aaxx = aa * xx;
        const float bbyy = bb * yy;

        float t0 = 0.0f;
        float t1 = max(sa, sb) * sqrt(xx + yy);
        float g0 = (aaxx) / ((aa + t0) * (aa + t0)) + (bbyy) / ((bb + t0) * (bb + t0)) - 1.0f;
        float g1 = (aaxx) / ((aa + t1) * (aa + t1)) + (bbyy) / ((bb + t1) * (bb + t1)) - 1.0f;

        // Scalar expansion coefficient of the interval
        const float expansion = 1.5f;
        while (g1 > 0)
        {
            t1 *= expansion;
            g1 = (aaxx) / ((aa + t1) * (aa + t1)) + (bbyy) / ((bb + t1) * (bb + t1)) - 1.0f;
        }

        // Perform binary search
        float t_k = (t0 + t1) * 0.5f;
        for (size_t k = 0; k < max_k; k++)
        {
            if (fabs(t1 - t0) < tol)
                break;

            const float aat = aa + t_k;
            const float bbt = bb + t_k;
            const float g_k = aaxx / (aat * aat) + bbyy / (bbt * bbt) - 1.0f;
            if (fabs(g_k) < tol)
                break;
            if (g_k > 0)
            {
                t0 = t_k;
                g0 = g_k;
            }
            else
            {
                t1 = t_k;
                g1 = g_k;
            }
            t_k = (t0 + t1) * 0.5;
        }

        // Convert root to unscaled problem
        t_k *= scale * scale;

        // Compute closest point
        x_k[1] = (a * a * z_s) / (a * a + t_k);
        x_k[2] = (b * b * z_t) / (b * b + t_k);
    }

    // lambda_n = prox_{R^+} (lambda_n - r (A * lambda_n + b))
    static inline void normal_solver(Eigen::Vector3f &z_k, Eigen::Vector3f &x_k)
    {
        // Max strategy
        x_k[0] = std::max(0.0f, z_k[0]);
    }

    // z_k = x_k - R_kk * ( wbar + b_k )
    static inline void compute_z_k(Eigen::Vector3f &x_k, Eigen::VectorXf &w, const Eigen::MatrixXf &R_k, const Eigen::Vector3f &b, Eigen::Vector3f &z_k, size_t k, Eigen::MatrixXf &J)
    {
        // DEBUG: Check component values
        //    printVectorXf("wbar", wbar);
        //    printMatrixXf("R_k", R_k);
        //    printVector3f("x_k", x_k);
        //    printVectorXf("b", b);
        //    printVector3f("z_k0", z_k);

        z_k = Eigen::Vector3f::Zero();

        //    Eigen::VectorXf wbar = Eigen::VectorXf::Zero(6);
        //    wbar[0] = w0[0];
        //    wbar[1] = w0[1];
        //    wbar[2] = w0[2];

        //    wbar[3] = w1[0];
        //    wbar[4] = w1[1];
        //    wbar[5] = w1[2];

        //    z_k = J * wbar;
        z_k = J * w;

        // Add b_k
        z_k += b;

        // temp += R_kk * z_k
        Eigen::Vector3f temp = Eigen::Vector3f::Zero();
        temp += R_k * z_k;

        // Substract x_k - temp
        z_k = x_k - temp;
    }

    static inline void initialize_R(const std::vector<Eigen::Matrix3f> &A, std::vector<Eigen::MatrixXf> &R, std::vector<Eigen::MatrixXf> &nu)
    {
        for (size_t i = 0; i < R.size(); i++)
        {
            // Set to zero
            R[i] = Eigen::Matrix3f::Zero();
            nu[i] = Eigen::Matrix3f::Zero();

            for (int u = 0; u < R[i].rows(); u++)
            {
                //            R[i](u,u) = 0.1f;
                //            nu[i](u,u) = 0.5f;

                R[i](u, u) = 1.0f / A[i](u, u);
                nu[i](u, u) = 0.9f;
            }
        }
    }

    static inline void multAndSub(const Eigen::MatrixXf &G, const Eigen::Vector3f &x, const Eigen::Vector3f &y, const float a, Eigen::Vector3f &b)
    {
        b += a * G.col(0) * x(0);
        b += a * G.col(1) * x(1);
        b += a * G.col(2) * x(2);
        b += a * G.col(3) * y(0);
        b += a * G.col(4) * y(1);
        b += a * G.col(5) * y(2);
    }

    // Computes the right-hand side vector of the Schur complement system:
    //      b = gamma*phi/h - J*vel - dt*JMinvJT*force
    //
    static inline void buildRHS(Contact *c, float h, Eigen::Vector3f &b)
    {
        const float gamma = h * c->k / (h * c->k + c->b); // error reduction parameter
        b = -gamma * c->phi / h;

        multAndSub(c->J0, c->body0->xdot, c->body0->omega, 1.0f, b);
        multAndSub(c->J1, c->body1->xdot, c->body1->omega, 1.0f, b);

        if (!c->body0->fixed)
        {
            multAndSub(c->J0Minv, c->body0->f, c->body0->tau, h, b);
        }
        if (!c->body1->fixed)
        {
            multAndSub(c->J1Minv, c->body1->f, c->body1->tau, h, b);
        }
    }

    // Loop over all other contacts for a body and compute modifications to the rhs vector b:
    //           x -= (JMinv*Jother^T) * lambda_other
    //
    static inline void accumulateCoupledContacts(Contact *c, const JBlock &JMinv, RigidBody *body, Eigen::VectorXf &b)
    {
        if (body->fixed)
            return;

        for (Contact *cc : body->contacts)
        {
            if (cc != c)
            {
                if (body == cc->body0)
                    b -= JMinv * (cc->J0.transpose() * cc->lambda);
                else
                    b -= JMinv * (cc->J1.transpose() * cc->lambda);
            }
        }
    }
}

// Inputs: TBD
void SolverPROX::solveContact(std::vector<Eigen::Matrix3f> &A, std::vector<Eigen::Vector3f> &b, std::vector<Eigen::Vector3f> &lambda, std::vector<Eigen::Vector4f> &coeffs, std::vector<Eigen::MatrixXf> &J, std::vector<Eigen::MatrixXf> &MinvJT, std::vector<Contact *> &contacts)
{
    std::vector<Eigen::Vector3f> x;
    x.resize(A.size());

    x = lambda; // warm starting

    int abs_conv_in_iteration = -1;
    int rel_conv_in_iteration = -1;
    size_t count_divergence = 0;

    float relative_threshold = 0.001;
    float absolute_threshold = 0.001;

    std::vector<Eigen::Vector3f> residual;
    residual.resize(A.size());

    std::vector<Eigen::MatrixXf> R;
    std::vector<Eigen::MatrixXf> nu;

    R.resize(A.size());
    nu.resize(A.size());

    float residual_norm = 0;

    initialize_R(A, R, nu);

    Eigen::Vector3f z_k = Eigen::Vector3f::Zero();
    Eigen::Vector3f x_k = Eigen::Vector3f::Zero();
    Eigen::Vector3f deltax = Eigen::Vector3f::Zero();

    float last_residual_norm = FLT_MAX;

    // GS Loop
    for (size_t iter = 0; iter < m_maxIter; iter++)
    {
        // Loop over contact points
        for (size_t k = 0; k < A.size(); k++)
        {
            x_k = x[k];
            deltax = x[k];

            // Compute z_k
            compute_z_k(x_k, contacts[k]->w, R[k], b[k], z_k, k, J[k]);

            // Project the normal component
            normal_solver(z_k, x_k);

            x_k[1] = z_k[1];
            x_k[2] = z_k[2];

            // Project the tangential component
            friction_solver(z_k, x_k, coeffs[k]);

            // delta_x = x_k_new - x_k_old
            deltax = x_k - deltax;

            // Update wbar
            Eigen::VectorXf temp;
            temp = Eigen::VectorXf::Zero(6);
            temp += MinvJT[k] * deltax;

            std::vector<size_t> alreadyVisistedContacts;
            for (size_t i = 0; i < contacts[k]->body0->contacts.size(); i++)
            {
                //                std:: cout << "Loop of body0" << std::endl;
                Contact *cc = contacts[k]->body0->contacts[i];
                if (!checkInVector(alreadyVisistedContacts, cc->tag))
                {
                    alreadyVisistedContacts.push_back(cc->tag);

                    // Check if contact has the same bodies as contact_k
                    if ((contacts[k]->body0->tag == cc->body0->tag && contacts[k]->body1->tag == cc->body1->tag) || (contacts[k]->body1->tag == cc->body0->tag && contacts[k]->body0->tag == cc->body1->tag) || (contacts[k]->body0->tag == cc->body1->tag && contacts[k]->body1->tag == cc->body0->tag))
                    {
                        //                        std::cout << "Contact: " << cc->tag << std::endl;
                        cc->w += temp;
                    }
                }
            }

            residual[k] = lambda[k] - x_k;

            x[k] = x_k;
        }

        residual_norm = infiniteNorm(residual);

        // Check if residual is under absolute threshold
        if (residual_norm < absolute_threshold)
        {
            // std::cout << "Absolute convergence in " << iter << " iterations \t |residual| = " << residual_norm << std::endl;
            abs_conv_in_iteration = iter;
            // outf << "abs_conv: " << abs_conv_in_iteration << " rel_conv: " << -1 << "divergence count: " << count_divergence << '\n';
            break;
        }

        // Check if residual norm difference is under relative threshold
        if (std::fabs(residual_norm - last_residual_norm) < relative_threshold * last_residual_norm)
        {
            // std::cout << "Absolute convergence in " << iter << " iterations \t |residual| = " << residual_norm << std::endl;
            rel_conv_in_iteration = iter;
            // outf << "abs_conv: " << -1 << " rel_conv: " << rel_conv_in_iteration << "divergence count: " << count_divergence << '\n';
            break;
        }

        if (residual_norm > last_residual_norm)
        {
            for (size_t j = 0; j < R.size(); j++)
            {
                R[j] = nu[j] * R[j];
            }
            x = lambda;
            ++count_divergence;
        }
        else
        {
            last_residual_norm = residual_norm;
            lambda = x;
        }
    }
}

SolverPROX::SolverPROX(RigidBodySystem* _rigidBodySystem) : Solver(_rigidBodySystem)
{

}

void SolverPROX::solve(float h)
{
    std::vector<Contact *> &contacts = m_rigidBodySystem->getContacts();
    std::vector<RigidBody *> &bodies = m_rigidBodySystem->getBodies();
    const int numContacts = contacts.size();

    // // Clear w from each body
    for (size_t i = 0; i < bodies.size(); i++)
    {
        bodies[i]->w = Eigen::Vector3f::Zero();
        bodies[i]->tag = i;
    }

    // Build diagonal matrices of contacts
    std::vector<Eigen::Matrix3f> Acontacts;
    std::vector<Eigen::MatrixXf> AJ;
    std::vector<Eigen::MatrixXf> MinvJT;
    std::vector<Eigen::Vector3f> lambdasContacts;
    std::vector<Eigen::Vector4f> coeffs;
    std::vector<Eigen::Vector3f> bContacts;
    if (numContacts > 0)
    {
        // Build diagonal matrices
        Acontacts.resize(numContacts);
        lambdasContacts.resize(numContacts);
        coeffs.resize(numContacts);
        AJ.resize(numContacts);
        MinvJT.resize(numContacts);

        for (int i = 0; i < numContacts; ++i)
        {
            Contact *c = contacts[i];
            c->tag = i;
            const float eps = 1.0f / (h * h * c->k + h * c->b); // constraint force mixing

            c->w = Eigen::VectorXf::Zero(6);

            // Compute the diagonal term : Aii = J0*Minv0*J0^T + J1*Minv1*J1^T
            //
            Acontacts[i] = Eigen::Matrix3f::Zero();
            Acontacts[i](0, 0) += eps;

            AJ[i] = Eigen::MatrixXf::Zero(3, 6);
            MinvJT[i] = Eigen::MatrixXf::Zero(6, 3);

            if (!c->body0->fixed)
            {
                JBlockTranspose JT = c->J0.transpose();
                Acontacts[i] += c->J0Minv * JT;
                AJ[i] += c->J0;
                MinvJT[i] += c->J0Minv.transpose();
            }
            if (!c->body1->fixed)
            {
                JBlockTranspose JT = c->J1.transpose();
                Acontacts[i] += c->J1Minv * JT;
                AJ[i] += c->J1;
                MinvJT[i] += c->J1Minv.transpose();
            }

            lambdasContacts[i] = c->lambda;
            coeffs[i][0] = c->mu;
            coeffs[i][1] = c->mu;
            coeffs[i][2] = c->mu;
            coeffs[i][3] = c->mu;
        }

        // Assemble b
        bContacts.resize(numContacts);
        for (int i = 0; i < numContacts; ++i)
        {
            Contact *c = contacts[i];
            buildRHS(c, h, bContacts[i]);
            c->lambda.setZero();
        }
    }

    if (numContacts > 0)
    {
        //        // DEBUGGING: Check size of lambdas before
        //        for (size_t i=0; i<lambdasContacts.size(); i++) {
        //            std::cout << i << " lambda: " << "(" << lambdasContacts[i][0] << ", " << lambdasContacts[i][1] << ", " << lambdasContacts[i][2] << ")" << std::endl;
        //        }

        solveContact(Acontacts, bContacts, lambdasContacts, coeffs, AJ, MinvJT, contacts);

        // DEBUGGING: Check size of lambdas after
        //        for (size_t i=0; i<lambdasContacts.size(); i++) {
        //            std::cout << i << " lambda: " << "(" << lambdasContacts[i][0] << ", " << lambdasContacts[i][1] << ", " << lambdasContacts[i][2] << ")" << std::endl;
        //        }

        // Set the new lambdas to the old contact lambdas
        for (size_t i = 0; i < lambdasContacts.size(); i++)
        {
            contacts[i]->lambda = lambdasContacts[i];
        }
        //        std::cout << "=============================================================================" << std::endl;
    }
}