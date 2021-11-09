#include "solvers/SolverPROX.h"

#include "contact/Contact.h"
#include "rigidbody/RigidBody.h"
#include "rigidbody/RigidBodySystem.h"
#include <set>

#include <fstream>


// Create file to gather data for LCP error
// Currently the format is: 
// #<timestep>
// per line lcp error 
std::ofstream prox_merit_file{ "prox_merit.csv" };
std::ofstream prox_lcp_energy_error_file{ "prox_lcp_energy_error.csv" };
std::ofstream conv_file{"convergence_prox.csv"};
std::set<int> timestep_samples{ 1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 100, 120, 140, 160, 180 };
bool ALL_SAMPLES = true;

namespace  {

static inline void print(std::string name, const Eigen::Vector3f& v) {
    std::cout << name << " : (" << std::scientific << v[0] << ", " << std::scientific << v[1] << ", " << std::scientific << v[2] << ")" << std::endl;
}

static inline void print(std::string name, const Eigen::VectorXf &v)
{
    std::cout << name << " : (";
    for (int i=0; i<v.size(); i++) {
        std::cout << std::scientific << v[i] << ", ";
    }
    std::cout << ")" << std::endl;
}

static inline void print(std::string name, const Eigen::Matrix3f &m)
{
    std::cout << name << " : ";
    for (int i=0; i<m.rows(); i++) {
        std::cout << std::scientific << m(i, 0) << ", " << std::scientific << m(i, 1) << ", " << std::scientific << m(i, 2) << std::endl;
    }
}

static inline void print(std::string name, const Eigen::MatrixXf& m) {
    std::cout << name << " : " << std::endl;
    for (int i=0; i<m.rows(); i++) {
        for (int j=0; j<m.cols(); j++) {
            std::cout << std::scientific << m(i, j) << ", ";
        }
        std::cout << std::endl;
    }
}

static inline void print(std::string name, const JBlock& m) {
    std::cout << name << " : " << std::endl;
    for (int i = 0; i < m.rows(); i++) {
        for (int j = 0; j < m.cols(); j++) {
            std::cout << std::scientific << m(i, j) << ", ";
        }
        std::cout << std::endl;
    }
}


// LCP error that is in units of Joules
float unit_consistent_energy_lcp_error(Eigen::Vector3f x_k, Eigen::Matrix3f Aii, Eigen::VectorXf b, const float mu)
{
    float err = 0.0f;
    for (int i = 0; i < x_k.size(); i++) {
        float lower = -mu * x_k(i);
        float upper = mu * x_k(i);
        float w = Aii(i, i) * x_k(i) + b(i);

        // Compute violation of the upperand lower bound
        float delta_x_u = std::max(x_k(i) - upper, 0.0f);
        float delta_x_l = std::max(lower - x_k(i), 0.0f);

        // Compute feasible component x0
        float x0 = std::max(lower, std::min(x_k(i), upper));

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
        float delta_e_w_plus = std::min(((1.0f / (2.0f * (Aii(i, i) + 1e-10f))) * (w_plus * w_plus)), (0.5f * Aii(i, i) * (sigma_l * sigma_l)));

        // Compute the negative velocity energy error
        float delta_e_w_minus = std::min(((1.0f / (2.0f * (Aii(i, i) + 1e-10f))) * (w_minus * w_minus)), (0.5f * Aii(i, i) * (sigma_u * sigma_u)));

        // Compute the energy error per constraint
        float delta_e = std::max(std::max(std::abs(delta_e_x_u), std::abs(delta_e_x_l)), std::max(std::abs(delta_e_w_plus), std::abs(delta_e_w_minus)));

        err += std::abs(delta_e);
    }

    return err;
}



// Loop over all other contacts involving _body
//      and compute : x -= (J0*Minv0*Jother^T) * lambda_other
//
static inline void accumulateCoupledContactsAndJoints(Contact* j, const JBlock& JMinv, RigidBody* body, Eigen::Vector3f& x)
{
    if( body->fixed )
        return;

    for(Contact* cc : body->contacts)
    {
        if( cc != j )
        {
            if( body == cc->body0 )
            {
                //JT = cc->J0.transpose();
                x.noalias() = x - JMinv * (cc->J0.transpose() * cc->lambda);
            }
            else
            {
                //JT = cc->J1.transpose();
                x.noalias() = x - JMinv * (cc->J1.transpose() * cc->lambda);
            }
        }
    }
}

// static inline void accumulateCoupledContactsAndJoints(const JBlock& JMinv, RigidBody* body, Eigen::VectorXf& x)
// {
//     if( body->fixed )
//         return;

//     for(Contact* cc : body->contacts)
//     {
//         if( body == cc->body0 )
//         {
//             //JT = cc->J0.transpose();
//             x.noalias() = x - JMinv * (cc->J0.transpose() * cc->lambda);
//         }
//         else
//         {
//             //JT = cc->J1.transpose();
//             x.noalias() = x - JMinv * (cc->J1.transpose() * cc->lambda);
//         }
//     }
// }

// Loop over all other contacts involving _body
//      and compute : x -= (J0*Minv0*Jother^T) * lambda_other
//
// static inline void accumulateCoupledContactsAndJoints(Contact* j, const JBlock& JMinv, RigidBody* body, Eigen::Vector3f& x, Eigen::Vector3f& b)
// {
//     if( body->fixed )
//         return;

//     for(Contact* cc : body->contacts)
//     {
//         if( cc != j )
//         {
//             Eigen::VectorXf wbar = Eigen::VectorXf::Zero(6);
//             wbar[0] = cc->body0->w[0];
//             wbar[1] = cc->body0->w[1];
//             wbar[2] = cc->body0->w[2];

//             wbar[3] = cc->body1->w[0];
//             wbar[4] = cc->body1->w[1];
//             wbar[5] = cc->body1->w[2];

//             if( body == cc->body0 )
//             {
//                 //JT = cc->J0.transpose();
//                 x.noalias() = b - cc->J0 * wbar;
//             }
//             else
//             {
//                 //JT = cc->J1.transpose();
//                 x.noalias() = b - cc->J1 * wbar;
//             }
//         }
//     }
// }


static inline void solveJoint(const Eigen::LLT<Eigen::MatrixXf>& LLT, const Eigen::VectorXf& b, Eigen::VectorXf& x)
{
    x = LLT.solve(b);
}

static inline void multAndSub(const Eigen::MatrixXf& G, const Eigen::Vector3f& x, const Eigen::Vector3f& y, const float a, Eigen::Vector3f& b)
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
static inline void buildRHS(Contact* c, float h, Eigen::Vector3f& b)
{
    const float gamma = h * c->k / (h * c->k + c->b);       // error reduction parameter
    b = gamma * c->phi / h;

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


//static inline void buildRHS(Contact* c, Eigen::Vector3f& b, const float h)
//{
//    // std::cout << "contact damping " << j->b << std::endl;
//    // std::cout << "contact stiffness " << j->k << std::endl;
//    // std::cout << "hinv " << hinv << std::endl;
//    // std::cout << "gamma " << gamma << std::endl;
//    // std::cout << "phi " << j->phi << std::endl;
//    // printVector3f("b ", b);
//
//    //const float gamma = h * c->k / (h * c->k + c->b);       // error reduction parameter
//    //b = -gamma * c->phi / h;
//
//    if( !c->body0->fixed )
//    {
//        multAndSub(c->J0Minv, c->body0->f, c->body0->tau, 1.0f, b);
//        multAndSub(c->J0, c->body0->xdot, c->body0->omega, h, b);
//    }
//    if( !c->body1->fixed )
//    {
//        multAndSub(c->J1Minv, c->body1->f, c->body1->tau, 1.0f, b);
//        multAndSub(c->J1, c->body1->xdot, c->body1->omega, h, b);
//    }
//}

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


static inline void solveFriction(const Eigen::Matrix3f& A, const Eigen::VectorXf& b, Eigen::Vector3f& x, const float mu)
{
    // Next, friction impulses are projected to [-mu * x(0), mu * x(1)]
    //
    x(1) = std::max(-mu * x(0), std::min(mu * x(0), (b(1) - A(1, 0) * x(0) - A(1, 2) * x(2)) / A(1, 1)));
    x(2) = std::max(-mu * x(0), std::min(mu * x(0), (b(2) - A(2, 0) * x(0) - A(2, 1) * x(1)) / A(2, 2)));
}

// z_k = x_k - R_kk * ( wbar + b_k )
static inline void compute_z_k(Eigen::Vector3f& x_k, Eigen::VectorXf& w0, Eigen::VectorXf& w1, const Eigen::MatrixXf& R_k, 
    const Eigen::Vector3f& b, Eigen::Vector3f& z_k, size_t k, JBlock& J0, JBlock& J1) {
    // DEBUG: Check component values

    //print("w0", w0);
    //print("w1", w1);
    //print("b", b);
    //print("R_k", R_k);

    z_k = Eigen::Vector3f::Zero();

//    z_k = J * wbar;
    z_k = J0 * w0 + J1 * w1;

    // Add b_k
    z_k += b;

    // temp += R_kk * z_k
    Eigen::Vector3f temp = Eigen::Vector3f::Zero();
    temp = R_k * z_k;

    // Substract x_k - temp
    z_k = x_k - temp;

    //print("z_k", z_k);
}

static inline void initialize_R(const std::vector<Eigen::Matrix3f>& A, std::vector<Eigen::MatrixXf>& R, std::vector<Eigen::MatrixXf>& nu) {
    for (size_t i=0; i<R.size(); i++) {
        // Set to zero
        R[i] = Eigen::Matrix3f::Zero();
        nu[i] = Eigen::Matrix3f::Zero();

        for (int u=0; u<R[i].rows(); u++) {
            //R[i](u,u) = 0.5f;
            //nu[i](u,u) = 0.7f;

            R[i](u,u) = 1.0f / A[i](u,u);
            nu[i](u,u) = 0.9f;
        }
    }
}

// lambda_n = prox_{R^+} (lambda_n - r (A * lambda_n + b))
static inline void normal_solver(Eigen::Vector3f& z_k, Eigen::Vector3f& x_k) {
    // Max strategy
    x_k[0] = std::max(0.0f, z_k[0]);
    x_k[1] = z_k[1];
    x_k[2] = z_k[2];
}

// lambda_f = prox_C (lambda_f - r (A * lambda_f + b))
static inline void friction_solver(Eigen::Vector3f& z_k, Eigen::Vector3f& x_k, Eigen::Vector4f& coeffs) {
    // 0: n
    // 1: s
    // 2: t
    // Set all directions to the same mu for now

    using std::sqrt;
    using std::max;
    using std::fabs;

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
    const float scale = max(1.0f, max( a, max( b, max(fabs( z_s ), fabs( z_t )) ) ) );

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
    float t1 = max(sa, sb) * sqrt(xx + yy);
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
        if (fabs(t1-t0) < tol) break;

        const float aat = aa + t_k;
        const float bbt = bb + t_k;
        const float g_k = aaxx/(aat*aat) + bbyy/(bbt*bbt) - 1.0f;
        if (fabs(g_k) < tol) break;
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

static inline bool checkInVector(std::vector<size_t>& v, size_t tag) {
    for (size_t i=0; i<v.size(); i++) {
        if (v[i] == tag) {
            return true;
        }
    }
    return false;
}

static inline void propagateW(RigidBody* body0, RigidBody* body1, float w0, float w1, float w2) {
    std::vector<size_t> alreadyVisitedTags;
    if (!body0->fixed) {
        for (Contact* cc : body0->contacts) {
            // Need to check if body already visited not to visit again

//                std::cout << "body0 loop: " << contacts[k]->body0->tag << std::endl;

            // Use the body that is not body0
            if (cc->body0->tag != body0->tag) {

                // Make sure it is also not the other body involved in the kth contact
                if (cc->body0->tag != body1->tag) {
                    if (!checkInVector(alreadyVisitedTags, cc->body0->tag)) {
                        alreadyVisitedTags.push_back(cc->body0->tag);

//                        std::cout << "cc->body0 tag: " << cc->body0->tag  << std::endl;

                        cc->body0->w[0] += w0;
                        cc->body0->w[1] += w1;
                        cc->body0->w[2] += w2;
                    }
                }
            } else {
                if (cc->body1->tag != body1->tag) {
                    if (!checkInVector(alreadyVisitedTags, cc->body1->tag)) {
                        alreadyVisitedTags.push_back(cc->body1->tag);

//                        std::cout << "cc->body1 tag: " << cc->body1->tag  << std::endl;

                        cc->body1->w[0] += w0;
                        cc->body1->w[1] += w1;
                        cc->body1->w[2] += w2;
                    }
                }
            }
        }
    }
}
static inline void accumulateCoupledContacts(Contact *c, const JBlock &JMinv, RigidBody *body, Eigen::Vector3f &b)
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

// Update the w of all the related bodies
//            propagateW(contacts[k]->body0, contacts[k]->body1, temp[3], temp[4], temp[5]);
//            propagateW(contacts[k]->body1, contacts[k]->body0, temp[0], temp[1], temp[2]);



// 2 body implementation
void SolverPROX::solveContact(std::vector<Eigen::Matrix3f>& A, std::vector<Eigen::Vector3f>& b, std::vector<Eigen::Vector3f>& lambda, 
    std::vector<Eigen::Vector4f>& coeffs, std::vector<Eigen::MatrixXf>& J, std::vector<Eigen::MatrixXf>& MinvJT, std::vector<Contact*>& contacts, 
    const float h, std::vector<Eigen::Matrix3f>& Aii) {
    std::vector<Eigen::Vector3f> x;
    x.resize(A.size());

    for (size_t i=0; i<x.size(); i++) {
        x[i] = Eigen::Vector3f::Zero();
    }

    x = lambda; // warm starting

    int abs_conv_in_iteration = -1;
    int rel_conv_in_iteration = -1;
    size_t count_divergence = 0;

    float relative_threshold = 0.0001;
    float absolute_threshold = 0.0001;

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

    const bool is_in = timestep_samples.find(m_rigidBodySystem->getFrameCount()) != timestep_samples.end();
    if (is_in || ALL_SAMPLES)
        prox_lcp_energy_error_file << "#" << m_rigidBodySystem->getFrameCount() << "\n";
        prox_merit_file << "#" << m_rigidBodySystem->getFrameCount() << "\n";

    // GS Loop
    for (size_t iter=0; iter<m_maxIter; iter++) {
        // Loop over contact points

        for (size_t i = 0; i < contacts.size(); i++) {
            contacts[i]->w = MinvJT[i] * lambda[i];
        }

        for (size_t k=0; k<A.size(); k++) {
            // accumulateCoupledContacts(contacts[k], contacts[k]->Minv0, contacts[k]->body0, b[k]);
            // accumulateCoupledContacts(contacts[k], contacts[k]->Minv1, contacts[k]->body1, b[k]);

            for (Contact *cc : contacts[k]->body0->contacts)
            {
                if (!contacts[k]->body0->fixed && cc != contacts[k])
                {
                    if (contacts[k]->body0 == cc->body0)
                    {
                        // contacts[k]->w += cc->w;
                        b[k] -= contacts[k]->J0Minv * (cc->J0.transpose() * cc->lambda);
                    }
                    else
                    {
                        // contacts[k]->w += cc->w;
                        b[k] -= contacts[k]->J1Minv * (cc->J1.transpose() * cc->lambda);
                    }
                }
            }

            // std::cout << "=============================================" << std::endl;
            x_k = x[k];
            deltax = x[k];

            // print("z_k before", z_k);
            // print("init x_k", x_k);
            // print("init deltax", deltax);
            // print("b_k", b[k]);
            // print("R_k", R[k]);
            // print("J_k", J[k]);
            // print("contact_k", contacts[k]->w);

            // Compute z_k
            //compute_z_k(x_k, contacts[k]->w, R[k], b[k], z_k, k, J[k]);

        //    print("z_k after", z_k);

            // Project the normal component
            normal_solver(z_k, x_k);

            x_k[1] = z_k[1];
            x_k[2] = z_k[2];

            // print("post normal x_k", x_k);
            // print("post normal deltax", deltax);

            // printVector3f("x_k before", x_k);

            // Project the tangential component
            friction_solver(z_k, x_k, coeffs[k]);

            // print("post friction x_k", x_k);
            // print("post friction deltax", deltax);

            // printVector3f("x_k after", x_k);

            // delta_x = x_k_new - x_k_old
            deltax = x_k - lambda[k];

            // Update wbar
            Eigen::VectorXf temp;
            temp = Eigen::VectorXf::Zero(6);
            temp = MinvJT[k] * deltax;

            // print("deltax", deltax);
            // print("MinvJT k", MinvJT[k]);

            // accumulateCoupledContacts(contacts[k], contacts[k]->J0Minv, contacts[k]->body0, b[k]);
            // accumulateCoupledContacts(contacts[k], contacts[k]->J1Minv, contacts[k]->body1, b[k]);

             //contacts[k]->w += temp;

            // print("kth contact w", contacts[k]->w);

                           std::vector<size_t> alreadyVisistedContacts;
                        //    std::cout << "Loop of body0: " << contacts[k]->body0->tag << " \twith contact: " << contacts[k]->tag << std::endl;
                           for (size_t i=0; i<contacts[k]->body0->contacts.size(); i++) {
               //                std:: cout << "Loop of body0" << std::endl;
                               Contact* cc = contacts[k]->body0->contacts[i];
                               if (!checkInVector(alreadyVisistedContacts, cc->tag)) {
                                   alreadyVisistedContacts.push_back(cc->tag);

                                   // Check if contact has the same bodies as contact_k
                                   if ((contacts[k]->body0->tag == cc->body0->tag && contacts[k]->body1->tag == cc->body1->tag)
                                           || (contacts[k]->body1->tag == cc->body0->tag && contacts[k]->body0->tag == cc->body1->tag)
                                           || (contacts[k]->body0->tag == cc->body1->tag && contacts[k]->body1->tag == cc->body0->tag)
                                           ) {
                                    //   std::cout << "Contact: " << cc->tag << std::endl;
                                       cc->w += temp;
                                    //    print("cc w", cc->w);
                                   }
                               }
                           }

            //            std::cout << "++++++++++++++" << std::endl;
            //            for (size_t i=0; i<contacts[k]->body1->contacts.size(); i++) {
            //                std:: cout << "Loop of body1" << std::endl;
            //                Contact* cc = contacts[k]->body1->contacts[i];
            //                if (!checkInVector(alreadyVisistedContacts, cc->tag)) {
            //                    alreadyVisistedContacts.push_back(cc->tag);

            //                    // Check if contact has the same bodies as contact_k
            //                    if ((contacts[k]->body0->tag == cc->body0->tag && contacts[k]->body1->tag == cc->body1->tag)
            //                        || (contacts[k]->body1->tag == cc->body0->tag && contacts[k]->body0->tag == cc->body1->tag)
            //                        || (contacts[k]->body0->tag == cc->body1->tag && contacts[k]->body1->tag == cc->body0->tag)
            //                        ) {
            //                        std::cout << "Contact: " << cc->tag << std::endl;
            //                        cc->w += temp;
            //                    }
            //                }
            //            }
            //            std::cout << "======================" << std::endl;
            //            contacts[k]->body0->w[0] += temp[0];
            //            contacts[k]->body0->w[1] += temp[1];
            //            contacts[k]->body0->w[2] += temp[2];

            //            contacts[k]->body1->w[0] += temp[3];
            //            contacts[k]->body1->w[1] += temp[4];
            //            contacts[k]->body1->w[2] += temp[5];

            //            propagateW(contacts[k]->body0, contacts[k]->body1, temp[3], temp[4], temp[5]);
            //            propagateW(contacts[k]->body1, contacts[k]->body0, temp[0], temp[1], temp[2]);

            // residual = lambda^k - lambda^(k+1)
            residual[k] = lambda[k] - x_k;

            x[k] = x_k;
            if (is_in || ALL_SAMPLES)
                prox_lcp_energy_error_file << unit_consistent_energy_lcp_error(x_k, Aii[k], b[k], contacts[k]->mu) << " ";
        }

        residual_norm = infiniteNorm(residual);
        if (is_in || ALL_SAMPLES)
            prox_merit_file << residual_norm << "\n";
            prox_lcp_energy_error_file << "\n";


//        outf2 << residual_norm << '\n';

        // Check if residual is under absolute threshold
        if (residual_norm < absolute_threshold) {
//            std::cout << "Absolute convergence in " << iter << " iterations \t |residual| = " << residual_norm << std::endl;
            abs_conv_in_iteration = iter;
            //conv_file << "abs_conv: " << abs_conv_in_iteration << " rel_conv: " << -1 << "divergence count: " << count_divergence << '\n';
            break;
        }

        // Check if residual norm difference is under relative threshold
        if (std::fabs(residual_norm-last_residual_norm) < relative_threshold*last_residual_norm) {
//            std::cout << "Absolute convergence in " << iter << " iterations \t |residual| = " << residual_norm << std::endl;
            rel_conv_in_iteration = iter;
            //conv_file << "abs_conv: " << -1 << " rel_conv: " << rel_conv_in_iteration << "divergence count: " << count_divergence << '\n';
            break;
        }

        if (residual_norm > last_residual_norm) {
            for (size_t j=0; j<R.size(); j++) {
                R[j] = nu[j] * R[j];
            }
            x = lambda;
            ++count_divergence;
        } else {
            last_residual_norm = residual_norm;
            lambda = x;
        }

    }

    lambda = x;
    // for (int i=0; i<lambda.size(); i++) {
    //     for (int j=0; j<lambda[i].size(); j++) {
    //         std::cout << lambda[i][j] << ", ";
    //     }
    //     std::cout << std::endl;
    // }
    // std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
    conv_file << abs_conv_in_iteration << "," << rel_conv_in_iteration << "," << count_divergence << "," << residual_norm  << '\n';
}

SolverPROX::SolverPROX(RigidBodySystem* _rigidBodySystem) : Solver(_rigidBodySystem)
{

}

void SolverPROX::solve(float h)
{
    std::vector<Contact*>& contacts = m_rigidBodySystem->getContacts();
    std::vector<RigidBody*>& bodies = m_rigidBodySystem->getBodies();
    const int numContacts = contacts.size();

    if (!prox_lcp_energy_error_file) {
        std::cerr << "uh oh, could not open/create data.csv" << std::endl;
    }

    // Clear w from each body
    for (size_t i=0; i<bodies.size(); i++) {
        bodies[i]->w = Eigen::Vector3f::Zero();
        bodies[i]->tag = i;
    }

    // Build diagonal matrices of contacts
    std::vector<Eigen::Matrix3f> Acontacts;
    std::vector<Eigen::Matrix3f> Acontactii;
    std::vector<Eigen::MatrixXf> AJ;
    //std::vector<Eigen::MatrixXf> AJ_conc;
    //std::vector<Eigen::MatrixXf> Minv_conc;
    std::vector<Eigen::MatrixXf> MinvJT;
    std::vector<Eigen::Vector3f> lambdasContacts;
    std::vector<Eigen::Vector4f> coeffs;
    std::vector<Eigen::Vector3f> bContacts;
    if( numContacts > 0)
    {
        // Build diagonal matrices
        Acontacts.resize(numContacts);
        lambdasContacts.resize(numContacts);
        coeffs.resize(numContacts);
        AJ.resize(numContacts);
        MinvJT.resize(numContacts);
        //AJ_conc.resize(numContacts);
        //Minv_conc.resize(numContacts);

        // Build diagonal matrices (For LCP energy error computation)
        Acontactii.resize(numContacts);
        for (int i = 0; i < numContacts; ++i)
        {
            Contact* c = contacts[i];
            const float eps = 1.0f / (h * h * c->k + h * c->b);    // constraint force mixing

            // Compute the diagonal term : Aii = J0*Minv0*J0^T + J1*Minv1*J1^T
            //
            Acontactii[i].setZero(3, 3);
            Acontactii[i](0, 0) += eps;

            if (!c->body0->fixed)
            {
                Acontactii[i] += c->J0Minv * c->J0.transpose();
            }
            if (!c->body1->fixed)
            {
                Acontactii[i] += c->J1Minv * c->J1.transpose();
            }
        }

        for(int i = 0; i < numContacts; ++i)
        {

            Contact* c = contacts[i];
            c->tag = i;
            contacts[i]->w = Eigen::Vector3f::Zero();

            contacts[i]->m_CFM = 1.0f / (h * h * c->k + h * c->b);    // constraint force mixing
            //contacts[i]->m_CFM = 0;

            c->w = Eigen::VectorXf::Zero(6);

            // Compute the diagonal term : Aii = J0*Minv0*J0^T + J1*Minv1*J1^T
            //
            Acontacts[i] = Eigen::Matrix3f::Zero();
            Acontacts[i](0,0) += contacts[i]->m_CFM;

            AJ[i] = Eigen::MatrixXf::Zero(3,6);
            MinvJT[i] = Eigen::MatrixXf::Zero(6,3);

            if( !c->body0->fixed )
            {
                JBlockTranspose JT = c->J0.transpose();
                Acontacts[i] += c->J0Minv * JT;
                AJ[i] += c->J0;
                MinvJT[i] += c->J0Minv.transpose();
            }
            if( !c->body1->fixed )
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

            // printMatrix3f("Acontacts[i]", Acontacts[i]);
            // printMatrixXf("AJ[i]", AJ[i]);
            // printMatrixXf("MinvJT[i]", MinvJT[i]);
            // printVector3f("lambdasContacts[i]", lambdasContacts[i]);
        }

        // Assemble AJ_conc and M_conc
        //for (int i = 0; i < numContacts; ++i) {
        //    Contact* c = contacts[i];
        //    AJ_conc[i] = Eigen::MatrixXf::Zero(3, 12);
        //    for (int j = 0; j < 6; j++) {
        //        AJ_conc[i](0, j) = c->J0(0, j);
        //        AJ_conc[i](1, j) = c->J0(1, j);
        //        AJ_conc[i](2, j) = c->J0(2, j);
        //    }
        //    for (int j = 6; j < 12; j++) {
        //        AJ_conc[i](0, j) = c->J1(0, j);
        //        AJ_conc[i](1, j) = c->J1(1, j);
        //        AJ_conc[i](2, j) = c->J1(2, j);
        //    }
        //    Minv_conc[i](0, 0) = c->Minv0(0, 0);
        //    Minv_conc[i](1, 1) = c->Minv0(1, 1);
        //    Minv_conc[i](2, 2) = c->Minv0(2, 2);

        //    Minv_conc[i](3, 3) = c->Minv1(0, 0);
        //    Minv_conc[i](4, 4) = c->Minv1(1, 1);
        //    Minv_conc[i](5, 5) = c->Minv1(2, 2);
        //}

        // Assemble b
        bContacts.resize(numContacts);
        for (int i=0; i<numContacts; ++i) {
            Contact* c = contacts[i];
            bContacts[i] = Eigen::Vector3f::Zero();
            const float hinv = 1.0f / h;
            contacts[i]->m_ERP = c->k / (c->k + hinv * c->b); // error reduction parameter
            //contacts[i]->m_ERP = 0;
            //bContacts[i] = hinv * contacts[i]->m_ERP * c->phi;
            buildRHS(c, h, bContacts[i]);

            //buildRHS(c, bContacts[i], hinv);
            //for (int j = 0; j < bContacts[i].size(); j++) {
            //    bContacts[i](j) = -bContacts[i](j);
            //}
            c->lambda.setZero();
        }
    }

    if (numContacts > 0) {
//        // DEBUGGING: Check size of lambdas before
//        for (size_t i=0; i<lambdasContacts.size(); i++) {
//            std::cout << i << " lambda: " << "(" << lambdasContacts[i][0] << ", " << lambdasContacts[i][1] << ", " << lambdasContacts[i][2] << ")" << std::endl;
//        }

        //solveContact(Acontacts, bContacts, lambdasContacts, coeffs, AJ, MinvJT, contacts, h, Acontactii);


        // Rewrite PROX, who knows...
        //std::cout << "GATE-00" << std::endl;
        std::vector<RigidBody*> bodies = m_rigidBodySystem->getBodies();
        for (int i = 0; i < bodies.size(); i++) {
            bodies[i]->w = Eigen::VectorXf::Zero(6);
        }
        //std::cout << "GATE-01" << std::endl;

        std::vector<Eigen::Vector3f> lambdaDummy;
        lambdaDummy.resize(Acontacts.size());
        for (int i=0; i<numContacts; ++i) {
            lambdaDummy[i] = Eigen::Vector3f::Zero();
        }
        
        std::vector<Eigen::MatrixXf> R;
        std::vector<Eigen::MatrixXf> nu;
        Eigen::Vector3f z_k = Eigen::Vector3f::Zero();

        std::vector<Eigen::Vector3f> residual;
        residual.resize(Acontacts.size());

        float residual_norm = 0;
        float last_residual_norm = FLT_MAX;

        R.resize(Acontacts.size());
        nu.resize(Acontacts.size());

        initialize_R(Acontacts, R, nu);

        for (int iter = 0; iter < m_maxIter; ++iter) {
            for (int i = 0; i < bodies.size(); i++) {
                bodies[i]->w = Eigen::VectorXf::Zero(6);
            }
            for (int i = 0; i < numContacts; ++i) {
                lambdaDummy[i] = Eigen::Vector3f::Zero();
                if (!contacts[i]->body0->fixed)
                    contacts[i]->body0->w += contacts[i]->Minv0 * contacts[i]->J0.transpose() * lambdasContacts[i];
                if (!contacts[i]->body1->fixed)
                    contacts[i]->body1->w += contacts[i]->Minv1 * contacts[i]->J1.transpose() * lambdasContacts[i];
            }
            for (int i = 0; i < numContacts; ++i) {
                compute_z_k(lambdasContacts[i], contacts[i]->body0->w, contacts[i]->body1->w, R[i], bContacts[i], z_k, i, contacts[i]->J0, contacts[i]->J1);
                normal_solver(z_k, lambdaDummy[i]);

                friction_solver(z_k, lambdaDummy[i], coeffs[i]);

                residual[i] = lambdaDummy[i] - lambdasContacts[i];
                if (!contacts[i]->body0->fixed)
                    contacts[i]->body0->w += contacts[i]->Minv0 * contacts[i]->J0.transpose() * (lambdaDummy[i] - lambdasContacts[i]);

                if (!contacts[i]->body1->fixed)
                    contacts[i]->body1->w += contacts[i]->Minv1 * contacts[i]->J1.transpose() * (lambdaDummy[i] - lambdasContacts[i]);

            }

            residual_norm = infiniteNorm(residual);
            if (residual_norm > last_residual_norm) {
                for (size_t j = 0; j < R.size(); j++) {
                    R[j] = nu[j] * R[j];
                }
            }
            else {
                last_residual_norm = residual_norm;
                for (int j = 0; j < numContacts; j++) {
                    lambdasContacts[j] = lambdaDummy[j];
                }
            }
        }


        // Set the new lambdas to the old contact lambdas
        for (size_t i=0; i<lambdasContacts.size(); i++) {
            contacts[i]->lambda = lambdasContacts[i];
        }

        //for (int i = 0; i < numContacts; i++) {
        //    //pgs_comp_error_file << get_complementarity_error(contacts[i], Acontactii[i], b[i]) << " ";
        //    prox_lcp_energy_error_file << unit_consistent_energy_lcp_error(contacts[i], Acontactii[i], bContacts[i]) << " ";
        //}
        //prox_lcp_energy_error_file << '\n';

        // DEBUGGING: Check size of lambdas after
    //     for (size_t i = 0; i < lambdasContacts.size(); i++)
    //     {
    //         std::cout << i << " lambda: "
    //                   << "(" << std::scientific << contacts[i]->lambda[0] << ", " << std::scientific << contacts[i]->lambda[1] << ", " << std::scientific << contacts[i]->lambda[2] << ")"
    //                   << "\tBody0: " << contacts[i]->body0->tag << "\tBody1: " << contacts[i]->body1->tag
    //                   << std::endl;
    //     }
    //    std::cout << "=============================================================================" << std::endl;
    }
}