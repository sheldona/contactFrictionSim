#include "solvers/SolverBoxBPP.h"

#include "rigidbody/RigidBodySystem.h"
#include "rigidbody/RigidBody.h"
#include "contact/Contact.h"

#include <Eigen/Dense>

namespace
{
    enum eIndexSet { kFree = 0, kLower, kUpper };


    static inline void multAndSub(const JBlock& G, const Eigen::Vector3f& x, const Eigen::Vector3f& y, const float a, Eigen::Ref<Eigen::VectorXf> b)
    {
            b -= a * G.col(0) * x(0);
            b -= a * G.col(1) * x(1);
            b -= a * G.col(2) * x(2);
            b -= a * G.col(3) * y(0);
            b -= a * G.col(4) * y(1);
            b -= a * G.col(5) * y(2);
    }

    static inline unsigned int initContacts(std::vector<Contact*>& contacts)
    {
        // First pass: count the number of rows/cols
        //
        unsigned int rows = 0;
        for(auto c : contacts)
        {
            c->index = rows;
            rows += c->J0.rows();
        }

        return rows;
    }

    static inline void updateBounds(std::vector<Contact*>& contacts, Eigen::VectorXf& lower, Eigen::VectorXf& upper)
    {
        for(auto c : contacts)
        {
            const unsigned int dim = c->J0.rows();
            // Non-interpenetration row.
            //
            lower(c->index) = 0.0f;
            upper(c->index) = std::numeric_limits<float>::max();

            // Friction rows.
            //
            for(unsigned int i = 1; i < dim; ++i)
            {
                lower(c->index+i) = -c->mu * c->lambda(0);
                upper(c->index+i) = c->mu * c->lambda(0);
            }
        }
    }

    static inline void buildRHS(const std::vector<Contact*>& contacts, float h, Eigen::VectorXf& b)
    {
        for(auto c : contacts)
        {
            const float hinv = 1.0f / h;
            const float gamma = 0.9f;//h * c->k / (h * c->k + c->b); // error reduction parameter
            const unsigned int dim = c->J0.rows();
            b.segment(c->index, dim) = hinv * gamma * c->phi;

            multAndSub(c->J0, c->body0->xdot, c->body0->omega, 1.0f, b.segment(c->index, dim));
            multAndSub(c->J1, c->body1->xdot, c->body1->omega, 1.0f, b.segment(c->index, dim));

            if( !c->body0->fixed )
            {
                multAndSub(c->J0Minv, c->body0->f, c->body0->tau, h, b.segment(c->index, dim));
            }
            if( !c->body1->fixed )
            {
                multAndSub(c->J1Minv, c->body1->f, c->body1->tau, h, b.segment(c->index, dim));
            }
        }
    }

    static inline void buildMatrix(const std::vector<Contact*>& contacts, float h, Eigen::MatrixXf& A)
    {
        for(auto c : contacts)
        {
            const unsigned int dim = c->J0.rows();
            const float eps = 1.0f / (h * h * c->k + h * c->b);    // constraint force mixing

            A.block(c->index, c->index, dim, dim) = 1e-4f * Eigen::MatrixXf::Identity(dim, dim);
            //A(c->index, c->index) += eps;

            if( !c->body0->fixed )
            {
                A.block(c->index, c->index, dim, dim) += c->J0Minv * c->J0.transpose();
                for(auto cc : c->body0->contacts)
                {
                    if( cc != c )
                    {
                        const int ddim = c->J0.rows();
                        if( cc->body0 == c->body0 )
                        {
                            A.block(c->index, cc->index, dim, ddim) += c->J0Minv * cc->J0.transpose();
                        }
                        else
                        {
                            A.block(c->index, cc->index, dim, ddim) += c->J0Minv * cc->J1.transpose();
                        }
                    }
                }
            }


            if( !c->body1->fixed )
            {
                A.block(c->index, c->index, dim, dim) += c->J1Minv*c->J1.transpose();

                for(auto cc : c->body1->contacts)
                {
                    if( cc != c )
                    {
                        const int ddim = c->J1.rows();
                        if( cc->body0 == c->body1 )
                        {
                            A.block(c->index, cc->index, dim, ddim) += c->J1Minv * cc->J0.transpose();
                        }
                        else
                        {
                            A.block(c->index, cc->index, dim, ddim) += c->J1Minv * cc->J1.transpose();
                        }
                    }
                }
            }
        }

    }

    static inline unsigned int pivot(Eigen::VectorXi& idx, Eigen::VectorXf& z, const Eigen::VectorXf& l, const Eigen::VectorXf& u, const Eigen::VectorXf& w)
    {
        static const float tol = 1e-4f;

        unsigned int numPivots = 0;
        const unsigned int n = idx.rows();
        for(unsigned int j = 0; j < n && (numPivots == 0); ++j)
        {
            int new_idx = kFree;

            if(idx[j] == kFree && z[j] <= l[j])
            {
                new_idx = kLower;
            }
            else if(idx[j] == kFree && z[j] >= u[j])
            {
                new_idx = kUpper;
            }
            else if(idx[j] == kLower && w[j] > -tol)
            {
                new_idx = kLower;
            }
            else if(idx[j] == kUpper && w[j] < tol)
            {
                new_idx = kUpper;
            }

            if(new_idx != idx[j])
            {
                ++numPivots;
                idx[j] = new_idx;
            }
        }
        return numPivots;
    }



    static inline unsigned int tightIndices(const Eigen::VectorXi& idx, std::vector<int>& tightIdx)
    {
        const unsigned int n = idx.rows();
        unsigned int numTight = 0;
        tightIdx.clear();
        for(unsigned int i = 0; i < n; ++i)
        {
            if( idx[i] == kLower || idx[i] == kUpper )
            {
                ++numTight;
                tightIdx.push_back(i);
            }
        }
        return numTight;
    }


    static inline unsigned int freeIndices(const Eigen::VectorXi& idx, std::vector<int>& freeIdx)
    {
        const unsigned int n = idx.rows();
        unsigned int numFree = 0;
        freeIdx.clear();
        for(unsigned int i = 0; i < n; ++i)
        {
            if( idx[i] == kFree )
            {
                ++numFree;
                freeIdx.push_back(i);
            }
        }
        return numFree;
    }

    static inline void solvePrincipalSubproblem(const Eigen::MatrixXf& A,
                                                const Eigen::VectorXf& b,
                                                const Eigen::VectorXi& idx,
                                                const Eigen::VectorXf& l,
                                                const Eigen::VectorXf& u,
                                                Eigen::VectorXf& x)
    {
        std::vector<int> freeIdx, tightIdx;
        const unsigned int numTight = tightIndices(idx, tightIdx);
        const unsigned int numFree = freeIndices(idx, freeIdx);

        if( numFree > 0 )
        {
            Eigen::MatrixXf Aff(numFree, numFree);
            Eigen::VectorXf bf(numFree);

            // Build sub-matrix and rhs vector using free indices
            //
            for(unsigned int j = 0; j < numFree; ++j)
            {
                for(unsigned int i = 0; i < numFree; ++i)
                {
                    Aff(i,j) = A(freeIdx[i], freeIdx[j]);
                }
                bf(j) = b(freeIdx[j]);
            }

            // Update rhs vector with tight indices
            //
            for(unsigned int j = 0; j < numTight; ++j)
            {
                for(unsigned int i = 0; i < numFree; ++i)
                {
                    if( idx[tightIdx[j]] == kLower )
                    {
                        bf(i) -= A(i, j) * l(tightIdx[j]);
                    }
                    else if( idx[tightIdx[j]] == kUpper )
                    {
                        bf(i) -= A(i, j) * u(tightIdx[j]);
                    }
                }
            }

            Eigen::LDLT<Eigen::MatrixXf> ldlt(Aff);
            const Eigen::VectorXf xf = ldlt.solve(bf);

            for(unsigned int i = 0; i < numFree; ++i)
            {
                x(freeIdx[i]) = xf(i);
            }

            for(unsigned int i = 0; i < numTight; ++i)
            {
                if( idx[tightIdx[i]] == kLower )
                {
                    x(tightIdx[i]) = l(tightIdx[i]);
                }
                else if( idx[tightIdx[i]] == kUpper )
                {
                    x(tightIdx[i]) = u(tightIdx[i]);
                }
            }
        }

    }

    static inline void updateContacts(const Eigen::VectorXf& x, std::vector<Contact*>& contacts)
    {
        // Distribute to the contacts
        for(auto c : contacts)
        {
            const unsigned int dim = c->J0.rows();
            c->lambda = x.segment(c->index, dim);
        }
    }

}


SolverBoxBPP::SolverBoxBPP(RigidBodySystem* _rigidBodySystem) :
    m_rigidBodySystem(_rigidBodySystem),
    m_maxIter(50)
{


}

void SolverBoxBPP::solve(float h)
{
    auto contacts = m_rigidBodySystem->getContacts();

    const unsigned int numContacts = contacts.size();
    if( numContacts > 0 )
    {

        const unsigned int dim = initContacts(contacts);

        Eigen::MatrixXf A = Eigen::MatrixXf::Zero(dim,dim);
        Eigen::VectorXf b = Eigen::VectorXf::Zero(dim);
        Eigen::VectorXf x = Eigen::VectorXf::Zero(dim);

        // Update box bounds.
        //
        Eigen::VectorXf lower = Eigen::VectorXf::Zero(dim);
        Eigen::VectorXf upper = Eigen::VectorXf::Zero(dim);
        updateBounds(contacts, lower, upper);

        // Initialize the index set.
        // All variables are initially set to 'free'.
        //
        Eigen::VectorXi idx = Eigen::VectorXi::Constant(dim, kFree);

        buildMatrix(contacts, h, A);
        buildRHS(contacts, h, b);

        for(int k = 0; k < 2; ++k)
        {
            idx.setConstant(kFree);
            for(unsigned int iter = 0; iter < m_maxIter; ++iter)
            {
                x.setZero();
                solvePrincipalSubproblem(A, b, idx, lower, upper, x);

                Eigen::VectorXf w = A*x - b;
                const int numPivots = pivot(idx, x, lower, upper, w);

                if( numPivots == 0 )
                {
                    break;      // Done
                }
            }

            updateContacts(x, contacts);
            updateBounds(contacts, lower, upper);
        }

    }
}

