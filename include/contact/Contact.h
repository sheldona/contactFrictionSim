#pragma once

#include <Eigen/Dense>

typedef Eigen::Matrix<float, -1, 6, Eigen::RowMajor> JBlock;
typedef Eigen::Matrix<float, 6, -1, Eigen::ColMajor> JBlockTranspose;

class RigidBody;

// Contact constraint with box friction.
// This class stores the contact normal @a n and contact point @a p,
// as well as the Jacobians @a J0 and @a J1 for each body,
// which are computed by computeJacobian().
//
class Contact
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

public:

    // Constructor with all parameters.
    Contact(RigidBody* _body0, RigidBody* _body1, const Eigen::Vector3f& _p, const Eigen::Vector3f& _n, float _phi);

    virtual ~Contact();

    Eigen::Vector3f p;          // The contact point.
    Eigen::Vector3f n;          // The contact normal.
    Eigen::Vector3f t1, t2;     // Tangent directions.
    float mu;                   // Coefficient of friction.
    float pene;                 // Penetration

    virtual void computeJacobian();

    // Computes a contact frame using the contact normal @a n
    // and the provided direction @a dir, which is aligned with the first tangent direction.
    // The resulting frame is stored in the basis vectors @a n, @a t1, and @a t2.
    //
    void computeContactFrame(const Eigen::Vector3f& dir)
    {
        // Compute first tangent direction.
        //
        t1 = dir - (dir.dot(n))*n;
        if( t1.norm() < 1e-5f )
        {
            // Fail-safe: use axis-aligned direction (0,0,-1) 
            t1 = -n.cross( Eigen::Vector3f(0,0,-1) );
        }
        t1.normalize();

        // Compute second tangent direction.
        //
        t2 = n.cross(t1);
        t2.normalize();
    }

    RigidBody* body0;           // The first body
    RigidBody* body1;           // The second boy
    JBlock J0;                  // The Jacobian of body0
    JBlock J1;                  // The Jacobian of body1
    JBlock J0Minv;              // The Jacobian post-multiplied by the mass inverse (body 0).
    JBlock J1Minv;              // The Jacobian post-multiplied by the mass inverse (body 1).
    Eigen::VectorXf phi;        // Contraint error
    Eigen::VectorXf lambda;     // Constraint impulse
    float k;                    // Contact stiffness (Baumgarte stabilization)
    float b;                    // Contact damping (Baumgarte stabilization)
    unsigned int index;         // Auxiliary variable used for global indexing

protected:

    // Default constructor.
    explicit Contact();
};
