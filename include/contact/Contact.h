#pragma once

#include <Eigen/Dense>

typedef Eigen::Matrix<float, -1, 6, Eigen::RowMajor> JBlock;
typedef Eigen::Matrix<float, 6, -1, Eigen::ColMajor> JBlockTranspose;

class RigidBody;

// Non-interpenetration contact constraint.
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

    void computeContactFrame(const Eigen::Vector3f& dir)
    {
        t1 = dir - (dir.dot(n))*n;
        // Compute tangent directions.
        //t2 = n.cross( dir );
        if( t1.norm() < 1e-5f )
        {
            t1 = -n.cross( Eigen::Vector3f(0,0,-1) );
        }
        t1.normalize();
        t2 = n.cross(t1);
        t2.normalize();
    }

    RigidBody* body0;           //< The first body
    RigidBody* body1;           //< The second boy
    JBlock J0;                  //< The Jacobian of body0
    JBlock J1;                  //< The Jacobian of body1
    JBlock J0Minv;
    JBlock J1Minv;
    Eigen::VectorXf phi;        //< Contraint error
    Eigen::VectorXf lambda;     //< Constraint impulse
    float k;                    //< Contact stiffness
    float b;
    unsigned int index;         // Used for global indexing (auxiliary)

protected:

    // Default constructor.
    explicit Contact();
};
