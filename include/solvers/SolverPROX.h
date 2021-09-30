#pragma once

#include <Eigen/Dense>
#include <vector>

#include "contact/Contact.h"
#include "solvers/Solver.h"
#include "rigidbody/RigidBody.h"
#include "rigidbody/RigidBodySystem.h"

// Proximal Algorithm (PROX) solver for the LCP.
// This implementation uses a matrix-free approach, where only
// the non-zero blocks of the lead matrix are assembled and stored.
//
class SolverPROX : public Solver
{
public:
    SolverPROX(RigidBodySystem* _rigidBodySystem);

    virtual void solve(float h);

protected:
    void solveContact(std::vector<Eigen::Matrix3f> &A, std::vector<Eigen::Vector3f> &b, std::vector<Eigen::Vector3f> &lambda, 
                        std::vector<Eigen::Vector4f> &coeffs, std::vector<Eigen::MatrixXf> &J, std::vector<Eigen::MatrixXf> &MinvJT, 
                        std::vector<Contact *> &contacts);
};