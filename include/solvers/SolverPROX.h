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

private:
    std::vector<Eigen::VectorXf> w;     // w = Minv * JT * lambda_k
};