#pragma once

#include "solvers/Solver.h"

// Proximal Algorithm (PROX) solver for the LCP.
// This implementation uses a matrix-free approach, where only
// the non-zero blocks of the lead matrix are assembled and stored.
//
class SolverPROX : public Solver
{
public:
    SolverPROX(RigidBodySystem* _rigidBodySystem);

    virtual void solve(float h);
};