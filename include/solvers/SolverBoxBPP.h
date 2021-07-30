#pragma once

#include "solvers/Solver.h"

class SolverBoxBPP : public Solver
{
public:
    SolverBoxBPP(RigidBodySystem* _rigidBodySystem);

    virtual void solve(float h);

};
