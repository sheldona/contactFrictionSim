#pragma once

class RigidBodySystem;

class SolverBoxBPP
{
public:
    SolverBoxBPP(RigidBodySystem* _rigidBodySystem);

    virtual void solve(float h);

protected:

    RigidBodySystem* m_rigidBodySystem;
    unsigned int m_maxIter;

};
