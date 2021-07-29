#pragma once

class RigidBodySystem;

// Projected Gauss-Seidel (PGS) Boxed LCP solver.
// This implementation uses a matrix-free approach, where only
// the non-zero blocks of the lead matrix are assembled and stored.
//
class SolverBoxPGS
{
public:

    SolverBoxPGS(RigidBodySystem* _rigidBodySystem);

    void setMaxIter(int _maxIter) { m_maxIter = _maxIter; }

    int getMaterIter() const { return m_maxIter; }

    // Implement solve method that solves for the constraint forces in @a m_rigidBodySystem.
    //
    virtual void solve(float h);

protected:

    RigidBodySystem* m_rigidBodySystem;
    int m_maxIter;

};
