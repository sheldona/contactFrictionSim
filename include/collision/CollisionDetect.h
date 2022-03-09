#pragma once

#include <Eigen/Dense>
#include <vector>

class Contact;
class RigidBody;
class RigidBodySystem;

//
// The main collision detection class,
// it will test the geometries of each pair of rigid bodies
// and populate an array with contacts.
//
class CollisionDetect
{
public:

    // Constructor.
    //
    CollisionDetect(RigidBodySystem* rigidBodySystem);

    // Tests for collisions between all pairs of bodies in the rigid body system
    // and generates contacts for any intersecting geometries.
    // The array of generated contacts can be retrieved by calling getContacts().
    //
    void detectCollisions();

    void clear();

    // Compute the Jacobians for contacts
    void computeContactJacobians();

    const std::vector<Contact*>& getContacts() const { return m_contacts; }

    std::vector<Contact*>& getContacts() { return m_contacts; }    
    
    const std::vector<Contact*>& getSubContacts() const { return m_subContacts; }

    std::vector<Contact*>& getSubContacts() { return m_subContacts; }

private:

    // Sphere-sphere collision test.
    // Assumes that both @a body0 and @a body1 have a sphere collision geometry.
    //
    void collisionDetectSphereSphere(RigidBody* body0, RigidBody* body1);

    // Sphere-box collision test.
    void collisionDetectSphereBox(RigidBody* body0, RigidBody* body1);

    // Box-plane collision test
    void collisionDetectBoxPlane(RigidBody* body0, RigidBody* body1);

    // Box-plane collision test - placing contact points in a regular grid
    void collisionDetectBoxPlaneGrid(RigidBody* body0, RigidBody* body1);

    // Box-plane collision test - placing contact points randomly on the contact surface
    void collisionDetectBoxPlaneRandom(RigidBody* body0, RigidBody* body1);

    // Box-SDF collision test
    void collisionDetectBoxSdf(RigidBody* body0, RigidBody* body1);

    // SDF-SDF collision test
    void collisionDetectSdfSdf(RigidBody* body0, RigidBody* body1);

    // Generate a random sample point given the coordinates of a triangle using barycentric coordinates
    Eigen::Vector3f generateSamplePoint(const std::array<Eigen::Vector3f, 3>& vertices);

private:

    RigidBodySystem* m_rigidBodySystem;
    std::vector<Contact*> m_contacts;
    std::vector<Contact*> m_subContacts;

};
