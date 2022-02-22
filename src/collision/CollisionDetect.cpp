#include "collision/CollisionDetect.h"

#include "contact/Contact.h"
#include "rigidbody/RigidBody.h"
#include "rigidbody/RigidBodySystem.h"

namespace
{
    // Plane-point collision test.
    //
    // Inputs:
    //   p - The point to test.
    //   plane_p - The plane origin.
    //   plane_n  - Direction perpendicular to the plane (the normal).
    // Outputs:
    //   phi - The penetration depth.
    // Returns:
    //   True if the point intersects the plane.
    //
    static inline bool collisionDetectPointPlane(const Eigen::Vector3f& p, const Eigen::Vector3f& plane_p, const Eigen::Vector3f& plane_n, float& phi)
    {
        const float dp = (p - plane_p).dot(plane_n);
        if (dp < 0.0f)
        {
            phi = std::min(0.0f, dp);
            return true;
        }
        return false;
    }

    // Point-SDF collision test.
    //
    static inline bool collisionDetectPointSDF(const Discregrid::CubicLagrangeDiscreteGrid& sdf, const Eigen::Vector3f& p, Eigen::Vector3f& n, float& phi)
    {
        const Eigen::Vector3d pdouble  = p.cast<double>();

        if( !sdf.domain().contains(pdouble) ) return false;

        // Interpolate to find the signed distance at p,
        // as well as the gradient.
        //
        Eigen::Vector3d grad;
        const double sdist = sdf.interpolate(0, pdouble, &grad);

        if ( sdist < 1e-5 )
        {
            n = grad.cast<float>();
            phi = (float) sdist;
            return true;
        }

        return false;
    }

    // SDF-Plane collision test.
    //
    static inline void collisionDetectSdfPlane(RigidBody* body0, RigidBody* body1, std::vector<Contact*>& contacts)
    {
        SDFGeometry* geomSDF = dynamic_cast<SDFGeometry*>(body0->geometry.get());
        Plane* plane = dynamic_cast<Plane*>(body1->geometry.get());

        float phi;
        // Check : body0 vertices vs. body 1 SDF
        for(int i = 0; i < geomSDF->mesh->nVertices(); ++i)
        {
            const Eigen::Vector3d& vert = geomSDF->mesh->vertex(i);
            const Eigen::Vector3f p0 = body0->R * vert.cast<float>() + body0->x;
            const Eigen::Vector3f pplane = body1->x;
            const Eigen::Vector3f nplane = body1->R * plane->normal;

            if( collisionDetectPointPlane(p0, pplane, nplane, phi) )
            {
                contacts.push_back( new Contact(body0, body1, p0, nplane, phi) );
            }
        }

    }
}

CollisionDetect::CollisionDetect(RigidBodySystem* rigidBodySystem) : m_rigidBodySystem(rigidBodySystem)
{

}

void CollisionDetect::detectCollisions()
{
    // First, clear any existing contacts.
    //
    clear();

    // Next, loop over all pairs of bodies and test for contacts.
    //
    auto bodies = m_rigidBodySystem->getBodies();
    for(unsigned int i = 0; i < bodies.size(); ++i)
    {
        for(unsigned int j = i+1; j < bodies.size(); ++j)
        {
            RigidBody* body0 = bodies[i];
            RigidBody* body1 = bodies[j];

            // Special case: skip tests for pairs of static bodies.
            //
            if (body0->fixed && body1->fixed) 
                continue;

            // Test for sphere-sphere collision.
            if( body0->geometry->getType() == kSphere &&
                body1->geometry->getType() == kSphere )
            {
                collisionDetectSphereSphere(body0, body1);
            }
            // Test for sphere-box collision
            else if( body0->geometry->getType() == kSphere &&
                     body1->geometry->getType() == kBox )
            {
                collisionDetectSphereBox(body0, body1);
            }
            // Test for box-sphere collision (order swap)
            else if( body1->geometry->getType() == kSphere &&
                     body0->geometry->getType() == kBox )
            {
                collisionDetectSphereBox(body1, body0);
            }
            // Test for plane-box collision
            else if( body1->geometry->getType() == kPlane &&
                     body0->geometry->getType() == kBox )
            {                
                if (m_rigidBodySystem->getSamplingType() == kCorners)
                {
                    collisionDetectBoxPlane(body0, body1);
                }
                else if (m_rigidBodySystem->getSamplingType() == kGrid)
                {
                    collisionDetectBoxPlaneGrid(body0, body1);
                }
            }
            // Test for plane-box collision (order swap)
            else if ( body0->geometry->getType() == kPlane &&
                      body1->geometry->getType() == kBox )
            {
                if (m_rigidBodySystem->getSamplingType() == kCorners)
                {
                    collisionDetectBoxPlane(body1, body0);
                }
                else if (m_rigidBodySystem->getSamplingType() == kGrid)
                {
                    collisionDetectBoxPlaneGrid(body1, body0);
                }
            }
            // Test for SDF-box collision
            else if ( body0->geometry->getType() == kSDF &&
                      body1->geometry->getType() == kBox )
            {
                collisionDetectBoxSdf(body1, body0);
            }
            // Test for SDF-box collision (order swap)
            else if ( body1->geometry->getType() == kSDF &&
                      body0->geometry->getType() == kBox )
            {
                collisionDetectBoxSdf(body0, body1);
            }
            // Test of SDF-plane collision
            else if( body0->geometry->getType() == kSDF &&
                     body1->geometry->getType() == kPlane )
            {
                collisionDetectSdfPlane(body0, body1, m_contacts);
            }
            // Test of SDF-plane collision (order swap)
            else if( body0->geometry->getType() == kPlane &&
                     body1->geometry->getType() == kSDF )
            {
                collisionDetectSdfPlane(body1, body0, m_contacts);
            }
            // Test for SDF-SDF collision
            else if( body0->geometry->getType() == kSDF &&
                     body1->geometry->getType() == kSDF )
            {
                collisionDetectSdfSdf(body0, body1);
            }

        }
    }
}

void CollisionDetect::computeContactJacobians()
{
    for(auto c : m_contacts)
    {
        c->computeContactFrame( Eigen::Vector3f(1,0,0) );
        c->computeJacobian();
    }
}

void CollisionDetect::clear()
{
    for(auto c : m_contacts)
    {
        delete c;
    }
    m_contacts.clear();

    auto bodies = m_rigidBodySystem->getBodies();
    for(auto b : bodies)
    {
        b->contacts.clear();
    }
}

void CollisionDetect::collisionDetectSphereSphere(RigidBody* body0, RigidBody* body1)
{
    Sphere* sphere0 = dynamic_cast<Sphere*>(body0->geometry.get());
    Sphere* sphere1 = dynamic_cast<Sphere*>(body1->geometry.get());

    // Implement sphere-sphere collision detection.
    // The function should check if a collision exists, and if it does
    // compute the contact normal, contact point, and penetration depth.
    //
    Eigen::Vector3f vec = body0->x - body1->x;

    const float rsum = (sphere0->radius + sphere1->radius);
    const float dist = vec.norm();
    if( dist < rsum )
    {
        const Eigen::Vector3f n = vec / dist;
        const Eigen::Vector3f p = 0.5f * ((body0->x - sphere0->radius*n) + (body1->x + sphere1->radius*n));
        const float phi = dist-rsum;

        m_contacts.push_back( new Contact(body0, body1, p, n, phi) );
    }
}

void CollisionDetect::collisionDetectSphereBox(RigidBody* body0, RigidBody* body1)
{
    Sphere* sphere = dynamic_cast<Sphere*>(body0->geometry.get());
    Box* box = dynamic_cast<Box*>(body1->geometry.get());

    const Eigen::Vector3f clocal = body1->R.transpose() * (body0->x - body1->x);

    Eigen::Vector3f q(0,0,0);
    for(unsigned int i = 0; i < 3; ++i)
    {
        q[i] = std::max(-box->dim[i]/2.0f, std::min(box->dim[i]/2.0f, clocal[i]));
    }

    const Eigen::Vector3f dx = clocal - q;
    const float dist = dx.norm();
    if( dist < sphere->radius )
    {
        const Eigen::Vector3f n = body1->R * (dx/dist);
        const Eigen::Vector3f p = body1->R * q + body1->x;
        const float phi = dist - sphere->radius;

        m_contacts.push_back( new Contact(body0, body1, p, n, phi) );
    }
}

void CollisionDetect::collisionDetectBoxPlane(RigidBody* body0, RigidBody* body1)
{
    Box* box = dynamic_cast<Box*>(body0->geometry.get());
    Plane* plane = dynamic_cast<Plane*>(body1->geometry.get());
    const Eigen::Vector3f pplane = body1->x;
    const Eigen::Vector3f nplane = body1->R * plane->normal;
    const Eigen::Vector3f plocal[8] = {
        0.5f*Eigen::Vector3f(-box->dim(0), -box->dim(1), -box->dim(2)),
        0.5f*Eigen::Vector3f(-box->dim(0), -box->dim(1),  box->dim(2)),
        0.5f*Eigen::Vector3f(-box->dim(0),  box->dim(1), -box->dim(2)),
        0.5f*Eigen::Vector3f(-box->dim(0),  box->dim(1),  box->dim(2)),
        0.5f*Eigen::Vector3f( box->dim(0), -box->dim(1), -box->dim(2)),
        0.5f*Eigen::Vector3f( box->dim(0), -box->dim(1),  box->dim(2)),
        0.5f*Eigen::Vector3f( box->dim(0),  box->dim(1), -box->dim(2)),
        0.5f*Eigen::Vector3f( box->dim(0),  box->dim(1),  box->dim(2))
    };

    // generate contacts at the corners of the box
    for (unsigned int i = 0; i < 8; ++i)
    {
        const Eigen::Vector3f pbox = body0->R * plocal[i] + body0->x;
        float phi;
        if  ( collisionDetectPointPlane(pbox, pplane, nplane, phi) )
        {
            m_contacts.push_back( new Contact(body0, body1, pbox, nplane, phi) );
        }
    }
}

void CollisionDetect::collisionDetectBoxPlaneGrid(RigidBody* body0, RigidBody* body1)
{
    Box* box = dynamic_cast<Box*>(body0->geometry.get());
    Plane* plane = dynamic_cast<Plane*>(body1->geometry.get());
    const Eigen::Vector3f pplane = body1->x;
    const Eigen::Vector3f nplane = body1->R * plane->normal;
    const unsigned int num = 4;
    const float dist[3] = { box->dim(0) / (num - 1), box->dim(1) / (num - 1),  box->dim(2) / (num - 1) };
    Eigen::Vector3f plocal[num][num][num];
    for (int i = 0; i < num; ++i)
    {
        const float startX = -0.5f * box->dim(0);
        for (int j = 0; j < num; ++j)
        {
            const float startY = -0.5f * box->dim(1);
            for (int k = 0; k < num; ++k)
            {
                const float startZ = -0.5f * box->dim(2);
                plocal[i][j][k] = Eigen::Vector3f(startX + i * dist[0], startY + j * dist[1] , startZ + k * dist[2]);
            }
        }
    }

    // generate contacts at the corners of the box
    unsigned int count = 0;
    for (int i = 0; i < num; ++i)
    {
        for (int j = 0; j < num; ++j)
        {
            for (int k = 0; k < num; ++k)
            {
                const Eigen::Vector3f pbox = body0->R * plocal[i][j][k] + body0->x;
                float phi;
                if (collisionDetectPointPlane(pbox, pplane, nplane, phi))
                {
                    m_contacts.push_back(new Contact(body0, body1, pbox, nplane, phi));
                    ++count;
                }
            }
        }
    }
}

void CollisionDetect::collisionDetectBoxSdf(RigidBody* body0, RigidBody* body1)
{
    SDFGeometry* geomSDF = dynamic_cast<SDFGeometry*>(body1->geometry.get());
    Box* box = dynamic_cast<Box*>(body0->geometry.get());
    const Eigen::Vector3f plocal[8] = {
        0.5f*Eigen::Vector3f(-box->dim(0), -box->dim(1), -box->dim(2)),
        0.5f*Eigen::Vector3f(-box->dim(0), -box->dim(1),  box->dim(2)),
        0.5f*Eigen::Vector3f(-box->dim(0),  box->dim(1), -box->dim(2)),
        0.5f*Eigen::Vector3f(-box->dim(0),  box->dim(1),  box->dim(2)),
        0.5f*Eigen::Vector3f( box->dim(0), -box->dim(1), -box->dim(2)),
        0.5f*Eigen::Vector3f( box->dim(0), -box->dim(1),  box->dim(2)),
        0.5f*Eigen::Vector3f( box->dim(0),  box->dim(1), -box->dim(2)),
        0.5f*Eigen::Vector3f( box->dim(0),  box->dim(1),  box->dim(2))
    };

    for (int i = 0; i < 8; ++i)
    {
        // Compute box point in world coordinates
        const Eigen::Vector3f pbox = body0->R * plocal[i] + body0->x;

        // Transform the box point into local SDF coordinate
        const Eigen::Vector3f psdf = body1->R.transpose() * (pbox - body1->x);

        Eigen::Vector3f n;
        float phi;
        if( collisionDetectPointSDF(*(geomSDF->sdf), psdf, n, phi) )
        {
            Eigen::Vector3f nworld = body1->R * n;
            nworld.normalize();
            m_contacts.push_back( new Contact(body0, body1, pbox, nworld, phi) );
        }
    }
}

void CollisionDetect::collisionDetectSdfSdf(RigidBody* body0, RigidBody* body1)
{
    SDFGeometry* geomSDF0 = dynamic_cast<SDFGeometry*>(body0->geometry.get());
    SDFGeometry* geomSDF1 = dynamic_cast<SDFGeometry*>(body1->geometry.get());

    Eigen::Vector3f n;
    float phi;
    // Check : body0 vertices vs. body 1 SDF
    for(unsigned int i = 0; i < geomSDF0->mesh->nVertices(); ++i)
    {
        const Eigen::Vector3d& vert0 = geomSDF0->mesh->vertex(i);
        const Eigen::Vector3f p0 = body0->R * vert0.cast<float>() + body0->x;
        const Eigen::Vector3f p0_1 = body1->R.transpose() * (p0 - body1->x);

        if( collisionDetectPointSDF(*(geomSDF1->sdf), p0_1, n, phi) )
        {
            Eigen::Vector3f nworld = body1->R * n;
            nworld.normalize();
            m_contacts.push_back( new Contact(body0, body1, p0, nworld, phi) );
        }
    }
    // Check : body1 vertices vs. body 0 SDF
    for(unsigned int i = 0; i < geomSDF1->mesh->nVertices(); ++i)
    {
        const Eigen::Vector3d& vert1 = geomSDF1->mesh->vertex(i);
        const Eigen::Vector3f p1 = body1->R * vert1.cast<float>() + body1->x;
        const Eigen::Vector3f p1_0 = body0->R.transpose() * (p1 - body0->x);

        if( collisionDetectPointSDF(*(geomSDF0->sdf), p1_0, n, phi) )
        {
            Eigen::Vector3f nworld = body0->R * n;
            nworld.normalize();
            m_contacts.push_back( new Contact(body1, body0, p1, nworld, phi) );
        }
    }
}
