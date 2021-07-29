#pragma once

#include <Eigen/Dense>
#include <Discregrid/All>
#include <memory>

// List of geometry type ids.
enum eGeometryType { kSphere, kBox, kPlane, kSDF };

// Generic geometry interface.
//
class Geometry
{
public:
    virtual Eigen::Matrix3f computeInertia(float _mass) = 0;

    virtual eGeometryType getType() const  = 0;

protected:
    Eigen::Matrix3f m_I;          // Inertia 3x3 matrix for this. Only used for local computations. (internal)
};

// Signed distance field (SDF) geometry.
// Our implementation uses Discregrid.
//
class SDFGeometry : public Geometry
{
private:
    Eigen::Matrix3f I0;         // Inertia matrix, computed from vertex positions

    // Compute the inertia of the mesh.
    // Assumes a manifold mesh and center of mass is the geometry origin.
    //
    // Based on code from: http://melax.github.io/volint.html (Stan Melax)
    //
    void computeInertia()
    {
        if( mesh )
        {
            Eigen::Matrix3f A;
            float vol = 0.0f;
            Eigen::Vector3f diag = Eigen::Vector3f::Zero();
            Eigen::Vector3f offd = Eigen::Vector3f::Zero();
            for(unsigned int i = 0; i < mesh->nFaces(); ++i)
            {
                const std::array<unsigned int, 3>& face = mesh->face(i);
                A.col(0) = mesh->vertex(face[0]).cast<float>();
                A.col(1) = mesh->vertex(face[1]).cast<float>();
                A.col(2) = mesh->vertex(face[2]).cast<float>();
                const float det = A.determinant();
                vol += det;

                for(int j = 0; j < 3; ++j)
                {
                    const int j1 = (j+1) % 3;
                    const int j2 = (j+2) % 3;
                    diag(j) += (A(j,0)*A(j,1) + A(j,1)*A(j,2) + A(j,2)*A(j,0) +
                                A(j,0)*A(j,0) + A(j,1)*A(j,1) + A(j,2)*A(j,2)  ) * det; // divide by 60.0f later
                    offd(j) += (A(j1,0)*A(j2,1) + A(j1,1)*A(j2,2) + A(j1,2)*A(j2,0)  +
                                A(j1,0)*A(j2,2) + A(j1,1)*A(j2,0) + A(j1,2)*A(j2,1)  +
                                2.0f*A(j1,0)*A(j2,0) + 2.0f*A(j1,1)*A(j2,1) + 2.0f*A(j1,2)*A(j2,2) ) * det; // divide by 120.0f later
                }
            }

            diag /= vol * (60.0f / 6.0f);  // divide by total volume (vol/6) since density= 1/volume
            offd /= vol * (120.0f / 6.0f);
            I0(0,0) = diag.y()+diag.z(); I0(0,1) = -offd.z(); I0(0,2) = -offd.y();
            I0(1,0) = -offd.z(); I0(1,1) = diag.x()+diag.z(); I0(1,2) = -offd.x();
            I0(2,0) = -offd.y(); I0(2,1) = -offd.x(); I0(2,2) = diag.x()+diag.y();
        }
        else
        {
            I0.setZero();
        }
    }

public:
    std::unique_ptr<Discregrid::TriangleMesh> mesh;
    std::unique_ptr<Discregrid::CubicLagrangeDiscreteGrid> sdf;

    SDFGeometry(const std::string& filename, const std::array<unsigned int, 3>& resolution)
        : mesh(), sdf()
    {
        Eigen::AlignedBox3d domain;
        domain.setEmpty();

        mesh = std::make_unique<Discregrid::TriangleMesh>(filename);

        Discregrid::MeshDistance md(*mesh);

        // Compute the bounding box of the mesh
        //
        for (auto const& x : mesh->vertices())
        {
            domain.extend(x);
        }

        // Enlarge the bounding box slightly
        //
        domain.max() += 1e-3f * domain.diagonal().norm() * Eigen::Vector3d::Ones();
        domain.min() -= 1e-3f * domain.diagonal().norm() * Eigen::Vector3d::Ones();

        sdf = std::make_unique<Discregrid::CubicLagrangeDiscreteGrid>(domain, resolution);

        // Create function that returns the SD
        auto func = Discregrid::DiscreteGrid::ContinuousFunction{};
        func = [&md](Eigen::Vector3d const& xi) { return md.signedDistanceCached(xi); };
        sdf->addFunction(func, true);

        computeInertia();
    }

    virtual ~SDFGeometry() {}

    // The mass is the total weight of the object
    // Use system of particle assumption when computing inertia.
    // All particles weight the same.
    virtual Eigen::Matrix3f computeInertia(float _mass) override
    {
        m_I = _mass * I0;
        return m_I;
    }

    virtual eGeometryType getType() const override { return kSDF; }

};

// Sphere geometry.
//
class Sphere : public Geometry
{
public:
    float radius;           // Sphere radius.

    Sphere(float _radius) : radius(_radius) {}
    virtual ~Sphere() {}

    virtual Eigen::Matrix3f computeInertia(float _mass) override
    {
        m_I.setZero();
        m_I(0,0) = m_I(1,1) = m_I(2,2) = (2.0f/5.0f) * _mass * radius * radius;
        return m_I;
    }

    virtual eGeometryType getType() const override { return kSphere; }

};

// Box geometry.
//
class Box : public Geometry
{
public:
    Eigen::Vector3f dim;        // Box dimensions.

    Box(const Eigen::Vector3f& _dim) : dim(_dim) {

    }
    virtual ~Box() {}

    virtual Eigen::Matrix3f computeInertia(float _mass) override
    {
        m_I.setZero();
        m_I(0,0) = (1.0f/12.0f)*_mass*(dim[1]*dim[1] + dim[2]*dim[2]);
        m_I(1,1) = (1.0f/12.0f)*_mass*(dim[0]*dim[0] + dim[2]*dim[2]);
        m_I(2,2) = (1.0f/12.0f)*_mass*(dim[0]*dim[0] + dim[1]*dim[1]);
        return m_I;
    }

    virtual eGeometryType getType() const override { return kBox; }

};

// Plane geometry.
//
// The infinite plane is defined by its perpendicular direction (normal),
// and we assume that the plane is centered at the body COM.
//
class Plane : public Geometry
{
public:
    Eigen::Vector3f normal;         // The plane normal.

    Plane(const Eigen::Vector3f& _normal)
        : normal(_normal) {}
    virtual ~Plane() {}

    virtual Eigen::Matrix3f computeInertia(float _mass) override
    {
        m_I.setIdentity();
        return _mass*m_I;
    }

    virtual eGeometryType getType() const override { return kPlane; }
};
