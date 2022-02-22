#pragma once

#include "util/MeshAssets.h"
#include "RigidBody.h"
#include "RigidBodySystem.h"

#include <Eigen/Dense>

class Scenarios
{
public:

    // Box on a plane
    //
    static void createBoxOnPlane(RigidBodySystem& rigidBodySystem)
    {
        rigidBodySystem.clear();

        RigidBody* body0 = new RigidBody(1.0f, new Plane(Eigen::Vector3f(0.0f, 1.0f, 0.0f)), "resources/plane.obj");
        body0->fixed = true;

        RigidBody* body1 = new RigidBody(1.0f, new Box(Eigen::Vector3f(1.0f, 1.0f, 1.0f)), "resources/box.obj");
        body1->x = Eigen::Vector3f(0.0f, 0.49f, 0.0f);

        RigidBody* body2 = new RigidBody(10.0f, new Box(Eigen::Vector3f(1.0f, 1.0f, 1.0f)), "resources/box.obj");
        body2->x = Eigen::Vector3f(2.0f, 0.49f, 0.0f);
        body2->xdot = Eigen::Vector3f(0.0f, 0.0f, 10.0f);

        rigidBodySystem.addBody(body0);
        rigidBodySystem.addBody(body1);
        rigidBodySystem.addBody(body2);
    }

    // Box on a plane
    //
    static void createBoxOnInclinedPlane(RigidBodySystem& rigidBodySystem)
    {
        rigidBodySystem.clear();

        auto angle = M_PI * 22 / 180;
        auto sinA = std::sin(angle / 2);
        auto cosA = std::cos(angle / 2);
        printf("sinA %f", sinA);
        printf("cosA %f", cosA);

        Eigen::Quaternionf q;
        q.x() = 0 * sinA;
        q.y() = 0 * sinA;
        q.z() = 1 * sinA;
        q.w() = cosA;
        q.normalize();

        RigidBody* body0 = new RigidBody(1.0f, new Plane(Eigen::Vector3f(0.0f, 1.0f, 0.0f)), "resources/plane.obj");
        body0->fixed = true; 
        body0->q = q;

        RigidBody* body1 = new RigidBody(1.0f, new Box(Eigen::Vector3f(1.0f, 1.0f, 1.0f)), "resources/box.obj");
        body1->q = q;
        body1->x = Eigen::Vector3f(0.0f, 0.52f, 0.0f);

        rigidBodySystem.addBody(body0);
        rigidBodySystem.addBody(body1);
    }

    // Stack of boxes and spheres.
    //
    static void createBoxBallStack(RigidBodySystem& rigidBodySystem)
    {
        rigidBodySystem.clear();

        RigidBody* body0 = new RigidBody(1.0f, new Plane(Eigen::Vector3f(0.0f, 1.0f, 0.0f)), "resources/plane.obj");
        body0->fixed = true;

        RigidBody* body1 = new RigidBody(1.0f, new Box(Eigen::Vector3f(1.0f, 1.0f, 1.0f)), "resources/box.obj");
        body1->x = Eigen::Vector3f(0.0f, 0.5f, 0.0f);
        RigidBody* body2 = new RigidBody(1.0f, new Sphere(0.5f), "resources/sphere.obj");
        body2->x = Eigen::Vector3f(0.0f, 1.5f, 0.0f);
        RigidBody* body3 = new RigidBody(1.0f, new Box(Eigen::Vector3f(1.0f, 1.0f, 1.0f)), "resources/box.obj");
        body3->x = Eigen::Vector3f(0.0f, 2.5f, 0.0f);
        RigidBody* body4 = new RigidBody(1.0f, new Sphere(0.5f), "resources/sphere.obj");
        body4->x = Eigen::Vector3f(0.0f, 3.5f, 0.0f);

        rigidBodySystem.addBody(body0);
        rigidBodySystem.addBody(body1);
        rigidBodySystem.addBody(body2);
        rigidBodySystem.addBody(body3);
        rigidBodySystem.addBody(body4);
    }

    // Box filled with balls.
    //
    static void createMarbleBox(RigidBodySystem& rigidBodySystem)
    {
        rigidBodySystem.clear();


        // Create two layers of "marbles", in a grid layout
        //
        for(int i = 0; i < 9; ++i)
        {
            for(int j = 0; j < 9; ++j)
            {
                RigidBody* body1 = new RigidBody(1.0f, new Sphere(0.5f), "resources/sphere.obj");
                body1->x.x() = -4.0f + (float)i*1.0f;
                body1->x.z() = -4.0f + (float)j*1.0f;
                body1->x.y() = 2.0f;
                rigidBodySystem.addBody(body1);
                RigidBody* body2 = new RigidBody(1.0f, new Sphere(0.5f), "resources/sphere.obj");
                body2->x.x() = -4.0f + (float)i*1.0f;
                body2->x.z() = -4.0f + (float)j*1.0f;
                body2->x.y() = 3.0f;
                rigidBodySystem.addBody(body2);
            }
        }

        RigidBody* body0 = new RigidBody(1.0f, new Box(Eigen::Vector3f(0.4f, 4.0f, 10.0f)), "resources/box_side.obj");
        RigidBody* body1 = new RigidBody(1.0f, new Box(Eigen::Vector3f(0.4f, 4.0f, 10.0f)), "resources/box_side.obj");
        RigidBody* body2 = new RigidBody(1.0f, new Box(Eigen::Vector3f(0.4f, 4.0f, 10.0f)), "resources/box_side.obj");
        RigidBody* body3 = new RigidBody(1.0f, new Box(Eigen::Vector3f(0.4f, 4.0f, 10.4f)), "resources/box_side.obj");
        RigidBody* body4 = new RigidBody(1.0f, new Box(Eigen::Vector3f(10.0f, 0.4f, 10.0f)), "resources/box_bot.obj");
        body0->fixed = true;
        body1->fixed = true;
        body2->fixed = true;
        body3->fixed = true;
        body4->fixed = true;
        body0->x.x() = 5.0f;
        body1->x.x() = -5.0f;
        body2->x.z() = 5.0f;
        body2->q = Eigen::AngleAxisf(1.57, Eigen::Vector3f(0, 1, 0));
        body3->x.z() = -5.0f;
        body3->q = Eigen::AngleAxisf(1.57, Eigen::Vector3f(0, 1, 0));
        body4->x.y() = -0.5f;

        rigidBodySystem.addBody(body0);
        rigidBodySystem.addBody(body1);
        rigidBodySystem.addBody(body2);
        rigidBodySystem.addBody(body3);
        rigidBodySystem.addBody(body4);
    }


    // Box on a plane
    //
    static void createBunnies(RigidBodySystem& rigidBodySystem)
    {
        rigidBodySystem.clear();

        RigidBody* body0 = new RigidBody(1.0f, new Plane(Eigen::Vector3f(0.0f, 1.0f, 0.0f)), "resources/plane.obj");
        body0->fixed = true;

        RigidBody* body1 = new RigidBody(1.0f, new SDFGeometry("resources/bunny.obj", {10, 10, 10}), "resources/bunny.obj");
        RigidBody* body2 = new RigidBody(1.0f, new SDFGeometry("resources/bunny.obj", {10, 10, 10}), "resources/bunny.obj");
        RigidBody* body3 = new RigidBody(1.0f, new SDFGeometry("resources/bunny.obj", {10, 10, 10}), "resources/bunny.obj");

        body1->x = Eigen::Vector3f(0.5f, 2.0f, -0.5f);
        body2->x = Eigen::Vector3f(0.5f, 4.0f, 0.5);
        body3->x = Eigen::Vector3f(-0.5f, 6.0f, 0.5f);

        rigidBodySystem.addBody(body0);
        rigidBodySystem.addBody(body1);
        rigidBodySystem.addBody(body2);
        rigidBodySystem.addBody(body3);
    }

};
