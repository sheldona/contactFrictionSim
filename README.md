# Rigid body contact and friction simulator

![Rigid body sim](https://siggraphcontact.github.io/assets/images/rigidbodysim_sshot1.png "Rigid body sim")

A simple C++ rigid body simulator that implements a variety of frictional contact models and solver methods.

The code roughly follows the [SIGGRAPH'21 Course](https://siggraphcontact.github.io/) on contact and friction simulation for computer graphics.

CMake is used to build the main application.  To build and run the code, first create a sub-directory: 

```    
> mkdir build
```

Then, configure the project for the specific build target:

```
> cd build
> cmake -DCMAKE_BUILD_TYPE=Debug ..
```

Finally, compile using your favorite build tool or IDE, e.g.
```
> make
```

The code has been compiled and tested on Windows (MS Visual C++ 2019), Ubuntu Linux (g++ 9.3 with VS Code), and MacOS (Clang)



## Dependencies

### Qt5
The main application uses [Qt5](https://doc.qt.io/qt-5/) for rendering the UI, handling mouse and keyboard calls, and managing calls to OpenGL and the shaders. For more details about installing the toolkit on your system, please see:
[https://doc.qt.io/qt-5/gettingstarted.html](https://doc.qt.io/qt-5/gettingstarted.html) 

Setting the environment variable *QTDIR* to the path of your installation should allow CMake to find it automatically.

### OpenGL
The shaders and rendering code use [OpenGL 4](https://www.opengl.org/). Hopefully you already have this installed through your operating system or video card drivers.

### Included dependencies
The following dependencies are included and compiled automatically (see the **3rdParty** directory):

 * [Eigen](https://eigen.tuxfamily.org/)
 * [Discregrid](https://github.com/InteractiveComputerGraphics/Discregrid)
 * [qglviewer](http://libqglviewer.com/)

## Source code

### solvers/

Implementations of the following contact solvers:
 * **SolverBoxBPP.h**: Block principal pivoting (BPP) for Boxed LCP
 * **SolverBoxPGS.h**: Matrix-free Projected Gauss-Seidel (PGS) for Boxed LCP

### collision/

The **CollisionDetect.h** class handles collision detection and contact generation for all pairs of simulation bodies. A variety of shapes and geometries are supported:
 * Spheres
 * Planes
 * Boxes

Contains collision tests are implemented for the following geometry pairs:
  * Sphere-sphere
  * Sphere-box
  * Plane-sphere
  * Plane-box
  * SDF-box
  * SDF-SDF
  * Box-box (TODO)
  
### contact/

Implementations of various contact and friction cone models:
 * **ContactBLCP.h**: Boxed friction cone
