cmake_minimum_required(VERSION 3.12)
project(contactFrictionSim)

set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(OpenGL REQUIRED)
find_package(OpenMP REQUIRED)
set(QTDIR "C:/Qt/5.15.1/msvc2019_64")
set(CMAKE_PREFIX_PATH "C:/Qt/5.15.1/msvc2019_64/")
set(CMAKE_LIBRARY_PATH "C:/Qt/5.15.1/msvc2019_64/lib")
find_package(Qt5 REQUIRED COMPONENTS Core Widgets Xml)

set(qglviewer_HEADERS 3rdParty/QGLViewer/camera.h
            3rdParty/QGLViewer/config.h
            3rdParty/QGLViewer/constraint.h
            3rdParty/QGLViewer/domUtils.h
            3rdParty/QGLViewer/frame.h
            3rdParty/QGLViewer/keyFrameInterpolator.h
            3rdParty/QGLViewer/manipulatedCameraFrame.h
            3rdParty/QGLViewer/manipulatedFrame.h
            3rdParty/QGLViewer/mouseGrabber.h
            3rdParty/QGLViewer/qglviewer.h
            3rdParty/QGLViewer/quaternion.h
            3rdParty/QGLViewer/vec.h )
set(qglviewer_SOURCE 3rdParty/QGLViewer/camera.cpp
        3rdParty/QGLViewer/constraint.cpp
        3rdParty/QGLViewer/frame.cpp
        3rdParty/QGLViewer/keyFrameInterpolator.cpp
        3rdParty/QGLViewer/manipulatedCameraFrame.cpp
        3rdParty/QGLViewer/manipulatedFrame.cpp
        3rdParty/QGLViewer/mouseGrabber.cpp
        3rdParty/QGLViewer/qglviewer.cpp
        3rdParty/QGLViewer/quaternion.cpp
        3rdParty/QGLViewer/saveSnapshot.cpp
        3rdParty/QGLViewer/vec.cpp )
set(qglviewer_FORMS 3rdParty/QGLViewer/ImageInterface.ui 
         3rdParty/QGLViewer/VRenderInterface.ui )
        
set(discregrid_SOURCE 3rdParty/Discregrid/src/acceleration/bounding_sphere_hierarchy.cpp 
           3rdParty/Discregrid/src/cubic_lagrange_discrete_grid.cpp 
           3rdParty/Discregrid/src/discrete_grid.cpp 
           3rdParty/Discregrid/src/geometry/mesh_distance.cpp 
           3rdParty/Discregrid/src/geometry/point_triangle_distance.cpp 
           3rdParty/Discregrid/src/mesh/entity_containers.cpp 
           3rdParty/Discregrid/src/mesh/entity_iterators.cpp 
           3rdParty/Discregrid/src/mesh/triangle_mesh.cpp 
           3rdParty/Discregrid/src/utility/timing.cpp )
set(discregrid_HEADERS 3rdParty/Discregrid/include/Discregrid/acceleration/bounding_sphere.hpp 
    3rdParty/Discregrid/include/Discregrid/acceleration/bounding_sphere_hierarchy.hpp 
    3rdParty/Discregrid/include/Discregrid/acceleration/kd_tree.hpp 
    3rdParty/Discregrid/include/Discregrid/acceleration/kd_tree.inl 
    3rdParty/Discregrid/include/Discregrid/cubic_lagrange_discrete_grid.hpp 
    3rdParty/Discregrid/include/Discregrid/discrete_grid.hpp 
    3rdParty/Discregrid/include/Discregrid/geometry/mesh_distance.hpp 
    3rdParty/Discregrid/include/Discregrid/mesh/entity_containers.hpp 
    3rdParty/Discregrid/include/Discregrid/mesh/entity_iterators.hpp 
    3rdParty/Discregrid/include/Discregrid/mesh/halfedge.hpp 
    3rdParty/Discregrid/include/Discregrid/mesh/triangle_mesh.hpp 
    3rdParty/Discregrid/include/Discregrid/utility/lru_cache.hpp 
    3rdParty/Discregrid/include/Discregrid/utility/serialize.hpp 
    3rdParty/Discregrid/src/data/z_sort_table.hpp 
    3rdParty/Discregrid/src/geometry/point_triangle_distance.hpp 
    3rdParty/Discregrid/src/utility/spinlock.hpp 
    3rdParty/Discregrid/src/utility/timing.hpp )   
        
set(collision_HEADERS include/collision/CollisionDetect.h 
                      include/collision/Geometry.h )
set(collision_SOURCE src/collision/CollisionDetect.cpp ) 

set(rigidbody_HEADERS include/rigidbody/RigidBody.h
                      include/rigidbody/RigidBodySystem.h 
                      include/rigidbody/Scenarios.h )
set(rigidbody_SOURCE src/rigidbody/RigidBody.cpp 
                     src/rigidbody/RigidBodySystem.cpp )

set(contact_HEADERS include/contact/Contact.h )
set(contact_SOURCE src/contact/Contact.cpp )

set(solver_HEADERS include/solvers/SolverBoxPGS.h 
                   include/solvers/SolverBoxBPP.h
                   include/solvers/SolverPROX.h )
set(solver_SOURCE src/solvers/SolverBoxPGS.cpp 
                   src/solvers/SolverBoxBPP.cpp
                   src/solvers/SolverPROX.cpp )

set(util_HEADERS include/util/MeshAssets.h
                 include/util/MeshUtils.h
                 include/util/OBJLoader.h )
set(util_SOURCE src/util/MeshAssets.cpp 
                src/util/OBJLoader.cpp ) 

set(viewer_HEADERS include/viewer/SimViewer.h
    include/viewer/RigidBodyRenderer.h
    include/viewer/ShaderVars.h)
set(viewer_SOURCE src/viewer/SimViewer.cpp
    src/viewer/RigidBodyRenderer.cpp )
set(viewer_FORMS mainwindow.ui )
         
add_executable (contactFrictionSim main.cpp MainWindow.h MainWindow.cpp
          ${contact_HEADERS} ${contact_SOURCE} 
          ${collision_HEADERS} ${collision_SOURCE}
          ${rigidbody_HEADERS} ${rigidbody_SOURCE}
          ${solver_HEADERS} ${solver_SOURCE}
          ${viewer_HEADERS} ${viewer_SOURCE} ${viewer_FORMS}
          ${util_HEADERS} ${util_SOURCE}
          ${discregrid_HEADERS} ${discregrid_SOURCE}
          ${qglviewer_HEADERS} ${qglviewer_SOURCE} ${qglviewer_FORMS} )

add_compile_definitions(QGLVIEWER_STATIC)
target_link_libraries(contactFrictionSim OpenGL::GL OpenMP::OpenMP_CXX Qt5::Core Qt5::Widgets Qt5::Xml)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include 3rdParty/QGLViewer 3rdParty/Eigen3/include/eigen3 3rdParty/Discregrid/include ${COMMON_INCLUDES})

# Setup source code groups
#
source_group(src FILES main.cpp MainWindow.cpp)
source_group(include FILES MainWindow.h)
source_group(ui FILES main.ui)
source_group(3rdParty/qglviewer/src FILES ${qglviewer_SOURCE})
source_group(3rdParty/qglviewer/include FILES ${qglviewer_HEADERS})
source_group(3rdParty/discregrid/src FILES ${discregrid_SOURCE})
source_group(3rdParty/discregrid/include FILES ${discregrid_HEADERS})
source_group(src/collision FILES ${collision_SOURCE})
source_group(include/collision FILES ${collision_HEADERS})
source_group(src/contact FILES ${contact_SOURCE})
source_group(include/contact FILES ${contact_HEADERS})
source_group(src/rigidbody FILES ${rigidbody_SOURCE})
source_group(include/rigidbody FILES ${rigidbody_HEADERS})
source_group(src/solvers FILES ${solver_SOURCE})
source_group(include/solvers FILES ${solver_HEADERS})
source_group(src/viewer FILES ${viewer_SOURCE})
source_group(include/viewer FILES ${viewer_HEADERS})
source_group(src/util FILES ${util_SOURCE})
source_group(include/util FILES ${util_HEADERS})

if(MSVC) 
	set_property(TARGET contactFrictionSim PROPERTY VS_DEBUGGER_WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
	
	# Discregrid compilation flags
	#
  	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /openmp")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4250")
endif()