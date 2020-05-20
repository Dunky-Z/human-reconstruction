
cmake_minimum_required (VERSION 3.8)

project ("BodyReconstruction")

#PMP
find_package(PMP)
#include_directories("D:/Library_Includes/glew/glew-2.1.0/include")
#include_directories("D:/Library_Includes/OpenGL/include")

#Eigen
find_package(Eigen3 REQUIRED)

#SURFACE MESH
include_directories("$ENV{Surface_Mesh_DIR}/include")

include_directories("./include")
LINK_DIRECTORIES("./Lib")
#include_directories("D:\Lib\OpenGL\include")

#find_package(GLFW REQUIRED)

add_executable(3D-HumanBodyshape src/testpmp.cpp src/binary_io.h src/main.cpp src/mesh_io.cpp src/mesh_io.h src/preprocessing.cpp src/preprocessing.h src/utils.h src/measure.cpp src/measure.h src/tri_mesh.h)
target_link_libraries(3D-HumanBodyshape pmp)
target_link_libraries(3D-HumanBodyshape Eigen3::Eigen)
target_link_libraries(3D-HumanBodyshape debug "$ENV{Surface_Mesh_DIR}/lib/surface_mesh_debug.lib")
target_link_libraries(3D-HumanBodyshape debug "D:/ITabc/ITabc/BodyReconstruction/BodyReconstruction/libs/python37_d.lib")

#target_link_libraries(3D-BodyShape debug "D:/Library_Includes/glew/glew-2.1.0/lib/Release/x64/glew32.lib")
#target_link_libraries(3D-BodyShape debug "D:/Library_Includes/OpenGL/Libs/glfw3.lib")
#target_link_libraries(3D-BodyShape debug "D:/Library_Includes/OpenGL/Libs/glut.lib")
#target_link_libraries(3D-BodyShape debug "D:/Library_Includes/OpenGL/Libs/glut32.lib")


#target_link_libraries(3D-BodyShape GLFW)