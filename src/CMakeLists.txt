
cmake_minimum_required (VERSION 3.8)

project ("human-reconstruction")

#PMP
#find_package(PMP)

#Eigen
find_package(Eigen3 REQUIRED)
include_directories("D:/ITabc/ITabc/BodyReconstruction/BodyReconstruction/include")
LINK_DIRECTORIES("D:/ITabc/ITabc/BodyReconstruction/BodyReconstruction/Lib")
#SURFACE MESH
include_directories("D:/Lib/pmp-library-debug/include")

#find_package(GLFW REQUIRED)
FILE(GLOB SC_FILES "*.cpp" "*.h")
FILE(GLOB ALGLIB_FILES "../alglib/cpp/src/*.cpp" "../alglib/cpp/src/*.h")
FILE(GLOB GEODESIC_FILES "../geodesic/*.cpp" "../geodesic/*.h")

source_group(alglib FILES ${ALGLIB_FILES})
source_group(geodesic FILES ${GEODESIC_FILES})

add_executable(${PROJECT_NAME} ${SC_FILES} ${ALGLIB_FILES} ${GEODESIC_FILES})

target_link_libraries(human-reconstruction debug "D:/Lib/pmp-library-debug/lib/pmp.lib")
target_link_libraries(human-reconstruction optimized "D:/Lib/pmp-library-release/lib/pmp.lib")
target_link_libraries(${PROJECT_NAME} Eigen3::Eigen)
#target_link_libraries(human-reconstruction debug "$ENV{Surface_Mesh_DIR}/lib/surface_mesh_debug.lib")
#target_link_libraries(3D-HumanBodyshape debug "D:/ITabc/ITabc/BodyReconstruction/BodyReconstruction/libs/python37_d.lib")
