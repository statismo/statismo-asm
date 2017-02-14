CMAKE_MINIMUM_REQUIRED(VERSION 3.1)
PROJECT(statismo-asm)

set(CMAKE_CXX_STANDARD 11) 
set(CMAKE_CXX_STANDARD_REQUIRED ON) 
set(CMAKE_CXX_EXTENSIONS OFF) 

FIND_PACKAGE(statismo REQUIRED)

FIND_PACKAGE(ITK REQUIRED)
include (${ITK_USE_FILE})



include_directories(${statismo_INCLUDE_DIRS}  ${CMAKE_CURRENT_SOURCE_DIR}/src)
link_directories( ${ITK_LIBRARY_DIRS} ${statismo_LIBRARY_DIRS})

add_executable (build-asm src/itkBuildASM.cxx )
target_link_libraries (build-asm ${ITK_LIBRARIES} ${statismo_LIBRARIES})

add_executable (fit-asm src/itkFitASM.cxx )
target_link_libraries (fit-asm ${ITK_LIBRARIES} ${statismo_LIBRARIES})

