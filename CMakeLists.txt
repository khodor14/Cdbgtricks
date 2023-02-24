cmake_minimum_required(VERSION 3.10)
add_compile_options(-O3)
add_compile_options(-march=native)
add_compile_options(-funroll-loops)
project(ccdbgupdater VERSION 0.1.0)
set(clear "khodor.hannoush@inria.fr")
#set to default type
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()
#adding the executable
add_executable(${PROJECT_NAME} src/main.cpp	src/ParseGFA.cpp src/unitigmodification.cpp  src/CommonUtils.cpp src/index_kmer.cpp)

#adding the directory
#add_subdirectory(src)
add_subdirectory(external/sparsehash)
target_include_directories(${PROJECT_NAME}
    #PUBLIC src
    PUBLIC external/sparsehash
    PUBLIC external/boost
)
#linking to directories
target_link_directories(${PROJECT_NAME}
    #PRIVATE src
    PRIVATE external/sparsehash
    PRIVATE external/boost
)

#linking to libraries
target_link_libraries(${PROJECT_NAME} boost_serialization)
