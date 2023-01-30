cmake_minimum_required(VERSION 3.14.0)
project(ccdbgupdater VERSION 0.1.0)
set(CONTACT "khodor.hannoush@inria.fr")

#set to default type
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()
#adding the executable
add_executable(${PROJECT_NAME} main.cpp	ParseGFA.cpp unitigmodification.cpp  CommonUtils.cpp index_kmer.cpp)

#adding the directory
add_subdirectory(src)
add_subdirectory(external/sparsehash)

#linking to directories
target_link_directories(${PROJECT_NAME}
    PRIVATE src
    PRIVATE external/sparsehash/src
)

#linking to libraries
target_link_libraries(${PROJECT_NAME} src sparsehash)