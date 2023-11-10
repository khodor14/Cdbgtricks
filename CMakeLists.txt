cmake_minimum_required(VERSION 3.10)
add_compile_options(-O3)
add_compile_options(-march=native)
add_compile_options(-funroll-loops)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0")
#set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS} -O0 -g -DNDEBUG")
#set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -g")
project(ccdbgupdater VERSION 0.1.0)
set(CONTACT "khodor.hannoush@inria.fr")
#set to default type

#adding the executable
add_executable(${PROJECT_NAME} src/main.cpp	src/ParseGFA.cpp src/unitigmodification.cpp  src/CommonUtils.cpp src/FileSerializer.hpp src/index_kmer.cpp src/unitig.cpp src/mapper.cpp src/zstr.hpp src/strict_fstream.hpp src/kseq.h src/index.cpp)

#adding the directory
add_subdirectory(external/sparsehash)
target_include_directories(${PROJECT_NAME}
    #PUBLIC src
    PUBLIC external/sparsehash
    PUBLIC external/boost
    PUBLIC external/unordered_dense
)
#linking to directories
target_link_directories(${PROJECT_NAME}
    #PRIVATE src
    PRIVATE external/sparsehash
    PRIVATE external/boost
    PRIVATE external/unordered_dense
)
find_package(ZLIB REQUIRED)
target_link_libraries(${PROJECT_NAME} ${ZLIB_LIBRARIES})

if (ZLIB_FOUND)
    include_directories(${ZLIB_INCLUDE_DIRS})
else()
    message(FATAL_ERROR "zlib not found. Required for to output files")
endif(ZLIB_FOUND)
find_package(zstd REQUIRED)
if (zstd_FOUND)
    target_link_libraries(${PROJECT_NAME} zstd)
    include_directories(${zstd_INCLUDE_DIRS})
else()
    message(FATAL_ERROR "Zstandard (zstd) not found.")
endif(zstd_FOUND)
target_link_libraries(${PROJECT_NAME} boost_serialization)
target_link_libraries(${PROJECT_NAME} pthread)
find_package(Boost REQUIRED COMPONENTS serialization)
add_library(boost::serialization INTERFACE IMPORTED)
set_target_properties(boost::serialization PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${Boost_INCLUDE_DIRS}"
    INTERFACE_LINK_LIBRARIES "${Boost_SERIALIZATION_LIBRARIES}")