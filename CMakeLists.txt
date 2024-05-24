cmake_minimum_required(VERSION 3.10)
add_compile_options(-O3)
add_compile_options(-march=native)
add_compile_options(-funroll-loops)
add_compile_options(-std=c++17)
project(cdbgtricks VERSION 0.1.1)
set(CONTACT "khodor.hannoush@inria.fr")
#set to default type
#adding the executable
add_executable(${PROJECT_NAME} src/main.cpp src/ParseGFA.cpp src/unitigmodification.cpp  src/QueryStreamer.cpp src/CommonUtils.cpp src/unitig.cpp src/mapper.cpp src/zstr.hpp src/strict_fstream.hpp src/kseq.h src/index.cpp)

#adding the directory
target_compile_definitions(${PROJECT_NAME} PRIVATE UNIT_TEST=1)
add_subdirectory(external/googletests)
target_include_directories(${PROJECT_NAME}
    #PUBLIC src
    PUBLIC external/unordered_dense
    PUBLIC external/pthash
    PUBLIC external/pthash/include
    PUBLIC external/pthash/external/xxHash
    PUBLIC external/pthash/external/essentials/include
)
#linking to directories
target_link_directories(${PROJECT_NAME}
    #PRIVATE src
    PRIVATE external/unordered_dense
    PRIVATE external/pthash
    PRIVATE external/pthash/include
    PRIVATE external/pthash/external/xxHash
    PRIVATE external/pthash/external/essentials/include
)
find_package(ZLIB REQUIRED)
target_link_libraries(${PROJECT_NAME} ${ZLIB_LIBRARIES})

if (ZLIB_FOUND)
    include_directories(${ZLIB_INCLUDE_DIRS})
else()
    message(FATAL_ERROR "zlib not found. Required for to output files")
endif(ZLIB_FOUND)

if (CompileTest)
    add_subdirectory(test)
endif()

# phthread
target_link_libraries(${PROJECT_NAME} pthread)
