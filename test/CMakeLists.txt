file(GLOB TEST_FILES "test_*.cpp")
file(GLOB CPP_FILES "../src/*.cpp")
list(FILTER CPP_FILES EXCLUDE REGEX "../src/main.cpp$")  # Exclude main.cpp
add_executable(ccdbgupdater_tests ${TEST_FILES} ${CPP_FILES} ../src/kseq.h ../src/zstr.hpp ../src/strict_fstream.hpp)

find_package(ZLIB REQUIRED)
target_link_libraries(ccdbgupdater_tests ${ZLIB_LIBRARIES})

find_package(zstd REQUIRED)
if (zstd_FOUND)
    target_link_libraries(ccdbgupdater_tests zstd)
    include_directories(${zstd_INCLUDE_DIRS})
else()
    message(FATAL_ERROR "Zstandard (zstd) not found.")
endif(zstd_FOUND)
target_link_libraries(ccdbgupdater_tests boost_serialization)
target_link_libraries(ccdbgupdater_tests pthread)
find_package(Boost REQUIRED COMPONENTS serialization)
add_library(boost::serialization INTERFACE IMPORTED)
set_target_properties(boost::serialization PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${Boost_INCLUDE_DIRS}"
    INTERFACE_LINK_LIBRARIES "${Boost_SERIALIZATION_LIBRARIES}")

target_link_libraries(ccdbgupdater_tests gtest gtest_main)