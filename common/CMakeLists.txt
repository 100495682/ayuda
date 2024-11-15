# Add to this list all files related to common library
file(GLOB_RECURSE COMMON_SOURCES "*.cpp" "*.hpp")

add_library(common ${COMMON_SOURCES}
        #info.cpp
        binaryo.hpp
        decompress.cpp)

# Use this line only if you have dependencies from this library to GSL
target_link_libraries(common PRIVATE Microsoft.GSL::GSL)