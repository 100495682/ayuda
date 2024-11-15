add_library(imgaos
        imgaos.cpp
        imgaos.hpp
        ../utest-imgaos/test_imgaos.cpp)
target_link_libraries(imgaos PUBLIC common PRIVATE Microsoft.GSL::GSL)