# Executable for all unit tests with list of sources
# For example, you may have one *_test.cpp for each *.cpp in img-aos
add_executable(utest-img-soa
        one_test.cpp
        other1_test.cpp
        other2_test.cpp)

# Library dependencies
target_link_libraries(utest-img-soa
        PRIVATE
        imgsoa
        GTest::gtest_main
        Microsoft.GSL::GSL)
