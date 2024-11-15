# Executable for all unit tests with list of sources
# For example, you may have one *_test.cpp for each *.cpp in img-aos
add_executable(utest-img-aos
)

# Library dependencies
target_link_libraries(utest-img-aos
        PRIVATE
        imgaos
        GTest::gtest_main
        Microsoft.GSL::GSL)
