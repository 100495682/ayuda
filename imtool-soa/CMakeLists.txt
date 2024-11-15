add_executable(imtool-soa main.cpp)
target_link_libraries(imtool-soa PRIVATE imgsoa common Microsoft.GSL::GSL)