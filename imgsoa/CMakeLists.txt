add_library(imgsoa
        imgsoa.cpp
        imgsoa.hpp)
target_link_libraries(imgsoa PUBLIC common PRIVATE Microsoft.GSL::GSL)