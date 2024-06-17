cmake_minimum_required(VERSION 3.14)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

add_library(qd STATIC
    src/bits.cpp
    src/c_dd.cpp
    src/c_qd.cpp 
    src/dd_const.cpp 
    src/dd_real.cpp 
    src/fpu.cpp 
    src/qd_const.cpp 
    src/qd_real.cpp 
    src/util.cpp 
)

SET_TARGET_PROPERTIES( qd PROPERTIES COMPILE_FLAGS -fPIC)

