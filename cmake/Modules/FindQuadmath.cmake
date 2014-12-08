# Module that checks whether the compiler supports the
# quadruple precision floating point math
#
# Sets the following variables:
# HAVE_QUAD
# QUADMATH_LIBRARIES
#
# perform tests
include(CheckCSourceCompiles)
include(CheckCXXSourceCompiles)
include(CMakePushCheckState)
include(CheckCXXCompilerFlag)

check_cxx_compiler_flag("-Werror -fext-numeric-literals" HAVE_EXTENDED_NUMERIC_LITERALS)
if (HAVE_EXTENDED_NUMERIC_LITERALS)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fext-numeric-literals")
endif()

cmake_push_check_state()
list(APPEND CMAKE_REQUIRED_LIBRARIES "quadmath")
CHECK_CXX_SOURCE_COMPILES("
#include <quadmath.h>

int main(void){
    __float128 foo = sqrtq(123.456);
    foo = FLT128_MIN;
}" HAVE_QUAD)
cmake_pop_check_state()

if (HAVE_QUAD)
  set(QUADMATH_LIBRARIES "quadmath")
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(QuadMath
  DEFAULT_MSG
  QUADMATH_LIBRARIES
  HAVE_QUAD
  )
