cmake_minimum_required (VERSION 2.6)
project (strassen)

set(CMAKE_BUILD_TYPE Debug)

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O4")

FILE (GLOB SOURCES "${PROJECT_SOURCE_DIR}/src/*.c"
	               "${PROJECT_SOURCE_DIR}/include/*.h")

include_directories("${PROJECT_SOURCE_DIR}/include")

SET(${LIB_REQUIRED_LIBRARIES} "")

if ("${CMAKE_SYSTEM_NAME}" STREQUAL "Linux")
  # If your system is a Linux system, then try to detect if librt is 
  # required, i.e.,  glibc < 2.17 || glibc = 2.29
  
  include(CheckCXXSourceCompiles)
  set(CLOCK_GETTIME_TEST
    "
    #include <time.h>
    int main() {
      clock_gettime(CLOCK_REALTIME, NULL);
      return 0;
    }
    "
  )
  check_cxx_source_compiles( "${CLOCK_GETTIME_TEST}" 
                             LIBRT_IS_NOT_REQUIRED )
  if (NOT LIBRT_IS_NOT_REQUIRED)
    # ch
    message(STATUS "libc doesn't contain clock_gettime()")
    set(CMAKE_REQUIRED_LIBS_OLD "${CMAKE_REQUIRED_LIBRARIES}")
    set(CMAKE_REQUIRED_LIBRARIES "${CMAKE_REQUIRED_LIBRARIES} rt")
	check_cxx_source_compiles( "${CLOCK_GETTIME_TEST}" 
                               LIBRT_CONTAINS_CLOCK_GETTIME)
    set(CMAKE_REQUIRED_LIBRARIES "${CMAKE_REQUIRED_LIBS_OLD}")
    if (LIBRT_CONTAINS_CLOCK_GETTIME)
      list(APPEND LIB_REQUIRED_LIBRARIES "rt")
    else()
      message(FATAL_ERROR "librt (required to get clock_gettime()) cannot be linked")
    endif()
  endif()
endif()


add_executable (strassen_test ${SOURCES})

