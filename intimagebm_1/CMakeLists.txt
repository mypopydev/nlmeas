cmake_minimum_required(VERSION 2.8)

# The build type "Release" adds some optimizations
if (NOT CMAKE_BUILD_TYPE)
     set (CMAKE_BUILD_TYPE "Release")
endif ()

# SETUP A MAKE CLEAN TO CLEAN EVERYTHING
SET_DIRECTORY_PROPERTIES( PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES 
      "CMakeCache.txt;cmake_install.cmake" )

# Enable C99
if (CMAKE_VERSION VERSION_LESS "3.1")
     set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=gnu99")
else ()
     set (CMAKE_C_STANDARD 99)
endif ()

# libmath
SET( LIBS m )

# IIO
if( NOT TARGET IIOLIB )
   INCLUDE_DIRECTORIES("iio")
   ADD_SUBDIRECTORY(iio)
ENDIF() 
SET( LIBS ${LIBS} IIOLIB )


# Link LibFFTW
find_path (FFTW_INCLUDE_DIR fftw3.h)
find_library (FFTWF_LIBRARIES NAMES fftw3)
include_directories (SYSTEM ${FFTW_INCLUDE_DIR})
set (LIBS ${LIBS} ${FFTWF_LIBRARIES})
if (NOT FFTW_INCLUDE_DIR)
     message (FATAL_ERROR "FFTW3 not found.")
endif ()

# Enable OpenMP
find_package (OpenMP)
if (OpenMP_CXX_FLAGS)
   set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
   set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
endif ()



#
ADD_EXECUTABLE( simplestereo main_stereo.c blockmatch_stereo.c blockmatch.c shear.c fimage2.c)
TARGET_LINK_LIBRARIES( simplestereo ${LIBS} )
SET_TARGET_PROPERTIES( simplestereo PROPERTIES 
   COMPILE_FLAGS "-DIIFLOAT='long double'")

## 
#ADD_EXECUTABLE (simplestereo_float main_stereo.c blockmatch_stereo.c blockmatch.c shear.c fimage2.c)
#TARGET_LINK_LIBRARIES(simplestereo_float ${LIBS})
#SET_TARGET_PROPERTIES( simplestereo_float PROPERTIES 
#   COMPILE_FLAGS "-DIIFLOAT='float'")
#
##
#ADD_EXECUTABLE (simplestereo_double main_stereo.c blockmatch_stereo.c blockmatch.c shear.c fimage2.c)
#TARGET_LINK_LIBRARIES(simplestereo_double ${LIBS})
#SET_TARGET_PROPERTIES( simplestereo_double PROPERTIES 
#   COMPILE_FLAGS "-DIIFLOAT='double'")

#
ADD_EXECUTABLE( simplebm main_blockmatch.c blockmatch_2d.c blockmatch_stereo.c blockmatch.c shear.c fimage2.c)
TARGET_LINK_LIBRARIES( simplebm ${LIBS})


