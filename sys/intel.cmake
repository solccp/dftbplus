#
# Toolchain file for
#
# Intel compiler, MKL library
#
# Notes:
#
#  * CMake format: Command line options (e.g. compiler flags) space separated, other kind
#    of lists semicolon separated.
#
#  * Variables containing library search paths are empty by default. The CMAKE_PREFIX_PATH
#    environment variable should be set up correctly, so that CMake can find those libraries
#    automatically. If that is not the case, override those variables to add search paths
#    manually
#
#  * Compiler flags specified via the environment variables FFLAGS and CFLAGS are *appended* to the
#    pre-configured flags. If you want to override all flags on the command line, use the
#    -DFortran_FLAGS=<flags> and -DC_FLAGS=<flags> options.
#


#
# Fortran compiler settings
#
set(Fortran_FLAGS "-standard-semantics"
  CACHE STRING "Additional general Fortran compiler flags")

set(Fortran_FLAGS_RELEASE "-ip -heap-arrays 10"
  CACHE STRING "Additional Fortran compiler flags for Release build")

set(Fortran_FLAGS_DEBUG "-warn all -stand f08 -check -diag-error-limit 1 -traceback"
  CACHE STRING "Additional Fortran compiler flags for Debug build")

set(FYPP_FLAGS "" CACHE STRING "Flags for the preprocessor")


#
# C compiler settings
#
set(C_FLAGS ""
  CACHE STRING "Additional general C compiler flags")

set(C_FLAGS_RELEASE "-ip"
  CACHE STRING  "Additional C compiler flags for Release build")

set(C_FLAGS_DEBUG "-Wall"
  CACHE STRING "Additional C compiler flags for Debug build")


#
# External libraries
#

# LAPACK and BLAS
if(WITH_OMP)
  set(LAPACK_LIBRARIES "mkl_intel_lp64;mkl_intel_thread;mkl_core" CACHE STRING
    "LAPACK and BLAS libraries to link")
else()
  set(LAPACK_LIBRARIES "mkl_intel_lp64;mkl_sequential;mkl_core" CACHE STRING
    "LAPACK and BLAS libraries to link")
endif()

#set(LAPACK_LIBRARY_DIRS "$ENV{MKLROOT}/lib/intel64" CACHE STRING
#  "Directories where LAPACK and BLAS libraries can be found")

# ARPACK -- only needed when built with ARPACK support
set(ARPACK_LIBRARIES "arpack" CACHE STRING "Arpack library")
set(ARPACK_LIBRARY_DIRS "" CACHE STRING "Directories where Arpack library can be found")

# ScaLAPACK -- only needed for MPI-parallel build
set(SCALAPACK_LIBRARIES "mkl_scalapack_lp64;mkl_blacs_intelmpi_lp64" CACHE STRING
  "Scalapack libraries to link")
set(SCALAPACK_LIBRARY_DIRS "" CACHE STRING
  "Directories where Scalapack libraries can be found")

# ELSI -- only needed when compiled with ELSI support
set(ELSI_ROOT "" CACHE STRING "Root directory of the ELSI installation")

set(ELSI_EXTERNAL_LIBRARIES "" CACHE STRING
  "Any EXTERNAL libraries ELSI needs apart of its own libraries (and scalapack)")
set(ELSI_EXTERNAL_LIBRARY_DIRS "" CACHE STRING
  "Directories where ELSI external libraries can be found")

# PEXSI -- only needed when ELSI was compiled with PEXSI support
# Note: PEXSI usually needs explicit linking of the standard C++ library. Make sure to
#     provide the library path to that C++ standard library, which was used to compile PEXSI.
set(PEXSI_EXTERNAL_LIBRARIES "" CACHE STRING
  "Any EXTERNAL libraries PEXSI needs apart of its own libraries")
set(PEXSI_EXTERNAL_LIBRARY_DIRS "" CACHE STRING "Directories with PEXSI external libraries")

# PLUMED -- only needed when compiled with PLUMED support
set(PLUMED_LIBRARIES "plumed;plumedKernel" CACHE STRING "Libraries to link for PLUMED support")
set(PLUMED_LIBRARY_DIRS "" CACHE STRING "Directories to scan for PLUMED libraries")

# MAGMA -- only needed when compiled with GPU support
set(MAGMA_LIBRARIES "magma" CACHE STRING "Magma library")
set(MAGMA_LIBRARY_DIRS "" CACHE STRING "Directories to scan for MAGMA library")
set(MAGMA_INCLUDE_DIRS "" CACHE STRING "Directories to scan for MAGMA include files")

# Any other library needed to be linked or considered as include
set(OTHER_LIBRARIES "" CACHE STRING "Other libraries to link")
set(OTHER_LIBRARY_DIRS "" CACHE STRING "Directories where the other libraries can be found")
set(OTHER_INCLUDE_DIRS "" CACHE STRING "Other include directories to consider")
