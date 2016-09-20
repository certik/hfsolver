#!/usr/bin/env bash

# Exit on error
set -e
# Echo each command
set -x

# Print CPU info for debugging purposes:
cat /proc/cpuinfo
gfortran -v
/usr/bin/f95 -v

# Install cmake
wget --no-check-certificate -O- http://www.cmake.org/files/v3.3/cmake-3.3.0-Linux-x86_64.tar.gz | tar xz
export PATH="`pwd`/cmake-3.3.0-Linux-x86_64/bin/:$PATH"

cmake --version
if [[ "${TEST_MPI}" == "yes" ]]; then
    export CC=mpicc CXX=mpicxx FC=mpif90
fi
cmake \
    -DCMAKE_INSTALL_PREFIX="$VIRTUAL_ENV" \
    -DWITH_LIBINT=${TEST_LIBINT} \
    -DLIBINT_DIR=$HOME/usr \
    -DWITH_FFTW=${TEST_FFTW} \
    -DCMAKE_BUILD_TYPE=${BUILD_TYPE} \
    -DWITH_PYTHON=${TEST_PYTHON} \
    -DWITH_OPENMP=${TEST_OPENMP} \
    -DWITH_MPI=${TEST_MPI} \
    .
make
make install
