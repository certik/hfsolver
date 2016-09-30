#!/usr/bin/env bash

# Exit on error
set -e
# Echo each command
set -x

if [[ "${TEST_PYTHON}" == "yes" ]]; then
    cd src/tests/bessel
    python test_modified_bessel.py
else
    if [[ "${TEST_SLOW}" == "yes" ]]; then
      ctest -V --output-on-failure -L slow
    elif [[ "${TEST_MPI}" == "yes" ]]; then
      export OMPI_MCA_mpi_yield_when_idle=1
      ctest --output-on-failure -L mpi
    else
      ctest --output-on-failure -LE slow
    fi
fi
