#!/usr/bin/python2
if __name__ == '__main__':
  import sys
  import os
  sys.path.insert(0, os.path.abspath('config'))
  import configure
  configure_options = [
    '--PETSC_ARCH=yambo_double_complex',
    '--prefix=/home/maurizia/D_PGM/D_exc/yambo-libs/gfortran/mpifort/double',
    '--with-blas-lapack-lib=-llapack -lblas',
    '--with-cc=mpicc',
    '--with-cxx=0',
    '--with-fc=mpifort',
    '--with-precision=double',
    '--with-scalar-type=complex',
    '--with-shared-libraries=0',
    '--with-ssl=0',
    '--with-x=0',
  ]
  configure.petsc_configure(configure_options)
