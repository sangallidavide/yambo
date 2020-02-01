BLAS_LIB="-L/usr/local/opt/openblas/lib -lopenblas"
LAPACK_LIB="-L/usr/local/opt/lapack/lib -llapack"
NETCDFF_LIB="-L/usr/local/lib -lnetcdff"
NETCDF_LIB="-L/usr/local/lib -lnetcdf"
HDF5_LIB="-L/usr/local/lib -lhdf5_fortran -lhdf5_hl -lhdf5"
#
NETCDF_INC="/usr/local/include"
NETCDFF_INC="/usr/local/include"
HDF5_INC="/usr/local/include"
LIBXC_LIBS="-L/usr/local/lib -lxc"
LIBXC_INC="/usr/local/include"
BLAS_INC="/usr/local/opt/openblas/include"
LAPACK_INC="/usr/local/opt/lapack/include"


./configure \
--with-extlibs-path="/Users/maurizia/Desktop/varie/D_PGM/D_Exc/yambo-libs" \
--with-libxc-libs="$LIBXC_LIBS" --with-libxc-includedir="$LIBXC_INC" \
--with-netcdf-libs="$NETCDF_LIB" --with-netcdf-includedir="$NETCDF_INC" \
--with-netcdff-libs="$NETCDFF_LIB" --with-netcdff-includedir="$NETCDFF_INC" \
--with-hdf5-libs="$HDF5_LIB" --with-hdf5-includedir="$HDF5_INC" \
--with-blas-libs="$BLAS_LIB" --with-lapack-libs="$LAPACK_LIB"  \
--enable-keep-src  \
--without-editor --enable-msgs-comps \
CC="gcc" \
FC="gfortran" \
F77="gfortran"
