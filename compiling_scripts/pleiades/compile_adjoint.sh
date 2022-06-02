#change var define first in geos_adjoint/model_code/adjoint/define_adj.h
module load comp-intel/2020.4.304 mpi-hpe/mpt.2.23 hdf4/4.2.12 hdf5/1.8.18_mpt netcdf/4.4.1.1_mpt
export ROOT_LIBRARY_DIR="/nasa/netcdf/4.4.1.1_mpt/lib"
export ROOT_LIBRARY_DIRF="/nasa/netcdf/4.4.1.1_mpt/lib"
export ROOT_LIBRARY_HDF5="/nasa/hdf5/1.12.0_mpt"
export GC_INCLUDE="$ROOT_LIBRARY_DIR/include"
export GC_F_INCLUDE="$ROOT_LIBRARY_DIRF/include"
export GC_F_BIN="$ROOT_LIBRARY_DIRF/bin"
export GC_BIN="$ROOT_LIBRARY_DIR/bin"
export LD_LIBRARY_PATH=$ROOT_LIBRARY_DIR/lib:$ROOT_LIBRARY_DIRF/lib:$ROOT_LIBRARY_HDF5/lib:$LD_LIBRARY_PATH
export CC=icc
export FC=ifort
export CXX=icpc
# for octave netcdf install
module load gcc/10.2
export NC_CONFIG=/home3/rlei/miniconda3/envs/bbb/bin/nc-config

export CC=gcc
export FC=gfortran
export CXX=g++
# for octave netcdf install

make clean
make SAT_NETCDF=yes USE_MKL=yes HDF=yes
