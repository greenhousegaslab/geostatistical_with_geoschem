ml intel/2022.0  
export ROOT_LIBRARY_DIR="/data/apps/linux-centos8-cascadelake/intel-19.1.2.254/netcdf-c-4.7.4-6y6obikbrvgtborvgu4fiqeakyqdbvlr"
export ROOT_LIBRARY_DIRF="/data/apps/linux-centos8-cascadelake/intel-19.1.2.254/netcdf-fortran-4.5.3-ao6brusl3l7xopfg27667ihk4mcseg5n"
export ROOT_LIBRARY_HDF5="/data/apps/linux-centos8-cascadelake/intel-19.1.2.254/hdf5-1.10.7-mscsqqv4iu6wc7ax75wkassr37dxzrlf"
export GC_INCLUDE="$ROOT_LIBRARY_DIR/include"
export GC_F_INCLUDE="$ROOT_LIBRARY_DIRF/include"
export GC_F_BIN="$ROOT_LIBRARY_DIRF/bin"
export GC_BIN="$ROOT_LIBRARY_DIR/bin"
export LD_LIBRARY_PATH=$ROOT_LIBRARY_DIR/lib:$ROOT_LIBRARY_DIRF/lib:$ROOT_LIBRARY_HDF5/lib:$LD_LIBRARY_PATH
make clean
make SAT_NETCDF=yes USE_MKL=yes HDF=yes

