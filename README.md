# geostatistical_with_geoschem

Code written by Ruixue Lei and Scot M. Miller

Run a geostatistical inverse model for CO2 using the GEOS-Chem adjoint model.

You need to compile fortran code and prepare .m files to run GEOS-Chem adjoint model.

For fortran code, enter model_code folder. 
cd geos_adjoint/model_code
Run chmod command to give run permissioin for two scripts first, then run two scripts.
compile_adjoint.sh will report an error of missing some libs, mx will fix it and genertate "geos". Move "geos" to your geoschem run directory.

chmod +x compile_adjoint.sh mx
./compile_adjoint.sh
./mx

As for .m files, they can be run by Matlab or Octave.

If you want to run it by Matlab, buy your own license or ask server admin for license.
Modify cost_gradient_fun.m in run dir.
Search all functions starting with "netcdf_" to "netcdf.". Then your code is already for matlab.


Octave is free software. If you want to run it by Octave, you may need to instal the lastest version by yourself. Below is the guidence:
1. Download and install miniconda.
Download page: https://docs.conda.io/en/latest/miniconda.html
Install guide: https://conda.io/projects/conda/en/latest/user-guide/install/linux.html

2. build a new python environment and install Octave.
conda create --yes --name octave
conda activate octave
conda install -c conda-forge octave_kernel netcdf

3. Set environment variables:
export CC=icc
export CXX=icpc
export FC=ifort
export F90=ifort
export F77=ifort
export NC_CONFIG=~/miniconda3/envs/octave/bin/nc_config

4. enter Octave gui mode
octave --gui

5. In Octave command window
pkg list
See if pkg command works. If no, exit octave and re-enter Octave gui mode. If yes, you can use "octave --no-gui" to enter no gui mode to speed up when you work on a server.

6. Install octave netcdf package, in Octave command window.
pkg install -forge netcdf
or you can download Octave netcdf package from https://octave.sourceforge.io/netcdf/index.html and install it mannually.

Then your code is ready for Octave.


