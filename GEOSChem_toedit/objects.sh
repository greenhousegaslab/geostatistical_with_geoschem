#!/bin/bash

INPUT_FOLDER=`cat INPUT_FOLDER`

cp Objects.default Objects.mkl

# Note this line gives an error but WORKS correctly
if  [ -n $(grep ' => offdiagonal                : F' $INPUT_FOLDER/input.gcadj) ]; then

find="geos_chem_adj_mod.o"
replace="covariance_mod.o 	      \
geos_chem_adj_mod.o"

sed -e "s/$find/$replace/g" Objects.mkl > Objects.off-diag
mv Objects.off-diag Objects.mkl
fi

# Note this line gives an error but WORKS correctly
if  [ -n $(grep 'Compute BFGS inverse Hessian   : F' $INPUT_FOLDER/input.gcadj) ]; then

find="inv_hessian_mod.o"

replace="inv_hessian_mod.o             \
inv_hessian_lbfgs_mod.o"

sed -e "s/$find/$replace/g" Objects.mkl > Objects.lbfgs
mv Objects.lbfgs Objects.mkl
fi


if [ $1 = "HDF" ]; then

find="getifsun.o"
replace="getifsun.o                    \
gvchsq.o "

sed -e "s/$find/$replace/g" Objects.mkl > output1

find="input_mod.o"
replace="input_mod.o            \
He4IncludeModule.o            \
He4ErrorModule.o              \
He4GridModule.o               \
He4SwathModule.o              \
findinv.o                     \
airsv5_mod.o                  \
airs_co_obs_mod.o             \
HdfIncludeModule.o            \
HdfSdModule.o                 \
HdfVdModule.o                 \
interp.o                      \
gaussj.o                      \
mopitt_obs_mod.o              \
omi_no2_obs_mod.o             \
omi_so2_obs_mod.o"

sed -e "s/$find/$replace/g" output1 > output
rm output1
mv output Objects.mk
fi

if [ $1 = "SAT_NETCDF" ]; then

find="rpmares_mod.o"
replace="rpmares_mod.o                 \
interpolate.o                 \
gosat_co2_mod.o               \
flask_co2_mod.o               \
tes_nh3_mod.o                 \
tes_o3_mod.o                  \
tes_o3_irk_mod.o"

sed -e "s/$find/$replace/g" Objects.mkl > output1

find="tes_ch4_mod.o"
replace="tes_ch4_mod.o                 \
scia_ch4_mod.o                         \
modis_aod_obs_mod.o"

sed -e "s/$find/$replace/g" output1 > output
rm output1
mv output Objects.mk
fi

if [ $1 = "LIDORT" ]; then

find="population_mod.o"
replace="population_mod.o              \
mie_mod.o                     \
lidort_mod.o"

sed -e "s/$find/$replace/g" Objects.mkl > output
mv output Objects.mk
fi

if [ $1 = "DEFAULT" ]; then
mv Objects.mkl Objects.mk
fi
