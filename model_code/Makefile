#==============================================================================
#
# GEOS-Chem Makefile for LINUX/ifort compiler
#
#==============================================================================

SHELL = /bin/sh

#==============================================================================
# Default settings
#==============================================================================

# OpenMP is turned on by default
ifndef OMP
OMP = yes
endif

# Turn on -traceback option by default
ifndef TRACEBACK
TRACEBACK=yes
endif

#==============================================================================
# Declare Options
#==============================================================================

# Pick compiler options for debug run or regular run
ifeq ($(DEBUG),yes)
FFLAGS   = -cpp -w -O0 -auto -noalign -convert big_endian -fp-model precise -CB -g -check all
else
FFLAGS   = -cpp -w -O3 -auto -noalign -convert big_endian -fp-model precise -vec-report0
endif

# Also add traceback option
ifeq ($(TRACEBACK),yes)
FFLAGS  += -traceback
endif

# Turn on OpenMP parallelization
ifeq ($(OMP),yes)
#FFLAGS  += -openmp -Dmultitask
FFLAGS  += -gopenmp -Dmultitask
endif

# Add special ifort optimization commands
ifeq ($(IPO),yes)
FFLAGS  += -ipo
endif

ifeq ($(STATIC),yes)
FFLAGS  += -static -parallel
endif

ifeq ($(AUTO_PARALLEL),yes)
FFLAGS  += -parallel
endif

# For arrays bigger than 2GB
#ifeq ($(BIG_ARRAYS),yes)
FFLAGS  += -mcmodel=medium -i-dynamic
#endif

F90 = ifort $(FFLAGS) $(INCLUDE)

# Library include path
INCLUDE   := -I$(GC_INCLUDE) -I$(GC_F_INCLUDE) -I$(ROOT_LIBRARY_HDF5)/include
# Library link path: first try to get the list of proper linking flags
# for this build of netCDF with nf-config and nc-config.
NCL       := $(shell $(GC_F_BIN)/nf-config --flibs)
NCL       += $(shell $(GC_BIN)/nc-config --libs)
NCL       := $(filter -l%,$(NCL)) 

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%% NOTE TO GEOS-CHEM USERS: If you do not have netCDF-4.2 installed
#%%%% Then you can add/modify the linking sequence here.  (This sequence
#%%%% is a guess, but is probably good enough for other netCDF builds.)
ifeq ($(NCL),)
NCL            :=-lnetcdf -lnetcdff -lhdf5_hl -lhdf5 -lz
endif
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Prepend the library directory path to the linking sequence
NCL       := -L$(GC_LIB) $(NCL)
LINK      := $(NCL)
ifeq ($(USE_MKL),yes)
LINK	+= -L$(MKLROOT)/lib/intel64_lin  $(MKLROOT)/lib/intel64_lin/libmkl_blas95_lp64.a  $(MKLROOT)/lib/intel64_lin/libmkl_lapack95_lp64.a -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -qopenmp -lpthread
LAPACK_BLAS_FFLAGS += -I$(MKLROOT)/include/intel64_lin/lp64/  -I$(MKLROOT)/include
endif

# Link to the HDF and HDF-EOS libraries
ifeq ($(HDF),yes)
HDFHOME=$(ROOT_LIBRARY_HDF5)
HDFINC=$(ROOT_LIBRARY_HDF5)/include
HDFLIB=$(ROOT_LIBRARY_HDF5)/lib

HDFEOS_HOME=$(ROOT_LIBRARY_HDF5)
HDFEOS_INC=$(ROOT_LIBRARY_HDF5)/include
HDFEOS_LIB=$(ROOT_LIBRARY_HDF5)/lib

HDF5HOME=$(ROOT_LIBRARY_HDF5)
HDF5INC=$(ROOT_LIBRARY_HDF5)/include
HDF5LIB=$(ROOT_LIBRARY_HDF5)/lib

FFLAGS +=  -I$(HDFEOS_INC) -I$(HDF5INC) -I$(HDFINC)
LINK +=  -L$(HDFEOS_LIB) -L$(HDF5LIB) -L$(HDFLIB) -lhdfeos -lGctp -lmfhdf -ldf -lz -lm -ljpeg -lsz -lhdf5 -lhdf5_hl -lhdf5hl_fortran -lhdf5_fortran -lhe5_hdfeos

endif

ifeq ($(SAT_NETCDF),yes)
LINK +=  -L /libmkl_solver_lp64.a -Wl,--start-group \
		-lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -Wl,--end-group \
		-qopenmp -lpthread
endif


#==============================================================================
# Include Objects
#==============================================================================

VPATH = ./modified ./adjoint ./new ./obs_operators ./NcdfUtil

ifeq ($(LIDORT),yes)
VPATH += ./lidort ./lidort/thread_sourcecode_MkII_F90
endif

include ./Objects.mk

# Add LIDORT Specific Code
ifeq ($(LIDORT),yes)
#====================
# LIDORT CODE
#====================
# dkh
#LIDORT_COMPILE  = ifort -c -warn all -check bounds -O3 -zero
LIDORT_COMPILE_FIXED = ifort -cpp -check bounds -O3 -zero  -noalign -fixed -qopenmp -Dmultitask
LIDORT_COMPILE       = ifort -c -check bounds -O3 -zero  -noalign -free -qopenmp -Dmultitask -traceback -CB -vec-report0
LAPACK_COMPILE  = ifort -c -warn all -check bounds -O3  -zero
LAPACK_NOPT_COMPILE = ifort -c -O3 -zero
FLINK = ifort

# Link definition
#################

#LINK.f90 = $(FLINK) -g -pg
LINK.f90 = $(FLINK)

# dkh
#LIDORT_PATH = ..
LIDORT_PATH = ./lidort

# other paths are relative

SPATH_S = $(LIDORT_PATH)/thread_sourcecode_MkII_F90/
OBJ     = $(LIDORT_PATH)/OBJECTS_F90

#  OBJECT MODULES

# LIDORT modules in directory sourcecode
MIE =                      \
RTS_mie_modules.o          \
RTS_mie_sourcecode.o       \
RTS_mie_sourcecode_plus.o  \
GC_forward_Mie.o           \
GC_adjoint_Mie.o

#  Masters set

OBJECTS_LIDORT_MASTERS      = $(OBJ)/lidort_masters_basic.o
OBJECTS_LIDORT_MASTERS_LCS  = $(OBJ)/lidort_masters_lcs.o
OBJECTS_LIDORT_MASTERS_LPS  = $(OBJ)/lidort_masters_lps.o

#  Basic set for Radiances

OBJECTS_LIDORT_BASIC     = $(OBJ)/lidort_solutions.o   \
                           $(OBJ)/lidort_bvproblem.o   \
                           $(OBJ)/lidort_intensity.o   \
                           $(OBJ)/lidort_corrections.o \
                           $(OBJ)/lidort_miscsetups.o  \
                           $(OBJ)/lidort_inputs.o      \
                           $(OBJ)/lidort_geometry.o

OBJECTS_LIDORT_AUX       = $(OBJ)/lidort_aux.o

OBJECTS_LIDORT_LA        = $(OBJ)/lidort_la_solutions.o   \
                           $(OBJ)/lidort_la_miscsetups.o

OBJECTS_LIDORT_LC        = $(OBJ)/lidort_lc_bvproblem.o   \
                           $(OBJ)/lidort_lc_wfatmos.o     \
                           $(OBJ)/lidort_lc_corrections.o \
                           $(OBJ)/lidort_lc_miscsetups.o

OBJECTS_LIDORT_LP        = $(OBJ)/lidort_lp_bvproblem.o   \
                           $(OBJ)/lidort_lp_wfatmos.o     \
                           $(OBJ)/lidort_lp_corrections.o \
                           $(OBJ)/lidort_lp_miscsetups.o

OBJECTS_LIDORT_LS        = $(OBJ)/lidort_ls_wfsurface.o   \
                           $(OBJ)/lidort_ls_corrections.o


# LIDORT environment & interface modules
#OBJECTS_LIDORT_3P5T_LPS_MT = $(OBJ)/lidort_mod.o

endif

#=============================================================================
#  Executables and Documentation
#=============================================================================

ifeq ($(LIDORT),yes)

geos:   $(MODS) $(OBJS) $(OBJSe) $(FJ)            \
                  $(OBJECTS_LIDORT_MASTERS_LPS)   \
                  $(OBJECTS_LIDORT_AUX)           \
                  $(OBJECTS_LIDORT_BASIC)         \
                  $(OBJECTS_LIDORT_LA)            \
                  $(OBJECTS_LIDORT_LP)            \
                  $(OBJECTS_LIDORT_LS)            \
                  $(MIE)
#	$(F90)  $(MODS) $(OBJS)  $(OBJSe) $(FJ) $(LIBS) -o geos
#	$(F90)  *.o -o geos
#  	$(F90)  $(MODS) $(OBJS)  $(OBJSe) $(FJ) -o geos
	$(F90)  $(MODS) $(OBJS)  $(OBJSe) $(FJ)   \
                  $(OBJECTS_LIDORT_MASTERS_LPS)   \
                  $(OBJECTS_LIDORT_AUX)           \
                  $(OBJECTS_LIDORT_BASIC)         \
                  $(OBJECTS_LIDORT_LA)            \
                  $(OBJECTS_LIDORT_LP)            \
                  $(OBJECTS_LIDORT_LS)            \
                  $(MIE) $(LINK) -o geos
else
geos:   $(MODS) $(OBJS) $(OBJSe) $(FJ)
	$(F90)  $(MODS) $(OBJS)  $(OBJSe) $(FJ) \
                $(LINK) -o geos

endif

# Build GEOS-Chem documenation w/ ProTeX
doc:
	@$(MAKE) -C doc all

# Remove all *.tex, *.ps, and *.pdf files from the doc subdirectory
docclean:
	@$(MAKE) -C doc clean

help:
	@echo '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
	@echo '%%%      GEOS-Chem Adjoint Help Screen     %%%'
	@echo '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
	@echo ''
	@echo 'Usage: make TARGET [ OPTIONAL-FLAGS ]'
	@echo ''
	@echo 'TARGET may be one of the following:'
	@echo 'geos			Builds GEOS-Chem Adjoint source code'
	@echo 'clean			Removes *.o, *.mod files and geos executable'
	@echo 'OPTIONAL-FLAGS may be:'
	@echo 'DEBUG=yes      Builds GEOS-Chem for a debugger (with -g -O0)'
	@echo 'HDF=yes       Enables writing diagnostic timeseries output to HDF files'
	@echo 'IPO=yes        Turns on optmization options -ipo -static (default is no)'
	@echo 'LIDORT=yes     Enables LIDORT Modules'
	@echo 'OMP=[yes|no]   Turns OpenMP parallelization on/off (default is yes)'
	@echo 'SAT_NETCDF=yes     Enables Satellite NetCDF'
	@echo 'TRACEBACK=yes  Turns on -traceback option (default is yes)'
	@echo ''
	@echo 'NOTE: This installation is set up to work with Intel Fortran Compilers only'

#==============================================================================
# Include Dependencies
#==============================================================================

include ./Dependencies.mk

ifeq ($(HDF),yes)
#====================
# MOPITT CODE
#====================
gvchsq.o             : gvchsq.f
	$(F90) -c -r8 obs_operators/gvchsq.f
HdfIncludeModule.o   : HdfIncludeModule.f90
	$(F90) -c -r8 obs_operators/HdfIncludeModule.f90
HdfSdModule.o        : HdfSdModule.f90
	$(F90) -c -r8 obs_operators/HdfSdModule.f90
HdfVdModule.o        : HdfVdModule.f90
	$(F90) -c -r8 obs_operators/HdfVdModule.f90
interp.o             : interp.f
	$(F90) -c -r8 obs_operators/interp.f
gaussj.o             : gaussj.f
	$(F90) -c -r8 obs_operators/gaussj.f
mopitt_obs_mod.o     : mopitt_obs_mod.f  CMN CMN_SIZE define.h define_adj.h
	$(F90) -c -r8 obs_operators/mopitt_obs_mod.f

#====================
# OMI NO2 CODE
#====================
omi_no2_obs_mod.o     : omi_no2_obs_mod.f90 CMN CMN_SIZE define.h define_adj.h
	$(F90) -c -r8 obs_operators/omi_no2_obs_mod.f90

#====================
# OMI L3 SO2
#====================
omi_so2_obs_mod.o   : omi_so2_obs_mod.f
	$(F90) -c -r8 obs_operators/omi_so2_obs_mod.f

#====================
# AIRS CODE
#====================
He4IncludeModule.o        : He4IncludeModule.f90
	$(F90) -c -r8 obs_operators/He4IncludeModule.f90
He4ErrorModule.o          : He4ErrorModule.f90
	$(F90) -c -r8 obs_operators/He4ErrorModule.f90
He4GridModule.o           : He4GridModule.f90
	$(F90) -c -r8 obs_operators/He4GridModule.f90
He4SwathModule.o          : He4SwathModule.f90
	$(F90) -c -r8 obs_operators/He4SwathModule.f90
airsv5_mod.o              : airsv5_mod.f90
	$(F90) -c -r8 obs_operators/airsv5_mod.f90
airs_co_obs_mod.o    : airs_co_obs_mod.f CMN_SIZE define.h
	$(F90) -c -r8 obs_operators/airs_co_obs_mod.f
findinv.o            : findinv.f
	$(F90) -c -r8 obs_operators/findinv.f
endif

ifeq ($(SAT_NETCDF),yes)
#====================
# SATELLITE CO2  CODE
#====================
interpolate.o      : interpolate.f90
	$(F90) -c -r8 obs_operators/interpolate.f90
gosat_co2_mod.o      : gosat_co2_mod.f
	$(F90) -c -r8 obs_operators/gosat_co2_mod.f
flask_co2_mod.o      : flask_co2_mod.f
	$(F90) -c -r8 obs_operators/flask_co2_mod.f
tes_nh3_mod.o        : tes_nh3_mod.f
	$(F90) -c -r8 obs_operators/tes_nh3_mod.f
tes_o3_mod.o         : tes_o3_mod.f
	$(F90) -c -r8 obs_operators/tes_o3_mod.f
tes_o3_irk_mod.o     : tes_o3_irk_mod.f
	$(F90) -c -r8 obs_operators/tes_o3_irk_mod.f

#====================
# MODIS AOD CODE (xxu, dkh, 01/09/12, adj32_011)
#====================
modis_aod_obs_mod.o  : modis_aod_obs_mod.f
	$(F90) -c -r8 obs_operators/modis_aod_obs_mod.f

#====================
# SCIA CODE
#====================
scia_ch4_mod.o       : scia_ch4_mod.f CMN_SIZE
	$(F90) -c -r8 $<
endif

ifeq ($(LIDORT),yes)
#====================
# LIDORT CODE
#====================
#--------------------------------------------------
#--------------------------Environment modules-----
#--------------------------------------------------

#lidort_mod.o: lidort_mod.f \
#	 $(SPATH_S)LIDORT.PARS_F90
#	$(F90) $(LIDORT_PATH)/lidort_mod.f
#	$(LIDORT_COMPILE_FIXED) $(LIDORT_PATH)/lidort_mod.f90
#$(LIDORT_COMPILE_FIXED) $(LIDORT_PATH)/lidort_mod.f90 -o lidort_mod.o
lidort_mod.o                   : lidort_mod.f LIDORT.PARS_F90
	$(F90) -c -r8 lidort/lidort_mod.f

mie_mod.o                      : mie_mod.f
	$(F90) -c -r8 lidort/mie_mod.f

RTS_mie_modules.o              : RTS_mie_modules.f90
	$(F90) -c -r8 lidort/RTS_mie_modules.f90

RTS_mie_sourcecode.o           : RTS_mie_sourcecode.f90
	$(F90) -c -r8 lidort/RTS_mie_sourcecode.f90

RTS_mie_sourcecode_plus.o      : RTS_mie_sourcecode_plus.f90
	$(F90) -c -r8 lidort/RTS_mie_sourcecode_plus.f90

GC_forward_Mie.o               : GC_forward_Mie.f90
	$(F90) -c -r8 lidort/GC_forward_Mie.f90

GC_adjoint_Mie.o               : GC_adjoint_Mie.f90
	$(F90) -c -r8 lidort/GC_adjoint_Mie.f90

#----------------------------------------------------
#----------------------LIDORT master modules --------
#----------------------------------------------------

$(OBJ)/lidort_masters_lps.o: $(SPATH_S)lidort_masters_lps.f90    \
	 $(SPATH_S)LIDORT.PARS_F90
	$(LIDORT_COMPILE) $(SPATH_S)lidort_masters_lps.f90 -o $(OBJ)/lidort_masters_lps.o

#----------------------------------------------------
#----------------------LIDORT Radiance modules ------
#----------------------------------------------------

$(OBJ)/lidort_solutions.o: $(SPATH_S)lidort_solutions.f90    \
	 $(SPATH_S)LIDORT.PARS_F90
	$(LIDORT_COMPILE) $(SPATH_S)lidort_solutions.f90 -o $(OBJ)/lidort_solutions.o

$(OBJ)/lidort_bvproblem.o: $(SPATH_S)lidort_bvproblem.f90    \
	 $(SPATH_S)LIDORT.PARS_F90
	$(LIDORT_COMPILE) $(SPATH_S)lidort_bvproblem.f90 -o $(OBJ)/lidort_bvproblem.o

$(OBJ)/lidort_geometry.o: $(SPATH_S)lidort_geometry.f90    \
	 $(SPATH_S)LIDORT.PARS_F90
	$(LIDORT_COMPILE) $(SPATH_S)lidort_geometry.f90 -o $(OBJ)/lidort_geometry.o

$(OBJ)/lidort_intensity.o: $(SPATH_S)lidort_intensity.f90    \
	 $(SPATH_S)LIDORT.PARS_F90
	$(LIDORT_COMPILE) $(SPATH_S)lidort_intensity.f90 -o $(OBJ)/lidort_intensity.o

$(OBJ)/lidort_miscsetups.o: $(SPATH_S)lidort_miscsetups.f90    \
	 $(SPATH_S)LIDORT.PARS_F90
	$(LIDORT_COMPILE) $(SPATH_S)lidort_miscsetups.f90 -o $(OBJ)/lidort_miscsetups.o

$(OBJ)/lidort_corrections.o: $(SPATH_S)lidort_corrections.f90    \
	 $(SPATH_S)LIDORT.PARS_F90
	$(LIDORT_COMPILE) $(SPATH_S)lidort_corrections.f90 -o $(OBJ)/lidort_corrections.o

$(OBJ)/lidort_inputs.o: $(SPATH_S)lidort_inputs.f90    \
	 $(SPATH_S)LIDORT.PARS_F90
	$(LIDORT_COMPILE) $(SPATH_S)lidort_inputs.f90 -o $(OBJ)/lidort_inputs.o

# ---------------------------------------------------
#----------------------LIDORT Auxiliary module ------
# ---------------------------------------------------

$(OBJ)/lidort_aux.o: $(SPATH_S)lidort_aux.f90    \
	 $(SPATH_S)LIDORT.PARS_F90
	$(LAPACK_COMPILE) $(SPATH_S)lidort_aux.f90 -o $(OBJ)/lidort_aux.o

# ---------------------------------------------------
#----------------------LIDORT Jacobian modules ------
# ---------------------------------------------------

#  General

$(OBJ)/lidort_la_solutions.o: $(SPATH_S)lidort_la_solutions.f90    \
	 $(SPATH_S)LIDORT.PARS_F90
	$(LIDORT_COMPILE) $(SPATH_S)lidort_la_solutions.f90 -o $(OBJ)/lidort_la_solutions.o

$(OBJ)/lidort_la_miscsetups.o: $(SPATH_S)lidort_la_miscsetups.f90    \
	 $(SPATH_S)LIDORT.PARS_F90
	$(LIDORT_COMPILE) $(SPATH_S)lidort_la_miscsetups.f90 -o $(OBJ)/lidort_la_miscsetups.o

#  Column specific

$(OBJ)/lidort_lc_bvproblem.o: $(SPATH_S)lidort_lc_bvproblem.f90    \
	 $(SPATH_S)LIDORT.PARS_F90
	$(LIDORT_COMPILE) $(SPATH_S)lidort_lc_bvproblem.f90 -o $(OBJ)/lidort_lc_bvproblem.o

$(OBJ)/lidort_lc_wfatmos.o: $(SPATH_S)lidort_lc_wfatmos.f90    \
	 $(SPATH_S)LIDORT.PARS_F90
	$(LIDORT_COMPILE) $(SPATH_S)lidort_lc_wfatmos.f90 -o $(OBJ)/lidort_lc_wfatmos.o

$(OBJ)/lidort_lc_corrections.o: $(SPATH_S)lidort_lc_corrections.f90    \
	 $(SPATH_S)LIDORT.PARS_F90
	$(LIDORT_COMPILE) $(SPATH_S)lidort_lc_corrections.f90 -o $(OBJ)/lidort_lc_corrections.o

$(OBJ)/lidort_lc_miscsetups.o: $(SPATH_S)lidort_lc_miscsetups.f90    \
	 $(SPATH_S)LIDORT.PARS_F90
	$(LIDORT_COMPILE) $(SPATH_S)lidort_lc_miscsetups.f90 -o $(OBJ)/lidort_lc_miscsetups.o

#  Profile specific

$(OBJ)/lidort_lp_bvproblem.o: $(SPATH_S)lidort_lp_bvproblem.f90    \
	 $(SPATH_S)LIDORT.PARS_F90
	$(LIDORT_COMPILE) $(SPATH_S)lidort_lp_bvproblem.f90 -o $(OBJ)/lidort_lp_bvproblem.o

$(OBJ)/lidort_lp_wfatmos.o: $(SPATH_S)lidort_lp_wfatmos.f90    \
	 $(SPATH_S)LIDORT.PARS_F90
	$(LIDORT_COMPILE) $(SPATH_S)lidort_lp_wfatmos.f90 -o $(OBJ)/lidort_lp_wfatmos.o

$(OBJ)/lidort_lp_corrections.o: $(SPATH_S)lidort_lp_corrections.f90    \
	 $(SPATH_S)LIDORT.PARS_F90
	$(LIDORT_COMPILE) $(SPATH_S)lidort_lp_corrections.f90 -o $(OBJ)/lidort_lp_corrections.o

$(OBJ)/lidort_lp_miscsetups.o: $(SPATH_S)lidort_lp_miscsetups.f90    \
	 $(SPATH_S)LIDORT.PARS_F90
	$(LIDORT_COMPILE) $(SPATH_S)lidort_lp_miscsetups.f90 -o $(OBJ)/lidort_lp_miscsetups.o

#  Surface

$(OBJ)/lidort_ls_wfsurface.o: $(SPATH_S)lidort_ls_wfsurface.f90  \
	 $(SPATH_S)LIDORT.PARS_F90
	$(LIDORT_COMPILE) $(SPATH_S)lidort_ls_wfsurface.f90 -o $(OBJ)/lidort_ls_wfsurface.o

$(OBJ)/lidort_ls_corrections.o: $(SPATH_S)lidort_ls_corrections.f90    \
	 $(SPATH_S)LIDORT.PARS_F90
	$(LIDORT_COMPILE) $(SPATH_S)lidort_ls_corrections.f90 -o $(OBJ)/lidort_ls_corrections.o

#----------------------------------------------------
#----------- BRDF supplement modules ----------------
#----------------------------------------------------

$(OBJ)/lidort_brdf_supplement.o: $(SPATH_S)lidort_brdf_supplement.f90    \
	 $(SPATH_S)LIDORT.PARS_F90
	$(LIDORT_COMPILE) $(SPATH_S)lidort_brdf_supplement.f90 -o $(OBJ)/lidort_brdf_supplement.o

$(OBJ)/lidort_brdf_kernels.o: $(SPATH_S)lidort_brdf_kernels.f90    \
	 $(SPATH_S)LIDORT.PARS_F90
	$(LIDORT_COMPILE) $(SPATH_S)lidort_brdf_kernels.f90 -o $(OBJ)/lidort_brdf_kernels.o

$(OBJ)/lidort_brdf_ls_supplement.o: $(SPATH_S)lidort_brdf_ls_supplement.f90    \
	 $(SPATH_S)LIDORT.PARS_F90
	$(LIDORT_COMPILE) $(SPATH_S)lidort_brdf_ls_supplement.f90 -o $(OBJ)/lidort_brdf_ls_supplement.o

$(OBJ)/lidort_brdf_ls_kernels.o: $(SPATH_S)lidort_brdf_ls_kernels.f90    \
	 $(SPATH_S)LIDORT.PARS_F90
	$(LIDORT_COMPILE) $(SPATH_S)lidort_brdf_ls_kernels.f90 -o $(OBJ)/lidort_brdf_ls_kernels.o
endif

#==============================================================================
#  Other compilation commands
#==============================================================================
ifort_errmsg.o                 : ifort_errmsg.f
linux_err.o                    : linux_err.c
	$(CC) -c linux_err.c

#=============================================================================
#  Other Makefile Commands
#=============================================================================
clean:
	rm -rf *.o *.mod ifc* geos rii_files

.PHONY: clean doc docclean

.SUFFIXES: .f .F .f90 .F90
.f.o:			; $(F90) -c $*.f
.F.o:			; $(F90) -c $*.F
.f90.o:                 ; $(F90) -c -free $*.f90
.F90.o:                 ; $(F90) -c -free $*.F90

%.o : %.mod

