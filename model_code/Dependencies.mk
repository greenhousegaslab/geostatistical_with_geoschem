#==============================================================================
#  Dependencies Listing
#==============================================================================
BLKSLV.o                       : BLKSLV.f                       jv_mie.h
CLDSRF.o                       : CLDSRF.f                       cmn_fj.h CMN_SIZE define.h jv_cmn.h
CO_strat_pl.o                  : CO_strat_pl.f                  CMN_SIZE define.h
critical_load_mod.o            : critical_load_mod.f            CMN_SIZE
EFOLD.o                        : EFOLD.f
FLINT.o                        : FLINT.f
GAUSSP.o                       : GAUSSP.f
GEN.o                          : GEN.f                          jv_mie.h
fjx_acet_mod.o                 : fjx_acet_mod.f                 cmn_fj.h jv_cmn.h
JRATET.o                       : JRATET.f                       cmn_fj.h CMN_SIZE define.h jv_cmn.h
JVALUE.o                       : JVALUE.f                       cmn_fj.h CMN_SIZE define.h jv_cmn.h
Kr85_mod.o                     : Kr85_mod.f                     CMN_DIAG CMN_O3 CMN_SIZE define.h define.h
LEGND0.o                       : LEGND0.f
MATIN4.o                       : MATIN4.f
MIESCT.o                       : MIESCT.f                       jv_mie.h
NOABS.o                        : NOABS.f
OPMIE.o                        : OPMIE.f                        cmn_fj.h CMN_SIZE define.h jv_cmn.h jv_mie.h
RD_TJPL.o                      : RD_TJPL.f                      cmn_fj.h CMN_SIZE define.h jv_cmn.h
RnPbBe_mod.o                   : RnPbBe_mod.f                   CMN_DEP CMN_DIAG CMN_SIZE define.h define.h
SPHERE.o                       : SPHERE.f                       cmn_fj.h CMN_SIZE define.h jv_cmn.h
XSEC1D.o                       : XSEC1D.f                       cmn_fj.h CMN_SIZE define.h jv_cmn.h
XSECO2.o                       : XSECO2.f                       cmn_fj.h CMN_SIZE define.h jv_cmn.h
XSECO3.o                       : XSECO3.f                       cmn_fj.h CMN_SIZE define.h jv_cmn.h
aerosol_mod.o                  : aerosol_mod.f                  CMN_DIAG CMN_SIZE define.h cmn_fj.h CMN_SIZE define.h comode.h jv_cmn.h
aircraft_nox_mod.o             : aircraft_nox_mod.f             CMN CMN_DIAG CMN_SIZE define.h
airmas.o                       : airmas.f
anthroems.o                    : anthroems.f                    CMN_O3 CMN_SIZE define.h comode.h
arctas_ship_emiss_mod.o        : arctas_ship_emiss_mod.f        CMN_SIZE define.h
arsl1k.o                       : arsl1k.f
backsub.o                      : backsub.f                      CMN_SIZE define.h comode.h
benchmark_mod.o                : benchmark_mod.f                CMN_SIZE define.h
biofit.o                       : biofit.f                       CMN_DEP CMN_SIZE define.h
biofuel_mod.o                  : biofuel_mod.f                  CMN_DIAG CMN_O3 CMN_SIZE define.h
biomass_mod.o                  : biomass_mod.f                  CMN_DIAG CMN_SIZE define.h
boxvl.o                        : boxvl.f
bravo_mod.o                    : bravo_mod.f                    CMN_O3 CMN_SIZE define.h
c2h6_mod.o                     : c2h6_mod.f                     CMN CMN_DIAG CMN_O3 CMN_SIZE define.h
cac_anthro_mod.o               : cac_anthro_mod.f               CMN_O3 CMN_SIZE define.h
ch3i_mod.o                     : ch3i_mod.f                     CMN_DEP CMN_DIAG CMN_SIZE define.h comode.h
charpak_mod.o                  : charpak_mod.f
cleanup.o                      : cleanup.f
comode_mod.o                   : comode_mod.f                   CMN_SIZE define.h comode.h
decomp.o                       : decomp.f                       CMN_SIZE define.h comode.h
diag03_mod.o                   : diag03_mod.f                   CMN_DIAG CMN_SIZE define.h
diag04_mod.o                   : diag04_mod.f                   CMN_DIAG CMN_SIZE define.h
diag1.o                        : diag1.f                        CMN_DIAG CMN_GCTM CMN_O3 CMN_SIZE define.h
diag3.o                        : diag3.f                        CMN CMN_DIAG CMN_O3 CMN_SIZE define.h comode.h
diag41_mod.o                   : diag41_mod.f                   CMN_DIAG CMN_SIZE define.h
diag42_mod.o                   : diag42_mod.f                   CMN_DIAG CMN_SIZE define.h
diag48_mod.o                   : diag48_mod.f                   CMN_GCTM CMN_O3 CMN_SIZE define.h cmn_fj.h CMN_SIZE define.h jv_cmn.h
diag49_mod.o                   : diag49_mod.f                   CMN_GCTM CMN_O3 CMN_SIZE define.h cmn_fj.h CMN_SIZE define.h jv_cmn.h
diag50_mod.o                   : diag50_mod.f                   CMN_GCTM CMN_O3 CMN_SIZE define.h cmn_fj.h CMN_SIZE define.h jv_cmn.h
diag51_mod.o                   : diag51_mod.f                   CMN_GCTM CMN_O3 CMN_SIZE define.h cmn_fj.h CMN_SIZE define.h jv_cmn.h
diag56_mod.o                   : diag56_mod.f                   CMN_DIAG CMN_SIZE define.h
diag59_mod.o                   : diag59_mod.f                   CMN_DIAG CMN_SIZE define.h
diag_2pm.o                     : diag_2pm.f                     CMN_DIAG CMN_SIZE define.h
diag_mod.o                     : diag_mod.f
diag_oh_mod.o                  : diag_oh_mod.f                  CMN_SIZE define.h comode.h
diag_pl_mod.o                  : diag_pl_mod.f                  CMN_DIAG CMN_SIZE define.h comode.h
diagoh.o                       : diagoh.f                       CMN_DIAG CMN_O3 CMN_SIZE define.h
directory_mod.o                : directory_mod.f
drydep_mod.o                   : drydep_mod.f                   CMN_DEP CMN_DIAG CMN_GCTM CMN_SIZE define.h CMN_VEL commsoil.h comode.h
dust_dead_mod.o                : dust_dead_mod.f                CMN_GCTM CMN_SIZE define.h
edgar_mod.o                    : edgar_mod.f                    CMN_SIZE define.h
emep_mod.o                     : emep_mod.f                     CMN_O3 CMN_SIZE define.h
nei2005_anthro_mod.o           : nei2005_anthro_mod.f           CMN_O3 CMN_SIZE define.h
nei2008_anthro_mod.o           : nei2008_anthro_mod.F90         CMN_O3 CMN_SIZE define.h
htap_mod.o                     : htap_mod.f90                   CMN_O3 CMN_SIZE define.h
emf_scale.o                    : emf_scale.f                    CMN_O3 CMN_SIZE define.h comode.h
emfossil.o                     : emfossil.f                     CMN_DIAG CMN_O3 CMN_SIZE define.h comode.h rcp_mod.o
emisop.o                       : emisop.f                       CMN_ISOP CMN_SIZE define.h CMN_VEL
emisop_grass.o                 : emisop_grass.f                 CMN_ISOP CMN_SIZE define.h CMN_VEL
emisop_mb.o                    : emisop_mb.f                    CMN_ISOP CMN_SIZE define.h CMN_VEL
emissdr.o                      : emissdr.f                      CMN CMN_DIAG CMN_MONOT CMN_NOX CMN_O3 CMN_SIZE define.h comode.h
emissions_mod.o                : emissions_mod.f                CMN_O3 CMN_SIZE define.h rcp_mod.o
emmonot.o                      : emmonot.f                      CMN_MONOT CMN_SIZE define.h CMN_VEL
epa_nei_mod.o                  : epa_nei_mod.f                  CMN_O3 CMN_SIZE define.h
error_mod.o                    : error_mod.f                    define.h
fast_j.o                       : fast_j.f                       cmn_fj.h CMN_SIZE define.h jv_cmn.h
fertadd.o                      : fertadd.f                      CMN_SIZE define.h commsoil.h
file_mod.o                     : file_mod.f                     define.h
findmon.o                      : findmon.f
fjfunc.o                       : fjfunc.f                       cmn_fj.h CMN_SIZE define.h
future_emissions_mod.o         : future_emissions_mod.f         CMN_SIZE define.h
fvdas_convect_mod.o            : fvdas_convect_mod.f            CMN_DIAG CMN_SIZE define.h
fcro2ho2.o                     : fcro2ho2.f
fyrno3.o                       : fyrno3.f
fyhoro.o                       : fyhoro.f
gc_biomass_mod.o               : gc_biomass_mod.f               CMN_SIZE define.h
gcap_convect_mod.o             : gcap_convect_mod.f             CMN_DIAG CMN_SIZE define.h
gcap_read_mod.o                : gcap_read_mod.f                CMN_DIAG CMN_GCTM CMN_SIZE define.h
geia_mod.o                     : geia_mod.f                     CMN_SIZE define.h
get_global_ch4.o               : get_global_ch4.f
getifsun.o                     : getifsun.f                     CMN_SIZE define.h comode.h
gfed2_biomass_mod.o            : gfed2_biomass_mod.f            CMN_SIZE define.h
gfed3_biomass_mod.o            : gfed3_biomass_mod.f            CMN_SIZE define.h
global_hno3_mod.o              : global_hno3_mod.f              CMN_SIZE define.h
global_no3_mod.o               : global_no3_mod.f               CMN_SIZE define.h
global_nox_mod.o               : global_nox_mod.f               CMN_SIZE define.h
global_o1d_mod.o               : global_o1d_mod.f               CMN_SIZE define.h
global_o3_mod.o                : global_o3_mod.f                CMN_SIZE define.h
h2_hd_mod.o                    : h2_hd_mod.f                    CMN_DEP CMN_DIAG CMN_O3 CMN_SIZE define.h
hcn_ch3cn_mod.o                : hcn_ch3cn_mod.f                CMN_DEP CMN_DIAG CMN_SIZE define.h
icoads_ship_mod.o              : icoads_ship_mod.f              CMN_O3 CMN_SIZE define.h
ifort_errmsg.o                 : ifort_errmsg.f
initialize.o                   : initialize.f                   CMN_DIAG CMN_SIZE define.h
inphot.o                       : inphot.f                       cmn_fj.h CMN_SIZE define.h jv_cmn.h
rd_aod.o                       : rd_aod.f                       cmn_fj.h CMN_SIZE define.h jv_cmn.h
isoropiaII_adj_mod.o           : isoropiaII_adj_mod.f           new/isrpia_adj.inc
	$(F90) -c -r8 new/isoropiaII_adj_mod.f
isoropiaIIcode_adj.o           : isoropiaIIcode_adj.f           new/isrpia_adj.inc
	$(F90) -c -r8 new/isoropiaIIcode_adj.f
adBuffer.o                     : adBuffer.c
	$(CC) -c  new/adBuffer.c
adStack.o                      : adStack.c
	$(CC)  -c new/adStack.c
inquire_Mod.o                  : inquireMod.F90
jsparse.o                      : jsparse.f                      CMN_SIZE define.h comode.h
jv_index.o                     : jv_index.f                     cmn_fj.h CMN_SIZE define.h comode.h
ksparse.o                      : ksparse.f                      CMN_SIZE define.h comode.h
lai_mod.o                      : lai_mod.f                      CMN_SIZE define.h
logical_mod.o                  : logical_mod.f
lump.o                         : lump.f                         CMN_SIZE define.h comode.h
main.o                         : main.f                         CMN_DIAG CMN_GCTM CMN_SIZE define.h
mercury_mod.o                  : mercury_mod.f                  CMN_DIAG CMN_GCTM CMN_SIZE define.h
mmran_16.o                     : mmran_16.f                     cmn_fj.h CMN_SIZE define.h jv_cmn.h
ndxx_setup.o                   : ndxx_setup.f                   CMN_DIAG CMN_SIZE define.h
ocean_mercury_mod.o            : ocean_mercury_mod.f            CMN_DEP CMN_SIZE define.h
ohsave.o                       : ohsave.f                       CMN_SIZE define.h comode.h
optdepth_mod.o                 : optdepth_mod.f                 CMN_DIAG CMN_SIZE define.h
pderiv.o                       : pderiv.f                       CMN_SIZE define.h comode.h
pjc_pfix_mod.o                 : pjc_pfix_mod.f                 CMN CMN_GCTM CMN_SIZE define.h
planeflight_mod.o              : planeflight_mod.f              CMN_DIAG CMN_SIZE define.h cmn_fj.h CMN_SIZE define.h comode.h jv_cmn.h
precipfrac.o                   : precipfrac.f                   CMN_SIZE define.h
pulsing.o                      : pulsing.f                      CMN_SIZE define.h commsoil.h
rcp_mod.o                      : rcp_mod.f                      diag_mod.o aircraft_nox_mod.o dao_mod.o
rd_js.o                        : rd_js.f                        cmn_fj.h CMN_SIZE define.h jv_cmn.h
rd_prof.o                      : rd_prof.f                      cmn_fj.h CMN_SIZE define.h jv_cmn.h
rdisopt.o                      : rdisopt.f                      CMN_SIZE define.h
rdlai.o                        : rdlai.f                        CMN_DEP CMN_SIZE define.h CMN_VEL
rdland.o                       : rdland.f                       CMN_DEP CMN_SIZE define.h CMN_VEL
rdlight.o                      : rdlight.f                      CMN_ISOP CMN_SIZE define.h
rdmonot.o                      : rdmonot.f                      CMN_SIZE define.h
rdsoil.o                       : rdsoil.f                       CMN_SIZE define.h commsoil.h
readchem.o                     : readchem.f                     CMN_SIZE define.h comode.h
reader.o                       : reader.f                       CMN_GCTM CMN_SIZE define.h comode.h
readlai.o                      : readlai.f                      CMN_DEP CMN_SIZE define.h CMN_VEL
regrid_1x1_mod.o               : regrid_1x1_mod.f               CMN_GCTM CMN_SIZE define.h
retro_mod.o                    : retro_mod.f                    CMN CMN_SIZE CMN_O3
regrid_a2a_mod.o               : regrid_a2a_mod.F90             CMN_GCTM CMN_SIZE
ruralbox.o                     : ruralbox.f                     CMN_SIZE define.h comode.h
scale_anthro_mod.o             : scale_anthro_mod.f             CMN_SIZE define.h
schem.o                        : schem.f                        CMN_SIZE define.h
seasalt_mod.o                  : seasalt_mod.f                  CMN_DIAG CMN_GCTM CMN_SIZE define.h
set_aer.o                      : set_aer.f                      cmn_fj.h CMN_SIZE define.h jv_cmn.h
setbase.o                      : setbase.f                      CMN CMN_ISOP CMN_MONOT CMN_SIZE define.h CMN_VEL
setemdep.o                     : setemdep.f                     CMN_SIZE define.h comode.h
setmodel.o                     : setmodel.f                     CMN_SIZE define.h comode.h
sfcwindsqr.o                   : sfcwindsqr.f                   CMN_SIZE define.h
smvgear.o                      : smvgear.f                      CMN_SIZE define.h comode.h
soilbase.o                     : soilbase.f                     CMN_SIZE define.h commsoil.h
soilcrf.o                      : soilcrf.f                      CMN_DEP CMN_SIZE define.h commsoil.h
soiltemp.o                     : soiltemp.f                     CMN_SIZE define.h commsoil.h
soiltype.o                     : soiltype.f                     CMN_SIZE define.h commsoil.h
streets_anthro_mod.o           : streets_anthro_mod.f           CMN_O3 CMN_SIZE define.h
subfun.o                       : subfun.f                       CMN_SIZE define.h comode.h
sunparam.o                     : sunparam.f
tcorr.o                        : tcorr.f
toms_mod.o                     : toms_mod.f                     CMN_SIZE define.h
tpcore_bc_mod.o                : tpcore_bc_mod.f                CMN CMN_SIZE define.h
tpcore_fvdas_mod.o             : tpcore_fvdas_mod.f90             CMN_GCTM
	$(F90) -c -r8 $*.f90
tpcore_mod.o                   : tpcore_mod.f                   CMN_DIAG CMN_GCTM CMN_SIZE define.h define.h
	$(F90) -c -r8 $*.f
tpcore_window_mod.o            : tpcore_window_mod.f            CMN_DIAG CMN_GCTM CMN_SIZE define.h define.h
	$(F90) -c -r8 $*.f
tracerid_mod.o                 : tracerid_mod.f                 CMN_SIZE define.h comode.h
transfer_mod.o                 : transfer_mod.f                 CMN_SIZE define.h
tropopause.o                   : tropopause.f                   CMN CMN_DIAG CMN_SIZE define.h
tropopause_mod.o               : tropopause_mod.f               CMN CMN_SIZE define.h comode.h
unix_cmds_mod.o                : unix_cmds_mod.f
update.o                       : update.f                       CMN_SIZE define.h comode.h
uvalbedo_mod.o                 : uvalbedo_mod.f                 CMN_SIZE define.h
vistas_anthro_mod.o            : vistas_anthro_mod.f            CMN_O3 CMN_SIZE define.h
xltmmp.o                       : xltmmp.f                       CMN_SIZE define.h
xtra_read_mod.o                : xtra_read_mod.f                CMN_DIAG CMN_SIZE define.h

#------------------------------------------------------------------------------
#	          MODIFIED FOR ORGANIZED DIRECTORY STRUCTURE
#------------------------------------------------------------------------------

#========================
# ADJOINT DIRECTORY FILES
#========================

adj_arrays_mod.o               : adj_arrays_mod.f               define_adj.h define.h CMN_SIZE comode.h gckpp_adj_Global.f90
	$(F90) -c -r8 adjoint/adj_arrays_mod.f
calcrate_adj.o                 : calcrate_adj.f                 CMN CMN_DIAG CMN_SIZE define.h comode.h
	$(F90) -c -r8 adjoint/calcrate_adj.f
carbon_adj_mod.o               : carbon_adj_mod.f               CMN CMN_DIAG CMN_GCTM CMN_O3 CMN_SIZE define.h comode.h
	$(F90) -c -r8 adjoint/carbon_adj_mod.f
checkpoint_mod.o               : checkpoint_mod.f               define.h
	$(F90) -c -r8 adjoint/checkpoint_mod.f
checkpt_mod.o                  : checkpt_mod.f                  CMN_SIZE comode.h CMN_VEL define.h CMN_DEP
	$(F90) -c -r8 adjoint/checkpt_mod.f
chemdr_adj.o                   : chemdr_adj.f                   CMN CMN_DEP CMN_DIAG CMN_NOX CMN_O3 CMN_SIZE define.h comode.h
	$(F90) -c -r8 adjoint/chemdr_adj.f
chemistry_adj_mod.o            : chemistry_adj_mod.f            gckpp_adj_Global.f90 gckpp_adj_Rates.f90 gckpp_adj_Integrator.f90 CMN_DIAG CMN_SIZE define.h comode.h
	$(F90) -c -r8 adjoint/chemistry_adj_mod.f
cleanup_adj.o                  : cleanup_adj.f
	$(F90) -c -r8 adjoint/cleanup_adj.f
co2_adj_mod.o                  : co2_adj_mod.f                  CMN_SIZE define.h
	$(F90) -c -r8 adjoint/co2_adj_mod.f
CO_strat_pl_adj.o              : CO_strat_pl_adj.f              CMN_SIZE define.h
	$(F90) -c -r8 adjoint/CO_strat_pl_adj.f
convection_adj_mod.o           : convection_adj_mod.f           CMN_DIAG CMN_SIZE define.h define.h
	$(F90) -c -r8 adjoint/convection_adj_mod.f
directory_adj_mod.o            : directory_adj_mod.f
	$(F90) -c -r8 adjoint/directory_adj_mod.f
dust_adj_mod.o                 : dust_adj_mod.f                 CMN_DIAG CMN_GCTM CMN_SIZE define.h define_adj.h
	$(F90) -c -r8 adjoint/dust_adj_mod.f
emissions_adj_mod.o            : emissions_adj_mod.f            CMN_O3 CMN_SIZE define.h
	$(F90) -c -r8 adjoint/emissions_adj_mod.f
fvdas_convect_adj_mod.o        : fvdas_convect_adj_mod.f        CMN_DIAG CMN_SIZE define.h
	$(F90) -c -r8 adjoint/fvdas_convect_adj_mod.f
gc_type_mod.o                  : gc_type_mod.F                  CMN_GCTM CMN_SIZE define.h
	$(F90) -c -r8 modified/gc_type_mod.F
global_ch4_adj_mod.o           : global_ch4_adj_mod.f           CMN CMN_GCTM CMN_DIAG CMN_SIZE define.h
	$(F90) -c -r8 $<
geos_chem_adj_mod.o            : geos_chem_adj_mod.f        CMN_O3 CMN_GCTM CMN_DIAG CMN_SIZE define.h comode.h define_adj.h
	$(F90) -c -r8 adjoint/geos_chem_adj_mod.f
covariance_mod.o	       : covariance_mod.f           CMN_SIZE CMN_DIAG CMN_GCTM define_adj.h
	$(F90) -c -r8 adjoint/covariance_mod.f
input_adj_mod.o                : input_adj_mod.f                CMN CMN_DIAG CMN_O3 CMN_SIZE define.h define_adj.h gckpp_adj_Global.f90
	$(F90) -c -r8 adjoint/input_adj_mod.f
inv_hessian_mod.o              : inv_hessian_mod.f              CMN_SIZE define_adj.h
	$(F90) -c -r8 adjoint/inv_hessian_mod.f
inv_hessian_lbfgs_mod.o        : inv_hessian_lbfgs_mod.f              CMN_SIZE define_adj.h
	$(F90) $(LAPACK_BLAS_FFLAGS) -c -r8 adjoint/inv_hessian_lbfgs_mod.f
inverse_driver.o               : inverse_driver.f               define_adj.h
	$(F90) -c -r8 adjoint/inverse_driver.f
inverse_mod.o                  : inverse_mod.f	      	   define_adj.h define.h
	$(F90) -c -r8 adjoint/inverse_mod.f
linoz_adj_mod.o                : linoz_adj_mod.f                CMN CMN_DIAG CMN_O3 CMN_SIZE define.h linoz.com
	$(F90) -c -r8 adjoint/linoz_adj_mod.f
logical_adj_mod.o              : logical_adj_mod.f
	$(F90) -c -r8 adjoint/logical_adj_mod.f
lump_adj.o                     : lump_adj.f                     CMN_SIZE define.h comode.h
	$(F90) -c -r8 adjoint/lump_adj.f
gckpp_adj_Precision.o	       : gckpp_adj_Precision.f90
	$(F90) -c -r8 adjoint/gckpp_adj_Precision.f90
gckpp_adj_Parameters.o	       : gckpp_adj_Parameters.f90 gckpp_adj_Precision.f90
	$(F90) -c -r8 adjoint/gckpp_adj_Parameters.f90
gckpp_adj_Global.o	       : gckpp_adj_Global.f90
	$(F90) -c -r8 adjoint/gckpp_adj_Global.f90
gckpp_adj_LinearAlgebra.o      : gckpp_adj_LinearAlgebra.f90 gckpp_adj_Parameters.f90 gckpp_adj_JacobianSP.f90
	$(F90) -c -r8 adjoint/gckpp_adj_LinearAlgebra.f90
gckpp_adj_Monitor.o	       : gckpp_adj_Monitor.f90
	$(F90) -c -r8 adjoint/gckpp_adj_Monitor.f90
gckpp_adj_Initialize.o         : gckpp_adj_Initialize.f90 gckpp_adj_Parameters.f90 gckpp_adj_Global.f90 gckpp_adj_Util.f90
	$(F90) -c -r8 adjoint/gckpp_adj_Initialize.f90
gckpp_adj_JacobianSP.o         : gckpp_adj_JacobianSP.f90
	$(F90) -c -r8 adjoint/gckpp_adj_JacobianSP.f90
gckpp_adj_Function.o           : gckpp_adj_Function.f90 gckpp_adj_Parameters.f90
	$(F90) -c -r8 adjoint/gckpp_adj_Function.f90
gckpp_adj_Jacobian.o           : gckpp_adj_Jacobian.f90 gckpp_adj_Parameters.f90 gckpp_adj_JacobianSP.f90
	$(F90) -c -r8 adjoint/gckpp_adj_Jacobian.f90
gckpp_adj_HessianSP.o          : gckpp_adj_HessianSP.f90
	$(F90) -c -r8 adjoint/gckpp_adj_HessianSP.f90
gckpp_adj_Hessian.o            : gckpp_adj_Hessian.f90 gckpp_adj_Parameters.f90 gckpp_adj_HessianSP.f90
	$(F90) -c -r8 adjoint/gckpp_adj_Hessian.f90
gckpp_adj_Util.o               : gckpp_adj_Util.f90 gckpp_adj_Parameters.f90 gckpp_adj_Global.f90
	$(F90) -c -r8 adjoint/gckpp_adj_Util.f90
gckpp_adj_StoichiomSP.o        : gckpp_adj_StoichiomSP.f90
	$(F90) -c -r8 adjoint/gckpp_adj_StoichiomSP.f90
gckpp_adj_Stoichiom.o	       : gckpp_adj_Stoichiom.f90 gckpp_adj_Parameters.f90 gckpp_adj_StoichiomSP.f90
	$(F90) -c -r8 adjoint/gckpp_adj_Stoichiom.f90
gckpp_adj_Rates.o	       : gckpp_adj_Rates.f90 gckpp_adj_Parameters.f90 gckpp_adj_Global.f90
	$(F90) -c -r8 adjoint/gckpp_adj_Rates.f90
gckpp_adj_Model.o	       : gckpp_adj_Model.f90
	$(F90) -c -r8 adjoint/gckpp_adj_Model.f90
gckpp_adj_Integrator.o	       : gckpp_adj_Integrator.f90 gckpp_adj_Parameters.f90 gckpp_adj_Precision.f90 gckpp_adj_Global.f90 gckpp_adj_LinearAlgebra.f90 gckpp_adj_Rates.f90 gckpp_adj_Function.f90 gckpp_adj_Jacobian.f90 gckpp_adj_Hessian.f90 gckpp_adj_Util.f90
	$(F90) -c -r8 adjoint/gckpp_adj_Integrator.f90
partition_adj.o                : partition_adj.f                CMN_SIZE define.h comode.h
	$(F90) -c -r8 adjoint/partition_adj.f
pbl_mix_adj_mod.o              : pbl_mix_adj_mod.f              CMN_DIAG CMN_GCTM CMN_SIZE define.h
	$(F90) -c -r8 adjoint/pbl_mix_adj_mod.f
rpmares_adj_mod.o              : rpmares_adj_mod.f              CMN_SIZE define.h
	$(F90) -c -extend_source -r8 adjoint/rpmares_adj_mod.f
schem_adj.o                    : schem_adj.f                    CMN_SIZE define.h
	$(F90) -c -r8 adjoint/schem_adj.f
setemis_adj.o                  : setemis_adj.f                  CMN_O3   CMN_NOX CMN_SIZE define.h comode.h
	$(F90) -c -r8 adjoint/setemis_adj.f
sulfate_adj_mod.o              : sulfate_adj_mod.f              CMN_GCTM CMN_O3 CMN_SIZE define.h cmn_fj.h CMN_SIZE define.h
	$(F90) -c -r8 adjoint/sulfate_adj_mod.f
strat_chem_adj_mod.o           : strat_chem_adj_mod.f           CMN_SIZE
	$(F90) -c -r8 adjoint/strat_chem_adj_mod.f
tagged_co_adj_mod.o            : tagged_co_adj_mod.f            CMN_DIAG CMN_O3 CMN_SIZE define.h
	$(F90) -c -r8 adjoint/tagged_co_adj_mod.f
tagged_ox_adj_mod.o            : tagged_ox_adj_mod.f            CMN_DIAG CMN_O3 CMN_SIZE define.h
	$(F90) -c -r8 adjoint/tagged_ox_adj_mod.f
upbdflx_adj_mod.o              : upbdflx_adj_mod.f              CMN_SIZE define.h
	$(F90) -c -r8 adjoint/upbdflx_adj_mod.f
wetscav_adj_mod.o              : wetscav_adj_mod.f              CMN_DIAG CMN_SIZE define.h define_adj.h
	$(F90) -c -r8 adjoint/wetscav_adj_mod.f
weak_constraint_mod.o          : weak_constraint_mod.f90        CMN_SIZE
	$(F90) -c -r8 adjoint/weak_constraint_mod.f90

#=========================
# MODIFIED DIRECTORY FILES
#=========================
a3_read_mod.o                  : a3_read_mod.f                  CMN_DIAG CMN_SIZE define.h
	$(F90) -c -r8 modified/a3_read_mod.f
a6_read_mod.o                  : a6_read_mod.f                  CMN_DIAG CMN_GCTM CMN_SIZE define.h
	$(F90) -c -r8 modified/a6_read_mod.f
acetone_mod.o                  : acetone_mod.f                  CMN_DEP CMN_DIAG CMN_MONOT CMN_SIZE define.h
	$(F90) -c -r8 modified/acetone_mod.f
bpch2_mod.o                    : bpch2_mod.f                    CMN_SIZE define.h define.h
	$(F90) -c -r8 modified/bpch2_mod.f
calcrate.o                     : calcrate.f                     CMN CMN_DIAG CMN_SIZE define.h comode.h
	$(F90) -c -r8 modified/calcrate.f
carbon_mod.o                   : carbon_mod.f                   CMN CMN_DIAG CMN_GCTM CMN_O3 CMN_SIZE define.h comode.h
	$(F90) -c -r8 modified/carbon_mod.f
chemdr.o                       : chemdr.f                       CMN CMN_DEP CMN_DIAG CMN_NOX CMN_O3 CMN_SIZE define.h comode.h define_adj.h
	$(F90) -c -r8 modified/chemdr.f
chemistry_mod.o                : chemistry_mod.f                gckpp_adj_Global.f90 gckpp_adj_Rates.f90 gckpp_adj_Integrator.f90 CMN_DIAG CMN_SIZE define.h comode.h
	$(F90) -c -r8 modified/chemistry_mod.f
co2_mod.o                      : co2_mod.f                      CMN_SIZE define.h
	$(F90) -c -r8 modified/co2_mod.f
comode_mod.o                   : comode_mod.f                   CMN_SIZE define.h comode.h
	$(F90) -c -r8 modified/comode_mod.f
convection_mod.o               : convection_mod.f               CMN_DIAG CMN_SIZE define.h
	$(F90) -c -r8 modified/convection_mod.f
dao_mod.o                      : dao_mod.f                      CMN_GCTM CMN_SIZE define.h
	$(F90) -c -r8 modified/dao_mod.f
dust_mod.o                     : dust_mod.f                     CMN_DIAG CMN_GCTM CMN_SIZE define.h cmn_fj.h comode.h jv_cmn.h
	$(F90) -c -r8 modified/dust_mod.f
geosfp_read_mod.o              : geosfp_read_mod.f              CMN_DIAG CMN_SIZE CMN_GCTM define.h
	$(F90) -c -r8 modified/geosfp_read_mod.f
geos_chem_mod.o                : geos_chem_mod.f                CMN_SIZE CMN_DIAG CMN_GCTM comode.h define.h define_adj.h
	$(F90) -c -r8 modified/geos_chem_mod.f
gamap_mod.o                    : gamap_mod.f                    CMN_DIAG CMN_SIZE define.h
	$(F90) -c     modified/gamap_mod.f
gasconc.o                      : gasconc.f                      CMN_SIZE define.h comode.h
	$(F90) -c -r8 modified/gasconc.f
global_oh_mod.o                : global_oh_mod.f                CMN_SIZE define.h
	$(F90) -c -r8 modified/global_oh_mod.f
global_ch4_mod.o               : global_ch4_mod.f               CMN CMN_DIAG CMN_SIZE define.h
	$(F90) -c -r8 modified/global_ch4_mod.f
grid_mod.o                     : grid_mod.f                     CMN_GCTM CMN_SIZE define.h
	$(F90) -c -r8 modified/grid_mod.f
gwet_read_mod.o                : gwet_read_mod.f                CMN_DIAG CMN_SIZE define.h
	$(F90) -c -r8 modified/gwet_read_mod.f
i6_read_mod.o                  : i6_read_mod.f                  CMN_DIAG CMN_SIZE define.h
	$(F90) -c -r8 modified/i6_read_mod.f
julday_mod.o                   : julday_mod.f
	$(F90) -c -r8 modified/julday_mod.f
input_mod.o                    : input_mod.f                    CMN CMN_DIAG CMN_O3 CMN_SIZE define.h diag59_mod.f
	$(F90) -c -r8 modified/input_mod.f
lightning_nox_mod.o            : lightning_nox_mod.f            CMN_DIAG CMN_GCTM CMN_SIZE define.h define.h
	$(F90) -c -r8 modified/lightning_nox_mod.f
megan_mod.o                    : megan_mod.f                    CMN_GCTM CMN_SIZE define.h
	$(F90) -c -r8 modified/megan_mod.f
merra2_read_mod.o              : merra2_read_mod.f              CMN_DIAG CMN_SIZE CMN_GCTM define.h
	$(F90) -c -r8 modified/merra2_read_mod.f
paranox_mod.o                  : paranox_mod.f                  CMN_SIZE define.h comode.h
	$(F90) -c -r8 paranox_mod.f
paranox_adj_mod.o                  : paranox_adj_mod.f          CMN_SIZE define.h comode.h
	$(F90) -c -r8 adjoint/paranox_adj_mod.f
partition.o                    : partition.f                    CMN_SIZE define.h comode.h
	$(F90) -c -r8 modified/partition.f
pbl_mix_mod.o                  : pbl_mix_mod.f                  CMN_DIAG CMN_GCTM CMN_SIZE define.h
	$(F90) -c -r8 modified/pbl_mix_mod.f
photoj.o                       : photoj.f                       define_adj.h cmn_fj.h CMN_SIZE define.h jv_cmn.h
	$(F90) -c -r8 modified/photoj.f
physproc.o                     : physproc.f                     CMN_SIZE define.h comode.h
	$(F90) -c -r8 modified/physproc.f
pjc_pfix_geos5_window_mod.o    : pjc_pfix_geos5_window_mod.f    CMN CMN_GCTM CMN_SIZE define.h
	$(F90) -c     modified/pjc_pfix_geos5_window_mod.f
pjc_pfix_geosfp_window_mod.o   : pjc_pfix_geosfp_window_mod.f    CMN CMN_GCTM CMN_SIZE define.h
	$(F90) -c     modified/pjc_pfix_geosfp_window_mod.f
pjc_pfix_merra2_window_mod.o   : pjc_pfix_merra2_window_mod.f    CMN CMN_GCTM CMN_SIZE define.h
	$(F90) -c     modified/pjc_pfix_merra2_window_mod.f
pressure_mod.o                 : pressure_mod.f                 CMN_SIZE define.h
	$(F90) -c -r8 modified/pressure_mod.f
rpmares_mod.o                  : rpmares_mod.f                  CMN_SIZE define.h
	$(F90) -c -r8 modified/rpmares_mod.f
restart_mod.o                  : restart_mod.f                  CMN_SIZE define.h
	$(F90) -c -r8 modified/restart_mod.f
set_prof.o                     : set_prof.f                     define_adj.h cmn_fj.h CMN_SIZE define.h jv_cmn.h
	$(F90) -c -r8 modified/set_prof.f
setemis.o                      : setemis.f                      CMN_DIAG CMN_NOX CMN_SIZE define.h comode.h
	$(F90) -c -r8 modified/setemis.f
soilnoxems.o                   : soilnoxems.f                   CMN_DEP CMN_DIAG CMN_NOX CMN_SIZE define.h commsoil.h
	$(F90) -c -r8 modified/soilnoxems.f
sulfate_mod.o                  : sulfate_mod.f                  CMN_DIAG CMN_GCTM CMN_O3 CMN_SIZE define.h cmn_fj.h CMN_SIZE define.h rcp_mod.o
	$(F90) -c -r8 modified/sulfate_mod.f
tagged_co_mod.o                : tagged_co_mod.f                CMN_DIAG CMN_O3 CMN_SIZE define.h
	$(F90) -c -r8 modified/tagged_co_mod.f
tagged_ox_mod.o                : tagged_ox_mod.f                CMN_DIAG CMN_GCTM CMN_SIZE define.h
	$(F90) -c -r8 modified/tagged_ox_mod.f
tpcore_geos5_window_mod.o      : tpcore_geos5_window_mod.f90
	$(F90) -c -r8 modified/tpcore_geos5_window_mod.f90
tpcore_geosfp_window_mod.o     : tpcore_geosfp_window_mod.f90
	$(F90) -c -r8 modified/tpcore_geosfp_window_mod.f90
tpcore_merra2_window_mod.o     : tpcore_merra2_window_mod.f90
	$(F90) -c -r8 modified/tpcore_merra2_window_mod.f90
time_mod.o                     : time_mod.f                     define.h
	$(F90) -c -r8 modified/time_mod.f
tracer_mod.o                   : tracer_mod.f                   CMN_SIZE define.h
	$(F90) -c -r8 modified/tracer_mod.f
transport_mod.o                : transport_mod.f                CMN CMN_DIAG CMN_GCTM CMN_SIZE define.h
	$(F90) -c -r8 modified/transport_mod.f
upbdflx_mod.o                  : upbdflx_mod.f                  CMN_GCTM CMN_SIZE define.h
	$(F90) -c -r8 modified/upbdflx_mod.f
wetscav_mod.o                  : wetscav_mod.f                  CMN_DIAG CMN_SIZE define.h
	$(F90) -c -r8 modified/wetscav_mod.f

#====================
# NEW DIRECTORY FILES
#====================

cgfam.o                        : cgfam.f
	$(F90) -c -r8 new/cgfam.f
cgsearch.o                     : cgsearch.f
	$(F90) -c -r8 new/cgsearch.f
linoz_mod.o                    : linoz_mod.f                    CMN CMN_DIAG CMN_O3 CMN_SIZE define.h linoz.com
	$(F90) -c -r8 new/linoz_mod.f
routines.o                     : routines.f
	$(F90) -c -r8 new/routines.f
blas.o                         : new/blas.f
	$(F90) -c -r8 new/blas.f
linpack.o                      : new/linpack.f
	$(F90) -c -r8 new/linpack.f
timer.o                        : new/timer.f
	$(F90) -c -r8 new/timer.f
netcdf_util_mod.o              : netcdf_util_mod.f
	$(F90) -c -r8 new/netcdf_util_mod.f
strat_chem_mod.o               : strat_chem_mod.f               CMN_DIAG CMN_SIZE define.h comode.h
	$(F90) -c -r8 new/strat_chem_mod.f

#====================
# NETCDF directory
#====================

m_do_err_out.o                : m_do_err_out.F90
	$(F90) -c -r8 NcdfUtil/m_do_err_out.F90
m_netcdf_io_checks.o		: m_netcdf_io_checks.F90
	$(F90) -c -r8 NcdfUtil/m_netcdf_io_checks.F90
m_netcdf_io_close.o		: m_netcdf_io_close.F90
	$(F90) -c -r8 NcdfUtil/m_netcdf_io_close.F90
m_netcdf_io_create.o     	: m_netcdf_io_create.F90
	$(F90) -c -r8 NcdfUtil/m_netcdf_io_create.F90
m_netcdf_io_define.o		: m_netcdf_io_define.F90
	$(F90) -c -r8 NcdfUtil/m_netcdf_io_define.F90
m_netcdf_io_get_dimlen.o 	: m_netcdf_io_get_dimlen.F90
	$(F90) -c -r8 NcdfUtil/m_netcdf_io_get_dimlen.F90
m_netcdf_io_handle_err.o	: m_netcdf_io_handle_err.F90
	$(F90) -c -r8 NcdfUtil/m_netcdf_io_handle_err.F90
m_netcdf_io_open.o       	: m_netcdf_io_open.F90
	$(F90) -c -r8 NcdfUtil/m_netcdf_io_open.F90
m_netcdf_io_read.o       	: m_netcdf_io_read.F90
	$(F90) -c -r8 NcdfUtil/m_netcdf_io_read.F90
m_netcdf_io_readattr.o		: m_netcdf_io_readattr.F90
	$(F90) -c -r8 NcdfUtil/m_netcdf_io_readattr.F90
m_netcdf_io_write.o      	: m_netcdf_io_write.F90
	$(F90) -c -r8 NcdfUtil/m_netcdf_io_write.F90


#====================
# SCIA CODE
#====================
ErrorModule.o        : ErrorModule.f90
	$(F90) -c -r8 obs_operators/ErrorModule.f90
sciabr_co_obs_mod.o  : sciabr_co_obs_mod.f CMN_SIZE
	$(F90) -c -r8 obs_operators/sciabr_co_obs_mod.f

#====================
# TES CODE
#====================
tes_ch4_mod.o        : tes_ch4_mod.f
	$(F90) -c -r8 $<

#====================
# Other CH4 obs operators
#====================
mem_ch4_mod.o        : mem_ch4_mod.f CMN_SIZE
	$(F90) -c -r8 $<
geocape_ch4_mod.o    : geocape_ch4_mod.f CMN_SIZE
	$(F90) -c -r8 $<
leo_ch4_mod.o        : leo_ch4_mod.f CMN_SIZE
	$(F90) -c -r8 $<

#====================
# Other
#====================
improve_bc_mod.o     : improve_bc_mod.f                        CMN_SIZE
	$(F90) -c -r8 obs_operators/improve_bc_mod.f
population_mod.o     : population_mod.f CMN_SIZE
	$(F90) -c -r8 obs_operators/population_mod.f
