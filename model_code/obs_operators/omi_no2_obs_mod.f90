MODULE OMI_NO2_OBS_MOD

  !
  ! Module OMI_NO2_OBS contains all subroutines and variables needed to assimilate OMI NO2 tropospheric column data
  !
  ! Module Routines:
  !
  ! (1) CALC_OMI_NO2_FORCE      : calculates adjoint forcing and cost function contribution for OMI tropospheric NO2 columns
  ! (2) TAI2UTC                 : converts TAI93 (seconds since 1.1.1993) to UTC
  ! (3) MAKE_OMI_BIAS_FILE_HDF5 : writes OMI satellite diagnostics in satellite diagnostic HDF5 file
  !

  IMPLICIT NONE

#include "CMN_SIZE"

  PRIVATE

  PUBLIC CALC_OMI_NO2_FORCE
  PUBLIC MAKE_OMI_BIAS_FILE_HDf5

  ! Module variables

  ! Arrays for diagnostic output

  TYPE FLEX_REAL                                ! Type to store information for "flexible" arrays. Think of this as a cruddy
     INTEGER :: CURRENT_N, MAX_N                ! implementation of some of the features of the C++ std::vector<T> container
     REAL*8,ALLOCATABLE :: DATA(:)              ! This only works in Fortran 2003, would have to use a pointer in Fortran 95
  ENDTYPE FLEX_REAL

  ! arrays to store diagnostic information
  REAL*4:: OMI_NO2_MEAN(IIPAR,JJPAR) = 0d0      ! Mean OMI columns
  REAL*4:: OMI_GEOS_NO2_MEAN(IIPAR,JJPAR) = 0d0 ! Mean GEOS-Chem columns
  REAL*4:: OMI_NO2_ERR_MEAN(IIPAR,JJPAR) = 0d0  ! Mean OMI observation errors
  REAL*4:: OMI_BIAS(IIPAR,JJPAR)=0d0            ! Model biases
  REAL*4:: OMI_VAR(IIPAR,JJPAR)=0d0             ! Model variances
  REAL*4:: OMI_DELTA=0d0                        ! temporary storage variable
  REAL*4:: OMI_BIAS_COUNT(IIPAR,JJPAR) = 0d0    ! counter for number of observations in grid box
  REAL*4:: OMI_CHISQUARED(IIPAR,JJPAR) = 0d0   ! Chi-squared values
  LOGICAL :: FIRST = .TRUE.
  TYPE(FLEX_REAL) :: FLEX_LON, FLEX_LAT, FLEX_TIME, FLEX_OMI_NO2, FLEX_GC_NO2 ! flex arrays to store satellite diagnostics sequentially

CONTAINS

  !-----------------------------------------------------------------------------!

  SUBROUTINE CALC_OMI_NO2_FORCE

    USE HDF5

    !!
    !! Subroutine CALC_OMI_NO2_FORCE computes the NO2 adjoint forcing and cost function contribution from OMI column data
    !!
    !! References:
    !!
    !! Bucsela2013:
    !! "A new stratospheric and tropospheric NO2 retrieval algorithm for nadir-viewing satellite instruments: applications to OMI"
    !! E.J. Bucsela et.al
    !! Atmos. Meas. Tech., 6, 2607-2626, 2013
    !! www.atmos-meas-tech.net/6/2607/2013/
    !! doi:10.5194/amt-6-2607-2013
    !!
    !! Chan83
    !! "Algorithms for Computing the Sample Variance: Analysis and Recommendations"
    !! Tony F. Chan, Gene H. Golub, Randall J. LeVeque
    !! The American Statistician
    !! Vol. 37, No. 3 (Aug. 1983), pp. 242-247
    !!

    USE COMODE_MOD,         ONLY : JLOP
    USE COMODE_MOD,         ONLY : CSPEC_AFTER_CHEM
    USE COMODE_MOD,         ONLY : CSPEC_AFTER_CHEM_ADJ
    USE GRID_MOD,           ONLY : GET_IJ
    USE GRID_MOD,           ONLY : GET_XMID, GET_YMID
    USE TIME_MOD,           ONLY : GET_HOUR, GET_DAY
    USE DAO_MOD,            ONLY : BXHEIGHT
    USE TROPOPAUSE_MOD,     ONLY : ITS_IN_THE_TROP
    USE ADJ_ARRAYS_MOD,     ONLY : ID2C
    USE TRACERID_MOD,       ONLY : IDNO2
    USE ADJ_ARRAYS_MOD,     ONLY : COST_FUNC
    USE PRESSURE_MOD,       ONLY : GET_PCENTER, GET_PEDGE
    USE TRACER_MOD,         ONlY : XNUMOLAIR
    USE DAO_MOD,            ONLY : T, AIRDEN
    USE LOGICAL_ADJ_MOD,    ONLY : LSAT_HDF_L2, LSAT_HDF_L3

    INTEGER :: I,J,L
    INTEGER :: I_OMI, J_OMI, K_OMI, JLOOP
    INTEGER :: IIJJ(2)
    INTEGER :: DAY
    CHARACTER(255) :: ORBIT_PATH,FILE_ORBIT
    CHARACTER(2) :: I_CHAR
    INTEGER :: IO_ORBIT_STATUS

    !! HDF5 variables for orbit OMI NO2 data

    INTEGER(HID_T) :: file_id_orbit
    INTEGER(HID_T) :: dset_id_orbit
    INTEGER(HID_T) :: dspace_id_orbit
    INTEGER :: error_orbit
    CHARACTER(LEN=255) :: filename_orbit, dsetname

    REAL*8, ALLOCATABLE :: LON_ORBIT(:,:), LAT_ORBIT(:,:)
    REAL*8, ALLOCATABLE :: TIME_ORBIT(:), AMF_TROP_ORBIT(:,:)
    REAL*8, ALLOCATABLE :: NO2_TROP(:,:), NO2_TROP_STD(:,:)
    REAL*8, ALLOCATABLE :: VIEW_ZENITH(:,:), SOLAR_ZENITH(:,:)
    REAL*8, ALLOCATABLE :: CLOUDFR(:,:), SCW_P(:)
    REAL*8, ALLOCATABLE :: SCATTERING_WEIGHTS(:,:,:)

    INTEGER :: rank_orbit
    INTEGER(HSIZE_T) :: dims_orbit(3), maxdims_orbit(3)
    INTEGER(HSIZE_T) :: data_dims_orbit(3)

    ! variables for time unit conversion
    REAL*8 :: tai93
    INTEGER :: iy,im,id,ih,imin
    REAL*8 :: sec
    INTEGER :: GC_HOUR, MIN_HOUR, MAX_HOUR

    ! variables for observation operator and adjoint thereof

    REAL*8 :: OMI_NO2_GC(IIPAR,JJPAR)
    REAL*8 :: SCW_GC(LLPAR), DP(LLPAR)
    REAL*8 :: AMF_GC
    REAL*8 :: GC_NO2(LLPAR)
    REAL*8 :: GC_NO2_COL
    REAL*8 :: DIFF
    REAL*8 :: OBS_ERROR

    ! arrays needed for superobservations
    LOGICAL :: SUPER_OBS = .TRUE.                       ! do super observations?
    REAL*8 ::  SOBS_COUNT(IIPAR,JJPAR)                   ! super observation count
    REAL*8 ::  SOBS_ADJ_FORCE(IIPAR,JJPAR,LLPAR)         ! super observation adjoint forcing
    REAL*8 ::  SOBS_COST_CONTRIBUTION(IIPAR,JJPAR)       ! super observation cost function contribution
    REAL*8 ::  SOBS_GC(IIPAR,JJPAR)
    REAL*8 ::  SOBS_OMI(IIPAR,JJPAR)
    REAL*8 ::  SOBS_BIAS(IIPAR,JJPAR)
    REAL*8 ::  SOBS_CHISQUARED(IIPAR,JJPAR)

    !=================================================================
    ! CALC_OMI_NO2_FORCE begins here!
    !=================================================================

    IF(LSAT_HDF_L2 .AND. FIRST) THEN

       ! initialize flexible arrays
       CALL INIT_FLEX_REAL(FLEX_LON)
       CALL INIT_FLEX_REAL(FLEX_LAT)
       CALL INIT_FLEX_REAL(FLEX_TIME)
       CALL INIT_FLEX_REAL(FLEX_OMI_NO2)
       CALL INIT_FLEX_REAL(FLEX_GC_NO2)
       
       FIRST = .FALSE.
       
    ENDIF

    ! initialize arrays

    GC_NO2 = 0d0
    GC_NO2_COL = 0d0
    OMI_NO2_GC = 0d0

    SOBS_COUNT = 0d0
    SOBS_ADJ_FORCE = 0d0
    SOBS_COST_CONTRIBUTION = 0d0
    SOBS_GC = 0d0
    SOBS_OMI = 0d0
    SOBS_BIAS = 0d0
    SOBS_CHISQUARED = 0d0

    ! Loop through data to find observations

    GC_HOUR = GET_HOUR()
    DAY = GET_DAY()

    ORBIT_PATH = '/met/gc/OMI_NO2/'

    WRITE(I_CHAR,'(I2.2)') DAY

    CALL SYSTEM("ls "//TRIM(ORBIT_PATH)//"OMI-Aura_L2-OMNO2_2011"//I_CHAR//"* > omi_file_list"//I_CHAR//".txt")

    CLOSE(65) ! ugly...

    OPEN(65,FILE="omi_file_list"//I_CHAR//".txt",ACTION="read",ACCESS="sequential",FORM="FORMATTED")

    DO ! loop over all available OMI NO2 files for the current day

       READ(65,'(A)',IOSTAT=IO_ORBIT_STATUS) FILE_ORBIT

       IF(IO_ORBIT_STATUS < 0) EXIT

       !FILE_ORBIT = TRIM(ORBIT_PATH) // FILE_ORBIT

       PRINT *,"Reading OMI file "//TRIM(FILE_ORBIT)

       !! open OMI orbit file

       CALL H5OPEN_F(error_orbit)

       CALL H5FOPEN_f (FILE_ORBIT, H5F_ACC_RDWR_F,file_id_orbit,error_orbit)

       PRINT *,"OMI file status: ",error_orbit

       ! Open an existing dataset.

       DSETNAME = '/HDFEOS/SWATHS/ColumnAmountNO2/Data Fields/ScatteringWeight'

       CALL H5DOPEN_F(FILE_ID_ORBIT, DSETNAME, DSET_ID_ORBIT, ERROR_ORBIT)

       ! open dataspace

       CALL H5DGET_SPACE_F(DSET_ID_ORBIT, DSPACE_ID_ORBIT, ERROR_ORBIT)

       CALL H5SGET_SIMPLE_EXTENT_NDIMS_F(DSPACE_ID_ORBIT, RANK_ORBIT,ERROR_ORBIT)

       CALL H5SGET_SIMPLE_EXTENT_DIMS_F(DSPACE_ID_ORBIT, DIMS_ORBIT, MAXDIMS_ORBIT, ERROR_ORBIT)

       CALL H5DCLOSE_F(DSET_ID_ORBIT,ERROR_ORBIT)

       !! allocate OMI time array

       PRINT *,"DIMS_ORBIT: ", DIMS_ORBIT

       ALLOCATE(TIME_ORBIT(DIMS_ORBIT(3)))

       !! read times

       dsetname = '/HDFEOS/SWATHS/ColumnAmountNO2/Geolocation Fields/Time'

       CALL H5DOPEN_F(FILE_ID_ORBIT, DSETNAME, DSET_ID_ORBIT, ERROR_ORBIT)

       CALL H5DREAD_F(DSET_ID_ORBIT, H5T_NATIVE_DOUBLE, TIME_ORBIT, (/DATA_DIMS_ORBIT(3),0/), ERROR_ORBIT)

       CALL H5DCLOSE_F(DSET_ID_ORBIT,ERROR_ORBIT)

       ! check if current hour is in dataset
       CALL TAI2UTC(REAL(TIME_ORBIT(1),8),IY,IM,ID,IH,IMIN,SEC)
       MIN_HOUR = IH

       CALL TAI2UTC(REAL(TIME_ORBIT(DIMS_ORBIT(3)),8),IY,IM,ID,IH,IMIN,SEC)
       MAX_HOUR = IH

       ! go to next dataset if current hour is not contained in dataset
       IF(GC_HOUR<MIN_HOUR .OR. GC_HOUR>MAX_HOUR) THEN
          CALL H5FCLOSE_F(FILE_ID_ORBIT,ERROR_ORBIT)
          CALL H5CLOSE_F(ERROR_ORBIT)
          DEALLOCATE(TIME_ORBIT)
          CYCLE
       ENDIF

       PRINT *,"Found matching OMI file for hour ", DAY, ",",GC_HOUR, ":", TRIM(FILE_ORBIT)

       ALLOCATE(LON_ORBIT(DIMS_ORBIT(2),DIMS_ORBIT(3)))
       ALLOCATE(LAT_ORBIT(DIMS_ORBIT(2),DIMS_ORBIT(3)))
       ALLOCATE(AMF_TROP_ORBIT(DIMS_ORBIT(2),DIMS_ORBIT(3)))
       ALLOCATE(NO2_TROP(DIMS_ORBIT(2),DIMS_ORBIT(3)))
       ALLOCATE(NO2_TROP_STD(DIMS_ORBIT(2),DIMS_ORBIT(3)))
       ALLOCATE(VIEW_ZENITH(DIMS_ORBIT(2),DIMS_ORBIT(3)))
       ALLOCATE(SOLAR_ZENITH(DIMS_ORBIT(2),DIMS_ORBIT(3)))
       ALLOCATE(CLOUDFR(DIMS_ORBIT(2),DIMS_ORBIT(3)))
       ALLOCATE(SCW_P(DIMS_ORBIT(1)))
       ALLOCATE(SCATTERING_WEIGHTS(DIMS_ORBIT(1),DIMS_ORBIT(2),DIMS_ORBIT(3)) )

       !! read in OMI data

       !! read tropospheric air mass factors

       DSETNAME = '/HDFEOS/SWATHS/ColumnAmountNO2/Data Fields/AmfTrop'

       CALL H5DOPEN_F(FILE_ID_ORBIT, DSETNAME, DSET_ID_ORBIT, ERROR_ORBIT)

       CALL H5DREAD_F(DSET_ID_ORBIT, H5T_NATIVE_DOUBLE, AMF_TROP_ORBIT, DATA_DIMS_ORBIT, ERROR_ORBIT)

       CALL H5DCLOSE_F(DSET_ID_ORBIT,ERROR_ORBIT)

       !! read tropospheric NO2 column

       DSETNAME = '/HDFEOS/SWATHS/ColumnAmountNO2/Data Fields/ColumnAmountNO2Trop'

       CALL H5DOPEN_F(FILE_ID_ORBIT, DSETNAME, DSET_ID_ORBIT, ERROR_ORBIT)

       CALL H5DREAD_F(DSET_ID_ORBIT, H5T_NATIVE_DOUBLE, NO2_TROP, DATA_DIMS_ORBIT(2:3), ERROR_ORBIT)

       CALL H5DCLOSE_F(DSET_ID_ORBIT,ERROR_ORBIT)

       !! read tropospheric NO2 column

       DSETNAME = '/HDFEOS/SWATHS/ColumnAmountNO2/Data Fields/ColumnAmountNO2TropStd'

       CALL H5DOPEN_F(FILE_ID_ORBIT, DSETNAME, DSET_ID_ORBIT, ERROR_ORBIT)

       CALL H5DREAD_F(DSET_ID_ORBIT, H5T_NATIVE_DOUBLE, NO2_TROP_STD, DATA_DIMS_ORBIT(2:3), ERROR_ORBIT)

       CALL H5DCLOSE_F(DSET_ID_ORBIT,ERROR_ORBIT)

       !! read longitudes

       DSETNAME = '/HDFEOS/SWATHS/ColumnAmountNO2/Geolocation Fields/Longitude'

       CALL H5DOPEN_F(FILE_ID_ORBIT, DSETNAME, DSET_ID_ORBIT, ERROR_ORBIT)

       CALL H5DREAD_F(DSET_ID_ORBIT, H5T_NATIVE_DOUBLE, LON_ORBIT, DATA_DIMS_ORBIT, ERROR_ORBIT)

       CALL H5DCLOSE_F(DSET_ID_ORBIT,ERROR_ORBIT)

       !! read latitudes

       DSETNAME = '/HDFEOS/SWATHS/ColumnAmountNO2/Geolocation Fields/Latitude'

       CALL H5DOPEN_F(FILE_ID_ORBIT, DSETNAME, DSET_ID_ORBIT, ERROR_ORBIT)

       CALL H5DREAD_F(DSET_ID_ORBIT, H5T_NATIVE_DOUBLE, LAT_ORBIT, DATA_DIMS_ORBIT(2:3), ERROR_ORBIT)

       CALL H5DCLOSE_F(DSET_ID_ORBIT,ERROR_ORBIT)

       !! read viewing zenith angles

       DSETNAME = '/HDFEOS/SWATHS/ColumnAmountNO2/Geolocation Fields/ViewingZenithAngle'

       CALL H5DOPEN_F(FILE_ID_ORBIT, DSETNAME, DSET_ID_ORBIT, ERROR_ORBIT)

       CALL H5DREAD_F(DSET_ID_ORBIT, H5T_NATIVE_DOUBLE, VIEW_ZENITH, DATA_DIMS_ORBIT(2:3), ERROR_ORBIT)

       CALL H5DCLOSE_F(DSET_ID_ORBIT,ERROR_ORBIT)

       !! read solar zenith angles

       DSETNAME = '/HDFEOS/SWATHS/ColumnAmountNO2/Geolocation Fields/SolarZenithAngle'

       CALL H5DOPEN_F(FILE_ID_ORBIT, DSETNAME, DSET_ID_ORBIT, ERROR_ORBIT)

       CALL H5DREAD_F(DSET_ID_ORBIT, H5T_NATIVE_DOUBLE, SOLAR_ZENITH, DATA_DIMS_ORBIT(2:3), ERROR_ORBIT)

       CALL H5DCLOSE_F(DSET_ID_ORBIT,ERROR_ORBIT)

       !! read cloud fraction

       DSETNAME = '/HDFEOS/SWATHS/ColumnAmountNO2/Data Fields/CloudFraction'

       CALL H5DOPEN_F(FILE_ID_ORBIT, DSETNAME, DSET_ID_ORBIT, ERROR_ORBIT)

       CALL H5DREAD_F(DSET_ID_ORBIT, H5T_NATIVE_DOUBLE, CLOUDFR, DATA_DIMS_ORBIT(2:3), ERROR_ORBIT)

       CALL H5DCLOSE_F(DSET_ID_ORBIT,ERROR_ORBIT)

       !! read scattering weight pressures

       DSETNAME = '/HDFEOS/SWATHS/ColumnAmountNO2/Data Fields/ScatteringWtPressure'

       CALL H5DOPEN_F(FILE_ID_ORBIT, DSETNAME, DSET_ID_ORBIT, ERROR_ORBIT)

       CALL H5DREAD_F(DSET_ID_ORBIT, H5T_NATIVE_DOUBLE, SCW_P, (/DATA_DIMS_ORBIT(1),0/), ERROR_ORBIT)

       CALL H5DCLOSE_F(DSET_ID_ORBIT,ERROR_ORBIT)

       !! read scattering weights

       DSETNAME = '/HDFEOS/SWATHS/ColumnAmountNO2/Data Fields/ScatteringWeight'

       CALL H5DOPEN_F(FILE_ID_ORBIT, DSETNAME, DSET_ID_ORBIT, ERROR_ORBIT)

       CALL H5DREAD_F(DSET_ID_ORBIT, H5T_NATIVE_DOUBLE, SCATTERING_WEIGHTS, DATA_DIMS_ORBIT, ERROR_ORBIT)

       CALL H5DCLOSE_F(DSET_ID_ORBIT,ERROR_ORBIT)

       !! close file

       CALL H5FCLOSE_F(FILE_ID_ORBIT,ERROR_ORBIT)

       !! loop over data

       DO I_OMI=1,DIMS_ORBIT(3)

          DO J_OMI=1,DIMS_ORBIT(2)

             IF(TIME_ORBIT(I_OMI)>0) THEN ! very basic quality check, most likely not needed anymore

                ! Convert TAI93 to UTC

                CALL TAI2UTC(TIME_ORBIT(I_OMI),IY,IM,ID,IH,IMIN,SEC)

                ! A number of conditions have to be met for OMI NO2 data to actually be assimilated

                IF ( ( GC_HOUR .EQ. ih                          ) .AND. &
#if defined(NESTED_NA) || defined(NESTED_CH)
                     ( LON_ORBIT(J_OMI,I_OMI) >= GET_XMID(1)    ) .AND. &
                     ( LON_ORBIT(J_OMI,I_OMI) <= GET_XMID(IIPAR)) .AND. &
                     ( LAT_ORBIT(J_OMI,I_OMI) >= GET_YMID(1)    ) .AND. &
                     ( LAT_ORBIT(J_OMI,I_OMI) <= GET_YMID(JJPAR)) .AND. &
#endif
                     ( LAT_ORBIT(J_OMI,I_OMI) < 60d0            ) .AND. &
                     ( NO2_TROP(J_OMI,I_OMI) > 0d0              ) .AND. &
                     ( NO2_TROP_STD(J_OMI,I_OMI) > 0d0          ) .AND. &
                     ( SOLAR_ZENITH(J_OMI,I_OMI) < 75d0         ) .AND. &
                     ( VIEW_ZENITH(J_OMI,I_OMI) < 65d0          ) .AND. &
                     ( AMF_TROP_ORBIT(J_OMI,I_OMI) > 0d0        ) .AND. &
                     ( CLOUDFR(J_OMI,I_OMI) < 0.2   ) ) THEN

                   ! Get model grid coordinate indices that correspond to the observation

                   IIJJ = GET_IJ(REAL(LON_ORBIT(J_OMI,I_OMI),4), REAL(LAT_ORBIT(J_OMI,I_OMI),4))

                   I = IIJJ(1)
                   J = IIJJ(2)

                   ! initialize variables & arrays

                   GC_NO2 = 0d0
                   GC_NO2_COL = 0d0
                   SCW_GC = 0d0
                   DP = 0d0

                   ! Get GEOS-CHEM NO2 values [#/cm3]

                   DO L = 1, LLPAR

                      IF( ITS_IN_THE_TROP(I,J,L) ) THEN

                         JLOOP = JLOP(I,J,L)
                         GC_NO2(L) = CSPEC_AFTER_CHEM(JLOOP,ID2C(IDNO2))

                      ENDIF

                   ENDDO

                   ! Compute tropospheric NO2 vertical column [#/cm2]

                   GC_NO2_COL = SUM(GC_NO2(:) * BXHEIGHT(I,J,:)*100d0)

                   ! interpolate scattering weights to GEOS-Chem grid to compute GEOS-Chem air mass factors
                   ! question: how do differences in surface pressures used in the retrieval and GEOS-Chem affect the computation below?

                   DO L=1,LLPAR
                      DO K_OMI = 2,dims_orbit(1)

                         IF( GET_PCENTER(I,J,L) < SCW_P(K_OMI-1) .AND. GET_PCENTER(I,J,L) > SCW_P(K_OMI) ) THEN

                            ! linearly interpolate scattering weights to GEOS-Chem center pressures

                            SCW_GC(L) = SCATTERING_WEIGHTS(K_OMI,J_OMI,I_OMI) + &
                                 ( SCATTERING_WEIGHTS(K_OMI-1,J_OMI,I_OMI) - SCATTERING_WEIGHTS(K_OMI,J_OMI,I_OMI) ) * &
                                 ( GET_PCENTER(I,J,L) - SCW_P(K_OMI) ) / ( SCW_P(K_OMI-1) - SCW_P(K_OMI) )

                            ! save pressure difference of edge pressures

                            DP(L) = GET_PEDGE(I,J,L) - GET_PEDGE(I,J,L+1)

                            ! apply temperature correction, as in Bucsela2013, eq. (4)

                            SCW_GC(L) = SCW_GC(L) * ( 1 - 0.003 * ( T(I,J,L) - 220 ) )

                            ! convert NO2 concentrations from number density to mixing ratio, as required for the calculation of the air mass factors from scattering weights

                            GC_NO2(L) = GC_NO2(L) *1d6 / ( AIRDEN(L,I,J)   * XNUMOLAIR )

                            EXIT ! exit loop over K_OMI to go to next cycle in loop over L

                         ENDIF

                      ENDDO
                   ENDDO

                   ! Use GEOS-Chem tropospheric air mass factors to convert vertical column to slant column
                   AMF_GC = SUM(GC_NO2 * DP * SCW_GC)/SUM(GC_NO2 * DP)
                   GC_NO2_COL = AMF_GC*GC_NO2_COL

                   ! The computation above is a little awkward, since the slant column can be computed directly from equation (2) in Bucsela2013 without
                   ! computing the airmass factors and NO2 column first.
                   ! I chose to compute the slant column from the computed air mass factors (which already included the computation of the slant column)
                   ! since the air mass factors might be of diagnostic interest and can be computed and saved
                   ! alongside other observation operator diagnostics. Furthermore, this formulation makes the adjoint of the observation operator somewhat simpler to handle.

                   ! compute slant column difference

                   DIFF = GC_NO2_COL - (NO2_TROP(J_OMI,I_OMI) * AMF_TROP_ORBIT(J_OMI, I_OMI))

                   ! compute slant column standard deviation

                   OBS_ERROR = NO2_TROP_STD(J_OMI,I_OMI) * AMF_TROP_ORBIT(J_OMI, I_OMI)

                   ! update adjoint NO2 concentration

                   DO L=1,LLPAR

                      IF(ITS_IN_THE_TROP(I,J,L)) THEN

                         ! question: how do errors in retrieved surface pressure impact the NO2 column values?
                         ! question: how do errors in simulated surface pressures impact the NO2 column values?

                         IF (SUPER_OBS) THEN
                            SOBS_ADJ_FORCE(I,J,L) = SOBS_ADJ_FORCE(I,J,L) + DIFF/(OBS_ERROR**2) * BXHEIGHT(I,J,L) * 100d0 * AMF_GC
                         ELSE
                            JLOOP = JLOP(I,J,L)
                            CSPEC_AFTER_CHEM_ADJ(JLOOP,ID2C(IDNO2)) = CSPEC_AFTER_CHEM_ADJ(JLOOP,ID2C(IDNO2)) &
                                 + DIFF/(OBS_ERROR**2) * BXHEIGHT(I,J,L) * 100d0 * AMF_GC
                         ENDIF

                      ENDIF

                   ENDDO

                   ! update cost function

                   IF(SUPER_OBS) THEN
                      SOBS_COST_CONTRIBUTION(I,J) = SOBS_COST_CONTRIBUTION(I,J) + 0.5 * (DIFF/OBS_ERROR)**2
                      SOBS_COUNT(I,J) = SOBS_COUNT(I,J) + 1d0
                   ELSE
                      COST_FUNC = COST_FUNC + 0.5 * (DIFF/OBS_ERROR)**2
                   ENDIF

                   ! update diagnostic arrays

                   IF( SUPER_OBS) THEN
                      SOBS_GC(I,J) = SOBS_GC(I,J) + GC_NO2_COL
                      SOBS_OMI(I,J) = SOBS_OMI(I,J) + NO2_TROP(J_OMI,I_OMI)
                      SOBS_BIAS(I,J) = SOBS_BIAS(I,J) + DIFF
                      SOBS_CHISQUARED(I,J) = SOBS_CHISQUARED(I,J) + 0.5 * (DIFF/OBS_ERROR)**2
                   ELSE

                      OMI_BIAS_COUNT(I,J) = OMI_BIAS_COUNT(I,J) + 1d0

                      OMI_NO2_MEAN(I,J) = OMI_NO2_MEAN(I,J) + NO2_TROP(J_OMI,I_OMI) * AMF_TROP_ORBIT(J_OMI, I_OMI)

                      OMI_GEOS_NO2_MEAN(I,J) = OMI_GEOS_NO2_MEAN(I,J) + GC_NO2_COL

                      OMI_NO2_ERR_MEAN(I,J) = OMI_NO2_ERR_MEAN(I,J) + OBS_ERROR

                      OMI_DELTA = DIFF - OMI_BIAS(I,J)

                      OMI_BIAS(I,J) = OMI_BIAS(I,J) + OMI_DELTA/OMI_BIAS_COUNT(I,J)

                      OMI_VAR(I,J) = OMI_VAR(I,J) + OMI_DELTA*(DIFF-OMI_BIAS(I,J))

                      OMI_CHISQUARED(I,J) = OMI_CHISQUARED(I,J) + ( DIFF/OBS_ERROR )**2

                   ENDIF

                   ! store current information in flexible arrays
                   
                   IF(LSAT_HDF_L2) THEN
                      CALL PUSH_FLEX_REAL(FLEX_LON, LON_ORBIT(J_OMI,I_OMI))
                      CALL PUSH_FLEX_REAL(FLEX_LAT, LAT_ORBIT(J_OMI,I_OMI))
                      CALL PUSH_FLEX_REAL(FLEX_TIME, REAL(TIME_ORBIT(I_OMI),8))
                      CALL PUSH_FLEX_REAL(FLEX_OMI_NO2, NO2_TROP(J_OMI,I_OMI))
                      CALL PUSH_FLEX_REAL(FLEX_GC_NO2, GC_NO2_COL)
                   ENDIF

                ENDIF ! data selection IF statement

             ENDIF ! OMI_TIME > 0

          ENDDO ! J
       ENDDO ! I

       IF(SUPER_OBS) THEN

          DO J=1,JJPAR
             DO I=1,IIPAR

                IF(SOBS_COUNT(I,J) > 0d0) THEN

                   DO L=1,LLPAR

                      IF(ITS_IN_THE_TROP(I,J,L)) THEN

                         JLOOP = JLOP(I,J,L)
                         CSPEC_AFTER_CHEM_ADJ(JLOOP,ID2C(IDNO2)) = CSPEC_AFTER_CHEM_ADJ(JLOOP,ID2C(IDNO2)) &
                              + SOBS_ADJ_FORCE(I,J,L)/SOBS_COUNT(I,J)

                      ENDIF

                   ENDDO

                   COST_FUNC = COST_FUNC + SOBS_COST_CONTRIBUTION(I,J)/SOBS_COUNT(I,J)

                   OMI_BIAS_COUNT(I,J) = OMI_BIAS_COUNT(I,J) + 1d0

                   OMI_NO2_MEAN(I,J) = OMI_NO2_MEAN(I,J) + SOBS_OMI(I,J)/SOBS_COUNT(I,J)

                   OMI_GEOS_NO2_MEAN(I,J) = OMI_GEOS_NO2_MEAN(I,J) + SOBS_GC(I,J)/SOBS_COUNT(I,J)

                   OMI_NO2_ERR_MEAN(I,J) = OMI_NO2_ERR_MEAN(I,J) + OBS_ERROR  ! mkeller: need to change this to reflect super observation error, but how?

                   ! calculate bias and variance of GC-OMI bias using numerically stable one-pass algorithm (Chan83)

                   OMI_DELTA = SOBS_BIAS(I,J)/SOBS_COUNT(I,J) - OMI_BIAS(I,J)

                   OMI_BIAS(I,J) = OMI_BIAS(I,J) + OMI_DELTA/OMI_BIAS_COUNT(I,J)

                   OMI_VAR(I,J) = OMI_VAR(I,J) + OMI_DELTA*(SOBS_BIAS(I,J)/SOBS_COUNT(I,J) - OMI_BIAS(I,J))

                   OMI_CHISQUARED(I,J) = OMI_CHISQUARED(I,J) + SOBS_CHISQUARED(I,J)/SOBS_COUNT(I,J)

                ENDIF

             ENDDO
          ENDDO

       ENDIF

       ! deallocate OMI arrays

       IF(ALLOCATED(LON_ORBIT)) DEALLOCATE(LON_ORBIT)
       IF(ALLOCATED(LAT_ORBIT)) DEALLOCATE(LAT_ORBIT)
       IF(ALLOCATED(TIME_ORBIT)) DEALLOCATE(TIME_ORBIT)
       IF(ALLOCATED(AMF_TROP_ORBIT)) DEALLOCATE(AMF_TROP_ORBIT)
       IF(ALLOCATED(NO2_TROP)) DEALLOCATE(NO2_TROP)
       IF(ALLOCATED(NO2_TROP_STD)) DEALLOCATE(NO2_TROP_STD)
       IF(ALLOCATED(VIEW_ZENITH)) DEALLOCATE(VIEW_ZENITH)
       IF(ALLOCATED(SOLAR_ZENITH)) DEALLOCATE(SOLAR_ZENITH)
       IF(ALLOCATED(CLOUDFR)) DEALLOCATE(CLOUDFR)
       IF(ALLOCATED(SCW_P)) DEALLOCATE(SCW_P)
       IF(ALLOCATED(SCATTERING_WEIGHTS)) DEALLOCATE(SCATTERING_WEIGHTS)

    ENDDO ! loop over OMI files

    CLOSE(65)

  END SUBROUTINE CALC_OMI_NO2_FORCE

  !-----------------------------------------------------------------------------!

  SUBROUTINE TAI2UTC(tai93,iy,im,id,ih,imin,sec)

    !!
    !! SUBROUTINE TAI2UTC converts TAI93 time (seconds since 1.1.1993) to UTC
    !!
    !! adapted from
    !! code.google.com/p/miyoshi/source/browse/trunk/common/common.f90
    !!

    IMPLICIT NONE

    INTEGER,PARAMETER :: N=7  ! number of leap seconds after Jan. 1, 1993
    INTEGER,PARAMETER :: LEAPSEC(N) = (/ 15638399, 47174400, 94608001, 141868802, 189302403, 410227204, 504921605/)
    REAL*8,INTENT(IN) :: TAI93
    INTEGER,INTENT(OUT) :: IY,IM,ID,IH,IMIN
    REAL*8,INTENT(OUT) :: SEC
    REAL*8,PARAMETER :: MINS = 60.0D0
    REAL*8,PARAMETER :: HOUR = 60.0D0*MINS
    REAL*8,PARAMETER :: DAY = 24.0D0*HOUR
    REAL*8,PARAMETER :: YEAR = 365.0D0*DAY
    INTEGER,PARAMETER :: MDAYS(12) = (/31,28,31,30,31,30,31,31,30,31,30,31/)
    REAL*8 :: WK,TAI
    INTEGER :: DAYS,I,LEAP

    TAI = TAI93
    SEC = 0.0D0

    DO I=1,N

       IF(FLOOR(TAI93) == LEAPSEC(I)+1) THEN
          SEC = 60.0D0 + TAI93-FLOOR(TAI93)
       ENDIF

       IF(FLOOR(TAI93) > LEAPSEC(I)) TAI = TAI -1.0D0

    END DO

    IY = 1993 + FLOOR(TAI /YEAR)
    WK = TAI - REAL(IY-1993)*YEAR - FLOOR(REAL(IY-1993)/4.0)*DAY

    IF(WK < 0.0D0) THEN
       IY = IY -1
       WK = TAI - REAL(IY-1993)*YEAR - FLOOR(REAL(IY-1993)/4.0)*DAY
    END IF

    DAYS = FLOOR(WK/DAY)
    WK = WK - REAL(DAYS)*DAY
    IM = 1

    DO I=1,12

       LEAP = 0
       IF(IM == 2 .AND. MOD(IY,4)==0) LEAP=1
       IF(IM == I .AND. DAYS >= MDAYS(I)+LEAP) THEN
          IM = IM + 1
          DAYS = DAYS - MDAYS(I)-LEAP
       END IF

    END DO

    ID = DAYS +1

    IH = FLOOR(WK/HOUR)
    WK = WK - REAL(IH)*HOUR
    IMIN = FLOOR(WK/MINS)

    IF(SEC < 60.0D0) SEC = WK - REAL(IMIN)*MINS

    RETURN

  END SUBROUTINE TAI2UTC

  !--------------------------------------------------------------------------------

  SUBROUTINE MAKE_OMI_BIAS_FILE_HDF5(FILE_ID)

    !!
    !! Subroutine MAKE_OMI_BIAS_FILE_HDF5
    !!
    !! Arguments:
    !! (1) FILE_ID: ID of the HDF5 diagnostic file
    !!

    USE HDF5

    USE GRID_MOD,            ONLY : GET_XMID, GET_YMID
    USE DIRECTORY_ADJ_MOD,   ONLY : OPTDATA_DIR
    USE ADJ_ARRAYS_MOD,      ONLY : N_CALC
    USE ADJ_ARRAYS_MOD,      ONLY : EXPAND_NAME
    USE LOGICAL_ADJ_MOD,     ONLY : LSAT_HDF_L2, LSAT_HDF_L3

    INTEGER(HID_T) :: FILE_ID

    INTEGER(HID_T) :: SPACE_LON, SPACE_LAT, SPACE_RAW
    INTEGER(HID_T) :: LON_ID, LAT_ID
    INTEGER(HID_T) :: LON_RAW_ID, LAT_RAW_ID, TIME_RAW_ID
    INTEGER(HID_T) :: SPACE_OMI
    INTEGER(HID_T) :: DSET_OMI_NO2_ID
    INTEGER(HID_T) :: DSET_OMI_GC_NO2_ID
    INTEGER(HID_T) :: DSET_OMI_NO2_RAW_ID
    INTEGER(HID_T) :: DSET_OMI_GC_NO2_RAW_ID
    INTEGER(HID_T) :: DSET_OMI_BIAS_ID
    INTEGER(HID_T) :: DSET_OMI_ERR_ID
    INTEGER(HID_T) :: DSET_OMI_CHISQUARED_ID
    INTEGER(HID_T) :: DSET_OMI_COUNT_ID

    INTEGER(HID_T) :: ASPACE_ID, ATYPE_ID, ATT_ID
    INTEGER(HSIZE_T) :: ADIMS(1)

    INTEGER(HID_T) :: OMI_GROUP_ID, LEVEL3_GROUP_ID, GRID_DATA_GROUP_ID, RAW_DATA_GROUP_ID, GRID_GROUP_ID

    INTEGER(HSIZE_T) :: DIMS(2), DIM_LON(1), DIM_LAT(1), DIM_RAW(1)

    INTEGER :: HDF_ERR
    INTEGER :: RANK = 2

    INTEGER :: I,J,L
    REAL*4 :: MISS_VAL = -999.9
    REAL*4 :: LON_VALS(IIPAR), LAT_VALS(JJPAR)

    ! populate lon & lat arrays

    DO I=1,IIPAR
       LON_VALS(I)=GET_XMID(I)
    ENDDO

    DO J=1,JJPAR
       LAT_VALS(J)=GET_YMID(J)
    ENDDO

    ! compute diagnostics

    DO J=1,JJPAR
       DO I=1,IIPAR

          IF(OMI_BIAS_COUNT(I,J)>0d0) THEN

             OMI_NO2_MEAN(I,J) = REAL(OMI_NO2_MEAN(I,J)/OMI_BIAS_COUNT(I,J))*1E-15
             OMI_GEOS_NO2_MEAN(I,J) = REAL(OMI_GEOS_NO2_MEAN(I,J)/OMI_BIAS_COUNT(I,J))*1E-15
             OMI_NO2_ERR_MEAN(I,J) = REAL(OMI_NO2_ERR_MEAN(I,J)/OMI_BIAS_COUNT(I,J))*1E-15
             OMI_BIAS(I,J) = OMI_BIAS(I,J)*1E-15
             OMI_VAR(I,J) = OMI_VAR(I,J)*1E-15

          ELSE

             OMI_NO2_MEAN(I,J) = MISS_VAL
             OMI_GEOS_NO2_MEAN(I,J) = MISS_VAL
             OMI_NO2_ERR_MEAN(I,J) = MISS_VAL
             OMI_BIAS(I,J) = MISS_VAL
             OMI_VAR(I,J) = MISS_VAL
             OMI_CHISQUARED(I,J) = MISS_VAL

          ENDIF

       ENDDO
    ENDDO

    ! define dimensions of diagnostic variables

    DIMS(1) = IIPAR
    DIMS(2) = JJPAR

    ADIMS(1) = 1

    DIM_LON = IIPAR
    DIM_LAT = JJPAR

    IF(LSAT_HDF_L2) DIM_RAW = FLEX_LON%CURRENT_N

    ! open HDF5 interface

    CALL H5OPEN_F(HDF_ERR)

    ! create group structure in file

    CALL H5GCREATE_F(FILE_ID,"OMI",OMI_GROUP_ID,HDF_ERR)
    
    IF(LSAT_HDF_L3) THEN
       CALL H5GCREATE_F(OMI_GROUP_ID,"Level3",LEVEL3_GROUP_ID,HDF_ERR)
       CALL H5GCREATE_F(LEVEL3_GROUP_ID,"Grid",GRID_GROUP_ID,HDF_ERR)
       CALL H5GCREATE_F(LEVEL3_GROUP_ID,"Data",GRID_DATA_GROUP_ID,HDF_ERR)
    ENDIF

    IF(LSAT_HDF_L2) THEN
       CALL H5GCREATE_F(OMI_GROUP_ID,"Level2",RAW_DATA_GROUP_ID,HDF_ERR)
    ENDIF

    IF(LSAT_HDF_L3) THEN

       ! write grid information for Level3 data
       
       CALL H5SCREATE_SIMPLE_F(1,DIM_LON,SPACE_LON,HDF_ERR)
       CALL H5SCREATE_SIMPLE_F(1,DIM_LAT,SPACE_LAT,HDF_ERR)
       
       CALL H5DCREATE_F(GRID_GROUP_ID, "/OMI/Level3/Grid/Longitude", H5T_IEEE_F32LE, SPACE_LON, LON_ID, HDF_ERR)
       CALL H5DCREATE_F(GRID_GROUP_ID, "/OMI/Level3/Grid/Latitude", H5T_IEEE_F32LE, SPACE_LAT, LAT_ID, HDF_ERR)
       
       CALL H5DWRITE_F(LON_ID, H5T_NATIVE_REAL, LON_VALS, DIM_LON, HDF_ERR)
       CALL H5DWRITE_F(LAT_ID, H5T_NATIVE_REAL, LAT_VALS, DIM_LAT, HDF_ERR)
       
       CALL WRITE_ATTRIBUTES(LON_ID,"Longitude","degrees")
       CALL WRITE_ATTRIBUTES(LAT_ID,"Latitude","degrees")
       
       CALL H5DCLOSE_F(LON_ID, HDF_ERR)
       CALL H5DCLOSE_F(LAT_ID, HDF_ERR)

       CALL H5SCLOSE_F(SPACE_LON, HDF_ERR)
       CALL H5SCLOSE_F(SPACE_LAT, HDF_ERR)

       ! create dataspace for Level3 OMI diagnostics
       
       CALL H5SCREATE_SIMPLE_F(RANK,DIMS,SPACE_OMI,HDF_ERR)
       
       ! note: all datasets are created  as little-endian 32 bit IEEE float
       ! write OMI NO2 concentrations
       
       CALL H5DCREATE_F(GRID_DATA_GROUP_ID,"/OMI/Level3/Data/OMI_NO2",H5T_IEEE_F32LE, SPACE_OMI, DSET_OMI_NO2_ID, HDF_ERR)
       
       CALL H5DWRITE_F(DSET_OMI_NO2_ID, H5T_NATIVE_REAL, OMI_NO2_MEAN, DIMS, HDF_ERR)
       
       CALL WRITE_ATTRIBUTES(DSET_OMI_NO2_ID,"Mean OMI NO2 Tropospheric Slant Columns", "#/cm^2 * 1e15")
       
       CALL H5DCLOSE_F(DSET_OMI_NO2_ID,HDF_ERR)
       
       ! write OMI_GC NO2 concentrations
       
       CALL H5DCREATE_F(GRID_DATA_GROUP_ID,"/OMI/Level3/Data/OMI_GC_NO2",H5T_IEEE_F32LE, SPACE_OMI, DSET_OMI_GC_NO2_ID, HDF_ERR)
       
       CALL H5DWRITE_F(DSET_OMI_GC_NO2_ID, H5T_NATIVE_REAL, OMI_GEOS_NO2_MEAN, ADIMS, HDF_ERR)
       
       CALL WRITE_ATTRIBUTES(DSET_OMI_GC_NO2_ID, "Mean OMI GEOS_Chem NO2 Tropospheric Slant Columns", "#/cm^2 * 1e15")
       
       CALL H5DCLOSE_F(DSET_OMI_GC_NO2_ID,HDF_ERR)
       
       ! write OMI_GC NO2 bias
       
       CALL H5DCREATE_F(GRID_DATA_GROUP_ID,"/OMI/Level3/Data/OMI_BIAS",H5T_IEEE_F32LE, SPACE_OMI, DSET_OMI_BIAS_ID, HDF_ERR)
       
       CALL H5DWRITE_F(DSET_OMI_BIAS_ID, H5T_NATIVE_REAL, OMI_BIAS, DIMS, HDF_ERR)
       
       CALL WRITE_ATTRIBUTES(DSET_OMI_BIAS_ID, "Mean OMI NO2 bias Tropospheric Slant Columns", "#/cm^2 * 1e15")
       
       CALL H5DCLOSE_F(DSET_OMI_BIAS_ID,HDF_ERR)
       
       ! write OMI_GC NO2 error
       
       CALL H5DCREATE_F(GRID_DATA_GROUP_ID,"/OMI/Level3/Data/OMI_ERROR",H5T_IEEE_F32LE, SPACE_OMI, DSET_OMI_ERR_ID, HDF_ERR)
       
       CALL H5DWRITE_F(DSET_OMI_ERR_ID, H5T_NATIVE_REAL, OMI_NO2_ERR_MEAN, DIMS, HDF_ERR)
       
       CALL WRITE_ATTRIBUTES(DSET_OMI_ERR_ID, "Mean OMI NO2 Observational Error", "#/cm^2 * 1e15")
       
       CALL H5DCLOSE_F(DSET_OMI_ERR_ID,HDF_ERR)
       
       ! write Chi-Squared
       
       CALL H5DCREATE_F(GRID_DATA_GROUP_ID, "/OMI/Level3/Data/Chi-Squared",H5T_IEEE_F32LE, SPACE_OMI, DSET_OMI_CHISQUARED_ID, HDF_ERR)
       
       CALL H5DWRITE_F(DSET_OMI_CHISQUARED_ID, H5T_NATIVE_REAL, OMI_CHISQUARED, DIMS, HDF_ERR)
       
       CALL WRITE_ATTRIBUTES(DSET_OMI_CHISQUARED_ID, "OMI NO2 Chi-Squared Value", "1")
       
       CALL H5DCLOSE_F(DSET_OMI_CHISQUARED_ID,HDF_ERR)
       
       ! write OMI_GC NO2 count
       
       CALL H5DCREATE_F(GRID_DATA_GROUP_ID, "/OMI/Level3/Data/OMI_COUNT",H5T_IEEE_F32LE, SPACE_OMI, DSET_OMI_COUNT_ID, HDF_ERR)
       
       CALL H5DWRITE_F(DSET_OMI_COUNT_ID, H5T_NATIVE_REAL, OMI_BIAS_COUNT, DIMS, HDF_ERR)
       
       CALL WRITE_ATTRIBUTES(DSET_OMI_COUNT_ID, "OMI data count", "1")
       
       CALL H5DCLOSE_F(DSET_OMI_COUNT_ID,HDF_ERR)
       
       !close data space & groups
       
       CALL H5SCLOSE_F(SPACE_OMI,HDF_ERR)
       CALL H5GCLOSE_F(GRID_DATA_GROUP_ID, HDF_ERR)
       CALL H5GCLOSE_F(GRID_GROUP_ID, HDF_ERR)

    ENDIF ! LSAT_HDF_L3

    IF(LSAT_HDF_L2) THEN

       ! create dataspace for Level2 OMI diagnostics
       
       CALL H5SCREATE_SIMPLE_F(1,DIM_RAW,SPACE_RAW,HDF_ERR)
       
       ! write raw longitudes
       
       CALL H5DCREATE_F(RAW_DATA_GROUP_ID,"/OMI/Level2/Longitude",H5T_IEEE_F32LE, SPACE_RAW, LON_RAW_ID, HDF_ERR)
       
       CALL H5DWRITE_F(LON_RAW_ID, H5T_NATIVE_REAL, REAL(FLEX_LON%DATA(1:FLEX_LON%CURRENT_N),4), DIM_RAW, HDF_ERR)
       
       CALL WRITE_ATTRIBUTES(LON_RAW_ID, "Longitude", "degrees")
       
       CALL H5DCLOSE_F(LON_RAW_ID,HDF_ERR)
       
       ! write raw latitudes
       
       CALL H5DCREATE_F(RAW_DATA_GROUP_ID, "/OMI/Level2/Latitude",H5T_IEEE_F32LE, SPACE_RAW, LAT_RAW_ID, HDF_ERR)
       
       CALL H5DWRITE_F(LAT_RAW_ID, H5T_NATIVE_REAL, REAL(FLEX_LAT%DATA(1:FLEX_LAT%CURRENT_N),4), DIM_RAW, HDF_ERR)
       
       CALL WRITE_ATTRIBUTES(LAT_RAW_ID, "Latitude", "degrees")
       
       CALL H5DCLOSE_F(LAT_RAW_ID,HDF_ERR)
       
       ! write raw times
       
       CALL H5DCREATE_F(RAW_DATA_GROUP_ID,"/OMI/Level2/Time",H5T_IEEE_F32LE, SPACE_RAW, TIME_RAW_ID, HDF_ERR)
       
       CALL H5DWRITE_F(TIME_RAW_ID, H5T_NATIVE_REAL, REAL(FLEX_TIME%DATA(1:FLEX_TIME%CURRENT_N),4), DIM_RAW, HDF_ERR)
       
       CALL WRITE_ATTRIBUTES(TIME_RAW_ID,"Time","seconds since 1.1.1993")
       
       CALL H5DCLOSE_F(TIME_RAW_ID,HDF_ERR)
       
       ! write raw OMI NO2 concentrations
       
       CALL H5DCREATE_F(RAW_DATA_GROUP_ID, "/OMI/Level2/OMI_NO2",H5T_IEEE_F32LE, SPACE_RAW, DSET_OMI_NO2_RAW_ID, HDF_ERR)
       
       CALL H5DWRITE_F(DSET_OMI_NO2_RAW_ID, H5T_NATIVE_REAL, REAL(FLEX_OMI_NO2%DATA(1:FLEX_OMI_NO2%CURRENT_N)*1e-15,4), DIM_RAW, HDF_ERR)
       
       CALL WRITE_ATTRIBUTES(DSET_OMI_NO2_RAW_ID,"OMI NO2 tropospheric columns","#/cm^2 * 1e15")
       
       CALL H5DCLOSE_F(DSET_OMI_NO2_RAW_ID,HDF_ERR)
       
       ! write raw OMI GC-NO2 concentrations
       
       CALL H5DCREATE_F(RAW_DATA_GROUP_ID, "/OMI/Level2/OMI_GC_NO2", H5T_IEEE_F32LE, SPACE_RAW, DSET_OMI_GC_NO2_RAW_ID, HDF_ERR)
       
       CALL H5DWRITE_F(DSET_OMI_GC_NO2_RAW_ID, H5T_NATIVE_REAL, REAL(FLEX_GC_NO2%DATA(1:FLEX_GC_NO2%CURRENT_N)*1e-15,4), DIM_RAW, HDF_ERR)
       
       CALL WRITE_ATTRIBUTES(DSET_OMI_GC_NO2_RAW_ID,"GEOS-Chem OMI NO2 tropospheric columns","#/cm^2 * 1e15")
       
       CALL H5DCLOSE_F(DSET_OMI_GC_NO2_RAW_ID,HDF_ERR)
       
       ! close data space & group
       
       CALL H5SCLOSE_F(SPACE_RAW,HDF_ERR)
       CALL H5GCLOSE_F(RAW_DATA_GROUP_ID, HDF_ERR)

       ! clear flexible arrays
       CALL CLEAR_FLEX_REAL(FLEX_LON)
       CALL CLEAR_FLEX_REAL(FLEX_LAT)
       CALL CLEAR_FLEX_REAL(FLEX_TIME)
       CALL CLEAR_FLEX_REAL(FLEX_OMI_NO2)
       CALL CLEAR_FLEX_REAL(FLEX_GC_NO2)
       
    ENDIF ! LSAT_HDF_L2

    CALL H5GCLOSE_F(OMI_GROUP_ID, HDF_ERR)

    ! close HDF5 interface

    CALL H5CLOSE_F(HDF_ERR)

  END SUBROUTINE MAKE_OMI_BIAS_FILE_HDF5

  !-------------------------------------------------------------------------------

  SUBROUTINE WRITE_ATTRIBUTES(DSET_ID,LONGNAME,UNIT)

    USE HDF5

    INTEGER(HID_T) :: DSET_ID
    CHARACTER(LEN=*) :: LONGNAME
    CHARACTER(LEN=*) :: UNIT

    INTEGER(HID_T) :: ASPACE_ID, ATYPE_ID, ATT_ID
    INTEGER(HSIZE_T) :: ADIMS(1)

    INTEGER :: HDF_ERR

    ADIMS(1) = 1

    ! create attribute "Long name"

    CALL H5SCREATE_SIMPLE_F(1,ADIMS,ASPACE_ID,HDF_ERR)

    CALL H5TCOPY_F(H5T_NATIVE_CHARACTER,ATYPE_ID,HDF_ERR)
    CALL H5TSET_SIZE_F(ATYPE_ID,LEN(LONGNAME),HDF_ERR)

    CALL H5ACREATE_F(DSET_ID,"Long name", ATYPE_ID,ASPACE_ID,ATT_ID,HDF_ERR)
    CALL H5AWRITE_F(ATT_ID,ATYPE_ID,LONGNAME, ADIMS,HDF_ERR)

    CALL H5ACLOSE_F(ATT_ID,HDF_ERR)
    CALL H5SCLOSE_F(ASPACE_ID,HDF_ERR)

    ! create attribute "Unit"

    CALL H5SCREATE_SIMPLE_F(1,ADIMS,ASPACE_ID,HDF_ERR)

    CALL H5TCOPY_F(H5T_NATIVE_CHARACTER,ATYPE_ID,HDF_ERR)
    CALL H5TSET_SIZE_F(ATYPE_ID,LEN(UNIT),HDF_ERR)

    CALL H5ACREATE_F(DSET_ID,"Unit", ATYPE_ID,ASPACE_ID,ATT_ID,HDF_ERR)
    CALL H5AWRITE_F(ATT_ID,ATYPE_ID,UNIT, ADIMS,HDF_ERR)

    CALL H5ACLOSE_F(ATT_ID,HDF_ERR)
    CALL H5SCLOSE_F(ASPACE_ID,HDF_ERR)

  END SUBROUTINE WRITE_ATTRIBUTES

  !--------------------------------------------------------------------------------

  !mkeller: helper routines for managing flexible arrays
  !         reinventing the wheel here, but hey...

  SUBROUTINE INIT_FLEX_REAL(INPUT)

    TYPE(FLEX_REAL):: INPUT
    INPUT%CURRENT_N = 0
    INPUT%MAX_N = 1000
    IF(ALLOCATED(INPUT%DATA)) DEALLOCATE(INPUT%DATA) ! safety first
    ALLOCATE(INPUT%DATA(INPUT%MAX_N))

  END SUBROUTINE INIT_FLEX_REAL

  SUBROUTINE GROW_FLEX_REAL(INPUT)

    TYPE(FLEX_REAL) :: INPUT
    REAL*8, ALLOCATABLE :: TEMP_ARRAY(:)
    ALLOCATE(TEMP_ARRAY(INPUT%MAX_N * 2))
    TEMP_ARRAY(1:INPUT%MAX_N) = INPUT%DATA
    DEALLOCATE(INPUT%DATA)
    ALLOCATE(INPUT%DATA(INPUT%MAX_N * 2))
    INPUT%DATA = TEMP_ARRAY
    DEALLOCATE(TEMP_ARRAY)
    INPUT%MAX_N = INPUT%MAX_N * 2

  END SUBROUTINE GROW_FLEX_REAL

  SUBROUTINE PUSH_FLEX_REAL(INPUT, NEW_VAL)

    TYPE(FLEX_REAL) :: INPUT
    REAL*8 :: NEW_VAL
    IF(INPUT%CURRENT_N == INPUT%MAX_N) THEN
       CALL GROW_FLEX_REAL(INPUT)
    ENDIF
    INPUT%CURRENT_N = INPUT%CURRENT_N + 1
    INPUT%DATA(INPUT%CURRENT_N) = NEW_VAL

  END SUBROUTINE PUSH_FLEX_REAL

  SUBROUTINE CLEAR_FLEX_REAL(INPUT)

    TYPE(FLEX_REAL) :: INPUT
    IF(ALLOCATED(INPUT%DATA)) DEALLOCATE(INPUT%DATA)

  END SUBROUTINE CLEAR_FLEX_REAL

END MODULE OMI_NO2_OBS_MOD
