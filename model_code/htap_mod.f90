!------------------------------------------------------------------------------
!              CU Boulder Adjoint Modeling Group (Henze Group)                !
!------------------------------------------------------------------------------
!
! !MODULE: htap_mod
!
! !DESCRIPTION: Module HTAP\_MOD contains variables and routines to
!  read the HTAP V2 anthropogenic emissions.
!\\
!\\
! !INTERFACE:
!
MODULE HTAP_MOD
  !
  ! !USES:
  !
  IMPLICIT NONE
#     include "define.h"
#     include "netcdf.inc"
  PRIVATE
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  !
  PUBLIC  :: CLEANUP_HTAP
  PUBLIC  :: EMISS_HTAP
  PUBLIC  :: GET_HTAP
  PUBLIC  :: SC_AIR, SC_SHIPS, SC_ENERGY, SC_INDUSTRY
  PUBLIC  :: SC_TRANSPORT, SC_RESIDENTIAL, SC_AGRICULTURE
  PUBLIC  :: SC_BC, SC_CO, SC_OC, SC_NH3, SC_NOX, SC_SO2
  PUBLIC  :: SC_CH2O, SC_PM25, htap_rcptr_mask
  PUBLIC  :: LOCN, LNAM, LEUR, LSAS, LEAS
  PUBLIC  :: LSEA, LPAN, LNAF, LSAF, LMDE, LMCA
  PUBLIC  :: LSAM, LRBU, LCAS, LNPO, LSPO
  PUBLIC  :: LOCN20, LOCN21, LOCN22, LOCN23, LOCN24
  PUBLIC  :: LOCN25, LOCN26, LOCN27, LOCN28
  PUBLIC  :: LNAM31, LNAM32, LNAM33, LNAM34, LNAM35
  PUBLIC  :: LNAM36, LEUR41, LEUR42, LEUR43, LEUR44
  PUBLIC  :: LSAS51, LSAS52, LSAS53, LEAS61, LEAS62
  PUBLIC  :: LEAS63, LEAS64, LEAS65, LEAS66, LSEA71
  PUBLIC  :: LSEA72, LPAN81, LPAN82, LPAN83, LNAF91
  PUBLIC  :: LNAF92, LNAF93
  PUBLIC  :: LSAF101, LSAF102, LSAF103, LMDE111, LMDE112
  PUBLIC  :: LMDE113, LMCA121, LMCA122, LMCA123, LMCA124
  PUBLIC  :: LSAM131, LSAM132, LSAM133, LSAM134, LRBU141
  PUBLIC  :: LRBU142, LRBU143, LCAS151, LNPO150, LSPO160
  PUBLIC  :: LSPO161
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  !
  PRIVATE :: INIT_HTAP
  PRIVATE :: TOTAL_HTAP_TG
  !
  ! !REMARKS:

  !   (1)
  !
  ! !REVISION HISTORY:
  !  Sep 2013 - Yanko Davila - initial version
  !------------------------------------------------------------------------------
  !
  ! !PRIVATE TYPES:
  !
  ! Arrays for emissions (lat/lon)
  REAL*8,  ALLOCATABLE, TARGET :: BC(:,:)
  REAL*8,  ALLOCATABLE, TARGET :: CO(:,:)
  REAL*8,  ALLOCATABLE, TARGET :: NH3(:,:)
  REAL*8,  ALLOCATABLE, TARGET :: NOX(:,:)
  REAL*8,  ALLOCATABLE, TARGET :: OC(:,:)
  REAL*8,  ALLOCATABLE, TARGET :: SO2(:,:)
  REAL*8,  ALLOCATABLE, TARGET :: PM25(:,:)
  REAL*8,  ALLOCATABLE, TARGET :: ALK4_HTAP(:,:)
  REAL*8,  ALLOCATABLE, TARGET :: ACET_HTAP(:,:)
  REAL*8,  ALLOCATABLE, TARGET :: MEK_HTAP (:,:)
  REAL*8,  ALLOCATABLE, TARGET :: ALD2_HTAP(:,:)
  REAL*8,  ALLOCATABLE, TARGET :: PRPE_HTAP(:,:)
  REAL*8,  ALLOCATABLE, TARGET :: C3H8_HTAP(:,:)
  REAL*8,  ALLOCATABLE, TARGET :: CH2O_HTAP(:,:)
  REAL*8,  ALLOCATABLE, TARGET :: C2H6_HTAP(:,:)

  ! Scaling Factors
  REAL*8                      :: SC_AIR, SC_SHIPS, SC_ENERGY, SC_INDUSTRY
  REAL*8                      :: SC_TRANSPORT, SC_RESIDENTIAL, SC_AGRICULTURE
  REAL*8                      :: SC_BC, SC_CO, SC_OC, SC_NH3, SC_NOX, SC_SO2
  REAL*8                      :: SC_CH2O, SC_PM25
  LOGICAL                     :: LOCN, LNAM, LEUR, LSAS, LEAS
  LOGICAL                     :: LSEA, LPAN, LNAF, LSAF, LMDE, LMCA
  LOGICAL                     :: LSAM, LRBU, LCAS, LNPO, LSPO
  LOGICAL                     :: LOCN20, LOCN21, LOCN22, LOCN23, LOCN24
  LOGICAL                     :: LOCN25, LOCN26, LOCN27, LOCN28
  LOGICAL                     :: LNAM31, LNAM32, LNAM33, LNAM34, LNAM35
  LOGICAL                     :: LNAM36, LEUR41, LEUR42, LEUR43, LEUR44
  LOGICAL                     :: LSAS51, LSAS52, LSAS53, LEAS61, LEAS62
  LOGICAL                     :: LEAS63, LEAS64, LEAS65, LEAS66, LSEA71
  LOGICAL                     :: LSEA72, LPAN81, LPAN82, LPAN83, LNAF91
  LOGICAL                     :: LNAF92, LNAF93
  LOGICAL                     :: LSAF101, LSAF102, LSAF103, LMDE111, LMDE112
  LOGICAL                     :: LMDE113, LMCA121, LMCA122, LMCA123, LMCA124
  LOGICAL                     :: LSAM131, LSAM132, LSAM133, LSAM134, LRBU141
  LOGICAL                     :: LRBU142, LRBU143, LCAS151, LNPO150, LSPO160
  LOGICAL                     :: LSPO161

  CHARACTER(LEN=255)          :: htap_rcptr_mask, htap_src_mask

  !
  ! !DEFINED PARAMETERS:
  !
  REAL*8,  PARAMETER   :: SEC_IN_YEAR  = 86400d0 * 365.25d0

CONTAINS
  !
  ! !IROUTINE: get_htap
  !
  ! !DESCRIPTION: Function GET\_HTAP returns the HTAP V2
  !  emission for GEOS-Chem grid box (I,J) and tracer N.
  !  Emissions ARE returned in units of [kg/cm2/s].
  !\\
  ! !INTERFACE:
  !
  FUNCTION GET_HTAP( I, J, N ) &
           RESULT( VALUE )
    !
    ! !USES:
    !
    USE TRACERID_MOD, ONLY : IDECO,   IDENOX
    USE TRACERID_MOD, ONLY : IDTSO2,  IDTNH3
    USE TRACERID_MOD, ONLY : IDTOCPO, IDTBCPO
    USE TRACERID_MOD, ONLY : IDEALK4, IDEACET
    USE TRACERID_MOD, ONLY : IDEALD2, IDEPRPE
    USE TRACERID_MOD, ONLY : IDEC3H8, IDECH2O
    USE TRACERID_MOD, ONLY : IDEC2H6, IDEMEK
    USE ERROR_MOD,    ONLY : ERROR_STOP


    !
    ! !INPUT PARAMETERS:
    !
    ! Longitude, latitude, hour, and tracer indices
    INTEGER, INTENT(IN)           :: I, J, N

    !
    ! !RETURN VALUE:
    !
    ! Emissions output
    REAL*8                        :: VALUE
    !
    ! !REVISION HISTORY:
    !  Sep 2013 - Yanko Davila - initial version
    !
    !------------------------------------------------------------------------------

    !=================================================================
    ! GET_HTAP begins here!
    !=================================================================

    IF  ( N==IDTBCPO ) THEN

          ! BC [kg/m2/s]
          VALUE = BC(I,J)

    ELSE IF ( N==IDECO ) THEN

          ! CO[kg/m2/s]
          VALUE = CO(I,J)

    ELSE IF ( N==IDTNH3 ) THEN

          ! NH3[kg/m2/s]
          VALUE = NH3(I,J)

    ELSE IF ( N==IDENOX ) THEN

          ! NOX[kg/m2/s]
          VALUE = NOX(I,J)

    ELSE IF ( N==IDTOCPO ) THEN

          ! OC [kg/m2/s]
          VALUE = OC(I,J)

    ELSE IF ( N==IDTSO2 ) THEN

          ! SO2 [kg/m2/s]
          VALUE = SO2(I,J)

    ELSE IF ( N==IDEALK4 ) THEN

          ! NMVOC[kg/m2/s]
          VALUE = ALK4_HTAP(I,J)

    ELSE IF ( N==IDEACET ) THEN

          ! NMVOC[kg/m2/s]
          VALUE = ACET_HTAP(I,J)

    ELSE IF ( N==IDEMEK ) THEN

          ! NMVOC[kg/m2/s]
          VALUE = MEK_HTAP(I,J)

    ELSE IF ( N==IDEALD2 ) THEN

          ! NMVOC[kg/m2/s]
          VALUE = ALD2_HTAP(I,J)

    ELSE IF ( N==IDEPRPE ) THEN

          ! NMVOC[kg/m2/s]
          VALUE = PRPE_HTAP(I,J)

    ELSE IF ( N==IDEC3H8 ) THEN

          ! NMVOC[kg/m2/s]
          VALUE = C3H8_HTAP(I,J)

    ELSE IF ( N==IDECH2O ) THEN

          ! NMVOC[kg/m2/s]
          VALUE = CH2O_HTAP(I,J)

    ELSE IF ( N==IDEC2H6 ) THEN

          ! NMVOC[kg/m2/s]
          VALUE = C2H6_HTAP(I,J)

    ELSE

          ! Otherwise stop simulation to indicate
          ! that there are no HTAP emissions for tracer N
          !VALUE = -1d0
          !RETURN
          CALL ERROR_STOP('HTAP emission missing for tracer' , 'GET_HTAP, htap_mod.f90'  )

    ENDIF


    ! Return to calling program
  END FUNCTION GET_HTAP
  !
  ! !IROUTINE: emiss_htap
  !
  ! !DESCRIPTION: Subroutine EMISS\_HTAP reads the HTAP v2
  !  emission fields at 0.1x0.1 resolution and regrids them to the
  !  current model resolution.
  !\\
  !\\
  ! !INTERFACE:
  !
  SUBROUTINE EMISS_HTAP
    !
    ! !USES:
    !
    USE DIRECTORY_MOD,     ONLY : DATA_DIR_1x1
    USE REGRID_A2A_MOD,    ONLY : DO_REGRID_A2A, DO_REGRID_DKH
    USE TIME_MOD,          ONLY : GET_YEAR, GET_MONTH, GET_DAY
    USE TIME_MOD,          ONLY : GET_DAY_OF_WEEK, GET_HOUR
    USE TRACER_MOD,        ONLY : ITS_A_FULLCHEM_SIM
    USE TRACERID_MOD,      ONLY : IDTCO
    USE TRACERID_MOD,      ONLY : IDTNOX, IDTOX
    USE TRACERID_MOD,      ONLY : IDTSO2, IDTNH3
    USE TRACERID_MOD,      ONLY : IDTOCPO, IDTBCPO
    USE TRACERID_MOD,      ONLY : IDEALK4, IDEACET
    USE TRACERID_MOD,      ONLY : IDEALD2, IDEPRPE
    USE TRACERID_MOD,      ONLY : IDEC3H8, IDECH2O
    USE TRACERID_MOD,      ONLY : IDEC2H6, IDEMEK
    USE TRACERID_MOD,      ONLY : NEMANTHRO
    USE m_netcdf_io_open
    USE m_netcdf_io_read
    USE m_netcdf_io_readattr
    USE m_netcdf_io_close
    USE m_netcdf_io_get_dimlen
    ! debug
    USE GRID_MOD,         ONLY : GET_AREA_M2

#     include "CMN_SIZE"          ! Size parameters
    !
    ! !REVISION HISTORY:
    !  Sep 2013 - Yanko Davila - initial version
    ! --------------------------------------------------------
    !
    ! !LOCAL VARIABLES:
    !
    LOGICAL, SAVE              :: FIRST = .TRUE.
    LOGICAL                    :: netCDF
    INTEGER, PARAMETER         :: I01x01 = 3600, J01x01 = 1800
    INTEGER                    :: I, J, IH, THISMONTH, THISYEAR
    INTEGER                    :: SNo, DAY_NUM, DOYT
    INTEGER                    :: HH, KLM, ID,  MN, N
    INTEGER                    :: OFFLINE_ID(15)
    INTEGER                    :: fId1, fId2
    INTEGER                    :: CAT
    INTEGER                    :: II, JJ, ALK
    REAL*8                     :: ARRAY(I01x01,J01x01),    MASK(I01x01,J01x01)
    REAL*8                     :: TMP_MASK(I01x01,J01x01), VOC_DATA(I01x01,J01x01,1)
    REAL*8                     :: ALK4_TMP(I01x01,J01x01)
    REAL*8                     :: C2H6_TMP (I01x01,J01x01)
    REAL*8                     :: C3H8_TMP (I01x01,J01x01)
    REAL*8                     :: PRPE_TMP (I01x01,J01x01)
    REAL*8                     :: CH2O_TMP (I01x01,J01x01)
    REAL*8                     :: ALD2_TMP (I01x01,J01x01)
    REAL*8                     :: ACET_TMP (I01x01,J01x01)
    REAL*8                     :: MEK_TMP  (I01x01,J01x01)
    REAL*8, TARGET             :: GEOS_01x01(I01x01,J01x01)
    CHARACTER(LEN=244)         :: DATA_DIR_HTAP
    CHARACTER(LEN=244)         :: FILENAME, VOC_FILE, ALK_MASK
    CHARACTER(LEN=4)           :: SYEAR
    CHARACTER(LEN=5)           :: SNAME, HTAP_YEAR_2,SId
    CHARACTER(LEN=6)           :: HTAP_YEAR
    CHARACTER(LEN=9)           :: VId
    CHARACTER(LEN=255)         :: LLFILENAME
    CHARACTER(LEN=24)          :: SPCLIST(8), EMIS_TYPE(7), VAR_NAME(8)
    CHARACTER(LEN=2)           :: MONTH
    CHARACTER(LEN=50)          :: ALK_NUM(3)
    REAL*8, POINTER            :: OUTGRID(:,:) => NULL()
    REAL*8, POINTER            :: INGRID(:,:) => NULL()
    ! Days per month
    REAL*8                     :: DAYS_IN_MONTH
    REAL*8                     :: DMON(12)

    !=================================================================
    ! HTAP_ANTHRO begins here!
    !=================================================================

    ! First-time initialization
    IF ( FIRST ) THEN
       CALL INIT_HTAP
       FIRST = .FALSE.
    ENDIF

    ! Get emissions year
    THISYEAR = GET_YEAR()

    ! Get month
    THISMONTH = GET_MONTH()

#if   defined( GEOS_5 ) || defined( MERRA ) || defined( GEOS_57 )
    SNAME = 'GEOS5'
#elif defined( GCAP  )
    SNAME = 'GEOS4'
#elif defined( GEOS_4 )
    SNAME = 'GEOS4'
#endif

    SPCLIST    = (/ 'BC','CO', 'NH3', 'NMVOC','NOx', 'OC', 'SO2', 'PM2.5' /)

    VAR_NAME   = (/ 'emi_bc','emi_co', 'emi_nh3', 'emi_nmvoc','emi_nox', 'emi_oc', 'emi_so2', 'emi_pm2.5' /)

    ALK_NUM    = (/ '_butanes','_pentanes', '_hexanes_and_higher_alkanes' /)

    IF ( THISYEAR  .eq. 2008 )  DMON = (/ 31d0, 29d0, 31d0, 30d0,  &
                                          31d0, 30d0, 31d0, 31d0,  &
                                          30d0, 31d0, 30d0, 31d0 /)

    IF ( THISYEAR  .eq. 2010 )  DMON = (/ 31d0, 28d0, 31d0, 30d0,  &
                                          31d0, 30d0, 31d0, 31d0,  &
                                          30d0, 31d0, 30d0, 31d0 /)

    ! Get days per month
    DAYS_IN_MONTH = DMON( GET_MONTH() )

    ! File with lat/lon edges for regridding
    LLFILENAME = TRIM( DATA_DIR_1x1) // &
         'MAP_A2A_Regrid_201203/MAP_HTAP.nc'

    DATA_DIR_HTAP = TRIM( DATA_DIR_1x1 ) // &
       'HTAP/edgar_HTAP_'


    htap_src_mask = TRIM( DATA_DIR_1x1 ) // &
                   'HTAP/MASKS/HTAP_Phase2_tier1NC01x01_v2.nc'

    ALK4_TMP  = 0d0
    C2H6_TMP  = 0d0
    C3H8_TMP  = 0d0
    PRPE_TMP  = 0d0
    CH2O_TMP  = 0d0
    ALD2_TMP  = 0d0
    ACET_TMP  = 0d0
    MEK_TMP   = 0d0


    ! Loop over species
    DO KLM = 1, SIZE( SPCLIST )

       ! Set GEOS_01x01
       GEOS_01x01 = 0d0

       IF ( ITS_A_FULLCHEM_SIM() ) THEN
          SId = SPCLIST( KLM )
          VId = VAR_NAME( KLM )
       ELSE
          SNo = OFFLINE_ID( KLM )
       ENDIF


       ! DataName for year
       IF (THISYEAR .le. 2009) THEN
          HTAP_YEAR = '_2008_'
          HTAP_YEAR_2 = '_2008'
       ELSE
          HTAP_YEAR = '_2010_'
          HTAP_YEAR_2 = '_2010'
       ENDIF

       SELECT CASE ( THISMONTH )
          CASE ( 1  )
            MONTH = '1'
          CASE ( 2  )
            MONTH = '2'
          CASE ( 3  )
            MONTH = '3'
          CASE ( 4  )
            MONTH = '4'
          CASE ( 5  )
            MONTH = '5'
          CASE ( 6  )
            MONTH = '6'
          CASE ( 7  )
            MONTH = '7'
          CASE ( 8  )
            MONTH = '8'
          CASE ( 9  )
            MONTH = '9'
          CASE ( 10 )
            MONTH = '10'
          CASE ( 11 )
            MONTH = '11'
          CASE ( 12 )
            MONTH = '12'
       END SELECT

       EMIS_TYPE = (/ 'AIR', 'ENERGY', 'INDUSTRY', 'RESIDENTIAL', &
                      'SHIPS', 'TRANSPORT', 'AGRICULTURE' /)

       DO CAT = 1, SIZE( EMIS_TYPE )

          ! Set Mask
          MASK = 0d0

          ! Open model_ready mask from netCDF file
          CALL Ncop_Rd(fId1, TRIM( htap_src_mask ))

          ! Read model_ready data from netCDF file
          CALL NcRd(TMP_MASK, fId1, 'region_code',   &
          (/ 1,  1 /),                               & !Start
          (/ I01x01, J01x01/) )                        !Count lon/lat

          ! Close netCDF file
          CALL NcCl( fId1 )


          ! Apply Source Mask Scaling
          DO I = 1, I01x01

          ! J on mask is N->S, but I on GEOS_01x01 is S->N
          JJ = J01x01

          DO J = 1, J01x01

             IF ( LOCN ) THEN
                IF ( TMP_MASK(I,JJ) .EQ. 2d0  ) MASK(I,J) = 1d0
             ENDIF

             IF ( LNAM ) THEN
                IF ( TMP_MASK(I,JJ) .EQ. 3d0  ) MASK(I,J) = 1d0
             ENDIF

             IF ( LEUR ) THEN
                IF ( TMP_MASK(I,JJ) .EQ. 4d0  ) MASK(I,J) = 1d0
             ENDIF

             IF ( LSAS ) THEN
                IF ( TMP_MASK(I,JJ) .EQ. 5d0  ) MASK(I,J) = 1d0
             ENDIF

             IF ( LEAS ) THEN
                IF ( TMP_MASK(I,JJ) .EQ. 6d0  ) MASK(I,J) = 1d0
             ENDIF

             IF ( LSEA ) THEN
                IF ( TMP_MASK(I,JJ) .EQ. 7d0  ) MASK(I,J) = 1d0
             ENDIF

             IF ( LPAN ) THEN
                IF ( TMP_MASK(I,JJ) .EQ. 8d0  ) MASK(I,J) = 1d0
             ENDIF

             IF ( LNAF ) THEN
                IF ( TMP_MASK(I,JJ) .EQ. 9d0  ) MASK(I,J) = 1d0
             ENDIF

             IF ( LSAF ) THEN
                IF ( TMP_MASK(I,JJ) .EQ. 10d0  ) MASK(I,J) = 1d0
             ENDIF

             IF ( LMDE ) THEN
                IF ( TMP_MASK(I,JJ) .EQ. 11d0  ) MASK(I,J) = 1d0
             ENDIF

             IF ( LMCA ) THEN
                IF ( TMP_MASK(I,JJ) .EQ. 12d0  ) MASK(I,J) = 1d0
             ENDIF

             IF ( LSAM ) THEN
                IF ( TMP_MASK(I,JJ) .EQ. 13d0  ) MASK(I,J) = 1d0
             ENDIF

             IF ( LRBU ) THEN
                IF ( TMP_MASK(I,JJ) .EQ. 14d0  ) MASK(I,J) = 1d0
             ENDIF

             IF ( LCAS ) THEN
                IF ( TMP_MASK(I,JJ) .EQ. 15d0  ) MASK(I,J) = 1d0
             ENDIF

             IF ( LNPO ) THEN
                IF ( TMP_MASK(I,JJ) .EQ. 16d0  ) MASK(I,J) = 1d0
             ENDIF

             IF ( LSPO ) THEN
                IF ( TMP_MASK(I,JJ) .EQ. 17d0  ) MASK(I,J) = 1d0
             ENDIF

             JJ = JJ - 1

          ENDDO
          ENDDO

          IF (.not. LOCN .and. .not. LNAM .and. .not. LEUR .and. &
              .not. LSAS .and. .not. LEAS .and. .not. LSEA .and. &
              .not. LPAN .and. .not. LNAF .and. .not. LSAF .and. &
              .not. LMDE .and. .not. LMCA .and. .not. LSAM .and. &
              .not. LRBU .and. .not. LCAS .and. .not. LNPO .and. &
              .not. LSPO ) MASK = 1d0

          ! Apply Sector Scaling Factor
          IF ( TRIM( EMIS_TYPE(CAT) ) == 'AIR'         ) &
               MASK = MASK * SC_AIR
          IF ( TRIM( EMIS_TYPE(CAT) ) == 'SHIPS'       ) &
               MASK = MASK * SC_SHIPS
          IF ( TRIM( EMIS_TYPE(CAT) ) == 'ENERGY'      ) &
               MASK = MASK * SC_ENERGY
          IF ( TRIM( EMIS_TYPE(CAT) ) == 'INDUSTRY'    ) &
               MASK = MASK * SC_INDUSTRY
          IF ( TRIM( EMIS_TYPE(CAT) ) == 'TRANSPORT'   ) &
               MASK = MASK * SC_TRANSPORT
          IF ( TRIM( EMIS_TYPE(CAT) ) == 'RESIDENTIAL' ) &
               MASK = MASK * SC_RESIDENTIAL
          IF ( TRIM( EMIS_TYPE(CAT) ) == 'AGRICULTURE' ) &
               MASK = MASK * SC_AGRICULTURE

          ! Apply Species Scaling Factor
          IF ( TRIM( SPCLIST(KLM) ) == 'BC'    )  MASK = MASK * SC_BC
          IF ( TRIM( SPCLIST(KLM) ) == 'CO'    )  MASK = MASK * SC_CO
          IF ( TRIM( SPCLIST(KLM) ) == 'OC'    )  MASK = MASK * SC_OC
          IF ( TRIM( SPCLIST(KLM) ) == 'NH3'   )  MASK = MASK * SC_NH3
          IF ( TRIM( SPCLIST(KLM) ) == 'NOx'   )  MASK = MASK * SC_NOX
          IF ( TRIM( SPCLIST(KLM) ) == 'SO2'   )  MASK = MASK * SC_SO2
          IF ( TRIM( SPCLIST(KLM) ) == 'NMVOC' )  MASK = MASK * SC_CH2O

          SELECT CASE ( EMIS_TYPE(CAT) )


             CASE ( 'AIR' )
             IF ( SPCLIST(KLM) == 'NH3' ) CYCLE
             FILENAME =  TRIM(DATA_DIR_HTAP) // TRIM(SPCLIST( KLM ))  &
                         // '_emi_' // TRIM(EMIS_TYPE( CAT ))         &
                         // TRIM(HTAP_YEAR_2)                         &
                         // '.0.1x0.1.nc'

             CASE ( 'SHIPS' )
             IF ( SPCLIST(KLM) == 'NH3' ) CYCLE
             FILENAME =  TRIM(DATA_DIR_HTAP) // TRIM(SPCLIST( KLM ))  &
                         // '_emi_' // TRIM(EMIS_TYPE( CAT ))         &
                         // TRIM(HTAP_YEAR_2)                         &
                         // '.0.1x0.1.nc'

             CASE ( 'AGRICULTURE' )
             IF ( SPCLIST(KLM) .NE. 'NH3' ) CYCLE
             IF ( SPCLIST(KLM) == 'NH3'  )                            &
             FILENAME =  TRIM(DATA_DIR_HTAP) // TRIM(SPCLIST( KLM ))  &
                         // '_emi_' // TRIM(EMIS_TYPE( CAT ))         &
                         // TRIM(HTAP_YEAR) // TRIM( MONTH )          &
                         // '.0.1x0.1.nc'

             CASE DEFAULT

             FILENAME =  TRIM(DATA_DIR_HTAP) // TRIM(SPCLIST( KLM ))  &
                         // '_emi_' // TRIM(EMIS_TYPE( CAT ))         &
                         // TRIM(HTAP_YEAR) // TRIM( MONTH )          &
                        // '.0.1x0.1.nc'
          END SELECT

          IF ( TRIM( SId ) .NE. 'NMVOC') THEN

             ! Echo info
             WRITE( 6, 100 )  TRIM( FILENAME )

             ! Open model_ready data from netCDF file
             CALL Ncop_Rd(fId1, TRIM(FILENAME))

             ! Read model_ready data from netCDF file
             CALL NcRd(ARRAY, fId1, TRIM(VId),   &
                  (/ 1,  1 /),                   & !Start
                  (/ I01x01, J01x01/) )            !Count lat/lon

             ! Close netCDF file
             CALL NcCl( fId1 )

          ENDIF

          IF (      TRIM( SId ) == 'NMVOC'           .and. &
              .not. TRIM( EMIS_TYPE(CAT) ) == 'AIR'  .and. &
              .not. TRIM( EMIS_TYPE(CAT) ) == 'SHIPS' ) THEN

             DO N = 1, NEMANTHRO

                IF ( N == IDEC2H6 ) VOC_FILE = TRIM( DATA_DIR_1x1 )  // &
                      'HTAP/HTAPv2_' // TRIM( EMIS_TYPE(CAT) ) // &
                      '_ethane' // TRIM(HTAP_YEAR) // '0.1x0.1.nc'

                IF ( N == IDEC3H8 ) VOC_FILE = TRIM( DATA_DIR_1x1 )  // &
                      'HTAP/HTAPv2_' // TRIM( EMIS_TYPE(CAT) ) // &
                      '_propane' // TRIM(HTAP_YEAR) // '0.1x0.1.nc'

                IF ( N == IDEPRPE ) VOC_FILE = TRIM( DATA_DIR_1x1 )  // &
                      'HTAP/HTAPv2_' // TRIM( EMIS_TYPE(CAT) ) // &
                      '_propene' // TRIM(HTAP_YEAR) // '0.1x0.1.nc'

                IF ( N == IDECH2O ) VOC_FILE = TRIM( DATA_DIR_1x1 )  // &
                      'HTAP/HTAPv2_' // TRIM( EMIS_TYPE(CAT) ) // &
                      '_formaldehyde' // TRIM(HTAP_YEAR) // '0.1x0.1.nc'

                IF ( N == IDEALD2 ) VOC_FILE = TRIM( DATA_DIR_1x1 )  // &
                      'HTAP/HTAPv2_' // TRIM( EMIS_TYPE(CAT) ) // &
                      '_other_alkanals' // TRIM(HTAP_YEAR) // '0.1x0.1.nc'

                IF ( N == IDEACET ) VOC_FILE = TRIM( DATA_DIR_1x1 )  // &
                      'HTAP/HTAPv2_' // TRIM( EMIS_TYPE(CAT) ) // &
                      '_ketones' // TRIM(HTAP_YEAR) // '0.1x0.1.nc'

                IF ( N == IDEMEK  ) VOC_FILE = TRIM( DATA_DIR_1x1 )  // &
                      'HTAP/HTAPv2_' // TRIM( EMIS_TYPE(CAT) ) // &
                      '_ketones' // TRIM(HTAP_YEAR) // '0.1x0.1.nc'

                IF ( N == IDEALK4 .or. N == IDEACET .or. &
                     N == IDEMEK  .or. N == IDEPRPE .or. &
                     N == IDEC3H8 .or. N == IDECH2O .or. &
                     N == IDEC2H6 .or. N == IDEALD2 ) THEN

                   ! Echo info
                   IF ( N .NE. IDEALK4 ) WRITE( 6, 100 )  TRIM( VOC_FILE )

                   ! Open model_ready data from netCDF file
                   CALL Ncop_Rd(fId1, TRIM(VOC_FILE))

                   SELECT CASE ( EMIS_TYPE(CAT) )

                      CASE ( 'ENERGY' )

                      ! Read model_ready data from netCDF file
                      CALL NcRd(VOC_DATA, fId1, 'emiss_ene',  &
                           (/ 1, 1, 1 /),                     & !Start
                           (/ I01x01, J01x01, 1/) )             !Count lat/lon

                      IF ( N == IDEALK4 ) THEN
                      DO ALK = 1, SIZE( ALK_NUM )

                         ALK_MASK = TRIM( DATA_DIR_1x1 )   // &
                                    'HTAP/HTAPv2_'   // &
                                    TRIM( EMIS_TYPE(CAT) ) // &
                                    TRIM( ALK_NUM( ALK ) ) // &
                                    TRIM(HTAP_YEAR) // '0.1x0.1.nc'

                         ! Echo info
                          WRITE( 6, 100 )  TRIM( ALK_MASK )

                         ! Open model_ready data from netCDF file
                         CALL Ncop_Rd(fId2, TRIM(ALK_MASK))

                         ! Read model_ready data from netCDF file
                         CALL NcRd(VOC_DATA, fId2, 'emiss_ene',  &
                              (/ 1, 1, 1 /),                     & !Start
                              (/ I01x01, J01x01, 1/) )             !Count lat/lon

                         DO I = 1, I01x01
                         DO J = 1, J01x01
                            IF (MASK(I,J) .GT. 0d0 ) THEN

                                 VOC_DATA(I,J,1) = VOC_DATA(I,J,1) * MASK(I,J)

                            ENDIF
                         ENDDO
                         ENDDO
                         ALK4_TMP = ALK4_TMP + VOC_DATA(:,:,1)

                         ! Close netCDF file
                         CALL NcCl( fId2 )

                     ENDDO
                     ENDIF


                      CASE ( 'TRANSPORT' )

                      ! Read model_ready data from netCDF file
                      CALL NcRd(VOC_DATA, fId1, 'emiss_tra',  &
                           (/ 1, 1, 1 /),                     & !Start
                           (/ I01x01, J01x01, 1/) )             !Count lat/lon

                      IF ( N == IDEALK4 ) THEN
                      DO ALK = 1, SIZE( ALK_NUM )

                          ALK_MASK = TRIM( DATA_DIR_1x1 )   // &
                                    'HTAP/HTAPv2_'   // &
                                    TRIM( EMIS_TYPE(CAT) ) // &
                                    TRIM( ALK_NUM( ALK ) ) // &
                                    TRIM(HTAP_YEAR) // '0.1x0.1.nc'

                         ! Echo info
                         IF ( N == IDEALK4 ) WRITE( 6, 100 )  TRIM( ALK_MASK )

                         ! Open model_ready data from netCDF file
                         CALL Ncop_Rd(fId2, TRIM(ALK_MASK))

                         ! Read model_ready data from netCDF file
                         CALL NcRd(VOC_DATA, fId2, 'emiss_tra',  &
                              (/ 1, 1, 1 /),                     & !Start
                              (/ I01x01, J01x01, 1/) )             !Count lat/lon

                         DO I = 1, I01x01
                         DO J = 1, J01x01
                            IF (MASK(I,J) .GT. 0d0 ) THEN

                                 VOC_DATA(I,J,1) = VOC_DATA(I,J,1) * MASK(I,J)

                            ENDIF
                         ENDDO
                         ENDDO
                         ALK4_TMP = ALK4_TMP + VOC_DATA(:,:,1)

                         ! Close netCDF file
                         CALL NcCl( fId2 )

                      ENDDO
                      ENDIF


                      CASE ( 'RESIDENTIAL' )

                      ! Read model_ready data from netCDF file
                      CALL NcRd(VOC_DATA, fId1, 'emiss_dom',  &
                           (/ 1, 1, 1 /),                     & !Start
                           (/ I01x01, J01x01, 1/) )             !Count lat/lon

                      IF ( N == IDEALK4 ) THEN
                      DO ALK = 1, SIZE( ALK_NUM )

                          ALK_MASK = TRIM( DATA_DIR_1x1 )   // &
                                    'HTAP/HTAPv2_'   // &
                                    TRIM( EMIS_TYPE(CAT) ) // &
                                    TRIM( ALK_NUM( ALK ) ) // &
                                    TRIM(HTAP_YEAR) // '0.1x0.1.nc'

                         ! Echo info
                         IF ( N == IDEALK4 ) WRITE( 6, 100 )  TRIM( ALK_MASK )

                         ! Open model_ready data from netCDF file
                         CALL Ncop_Rd(fId2, TRIM(ALK_MASK))

                         ! Read model_ready data from netCDF file
                         CALL NcRd(VOC_DATA, fId2, 'emiss_dom',  &
                              (/ 1, 1, 1 /),                     & !Start
                              (/ I01x01, J01x01, 1/) )             !Count lat/lon

                         DO I = 1, I01x01
                         DO J = 1, J01x01
                            IF (MASK(I,J) .GT. 0d0 ) THEN

                                 VOC_DATA(I,J,1) = VOC_DATA(I,J,1) * MASK(I,J)

                            ENDIF
                         ENDDO
                         ENDDO
                         ALK4_TMP = ALK4_TMP + VOC_DATA(:,:,1)

                         ! Close netCDF file
                         CALL NcCl( fId2 )

                      ENDDO
                      ENDIF


                      CASE ( 'INDUSTRY'  )

                      ! Read model_ready data from netCDF file
                      CALL NcRd(VOC_DATA, fId1, 'emiss_ind',  &
                           (/ 1, 1, 1 /),                     & !Start
                           (/ I01x01, J01x01, 1/) )             !Count lat/lon


                      IF ( N == IDEALK4 ) THEN
                         DO ALK = 1, SIZE( ALK_NUM )

                            ALK_MASK = TRIM( DATA_DIR_1x1 )   // &
                                    'HTAP/HTAPv2_'   // &
                                    TRIM( EMIS_TYPE(CAT) ) // &
                                    TRIM( ALK_NUM( ALK ) ) // &
                                    TRIM(HTAP_YEAR) // '0.1x0.1.nc'

                           ! Echo info
                           IF ( N == IDEALK4 ) WRITE( 6, 100 )  TRIM( ALK_MASK )

                           ! Open model_ready data from netCDF file
                           CALL Ncop_Rd(fId2, TRIM(ALK_MASK))

                           ! Read model_ready data from netCDF file
                           CALL NcRd(VOC_DATA, fId2, 'emiss_ind',  &
                                (/ 1, 1, 1 /),                     & !Start
                                (/ I01x01, J01x01, 1/) )             !Count lat/lon

                           DO I = 1, I01x01
                           DO J = 1, J01x01
                              IF (MASK(I,J) .GT. 0d0 ) THEN

                                 VOC_DATA(I,J,1) = VOC_DATA(I,J,1) * MASK(I,J)

                              ENDIF
                           ENDDO
                           ENDDO
                           ALK4_TMP = ALK4_TMP + VOC_DATA(:,:,1)

                           ! Close netCDF file
                           CALL NcCl( fId2 )

                         ENDDO
                      ENDIF

                      CASE DEFAULT

                      !Do nothing

                   END SELECT

                   ! Close netCDF file
                   CALL NcCl( fId1 )

                ENDIF

                DO I = 1, I01x01
                DO J = 1, J01x01
                   IF (MASK(I,J) .GT. 0d0 ) THEN

                      VOC_DATA(I,J,1) = VOC_DATA(I,J,1) * MASK(I,J)

                   ENDIF
                ENDDO
                ENDDO

                !Apply VOCs Speciation
                IF ( N == IDEC2H6 ) C2H6_TMP  = C2H6_TMP + VOC_DATA(:,:,1)
                IF ( N == IDEC3H8 ) C3H8_TMP  = C3H8_TMP + VOC_DATA(:,:,1)
                IF ( N == IDEPRPE ) PRPE_TMP  = PRPE_TMP + VOC_DATA(:,:,1)
                IF ( N == IDECH2O ) CH2O_TMP  = CH2O_TMP + VOC_DATA(:,:,1)
                IF ( N == IDEALD2 ) ALD2_TMP  = ALD2_TMP + VOC_DATA(:,:,1)
                IF ( N == IDEACET ) ACET_TMP  = ACET_TMP + VOC_DATA(:,:,1)
                IF ( N == IDEMEK  ) MEK_TMP   = MEK_TMP  + VOC_DATA(:,:,1)

             ENDDO

          ENDIF

          ! Apply Regional Mask if defined
          DO I = 1, I01x01
          DO J = 1, J01x01
             IF (MASK(I,J) .GT. 0d0 ) THEN

                ARRAY(I,J)    = ARRAY(I,J)    * MASK(I,J)

             ENDIF
          ENDDO
          ENDDO

          ! Add sectors before regridding
          GEOS_01x01(:,:) = GEOS_01x01(:,:) + ARRAY(:,:)

       ENDDO

100    FORMAT( '     - EMISS_HTAP_0.1x0.1:  &
            Reading : ', a )

       ! Regrid from GEOS 0.1x0.1 --> current model resolution
       SELECT CASE ( SId )

          CASE ( 'BC' )

                !-----------------
                ! BCPO
                !-----------------

                INGRID  => GEOS_01x01(:,:)
                OUTGRID => BC(:,:)

                ! Regrid
                CALL DO_REGRID_DKH( LLFILENAME, I01x01,    J01x01,  &
                                    INGRID,     OUTGRID, IS_MASS=0, &
                                    netCDF=.TRUE.)

                ! Free pointers
                NULLIFY( INGRID, OUTGRID )

          CASE ( 'CO' )

                !-----------------
                ! CO
                !-----------------
                ! Point to array slices
                INGRID  => GEOS_01x01(:,:)
                OUTGRID => CO(:,:)

                ! Regrid
                CALL DO_REGRID_DKH( LLFILENAME, I01x01,    J01x01, &
                                    INGRID,     OUTGRID, IS_MASS=0, &
                                    netCDF=.TRUE. )

                ! Free pointers
                 NULLIFY( INGRID, OUTGRID )

          CASE ( 'NH3' )

                !-----------------
                ! NH3
                !-----------------

                ! Point to array slices
                INGRID  => GEOS_01x01(:,:)
                OUTGRID => NH3(:,:)

                ! Regrid
                CALL DO_REGRID_DKH( LLFILENAME, I01x01,    J01x01, &
                                    INGRID,     OUTGRID, IS_MASS=0, &
                                    netCDF=.TRUE. )

                ! Free pointers
                NULLIFY( INGRID, OUTGRID )

          CASE ( 'NMVOC' )

                !-----------------
                ! VOC
                !-----------------

                ! Regrid
                CALL DO_REGRID_DKH( LLFILENAME, I01x01,    J01x01, &
                                    ALK4_TMP,     ALK4_HTAP, IS_MASS=0, &
                                    netCDF=.TRUE. )

                ! Regrid
                CALL DO_REGRID_DKH( LLFILENAME, I01x01,    J01x01, &
                                    ACET_TMP,     ACET_HTAP, IS_MASS=0, &
                                    netCDF=.TRUE. )

                ! Regrid
                CALL DO_REGRID_DKH( LLFILENAME, I01x01,    J01x01, &
                                    MEK_TMP,      MEK_HTAP,  IS_MASS=0, &
                                    netCDF=.TRUE. )

                ! Regrid
                CALL DO_REGRID_DKH( LLFILENAME, I01x01,    J01x01, &
                                    ALD2_TMP,     ALD2_HTAP, IS_MASS=0, &
                                    netCDF=.TRUE. )

                ! Regrid
                CALL DO_REGRID_DKH( LLFILENAME, I01x01,    J01x01, &
                                    PRPE_TMP,     PRPE_HTAP, IS_MASS=0, &
                                    netCDF=.TRUE. )

                ! Regrid
                CALL DO_REGRID_DKH( LLFILENAME, I01x01,    J01x01, &
                                    C3H8_TMP,     C3H8_HTAP, IS_MASS=0, &
                                    netCDF=.TRUE. )

                ! Regrid
                CALL DO_REGRID_DKH( LLFILENAME, I01x01,    J01x01, &
                                    CH2O_TMP,     CH2O_HTAP, IS_MASS=0, &
                                    netCDF=.TRUE. )

                ! Regrid
                CALL DO_REGRID_DKH( LLFILENAME, I01x01,    J01x01, &
                                    C2H6_TMP,     C2H6_HTAP, IS_MASS=0, &
                                    netCDF=.TRUE. )

          CASE ( 'NOx' )

                !-----------------
                ! NOX
                !-----------------
                ! Point to array slices
                INGRID  => GEOS_01x01(:,:)
                OUTGRID => NOX(:,:)

                ! Regrid
                CALL DO_REGRID_DKH( LLFILENAME, I01x01,    J01x01,  &
                                    INGRID,     OUTGRID, IS_MASS=0, &
                                    netCDF=.TRUE. )

                ! Free pointers
                NULLIFY( INGRID, OUTGRID )

          CASE ( 'OC' )

                !-----------------
                ! OCPO
                !-----------------

                INGRID  => GEOS_01x01(:,:)
                OUTGRID => OC(:,:)

                ! Regrid
                CALL DO_REGRID_DKH( LLFILENAME, I01x01,    J01x01, &
                                    INGRID,     OUTGRID, IS_MASS=0, &
                                    netCDF=.TRUE. )

                ! Free pointers
                NULLIFY( INGRID, OUTGRID )
                ! Reset GEOS_01x01
                GEOS_01x01 = 0d0

          CASE ( 'SO2' )

                !-----------------
                ! SO2
                !-----------------

                ! Point to array slices
                INGRID  => GEOS_01x01(:,:)
                OUTGRID => SO2(:,:)

                ! Regrid
                CALL DO_REGRID_DKH( LLFILENAME, I01x01,    J01x01,  &
                                    INGRID,     OUTGRID, IS_MASS=0, &
                                    netCDF=.TRUE. )

                ! Free pointers
                NULLIFY( INGRID, OUTGRID )

          CASE ( 'PM2.5' )

                !-----------------
                ! PM2.5
                !-----------------

                ! Point to array slices
                INGRID  => GEOS_01x01(:,:)
                OUTGRID => PM25(:,:)

                ! Regrid
                CALL DO_REGRID_DKH( LLFILENAME, I01x01,    J01x01,  &
                                    INGRID,     OUTGRID, IS_MASS=0, &
                                    netCDF=.TRUE. )

                ! Free pointers
                NULLIFY( INGRID, OUTGRID )


          END SELECT

    ENDDO

    !--------------------------
    ! Print emission totals for the day
    !--------------------------
    CALL TOTAL_HTAP_Tg( DOYT )

    ! Return to calling program
  END SUBROUTINE EMISS_HTAP
  !
  ! !IROUTINE: total_htap_Tg
  !
  ! !DESCRIPTION: Subroutine TOTAL\_ANTHRO\_TG prints the totals for the
  !  anthropogenic emissions of NOx, CO, SO2 and NH3.
  !\\
  !\\
  ! !INTERFACE:
  !
  SUBROUTINE TOTAL_HTAP_TG( DAY )
    !
    ! !USES:
    !
    USE GRID_MOD,     ONLY : GET_AREA_M2
    USE TIME_MOD,     ONLY : GET_YEAR, GET_MONTH
    USE TRACER_MOD,   ONLY : XNUMOL
    USE TRACERID_MOD, ONLY : IDTCO, IDTCH2O
    USE TRACERID_MOD, ONLY : IDTSO2, IDTNH3, IDTNOx
    USE TRACERID_MOD, ONLY : IDTOCPO,  IDTBCPO

#     include "CMN_SIZE"          ! Size parameters

    ! !INPUT PARAMETERS:
    !
    INTEGER, INTENT(IN) :: DAY   ! Day of data to compute totals
    !
    ! !REVISION HISTORY:
    !   Sept 2013 - Yanko Davila - initial version
    !EOP
    !------------------------------------------------------------------------------
    !BOC
    !
    ! !LOCAL VARIABLES:
    !
    INTEGER             :: II, JJ, IH, LL
    REAL*8              :: T_CO,   T_NOx, T_SO2, T_NH3
    REAL*8              :: T_BC,  T_OC
    REAL*8              :: T_PM25
    REAL*8              :: T_ALK4, T_ACET, T_MEK,  T_ALD2
    REAL*8              :: T_PRPE, T_C3H8, T_CH2O, T_C2H6

    REAL*8              :: tmpArea(IIPAR, JJPAR)
    CHARACTER(LEN=3)    :: UNIT
    REAL*8,  PARAMETER  :: SEC_IN_DAY  = 86400d0
    ! Days per month
    REAL*8                :: DAYS_IN_MONTH
    REAL*8                :: DMON(12)

    !=================================================================
    ! TOTAL_HTAP_TG begins here!
    !=================================================================

    ! Fancy output
    WRITE( 6, '(a)' ) REPEAT( '=', 79 )
    WRITE( 6, 100  )
100 FORMAT( 'H. T. A. P.  E M I S S I O N S', / )

    DO II = 1, IIPAR
       DO JJ = 1, JJPAR
          tmpArea(II,JJ) = GET_AREA_M2(JJ)
       ENDDO
    ENDDO

    IF ( GET_YEAR()  .le. 2009 )  DMON = (/ 31d0, 29d0, 31d0, 30d0,  &
                                            31d0, 30d0, 31d0, 31d0,  &
                                            30d0, 31d0, 30d0, 31d0 /)

    IF ( GET_YEAR()  .ge. 2010 )  DMON = (/ 31d0, 28d0, 31d0, 30d0,  &
                                            31d0, 30d0, 31d0, 31d0,  &
                                            30d0, 31d0, 30d0, 31d0 /)

    ! Get days per month
    DAYS_IN_MONTH = DMON( GET_MONTH() )

    ! Total CO  [Tg CO]
    T_CO  = SUM( CO * tmpArea )  * &
         SEC_IN_DAY * DAYS_IN_MONTH * 1d-9

    ! Total NOX [Tg NO2]
    T_NOx = SUM( NOX * tmpArea) * &
         SEC_IN_DAY * DAYS_IN_MONTH * 1d-9

    ! Total SO2 [Tg SO2]
    T_SO2 = SUM( SO2 * tmpArea ) * &
         SEC_IN_DAY * DAYS_IN_MONTH * 1d-9

    ! Total NH3 [Tg NH3]
    T_NH3 = SUM( NH3 * tmpArea ) * &
         SEC_IN_DAY * DAYS_IN_MONTH * 1d-9

    ! Total BC [Tg C]
    T_BC = SUM( BC  * tmpArea ) * &
         SEC_IN_DAY * DAYS_IN_MONTH * 1d-9

    ! Total OC [Tg C]
    T_OC = SUM( OC * tmpArea ) * &
         SEC_IN_DAY * DAYS_IN_MONTH * 1d-9

    ! Total VOC [Tg C]
    T_ALK4 = SUM( ALK4_HTAP * tmpArea ) * &
         SEC_IN_DAY * DAYS_IN_MONTH * 1d-9

    ! Total VOC [Tg C]
    T_ACET = SUM( ACET_HTAP * tmpArea ) * &
         SEC_IN_DAY * DAYS_IN_MONTH * 1d-9 * 0.75d0

    ! Total VOC [Tg C]
    T_MEK  = SUM( MEK_HTAP  * tmpArea ) * &
         SEC_IN_DAY * DAYS_IN_MONTH * 1d-9 * 0.25d0

    ! Total VOC [Tg C]
    T_ALD2 = SUM( ALD2_HTAP * tmpArea ) * &
         SEC_IN_DAY * DAYS_IN_MONTH * 1d-9

    ! Total VOC [Tg C]
    T_PRPE = SUM( PRPE_HTAP * tmpArea ) * &
         SEC_IN_DAY * DAYS_IN_MONTH * 1d-9

    ! Total VOC [Tg C]
    T_C3H8 = SUM( C3H8_HTAP * tmpArea ) * &
         SEC_IN_DAY * DAYS_IN_MONTH * 1d-9

    ! Total VOC [Tg C]
    T_CH2O = SUM( CH2O_HTAP * tmpArea ) * &
         SEC_IN_DAY * DAYS_IN_MONTH * 1d-9

    ! Total VOC [Tg C]
    T_C2H6 = SUM( C2H6_HTAP * tmpArea ) * &
         SEC_IN_DAY * DAYS_IN_MONTH * 1d-9

    ! Total PM25 [Tg C]
    T_PM25 = SUM( PM25 * tmpArea ) * &
         SEC_IN_DAY * DAYS_IN_MONTH * 1d-9


    ! Print totals in [Tg]
    WRITE( 6, 110 ) 'CO    ', T_CO,   '[Tg CO]'
    WRITE( 6, 110 ) 'NOx   ', T_NOx,  '[Tg NO2]'
    WRITE( 6, 110 ) 'SO2   ', T_SO2,  '[Tg SO2]'
    WRITE( 6, 110 ) 'NH3   ', T_NH3,  '[Tg NH3]'
    WRITE( 6, 110 ) 'BC    ', T_BC,   '[Tg C]'
    WRITE( 6, 110 ) 'OC    ', T_OC,   '[Tg C]'
    WRITE( 6, 110 ) 'PM2.5 ', T_PM25, '[Tg C]'
    WRITE( 6, 110 ) 'ALK4  ', T_ALK4, '[Tg Molec]'
    WRITE( 6, 110 ) 'ACET  ', T_ACET, '[Tg Molec]'
    WRITE( 6, 110 ) 'MEK   ', T_MEK , '[Tg Molec]'
    WRITE( 6, 110 ) 'ALD2  ', T_ALD2, '[Tg Molec]'
    WRITE( 6, 110 ) 'PRPE  ', T_PRPE, '[Tg Molec]'
    WRITE( 6, 110 ) 'C3H8  ', T_C3H8, '[Tg Molec]'
    WRITE( 6, 110 ) 'CH2O  ', T_CH2O, '[Tg Molec]'
    WRITE( 6, 110 ) 'C2H6  ', T_C2H6, '[Tg Molec]'

    ! Format statement
110 FORMAT( 'HTAP Emissions ', a6, &
            ': ', f20.8, 1x, a10 )

    ! Fancy output
    WRITE( 6, '(a)' ) REPEAT( '=', 79 )

    ! Return to calling program
  END SUBROUTINE TOTAL_HTAP_Tg
  !
  ! !IROUTINE: init_htap
  !
  ! !DESCRIPTION: Subroutine INIT\_HTAP allocates and zeroes all
  !  module arrays.
  !\\
  !\\
  ! !INTERFACE:
  !
  SUBROUTINE INIT_HTAP
    !
    ! !USES:
    !
    USE ERROR_MOD,       ONLY : ALLOC_ERR
    USE LOGICAL_MOD,     ONLY : LHTAP

#     include "CMN_SIZE"          ! Size parameters!
    !EOP
    !------------------------------------------------------------------------------
    !BOC
    !
    ! !LOCAL VARIABLES:
    !
    INTEGER              :: AS, J

    !=================================================================
    ! INIT_NEI2008_ANTHRO begins here!
    !=================================================================

    ! Return if LHTAP is false
    IF ( .not. LHTAP ) RETURN

    !--------------------------------------------------
    ! Allocate and zero arrays for emissions
    !--------------------------------------------------

    ALLOCATE( BC( IIPAR, JJPAR ), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'BC' )
    BC = 0d0

    ALLOCATE( CO( IIPAR, JJPAR ), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'CO' )
    CO = 0d0

    ALLOCATE( NH3( IIPAR, JJPAR ), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'NH3' )
    NH3 = 0d0

    ALLOCATE( NOX( IIPAR, JJPAR ), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'NOX' )
    NOX = 0d0

    ALLOCATE( OC( IIPAR, JJPAR ), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'OC' )
    OC = 0d0

    ALLOCATE( SO2( IIPAR, JJPAR ), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'SO2' )
    SO2 = 0d0

    ALLOCATE( PM25( IIPAR, JJPAR ), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'PM25' )
    PM25 = 0d0

    ALLOCATE( ALK4_HTAP( IIPAR, JJPAR ), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'ALK4' )
    ALK4_HTAP = 0d0

    ALLOCATE( ACET_HTAP( IIPAR, JJPAR ), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'ACET' )
    ACET_HTAP = 0d0

    ALLOCATE( MEK_HTAP ( IIPAR, JJPAR ), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'MEK'  )
    MEK_HTAP = 0d0

    ALLOCATE( ALD2_HTAP( IIPAR, JJPAR ), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'ALD2' )
    ALD2_HTAP = 0d0

    ALLOCATE( PRPE_HTAP( IIPAR, JJPAR ), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'PRPE' )
    PRPE_HTAP = 0d0

    ALLOCATE( C3H8_HTAP( IIPAR, JJPAR ), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'C3H8' )
    C3H8_HTAP = 0d0

    ALLOCATE( CH2O_HTAP( IIPAR, JJPAR ), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'CH2O' )
    CH2O_HTAP = 0d0

    ALLOCATE( C2H6_HTAP( IIPAR, JJPAR ), STAT=AS )
    IF ( AS /= 0 ) CALL ALLOC_ERR( 'C2H6' )
    C2H6_HTAP = 0d0


    ! Return to calling program
  END SUBROUTINE INIT_HTAP
  !
  ! !IROUTINE: cleanup_htap
  !
  ! !DESCRIPTION: Subroutine CLEANUP\_HTAP deallocates all module
  !  arrays.
  !\\
  !\\
  ! !INTERFACE:
  !
  SUBROUTINE CLEANUP_HTAP
    !
    ! !REVISION HISTORY:
    !  Sept 2013 - Yanko Davila - Initial revision
    !EOP
    !------------------------------------------------------------------------------
    !BOC
    !=================================================================
    ! CLEANUP_NEIO2008_ANTHRO begins here!
    !=================================================================

    IF ( ALLOCATED( BC     ) ) DEALLOCATE( BC    )
    IF ( ALLOCATED( CO     ) ) DEALLOCATE( CO    )
    IF ( ALLOCATED( NH3    ) ) DEALLOCATE( NH3   )
    IF ( ALLOCATED( NOX    ) ) DEALLOCATE( NOX   )
    IF ( ALLOCATED( OC     ) ) DEALLOCATE( OC    )
    IF ( ALLOCATED( SO2    ) ) DEALLOCATE( SO2   )
    IF ( ALLOCATED( PM25   ) ) DEALLOCATE( PM25  )
    IF ( ALLOCATED( ALK4_HTAP   ) ) DEALLOCATE( ALK4_HTAP  )
    IF ( ALLOCATED( ACET_HTAP   ) ) DEALLOCATE( ACET_HTAP  )
    IF ( ALLOCATED( MEK_HTAP    ) ) DEALLOCATE( MEK_HTAP   )
    IF ( ALLOCATED( ALD2_HTAP   ) ) DEALLOCATE( ALD2_HTAP  )
    IF ( ALLOCATED( PRPE_HTAP   ) ) DEALLOCATE( PRPE_HTAP  )
    IF ( ALLOCATED( C3H8_HTAP   ) ) DEALLOCATE( C3H8_HTAP  )
    IF ( ALLOCATED( CH2O_HTAP   ) ) DEALLOCATE( CH2O_HTAP  )
    IF ( ALLOCATED( C2H6_HTAP   ) ) DEALLOCATE( C2H6_HTAP  )


  END SUBROUTINE CLEANUP_HTAP

  !EOC
END MODULE HTAP_MOD

