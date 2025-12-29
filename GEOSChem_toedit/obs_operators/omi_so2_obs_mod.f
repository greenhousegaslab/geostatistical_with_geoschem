! $Id: omi_so2_obs_mod.f
      MODULE OMI_SO2_OBS_MOD
!
!*****************************************************************************
! MODULE OMI_SO2_OBS_MOD contians subroutines necessary to
! 1. Read OMI L3 SO2 observations
! 2. Compute OMI-GEOS-Chem difference, cost function, and adjoint
! forcing
!
!    (ywang (yi.wang@huskers.unl.edu), 07/16/04)
!
!  Module Variables:
!  ===========================================================================
!
!  Module Routines:
!  ===========================================================================
!  ( 1) CALC_OMI_SO2_FORCE
!  ( 2) CHECK                : Check status for calling netCDF
!  ( 1) INIT_OMI_SO2_OBS     : Initialize OMI SO2 observation operator
!  ( 3) WRITE_GC_OMI_SO2_OBS
!  ( 4) READ_GC_OMI_SO2_OBS
!  ( 5) MAKE_CURRENT_OMI_SO2 :
!  (  ) MAKE_AVERAGE_OMI_SO2 :
!
!  ===========================================================================
!  NOTES:
!  ( 1) The original OMI L3 SO2 HDFEOS5 data are preprocessed by
!  IDL code OMI_SO2_L3_preprocess.pro. Quality control is done through
!  the IDL code.
!  ( 2) OMI L3 SO2 contains ColumnAmountSO2_PBL (Center Mass of
!  Alitude = 0.9km), ColumnAmountSO2_TRL (CMA = 2.5km), ColumnAmountSO2_TRM
!  (CMA = 7.5km), and ColumnAmountSO2_STL (CMA = 17km). Only the first
!  three are used in the observation operatpor. The prior GEOS-Chem CMA
!  (where SO2 peaks in the vertical column) is calculated. We choose the
!  observation whose CMA is closest to GEOS-Chem CMA to be assimilated.
!
!*****************************************************************************
!

      IMPLICIT NONE

      ! Header files
#     include "define.h"
#     include "CMN_SIZE"       ! Size parameters

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep centain internal variables
      ! and routines from being seen outside "omi_so2_obs_mod.f"
      !=================================================================

      ! Make everything PRIVATE...
      PRIVATE

      ! ... except these routines
      PUBLIC :: CALC_OMI_SO2_FORCE

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================

      ! Parameters
      CHARACTER(LEN=255) :: FILENAME = 'OMI_L3_SO2_grid_YYYYMMDD.nc'
      INTEGER, PARAMETER :: TIME_WINDOW = 60 ! (units: min)

      ! Variables

      ! Recored to store OMI SO2 observations
      TYPE OMI_SO2_OBS
         REAL    :: LON           ! Latitude (units: degree)
         REAL    :: LAT           ! Longitude (units: dgeree)
         REAL*8  :: TIME          ! Time at start of scan (TAI93) (units: s)
         REAL*8  :: SO2_PBL       ! ColumnAmountSO2_PBL (units: DU)
         REAL*8  :: SO2_TRL       ! ColumnAmountSO2_TRL (units: DU)
         REAL*8  :: SO2_TRM       ! ColumnAmountSO2_TRM (units: DU)
         REAL*8  :: SO2           ! OMI SO2 whose CMA is closest to GEOS-Chem CMA (units: DU)
         REAL*8  :: ERR_PBL       ! SO2_PBL error (units: DU)
         REAL*8  :: ERR_TRL       ! SO2_TRL error (units: DU)
         REAL*8  :: ERR_TRM       ! SO2_TRM error (units: DU)
         REAL*8  :: ERR           ! SO2 error (units: DU)
!         INTEGER :: MARK          ! Mark of which observation is expected. 0: No observation, 1: PBL, 2: TRL, 3: TRM
         INTEGER :: OPT           ! Mark of which observation whose CMA is close to GEOS-Chem CMA, though the observation may not exists
         INTEGER :: FLAG          ! Does the observation expected really exits?
      ENDTYPE OMI_SO2_OBS

      TYPE(OMI_SO2_OBS), ALLOCATABLE :: OMI_SO2(:)

      LOGICAL,           ALLOCATABLE :: FLAGS(:)

      REAL*8                      :: CURR_GC_SO2(IIPAR, JJPAR)
      REAL*8                      :: CURR_PBL_SO2(IIPAR, JJPAR)
      REAL*8                      :: CURR_TRL_SO2(IIPAR, JJPAR)
      REAL*8                      :: CURR_TRM_SO2(IIPAR, JJPAR)
      REAL*8                      :: CURR_SO2(IIPAR, JJPAR)
      REAL*8                      :: CURR_DIFF_SO2(IIPAR, JJPAR)
      REAL*8                      :: CURR_FORCING(IIPAR, JJPAR)
      REAL*8                      :: CURR_COST(IIPAR, JJPAR)
      REAL*8                      :: CURR_CMA(IIPAR, JJPAR)
      INTEGER                     :: CURR_COUNT(IIPAR, JJPAR) ! number of observations in current time window
      INTEGER                     :: CURR_PBL_C(IIPAR, JJPAR) ! number of PBL observations used in current time window
      INTEGER                     :: CURR_TRL_C(IIPAR, JJPAR) ! number of TRL observations used in current time window
      INTEGER                     :: CURR_TRM_C(IIPAR, JJPAR) ! numberof TRM observations used in current time window
      INTEGER                     :: CURR_OPT(IIPAR, JJPAR)   !

      REAL*8                      :: ALL_GC_SO2(IIPAR, JJPAR)   = 0D0
      REAL*8                      :: ALL_PBL_SO2(IIPAR, JJPAR)  = 0D0
      REAL*8                      :: ALL_TRL_SO2(IIPAR, JJPAR)  = 0D0
      REAL*8                      :: ALL_TRM_SO2(IIPAR, JJPAR)  = 0D0
      REAL*8                      :: ALL_SO2(IIPAR, JJPAR)      = 0D0
      REAL*8                      :: ALL_DIFF_SO2(IIPAR, JJPAR) = 0D0
      REAL*8                      :: ALL_FORCING(IIPAR, JJPAR)  = 0D0
      REAL*8                      :: ALL_COST(IIPAR, JJPAR)     = 0D0
      REAL*8                      :: ALL_CMA(IIPAR, JJPAR)      = 0D0
      INTEGER                     :: ALL_COUNT(IIPAR, JJPAR)    = 0 ! number of observations in simulation time
      INTEGER                     :: ALL_PBL_C(IIPAR, JJPAR)    = 0
      INTEGER                     :: ALL_TRL_C(IIPAR, JJPAR)    = 0
      INTEGER                     :: ALL_TRM_C(IIPAR, JJPAR)    = 0
      INTEGER                     :: ALL_OPT(IIPAR, JJPAR)      = 0 ! number of "best" observations used

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAIN" statement
      !=================================================================
      CONTAINS
!
!-----------------------------------------------------------------------------
!
      SUBROUTINE INIT_OMI_SO2_OBS
!
!*****************************************************************************
!  Subroutine INIT_OMI_SO2_OBS initialize the OMI SO2 observation
!   ( 1) Check if SO2 observation is specified by the adjoint input
!   files
!  (ywang, 07/16/14)
!*****************************************************************************
!
      ! Reference to f90 modules
      USE TRACER_MOD,         ONLY : TRACER_NAME
      USE TRACERID_MOD,       ONLY : IDTSO2
      USE ADJ_ARRAYS_MOD,     ONLY : OBS_THIS_TRACER
      USE ERROR_MOD,          ONLY : ERROR_STOP

      !=================================================================
      ! INIT_OMI_SO2_OBS begins here!
      !=================================================================

      IF ( OBS_THIS_TRACER(IDTSO2) ) THEN

         WRITE( 6, 100 ) IDTSO2, TRACER_NAME(IDTSO2)

      ELSE

         CALL ERROR_STOP( 'Error: Obs SO2 tracer is not specified in
     &OBSERVATION MENU in input.gcadj',
     &                    'INIT_OMI_SO2_OBS (omi_so2_obs_mod.f)' )

      END IF

 100  FORMAT( 3X, 'Tracer ID: ', I4, 'Use OMI L3 ', A6 )

      ! Return to the calling routine
      END SUBROUTINE INIT_OMI_SO2_OBS
!
!-----------------------------------------------------------------------------
!
      SUBROUTINE READ_GC_OMI_SO2_OBS( YYYYMMDD, N_SO2 )
!
!*****************************************************************************
!  Subroutine READ_OMI_SO2_OBS read OMI SO2 data preprocessed by
!  OMI_SO2_L3_preprocess.pro if N_CALC > 1 (ref: note 1) or
!  observations whose CMA are close to GEOS-Chem CMA if N_CALC == 1
!  (ref: note 2)
!  (ywang, 07/16/14)
!
!  Arguements as Input:
!  ===========================================================================
!  ( 1) YYYYMMDD    (INTEGER) : Current year-month-day
!  Arguements as Output:
!  ===========================================================================
!  ( 1) N_SO2       (INTEGER) : Number of OMI SO2 observations for
!  current day
!
!  Module variable as Output:
!  ===========================================================================
!  ( 1) OMI_SO2 (OMI_SO2_OBS) : OMI SO2 observations
!
!*****************************************************************************
!
      ! Reference to f90 modules
      USE ADJ_ARRAYS_MOD,   ONLY : N_CALC
      USE NETCDF
      USE TIME_MOD,         ONLY : EXPAND_DATE, GET_NYMD

      ! Arguements
      INTEGER, INTENT( IN)   :: YYYYMMDD
      INTEGER, INTENT(OUT)   :: N_SO2

      ! Local variables
      INTEGER                :: FID,        N_ID
      INTEGER                :: LON_ID,     LAT_ID,     TIME_ID
      INTEGER                :: SO2_PBL_ID, SO2_TRL_ID, SO2_TRM_ID
      INTEGER                :: SO2_ID
      INTEGER                :: ERR_PBL_ID, ERR_TRL_ID, ERR_TRM_ID
      INTEGER                :: ERR_ID
!      INTEGER                :: MARK_ID
      INTEGER                :: OPT_ID  ! (ywang, 09/12/14)
      INTEGER                :: FLAG_ID ! (ywang, 09/12/14)
      CHARACTER(LEN=255)     :: DIR
      CHARACTER(LEN=255)     :: READ_FILENAME
      CHARACTER(LEN=4)       :: TMP
      INTEGER, ALLOCATABLE   :: TMPINT(:)
      REAL*4,  ALLOCATABLE   :: TMP4(:)
      REAL*8,  ALLOCATABLE   :: TMP8(:)
      LOGICAL                :: LF

      !=================================================================
      ! READ_OMI_SO2_OBS begins here!
      !=================================================================

      ! Filename root
      IF (N_CALC ==1 ) THEN

         READ_FILENAME = TRIM( FILENAME )

      ELSE IF (N_CALC > 1) THEN

         READ_FILENAME = 'GC_' // TRIM( FILENAME )

      ELSE

         STOP

      END IF

      ! Expand date tokens in filename
      CALL EXPAND_DATE( READ_FILENAME, YYYYMMDD, 9999 )

      ! Construct complete filename
      DIR           = './data/OMI_SO2/'
      READ_FILENAME = TRIM( DIR ) // TRIM( READ_FILENAME )

      ! Does data file exist? If not, it means no data in the day.
      ! (ywang, 09/23/2014)
      INQUIRE( FILE = TRIM( READ_FILENAME ), EXIST = LF )
      IF ( .NOT. LF ) THEN

         ! No data
         N_SO2 = 0

         PRINT*, ' - READ_GC_OMI_SO2_OBS: No data file (warning)'

         RETURN

      END IF



      ! Print to screen
      WRITE(6, 100) TRIM( READ_FILENAME )
 100  FORMAT(' - READ_GC_OMI_SO2_OBS: reading file: ', A)

!      PRINT*, 'N_CALC:', N_CALC

      ! Open file and assign file id (FID)
      CALL CHECK( NF90_OPEN( READ_FILENAME, NF90_NOWRITE, FID ), 0 )

      !-----------------------------------
      ! Get data record IDs
      !-----------------------------------
      CALL CHECK( NF90_INQ_DIMID( FID, "time",    N_ID       ), 100 )
      CALL CHECK( NF90_INQ_VARID( FID, "lon",     LON_ID     ), 101 )
      CALL CHECK( NF90_INQ_VARID( FID, "lat",     LAT_ID     ), 102 )
      CALL CHECK( NF90_INQ_VARID( FID, "time",    TIME_ID    ), 103 )
      CALL CHECK( NF90_INQ_VARID( FID, "SO2_PBL", SO2_PBL_ID ), 104 )
      CALL CHECK( NF90_INQ_VARID( FID, "SO2_TRL", SO2_TRL_ID ), 105 )
      CALL CHECK( NF90_INQ_VARID( FID, "SO2_TRM", SO2_TRM_ID ), 106 )
      CALL CHECK( NF90_INQ_VARID( FID, "ERR_PBL", ERR_PBL_ID ), 108 )
      CALL CHECK( NF90_INQ_VARID( FID, "ERR_TRL", ERR_TRL_ID ), 109 )
      CALL CHECK( NF90_INQ_VARID( FID, "ERR_TRM", ERR_TRM_ID ), 110 )

      IF ( N_CALC > 1 ) THEN

         CALL CHECK( NF90_INQ_VARID( FID, "SO2",  SO2_ID  ), 107 )
         CALL CHECK( NF90_INQ_VARID( FID, "ERR",  ERR_ID  ), 111 )
!         CALL CHECK( NF90_INQ_VARID( FID, "MARK", MARK_ID ), 112 )
         CALL CHECK( NF90_INQ_VARID( FID, "OPT",  OPT_ID  ), 113 )
         ! (ywang, 09/12/14)
         CALL CHECK( NF90_INQ_VARID( FID, "FLAG", FLAG_ID ), 114 )

      END IF

      !------------------------------------
      ! Read dimensions
      !------------------------------------

      ! Read number of observations, N_SO2
      CALL CHECK( NF90_INQUIRE_DIMENSION( FID, N_ID, TMP, N_SO2 ), 200 )

      ! Print to screen
      WRITE(6, 110) N_SO2, GET_NYMD()
 110  FORMAT('      Number of OMI SO2 observations: ' I10, ' in ' I10)

      !-------------------------------------
      ! Read 1D data
      !-------------------------------------

      ! Allocate temporal arrays for 1D data
      ALLOCATE( TMPINT(N_SO2) )
      ALLOCATE( TMP4(N_SO2)   )
      ALLOCATE( TMP8(N_SO2)   )
      TMPINT = 0
      TMP4   = 0E0
      TMP8   = 0D0

      ! Allocate OMI SO2 observations array
      IF ( ALLOCATED( OMI_SO2 ) ) DEALLOCATE( OMI_SO2 )
      ALLOCATE( OMI_SO2(N_SO2) )

      IF ( ALLOCATED( FLAGS) ) DEALLOCATE( FLAGS )
      ALLOCATE( FLAGS(N_SO2) )


      ! Read longitude
      CALL CHECK( NF90_GET_VAR( FID, LON_ID, TMP4 ), 301 )
      OMI_SO2(1:N_SO2)%LON = TMP4(1:N_SO2)

      ! Read latitude
      CALL CHECK( NF90_GET_VAR( FID, LAT_ID, TMP4 ), 302 )
      OMI_SO2(1:N_SO2)%LAT = TMP4(1:N_SO2)

      ! Read time
      CALL CHECK( NF90_GET_VAR( FID, TIME_ID, TMP8 ), 303 )
      OMI_SO2(1:N_SO2)%TIME = TMP8(1:N_SO2)

      ! Read ColumnAmountSO2_PBL
      CALL CHECK( NF90_GET_VAR( FID, SO2_PBL_ID, TMP4 ), 304 )
      OMI_SO2(1:N_SO2)%SO2_PBL = TMP4(1:N_SO2)

      ! Read ColumnAmountSO2_TRL
      CALL CHECK( NF90_GET_VAR( FID, SO2_TRL_ID, TMP4 ), 305 )
      OMI_SO2(1:N_SO2)%SO2_TRL = TMP4(1:N_SO2)

      ! Read ColumnAmountSO2_TRM
      CALL CHECK( NF90_GET_VAR( FID, SO2_TRM_ID, TMP4 ), 306 )
      OMI_SO2(1:N_SO2)%SO2_TRM = TMP4(1:N_SO2)

      ! Read ColumnAmountSO2_PBL error
      CALL CHECK( NF90_GET_VAR( FID, ERR_PBL_ID, TMP4 ), 308 )
      OMI_SO2(1:N_SO2)%ERR_PBL = TMP4(1:N_SO2)

      ! Read ColumnAmountSO2_TRL error
      CALL CHECK( NF90_GET_VAR( FID, ERR_TRL_ID, TMP4 ), 309 )
      OMI_SO2(1:N_SO2)%ERR_TRL = TMP4(1:N_SO2)

      ! Read ColumnAmountSO2_TRM error
      CALL CHECK( NF90_GET_VAR( FID, ERR_TRM_ID, TMP4 ), 310 )
      OMI_SO2(1:N_SO2)%ERR_TRM = TMP4(1:N_SO2)

      IF ( N_CALC > 1 ) THEN

         ! Read SO2 Whose CMA is closest to GEOS-Chem CMA
         CALL CHECK( NF90_GET_VAR( FID, SO2_ID,  TMP4   ), 307 )
         OMI_SO2(1:N_SO2)%SO2  = TMP4(1:N_SO2)

         ! Read SO2 error
         CALL CHECK( NF90_GET_VAR( FID, ERR_ID,  TMP4   ), 311 )
         OMI_SO2(1:N_SO2)%ERR  = TMP4(1:N_SO2)

!         ! Read MARK
!         CALL CHECK( NF90_GET_VAR( FID, MARK_ID, TMPINT ), 312 )
!         OMI_SO2(1:N_SO2)%MARK = TMPINT(1:N_SO2)

         ! Read OPT
         CALL CHECK( NF90_GET_VAR( FID, OPT_ID,  TMPINT ), 313 )
         OMI_SO2(1:N_SO2)%OPT  = TMPINT(1:N_SO2)

         ! Read FLAG
         CALL CHECK( NF90_GET_VAR( FID, FLAG_ID, TMPINT ), 314 )
         OMI_SO2(1:N_SO2)%FLAG = TMPINT(1:N_SO2)

      END IF

      ! Close the file
      CALL CHECK( NF90_CLOSE( FID ), 9999 )

      DEALLOCATE( TMPINT )
      DEALLOCATE( TMP4   )
      DEALLOCATE( TMP8   )

!      IF ( N_CALC > 1 ) THEN
!      WRITE(*,*) '-----------SO2-----------'
!      WRITE(*,*) OMI_SO2(1:N_SO2)%SO2
!      WRITE(*,*) '-----------SO2-----------'
!      WRITE(*,*) '-----------ERR-----------'
!      WRITE(*,*) OMI_SO2(1:N_SO2)%ERR
!      WRITE(*,*) '-----------ERR-----------'
!      WRITE(*,*) '-----------MARK-----------'
!      WRITE(*,*) OMI_SO2(1:N_SO2)%MARK
!      WRITE(*,*) '-----------MARK-----------'
!      WRITE(*,*) '-----------OPT-----------'
!      WRITE(*,*) OMI_SO2(1:N_SO2)%OPT
!      WRITE(*,*) '-----------OPT-----------'
!      END IF

      ! Return to the calling routines
      END SUBROUTINE READ_GC_OMI_SO2_OBS
!
!-----------------------------------------------------------------------------
!
      SUBROUTINE WRITE_GC_OMI_SO2_OBS
!
!*****************************************************************************
!  Subroutine WRITE_GC_OMI_SO2_OBS write observations whose CMA are
!  close to GEOS-Chem CMA
!  (ywang, 07/19/14)
!
!*****************************************************************************
!
      ! Reference to f90 modules
      USE NETCDF
      USE TIME_MOD,          ONLY : EXPAND_DATE, GET_NYMD

      ! Local variables
      INTEGER                :: FID,        N_ID
      INTEGER                :: LON_ID,     LAT_ID,     TIME_ID
      INTEGER                :: SO2_PBL_ID, SO2_TRL_ID, SO2_TRM_ID
      INTEGER                :: SO2_ID
      INTEGER                :: ERR_PBL_ID, ERR_TRL_ID, ERR_TRM_ID
      INTEGER                :: ERR_ID
!      INTEGER                :: MARK_ID
      INTEGER                :: OPT_ID
      INTEGER                :: FLAG_ID
      CHARACTER(LEN=255)     :: DIR
      CHARACTER(LEN=255)     :: WRITE_FILENAME
      INTEGER                :: ICOUNT
      REAL*4,  ALLOCATABLE   :: TMP4(:)
      REAL*8,  ALLOCATABLE   :: TMP8(:)

      !=================================================================
      ! WRITE_GC_OMI_SO2_OBS begins here!
      !=================================================================

      ! Filename root
      WRITE_FILENAME = 'GC_' // TRIM( FILENAME )

      ! Expand date tokens in filename
      CALL EXPAND_DATE( WRITE_FILENAME, GET_NYMD(), 9999 )

      ! Construct complete filename
      DIR            = './data/OMI_SO2/'
      WRITE_FILENAME = TRIM( DIR ) // TRIM( WRITE_FILENAME )

      ! Print to screen
      WRITE(6, 100) TRIM( WRITE_FILENAME )
 100  FORMAT(' - WRITE_GC_OMI_SO2_OBS: reading file: ', A)

      ! Create the netCDF file
      CALL CHECK( NF90_CREATE( WRITE_FILENAME, 0, FID ), 0 )

      ! Define time dimension
      CALL CHECK( NF90_DEF_DIM( FID, 'time', NF90_UNLIMITED, N_ID ),
     &            100 )

      CALL CHECK( NF90_ENDDEF( FID ), 99 )

      ICOUNT = SIZE( OMI_SO2 )

      ALLOCATE( TMP4(ICOUNT)   )
      ALLOCATE( TMP8(ICOUNT)   )

      ! put variable longitude
      TMP4(:) = OMI_SO2(:)%LON
      CALL NCIO_1D( FID, TMP4, 'lon', 'longitude',
     &              'degree', N_ID, ICOUNT, 101 )

      ! put variable latitude
      TMP4(:) = OMI_SO2(:)%LAT
      CALL NCIO_1D( FID, TMP4, 'lat', 'latitude',
     &              'degree', N_ID, ICOUNT, 102 )

      ! put variable time
      TMP8(:) = OMI_SO2(:)%TIME
      CALL NCIO_1D_DBL( FID, TMP8, 'time',
     &                  'time at start of scan (TAI93)',
     &                  's', N_ID, ICOUNT, 103 )

      ! put variable ColumnAmountSO2_PBL
      TMP4(:) = OMI_SO2(:)%SO2_PBL
      CALL NCIO_1D( FID, TMP4, 'SO2_PBL',
     &              'ColumnAmountSO2_PBL',
     &              'DU', N_ID, ICOUNT, 104 )

      ! put variable ColumnAmountSO2_TRL
      TMP4(:) = OMI_SO2(:)%SO2_TRL
      CALL NCIO_1D( FID, TMP4, 'SO2_TRL',
     &              'ColumnAmountSO2_TRL',
     &              'DU', N_ID, ICOUNT, 105 )

      ! put variable ColumnAmountSO2_TRM
      TMP4(:) = OMI_SO2(:)%SO2_TRM
      CALL NCIO_1D( FID, TMP4, 'SO2_TRM',
     &              'ColumnAmountSO2_TRM',
     &              'DU', N_ID, ICOUNT, 106 )

      ! put variable SO2 Whose CMA is closest to GEOS-Chem CMA
      TMP4(:) = OMI_SO2(:)%SO2
      CALL NCIO_1D( FID, TMP4, 'SO2',
     &              'GC_OMI_SO2',
     &              'DU', N_ID, ICOUNT, 107 )

      ! put variable ColumnAmountSO2_PBL error
      TMP4(:) = OMI_SO2(:)%ERR_PBL
      CALL NCIO_1D( FID, TMP4, 'ERR_PBL',
     &              'ColumnAmountSO2_PBL error',
     &              'DU', N_ID, ICOUNT, 108 )

      ! put variable ColumnAmountSO2_TRL error
      TMP4(:) = OMI_SO2(:)%ERR_TRL
      CALL NCIO_1D( FID, TMP4, 'ERR_TRL',
     &              'ColumnAmountSO2_TRL error',
     &              'DU', N_ID, ICOUNT, 109 )

      ! put variable ColumnAmountSO2_TRM error
      TMP4(:) = OMI_SO2(:)%ERR_TRM
      CALL NCIO_1D( FID, TMP4, 'ERR_TRM',
     &              'ColumnAmountSO2_TRM error',
     &              'DU', N_ID, ICOUNT, 110 )

      ! put variable GC_OMI_SO2 error
      TMP4(:) = OMI_SO2(:)%ERR
      CALL NCIO_1D( FID, TMP4, 'ERR',
     &              'GC_OMI_SO2 error',
     &              'DU', N_ID, ICOUNT, 111 )

!      ! put variable MARK
!      CALL NCIO_1D_INT( FID, OMI_SO2(:)%MARK, 'MARK',
!     &                  'MARK',
!     &                  'unitless', N_ID, ICOUNT, 112 )

      ! put variable OPT
      CALL NCIO_1D_INT( FID, OMI_SO2(:)%OPT, 'OPT',
     &                  'OPT',
     &                  'unitless', N_ID, ICOUNT, 113 )

      ! put variable FLAG (ywang, 09/12/14)
      CALL NCIO_1D_INT( FID, OMI_SO2(:)%FLAG, 'FLAG',
     &                  'FLAG',
     &                  'unitless', N_ID, ICOUNT, 114 )

      ! close netCDF file
      CALL CHECK( NF90_CLOSE( FID ), 999 )

      DEALLOCATE( TMP4   )
      DEALLOCATE( TMP8   )

      ! Return to the calling routines
      END SUBROUTINE WRITE_GC_OMI_SO2_OBS
!
!-----------------------------------------------------------------------------
!
      SUBROUTINE CALC_OMI_SO2_FORCE( COST_FUNC )
!
!*****************************************************************************
!  Subroutine CALC_OMI_SO2_FORCE calculate the adjoint forcing from OMI
!  L3 SO2 observation and updates the cost function.
!  (ywang, 07/16/14)
!
!  Arguments as Input/Output:
!  ===========================================================================
!  ( 1) COST_FUNC (REAL*8) : Cost funtion                     [unitless]
!
!  NOTES:
!*****************************************************************************
!
      ! Reference to f90 modules
      USE ADJ_ARRAYS_MOD,     ONLY : STT_ADJ
      USE ADJ_ARRAYS_MOD,     ONLY : OBS_FREQ, N_CALC
      USE CHECKPT_MOD,        ONLY : CHK_STT
      USE DAO_MOD,            ONLY : AIRVOL
      USE DAO_MOD,            ONLY : BXHEIGHT
      USE DIRECTORY_ADJ_MOD,  ONLY : DIAGADJ_DIR
      USE GRID_MOD,           ONLY : GET_IJ
      USE TIME_MOD,           ONLY : GET_NYMD, GET_NHMS
      USE TIME_MOD,           ONLY : GET_TAUb, GET_TAU, GET_TAUe
      USE TRACER_MOD,         ONLY : XNUMOL
      USE TRACERID_MOD,       ONLY : IDTSO2
      USE TROPOPAUSE_MOD,     ONLY : ITS_IN_THE_TROP

      ! Arguements
      REAL*8, INTENT(INOUT)       :: COST_FUNC

      ! Parameters
      REAL*8, PARAMETER           :: M2DU = 2.69D20 ! 1 DU = 2.69D20 molec/m2

      ! Local variables
      INTEGER, SAVE               :: N_SO2 = 0
      INTEGER                     :: N_CURR
      INTEGER                     :: NT
      INTEGER                     :: NC
      INTEGER                     :: IIJJ(2), I, J, L

      REAL*8                      :: OLD_COST
      REAL*8                      :: TMP_COST
      REAL*8, ALLOCATABLE         :: NEW_COST(:)

      REAL*8                      :: GC_CMA              ! GEOS-Chem Center Mass of Alitude. (units: km)
      REAL*8                      :: GC_SO2_CONC(LLPAR)  ! GEOS-Chem SO2 at each layer (units: kg/m3)
      REAL*8                      :: LAYER_GC_SO2(LLPAR) ! GEOS-Chem SO2 at each layer (units: DU)
      REAL*8                      :: COLUMN_GC_SO2       ! GEOS-Chem column SO2 (units: DU)
      REAL*8                      :: DIFF
      REAL*8                      :: FORCING



      !=================================================================
      ! CALC_OMI_SO2_FORCE begins here!
      !=================================================================

      PRINT*, '    - CALC_OMI_SO2_FORCE: OMI SO2 forcing '

      ! Initialize
      CALL INIT_OMI_SO2_OBS

      ! Save a value of the cost function first
      OLD_COST = COST_FUNC
!      PRINT*, 'OLD_COST', OLD_COST, COST_FUNC

      ! Check if it is the last time of a day
      IF ( OBS_FREQ > 60 ) THEN
         PRINT*, '236000 - OBS_FREQ * 100 is not valid'
         STOP
      END IF
      IF ( GET_NHMS() == 236000 - OBS_FREQ * 100 ) THEN
!      IF ( (GET_NHMS() == 236000 - OBS_FREQ * 100) .OR.
!     &     ( ABS(GET_TAU() -GET_TAUe()) < 1e-6  ) ) THEN

         CALL READ_GC_OMI_SO2_OBS( GET_NYMD(), N_SO2 )

      END IF

      ! No observations for current day
      IF ( N_SO2 == 0 ) THEN

         PRINT*, '    - CALC_OMI_SO2_FORCE: No OMI SO2 obsevations for
     &current day'

         IF ( ABS( GET_TAUb() - GET_TAU() ) < 1E-6) THEN

            CALL MAKE_AVERAGE_OMI_SO2

         END IF

         ! The start time of one day
         IF ( GET_NHMS() == 0 ) THEN

            CALL WRITE_GC_OMI_SO2_OBS

         END IF


         RETURN

      END IF


      ! GET observations in time window
      CALL GET_OBS( N_SO2, N_CURR )

      ! No observations for time window
      IF ( N_CURR == 0 ) THEN

         PRINT*, '    - CALC_OMI_SO2_FORCE: No OMI SO2 obsevations for
     &current time window'

         IF ( ABS( GET_TAUb() - GET_TAU() ) < 1E-6) THEN

            CALL MAKE_AVERAGE_OMI_SO2

         END IF

         ! The start time of one day
         IF ( GET_NHMS() == 0 ) THEN

            CALL WRITE_GC_OMI_SO2_OBS

         END IF

         RETURN

      END IF

      ! Reset
      CURR_GC_SO2   = 0D0
      CURR_PBL_SO2  = 0D0
      CURR_TRL_SO2  = 0D0
      CURR_TRM_SO2  = 0D0
      CURR_SO2      = 0D0
      CURR_DIFF_SO2 = 0D0
      CURR_FORCING  = 0D0
      CURR_COST     = 0D0
      CURR_CMA      = 0D0
      CURR_COUNT    = 0
      CURR_PBL_C    = 0
      CURR_TRL_C    = 0
      CURR_TRM_C    = 0
      CURR_OPT     = 0

      ! Reset
      IF ( ALLOCATED( NEW_COST ) ) DEALLOCATE( NEW_COST )
      ALLOCATE ( NEW_COST(N_CURR) )
      NEW_COST = 0D0

      NC = 0
      ! Loop for all observations
      DO NT = 1, N_SO2, 1

         ! Observations in time window and simulation area
         IF ( FLAGS(NT) ) THEN

            ! Get grid box of current record
            IIJJ = GET_IJ( REAL(OMI_SO2(NT)%LON ,4),
     &                     REAL(OMI_SO2(NT)%LAT ,4) )

            I    = IIJJ(1)
            J    = IIJJ(2)

! These code is moved below
!            ! Count observations in gridbox
!            CURR_COUNT(I,J) = CURR_COUNT(I,J) + 1
!            ALL_COUNT(I,J)  = ALL_COUNT(I,J)  + 1

            ! SO2 outside troposphere is set to 0
            GC_SO2_CONC  = 0D0
            LAYER_GC_SO2 = 0D0
            DO L = 1, LLPAR, 1

               IF ( ITS_IN_THE_TROP(I,J,L) ) THEN

                  ! Units conversion [ kg/gridbox => DU ]
                  ! Units of LAYER_GC_SO2: DU
                  !           GC_SO2_CONC: kg/m3
                  !              CHK_STT : kg/gridbox
                  !               AIRVOL : m3/gridbox
                  !             BXHEIGHT : m
                  !               XNUMOL : molec/kg
                  ! M2DU = 2.69D20 , 1 DU = 2.69D20 molec/m2

                  GC_SO2_CONC(L)  = CHK_STT(I,J,L,IDTSO2) /
     &                              AIRVOL(I,J,L)

                  LAYER_GC_SO2(L) = GC_SO2_CONC(L)        *
     &                              BXHEIGHT(I,J,L)       *
     &                              XNUMOL(IDTSO2)        /
     &                              M2DU

               END IF

            END DO

            COLUMN_GC_SO2    = SUM( LAYER_GC_SO2 )

! These code is moved below
!            CURR_GC_SO2(I,J) = CURR_GC_SO2(I,J) + COLUMN_GC_SO2
!            ALL_GC_SO2(I,J)  = ALL_GC_SO2(I,J)  + COLUMN_GC_SO2

            ! Get GEOS-Chem CMA
            CALL CALC_CMA( I, J, GC_SO2_CONC, GC_CMA )

! These code is moved below
!            CURR_CMA(I,J) = CURR_CMA(I,J) + GC_CMA
!            ALL_CMA(I,J)  = ALL_CMA(I,J)  + GC_CMA

            ! Determine SO2 observations whoes CMA is closest to GEOS-Chem CMA
            ! DETERMINE_OBS is modified (ywang, 09/12/14)
            IF ( N_CALC == 1 ) CALL DETERMINE_OBS(GC_CMA, NT)


            IF ( OMI_SO2(NT)%FLAG == 1 ) THEN ! expected observation exist

               ! Count observations in gridbox
               CURR_COUNT(I,J) = CURR_COUNT(I,J) + 1
               ALL_COUNT(I,J)  = ALL_COUNT(I,J)  + 1

               CURR_GC_SO2(I,J) = CURR_GC_SO2(I,J) + COLUMN_GC_SO2
               ALL_GC_SO2(I,J)  = ALL_GC_SO2(I,J)  + COLUMN_GC_SO2

               CURR_CMA(I,J) = CURR_CMA(I,J) + GC_CMA
               ALL_CMA(I,J)  = ALL_CMA(I,J)  + GC_CMA

               CURR_OPT(I,J) = CURR_OPT(I,J) + 1
               ALL_OPT(I,J)  = ALL_OPT(I,J) + 1

               ! Calculate sum here. In the end, average will be calculated
               CURR_SO2(I,J) = CURR_SO2(I,J) + OMI_SO2(NT)%SO2
               ALL_SO2(I,J)  = ALL_SO2(I,J)  + OMI_SO2(NT)%SO2

! The code is moved below
!            IF ( OMI_SO2(NT)%SO2_PBL > 0.0  ) THEN
!
!               CURR_PBL_SO2(I,J) = CURR_PBL_SO2(I,J) +
!     &                             OMI_SO2(NT)%SO2_PBL
!               CURR_PBL_C(I,J)   = CURR_PBL_C(I,J)   + 1
!
!               ALL_PBL_SO2(I,J)  = ALL_PBL_SO2(I,J)  +
!     &                             OMI_SO2(NT)%SO2_PBL
!               ALL_PBL_C(I,J)    = ALL_PBL_C(I,J)    + 1
!
!            END IF
!
!            IF ( OMI_SO2(NT)%SO2_TRL > 0.0  ) THEN
!
!               CURR_TRL_SO2(I,J) = CURR_TRL_SO2(I,J) +
!     &                             OMI_SO2(NT)%SO2_TRL
!               CURR_TRL_C(I,J)   = CURR_TRL_C(I,J)   + 1
!
!               ALL_TRL_SO2(I,J)  = ALL_TRL_SO2(I,J)  +
!     &                             OMI_SO2(NT)%SO2_TRL
!               ALL_TRL_C(I,J)    = ALL_TRL_C(I,J)    + 1
!
!            END IF
!
!            IF ( OMI_SO2(NT)%SO2_TRM > 0.0  ) THEN
!
!               CURR_TRM_SO2(I,J) = CURR_TRM_SO2(I,J) +
!     &                             OMI_SO2(NT)%SO2_TRM
!               CURR_TRM_C(I,J)   = CURR_TRM_C(I,J)   + 1
!
!               ALL_TRM_SO2(I,J)  = ALL_TRM_SO2(I,J)  +
!     &                             OMI_SO2(NT)%SO2_TRM
!               ALL_TRM_C(I,J)    = ALL_TRM_C(I,J)    + 1
!
!            END IF

               !-----------------------------
               ! Calculate adjoint forcing
               !-----------------------------

               ! The difference between GEOS-Chem SO2 and OMI SO2
               DIFF               = COLUMN_GC_SO2      - OMI_SO2(NT)%SO2
               CURR_DIFF_SO2(I,J) = CURR_DIFF_SO2(I,J) + DIFF
               ALL_DIFF_SO2(I,J)  = ALL_DIFF_SO2(I,J)  + DIFF

               ! S_{obs}^{-1} * DIFF
               FORCING           = DIFF / (OMI_SO2(NT)%ERR ** 2)
               CURR_FORCING(I,J) = CURR_FORCING(I,J) + FORCING
               ALL_FORCING(I,J)  = ALL_FORCING(I,J)  + FORCING

               ! Contribution to the cost function
               TMP_COST       = 0.5D0 * DIFF * FORCING
               NC             = NC + 1
               NEW_COST(NC)   = NEW_COST(NC)   + TMP_COST
               CURR_COST(I,J) = CURR_COST(I,J) + TMP_COST
               ALL_COST(I,J)  = ALL_COST(I,J)  + TMP_COST


               ! Now pass the adjoint back to the adjoint tracer array
               DO L = 1, LLPAR, 1

                  IF ( ITS_IN_THE_TROP(I,J,L) ) THEN

                     STT_ADJ(I,J,L,IDTSO2) = STT_ADJ(I,J,L,IDTSO2)   +
     &                                       FORCING / AIRVOL(I,J,L) *
     &                                       BXHEIGHT(I,J,L)         *
     &                                       XNUMOL(IDTSO2) / M2DU

                  END IF

               END DO

            END IF !OMI_SO2(NT)%FLAG == 1

            IF ( OMI_SO2(NT)%SO2_PBL > 0.0  ) THEN

               CURR_PBL_SO2(I,J) = CURR_PBL_SO2(I,J) +
     &                             OMI_SO2(NT)%SO2_PBL
               CURR_PBL_C(I,J)   = CURR_PBL_C(I,J)   + 1

               ALL_PBL_SO2(I,J)  = ALL_PBL_SO2(I,J)  +
     &                             OMI_SO2(NT)%SO2_PBL
               ALL_PBL_C(I,J)    = ALL_PBL_C(I,J)    + 1

            END IF

            IF ( OMI_SO2(NT)%SO2_TRL > 0.0  ) THEN

               CURR_TRL_SO2(I,J) = CURR_TRL_SO2(I,J) +
     &                             OMI_SO2(NT)%SO2_TRL
               CURR_TRL_C(I,J)   = CURR_TRL_C(I,J)   + 1

               ALL_TRL_SO2(I,J)  = ALL_TRL_SO2(I,J)  +
     &                             OMI_SO2(NT)%SO2_TRL
               ALL_TRL_C(I,J)    = ALL_TRL_C(I,J)    + 1

            END IF

            IF ( OMI_SO2(NT)%SO2_TRM > 0.0  ) THEN

               CURR_TRM_SO2(I,J) = CURR_TRM_SO2(I,J) +
     &                             OMI_SO2(NT)%SO2_TRM
               CURR_TRM_C(I,J)   = CURR_TRM_C(I,J)   + 1

               ALL_TRM_SO2(I,J)  = ALL_TRM_SO2(I,J)  +
     &                             OMI_SO2(NT)%SO2_TRM
               ALL_TRM_C(I,J)    = ALL_TRM_C(I,J)    + 1

            END IF

         END IF ! FLAGS(NT)

      END DO ! NT


      ! Update cost function
      COST_FUNC = COST_FUNC + SUM( NEW_COST )
      PRINT*, ' Update value of COST_FUNC = ', COST_FUNC
      PRINT*, ' OMI SO2 contribution      = ', COST_FUNC - OLD_COST

      PRINT*, ' MIN/MAX STT_ADJ  = ', MINVAL(STT_ADJ), MAXVAL(STT_ADJ)
      PRINT*, ' MIN/MAX in       = ', MINLOC(STT_ADJ), MAXLOC(STT_ADJ)
      PRINT*, ' MIN/MAX NEW_COST = ', MINVAL(NEW_COST), MAXVAL(NEW_COST)
      PRINT*, ' MIN/MAX cost in  = ', MINLOC(NEW_COST),MAXLOC(NEW_COST)


      CALL MAKE_CURRENT_OMI_SO2


      IF ( ABS( GET_TAUb() - GET_TAU() ) < 1E-6) THEN

         CALL MAKE_AVERAGE_OMI_SO2

      END IF

      ! The start time of one day
      IF ( ( GET_NHMS() == 0 ) .AND. ( N_CALC == 1 ) ) THEN

         CALL WRITE_GC_OMI_SO2_OBS

      END IF

      ! Return to the calling routines
      END SUBROUTINE CALC_OMI_SO2_FORCE
!
!-----------------------------------------------------------------------------
!
      SUBROUTINE MAKE_CURRENT_OMI_SO2
!
!*****************************************************************************
!  Subroutine MAKE_CURRENT_OMI_SO2 output some dignostic data for
!  current assimilation time window
!  (ywang, 07/19/14)
!
!*****************************************************************************
!
      ! Reference to f90 module
      USE ADJ_ARRAYS_MOD,    ONLY : N_CALC
      USE ADJ_ARRAYS_MOD,    ONLY : EXPAND_NAME
      USE BPCH2_MOD
      USE DIRECTORY_ADJ_MOD, ONLY : DIAGADJ_DIR
      USE FILE_MOD,          ONLY : IU_RST,      IOERROR
      USE GRID_MOD,          ONLY : GET_XOFFSET, GET_YOFFSET
      USE TIME_MOD,          ONLY : EXPAND_DATE
      USE TIME_MOD,          ONLY : GET_NYMD, GET_NHMS
      USE TIME_MOD,          ONLY : GET_TAU

      REAL*4, PARAMETER    :: UNDEF = -999.0

      ! Local Variables
      INTEGER              :: I,    I0, J,  J0, L
      REAL*4               :: SO2(IIPAR,JJPAR,14)
      CHARACTER(LEN=255)   :: FILENAME

      ! For binary punch file, version 2.0
      REAL*4               :: LONRES, LATRES
      INTEGER, PARAMETER   :: HALFPOLAR = 1
      INTEGER, PARAMETER   :: CENTER180 = 1

      CHARACTER(LEN=255)   :: OUTPUT_FILE
      CHARACTER(LEN=20)    :: MODELNAME
      CHARACTER(LEN=40)    :: CATEGORY
      CHARACTER(LEN=40)    :: UNIT
      CHARACTER(LEN=40)    :: RESERVED = ''
      CHARACTER(LEN=80)    :: TITLE

      !=================================================================
      ! MAKE_CURRENT_OMI_SO2 begins here!
      !=================================================================

      ! Hardwire output file for now
      OUTPUT_FILE = 'gctm.omi.so2.YYYYMMDD.hhmm.NN'

      ! Define variables for BINARY PUNCH FILE OUTPUT
      TITLE    = 'GEOS-CHEM diag File: ' //
     &           'OMI SO2'
      UNIT     = 'DU'
      CATEGORY = 'IJ-AVG-$'
      LONRES   = DISIZE
      LATRES   = DJSIZE

      ! Call GET_MODELNAME to return the proper model name for
      ! the given met data being used
      MODELNAME = GET_MODELNAME()

      ! Get the nested-grid offsets
      I0 = GET_XOFFSET( GLOBAL=.TRUE. )
      J0 = GET_YOFFSET( GLOBAL=.TRUE. )

      ! Copy the output observation file name into a local variable
      FILENAME = TRIM( OUTPUT_FILE )

      ! Replace YYMMDD and hhmmss token w/ actucl value
      CALL EXPAND_DATE( FILENAME, GET_NYMD(), GET_NHMS() )

      ! Replace NN token w/ actual value
      CALL EXPAND_NAME( FILENAME, N_CALC )

      ! Add DIAGADJ_DIR prefix to FILENAME
      FILENAME = TRIM( DIAGADJ_DIR ) // TRIM( FILENAME )

      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - MAKE_CURRENT_OMI_SO2:  Writing ', a )

      ! Open file for output
      CALL OPEN_BPCH2_FOR_WRITE( IU_RST, FILENAME, TITLE )

      SO2 = 0.0

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         IF ( CURR_COUNT(I,J) > 0 ) THEN
            SO2(I,J, 1) = REAL(CURR_GC_SO2(I,J))   / CURR_COUNT(I,J)
            SO2(I,J, 5) = REAL(CURR_SO2(I,J))      / CURR_COUNT(I,J)
            SO2(I,J, 6) = REAL(CURR_DIFF_SO2(I,J)) / CURR_COUNT(I,J)
            SO2(I,J, 7) = REAL(CURR_FORCING(I,J))
            SO2(I,J, 8) = REAL(CURR_COST(I,J))
            SO2(I,J, 9) = REAL(CURR_CMA(I,J))      / CURR_COUNT(I,J)
            SO2(I,J,10) = REAL(CURR_COUNT(I,J))
            SO2(I,J,11) = REAL(CURR_PBL_C(I,J))
            SO2(I,J,12) = REAL(CURR_TRL_C(I,J))
            SO2(I,J,13) = REAL(CURR_TRM_C(I,J))
            SO2(I,J,14) = REAL(CURR_OPT(I,J))
         END IF

         IF ( CURR_PBL_C(I,J) > 0 ) THEN
            SO2(I,J,2) = REAL(CURR_PBL_SO2(I,J)) / CURR_PBL_C(I,J)
         END IF

         IF ( CURR_TRL_C(I,J) > 0 ) THEN
            SO2(I,J,3) = REAL(CURR_TRL_SO2(I,J)) / CURR_TRL_C(I,J)
         END IF

         IF ( CURR_TRM_C(I,J) > 0 ) THEN
            SO2(I,J,4) = REAL(CURR_TRM_SO2(I,J)) / CURR_TRM_C(I,J)
         END IF

      ENDDO
      ENDDO

!$OMP END PARALLEL DO

      CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,
     &            HALFPOLAR, CENTER180, CATEGORY,  26,
     &            UNIT,      GET_TAU(), GET_TAU(), RESERVED,
     &            IIPAR,     JJPAR,     14,        I0+1,
     &            J0+1,      1,         SO2 )

      ! Close file
      CLOSE( IU_RST )

      ! Return to the calling routines
      END SUBROUTINE MAKE_CURRENT_OMI_SO2
!
!-----------------------------------------------------------------------------
!
      SUBROUTINE MAKE_AVERAGE_OMI_SO2
!
!*****************************************************************************
!  Subroutine MAKE_AVERAGE_OMI_SO2 output some dignostic data for
!  simulation time
!  (ywang, 07/19/14)
!
!*****************************************************************************
!
      ! Reference to f90 module
      USE ADJ_ARRAYS_MOD,    ONLY : N_CALC
      USE ADJ_ARRAYS_MOD,    ONLY : EXPAND_NAME
      USE BPCH2_MOD
      USE DIRECTORY_ADJ_MOD, ONLY : DIAGADJ_DIR
      USE FILE_MOD,          ONLY : IU_RST,      IOERROR
      USE GRID_MOD,          ONLY : GET_XOFFSET, GET_YOFFSET
      USE TIME_MOD,          ONLY : EXPAND_DATE
      USE TIME_MOD,          ONLY : GET_NYMD
      USE TIME_MOD,          ONLY : GET_TAUb, GET_TAUe

      REAL*4, PARAMETER    :: UNDEF = -999.0

      ! Local Variables
      INTEGER              :: I,    I0, J,  J0, L
      REAL*4               :: SO2(IIPAR,JJPAR,14)
      CHARACTER(LEN=255)   :: FILENAME

      ! For binary punch file, version 2.0
      REAL*4               :: LONRES, LATRES
      INTEGER, PARAMETER   :: HALFPOLAR = 1
      INTEGER, PARAMETER   :: CENTER180 = 1

      CHARACTER(LEN=255)   :: OUTPUT_FILE
      CHARACTER(LEN=20)    :: MODELNAME
      CHARACTER(LEN=40)    :: CATEGORY
      CHARACTER(LEN=40)    :: UNIT
      CHARACTER(LEN=40)    :: RESERVED = ''
      CHARACTER(LEN=80)    :: TITLE

      !=================================================================
      ! MAKE_AVERAGE_OMI_SO2 begins here!
      !=================================================================

      ! Hardwire output file for now
      OUTPUT_FILE = 'gctm.omi.so2.ave.YYYYMMDD.NN'

      ! Define variables for BINARY PUNCH FILE OUTPUT
      TITLE    = 'GEOS-CHEM diag File: ' //
     &           'OMI SO2'
      UNIT     = 'DU'
      CATEGORY = 'IJ-AVG-$'
      LONRES   = DISIZE
      LATRES   = DJSIZE

      ! Call GET_MODELNAME to return the proper model name for
      ! the given met data being used
      MODELNAME = GET_MODELNAME()

      ! Get the nested-grid offsets
      I0 = GET_XOFFSET( GLOBAL=.TRUE. )
      J0 = GET_YOFFSET( GLOBAL=.TRUE. )

      ! Copy the output observation file name into a local variable
      FILENAME = TRIM( OUTPUT_FILE )

      ! Replace YYMMDD token w/ actucl value
      CALL EXPAND_DATE( FILENAME, GET_NYMD(), 9999 )

      ! Replace NN token w/ actual value
      CALL EXPAND_NAME( FILENAME, N_CALC )

      ! Add DIAGADJ_DIR prefix to FILENAME
      FILENAME = TRIM( DIAGADJ_DIR ) // TRIM( FILENAME )

      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - MAKE_AVERAGE_OMI_SO2:  Writing ', a )

      ! Open file for output
      CALL OPEN_BPCH2_FOR_WRITE( IU_RST, FILENAME, TITLE )

      SO2 = 0.0

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         IF ( ALL_COUNT(I,J) > 0 ) THEN
            SO2(I,J, 1) = REAL(ALL_GC_SO2(I,J))   / ALL_COUNT(I,J)
            SO2(I,J, 5) = REAL(ALL_SO2(I,J))      / ALL_COUNT(I,J)
            SO2(I,J, 6) = REAL(ALL_DIFF_SO2(I,J)) / ALL_COUNT(I,J)
            SO2(I,J, 7) = REAL(ALL_FORCING(I,J))
            SO2(I,J, 8) = REAL(ALL_COST(I,J))
            SO2(I,J, 9) = REAL(ALL_CMA(I,J))      / ALL_COUNT(I,J)
            SO2(I,J,10) = REAL(ALL_COUNT(I,J))
            SO2(I,J,11) = REAL(ALL_PBL_C(I,J))
            SO2(I,J,12) = REAL(ALL_TRL_C(I,J))
            SO2(I,J,13) = REAL(ALL_TRM_C(I,J))
            SO2(I,J,14) = REAL(ALL_OPT(I,J))
         END IF

         IF ( ALL_PBL_C(I,J) > 0 ) THEN
            SO2(I,J,2) = REAL(ALL_PBL_SO2(I,J)) / ALL_PBL_C(I,J)
         END IF

         IF ( ALL_TRL_C(I,J) > 0 ) THEN
            SO2(I,J,3) = REAL(ALL_TRL_SO2(I,J)) / ALL_TRL_C(I,J)
         END IF

         IF ( ALL_TRM_C(I,J) > 0 ) THEN
            SO2(I,J,4) = REAL(ALL_TRM_SO2(I,J)) / ALL_TRM_C(I,J)
         END IF

      ENDDO
      ENDDO

!$OMP END PARALLEL DO

      CALL BPCH2( IU_RST,    MODELNAME,  LONRES,     LATRES,
     &            HALFPOLAR, CENTER180,  CATEGORY,   26,
     &            UNIT,      GET_TAUb(), GET_TAUe(), RESERVED,
     &            IIPAR,     JJPAR,      14,         I0+1,
     &            J0+1,      1,          SO2 )

      ! Close file
      CLOSE( IU_RST )

      ! Return to the calling routines
      END SUBROUTINE MAKE_AVERAGE_OMI_SO2
!
!-----------------------------------------------------------------------------
!
      SUBROUTINE GET_OBS( N_SO2, N_CURR )
!
!*****************************************************************************
!  Subroutine GET_OBS finds all obsevations in the current time window
!  and simulation area.
!  (ywang, 07/17/14)
!
!  Arguements as Input:
!  ===========================================================================
!  ( 1) N_SO2     (INTEGER) : Number of observation in current day
!  Arguements as Output:
!  ===========================================================================
!  ( 2) N_CURR    (INTEGER) : Number of obsevations in the current time
!  window and simulation area.
!
!  Module variable as Output:
!  ===========================================================================
!  ( 1) FLAGS     (LOGICAL) : Whether or not a speicific obsevation is
!  in current time window and simulation area.
!
!*****************************************************************************
!
      ! Reference to f90 module
      USE GRID_MOD,      ONLY : GET_XEDGE, GET_YEDGE
      USE TIME_MOD,      ONLY : GET_JD,   GET_TAU
      USE TIME_MOD,      ONLY : GET_NYMD, GET_NHMS

      ! Arguemens
      INTEGER, INTENT( IN)   :: N_SO2
      INTEGER, INTENT(OUT)   :: N_CURR

      ! Local variables
      REAL*8                 :: HALF_TIME_WINDOW
      REAL*8                 :: WINDOW_BEGIN,    WINDOW_END
      REAL*8                 :: JD85,            JD93,      JD93_85
      REAL*8                 :: CURRENT_TAU
      INTEGER                :: NT

#if defined( NESTED_CH ) || defined( NESTED_NA ) || defined( NESTED_SD )
      REAL*8, SAVE           :: XEDGE_MIN, XEDGE_MAX
      REAL*8, SAVE           :: YEDGE_MIN, YEDGE_MAX
#endif

      !=================================================================
      ! GET_OBS begins here!
      !=================================================================

      ! Get the difference between JD93 and JD85
      ! In GEOS-Chem, it is since 1/1/1985, while in OMI, it is since
      ! 1/1/1993
      JD85    = GET_JD( 19850000, 000000 )
      JD93    = GET_JD( 19930000, 000000 )
      JD93_85 = JD93 - JD85
      JD93_85 = JD93_85 * 24D0 ! days => hours

!      write(*,  '(A15, f30.15)' ) 'JD85', JD85
!      write(*, '(A15, f30.15)' ) 'JD93', JD93
!      write(*, '(A15, f30.15)' ) 'JD93_85', JD93_85

      ! Get current GEOS-Chem TAU
      CURRENT_TAU = GET_TAU()
!      write(*,'(A15, f30.15)') 'CURRENT_TAU',CURRENT_TAU

      ! Change current GEOS-Chem TAU into OMI TAU
      CURRENT_TAU = CURRENT_TAU - JD93_85
!      write(*,'(A15, f30.15)') 'CURRENT_TAU',CURRENT_TAU

      ! Change TAU units ( hours => second )
      CURRENT_TAU = CURRENT_TAU * 3600D0
!      write(*,'(A15, f30.15)') 'CURRENT_TAU',CURRENT_TAU

      ! Get half time window
      HALF_TIME_WINDOW = TIME_WINDOW / 2D0
      HALF_TIME_WINDOW = HALF_TIME_WINDOW * 60D0 ! ( minute => second )
!      write(*, '(A15, f30.15)') 'HALF_TIME_WINDOW',HALF_TIME_WINDOW

      ! Get current time window
      WINDOW_BEGIN = CURRENT_TAU - HALF_TIME_WINDOW
      WINDOW_END   = CURRENT_TAU + HALF_TIME_WINDOW

!      write(*, '(A15,f30.15)') 'WINDOW_BEGIN',WINDOW_BEGIN
!      write(*, '(A15,f30.15)') 'WINDOW_END',WINDOW_END

#if defined( NESTED_CH ) || defined( NESTED_NA ) || defined( NESTED_SD )
         XEDGE_MIN = GET_XEDGE( 1 )
         XEDGE_MAX = GET_XEDGE( IIPAR+1 )
         YEDGE_MIN = GET_YEDGE( 1 )
         YEDGE_MAX = GET_YEDGE( JJPAR+1 )
         PRINT*, 'Nested region edge limit'
         PRINT*, 'XEDGE_MIN: ', XEDGE_MIN
         PRINT*, 'XEDGE_MAX: ', XEDGE_MAX
         PRINT*, 'YEDGE_MIN: ', YEDGE_MIN
         PRINT*, 'YEDGE_MAX: ', YEDGE_MAX
#endif
!      Write(*, '(f30.15)') OMI_SO2(1:10)%TIME
       N_CURR = 0
      ! Find observations in current time window and simulation area
      DO NT = 1, N_SO2, 1

         IF (      ( OMI_SO2(NT)%TIME >= WINDOW_BEGIN )
     &       .AND. ( OMI_SO2(NT)%TIME <  WINDOW_END   )
#if defined( NESTED_CH ) || defined( NESTED_NA ) || defined( NESTED_SD )
     &       .AND. ( OMI_SO2(NT)%LON  >= XEDGE_MIN    )
     &       .AND. ( OMI_SO2(NT)%LON  <= XEDGE_MAX    )
     &       .AND. ( OMI_SO2(NT)%LAT  >= YEDGE_MIN    )
     &       .AND. ( OMI_SO2(NT)%LAT  <= YEDGE_MAX    )
#endif
     &                                                  ) THEN


            FLAGS(NT) = .TRUE.

            N_CURR    = N_CURR + 1

         ELSE

            FLAGS(NT) = .FALSE.

         END IF

      END DO

      WRITE(6, 100) N_CURR, GET_NHMS()
 100  FORMAT('      Number of OMI SO2 observations: ' I10 ' at ' I10.6)

      ! Return to calling program
      END SUBROUTINE GET_OBS
!
!-----------------------------------------------------------------------------
!
      SUBROUTINE CALC_CMA( I, J, GC_SO2_CONC, GC_CMA )
!
!*****************************************************************************
!  Subroutine CALC_CMA calculate GEOS-Chem Center Mass of Altitude (
!  The height where SO2 reach its maximum value )
!  (ywang, 07/17/14)
!
!  Arguements as Input:
!  ===========================================================================
!  ( 1) I                  (INTEGER) : LON INDEX
!  ( 2) J                  (INTEGER) : LAT INDEX
!  ( 3) GC_SO2_CONC(LLPAR)  (REAL*8) : SO2 concentration at each layer
!  Arguements as Output:
!  ===========================================================================
!  ( 1) GC_CMA              (REAL*8) : Center Mass of Altitude (km)
!
!*****************************************************************************
!
      ! Reference to f90 modules
      USE DAO_MOD,            ONLY : BXHEIGHT

      ! Arguements
      INTEGER, INTENT( IN)        :: I, J
      REAL*8,  INTENT( IN)        :: GC_SO2_CONC(LLPAR)
      REAL*8,  INTENT(OUT)        :: GC_CMA

      ! Local variabls
      INTEGER                     :: L, MAX_SO2_L

      !=================================================================
      ! CALC_CMA begins here!
      !=================================================================

      MAX_SO2_L = 1

      ! Find the layer SO2 concentration reaches its maximum value
      DO L = 2, LLPAR, 1

         IF ( GC_SO2_CONC(L) > GC_SO2_CONC(MAX_SO2_L) ) THEN

            MAX_SO2_L = L

         END IF

      END DO
!      MAX_SO2_L = MAXLOC( GC_SO2_CONC )

      ! Calculate CMA
      GC_CMA = 0D0
      IF ( MAX_SO2_L == 1) THEN

         GC_CMA = BXHEIGHT(I,J,MAX_SO2_L) / 2D0

      ELSE

         DO L = 1, MAX_SO2_L-1, 1

            GC_CMA = GC_CMA + BXHEIGHT(I,J,L)

         END DO

         GC_CMA = GC_CMA + BXHEIGHT(I,J,MAX_SO2_L) / 2D0

      END IF

      ! m => km
      GC_CMA = GC_CMA / 1000D0

      ! Return to calling program
      END SUBROUTINE CALC_CMA
!
!-----------------------------------------------------------------------------
!
      SUBROUTINE DETERMINE_OBS( GC_CMA, NT )
!
!*****************************************************************************
!  Subroutine DETERMINE_OBS finds the OMI SO2 observation whose CMA is
!  closest to GEOS-Chem CMA
!  (ywang, 07/18/14)
!
!  Arguements as Input:
!  ===========================================================================
!  ( 1) GC_CMA           (REAL*8 ) : GEOS-Chem Center Mass of Altitude (km)
!  ( 2) NT               (INTEGER) : OMI_SO2 index
!                                    1: PBL, 2: TRL, 3: TRM
!
!  Moudule variables as Output:
!  ( 1) OMI_SO2(NT)%SO2   (REAL*8) : OMI SO2 will be used
!  ( 2) OMI_SO2(NT)%ERR   (REAL*8) : OMI SO2 error will be used
!!  ( 3) OMI_SO2(NT)%MARK (INTEGER) : Marks of which observation is
!       expected
!  ( 4) OMI_SO2(NT)%OPT  (INTEGER) : Mark of which observation whose CMA
!       is close to GEOS-Chem CMA, though the observation may not exists
!  ( 5) OMI_SO2(NT)%FLAG (INTEGER) : Does the observartion expected
!       really exist? 0: not exist, 1: exist
!
!  Notes:
!  ( 1) The alogrithm is changed. Now if the expected obsevation does
!  not exit, we shall not use other observation to substutite it.
!
!*****************************************************************************
!
      ! Parameters
!      ! PBL_CMA, TRL_CMA, and TRM_CMA references to moudle note 2
!      REAL*8, PARAMETER      :: PBL_CMA    = 0.9D0 ! units: km
!      REAL*8, PARAMETER      :: TRL_CMA    = 2.5D0 ! units: km
!      REAL*8, PARAMETER      :: TRM_CMA    = 7.5D0 ! units: km
      REAL*8, PARAMETER      :: MIN_PBL = 0.0D0   ! units: km
      REAL*8, PARAMETER      :: MAX_PBL = 1.7D0   ! units: km
      REAL*8, PARAMETER      :: MIN_TRL = MAX_PBL
      REAL*8, PARAMETER      :: MAX_TRL = 5.0D0   ! units: km
      REAL*8, PARAMETER      :: MIN_TRM = MAX_TRL
      REAL*8, PARAMETER      :: MAX_TRM = 10.0D0  ! units: km

      ! minimum valid observations
      REAL*8, PARAMETER      :: MIN_VAL = -10.0   ! units: DU

      REAL*8, PARAMETER      :: UNDEF = -999.0D0

      ! Arguments
      REAL*8,  INTENT( IN)   :: GC_CMA
      INTEGER, INTENT( IN)   :: NT

!      ! Local variables
!      REAL*8                 :: DIST1(3), DIST2(3)
!      INTEGER                :: I(1)

      !=================================================================
      ! DETERMINE_OBS begins here!
      !=================================================================

!      ! for OPT
!      DIST1(1) = ABS( GC_CMA-PBL_CMA )
!      DIST1(2) = ABS( GC_CMA-TRL_CMA )
!      DIST1(3) = ABS( GC_CMA-TRM_CMA )

!      I = MINLOC( DIST1 )
!
!      SELECT CASE ( I(1) )
!
!         CASE ( 1 )
!            OMI_SO2(NT)%SO2 = OMI_SO2(NT)%SO2_PBL
!            OMI_SO2(NT)%ERR = OMI_SO2(NT)%ERR_PBL
!
!         CASE ( 2 )
!            OMI_SO2(NT)%SO2 = OMI_SO2(NT)%SO2_TRL
!            OMI_SO2(NT)%ERR = OMI_SO2(NT)%ERR_TRL
!         CASE ( 3 )
!            OMI_SO2(NT)%SO2 = OMI_SO2(NT)%SO2_TRM
!            OMI_SO2(NT)%ERR = OMI_SO2(NT)%ERR_TRM
!
!         CASE DEFAULT
!            PRINT*, 'DETERMINE_OBS ERROR'
!            STOP
!
!      END SELECT
!
!      OMI_SO2(NT)%OPT = I(1)

!      ! for MARK
!      DIST2 = 1D10
!      IF ( OMI_SO2(NT)%SO2_PBL > 0.0 ) DIST2(1) = ABS( GC_CMA-PBL_CMA )
!      IF ( OMI_SO2(NT)%SO2_TRL > 0.0 ) DIST2(2) = ABS( GC_CMA-TRL_CMA )
!      IF ( OMI_SO2(NT)%SO2_TRM > 0.0 ) DIST2(3) = ABS( GC_CMA-TRM_CMA )
!
!      I = MINLOC( DIST2 )
!
!      SELECT CASE ( I(1) )
!
!         CASE ( 1 )
!            OMI_SO2(NT)%SO2 = OMI_SO2(NT)%SO2_PBL
!            OMI_SO2(NT)%ERR = OMI_SO2(NT)%ERR_PBL
!
!         CASE ( 2 )
!            OMI_SO2(NT)%SO2 = OMI_SO2(NT)%SO2_TRL
!            OMI_SO2(NT)%ERR = OMI_SO2(NT)%ERR_TRL
!         CASE ( 3 )
!            OMI_SO2(NT)%SO2 = OMI_SO2(NT)%SO2_TRM
!            OMI_SO2(NT)%ERR = OMI_SO2(NT)%ERR_TRM
!
!         CASE DEFAULT
!            PRINT*, 'DETERMINE_OBS ERROR'
!            STOP
!
!      END SELECT
!
!      OMI_SO2(NT)%MARK = I(1)

      OMI_SO2(NT)%SO2 = UNDEF
      OMI_SO2(NT)%ERR = UNDEF

      IF ( (GC_CMA >= MIN_PBL) .AND. (GC_CMA < MAX_PBL) ) THEN

         OMI_SO2(NT)%OPT = 1

         IF ( OMI_SO2(NT)%SO2_PBL >= MIN_VAL ) THEN
            OMI_SO2(NT)%SO2 = OMI_SO2(NT)%SO2_PBL
            OMI_SO2(NT)%ERR = OMI_SO2(NT)%ERR_PBL
         END IF

      ELSE IF ( (GC_CMA >= MIN_TRL) .AND. (GC_CMA < MAX_TRL) ) THEN

         OMI_SO2(NT)%OPT = 2

         IF ( OMI_SO2(NT)%SO2_TRL >= MIN_VAL ) THEN
            OMI_SO2(NT)%SO2 = OMI_SO2(NT)%SO2_TRL
            OMI_SO2(NT)%ERR = OMI_SO2(NT)%ERR_TRL
         END IF

      ELSE IF ( (GC_CMA >= MIN_TRM) .AND. (GC_CMA < MAX_TRM) ) THEN

         OMI_SO2(NT)%OPT = 3

         IF ( OMI_SO2(NT)%SO2_TRM >= MIN_VAL ) THEN
            OMI_SO2(NT)%SO2 = OMI_SO2(NT)%SO2_TRM
            OMI_SO2(NT)%ERR = OMI_SO2(NT)%ERR_TRM
         END IF

      ELSE

         OMI_SO2(NT)%OPT = 0

      END IF

      IF ( OMI_SO2(NT)%SO2 > MIN_VAL ) THEN

         OMI_SO2(NT)%FLAG = 1

      ELSE

         OMI_SO2(NT)%FLAG = 0

      END IF


      ! Return to calling program
      END SUBROUTINE DETERMINE_OBS
!
!-----------------------------------------------------------------------------
!
      SUBROUTINE CHECK( STATUS, LOCATION )
!
!*****************************************************************************
!  Subroutine CHECK checks the status of calls to netCDF libraries
!  routines
!  (dkh, 02/15/09)
!
!  Arguments as Input:
!  ===========================================================================
!  (1 ) STATUS    (INTEGER) : Completion status of netCDF library call
!  (2 ) LOCATION  (INTEGER) : Location at which netCDF library call was
!  made
!
!  NOTES:
!
!*****************************************************************************
!
      ! Reference to f90 modules
      USE ERROR_MOD,    ONLY  : ERROR_STOP
      USE NETCDF

      ! Arguments
      INTEGER, INTENT(IN)    :: STATUS
      INTEGER, INTENT(IN)    :: LOCATION

      !=================================================================
      ! CHECK begins here!
      !=================================================================

      IF ( STATUS /= NF90_NOERR ) THEN
        WRITE(6,*) TRIM( NF90_STRERROR( STATUS ) )
        WRITE(6,*) 'At location = ', LOCATION
        CALL ERROR_STOP('netCDF error', 'omi_so2_mod')
      ENDIF

      ! Return to calling program
      END SUBROUTINE CHECK
!
!-----------------------------------------------------------------------------
!
      SUBROUTINE NCIO_1D (NCID, VAR1D, VARNAME, LONGNAME, VARUNIT,
     &                    DIMID, DIMV, INDEX)

      ! References to F90 modules
      USE NETCDF

      ! Arguments
      INTEGER, INTENT(IN)          :: NCID, DIMID, DIMV, INDEX
      REAL*4,  INTENT(IN)          :: VAR1D(DIMV)
      CHARACTER(LEN=*), INTENT(IN) :: VARNAME, LONGNAME, VARUNIT

      ! Local variable
      INTEGER                     :: DIMS(1)
      INTEGER                     :: VAR_ID

      CALL CHECK( NF90_REDEF(NCID), INDEX )

      DIMS(1) = DIMID
      CALL CHECK( NF90_DEF_VAR(NCID, TRIM(VARNAME),
     &            NF90_FLOAT, DIMS, VAR_ID), INDEX)
      CALL CHECK( NF90_PUT_ATT(NCID,VAR_ID,"name",TRIM(LONGNAME)),INDEX)
      CALL CHECK( NF90_PUT_ATT(NCID,VAR_ID,"units",TRIM(VARUNIT)),INDEX)
      CALL CHECK( NF90_ENDDEF(NCID), INDEX )
      CALL CHECK( NF90_PUT_VAR(NCID, VAR_ID, VAR1D ), INDEX )

      END SUBROUTINE NCIO_1D
!
!-----------------------------------------------------------------------------
!
      SUBROUTINE NCIO_1D_DBL (NCID, VAR1D, VARNAME, LONGNAME, VARUNIT,
     &                        DIMID, DIMV, INDEX)

      ! References to F90 modules
      USE NETCDF

      ! Arguments
      INTEGER, INTENT(IN)          :: NCID, DIMID, DIMV, INDEX
      REAL*8,  INTENT(IN)          :: VAR1D(DIMV)
      CHARACTER(LEN=*), INTENT(IN) :: VARNAME, LONGNAME, VARUNIT

      ! Local variable
      INTEGER                     :: DIMS(1)
      INTEGER                     :: VAR_ID

      CALL CHECK( NF90_REDEF(NCID), INDEX )

      DIMS(1) = DIMID
      CALL CHECK( NF90_DEF_VAR(NCID, TRIM(VARNAME),
     &            NF90_DOUBLE, DIMS, VAR_ID), INDEX)
      CALL CHECK( NF90_PUT_ATT(NCID,VAR_ID,"name",TRIM(LONGNAME)),INDEX)
      CALL CHECK( NF90_PUT_ATT(NCID,VAR_ID,"units",TRIM(VARUNIT)),INDEX)
      CALL CHECK( NF90_ENDDEF(NCID), INDEX )
      CALL CHECK( NF90_PUT_VAR(NCID, VAR_ID, VAR1D ), INDEX )

      END SUBROUTINE NCIO_1D_DBL
!
!-----------------------------------------------------------------------------
!
      SUBROUTINE NCIO_1D_INT (NCID, VAR1D, VARNAME, LONGNAME, VARUNIT,
     &                        DIMID, DIMV, INDEX)

      ! References to F90 modules
      USE NETCDF

      ! Arguments
      INTEGER, INTENT(IN)          :: NCID, DIMID, DIMV, INDEX
      INTEGER, INTENT(IN)          :: VAR1D(DIMV)
      CHARACTER(LEN=*), INTENT(IN) :: VARNAME, LONGNAME, VARUNIT

      ! Local variable
      INTEGER                     :: DIMS(1)
      INTEGER                     :: VAR_ID

      CALL CHECK( NF90_REDEF(NCID), INDEX )

      DIMS(1) = DIMID
      CALL CHECK( NF90_DEF_VAR(NCID, TRIM(VARNAME),
     &            NF90_INT, DIMS, VAR_ID), INDEX)
      CALL CHECK( NF90_PUT_ATT(NCID,VAR_ID,"name",TRIM(LONGNAME)),INDEX)
      CALL CHECK( NF90_PUT_ATT(NCID,VAR_ID,"units",TRIM(VARUNIT)),INDEX)
      CALL CHECK( NF90_ENDDEF(NCID), INDEX )
      CALL CHECK( NF90_PUT_VAR(NCID, VAR_ID, VAR1D ), INDEX )

      END SUBROUTINE NCIO_1D_INT
!
!-----------------------------------------------------------------------------
!
      END MODULE OMI_SO2_OBS_MOD
