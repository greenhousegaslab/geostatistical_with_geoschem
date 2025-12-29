      MODULE MOPITT_OBS_MOD
!
!*****************************************************************************
!  Module MOPITT_OBS_MOD contains all the subroutines for the using of MOPITT
!  observation (version 3 and version 4).(zhe 1/19/11)
!  Remove the support to MOPITT v3 and v4. Now support v5 and v6. (Zhe 1/20/14)
!  Module Routines:
!  ============================================================================
!  (1 ) READ_MOPITT_FILE       : Read MOPITT hdf file
!  (2 ) CALC_MOPITT_FORCE      : Calculates cost function and STT_ADJ increments
!  (3 ) CALC_AVGKER            : Construct the averging kernel matrix
!  (4 ) BIN_DATA               : Interpolation between different vertical resolutions
!  (5 ) INIT_DOMAIN            : Define the observation window
!  (6 ) CALC_OBS_HOUR          : Calculated hour of morning obs
!  (7 ) ITS_TIME_FOR_MOPITT_OBS: FUNCTION that checks time vs. OBS_HOUR array
!  (8 ) READ_MOP02             : Reads MOPITT data fields from the HDF-EOS file
!  (9 ) APRIORI_MOP02          : Read A priori field for MOPITT version 3
!  (10) INFO_MOP02             : Prints name, dims, type, etc. of MOPITT data fields
!  (11) CLEANUP_MOP02          : Deallocates all module arrays
!  =============================================================================

      IMPLICIT NONE

#     include "CMN_SIZE"
#     include "../adjoint/define_adj.h"

      PRIVATE

      PUBLIC OBS_HOUR_MOPITT
      PUBLIC COUNT_TOTAL
      PUBLIC ITS_TIME_FOR_MOPITT_OBS
      PUBLIC READ_MOPITT_FILE
      PUBLIC CALC_MOPITT_FORCE
      PUBLIC MAKE_MOPITT_BIAS_FILE_HDF5

      !=============================================================================
      ! MODULE VARIABLES
      !=============================================================================

      INTEGER   :: OBS_HOUR_MOPITT(IIPAR,JJPAR)
      INTEGER   :: DOMAIN_OBS(IIPAR,JJPAR)
      REAL*8    :: COUNT_TOTAL

      REAL*4    :: ERR_PERCENT(IIPAR,JJPAR)
      REAL*4,  ALLOCATABLE :: A(:,:)
      REAL*4,  ALLOCATABLE :: T(:)
      REAL*4,  ALLOCATABLE :: XA(:)
      REAL*8,  ALLOCATABLE :: AC(:)

      ! MOPITT dimension fields
      INTEGER   :: T_DIM, Z_DIM
      REAL*4, ALLOCATABLE :: LATITUDE(:)
      REAL*4, ALLOCATABLE :: LONGITUDE(:)
      REAL*4, ALLOCATABLE :: PRESSURE(:)
      REAL*4, ALLOCATABLE :: SECONDS_IN_DAY(:)
      REAL*4, ALLOCATABLE :: MOPITT_GMT(:)
      REAL*8, ALLOCATABLE :: TAU(:)

      ! MOPITT data quantities
      REAL*4, ALLOCATABLE :: BOTTOM_PRESSURE(:)
      REAL*4, ALLOCATABLE :: CO_MIXING_RATIO(:,:,:)
      REAL*4, ALLOCATABLE :: CO_RET_BOT_MIXING_RATIO(:,:)
      REAL*4, ALLOCATABLE :: CO_TOTAL_COLUMN(:,:)
      REAL*4, ALLOCATABLE :: COVARIANCE(:,:,:)
      REAL*4, ALLOCATABLE :: AVGKER(:,:,:)
      INTEGER, ALLOCATABLE :: CLOUD_DES(:)
      INTEGER, ALLOCATABLE :: SURFACE_INDEX(:)

      ! MOPITT a priori
      INTEGER :: NLEV_AP
      REAL*4, ALLOCATABLE :: PLEV_AP(:)
      REAL*4, ALLOCATABLE :: CH4_MR_AP(:)
      REAL*4, ALLOCATABLE :: CO_MR_AP(:,:,:)
      REAL*4, ALLOCATABLE :: CO_MR_AP_BOTTOM(:,:)
      REAL*4, ALLOCATABLE :: COV_CO_AP(:,:)

      ! mkeller:
      REAL*4 :: MOPITT_CO_MEAN(IIPAR,JJPAR,10) = 0d0
      REAL*4 :: MOPITT_GEOS_CO_MEAN(IIPAR,JJPAR,10) = 0d0
      REAL*4 :: MOPITT_BIAS(IIPAR,JJPAR,10) = 0d0
      REAL*4 :: MOPITT_BIAS_COUNT(IIPAR,JJPAR,10) = 0d0
      REAL*4 :: MOPITT_CHI_SQUARED(IIPAR,JJPAR,10) = 0d0

      REAL*4 :: MOPITT_GC_CO_SOBS(IIPAR,JJPAR,10) = 0d0
      REAL*4 :: MOPITT_MOP_CO_SOBS(IIPAR,JJPAR,10) = 0d0
      REAL*4 :: MOPITT_BIAS_SOBS(IIPAR,JJPAR,10) = 0d0
      REAL*4 :: MOPITT_COUNT_SOBS(IIPAR,JJPAR,10) = 0d0

      REAL*4 :: MOPITT_GC_CO_TMP(IIPAR,JJPAR,10) = 0d0
      REAL*4 :: MOPITT_MOP_CO_TMP(IIPAR,JJPAR,10) = 0d0
      REAL*4 :: MOPITT_BIAS_TMP(IIPAR,JJPAR,10) = 0d0
      REAL*4 :: MOPITT_COUNT_TMP(IIPAR,JJPAR,10) = 0d0

      INTEGER :: MAXLEV = 10
      LOGICAL :: FIRST_FLEX = .TRUE.
      
      !mkeller: arrays for saving diagnostics

      TYPE FLEX_REAL_1D            
      INTEGER :: CURRENT_N, MAX_N 
      REAL*8,ALLOCATABLE :: DATA(:) 
      ENDTYPE FLEX_REAL_1D

      TYPE FLEX_REAL_2D            
      INTEGER :: CURRENT_N, MAX_N 
      REAL*8,ALLOCATABLE :: DATA(:,:) 
      ENDTYPE FLEX_REAL_2D

      TYPE(FLEX_REAL_1D) :: FLEX_LON, FLEX_LAT, FLEX_TIME
      TYPE(FLEX_REAL_2D) :: FLEX_MOP_CO, FLEX_GC_CO 


      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE READ_MOPITT_FILE( YYYYMMDD, HHMMSS )

!******************************************************************************
!  Subroutine READ_MOPITT_FILE reads the MOPITT hdf file.
!  (mak, 7/12/07, zhe 1/19/11)
!  Arguments as input:
!  ============================================================================
!  (1 ) YYYYMMDD : Year-Month-Day
!  (2 ) HHMMSS   :  and Hour-Min-Sec for which to read restart file
!******************************************************************************

      USE ERROR_MOD, ONLY : ALLOC_ERR
      USE TIME_MOD,  ONLY : EXPAND_DATE, GET_MONTH, GET_YEAR

#     include "CMN_SIZE"   ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN) :: YYYYMMDD, HHMMSS

      ! Local Variables
      INTEGER             :: I, IOS, J, L, N, AS
      CHARACTER(LEN=255)  :: DIR_MOPITT
      CHARACTER(LEN=255)  :: DIR_MONTH
      CHARACTER(LEN=255)  :: FILENAMEM
      CHARACTER(LEN=255)  :: FILENAME2
      LOGICAL             :: IT_EXISTS
      LOGICAL, SAVE       :: FIRST = .TRUE.

      !=================================================================
      ! READ_MOPITT_FILE begins here!
      !=================================================================
#if defined( MOPITT_V5_CO_OBS )
      DIR_MOPITT = '/met/gc/MOPITT_CO/'
      DIR_MONTH = 'v5/YYYY/MM/'
      FILENAMEM = 'MOP02J-YYYYMMDD-L2V10.1.3.beta.hdf'
#endif
#if defined( MOPITT_V6_CO_OBS )
      DIR_MOPITT = '/met/gc/MOPITT_CO/'
      DIR_MONTH = 'v6/YYYY/MM/'
      FILENAMEM = 'MOP02J-YYYYMMDD-L2V16.2.3.he5'
#endif

      IF ( FIRST ) THEN
         ERR_PERCENT(:,:) = 0.0
         COUNT_TOTAL = 0
         FIRST            = .FALSE.
      ENDIF

      OBS_HOUR_MOPITT(:,:) = -99

      CALL EXPAND_DATE( FILENAMEM, YYYYMMDD, 0 )
      CALL EXPAND_DATE( DIR_MONTH, YYYYMMDD, 0 )

      FILENAME2 = TRIM( DIR_MOPITT ) // TRIM( DIR_MONTH ) // FILENAMEM
      PRINT*, '=== Reading ===:', TRIM( FILENAME2 )

      INQUIRE( FILE = FILENAME2, EXIST = IT_EXISTS )
      IF (IT_EXISTS) THEN

         !CALL INFO_MOP02(FILENAME2)

         CALL READ_MOP02( FILENAME2 )

         CALL INIT_DOMAIN

         ! Calculate hour of day when obs should be compared to model
         CALL CALC_OBS_HOUR

      ENDIF
      !CALL READ_ERROR_VARIANCE

      END SUBROUTINE READ_MOPITT_FILE
!-------------------------------------------------------------------------------------------------

      SUBROUTINE CALC_MOPITT_FORCE

!******************************************************************************
! CALC_MOPITT_FORCE calculate cost function and STT_ADJ increments
! Now support v5 and v6. (zhe 1/20/14)
!******************************************************************************

      USE PRESSURE_MOD, ONLY : GET_PCENTER, GET_AP, GET_BP
      USE BPCH2_MOD,    ONLY : GET_TAU0
      USE TIME_MOD,     ONLY : GET_MONTH, GET_DAY, GET_YEAR
      USE TIME_MOD,     ONLY : GET_HOUR
      USE CHECKPT_MOD,  ONLY : CHK_STT
      USE TRACER_MOD,   ONLY : TCVV
      USE TRACERID_MOD, ONLY : IDTCO
      USE DAO_MOD,      ONLY : AD, IS_LAND
      USE ADJ_ARRAYS_MOD,   ONLY : SET_FORCING, SET_MOP_MOD_DIFF,
     &                  SET_MODEL_BIAS, SET_MODEL, SET_OBS,
     &                  COST_ARRAY, DAY_OF_SIM, IFD, JFD, LFD, NFD,
     &                  COST_FUNC, ADJ_FORCE, STT_ADJ
      USE LOGICAL_ADJ_MOD,  ONLY : LPRINTFD, LDCOSAT
      USE ERROR_MOD,        ONLY : IT_IS_NAN, ERROR_STOP
      USE GRID_MOD,         ONLY : GET_IJ

      ! Local Variables
      INTEGER ::  W, I, J, Z, ZZ, L,LL
      INTEGER ::  LON15, IIJJ(2)
      INTEGER ::  NLEV_RET
      INTEGER ::  SKIP_LEVELS

      REAL*4  ::  RETLEV(Z_DIM+1)
      REAL*8  ::  P_EDGE(Z_DIM+2), MODEL_COL, MOPITT_COL
      REAL*8  ::  UTC, TAU0
      REAL*8  ::  MODEL_P(LLPAR), MODEL_CO_MR(LLPAR)
      REAL*8  ::  COUNT_GRID(IIPAR,JJPAR)
      REAL*8  ::  COUNT(IIPAR,JJPAR)
      REAL*8  ::  MOP_COL_GRID(IIPAR,JJPAR)
      REAL*8  ::  MODEL_COL_GRID(IIPAR,JJPAR)
      REAL*8  ::  NEW_COST(IIPAR,JJPAR)
      REAL*8  ::  ADJ_F(LLPAR)
      REAL*8  ::  MODEL_P_EDGE(LLPAR+1)

      REAL*8, ALLOCATABLE :: GEOS_RAW(:)
      REAL*8, ALLOCATABLE :: MOP_CO(:)
      REAL*8, ALLOCATABLE :: DIFF_ADJ(:)
      REAL*8, ALLOCATABLE :: GEOS_CO(:)
      REAL*8, ALLOCATABLE :: DIFF_COST(:)
      REAL*8, ALLOCATABLE :: FORCING(:)
      REAL*4, ALLOCATABLE :: INV_S(:,:)

      !mkeller: initialize NLEV_RET to zero
      !         necessary to avoid avoid problems when aggregating super observations
      !         in the absense of data
      NLEV_RET = 0d0

      IF( FIRST_FLEX ) THEN
         CALL INIT_FLEX_REAL_1D(FLEX_LON)
         CALL INIT_FLEX_REAL_1D(FLEX_LAT)
         CALL INIT_FLEX_REAL_1D(FLEX_TIME)
         CALL INIT_FLEX_REAL_2D(FLEX_MOP_CO)
         CALL INIT_FLEX_REAL_2D(FLEX_GC_CO)
         FIRST_FLEX = .FALSE.
      ENDIF

      !=================================================================
      ! CALC_MOPITT_FORCE begins here!
      !=================================================================

      TAU0 = GET_TAU0( GET_MONTH(), GET_DAY(), GET_YEAR() )

      COUNT_GRID(:,:)     = 0d0
      COUNT(:,:)          = 0d0
      MOP_COL_GRID(:,:)   = -999.0
      MODEL_COL_GRID(:,:) = -999.0
      ADJ_FORCE(:,:,:,:)  = 0d0
      NEW_COST(:,:)       = 0d0

      !=================================================================
      ! Loop over MOPITT data
      !=================================================================
      DO W = 1, T_DIM

         ! Compute local time:
         ! Local TIME = GMT + ( LONGITUDE / 15 ) since each hour of time
         ! corresponds to 15 degrees of LONGITUDE on the globe
         LON15 = LONGITUDE(W) / 15.
         UTC   = TAU(W) - TAU0 + LON15
         IF (UTC < 0. )  UTC = UTC + 24
         IF (UTC > 24.)  UTC = UTC - 24


         !Only consider day time MOPITT measurements
         ! am = 12 hrs centered on 10:30am local time (so 4:30am-4:30pm)
#if defined( GRID05x0666 ) && defined( NESTED_CH ) && !defined( NESTED_SD )
         IF ( UTC >= 4.5 .and. UTC <= 16.5
     &    .and. LONGITUDE(W) > 70
     &    .and. LONGITUDE(W) < 150
     &    .and. LATITUDE(W)  > -11
     &    .and. LATITUDE(W)  < 55 ) THEN
#elif defined( GRID05x0666 ) && defined( NESTED_NA ) && !defined( NESTED_SD )
         IF ( UTC >= 4.5 .and. UTC <= 16.5
     &    .and. LONGITUDE(W) > -140
     &    .and. LONGITUDE(W) < -40
     &    .and. LATITUDE(W)  > 10
     &    .and. LATITUDE(W)  < 70 ) THEN
#elif defined( GRID05x0666 ) && (defined( NESTED_NA ) || defined( NESTED_CH )) && defined( NESTED_SD )
         IF ( UTC >= 4.5 .and. UTC <= 16.5
     &    .and. LONGITUDE(W) > -126
     &    .and. LONGITUDE(W) < -66
     &    .and. LATITUDE(W)  > 13
     &    .and. LATITUDE(W)  < 57 ) THEN
#else
         IF ( UTC >= 4.5 .and. UTC <= 16.5 ) THEN
#endif

         ! Get grid box
         IIJJ  = GET_IJ( LONGITUDE(W), LATITUDE(W))
         I = IIJJ(1)
         J = IIJJ(2)

         !=================================================================
         ! Data selection
         !=================================================================
         IF( GET_HOUR() == OBS_HOUR_MOPITT(I,J) .and.
     &      CLOUD_DES(W) == 2.0 .and.
     &      CO_TOTAL_COLUMN(1,W) > 5E17 .and.
     &      DOMAIN_OBS(I,J) == 1 ) THEN

!         IF ( (IS_LAND(I,J) .AND.
!     &         LATITUDE(W) .GE. -52 .AND. LATITUDE(W) .LE. 52 ) .OR.   !52S-52N
!     &        (LATITUDE(W) .GE. -40 .AND. LATITUDE(W) .LE. 40) ) THEN  !40S-40N

         IF (LATITUDE(W) .GE. -60 .AND. LATITUDE(W) .LE. 60) THEN  !40S-40N

            RETLEV(:) = -999.0
            MODEL_COL = 0D0
            MOPITT_COL = 0D0

            ! Create pressure profile
            RETLEV(1) = BOTTOM_PRESSURE(W)
            ZZ = 0
            ! Loop over Mopitt levels
            DO Z = 1, Z_DIM
            ! Always start from the bottom pressure,
            ! even if it means skipping a MOPITT pressure level
               IF ( PRESSURE(Z) >= RETLEV(1) ) THEN
                  ZZ = ZZ + 1
                  CYCLE
               ENDIF
               ! Save into profile
               RETLEV(Z+1-ZZ) = PRESSURE(Z)
            ENDDO

            NLEV_RET = Z_DIM+1 - ZZ

            DO L = 1, NLEV_RET
               P_EDGE(L) = RETLEV(L)
            ENDDO
            P_EDGE(NLEV_RET+1) = 36


            ALLOCATE( XA( NLEV_RET ) )
            ALLOCATE( T( NLEV_RET ) )
            ALLOCATE( A( NLEV_RET,NLEV_RET ) )
            ALLOCATE( AC( NLEV_RET ) )
            ALLOCATE( MOP_CO( NLEV_RET ) )
            ALLOCATE( GEOS_RAW( NLEV_RET ) )
            ALLOCATE( DIFF_ADJ( NLEV_RET ) )
            ALLOCATE( GEOS_CO( NLEV_RET ) )
            ALLOCATE( DIFF_COST( NLEV_RET ) )
            ALLOCATE( FORCING( NLEV_RET ) )
            ALLOCATE( INV_S( NLEV_RET,NLEV_RET ) )


            ! MOPITT CO vertical profile
            MOP_CO(1) = CO_RET_BOT_MIXING_RATIO(1,W)
            MOP_CO(2:NLEV_RET) = CO_MIXING_RATIO(1,11-NLEV_RET:9,W)
            MOP_CO = MOP_CO * 1E-9

            ! COMPUTE AVERAGING KERNEL
            CALL CALC_AVGKER(NLEV_RET, W, RETLEV, MOP_CO)

            ! COMPUTE INVERSE OF CONVARIANCE
            CALL CALL_INV_S(NLEV_RET, W, INV_S, 0.0)

            !USE MOPITT SURFACE PRESSURE
            !DO L=1, LLPAR + 1
            !   MODEL_P_EDGE(L) = GET_AP(L) + GET_BP(L) * RETLEV(1)
            !ENDDO

            DO L = 1, LLPAR
               !MOPITT PRESSURE LEVEL
               !MODEL_P(L) = (MODEL_P_EDGE(L) + MODEL_P_EDGE(L+1)) / 2

               ! Get GC pressure levels (mbar)
               MODEL_P(L) = GET_PCENTER(I,J,L)

               ! Obtain archieved forward model results
               ! kg -> v/v
               MODEL_CO_MR(L) = CHK_STT(I,J,L,IDTCO) *
     &                            TCVV(IDTCO) / AD(I,J,L)
            ENDDO

            ! Interplote the model to MOPITT vertical grids
            CALL BIN_DATA(MODEL_P, P_EDGE, MODEL_CO_MR(:),
     &            GEOS_RAW, NLEV_RET, 1)

            !=================================================================
            ! Apply MOPITT observation operator
            !=================================================================

            ! Total Column: C = T * XA  + AC * ( Xm - XA )
            ! Stratosphere Levels are removed
            DO L = 1, NLEV_RET - 1
               MODEL_COL = MODEL_COL
     &                   + T(L) * XA(L)
     &                   + AC(L) * (LOG10(GEOS_RAW(L))
     &                   - LOG10(XA(L)))
               MOPITT_COL = MOPITT_COL + T(L) * MOP_CO(L)
            ENDDO

            GEOS_CO(:) = 0d0
            ! Smoothed Profile: X_hat = XA  + A * ( Xm - XA )
            DO L = 1, NLEV_RET
               DO LL = 1, NLEV_RET
                  GEOS_CO(L) = GEOS_CO(L)
     &                       + A(L,LL)
     &                       * (LOG10( GEOS_RAW(LL) ) - LOG10( XA(LL) ))
               ENDDO
               GEOS_CO(L) = LOG10( XA(L) ) + GEOS_CO(L)
            ENDDO

            !=================================================================
            ! COST FUNCTION
            !=================================================================
            DIFF_COST(:) = 0D0
            DO L = 1, NLEV_RET - 1
               DIFF_COST(L)  = GEOS_CO(L) - LOG10( MOP_CO(L) )
            ENDDO

            FORCING(:) = 0d0
            DO L = 1, NLEV_RET - 1

               DO LL = 1, NLEV_RET - 1
                  FORCING(L) = FORCING(L) + INV_S(L,LL) * DIFF_COST(LL)
               ENDDO

               NEW_COST(I,J) = NEW_COST(I,J)
     &                       + 0.5d0 * DIFF_COST(L) * FORCING(L)
               COUNT(I,J) = COUNT(I,J) + 1

               !mkeller: populate a number of arrays for diagnostics
               !         only use super observations for now
!               MOPITT_CO_MEAN(I,J,L) = MOPITT_CO_MEAN(I,J,L) +
!     &              MOP_CO(L)
!               MOPITT_GEOS_CO_MEAN(I,J,L) = MOPITT_GEOS_CO_MEAN(I,J,L) +
!     &              10**(GEOS_CO(L))
!               MOPITT_BIAS(I,J,L) = MOPITT_BIAS(I,J,L) +
!     &              10**(GEOS_CO(L)) - MOP_CO(L)
!               MOPITT_CHI_SQUARED(I,J,L) = MOPITT_CHI_SQUARED(I,J,L) +
!     &              ( (GEOS_CO(L) - LOG10(MOP_CO(L)))**2 )/SY
!               MOPITT_BIAS_COUNT(I,J,L) = MOPITT_BIAS_COUNT(I,J,L) + 1d0
!
               MOPITT_MOP_CO_TMP(I,J,L) = MOPITT_MOP_CO_TMP(I,J,L) +
     &              MOP_CO(L)
               MOPITT_GC_CO_TMP(I,J,L) = MOPITT_GC_CO_TMP(I,J,L) +
     &              10**(GEOS_CO(L))
               MOPITT_BIAS_TMP(I,J,L) = MOPITT_BIAS_TMP(I,J,L) +
     &              10**(GEOS_CO(L)) - MOP_CO(L)
               MOPITT_COUNT_TMP(I,J,L) = MOPITT_COUNT_TMP(I,J,L) + 1d0

            ENDDO

            CALL PUSH_FLEX_REAL_1D(FLEX_LON, REAL(LONGITUDE(W),8))
            CALL PUSH_FLEX_REAL_1D(FLEX_LAT, REAL(LATITUDE(W),8))
            CALL PUSH_FLEX_REAL_1D(FLEX_TIME, TAU(W))
      
            CALL PUSH_FLEX_REAL_2D(FLEX_MOP_CO, MOP_CO, 10)
            CALL PUSH_FLEX_REAL_2D(FLEX_GC_CO, 10**(GEOS_CO),10)
            

            !=================================================================
            ! adjoint operator
            !=================================================================
            DIFF_ADJ(:) = 0D0
            DO L = 1, NLEV_RET
               DO LL = 1, NLEV_RET
                  DIFF_ADJ(L) = DIFF_ADJ(L) + A(LL,L) * FORCING(LL)
               ENDDO
               ! fwd code: LOG(GEOS_RAW) - LOG(XA)
               DIFF_ADJ(L) = DIFF_ADJ(L) / (GEOS_RAW(L) * LOG(10.0))
            ENDDO

            CALL BIN_DATA( MODEL_P,  P_EDGE, ADJ_F,
     &                        DIFF_ADJ, NLEV_RET, -1   )

            ! adjoint FORCE
            DO L = 1, LLPAR

               !v/v->kg
               ADJ_FORCE(I,J,L,IDTCO) = ADJ_FORCE(I,J,L,IDTCO)
     &                            + ADJ_F(L) * TCVV(IDTCO)
     &                            / AD(I,J,L)

            ENDDO

            COUNT_GRID(I,J)     = COUNT_GRID(I,J) + 1.d0
            MOP_COL_GRID(I,J)   = MOP_COL_GRID(I,J) + MOPITT_COL
            MODEL_COL_GRID(I,J) = MODEL_COL_GRID(I,J) + MODEL_COL

            IF ( ALLOCATED( GEOS_RAW ) ) DEALLOCATE( GEOS_RAW )
            IF ( ALLOCATED( MOP_CO   ) ) DEALLOCATE( MOP_CO   )
            IF ( ALLOCATED( DIFF_ADJ ) ) DEALLOCATE( DIFF_ADJ )
            IF ( ALLOCATED( A        ) ) DEALLOCATE( A        )
            IF ( ALLOCATED( AC       ) ) DEALLOCATE( AC       )
            IF ( ALLOCATED( T        ) ) DEALLOCATE( T        )
            IF ( ALLOCATED( XA       ) ) DEALLOCATE( XA       )
            IF ( ALLOCATED( GEOS_CO  ) ) DEALLOCATE( GEOS_CO  )
            IF ( ALLOCATED( DIFF_COST) ) DEALLOCATE( DIFF_COST)
            IF ( ALLOCATED( FORCING  ) ) DEALLOCATE( FORCING  )
            IF ( ALLOCATED( INV_S    ) ) DEALLOCATE( INV_S    )

         ENDIF !IS_LAND

         ENDIF !OBS_HOUR

         ENDIF !local time



      ENDDO  !loop over MOPITT data

      !=================================================================
      ! BIN OUTPUT INFO INTO MODEL GRID BOXES
      !=================================================================
      DO I = 1, IIPAR
         DO J = 1, JJPAR

            IF ( COUNT_GRID(I,J) > 0d0 ) THEN

               !The mean value in the grid
               MOP_COL_GRID(I,J)   = MOP_COL_GRID(I,J)
     &                             / COUNT_GRID(I,J)
               MODEL_COL_GRID(I,J) = MODEL_COL_GRID(I,J)
     &                             / COUNT_GRID(I,J)
               ADJ_FORCE(I,J,:,IDTCO)  = ADJ_FORCE(I,J,:,IDTCO)
     &                             / COUNT_GRID(I,J)
               NEW_COST(I,J)       = NEW_COST(I,J) / COUNT_GRID(I,J)
               COUNT(I,J)          = COUNT(I,J) / COUNT_GRID(I,J)

               !Update adjoint tracer
               STT_ADJ(I,J,:,IDTCO) = STT_ADJ(I,J,:,IDTCO) +
     &                                   ADJ_FORCE(I,J,:,IDTCO)

               ! Diagnostic stuff: FORCING, MOP_MOD_DIFF, MODEL_BIAS
               IF( LDCOSAT )THEN

                  CALL SET_FORCING( I, J, DAY_OF_SIM,
     &                              ADJ_FORCE(I,J,1,IDTCO) )
                  CALL SET_MOP_MOD_DIFF( I, J, DAY_OF_SIM,
     &               MODEL_COL_GRID(I,J) - MOP_COL_GRID(I,J) )

                  CALL SET_MODEL_BIAS( I, J, DAY_OF_SIM, 1,
     &              ( MODEL_COL_GRID(I,J) - MOP_COL_GRID(I,J) ) /
     &                                 MOP_COL_GRID(I,J)           )
                  CALL SET_MODEL     ( I, J, DAY_OF_SIM, 1,
     &                                 MODEL_COL_GRID(I,J)         )
                  CALL SET_OBS       ( I, J, DAY_OF_SIM, 1,
     &                                 MOP_COL_GRID(I,J)           )

                  COST_ARRAY(I,J,DAY_OF_SIM) =
     &               COST_ARRAY(I,J,DAY_OF_SIM) + NEW_COST(I,J)

               ENDIF

               IF ( IT_IS_NAN( NEW_COST(I,J) ) ) THEN
                  PRINT*, 'I=', I, 'J=', J
                  CALL ERROR_STOP( 'NEW_COST is NaN',
     &                             'CALC_MOPITT_FORCE')
               ENDIF

            ENDIF !COUNT_GRID

            !mkeller: populate super observation arrays

            IF(NLEV_RET > 0d0) THEN

            DO L=1,NLEV_RET

               IF(MOPITT_COUNT_TMP(I,J,L) > 0d0) THEN

                  MOPITT_GC_CO_SOBS(I,J,L) =
     &                 MOPITT_GC_CO_SOBS(I,J,L) +
     &                 MOPITT_GC_CO_TMP(I,J,L)/MOPITT_COUNT_TMP(I,J,L)

                  MOPITT_MOP_CO_SOBS(I,J,L) =
     &                 MOPITT_MOP_CO_SOBS(I,J,L) +
     &                 MOPITT_MOP_CO_TMP(I,J,L)/MOPITT_COUNT_TMP(I,J,L)

                  MOPITT_BIAS_SOBS(I,J,L) =
     &                 MOPITT_BIAS_SOBS(I,J,L) +
     &                 MOPITT_BIAS_TMP(I,J,L)/MOPITT_COUNT_TMP(I,J,L)

                  MOPITT_COUNT_SOBS(I,J,L) =
     &                 MOPITT_COUNT_SOBS(I,J,L) + 1d0

               ENDIF

            ENDDO

            ENDIF

         ENDDO
      ENDDO      
      
      ! mkeller: reset temporary arrays

      MOPITT_MOP_CO_TMP = 0d0
      MOPITT_GC_CO_TMP = 0d0
      MOPITT_BIAS_TMP = 0d0
      MOPITT_COUNT_TMP = 0d0

      IF (LPRINTFD)  THEN
         PRINT*, 'IFD, JFD= ', IFD, JFD
         PRINT*, 'MODEL_STT:', MODEL_COL_GRID(IFD,JFD)
         PRINT*, 'OBS_STT:', MOP_COL_GRID(IFD,JFD)
         PRINT*, 'NEW_COST', NEW_COST(IFD,JFD)
         PRINT*, 'ADJ_FORCE:', ADJ_FORCE(IFD,JFD,:,IDTCO)
         PRINT*, 'STT_ADJ:', STT_ADJ(IFD,JFD,:,IDTCO)
      ENDIF

      ! Update cost function
      PRINT*, 'TOTAL NEW_COST = ', SUM(NEW_COST)
      PRINT*, 'COST_FUNC BEFORE ADDING NEW_COST=', COST_FUNC
      COST_FUNC   = COST_FUNC   + SUM ( NEW_COST )
      COUNT_TOTAL = COUNT_TOTAL + SUM ( COUNT    )
      PRINT*, 'Total observation number:', COUNT_TOTAL

      ! Return to calling program
      END SUBROUTINE CALC_MOPITT_FORCE
!--------------------------------------------------------------------------------------------

      SUBROUTINE CALC_AVGKER( NLEV_RET, W, RETLEV, MOP_CO )

!******************************************************************************
! SUBROUTINE CALC_AVGKER construct the averging kernel matrix
! (zhe 1/19/11)
!******************************************************************************

      INTEGER :: ILEV, JLEV, ILEV2, JLEV2, Z, W
      INTEGER :: NLEV_RET
      REAL*4  :: DELP(NLEV_RET)
      REAL*4  :: RETLEV(NLEV_RET)
      REAL*8  :: MOP_CO(NLEV_RET)

      REAL*8 :: AVGKER_RET(NLEV_RET, NLEV_RET)
      REAL*8, PARAMETER :: log10e = LOG10(2.71828183)

      !=================================================================
      ! CALC_AVGKER begins here!
      !=================================================================

      A(:,:) = 0d0
      AC(:)  = 0d0

      XA(1) = CO_MR_AP_BOTTOM(1, W)
      XA(2:NLEV_RET) = CO_MR_AP(1,11-NLEV_RET:9,W)
      XA = XA * 1E-9

      !Remove bad levels from averging kernel matrix
      IF ( NLEV_RET < 10 ) THEN
         DO ILEV = 1, NLEV_RET
            ILEV2 = ILEV + ( 10 - NLEV_RET )
            DO JLEV =1, NLEV_RET
               JLEV2 = JLEV + ( 10 - NLEV_RET)
               A(ILEV,JLEV) =
     &              AVGKER(ILEV2,JLEV2,W)
            ENDDO
         ENDDO
      ELSE
         A(:,:) = AVGKER(:,:,W)
      ENDIF

      DELP(1) = RETLEV(1) - RETLEV(2)
      DELP(2:NLEV_RET-1) = 100D0
      DELP(NLEV_RET) = 74D0

      ! transfer function [v/v -> molec/cm2]
      T = 2.12E+22 * DELP

      ! Convert to column averaging kernel
      DO JLEV = 1, NLEV_RET
         DO ILEV = 1, NLEV_RET
            AC(JLEV) = AC(JLEV) + DELP(ILEV) * MOP_CO(ILEV)
     &               * A(ILEV,JLEV)
         ENDDO
         AC(JLEV) = (2.12E+22 / log10e ) * AC(JLEV)
      ENDDO

      END SUBROUTINE CALC_AVGKER
!------------------------------------------------------------------------------------

      SUBROUTINE BIN_DATA( P_MODEL,  P_EDGE, DATA_MODEL, DATA_MOP,
     &                        NLEV_RET, FB                            )

!******************************************************************************
!Based on the code from Monika.  (zhe 1/19/11)
!FB = 1 for forward
!FB = -1 for adjoint
!******************************************************************************

      INTEGER :: L, LL, FB
      INTEGER :: NLEV_RET, NB
      REAL*8  :: P_MODEL(LLPAR)
      REAL*8  :: DATA_MODEL(LLPAR), DATA_MOP(NLEV_RET), DATA_TEM
      REAL*8  :: P_EDGE(NLEV_RET+1)

      !=================================================================
      ! BIN_DATA begins here!
      !=================================================================

      IF (FB > 0) THEN

         DO L = 1, NLEV_RET
            DO LL = 1, LLPAR
               IF ( P_MODEL(LL) <= P_EDGE(L) ) THEN
                  DATA_MOP(L) = DATA_MODEL(LL)
                  EXIT
               ENDIF
            ENDDO
         ENDDO

         DO L = 1, NLEV_RET
            NB = 0
            DATA_TEM = 0
            DO LL = 1, LLPAR
               IF ( ( P_MODEL(LL) <= P_EDGE(L)) .and.
     &              ( P_MODEL(LL) > P_EDGE(L+1)) ) THEN
                  DATA_TEM = DATA_TEM + DATA_MODEL(LL)
                  NB = NB + 1
               ENDIF
            ENDDO
            IF (NB > 0) DATA_MOP(L) = DATA_TEM / NB
         ENDDO

      ELSE

         DATA_MODEL(:) = 0.
         DO L = 1, LLPAR
            DO LL = 1, NLEV_RET
               IF ( ( P_MODEL(L) <= P_EDGE(LL)) .and.
     &              ( P_MODEL(L) > P_EDGE(LL+1)) ) THEN
                  DATA_MODEL(L) = DATA_MOP(LL)
               ENDIF
            ENDDO
         ENDDO

      ENDIF


      ! Return to calling program
      END SUBROUTINE BIN_DATA
!-----------------------------------------------------------------------------------

      SUBROUTINE INIT_DOMAIN

!******************************************************************************
!Define the observatio region
!******************************************************************************
#     include "CMN_SIZE"   ! Size parameters

      !local variables
      INTEGER :: I, J

      !=================================================================
      ! INIT_DOMAIN begins here!
      !=================================================================

      DOMAIN_OBS(:,:) = 0d0

      DO J = 1, JJPAR
      DO I = 1, IIPAR

#if   defined( GRID05x0666 )
!     The surrounding region is used as cushion
!     (zhe 11/28/10)
         IF ( J >= 8 .and. J <= JJPAR-7 .and.
     &        I >= 7 .and. I <= IIPAR-6
#elif defined( GRID2x25 )
         IF ( J >= 16 .and. J <= 76   !60S-60N
#elif defined( GRID4x5 )
         IF ( J >= 9 .and. J <= 39    !60S-60N
#endif
     &       ) DOMAIN_OBS(I,J) = 1d0

      ENDDO
      ENDDO

      PRINT*, sum(DOMAIN_obs), 'MAX observations today'

      END SUBROUTINE INIT_DOMAIN

!-----------------------------------------------------------------------------

      SUBROUTINE CALC_OBS_HOUR

!***************************************************************************
! Subroutine CALC_OBS_HOUR computes an array of hours for each day of obs.
! If there is an obs in a particular gridbox on that day, it assigns the
! hour (0..23). If there isn't, OBS_HOUR stays initialized to -1.
! (mak, 12/14/05)
!***************************************************************************

      USE BPCH2_MOD,    ONLY : GET_TAU0
      USE TIME_MOD,     ONLY : GET_MONTH, GET_DAY,
     &                         GET_YEAR, GET_HOUR
      USE GRID_MOD,     ONLY : GET_IJ

#     include "CMN_SIZE"

      REAL*4   :: OBS_HOUR(IIPAR,JJPAR)
      REAL*8   :: TAU0, UTC
      INTEGER  :: W, I, J
      INTEGER  :: LON15, IIJJ(2)
      INTEGER  :: COUNT_GRID(IIPAR,JJPAR)

      !=================================================================
      ! CALC_OBS_HOUR begins here!
      !=================================================================

      ! Get TAU0 from the date (at 0GMT)
      TAU0 = GET_TAU0(GET_MONTH(), GET_DAY(), GET_YEAR())

      OBS_HOUR_MOPITT(:,:) = -1
      OBS_HOUR(:,:)        = 0
      COUNT_GRID(:,:)      = 0

      DO W = 1, T_DIM

         ! Compute local time:
         ! Local TIME = GMT + ( LONGITUDE / 15 ) since each hour of time
         ! corresponds to 15 degrees of LONGITUDE on the globe
         !============================================================
         LON15 = LONGITUDE(W) / 15d0
         UTC   = TAU(W) - TAU0 + LON15
         IF ( UTC < 0d0  )  UTC = UTC + 24
         IF ( UTC > 24d0 )  UTC = UTC - 24

         !Only consider day time MOPITT measurements
         !am = 12 hrs centered on 10:30am local time (so 4:30am-4:30pm)

#if defined( GRID05x0666 ) && defined( NESTED_CH ) && !defined( NESTED_SD )
         IF ( UTC >= 4.5 .and. UTC <= 16.5
     &    .and. LONGITUDE(W) > 70
     &    .and. LONGITUDE(W) < 150
     &    .and. LATITUDE(W)  > -11
     &    .and. LATITUDE(W)  < 55 ) THEN
#elif defined( GRID05x0666 ) && defined( NESTED_NA ) && !defined( NESTED_SD )
         IF ( UTC >= 4.5 .and. UTC <= 16.5
     &    .and. LONGITUDE(W) > -140
     &    .and. LONGITUDE(W) < -40
     &    .and. LATITUDE(W)  > 10
     &    .and. LATITUDE(W)  < 70 ) THEN
#elif defined( GRID05x0666 ) && (defined( NESTED_NA ) || defined( NESTED_CH )) && defined( NESTED_SD )
         IF ( UTC >= 4.5 .and. UTC <= 16.5
     &    .and. LONGITUDE(W) > -126
     &    .and. LONGITUDE(W) < -66
     &    .and. LATITUDE(W)  > 13
     &    .and. LATITUDE(W)  < 57 ) THEN
#else
         IF ( UTC >= 4.5 .and. UTC <= 16.5 ) THEN
#endif

         ! Get grid box of current record
            IIJJ  = GET_IJ( LONGITUDE(W), LATITUDE(W))
            I = IIJJ(1)
            J = IIJJ(2)

            ! If there's an obs, calculate the time
            IF ( CO_TOTAL_COLUMN(1,W) > 0d0 ) THEN

               COUNT_GRID(I,J) = COUNT_GRID(I,J) + 1d0
               !Add the time of obs, to be averaged and floored later
               OBS_HOUR(I,J) = OBS_HOUR(I,J) + MOPITT_GMT(W)

            ENDIF
         ENDIF
      ENDDO

      ! average obs_hour on the grid
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         IF ( COUNT_GRID(I,J) > 0d0 ) THEN

            OBS_HOUR_MOPITT(I,J) =
     &         FLOOR( OBS_HOUR(I,J) / COUNT_GRID(I,J) )

         ENDIF
      ENDDO
      ENDDO

      END SUBROUTINE CALC_OBS_HOUR

!----------------------------------------------------------------------------------------------

      FUNCTION ITS_TIME_FOR_MOPITT_OBS( ) RESULT( FLAG )

!******************************************************************************
!  Function ITS_TIME_FOR_MOPITT_OBS returns TRUE if there are observations
!  available for particular time (hour of a particular day) based on
!  the OBS_HOUR_MOPITT array which holds the hour of obs in each gridbox
!  (computed when file read in mop02_mod.f) (mak, 7/12/07)
!******************************************************************************

      USE TIME_MOD, ONLY : GET_HOUR, GET_MINUTE

#     include "CMN_SIZE"  ! Size params

      ! Function value
      LOGICAL :: FLAG

      INTEGER :: I,J

      !=================================================================
      ! ITS_TIME_FOR_MOPITT_OBS begins here!
      !=================================================================

      ! Default to false
      FLAG = .FALSE.

      DO J = 1,JJPAR
      DO I = 1,IIPAR
         IF( GET_HOUR()   == OBS_HOUR_MOPITT(I,J) .and.
     &       GET_MINUTE() == 0                          ) THEN

               PRINT*, 'obs_hour was', get_hour(), 'in box', I, J
               FLAG = .TRUE.

               !GOTO 11
               RETURN

         ENDIF
      ENDDO
      ENDDO

      END FUNCTION ITS_TIME_FOR_MOPITT_OBS

!----------------------------------------------------------------------------

      SUBROUTINE READ_MOP02( FILENAME )

!******************************************************************************
!  Subroutine READ_MOP02 allocates all module arrays and reads data into
!  them from the HDF file. (bmy, 7/2/03, zhe 1/19/11)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) FILENAME (CHARACTER) : Name of MOPITT file to read
!
!  NOTES:
!******************************************************************************

      ! References to F90 modules
#if defined( MOPITT_V5_CO_OBS )
      USE HdfSdModule
      USE HdfVdModule
#endif
      USE BPCH2_MOD,  ONLY : GET_TAU0
      USE ERROR_MOD,  ONLY : ALLOC_ERR

      ! Arguments
      CHARACTER(LEN=*), INTENT(IN) :: FILENAME

      ! Local variables

      INTEGER  :: as, i, year, month, day
      REAL*8   :: TAU0

#if defined( MOPITT_V5_CO_OBS )
      INTEGER  :: sId, vId, vSize, nDims, dims(4)
#endif
#if defined( MOPITT_V6_CO_OBS )

      INTEGER  :: he5_swopen, he5_swattach, he5_swfldinfo,
     &            he5_swrdfld, he5_swdetach, he5_swclose

      INTEGER    :: N, fId, swathid, rank
      INTEGER    :: ntype(4)
      INTEGER*8  :: dims8(4)
      INTEGER*8  :: START1(1), STRIDE1(1), EDGE1(1)
      INTEGER*8  :: START2(2), STRIDE2(2), EDGE2(2)
      INTEGER*8  :: START3(3), STRIDE3(3), EDGE3(3)
      INTEGER, PARAMETER  :: HE5F_ACC_RDONLY=101
      character*72 dimlist, maxdimlist

#endif

      !=================================================================
      ! Mop02Read begins here!
      !=================================================================

      ! Deallocate arrays
      CALL CLEANUP_MOP02

      ! Get date from filename (next to the '-' character)
      i    = INDEX( FILENAME, '-' )
      READ( FILENAME(i+1:i+4), '(i4)' ) year
      READ( FILENAME(i+5:i+6), '(i2)' ) month
      READ( FILENAME(i+7:i+8), '(i2)' ) day

      ! Get TAU0 from the date (at 0GMT)
      TAU0 = GET_TAU0( month, day, year )

#if defined( MOPITT_V6_CO_OBS )

      ! Opening an HDF-EOS5 swath file
      fId = he5_swopen(FILENAME, HE5F_ACC_RDONLY)

      ! Attaching to a swath object
      swathid = he5_swattach(fId, 'MOP02' )

      !=================================================================
      ! Seconds in day (1-D)
      !=================================================================

      ! Retrieve information about a specific geolocation or data field in the swath.
      as = he5_swfldinfo(swathid, "SecondsinDay", rank, dims8,
     &                   ntype, dimlist, maxdimlist)

      ! Allocate arrays
      ALLOCATE( SECONDS_IN_DAY( dims8(1) ), stat=as )
      IF ( as /= 0 ) CALL ALLOC_ERR( 'SECONDS_IN_DAY' )

      ALLOCATE( TAU( dims8(1) ), stat=as )
      IF ( as /= 0 ) CALL ALLOC_ERR( 'TAU' )

      ALLOCATE( MOPITT_GMT( dims8(1) ), stat=as )
      IF ( as /= 0 ) CALL ALLOC_ERR( 'MOPITT_GMT' )

      START1  = (/0/)
      STRIDE1 = (/1/)
      EDGE1   = (/dims8(1)/)

      ! Reading data from a data field
      as = he5_swrdfld(swathid, 'SecondsinDay',
     &     START1, STRIDE1, EDGE1, SECONDS_IN_DAY)

      ! Compute GMT of MOPITT observations
      MOPITT_GMT = ( DBLE( SECONDS_IN_DAY ) / 3600d0 )

      ! Compute TAU values for GAMAP from SECONDS_IN_DAY
      TAU        = MOPITT_GMT + TAU0

      ! Save time dimension in T_DIM
      T_DIM      = dims8(1)

      !=================================================================
      ! LONGITUDE (1-D)
      !=================================================================

      ! Retrieve information about a specific geolocation or data field in the swath.
      as = he5_swfldinfo(swathid, "Longitude", rank, dims8,
     &                   ntype, dimlist, maxdimlist)

      ! Allocate array
      ALLOCATE( LONGITUDE( dims8(1) ), stat=as )
      IF ( as /= 0 ) CALL ALLOC_ERR( 'LONGITUDE' )

      START1  = (/0/)
      STRIDE1 = (/1/)
      EDGE1   = (/dims8(1)/)

      ! Reading data from a data field
      as = he5_swrdfld(swathid, 'Longitude',
     &     START1, STRIDE1, EDGE1, LONGITUDE)

      !=================================================================
      ! LATITUDE (1-D)
      !=================================================================

      ! Retrieve information about a specific geolocation or data field in the swath.
      as = he5_swfldinfo(swathid, "Latitude", rank, dims8,
     &                   ntype, dimlist, maxdimlist)

      ! Allocate array
      ALLOCATE( LATITUDE( dims8(1) ), stat=as )
      IF ( as /= 0 ) CALL ALLOC_ERR( 'LATITUDE' )

      START1  = (/0/)
      STRIDE1 = (/1/)
      EDGE1   = (/dims8(1)/)

      ! Reading data from a data field
      as = he5_swrdfld(swathid, 'Latitude',
     &     START1, STRIDE1, EDGE1, LATITUDE)

      !=================================================================
      ! PRESSURE (1-D)
      !=================================================================

      ! Retrieve information about a specific geolocation or data field in the swath.
      as = he5_swfldinfo(swathid, "Pressure", rank, dims8,
     &                   ntype, dimlist, maxdimlist)

      ! Allocate array
      ALLOCATE(PRESSURE( dims8(1) ), stat=as )
      IF ( as /= 0 ) CALL ALLOC_ERR( 'PRESSURE' )

      START1  = (/0/)
      STRIDE1 = (/1/)
      EDGE1   = (/dims8(1)/)

      ! Reading data from a data field
      as = he5_swrdfld(swathid, 'Pressure',
     &     START1, STRIDE1, EDGE1, PRESSURE)

      ! Save PRESSURE dimension in Z_DIM
      Z_DIM = dims8(1)

      !=================================================================
      ! Cloud Description (1-D)
      !=================================================================

      ! Retrieve information about a specific geolocation or data field in the swath.
      as = he5_swfldinfo(swathid, "CloudDescription", rank, dims8,
     &                   ntype, dimlist, maxdimlist)

      ! Allocate array
      ALLOCATE( CLOUD_DES( dims8(1) ), stat=as )
      IF ( as /= 0 ) CALL ALLOC_ERR( 'CLOUD_DES' )

      START1  = (/0/)
      STRIDE1 = (/1/)
      EDGE1   = (/dims8(1)/)

      ! Reading data from a data field
      as = he5_swrdfld(swathid, 'CloudDescription',
     &     START1, STRIDE1, EDGE1, CLOUD_DES)

      !=================================================================
      ! Surface Index (1-D)
      !=================================================================

      ! Retrieve information about a specific geolocation or data field in the swath.
      as = he5_swfldinfo(swathid, "SurfaceIndex", rank, dims8,
     &                   ntype, dimlist, maxdimlist)

      ! Allocate array
      ALLOCATE(SURFACE_INDEX( dims8(1) ), stat=as )
      IF ( as /= 0 ) CALL ALLOC_ERR( 'SURFACE_INDEX' )

      START1  = (/0/)
      STRIDE1 = (/1/)
      EDGE1   = (/dims8(1)/)

      ! Reading data from a data field
      as = he5_swrdfld(swathid, 'SurfaceIndex',
     &     START1, STRIDE1, EDGE1, SURFACE_INDEX)

      !=================================================================
      ! Retrieval Bottom Pressure (1-D)
      !=================================================================

      ! Retrieve information about a specific geolocation or data field in the swath.
      as = he5_swfldinfo(swathid, "SurfacePressure", rank, dims8,
     &                   ntype, dimlist, maxdimlist)

      ! Allocate array
      ALLOCATE(BOTTOM_PRESSURE( dims8(1) ), stat=as )
      IF ( as /= 0 ) CALL ALLOC_ERR( 'BOTTOM_PRESSURE' )

      START1  = (/0/)
      STRIDE1 = (/1/)
      EDGE1   = (/dims8(1)/)

      ! Reading data from a data field
      as = he5_swrdfld(swathid, 'SurfacePressure',
     &     START1, STRIDE1, EDGE1, BOTTOM_PRESSURE)

      !=================================================================
      ! CO Mixing Ratio (3-D)
      !=================================================================

      ! Retrieve information about a specific geolocation or data field in the swath.
      as = he5_swfldinfo(swathid, "RetrievedCOMixingRatioProfile",
     &                   rank, dims8, ntype, dimlist, maxdimlist)

      ! Allocate array
      ALLOCATE( CO_MIXING_RATIO( dims8(1), dims8(2), dims8(3) ),
     &          stat=as )
      IF ( as /= 0 ) CALL ALLOC_ERR( 'CO_MIXING_RATIO' )

      START3  = (/0, 0, 0/)
      STRIDE3 = (/1, 1, 1/)
      EDGE3   = (/dims8(1), dims8(2), dims8(3)/)

      ! Reading data from a data field
      as = he5_swrdfld(swathid, 'RetrievedCOMixingRatioProfile',
     &     START3, STRIDE3, EDGE3, CO_MIXING_RATIO)

      !=================================================================
      ! SDATA field: CO Retrieval Bottom Mixing Ratio (2-D)
      !=================================================================

      ! Retrieve information about a specific geolocation or data field in the swath.
      as = he5_swfldinfo(swathid, "RetrievedCOSurfaceMixingRatio",
     &                   rank, dims8, ntype, dimlist, maxdimlist)

      ! Allocate array
      ALLOCATE( CO_RET_BOT_MIXING_RATIO( dims8(1), dims8(2) ), stat=as )
      IF ( as /= 0 ) CALL ALLOC_ERR( 'CO_RET_BOT_MIXING_RATIO' )

      START2  = (/0, 0/)
      STRIDE2 = (/1, 1/)
      EDGE2   = (/dims8(1), dims8(2)/)

      ! Reading data from a data field
      as = he5_swrdfld(swathid, 'RetrievedCOSurfaceMixingRatio',
     &     START2, STRIDE2, EDGE2, CO_RET_BOT_MIXING_RATIO)

      !=================================================================
      ! CO Total Column (2-D)
      !=================================================================

      ! Retrieve information about a specific geolocation or data field in the swath.
      as = he5_swfldinfo(swathid, "RetrievedCOTotalColumn",
     &                   rank, dims8, ntype, dimlist, maxdimlist)

      ! Allocate array
      ALLOCATE( CO_TOTAL_COLUMN( dims8(1), dims8(2) ), stat=as )
      IF ( as /= 0 ) CALL ALLOC_ERR( 'CO_TOTAL_COLUMN' )

      START2  = (/0, 0/)
      STRIDE2 = (/1, 1/)
      EDGE2   = (/dims8(1), dims8(2)/)

      ! Reading data from a data field
      as = he5_swrdfld(swathid, 'RetrievedCOTotalColumn',
     &     START2, STRIDE2, EDGE2, CO_TOTAL_COLUMN)

      !=================================================================
      ! Retrieval Averaging Kernel Matrix (3-D)
      !=================================================================

      ! Retrieve information about a specific geolocation or data field in the swath.
      as = he5_swfldinfo(swathid, "RetrievalAveragingKernelMatrix",
     &                   rank, dims8, ntype, dimlist, maxdimlist)

      ! Allocate array
      ALLOCATE( AVGKER( dims8(1), dims8(2), dims8(3) ), stat=as )
      IF ( as /= 0 ) CALL ALLOC_ERR( 'AVGKER' )

      START3  = (/0, 0, 0/)
      STRIDE3 = (/1, 1, 1/)
      EDGE3   = (/dims8(1), dims8(2), dims8(3)/)

      ! Reading data from a data field
      as = he5_swrdfld(swathid, 'RetrievalAveragingKernelMatrix',
     &     START3, STRIDE3, EDGE3, AVGKER)

      !=================================================================
      ! A Priori CO Mixing Ratio Profile (3-D)
      !=================================================================

      ! Retrieve information about a specific geolocation or data field in the swath.
      as = he5_swfldinfo(swathid, "APrioriCOMixingRatioProfile",
     &                   rank, dims8, ntype, dimlist, maxdimlist)

      ! Allocate array
      ALLOCATE( CO_MR_AP( dims8(1), dims8(2), dims8(3) ), stat=as )
      IF ( as /= 0 ) CALL ALLOC_ERR( 'CO_MR_AP' )

      START3  = (/0, 0, 0/)
      STRIDE3 = (/1, 1, 1/)
      EDGE3   = (/dims8(1), dims8(2), dims8(3)/)

      ! Reading data from a data field
      as = he5_swrdfld(swathid, 'APrioriCOMixingRatioProfile',
     &     START3, STRIDE3, EDGE3, CO_MR_AP)

      !=================================================================
      ! A Priori CO Surface Mixing Ratio (2-D)
      !=================================================================

      ! Retrieve information about a specific geolocation or data field in the swath.
      as = he5_swfldinfo(swathid, "APrioriCOSurfaceMixingRatio",
     &                   rank, dims8, ntype, dimlist, maxdimlist)

      ! Allocate array
      ALLOCATE( CO_MR_AP_BOTTOM( dims8(1), dims8(2)), stat=as )
      IF ( as /= 0 ) CALL ALLOC_ERR( 'CO_MR_AP_BOTTOM' )

      START2  = (/0, 0/)
      STRIDE2 = (/1, 1/)
      EDGE2   = (/dims8(1), dims8(2)/)

      ! Reading data from a data field
      as = he5_swrdfld(swathid, 'APrioriCOSurfaceMixingRatio',
     &     START2, STRIDE2, EDGE2, CO_MR_AP_BOTTOM)

      !=================================================================
      ! Retrieval Error Covariance Matrix (3-D)
      !=================================================================

      ! Retrieve information about a specific geolocation or data field in the swath.
      as = he5_swfldinfo(swathid, "RetrievalErrorCovarianceMatrix",
     &                   rank, dims8, ntype, dimlist, maxdimlist)

      ! Allocate array
      ALLOCATE( COVARIANCE( dims8(1), dims8(2), dims8(3) ), stat=as )
      IF ( as /= 0 ) CALL ALLOC_ERR( 'COVARIANCE' )

      START3  = (/0, 0, 0/)
      STRIDE3 = (/1, 1, 1/)
      EDGE3   = (/dims8(1), dims8(2), dims8(3)/)

      ! Reading data from a data field
      as = he5_swrdfld(swathid, 'RetrievalErrorCovarianceMatrix',
     &     START3, STRIDE3, EDGE3, COVARIANCE)

      ! Detaching from the swath object
      as = he5_swdetach(swathid)

      ! Closing the file
      as = he5_swclose(fId)


#endif  !MOPITT v6

#if defined( MOPITT_V5_CO_OBS )

      ! Open file for HDF-VDATA interface
      CALL vdOpen( FILENAME )

      !=================================================================
      ! VDATA field: Time (1-D)
      !=================================================================

      ! Open field for reading
      CALL vdOpenField( 'Seconds in Day', vId )

      ! Get size of field
      CALL vdGetFieldDim( vId, vSize )

      ! Allocate arrays
      ALLOCATE( SECONDS_IN_DAY( vSize ), stat=as )
      IF ( as /= 0 ) CALL ALLOC_ERR( 'SECONDS_IN_DAY' )

      ALLOCATE( TAU( vSize ), stat=as )
      IF ( as /= 0 ) CALL ALLOC_ERR( 'TAU' )

      ALLOCATE( MOPITT_GMT( vSize ), stat=as )
      IF ( as /= 0 ) CALL ALLOC_ERR( 'MOPITT_GMT' )

      ! Read data
      CALL vdGetData( vId, vSize, SECONDS_IN_DAY )

      ! Close field
      CALL vdCloseField( vId )

      ! Compute GMT of MOPITT observations
      MOPITT_GMT = ( DBLE( SECONDS_IN_DAY ) / 3600d0 )

      ! Compute TAU values for GAMAP from SECONDS_IN_DAY
      TAU       = MOPITT_GMT + TAU0

      ! Save time dimension in T_DIM
      T_DIM      = vSize

      !=================================================================
      ! VDATA field: LONGITUDE (1-D)
      !=================================================================

      ! Open field for reading
      CALL vdOpenField( 'Longitude', vId )

      ! Get size of field
      CALL vdGetFieldDim( vId, vSize )

      ! Allocate array
      ALLOCATE( LONGITUDE( vSize ), stat=as )
      IF ( as /= 0 ) CALL ALLOC_ERR( 'LONGITUDE' )

      ! Read data
      CALL vdGetData( vId, vSize, LONGITUDE )

      ! Close field
      CALL vdCloseField( vId )

      !=================================================================
      ! VDATA field: LATITUDE (1-D)
      !=================================================================

      ! Open field for reading
      CALL vdOpenField( 'Latitude', vId )

      ! Get size of field
      CALL vdGetFieldDim( vId, vSize )

      ! Allocate array
      ALLOCATE( LATITUDE( vSize ), stat=as )
      IF ( as /= 0 ) CALL ALLOC_ERR( 'LATITUDE' )

      ! Read data
      CALL vdGetData( vId, vSize, LATITUDE )

      ! Close field
      CALL vdCloseField( vId )

      !=================================================================
      ! VDATA field: Cloud Description (1-D)
      !=================================================================

      ! Open field for reading
      CALL vdOpenField( 'Cloud Description', vId )

      ! Get size of field
      CALL vdGetFieldDim( vId, vSize )

      ! Allocate array
      ALLOCATE( CLOUD_DES( vSize ), stat=as )
      IF ( as /= 0 ) CALL ALLOC_ERR( 'CLOUD_DES' )

      ! Read data
      CALL vdGetData( vId, vSize, CLOUD_DES )

      ! Close field
      CALL vdCloseField( vId )

      !=================================================================
      ! VDATA field: Surface Index (1-D)
      !=================================================================

      ! Open field for reading
      CALL vdOpenField( 'Surface Index', vId )

      ! Get size of field
      CALL vdGetFieldDim( vId, vSize )

      ! Allocate array
      ALLOCATE( SURFACE_INDEX( vSize ), stat=as )
      IF ( as /= 0 ) CALL ALLOC_ERR( 'SURFACE_INDEX' )

      ! Read data
      CALL vdGetData( vId, vSize, SURFACE_INDEX )

      ! Close field
      CALL vdCloseField( vId )

      !=================================================================
      ! VDATA field: PRESSURE (1-D)
      !=================================================================

      ! Open field for reading
      CALL vdOpenField( 'Pressure Grid', vId )

      ! Get size of field
      CALL vdGetFieldDim( vId, vSize )

      ! Allocate array
      ALLOCATE( PRESSURE( vSize ), stat=as )
      IF ( as /= 0 ) CALL ALLOC_ERR( 'PRESSURE' )

      ! Read data
      CALL vdGetData( vId, vSize, PRESSURE )

      ! Close field
      CALL vdCloseField( vId )

      ! Save PRESSURE dimension in Z_DIM
      Z_DIM = vSize

      !=================================================================
      ! VDATA field: Retrieval Bottom Pressure (1-D)
      !=================================================================

      ! Open field for reading
      CALL vdOpenField( 'Surface Pressure', vId )

      ! Get size of field
      CALL vdGetFieldDim( vId, vSize )

      ! Allocate array
      ALLOCATE( BOTTOM_PRESSURE( vSize ), stat=as )
      IF ( as /= 0 ) CALL ALLOC_ERR( 'BOTTOM_PRESSURE' )

      ! Read data
      CALL vdGetData( vId, vSize, BOTTOM_PRESSURE )

      ! Close field
      CALL vdCloseField( vId )

      ! Close HDF-VDATA interface
      CALL vdClose( FILENAME )


      ! Open file for HDF-SDATA interface
      CALL sdOpen( FILENAME )

      !=================================================================
      ! SDATA field: CO Mixing Ratio (3-D)
      !=================================================================

      ! Open field
      CALL sdOpenFieldByName( 'Retrieved CO Mixing Ratio Profile', sId )

      ! Get # of dimensions
      CALL sdGetFieldDims( sId, nDims, dims )

      ! Allocate array
      ALLOCATE( CO_MIXING_RATIO( dims(1), dims(2), dims(3) ), stat=as )
      IF ( as /= 0 ) CALL ALLOC_ERR( 'CO_MIXING_RATIO' )

      ! Read data
      CALL sdGetData( sId, dims(1), dims(2), dims(3), CO_MIXING_RATIO )

      ! Close field
      CALL sdCloseField( sId )

      !=================================================================
      ! SDATA field: CO Retrieval Bottom Mixing Ratio (2-D)
      !=================================================================

      ! Open field
      CALL sdOpenFieldByName( 'Retrieved CO Surface Mixing Ratio', sId )

      ! Get # of dimensions
      CALL sdGetFieldDims( sId, nDims, dims )

      ! Allocate array
      ALLOCATE( CO_RET_BOT_MIXING_RATIO( dims(1), dims(2) ), stat=as )
      IF ( as /= 0 ) CALL ALLOC_ERR( 'CO_RET_BOT_MIXING_RATIO' )

      ! Read data
      CALL sdGetData( sId, dims(1), dims(2), CO_RET_BOT_MIXING_RATIO )

      ! Close field
      CALL sdCloseField( sId )

      !=================================================================
      ! SDATA field: CO Total Column (2-D)
      !=================================================================

      ! Open field

      CALL sdOpenFieldByName( 'Retrieved CO Total Column', sId )

      ! Get # of dimensions
      CALL sdGetFieldDims( sId, nDims, dims )

      ! Allocate array
      ALLOCATE( CO_TOTAL_COLUMN( dims(1), dims(2) ), stat=as )
      IF ( as /= 0 ) CALL ALLOC_ERR( 'CO_TOTAL_COLUMN' )

      ! Read data
      CALL sdGetData( sId, dims(1), dims(2), CO_TOTAL_COLUMN )

      ! Close field
      CALL sdCloseField( sId )

      !=================================================================
      ! SDATA field: Retrieval Averaging Kernel Matrix (3-D)
      !=================================================================

      ! Open field
      CALL sdOpenFieldByName( 'Retrieval Averaging Kernel Matrix', sId )

      ! Get # of dimensions
      CALL sdGetFieldDims( sId, nDims, dims )

      ! Allocate array
      ALLOCATE( AVGKER( dims(1), dims(2), dims(3) ), stat=as )
      IF ( as /= 0 ) CALL ALLOC_ERR( 'AVGKER' )

      ! Read data
      CALL sdGetData( sId, dims(1), dims(2), dims(3), AVGKER )

      ! Close field
      CALL sdCloseField( sId )

      !=================================================================
      ! SDATA field: A Priori CO Mixing Ratio Profile (3-D)
      !=================================================================

      ! Open field
      CALL sdOpenFieldByName( 'A Priori CO Mixing Ratio Profile', sId )

      ! Get # of dimensions
      CALL sdGetFieldDims( sId, nDims, dims )

      ! Allocate array
      ALLOCATE( CO_MR_AP( dims(1), dims(2), dims(3) ), stat=as )
      IF ( as /= 0 ) CALL ALLOC_ERR( 'CO_MR_AP' )

      ! Read data
      CALL sdGetData( sId, dims(1), dims(2), dims(3), CO_MR_AP )

      ! Close field
      CALL sdCloseField( sId )

      !=================================================================
      ! SDATA field: A Priori CO Surface Mixing Ratio (2-D)
      !=================================================================

      ! Open field
      CALL sdOpenFieldByName( 'A Priori CO Surface Mixing Ratio', sId )

      ! Get # of dimensions
      CALL sdGetFieldDims( sId, nDims, dims )

      ! Allocate array
      ALLOCATE( CO_MR_AP_BOTTOM( dims(1), dims(2)), stat=as )
      IF ( as /= 0 ) CALL ALLOC_ERR( 'CO_MR_AP_BOTTOM' )

      ! Read data
      CALL sdGetData( sId, dims(1), dims(2), CO_MR_AP_BOTTOM )

      ! Close field
      CALL sdCloseField( sId )

      !=================================================================
      ! SDATA field: Retrieval Error Covariance Matrix (3-D)
      !=================================================================

      ! Open field
      CALL sdOpenFieldByName( 'Retrieval Error Covariance Matrix', sId )

      ! Get # of dimensions
      CALL sdGetFieldDims( sId, nDims, dims )

      ! Allocate array
      ALLOCATE( COVARIANCE ( dims(1), dims(2), dims(3) ), stat=as )
      IF ( as /= 0 ) CALL ALLOC_ERR( 'COVARIANCE' )

      ! Read data
      CALL sdGetData( sId, dims(1), dims(2), dims(3), COVARIANCE )

      ! Close field
      CALL sdCloseField( sId )

      ! Close file and quit
      CALL sdClose( FILENAME )

#endif  !MOPITT v5

      ! Return to calling program
      END SUBROUTINE READ_MOP02

!------------------------------------------------------------------------------------

      SUBROUTINE APRIORI_MOP02

!******************************************************************************
! Calculates averaging kernels
! from a MOPITT retrieved covariance matrix on 7-levels,
! using apriori profile and covariance matrix for v3 MOPITT retrievals
!******************************************************************************

      CHARACTER(LEN=255) AP_FILENAME1
      CHARACTER(LEN=255) AP_FILENAME
      CHARACTER(LEN=80) LINE
      INTEGER :: ILEV,JLEV
      REAL*4, ALLOCATABLE :: DATA1(:)

      ! Read CO a priori profile and covariance matrix
      AP_FILENAME1 = 'mopitt_v3_apriori.dat'
      AP_FILENAME = TRIM(AP_FILENAME1)

      OPEN ( FILE = TRIM(AP_FILENAME),
     &     UNIT = 20)
!     &     ,form='unformatted')
      READ(20,101) LINE
      READ(20,101) LINE
      READ(20,*) NLEV_AP
!      print*,'NLEV_AP=',NLEV_AP
      ALLOCATE( PLEV_AP ( NLEV_AP ) )
      ALLOCATE( CH4_MR_AP ( NLEV_AP ) )
      ALLOCATE( COV_CO_AP ( NLEV_AP , NLEV_AP) )
      READ(20,101) LINE
      READ(20,*) PLEV_AP
      !PRINT*, PLEV_AP
      READ(20,101) LINE
      READ(20,*) CO_MR_AP
      READ(20,101) LINE
      READ(20,*) CH4_MR_AP
      READ(20,101) LINE
      DO ILEV= 1, NLEV_AP
         ALLOCATE( DATA1 ( ILEV ) )
         READ(20,101) LINE
         READ(20,*) DATA1
         COV_CO_AP(ILEV,1:ILEV) = DATA1
         COV_CO_AP(1:ILEV,ILEV) = DATA1
         IF ( ALLOCATED( DATA1 ) ) DEALLOCATE( DATA1 )
      ENDDO

      CLOSE(20)

 101  FORMAT(A80)

      END SUBROUTINE APRIORI_MOP02
!-----------------------------------------------------------------------------------------

      SUBROUTINE READ_ERROR_VARIANCE
!
!******************************************************************************
!  Subroutine READ_ERROR_VARIANCE reads observation error from binary punch files
!  (zhe 4/20/11)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD
      USE TIME_MOD,   ONLY : GET_TAUb

      IMPLICIT NONE

#     include "CMN_SIZE"    ! Size parameters

      ! Local variables
      CHARACTER(LEN=255)   :: FILENAME

      !=================================================================
      ! READ_ERROR_VARIANCE begins here!
      !=================================================================

      ! Filename
        FILENAME = TRIM( 'OBS_ERR_' ) // GET_RES_EXT()

      ! Echo some information to the standard output
        WRITE( 6, 110 ) TRIM( FILENAME )
 110    FORMAT( '     - READ_ERROR_VARIANCE: Reading ERR_PERCENT
     &                from: ', a )

      ! Read data from the binary punch file
        CALL READ_BPCH2( FILENAME, 'IJ-AVG-$', 1,
     &           GET_TAUb(),    IGLOB,     JGLOB,
     &           1,  ERR_PERCENT,  QUIET=.TRUE. )

      ! Return to calling program
      END SUBROUTINE READ_ERROR_VARIANCE

!------------------------------------------------------------------------------

      SUBROUTINE INFO_MOP02( FILENAME )
!
!******************************************************************************
!  Subroutine INFO_MOP02 Info prints info about all VDATA and SDATA fields
!  contained within the MOPITT HDF file. (bmy, 7/3/03, 4/27/05)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) FILENAME (CHARACTER) : Name of MOPITT file to read
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE HdfSdModule
      USE HdfVdModule

      ! Arguments
      CHARACTER(LEN=*), INTENT(IN) :: FILENAME

      !=================================================================
      ! INFO_MOP02 begins here!
      !=================================================================

      ! Print HDF-VDATA variables
      CALL vdOpen( FILENAME )
      CALL vdPrintInfo
      CALL vdClose( FILENAME )

      ! Print HDF-SDATA variables
      CALL sdOpen( FILENAME )
      CALL sdPrintInfo
      CALL sdClose( FILENAME )

      ! Return to calling program
      END SUBROUTINE INFO_MOP02

!-----------------------------------------------------------------------------

      SUBROUTINE CALL_INV_S(NLEV_RET, W, INV_S, SY)

!***************************************************************************
! Subroutine CALL_INV_S compute the inverse of convariance matrix
! Uniform error matrix, if SY > 0
! MOPITT error covariance matrix, if SY <= 0
!***************************************************************************

!     SY = (SY / LOG(10.0)) **2

      REAL*8   :: OBS_S(NLEV_RET,NLEV_RET)
      REAL*8   :: U(NLEV_RET,NLEV_RET)
      REAL*8   :: VT(NLEV_RET,NLEV_RET)
      REAL*8   :: S(NLEV_RET), SY, TMP
      REAL*8   :: TEST(NLEV_RET,NLEV_RET)
      REAL*4   :: INV_S(NLEV_RET,NLEV_RET)
      INTEGER  :: W, I, J, K, NLEV_RET
      INTEGER  :: ILEV, ILEV2, JLEV, JLEV2



      !=================================================================
      ! Unifrom observation error, if SY is defined
      !=================================================================
      IF (SY > 0) THEN

         INV_S(:,:) = 0.0
         TMP = 1.0 / ((SY / LOG(10.0))**2)

         DO I = 1, NLEV_RET
            INV_S(I,I) = TMP
         ENDDO

         RETURN

      ENDIF

      !Remove bad levels from covariance matrix
      IF ( NLEV_RET < 10 ) THEN
         DO ILEV = 1, NLEV_RET
            ILEV2 = ILEV + ( 10 - NLEV_RET )
            DO JLEV =1, NLEV_RET
               JLEV2 = JLEV + ( 10 - NLEV_RET)
               OBS_S(ILEV,JLEV) =
     &              COVARIANCE(ILEV2,JLEV2,W)
            ENDDO
         ENDDO
      ELSE
         OBS_S(:,:) = COVARIANCE(:,:,W)
      ENDIF

      ! Add a bit to the diagonal to regularize the inversion
      DO I = 1, NLEV_RET
         OBS_S(I,I) = OBS_S(I,I) + 0.001D0
      ENDDO

      CALL SVD( OBS_S, NLEV_RET, U, S, VT)

      ! U = S^-1 * U^T
      DO I = 1, NLEV_RET
         DO J = 1, NLEV_RET
            TEST(I,J) = U(J,I) / S(I)
         ENDDO
      ENDDO
      U = TEST

      ! S_OER_INV = V * S^-1 * U^T
      DO I = 1, NLEV_RET
         DO J = 1, NLEV_RET
            TMP = 0d0
            DO K = 1, NLEV_RET
               TMP = TMP + VT(K,I) * U(K,J)
            ENDDO
            INV_S(I,J) = TMP
         ENDDO
      ENDDO


      END SUBROUTINE CALL_INV_S

!-----------------------------------------------------------------------------

      SUBROUTINE SVD(A,N,U,S,VT)
!
!******************************************************************************
!  Subroutine SVD is a driver for the LAPACK SVD routine DGESVD. (dkh, 05/04/10)
!
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) A   (REAL*8) :  N x N matrix to decompose
!  (2 ) N  (INTEGER) :  N is dimension of A
!
!  Arguments as Output:
!  ============================================================================
!  (1 ) U   (REAL*8) :  Array of left singular vectors
!  (2 ) S   (REAL*8) :  Vector of singular values
!  (3 ) VT  (REAL*8) :  Array of right singular vectors, TRANSPOSED
!
!
!  NOTES:
!
*  Copyright (C) 2009-2010 Intel Corporation. All Rights Reserved.
*  The information and material ("Material") provided below is owned by Intel
*  Corporation or its suppliers or licensors, and title to such Material remains
*  with Intel Corporation or its suppliers or licensors. The Material contains
*  proprietary information of Intel or its suppliers and licensors. The Material
*  is protected by worldwide copyright laws and treaty provisions. No part of
*  the Material may be copied, reproduced, published, uploaded, posted,
*  transmitted, or distributed in any way without Intel's prior express written
*  permission. No license under any patent, copyright or other intellectual
*  property rights in the Material is granted to or conferred upon you, either
*  expressly, by implication, inducement, estoppel or otherwise. Any license
*  under such intellectual property rights must be express and approved by Intel
*  in writing.
*  =============================================================================
*
*  DGESVD Example.
*  ==============
*
*  Program computes the singular value decomposition of a general
*  rectangular matrix A:
*
*    8.79   9.93   9.83   5.45   3.16
*    6.11   6.91   5.04  -0.27   7.98
*   -9.15  -7.93   4.86   4.85   3.01
*    9.57   1.64   8.83   0.74   5.80
*   -3.49   4.02   9.80  10.00   4.27
*    9.84   0.15  -8.99  -6.02  -5.31
*
*  Description.
*  ============
*
*  The routine computes the singular value decomposition (SVD) of a real
*  m-by-n matrix A, optionally computing the left and/or right singular
*  vectors. The SVD is written as
*
*  A = U*SIGMA*VT
*
*  where SIGMA is an m-by-n matrix which is zero except for its min(m,n)
*  diagonal elements, U is an m-by-m orthogonal matrix and VT (V transposed)
*  is an n-by-n orthogonal matrix. The diagonal elements of SIGMA
*  are the singular values of A; they are real and non-negative, and are
*  returned in descending order. The first min(m, n) columns of U and V are
*  the left and right singular vectors of A.
*
*  Note that the routine returns VT, not V.
*
*  Example Program Results.
*  ========================
*
* DGESVD Example Program Results
*
* Singular values
*  27.47  22.64   8.56   5.99   2.01
*
* Left singular vectors (stored columnwise)
*  -0.59   0.26   0.36   0.31   0.23
*  -0.40   0.24  -0.22  -0.75  -0.36
*  -0.03  -0.60  -0.45   0.23  -0.31
*  -0.43   0.24  -0.69   0.33   0.16
*  -0.47  -0.35   0.39   0.16  -0.52
*   0.29   0.58  -0.02   0.38  -0.65
*
* Right singular vectors (stored rowwise)
*  -0.25  -0.40  -0.69  -0.37  -0.41
*   0.81   0.36  -0.25  -0.37  -0.10
*  -0.26   0.70  -0.22   0.39  -0.49
*   0.40  -0.45   0.25   0.43  -0.62
*  -0.22   0.14   0.59  -0.63  -0.44
*  =============================================================================
!******************************************************************************
!
      ! Arguements
      INTEGER,INTENT(IN)     :: N
      REAL*8, INTENT(IN)     :: A(N,N)
      REAL*8, INTENT(OUT)    :: U(N,N)
      REAL*8, INTENT(OUT)    :: S(N)
      REAL*8, INTENT(OUT)    :: VT(N,N)

      ! Local variables
      INTEGER, PARAMETER     :: LWMAX = 350
      INTEGER                :: INFO, LWORK
      DOUBLE PRECISION       :: WORK( LWMAX )

*     .. External Subroutines ..
      EXTERNAL               :: DGESVD

*     .. Intrinsic Functions ..
      INTRINSIC              :: INT, MIN

      !=================================================================
      ! SVD begins here!
      !=================================================================

*     .. Executable Statements ..
      !WRITE(*,*)'DGESVD Example Program Results'
*
*     Query the optimal workspace.
*
      LWORK = -1
      CALL DGESVD( 'All', 'All', N, N, A, N, S, U, N, VT, N,
     $             WORK, LWORK, INFO )
      LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
*
*     Compute SVD.
*
      CALL DGESVD( 'All', 'All', N, N, A, N, S, U, N, VT, N,
     $             WORK, LWORK, INFO )
*
*     Check for convergence.
*
      IF( INFO.GT.0 ) THEN
         WRITE(*,*)'The algorithm computing SVD failed to converge.'
         STOP
      END IF

!  Uncomment the following to print out singlular values, vectors (dkh, 05/04/10)
!!
!!     Print singular values.
!!
!      CALL PRINT_MATRIX( 'Singular values', 1, N, S, 1 )
!!
!!     Print left singular vectors.
!!
!      CALL PRINT_MATRIX( 'Left singular vectors (stored columnwise)',
!     $                   N, N, U, N   )
!!
!!     Print right singular vectors.
!!
!      CALL PRINT_MATRIX( 'Right singular vectors (stored rowwise)',
!     $                   N, N, VT, N    )

      ! Return to calling program
      END SUBROUTINE SVD

!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_MOP02
!
!******************************************************************************
!  Subroutine CLEANUP_MOP02 deallocates all module arrays (bmy, 4/27/05)
!
!  NOTES:
!******************************************************************************
!
      !=================================================================
      ! CLEANUP_MOP02 begins here!
      !=================================================================
      IF ( ALLOCATED( LATITUDE        ) ) DEALLOCATE( LATITUDE        )
      IF ( ALLOCATED( LONGITUDE       ) ) DEALLOCATE( LONGITUDE       )
      IF ( ALLOCATED( PRESSURE        ) ) DEALLOCATE( PRESSURE        )
      IF ( ALLOCATED( CLOUD_DES       ) ) DEALLOCATE( CLOUD_DES       )
      IF ( ALLOCATED( SURFACE_INDEX   ) ) DEALLOCATE( SURFACE_INDEX   )
      IF ( ALLOCATED( TAU             ) ) DEALLOCATE( TAU             )
      IF ( ALLOCATED( SECONDS_IN_DAY  ) ) DEALLOCATE( SECONDS_IN_DAY  )
      IF ( ALLOCATED( MOPITT_GMT      ) ) DEALLOCATE( MOPITT_GMT      )
      IF ( ALLOCATED( BOTTOM_PRESSURE ) ) DEALLOCATE( BOTTOM_PRESSURE )
      IF ( ALLOCATED( CO_MIXING_RATIO ) ) DEALLOCATE( CO_MIXING_RATIO )
      IF ( ALLOCATED( COVARIANCE      ) ) DEALLOCATE( COVARIANCE      )

      IF ( ALLOCATED( CO_RET_BOT_MIXING_RATIO)) THEN
         DEALLOCATE( CO_RET_BOT_MIXING_RATIO )
      ENDIF

      IF ( ALLOCATED( CO_TOTAL_COLUMN ) ) DEALLOCATE( CO_TOTAL_COLUMN )
      IF ( ALLOCATED( AVGKER          ) ) DEALLOCATE( AVGKER          )
      IF ( ALLOCATED( PLEV_AP         ) ) DEALLOCATE( PLEV_AP         )
      IF ( ALLOCATED( CO_MR_AP        ) ) DEALLOCATE( CO_MR_AP        )
      IF ( ALLOCATED( CH4_MR_AP       ) ) DEALLOCATE( CH4_MR_AP       )
      IF ( ALLOCATED( CO_MR_AP        ) ) DEALLOCATE( CO_MR_AP        )
      IF ( ALLOCATED( CO_MR_AP_BOTTOM ) ) DEALLOCATE( CO_MR_AP_BOTTOM )
      IF ( ALLOCATED( COV_CO_AP       ) ) DEALLOCATE( COV_CO_AP       )

      ! Return to calling program
      END SUBROUTINE CLEANUP_MOP02

!---------------------------------------------------------------------------------------------------

      SUBROUTINE MAKE_MOPITT_BIAS_FILE_HDF5(FILE_ID)

      USE HDF5

      USE GRID_MOD,            ONLY : GET_XMID, GET_YMID
      USE DIRECTORY_ADJ_MOD,   ONLY : OPTDATA_DIR
      USE ADJ_ARRAYS_MOD,      ONLY : N_CALC
      USE ADJ_ARRAYS_MOD,      ONLY : EXPAND_NAME
      USE LOGICAL_ADJ_MOD,     ONLY : LSAT_HDF_L2, LSAT_HDF_L3

      INTEGER(HID_T) :: FILE_ID

      CHARACTER(LEN=255) :: LON_NAME, LAT_NAME, LEV_NAME
      CHARACTER(LEN=255) :: MOPITT_CO_NAME
      CHARACTER(LEN=255) :: MOPITT_GC_CO_NAME
      CHARACTER(LEN=255) :: MOPITT_BIAS_NAME
      CHARACTER(LEN=255) :: MOPITT_COUNT_NAME

      CHARACTER(LEN=255) :: MOPITT_CO_LONGNAME
      CHARACTER(LEN=255) :: MOPITT_GC_CO_LONGNAME
      CHARACTER(LEN=255) :: MOPITT_BIAS_LONGNAME
      CHARACTER(LEN=255) :: MOPITT_COUNT_LONGNAME

      CHARACTER(LEN=255) :: MOPITT_CO_UNIT
      CHARACTER(LEN=255) :: MOPITT_GC_CO_UNIT
      CHARACTER(LEN=255) :: MOPITT_BIAS_UNIT
      CHARACTER(LEN=255) :: MOPITT_COUNT_UNIT

      CHARACTER(LEN=255) :: LON_LONGNAME, LAT_LONGNAME, LEV_LONGNAME
      CHARACTER(LEN=255) :: LON_UNIT, LAT_UNIT, LEV_UNIT

      INTEGER(HID_T) :: SPACE_LON, SPACE_LAT, SPACE_LEV
      INTEGER(HID_T) :: SPACE_RAW_1D, SPACE_RAW_2D
      INTEGER(HID_T) :: LON_ID, LAT_ID, LEV_ID
      INTEGER(HID_T) :: SPACE_MOPITT, DSET_MOPITT_CO_ID
      INTEGER(HID_T) :: DSET_MOPITT_GC_CO_ID
      INTEGER(HID_T) :: DSET_MOPITT_BIAS_ID
      INTEGER(HID_T) :: DSET_MOPITT_COUNT_ID

      INTEGER(HID_T) :: LON_RAW_ID, LAT_RAW_ID, TIME_RAW_ID
      INTEGER(HID_T) :: DSET_MOP_CO_RAW_ID, DSET_GC_CO_RAW_ID

      INTEGER(HID_T) :: ASPACE_ID, ATYPE_ID, ATT_ID
      INTEGER(HSIZE_T) :: ADIMS(1)

      INTEGER(HID_T) :: LEVEL3_GROUP_ID, LEVEL2_GROUP_ID
      INTEGER(HID_T) :: MOPITT_GROUP_ID, DATA_GROUP_ID, GRID_GROUP_ID

      INTEGER(HSIZE_T) :: DIMS(3), DIM_LON(1), DIM_LAT(1), DIM_LEV(1)
      INTEGER(HSIZE_T) :: DIM_RAW_1D(1), DIM_RAW_2D(2)

      INTEGER :: HDF_ERR
      INTEGER :: RANK = 3

      INTEGER :: I,J,L
      REAL*4 :: MISS_VAL = -999.9
      REAL*4 :: LON_VALS(IIPAR), LAT_VALS(JJPAR), LEV_VALS(10)

      ! populate lon & lat arrays

      DO I=1,IIPAR
         LON_VALS(I)=GET_XMID(I)
      ENDDO

      DO J=1,JJPAR
         LAT_VALS(J)=GET_YMID(J)
      ENDDO

      DO I=1,10
         LEV_VALS(I)=1100-I*100 ! mkeller: needs to be handled better
      ENDDO

      DO L=1,10
         DO J=1,JJPAR
            DO I=1,IIPAR

               IF(MOPITT_COUNT_SOBS(I,J,L)>0d0) THEN

                  MOPITT_MOP_CO_SOBS(I,J,L) =
     & REAL(MOPITT_MOP_CO_SOBS(I,J,L)/MOPITT_COUNT_SOBS(I,J,L))*1E9
                  MOPITT_GC_CO_SOBS(I,J,L) =
     & REAL(MOPITT_GC_CO_SOBS(I,J,L)/MOPITT_COUNT_SOBS(I,J,L))*1E9
                  MOPITT_BIAS_SOBS(I,J,L) =
     & REAL(MOPITT_BIAS_SOBS(I,J,L)/MOPITT_COUNT_SOBS(I,J,L))*1E9

               ELSE

                  MOPITT_MOP_CO_SOBS(I,J,L) = MISS_VAL
                  MOPITT_GC_CO_SOBS(I,J,L) = MISS_VAL
                  MOPITT_BIAS_SOBS(I,J,L) = MISS_VAL
                  !MOPITT_CHI_SQUARED(I,J,L) = MISS_VAL

               ENDIF

            ENDDO
         ENDDO
      ENDDO

      DIMS(1) = IIPAR
      DIMS(2) = JJPAR
      DIMS(3) = 10

      ADIMS(1) = 1

      DIM_LON = IIPAR
      DIM_LAT = JJPAR
      DIM_LEV = 10

      DIM_RAW_1D = FLEX_LON%CURRENT_N      
      DIM_RAW_2D(1) = MAXLEV
      DIM_RAW_2D(2) = FLEX_LON%CURRENT_N

      LON_NAME = "/MOPITT/Level3/Grid/Longitude"
      LAT_NAME = "/MOPITT/Level3/Grid/Latitude"
      LEV_NAME = "/MOPITT/Level3/Grid/Level"

      MOPITT_CO_NAME = "/MOPITT/Level3/Data/MOPITT_CO"
      MOPITT_GC_CO_NAME = "/MOPITT/Level3/Data/MOPITT_GC_CO"
      MOPITT_BIAS_NAME = "/MOPITT/Level3/Data/MOPITT_BIAS"
      MOPITT_COUNT_NAME = "/MOPITT/Level3/Data/MOPITT_COUNT"

      MOPITT_CO_LONGNAME = "Mean MOPITT CO Profiles"
      MOPITT_GC_CO_LONGNAME = "Mean MOPITT GEOS_Chem CO Profiles"
      MOPITT_BIAS_LONGNAME = "Mean MOPITT CO bias Profiles"
      MOPITT_COUNT_LONGNAME = "MOPITT data count"

      MOPITT_CO_UNIT = "ppbv"
      MOPITT_GC_CO_UNIT = "ppbv"
      MOPITT_BIAS_UNIT = "ppbv"
      MOPITT_COUNT_UNIT = "1"

      LON_LONGNAME = "Longitudes"
      LAT_LONGNAME = "Latitudes"
      LEV_LONGNAME = "Vertical Levels. "
     &     //"1000hPa is a placeholder for the surface level here."

      LON_UNIT = "degrees"
      LAT_UNIT = "degrees"
      LEV_UNIT = "hPa"

      ! open HDF5 interface

      CALL H5OPEN_F(HDF_ERR)

      ! create group structure in file

      CALL H5GCREATE_F(FILE_ID,"MOPITT",MOPITT_GROUP_ID,HDF_ERR)
      
      IF(LSAT_HDF_L2) THEN
         CALL H5GCREATE_F(MOPITT_GROUP_ID,"Level2",
     &        LEVEL2_GROUP_ID,HDF_ERR)
      ENDIF
      
      IF(LSAT_HDF_L3) THEN
         CALL H5GCREATE_F(MOPITT_GROUP_ID,"Level3",
     &        LEVEL3_GROUP_ID,HDF_ERR)
         CALL H5GCREATE_F(LEVEL3_GROUP_ID,"Data",DATA_GROUP_ID,HDF_ERR)
         CALL H5GCREATE_F(LEVEL3_GROUP_ID,"Grid",GRID_GROUP_ID,HDF_ERR)

         ! write grid information
   
         CALL H5SCREATE_SIMPLE_F(1,DIM_LON,SPACE_LON,HDF_ERR)
         CALL H5SCREATE_SIMPLE_F(1,DIM_LAT,SPACE_LAT,HDF_ERR)
         CALL H5SCREATE_SIMPLE_F(1,DIM_LEV,SPACE_LEV,HDF_ERR)
   
         CALL H5DCREATE_F(GRID_GROUP_ID, LON_NAME, H5T_IEEE_F32LE,
     &        SPACE_LON, LON_ID, HDF_ERR)
         CALL H5DCREATE_F(GRID_GROUP_ID, LAT_NAME, H5T_IEEE_F32LE,
     &        SPACE_LAT, LAT_ID, HDF_ERR)
         CALL H5DCREATE_F(GRID_GROUP_ID, LEV_NAME, H5T_IEEE_F32LE,
     &        SPACE_LEV, LEV_ID, HDF_ERR)
   
         CALL H5DWRITE_F(LON_ID, H5T_NATIVE_REAL, LON_VALS,
     &        DIM_LON, HDF_ERR)
         CALL H5DWRITE_F(LAT_ID, H5T_NATIVE_REAL, LAT_VALS,
     &        DIM_LAT, HDF_ERR)
         CALL H5DWRITE_F(LEV_ID, H5T_NATIVE_REAL, LEV_VALS,
     &        DIM_LEV, HDF_ERR)
   
         CALL WRITE_ATTRIBUTES(LON_ID,LON_LONGNAME,LON_UNIT)
         CALL WRITE_ATTRIBUTES(LAT_ID,LAT_LONGNAME,LAT_UNIT)
         CALL WRITE_ATTRIBUTES(LEV_ID,LEV_LONGNAME,LEV_UNIT)
   
         CALL H5DCLOSE_F(LON_ID, HDF_ERR)
         CALL H5DCLOSE_F(LAT_ID, HDF_ERR)
         CALL H5DCLOSE_F(LEV_ID, HDF_ERR)
   
         CALL H5SCLOSE_F(SPACE_LON, HDF_ERR)
         CALL H5SCLOSE_F(SPACE_LAT, HDF_ERR)
         CALL H5SCLOSE_F(SPACE_LEV, HDF_ERR)
   
         ! create dataspace for MOPITT diagnostics
   
         CALL H5SCREATE_SIMPLE_F(RANK,DIMS,SPACE_MOPITT,HDF_ERR)
   
         ! create all datasets as little-endian 32 bit IEEE float
   
         ! write MOPITT CO concentrations
   
         CALL H5DCREATE_F(DATA_GROUP_ID,MOPITT_CO_NAME,H5T_IEEE_F32LE,
     &        SPACE_MOPITT, DSET_MOPITT_CO_ID, HDF_ERR)
   
         CALL H5DWRITE_F(DSET_MOPITT_CO_ID, H5T_NATIVE_REAL,
     &        REAL(MOPITT_MOP_CO_SOBS), DIMS, HDF_ERR)
   
         CALL WRITE_ATTRIBUTES(DSET_MOPITT_CO_ID,MOPITT_CO_LONGNAME,
     &        MOPITT_CO_UNIT)
   
         CALL H5DCLOSE_F(DSET_MOPITT_CO_ID,HDF_ERR)
   
         ! write MOPITT_GC CO concentrations
   
         CALL H5DCREATE_F(DATA_GROUP_ID,MOPITT_GC_CO_NAME,H5T_IEEE_F32LE
     &        ,SPACE_MOPITT, DSET_MOPITT_GC_CO_ID, HDF_ERR)
   
         CALL H5DWRITE_F(DSET_MOPITT_GC_CO_ID, H5T_NATIVE_REAL,
     &        REAL(MOPITT_GC_CO_SOBS), ADIMS, HDF_ERR)
   
         CALL WRITE_ATTRIBUTES(DSET_MOPITT_GC_CO_ID,
     &        MOPITT_GC_CO_LONGNAME, MOPITT_GC_CO_UNIT)
   
         CALL H5DCLOSE_F(DSET_MOPITT_GC_CO_ID,HDF_ERR)
   
         ! write MOPITT_GC CO bias
   
         CALL H5DCREATE_F(DATA_GROUP_ID,MOPITT_BIAS_NAME,H5T_IEEE_F32LE,
     &        SPACE_MOPITT, DSET_MOPITT_BIAS_ID, HDF_ERR)
   
         CALL H5DWRITE_F(DSET_MOPITT_BIAS_ID, H5T_NATIVE_REAL,
     &        REAL(MOPITT_BIAS_SOBS), DIMS, HDF_ERR)
   
         CALL WRITE_ATTRIBUTES(DSET_MOPITT_BIAS_ID,MOPITT_BIAS_LONGNAME,
     &        MOPITT_BIAS_UNIT)
   
         CALL H5DCLOSE_F(DSET_MOPITT_BIAS_ID,HDF_ERR)
   
         ! write MOPITT_GC CO count
   
         CALL H5DCREATE_F(DATA_GROUP_ID,MOPITT_COUNT_NAME,H5T_IEEE_F32LE
     &       ,SPACE_MOPITT, DSET_MOPITT_COUNT_ID, HDF_ERR)
   
         CALL H5DWRITE_F(DSET_MOPITT_COUNT_ID, H5T_NATIVE_REAL,
     &        REAL(MOPITT_COUNT_SOBS), DIMS, HDF_ERR)
   
         CALL WRITE_ATTRIBUTES(DSET_MOPITT_COUNT_ID,
     &        MOPITT_COUNT_LONGNAME,MOPITT_COUNT_UNIT)
   
         CALL H5DCLOSE_F(DSET_MOPITT_COUNT_ID,HDF_ERR)

      ENDIF
      !------------------------------------------------------------------
      IF(LSAT_HDF_L2) THEN
         ! create dataspace for raw 1D (Level2) diagnostics
   
         CALL H5SCREATE_SIMPLE_F(1,DIM_RAW_1D,SPACE_RAW_1D,HDF_ERR)
   
         ! write raw longitudes
   
         CALL H5DCREATE_F(LEVEL2_GROUP_ID,"/MOPITT/Level2/Longitude",
     &        H5T_IEEE_F32LE, SPACE_RAW_1D, LON_RAW_ID, HDF_ERR)
         
         CALL H5DWRITE_F(LON_RAW_ID, H5T_NATIVE_REAL, 
     &        REAL(FLEX_LON%DATA(1:FLEX_LON%CURRENT_N),4), 
     &        DIM_RAW_1D, HDF_ERR)
         
         CALL WRITE_ATTRIBUTES(LON_RAW_ID,"Longitude", "degrees")
         
         CALL H5DCLOSE_F(LON_RAW_ID,HDF_ERR)
   
         ! write raw latitudes
   
         CALL H5DCREATE_F(LEVEL2_GROUP_ID,"/MOPITT/Level2/Latitude",
     &        H5T_IEEE_F32LE, SPACE_RAW_1D, LAT_RAW_ID, HDF_ERR)
         
         CALL H5DWRITE_F(LAT_RAW_ID, H5T_NATIVE_REAL, 
     &        REAL(FLEX_LAT%DATA(1:FLEX_LAT%CURRENT_N),4), 
     &        DIM_RAW_1D, HDF_ERR)
         
         CALL WRITE_ATTRIBUTES(LAT_RAW_ID,"Latitude", "degrees")
         
         CALL H5DCLOSE_F(LAT_RAW_ID,HDF_ERR)
   
         ! write raw times
   
         CALL H5DCREATE_F(LEVEL2_GROUP_ID,"/MOPITT/Level2/Time",
     &        H5T_IEEE_F64LE, SPACE_RAW_1D, TIME_RAW_ID, HDF_ERR)
         
         CALL H5DWRITE_F(TIME_RAW_ID, H5T_NATIVE_DOUBLE, 
     &        FLEX_TIME%DATA(1:FLEX_TIME%CURRENT_N), 
     &        DIM_RAW_1D, HDF_ERR)
         
         CALL WRITE_ATTRIBUTES(TIME_RAW_ID,"Time","???")
         
         CALL H5DCLOSE_F(TIME_RAW_ID,HDF_ERR)
   
         ! create dataspace for raw 2D diagnostics
   
         CALL H5SCREATE_SIMPLE_F(2,DIM_RAW_2D,SPACE_RAW_2D,HDF_ERR)
   
         ! write raw MOPITT CO profiles
   
         CALL H5DCREATE_F(LEVEL2_GROUP_ID,"/MOPITT/Level2/MOP_CO",
     &        H5T_IEEE_F32LE, SPACE_RAW_2D, DSET_MOP_CO_RAW_ID, HDF_ERR)
         
         CALL H5DWRITE_F(DSET_MOP_CO_RAW_ID, H5T_NATIVE_REAL, 
     &        REAL(FLEX_MOP_CO%DATA(:,1:FLEX_TIME%CURRENT_N)*1e9,4), 
     &        DIM_RAW_2D, HDF_ERR)
         
         CALL WRITE_ATTRIBUTES(DSET_MOP_CO_RAW_ID,"MOPITT CO profiles",
     &        "ppbv")
         
         CALL H5DCLOSE_F(DSET_MOP_CO_RAW_ID,HDF_ERR)
   
         ! write raw GC CO profiles as observed by GEOS-Chem
   
         CALL H5DCREATE_F(LEVEL2_GROUP_ID,"/MOPITT/Level2/MOP_GC_CO",
     &        H5T_IEEE_F32LE, SPACE_RAW_2D, DSET_GC_CO_RAW_ID, HDF_ERR)
         
         CALL H5DWRITE_F(DSET_GC_CO_RAW_ID, H5T_NATIVE_REAL, 
     &        REAL(FLEX_GC_CO%DATA(:,1:FLEX_TIME%CURRENT_N)*1e9,4), 
     &        DIM_RAW_2D, HDF_ERR)
         
         CALL WRITE_ATTRIBUTES(DSET_GC_CO_RAW_ID,
     &       "GEOS-Chem CO profiles in MOPITT observation space","ppbv")
         
         CALL H5DCLOSE_F(DSET_GC_CO_RAW_ID,HDF_ERR)
   
         !close file
   
         CALL H5SCLOSE_F(SPACE_MOPITT,HDF_ERR)
   
         CALL H5GCLOSE_F(DATA_GROUP_ID, HDF_ERR)
         CALL H5GCLOSE_F(GRID_GROUP_ID, HDF_ERR)
         CALL H5GCLOSE_F(MOPITT_GROUP_ID, HDF_ERR)

      ENDIF
      !CALL H5FCLOSE_F(FILE_ID,HDF_ERR)

      ! close HDF5 interface

      CALL H5CLOSE_F(HDF_ERR)

      CALL H5EPRINT_F(HDF_ERR,"hdf_error")

      ! clear flexible arrays

      CALL CLEAR_FLEX_REAL_1D(FLEX_LON)
      CALL CLEAR_FLEX_REAL_1D(FLEX_LAT)
      CALL CLEAR_FLEX_REAL_1D(FLEX_TIME)

      CALL CLEAR_FLEX_REAL_2D(FLEX_MOP_CO)
      CALL CLEAR_FLEX_REAL_2D(FLEX_GC_CO)

      END SUBROUTINE MAKE_MOPITT_BIAS_FILE_HDF5


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
      CALL H5TSET_SIZE_F(ATYPE_ID,80,HDF_ERR)

      CALL H5ACREATE_F(DSET_ID,"Long name",
     &     ATYPE_ID,ASPACE_ID,ATT_ID,HDF_ERR)
      CALL H5AWRITE_F(ATT_ID,ATYPE_ID,LONGNAME,
     &     ADIMS,HDF_ERR)

      CALL H5ACLOSE_F(ATT_ID,HDF_ERR)
      CALL H5SCLOSE_F(ASPACE_ID,HDF_ERR)

      ! create attribute "Unit"

      CALL H5SCREATE_SIMPLE_F(1,ADIMS,ASPACE_ID,HDF_ERR)

      CALL H5TCOPY_F(H5T_NATIVE_CHARACTER,ATYPE_ID,HDF_ERR)
      CALL H5TSET_SIZE_F(ATYPE_ID,80,HDF_ERR)

      CALL H5ACREATE_F(DSET_ID,"Unit",
     &     ATYPE_ID,ASPACE_ID,ATT_ID,HDF_ERR)
      CALL H5AWRITE_F(ATT_ID,ATYPE_ID,UNIT,
     &     ADIMS,HDF_ERR)

      CALL H5ACLOSE_F(ATT_ID,HDF_ERR)
      CALL H5SCLOSE_F(ASPACE_ID,HDF_ERR)

      END SUBROUTINE WRITE_ATTRIBUTES

! mkeller: helper routines for managing flexible arrays
!          reinventing the wheel here, but hey...

      SUBROUTINE INIT_FLEX_REAL_1D(INPUT)
      
      TYPE(FLEX_REAL_1D):: INPUT
      INPUT%CURRENT_N = 0
      INPUT%MAX_N = 1000
      IF(ALLOCATED(INPUT%DATA)) DEALLOCATE(INPUT%DATA) ! safety first
      ALLOCATE(INPUT%DATA(INPUT%MAX_N))
      
      END SUBROUTINE INIT_FLEX_REAL_1D
      
      SUBROUTINE GROW_FLEX_REAL_1D(INPUT)
      
      TYPE(FLEX_REAL_1D) :: INPUT
      REAL*8, ALLOCATABLE :: TEMP_ARRAY(:)
      ALLOCATE(TEMP_ARRAY(INPUT%MAX_N * 2))
      TEMP_ARRAY(1:INPUT%MAX_N) = INPUT%DATA
      DEALLOCATE(INPUT%DATA)
      ALLOCATE(INPUT%DATA(INPUT%MAX_N * 2))
      INPUT%DATA = TEMP_ARRAY
      DEALLOCATE(TEMP_ARRAY)
      INPUT%MAX_N = INPUT%MAX_N * 2
      
      END SUBROUTINE GROW_FLEX_REAL_1D
      
      SUBROUTINE PUSH_FLEX_REAL_1D(INPUT, NEW_VAL)
      
      TYPE(FLEX_REAL_1D) :: INPUT
      REAL*8 :: NEW_VAL
      IF(INPUT%CURRENT_N == INPUT%MAX_N) THEN
         CALL GROW_FLEX_REAL_1D(INPUT)
      ENDIF
      INPUT%CURRENT_N = INPUT%CURRENT_N + 1
      INPUT%DATA(INPUT%CURRENT_N) = NEW_VAL
      
      END SUBROUTINE PUSH_FLEX_REAL_1D
      
      SUBROUTINE CLEAR_FLEX_REAL_1D(INPUT)
      
      TYPE(FLEX_REAL_1D) :: INPUT
      IF(ALLOCATED(INPUT%DATA)) DEALLOCATE(INPUT%DATA)
      
      END SUBROUTINE CLEAR_FLEX_REAL_1D

!--------------------------------------------------------------------------------
      
      SUBROUTINE INIT_FLEX_REAL_2D(INPUT)
      
      TYPE(FLEX_REAL_2D):: INPUT
      INPUT%CURRENT_N = 0
      INPUT%MAX_N = 1000
      IF(ALLOCATED(INPUT%DATA)) DEALLOCATE(INPUT%DATA) ! safety first
      ALLOCATE(INPUT%DATA(MAXLEV,INPUT%MAX_N))
      
      END SUBROUTINE INIT_FLEX_REAL_2D
      
      SUBROUTINE GROW_FLEX_REAL_2D(INPUT)
      
      TYPE(FLEX_REAL_2D) :: INPUT
      REAL*8, ALLOCATABLE :: TEMP_ARRAY(:,:)
      ALLOCATE(TEMP_ARRAY(MAXLEV,INPUT%MAX_N * 2))
      TEMP_ARRAY(:,1:INPUT%MAX_N) = INPUT%DATA
      DEALLOCATE(INPUT%DATA)
      ALLOCATE(INPUT%DATA(MAXLEV,INPUT%MAX_N * 2))
      INPUT%DATA = TEMP_ARRAY
      DEALLOCATE(TEMP_ARRAY)
      INPUT%MAX_N = INPUT%MAX_N * 2
      
      END SUBROUTINE GROW_FLEX_REAL_2D
      
      SUBROUTINE PUSH_FLEX_REAL_2D(INPUT, NEW_VAL, NLEV)
      
      TYPE(FLEX_REAL_2D) :: INPUT
      REAL*8 :: NEW_VAL(MAXLEV)
      INTEGER :: NLEV
      IF(INPUT%CURRENT_N == INPUT%MAX_N) THEN
         CALL GROW_FLEX_REAL_2D(INPUT)
      ENDIF
      INPUT%CURRENT_N = INPUT%CURRENT_N + 1
      INPUT%DATA(MAXLEV-NLEV+1:MAXLEV,INPUT%CURRENT_N) = NEW_VAL(1:NLEV)
      
      END SUBROUTINE PUSH_FLEX_REAL_2D
      
      SUBROUTINE CLEAR_FLEX_REAL_2D(INPUT)
      
      TYPE(FLEX_REAL_2D) :: INPUT
      IF(ALLOCATED(INPUT%DATA)) DEALLOCATE(INPUT%DATA)
      
      END SUBROUTINE CLEAR_FLEX_REAL_2D

      END MODULE MOPITT_OBS_MOD
