      MODULE WEAK_CONSTRAINT_MOD

!!
!! ********************************************************************************
!!
!! Module WEAK_CONSTRAINT contains subroutines and arrays needed for weak-constraint
!! 4D-VAR.
!!
!! mkeller, on a cold september morning
!!
!!
!! Module Variables:
!! ********************************************************************************
!!
!! ( 1) FORCE_U_FULLGRID (REAL*8): values of estimated forcing on full model grid
!! ( 2) FORCE_U_SUBGRID (REAL*8): values of forcing on subgrid
!! ( 3) GRADNT_U (REAL*8): gradient of cost function with respect to forcing terms
!! ( 4) STT_ADJ_SUBGRID (REAL*8): values of adjoint variables on forcing subgrid
!! ( 5) X_U (REAL*8): values of forcing vector for optimization routine
!! ( 6) X_GRADNT_U (REAL*8): gradient of cost function for optimization routine
!! ( 7) MIN_LON_U (REAL*8): minimal longitude of forcing subgrid
!! ( 8) MAX_LON_U (REAL*8): maximal longitude of forcing subgrid
!! ( 9) MIN_LAT_U (REAL*8): minimal latitude of forcing subgrid
!! (10) MAX_LAT_U (REAL*8): maximal latitude of forcing subgrid
!! (11) SIGMA_U (REAL*8): standard deviation of forcing covariance matrix
!! (12) MIN_LON_U_INDEX (INTEGER): index of smallest subgrid longitude on full grid
!! (13) MAX_LON_U_INDEX (INTEGER): index of largest subgrid longitude on full grid
!! (14) MIN_LAT_U_INDEX (INTEGER): index of smallest subgrid latitude on full grid
!! (15) MAX_LAT_U_INDEX (INTEGER): index of largest subgrid latitude on full grid
!! (16) MIN_LEV_U_INDEX (INTEGER): index of smallest subgrid vertical level on full grid
!! (17) MAX_LEV_U_INDEX (INTEGER): index of largest subgrid vertical level on full grid
!! (18) LEN_LON_U (INTEGER): number of longitudes on subgrid
!! (19) LEN_LAT_U (INTEGER): number of latitudes on subgrid
!! (20) LEN_LEV_U (INTEGER): number of vertical levels on subgrid
!! (21) N_TIMESTEPS_U (INTEGER): number of forcing timesteps
!! (22) LEN_SUBWINDOW_U (INTEGER): number of subwindow timesteps (i.e. timesteps with the same forcing)
!! (23) CT_SUB_U (INTEGER): timestep in subwindow
!! (24) CT_MAIN_U (INTEGER): main forcing timestep
!!
!! Module Routines:
!! ********************************************************************************
!!
!! ( 1) INITIALIZE_GRID_INDICES_U : initializes grid indices
!! ( 2) INIT_SETULB_U: allocates & zeroes arrays for optimization
!! ( 3) CLEAN_SETULB_U: deallocates arrays for optimization
!! ( 4) INIT_WEAK_CONSTRAINT: allocates & zeroes module arrays
!! ( 5) CLEAN_WEAK_CONSTRAINT: deallocates module arrays
!! ( 6) FORCE_SUBGRID_TO_FULLGRID_U: maps the estimated force from the subgrid to the full grid
!! ( 7) STT_ADJ_FULLGRID_TO_SUBGRID: maps the adjoint variables from the full grid to the subgrid
!! ( 8) SET_CT_U: increments subwindow timestep
!! ( 9) SET_CT_MAIN_U: increments main forcing timestep
!! (10) ITS_TIME_FOR_U: function that checks whether its time to update the forcing estimate
!! (11) CALC_GRADNT_U: update gradient of cost function with respect to u and call optimization routines if necessary
!! (12) GET_X_U_FROM_FORCE_U: get X_U from FORCE_U_SUBGRID
!! (12) GET_FORCE_U_FROM_X_U: get FORCE_U_SUBGRID from X_U
!! (13) WRITE_WA_U: writes the WA_U array (needed for optimization) to disk
!! (14) READ_WA_U: reads the WA_U array (needed for optimization) from disk
!! (15) MAKE_FORCE_U_FILE: write forcing values to disk
!! (16) READ_FORCE_U_FILE: read forcing values from disk


      IMPLICIT NONE

  !! Define arrays needed for the computation of model forcing terms

      PUBLIC

      REAL*8, ALLOCATABLE :: FORCE_U_FULLGRID(:,:,:) !! Model forcing on full grid
      REAL*8, ALLOCATABLE :: X_U(:,:) !! forcing array for optimization routines
      REAL*8, ALLOCATABLE :: X_U_TEMP(:,:) !! temporary array needed during the computation of the gradient of the cost function
      REAL*8, ALLOCATABLE :: X_GRADNT_U(:,:) !! forcing gradient array for optimization routines

      REAL*8 :: MIN_LON_U
      REAL*8 :: MAX_LON_U
      REAL*8 :: MIN_LAT_U
      REAL*8 :: MAX_LAT_U
      REAL*8 :: SCALE_FACTOR_U
      REAL*8 :: SIGMA_U

      INTEGER :: N_TRACER_U
      INTEGER :: MIN_LON_U_INDEX
      INTEGER :: MAX_LON_U_INDEX
      INTEGER :: MIN_LAT_U_INDEX
      INTEGER :: MAX_LAT_U_INDEX
      INTEGER :: MIN_LEV_U_INDEX
      INTEGER :: MAX_LEV_U_INDEX
      INTEGER :: LEN_LON_U
      INTEGER :: LEN_LAT_U
      INTEGER :: LEN_LEV_U
      INTEGER :: N_TIMESTEPS_U
      INTEGER :: LEN_SUBWINDOW_U
      INTEGER :: CT_SUB_U
      INTEGER :: CT_MAIN_U

      LOGICAL :: DO_WEAK_CONSTRAINT
      LOGICAL :: PERTURB_OBS_U
      INTEGER :: NOPT_U         !! total number of gridpoints where forcing terms are to be estimated

      REAL*8 :: MOP_COST, FORCE_COST, BG_COST

      CONTAINS

      SUBROUTINE INITIALIZE_GRID_INDICES_U

      USE GRID_MOD, ONLY : GET_XMID, GET_YMID

                                ! local variables
      INTEGER                    :: I
      INTEGER                    :: J
#     include "CMN_SIZE"        ! Size params

                                !! Calculate longitude array indices for forcing subgrid

      IF ( GET_XMID(1) >= MIN_LON_U ) MIN_LON_U_INDEX = 1

      DO I = 2, IIPAR

      IF ( ( GET_XMID(I-1) < MIN_LON_U) .AND.
     & ( GET_XMID(I) >= MIN_LON_U ) ) MIN_LON_U_INDEX = I

      IF ( ( GET_XMID(I-1) <= MAX_LON_U) .AND.
     & ( GET_XMID(I) > MAX_LON_U ) ) MAX_LON_U_INDEX = I-1


      ENDDO

      IF ( GET_XMID(IIPAR) <= MAX_LON_U ) MAX_LON_U_INDEX = IIPAR

      LEN_LON_U = MAX_LON_U_INDEX - MIN_LON_U_INDEX + 1

                                !! Calculate latitude array indices for forcing subgrid

      IF ( GET_YMID(1) >= MIN_LAT_U ) MIN_LAT_U_INDEX = 1

      DO J = 2, JJPAR

      IF ( ( GET_YMID(J-1) < MIN_LAT_U) .AND.
     &( GET_YMID(J) >= MIN_LAT_U ) ) MIN_LAT_U_INDEX = J
      IF ( ( GET_YMID(J-1) <= MAX_LAT_U) .AND.
     &( GET_YMID(J) > MAX_LAT_U ) ) MAX_LAT_U_INDEX = J-1

      ENDDO

      IF ( GET_YMID(JJPAR) <= MAX_LAT_U ) MAX_LAT_U_INDEX = JJPAR

      LEN_LAT_U = MAX_LAT_U_INDEX - MIN_LAT_U_INDEX + 1

      LEN_LEV_U = MAX_LEV_U_INDEX - MIN_LEV_U_INDEX + 1

                                !! Calculate total number of forcing gridpoints to be optimized

      NOPT_U = LEN_LON_U * LEN_LAT_U * LEN_LEV_U

      PRINT *,'WEAK_CONSTRAINT: MIN_LON_U' , MIN_LON_U,
     &     GET_XMID(MIN_LON_U_INDEX)
      PRINT *,'WEAK_CONSTRAINT: MIN_LAT_U' , GET_YMID(MIN_LAT_U_INDEX)
      PRINT *,'WEAK_CONSTRAINT: MAX_LON_U' , MAX_LON_U,
     &     GET_XMID(MAX_LON_U_INDEX)
      PRINT *,'WEAK_CONSTRAINT: MAX_LAT_U' , GET_YMID(MAX_LAT_U_INDEX)


      END SUBROUTINE INITIALIZE_GRID_INDICES_U

      SUBROUTINE INIT_WEAK_CONSTRAINT
                                !!
                                !! **********************************************************
                                !! Subroutine INIT_WEAK_CONSTRAINT initializes and zeroes all module arrays.
                                !! mkeller (15/09/2011)
                                !! **********************************************************
                                !!

      USE ERROR_MOD, ONLY : ALLOC_ERR
      USE TIME_MOD, ONLY : CALC_RUN_DAYS
      USE TIME_MOD, ONLY : GET_TS_DYN
      USE TRACER_MOD,      ONLY : N_TRACERS

#       include "CMN_SIZE"      ! Size parameters

      INTEGER :: AS, DAYS, TS_DYN

      SCALE_FACTOR_U = 0.5 * 10.0 ** (-10)
      SIGMA_U = 0.5 * 10.0**(-10)

                   !! Compute grid indices
      CALL INITIALIZE_GRID_INDICES_U

      DAYS = CALC_RUN_DAYS()
      TS_DYN = GET_TS_DYN()

                                !! Allow for the force to be estimated at a subset of timesteps
      N_TIMESTEPS_U = (DAYS * 24 * 60 / TS_DYN) / (LEN_SUBWINDOW_U) + 1

                                !! Initialize time counters. Note that in the adjoint calculation, time runs backwards
      CT_SUB_U = -1
      CT_MAIN_U = 1
                                !! Allocate all required arrays
      ALLOCATE( FORCE_U_FULLGRID(IIPAR,JJPAR,LLPAR)
     &, STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'FORCE_U_FULLGRID' )
      FORCE_U_FULLGRID = 0d0

      ALLOCATE( X_U ( NOPT_U, N_TIMESTEPS_U ) )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'X_U' )
      X_U = 0d0

      ALLOCATE( X_U_TEMP ( NOPT_U, N_TIMESTEPS_U ) )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'X_U_TEMP' )
      X_U = 0d0

      ALLOCATE( X_GRADNT_U ( NOPT_U, N_TIMESTEPS_U ) )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'X_GRADNT_U' )
      X_GRADNT_U = 0d0

      PRINT *,'WEAK_CONSTRAINT: NOPT_U: ',NOPT_U

      FORCE_COST = 0d0

      END SUBROUTINE INIT_WEAK_CONSTRAINT

      SUBROUTINE CLEAN_WEAK_CONSTRAINT

                                !!
                                !! **********************************************************
                                !! Subroutine CLEAN_SETULB_U deallocates the arrays needed for weak-constraint 4D-VAR
                                !! mkeller (15/09/2011)
                                !! **********************************************************
                                !!

      IF ( ALLOCATED ( FORCE_U_FULLGRID ) )DEALLOCATE (FORCE_U_FULLGRID)
      IF ( ALLOCATED ( X_U              ) )DEALLOCATE ( X_U )
      IF ( ALLOCATED ( X_U_TEMP         ) )DEALLOCATE ( X_U_TEMP )
      IF ( ALLOCATED ( X_GRADNT_U       ) )DEALLOCATE ( X_GRADNT_U )


      END SUBROUTINE CLEAN_WEAK_CONSTRAINT


      SUBROUTINE SET_CT_U ( INCREASE, RESET , FLIP )

                                !!
                                !! **********************************************************
                                !! Subroutine SET_CT_U increments the time counter for the forcing subwindow.
                                !! **********************************************************
                                !!

      LOGICAL, INTENT(IN), OPTIONAL :: INCREASE
      LOGICAL, INTENT(IN), OPTIONAL :: RESET
      LOGICAL, INTENT(IN), OPTIONAL :: FLIP

      IF ( PRESENT ( RESET ) ) THEN
         CT_SUB_U = 0
      ENDIF

      IF( PRESENT ( FLIP ) ) THEN
         PRINT *,'WEAK_CONSTRAINT: BEFORE',CT_SUB_U
         IF(CT_SUB_U == 0) THEN
            CT_SUB_U=0
         ELSE
            CT_SUB_U = LEN_SUBWINDOW_U - CT_SUB_U
         ENDIF
         PRINT *,'WEAK_CONSTRAINT: AFTER',CT_SUB_U
      ENDIF

      IF( PRESENT ( INCREASE ) ) THEN
         IF ( INCREASE) THEN
            CT_SUB_U = CT_SUB_U + 1
         ELSE
            CT_SUB_U = CT_SUB_U - 1
         ENDIF
      ENDIF
                                !! return to calling program

      END SUBROUTINE SET_CT_U

      SUBROUTINE SET_CT_MAIN_U( INCREASE, RESET )

                                !!
                                !! **********************************************************
                                !! Subroutine SET_CT_MAIN_U increases/decreases the time counter for the main forcing window.
                                !! **********************************************************
                                !!


      LOGICAL, INTENT(IN), OPTIONAL :: INCREASE
      LOGICAL, INTENT(IN), OPTIONAL :: RESET

      IF ( PRESENT ( RESET ) ) THEN
         CT_MAIN_U = 1
         RETURN
      END IF
      IF ( PRESENT ( INCREASE ) ) THEN
         IF ( INCREASE ) THEN
            CT_MAIN_U = CT_MAIN_U + 1
         ELSE
            CT_MAIN_U = CT_MAIN_U - 1
         END IF
      ELSE
         CT_MAIN_U = CT_MAIN_U + 1
      END IF

                                !! return to calling program

      END SUBROUTINE SET_CT_MAIN_U

      FUNCTION ITS_TIME_FOR_U () RESULT ( FLAG )

                                !!
                                !! **********************************************************
                                !! Subroutine ITS_TIME_FOR_U returns true if it's time to estimate the new forcing terms
                                !! **********************************************************
                                !!

      LOGICAL :: FLAG

      FLAG = .FALSE.

      IF (CT_SUB_U == LEN_SUBWINDOW_U) FLAG = .TRUE.

                                !! return to calling program

      END FUNCTION ITS_TIME_FOR_U

      SUBROUTINE CALC_GRADNT_U(YYYYMMDD,HHMMSS)

                                !!
                                !! **********************************************************
                                !! Subroutine CALC_GRADNT_U updates the gradient of the cost function with respect to u
                                !! and calculates the new estimates once all the required data have been gathered
                                !! **********************************************************
                                !!

      USE ADJ_ARRAYS_MOD,       ONLY : N_CALC, COST_FUNC

      INTEGER, INTENT(IN) :: YYYYMMDD, HHMMSS

      X_GRADNT_U(:,CT_MAIN_U) = X_GRADNT_U(:,CT_MAIN_U) + GET_STT_X_U()

      PRINT *,'WEAK_CONSTRAINT: IN CALC_GRADNT_U'
      PRINT *,'WEAK_CONSTARAINT:',YYYYMMDD,HHMMSS
      PRINT *,'WEAK_CONSTRAINT:',CT_SUB_U,CT_MAIN_U

      IF ( ITS_TIME_FOR_U() ) THEN

      CALL GET_FORCE_U_FROM_X_U

      CALL COMPUTE_APRIORI_U

      X_GRADNT_U(:,CT_MAIN_U) = X_GRADNT_U(:,CT_MAIN_U) +
     &     X_U_TEMP(:,CT_MAIN_U)

      X_GRADNT_U(:,CT_MAIN_U) = X_GRADNT_U(:,CT_MAIN_U) * 2.0

      CALL SET_CT_U ( RESET = .TRUE. )

      CALL SET_CT_MAIN_U ( INCREASE = .FALSE. )

      ENDIF

                                !! return to calling program

      END SUBROUTINE CALC_GRADNT_U


      SUBROUTINE GET_X_U_FROM_FORCE_U

                                !!
                                !! **********************************************************
                                !! SUBROUTINE GET_FORCE_U FROM_X_U obtains the array X_U from the current forcing estimate
                                !! **********************************************************
                                !!

      USE TRACER_MOD, ONLY : N_TRACERS

                                ! Local variables
      INTEGER                    :: I, J, L, M, N
      INTEGER                    :: I_DUM

      !DO N = 1, N_TRACERS
      DO L = MIN_LEV_U_INDEX, MAX_LEV_U_INDEX
      DO J = MIN_LAT_U_INDEX, MAX_LAT_U_INDEX
      DO I = MIN_LON_U_INDEX, MAX_LON_U_INDEX

      I_DUM    = (I - MIN_LON_U_INDEX + 1) +
     &(  LEN_LON_U * ( J - MIN_LAT_U_INDEX )  )
     &+( LEN_LON_U * LEN_LAT_U * ( L - MIN_LEV_U_INDEX ))
      !&+( LEN_LON_U * LEN_LAT_U * LEN_LEV_U * ( N - 1 )  )

      X_U( I_DUM , CT_MAIN_U) = FORCE_U_FULLGRID( I, J, L)!, N )

      ENDDO
      ENDDO
      ENDDO
      !ENDDO

                                !! return to calling program

      END SUBROUTINE GET_X_U_FROM_FORCE_U

      SUBROUTINE COMPUTE_APRIORI_U

      !!
      !!****************************************************************************
      !! SUBROUTINE COMPUTE_APRIORI_U computes the apriori-term (i.e. the Q-matrix term)
      !! to the gradient of the cost function with respect to u_n
      !!****************************************************************************
      !!

      USE TRACER_MOD, ONLY : N_TRACERS
      USE ADJ_ARRAYS_MOD,       ONLY : N_CALC, COST_FUNC

                                ! Local variables
      INTEGER                    :: I, J, L, M, N
      INTEGER                    :: I_DUM
      REAL*8                     :: COST_TEMP

      PRINT *,"MKDB: COST BEFORE FORCE: ",COST_FUNC
      COST_TEMP = COST_FUNC

      PRINT *,"MKDB: MAX FORCE", MAXVAL(FORCE_U_FULLGRID)

!      DO N = 1, N_TRACERS
      DO L = MIN_LEV_U_INDEX, MAX_LEV_U_INDEX
      DO J = MIN_LAT_U_INDEX, MAX_LAT_U_INDEX
      DO I = MIN_LON_U_INDEX, MAX_LON_U_INDEX

      I_DUM    = (I - MIN_LON_U_INDEX + 1) +
     &(  LEN_LON_U * ( J - MIN_LAT_U_INDEX )  )
     &+( LEN_LON_U * LEN_LAT_U * ( L - MIN_LEV_U_INDEX ))
!     &+( LEN_LON_U * LEN_LAT_U * LEN_LEV_U * ( N - 1 )  )

      IF ( L< ( MIN_LEV_U_INDEX + 3 ) ) THEN
         X_U_TEMP( I_DUM , CT_MAIN_U) = FORCE_U_FULLGRID( I, J, L ) /
     &       ( (SIGMA_U*1/(1+MIN_LEV_U_INDEX + 3-L) )**2)

         COST_FUNC = COST_FUNC +  FORCE_U_FULLGRID( I, J, L)**2 /
     &        ( (SIGMA_U *1/(1+MAX_LEV_U_INDEX + 3 - L) )**2)

      ELSE IF ( L > ( MAX_LEV_U_INDEX - 3 ) ) THEN
         X_U_TEMP( I_DUM , CT_MAIN_U) = FORCE_U_FULLGRID( I, J, L ) /
     &       ( (SIGMA_U*1/(1+MAX_LEV_U_INDEX - 3 + L) )**2)

         COST_FUNC = COST_FUNC +  FORCE_U_FULLGRID( I, J, L )**2 /
     &        ( (SIGMA_U *1/(1+MAX_LEV_U_INDEX - 3 + L) )**2)

      ELSE
         X_U_TEMP( I_DUM , CT_MAIN_U) = FORCE_U_FULLGRID( I, J, L ) /
     &       ( SIGMA_U **2)

         COST_FUNC = COST_FUNC +  FORCE_U_FULLGRID( I, J, L )**2 /
     &       ( SIGMA_U **2)

      ENDIF



      ENDDO
      ENDDO
      ENDDO
!      ENDDO

      FORCE_COST = COST_FUNC - COST_TEMP + FORCE_COST

                                !! return to calling program

      END SUBROUTINE COMPUTE_APRIORI_U



      FUNCTION GET_STT_X_U() RESULT ( STT_X_U )

                                !!
                                !! **********************************************************
                                !! FUNCTION GET_STT_X_U reads the value of STT_ADJ into a 1D-array
                                !! **********************************************************
                                !!

      USE TRACER_MOD, ONLY     : N_TRACERS
      USE ADJ_ARRAYS_MOD, ONLY : STT_ADJ

      ! Local variables
      INTEGER                    :: I, J, L, M, N
      INTEGER                    :: I_DUM
      REAL*8                     :: STT_X_U ( NOPT_U )

      PRINT *,'WEAK_CONSTRAINT: MAX_STT_ADJ',MAXVAL(STT_ADJ)

!      DO N = 1, N_TRACERS
      DO L = MIN_LEV_U_INDEX, MAX_LEV_U_INDEX
      DO J = MIN_LAT_U_INDEX, MAX_LAT_U_INDEX
      DO I = MIN_LON_U_INDEX, MAX_LON_U_INDEX

      I_DUM    = (I - MIN_LON_U_INDEX + 1) +
     &(  LEN_LON_U * ( J - MIN_LAT_U_INDEX )  )
     &+( LEN_LON_U * LEN_LAT_U * ( L - MIN_LEV_U_INDEX ))
!     &+( LEN_LON_U * LEN_LAT_U * LEN_LEV_U * ( N - 1 )  )

         !     I_DUM    = (I - MIN_LON_U_INDEX + 1) + (  LEN_LON_U * ( J - MIN_LAT_U_INDEX )  ) &
         !     + (  LEN_LON_U * LEN_LAT_U * ( L - MIN_LEV_U_INDEX )  ) &
         !     + (  LEN_LON_U * LEN_LAT_U * LEN_LEV_U * ( N - 1 )  )

      STT_X_U( I_DUM ) = STT_ADJ( I, J, L, N_TRACER_U )

      ENDDO
      ENDDO
      ENDDO
!      ENDDO

                                !! return to calling program

      END FUNCTION GET_STT_X_U


      SUBROUTINE MAKE_GDT_U_FILE( )
!
!******************************************************************************
!  Subroutine MAKE_GDT_FILE creates a binary file of forcing gradients
!  mkeller
!
!  Module Variable as Input:
!  ============================================================================
!  (1 ) N_CALC     : Current iteration number
!******************************************************************************
!
      ! References to F90 modules

      USE ADJ_ARRAYS_MOD,      ONLY : N_CALC
      USE ADJ_ARRAYS_MOD,      ONLY : EXPAND_NAME
      USE DIRECTORY_ADJ_MOD,   ONLY : OPTDATA_DIR
      USE ERROR_MOD,           ONLY : DEBUG_MSG, ERROR_STOP
      USE FILE_MOD,            ONLY : IU_RST,      IOERROR
      USE TIME_MOD,            ONLY : EXPAND_DATE, GET_TAU
      USE TRACER_MOD,          ONLY : N_TRACERS

      ! Local Variables

      CHARACTER(LEN=255)   :: FILENAME

      INTEGER :: I, J

      CHARACTER(LEN=20)    :: OUTPUT_GDT_FILE
      CHARACTER(LEN=40)    :: UNIT



      ! Hardwire output file for now
      OUTPUT_GDT_FILE = 'gctm.gdt.forcing.NN'

      ! Copy the output observation file name into a local variable
      FILENAME = TRIM( OUTPUT_GDT_FILE )

      ! Append the iteration number suffix to the file name
      CALL EXPAND_NAME( FILENAME, N_CALC )

      ! Add the OPTDATA_DIR prefix to the file name
      FILENAME = TRIM( OPTDATA_DIR ) //  TRIM( FILENAME )

      OPEN( UNIT = IU_RST, FILE = FILENAME, FORM = 'UNFORMATTED' )

      DO I=1,NOPT_U
         DO J=1,N_TIMESTEPS_U

            WRITE( IU_RST ) X_GRADNT_U(I,J)

         ENDDO
      ENDDO

      CLOSE ( IU_RST )

      ! Return to calling program
      END SUBROUTINE MAKE_GDT_U_FILE

      SUBROUTINE READ_GDT_U_FILE( )
!
!******************************************************************************
!  Subroutine MAKE_GDT_FILE creates a binary file of forcing gradients
!  mkeller
!
!  Module Variable as Input:
!  ============================================================================
!  (1 ) N_CALC     : Current iteration number
!******************************************************************************
!
      ! References to F90 modules

      USE ADJ_ARRAYS_MOD,      ONLY : N_CALC
      USE ADJ_ARRAYS_MOD,      ONLY : EXPAND_NAME
      USE DIRECTORY_ADJ_MOD,   ONLY : OPTDATA_DIR
      USE ERROR_MOD,           ONLY : DEBUG_MSG, ERROR_STOP
      USE FILE_MOD,            ONLY : IU_RST,      IOERROR
      USE TIME_MOD,            ONLY : EXPAND_DATE, GET_TAU
      USE TRACER_MOD,          ONLY : N_TRACERS

      ! Local Variables

      CHARACTER(LEN=255)   :: FILENAME

      INTEGER :: I, J

      CHARACTER(LEN=20)    :: OUTPUT_GDT_FILE
      CHARACTER(LEN=40)    :: UNIT


      ! Hardwire output file for now
      OUTPUT_GDT_FILE = 'gctm.gdt.forcing.NN'

      ! Copy the output observation file name into a local variable
      FILENAME = TRIM( OUTPUT_GDT_FILE )

      ! Append the iteration number suffix to the file name
      CALL EXPAND_NAME( FILENAME, N_CALC )

      ! Add the OPTDATA_DIR prefix to the file name
      FILENAME = TRIM( OPTDATA_DIR ) //  TRIM( FILENAME )

      OPEN( UNIT = IU_RST, FILE = FILENAME, FORM = 'UNFORMATTED' )

      DO I=1,NOPT_U
         DO J=1,N_TIMESTEPS_U

            READ( IU_RST ) X_GRADNT_U(I,J)

         ENDDO
      ENDDO

      CLOSE ( IU_RST )

      ! Return to calling program
      END SUBROUTINE READ_GDT_U_FILE

      SUBROUTINE GET_FORCE_U_FROM_X_U

                                !!
                                !! **********************************************************
                                !! SUBROUTINE GET_FORCE_U FROM_X_U obtains the current forcing estimate from the X_U array
                                !! **********************************************************
                                !!

      USE TRACER_MOD, ONLY : N_TRACERS
#     include "CMN_SIZE"      ! Size parameters

                                ! Local variables
      INTEGER                    :: I, J, L, M, N
      INTEGER                    :: I_DUM

      I_DUM = 0
!      DO N = 1, N_TRACERS
      DO L = MIN_LEV_U_INDEX, MAX_LEV_U_INDEX
      DO J = MIN_LAT_U_INDEX, MAX_LAT_U_INDEX
      DO I = MIN_LON_U_INDEX, MAX_LON_U_INDEX

      I_DUM    = (I - MIN_LON_U_INDEX + 1) +
     &(  LEN_LON_U * ( J - MIN_LAT_U_INDEX )  )
     &+( LEN_LON_U * LEN_LAT_U * ( L - MIN_LEV_U_INDEX ))
 !    &+( LEN_LON_U * LEN_LAT_U * LEN_LEV_U * ( N - 1 )  )

      FORCE_U_FULLGRID( I, J, L ) = X_U( I_DUM, CT_MAIN_U )

      ENDDO
      ENDDO
      ENDDO
!      ENDDO

                                !! return to calling program

      END SUBROUTINE GET_FORCE_U_FROM_X_U


      SUBROUTINE MAKE_FORCE_U_FILE( YYYYMMDD, HHMMSS )
      !
      !******************************************************************************
      !  Subroutine MAKE_ADJ_FILE creates a binary file of STT_ADJ
      !  (dkh, 10/03/04)exit

      !
      !  Arguments as Input:
      !  ============================================================================
      !  (1 ) YYYYMMDD : Year-Month-Date
      !  (2 ) HHMMSS   :  and Hour-Min-Sec for which to create an adjoint file
      !

      ! References to F90 modules
      USE ADJ_ARRAYS_MOD,    ONLY : STT_ADJ
      USE BPCH2_MOD
      USE DIRECTORY_ADJ_MOD, ONLY : DIAGADJ_DIR
      USE DIRECTORY_ADJ_MOD, ONLY : OPTDATA_DIR
      USE ERROR_MOD,         ONLY : DEBUG_MSG
      USE FILE_MOD,          ONLY : IU_RST,      IOERROR
      USE GRID_MOD,          ONLY : GET_XOFFSET, GET_YOFFSET
      USE LOGICAL_MOD,       ONLY : LPRT
      USE LOGICAL_ADJ_MOD,   ONLY : LTRAJ_SCALE
      USE TIME_MOD,          ONLY : EXPAND_DATE, GET_TAU
      USE TRACER_MOD,        ONLY : N_TRACERS
      USE TRACER_MOD,        ONLY : STT

#     include "CMN_SIZE"        ! Size parameters

                                ! Arguments
      INTEGER, INTENT(IN)     :: YYYYMMDD, HHMMSS

      ! Local Variables
      INTEGER                 :: I,    I0, IOS, J,  J0, L, N
      INTEGER                 :: YYYY, MM, DD,  HH, SS
      REAL*4                  :: TRACER(IIPAR,JJPAR,LLPAR)
      CHARACTER(LEN=255)      :: FILENAME

                                ! For binary punch file, version 2.0
      REAL*4                  :: LONRES, LATRES
      INTEGER, PARAMETER      :: HALFPOLAR = 1
      INTEGER, PARAMETER      :: CENTER180 = 1

      CHARACTER(LEN=40)       :: OUTPUT_FORCE_U_FILE
      CHARACTER(LEN=20)       :: MODELNAME
      CHARACTER(LEN=40)       :: CATEGORY
      CHARACTER(LEN=40)       :: UNIT
      CHARACTER(LEN=40)       :: RESERVED = ''
      CHARACTER(LEN=80)       :: TITLE

      ! Should make these user defined in input.gcadj
      !! Parameter
      INTEGER, PARAMETER      :: LLADJKEEP   = LLPAR
      !INTEGER, PARAMETER      :: NNADJKEEP   = N_TRACERS
      ! Now specify this input.gcadj
      !LOGICAL, PARAMETER      :: LTRAJ_SCALE = .TRUE.

      !=================================================================
      ! MAKE_FORCE_U_FILE begins here!
      !=================================================================


      ! Hardwire output file for now
      OUTPUT_FORCE_U_FILE = 'gctm.forcing.YYYYMMDD.hhmm'

      ! Define variables for BINARY PUNCH FILE OUTPUT
      TITLE    = 'GEOS-CHEM Forcing File: ' //
     &      'Instantaneous Forcing Values '
      CATEGORY = 'IJ-ADJ-$'
      LONRES   = DISIZE
      LATRES   = DJSIZE

      ! Call GET_MODELNAME to return the proper model name for
      ! the given met data being used (bmy, 6/22/00)
      MODELNAME = GET_MODELNAME()

      ! Get the nested-grid offsets
      I0 = GET_XOFFSET( GLOBAL=.TRUE. )
      J0 = GET_YOFFSET( GLOBAL=.TRUE. )

      !=================================================================
      ! Open the force file for output -- binary punch format
      !=================================================================

      ! Copy the output observation file name into a local variable
      FILENAME = TRIM( OUTPUT_FORCE_U_FILE )

      ! Append the iteration number suffix to the file name
      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )

      ! Add the ADJ_DIR prefix to the file name
      FILENAME = TRIM( OPTDATA_DIR ) //  TRIM( FILENAME )

      WRITE( 6, 100 ) TRIM( FILENAME )
100   FORMAT( '     - MAKE_ADJ_FILE: Writing ', a )

      ! Open checkpoint file for output
      CALL OPEN_BPCH2_FOR_WRITE( IU_RST, FILENAME, TITLE )

      !=================================================================
      ! Write forcing for each observed quantity to the forcing file
      !=================================================================
      !DO N = 1, N_TRACERS

      UNIT     = 'J'

      !Temporarily store quantities in the TRACER array


      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

      TRACER(I,J,L) = FORCE_U_FULLGRID(I,J,L)

      ENDDO
      ENDDO
      ENDDO


      CALL BPCH2( IU_RST,    MODELNAME, LONRES, LATRES,
     &         HALFPOLAR, CENTER180, CATEGORY,  1,
     &         UNIT,      GET_TAU(), GET_TAU(), RESERVED,
     &         IIPAR,     JJPAR,     LLPAR, I0+1,
     &         J0+1,      1,         TRACER  )

!      ENDDO

 !    TRACER = STT(:,:,:,1)

!      CALL BPCH2( IU_RST,    MODELNAME, LONRES, LATRES,
!     &         HALFPOLAR, CENTER180, CATEGORY,  N_TRACERS+1,
!     &         UNIT,      GET_TAU(), GET_TAU(), RESERVED,
!     &         IIPAR,     JJPAR,     LLPAR, I0+1,
!     &         J0+1,      1,         TRACER  )

      ! Close file
      CLOSE( IU_RST )

      ! Return to calling program

      END SUBROUTINE MAKE_FORCE_U_FILE

      SUBROUTINE READ_FORCE_U_FILE( YYYYMMDD, HHMMSS )
      !
      !******************************************************************************
      !  Subroutine READ_FORCE_U_FILE reads forcing values
      !  from a checkpoint file (binary punch file format)
      !  (dkh, 8/30/04)
      !
      !  Arguments as input:
      !  ============================================================================
      !  (1 ) YYYYMMDD : Year-Month-Day
      !  (2 ) HHMMSS   :  and Hour-Min-Sec for which to read restart file
      !
      !  Passed via CMN:
      !  ============================================================================
      !  (1 ) TAU    : TAU value (elapsed hours) at start of diagnostic interval
      !

      ! References to F90 modules
      USE BPCH2_MOD,         ONLY : OPEN_BPCH2_FOR_READ
      USE COMODE_MOD,        ONLY : CHK_CSPEC, JLOP
      USE DIRECTORY_ADJ_MOD, ONLY : ADJTMP_DIR
      USE DIRECTORY_ADJ_MOD, ONLY : DIAGADJ_DIR
      USE DIRECTORY_ADJ_MOD, ONLY : OPTDATA_DIR
      USE ERROR_MOD,         ONLY : DEBUG_MSG, ERROR_STOP
      USE FILE_MOD,          ONLY : IU_RST, IOERROR
      USE GCKPP_ADJ_GLOBAL,  ONLY : NTT
      USE LOGICAL_MOD,       ONLY : LCHEM , LSULF
      USE LOGICAL_MOD,       ONLY : LSOILNOX, LLIGHTNOX
      USE LOGICAL_MOD,       ONLY : LPRT
      USE LOGICAL_ADJ_MOD,   ONLY : LAERO_THERM
      USE LOGICAL_ADJ_MOD,   ONLY : LDEL_CHKPT
      USE RESTART_MOD,       ONLY : CHECK_DIMENSIONS
      USE TIME_MOD,          ONLY : EXPAND_DATE
      USE TRACER_MOD,        ONLY : ITS_A_FULLCHEM_SIM
      USE TRACER_MOD,        ONLY : ITS_AN_AEROSOL_SIM
      USE TRACER_MOD,        ONLY : N_TRACERS
      USE UNIX_CMDS_MOD,     ONLY : REMOVE_CMD

#     include "CMN_SIZE"   ! Size parameters
#     include "comode.h"   ! ITLOOP, IGAS
#     include "CMN_VEL"    ! XYLAI

      ! Arguments
      INTEGER, INTENT(IN) :: YYYYMMDD, HHMMSS

      ! Local Variables
      INTEGER             :: I, IOS, J, L, N, JLOOP, NN, NTL
      INTEGER             :: NCOUNT(NNPAR)
      REAL*4              :: TRACER(IIPAR,JJPAR,LLPAR)
      ! Remove these since we always recompute instead
      ! of checkpointing (dkh, 06/11/09)
      !      REAL*4		  :: CHECK1(IIPAR,JJPAR,LLPAR)
      !      REAL*4		  :: CHECK2(IIPAR,JJPAR,LLPAR)
      !      REAL*4		  :: CHECK3(IIPAR,JJPAR,LLPAR)
      !      REAL*4		  :: CHECK4(IIPAR,JJPAR,LLPAR)
      !      REAL*4		  :: CHECK5(IIPAR,JJPAR,LLPAR)
      !      REAL*4		  :: CHECK6(IIPAR,JJPAR,LLPAR)
      !      REAL*4		  :: CHECK7(IIPAR,JJPAR,LLPAR)
      !      REAL*4		  :: CHECK8(IIPAR,JJPAR,LLPAR)
      !      REAL*4		  :: CHECK9(IIPAR,JJPAR,LLPAR)
      REAL*4		  :: SMVGARRAY(ITLOOP,IGAS)

      !>>>
      ! Now include adjoint of F (dkh, 10/03/08)
      INTEGER             :: NS
      !<<<

      REAL*8              :: SUMTC
      CHARACTER(LEN=255)  :: FILENAME
      CHARACTER(LEN=255)  :: UNZIP_FILE_CMD
      CHARACTER(LEN=255)  :: REMOVE_CHK_FILE_CMD


      ! For binary punch file, version 2.0
      INTEGER             :: NI,     NJ,     NL,    NV
      INTEGER             :: IJLOOP
      INTEGER             :: IFIRST, JFIRST, LFIRST
      INTEGER             :: NTRACER,   NSKIP
      INTEGER             :: HALFPOLAR, CENTER180
      REAL*4              :: LONRES,    LATRES
      REAL*8              :: ZTAU0,     ZTAU1
      CHARACTER(LEN=20)   :: MODELNAME
      CHARACTER(LEN=40)   :: CATEGORY
      CHARACTER(LEN=40)   :: UNIT
      CHARACTER(LEN=40)   :: RESERVED

      ! added by Martin Keller

      CHARACTER(LEN=255)  :: INPUT_FILE
      LOGICAL             :: EX

      !=================================================================
      ! READ_FORCE_U_FILE begins here!
      !=================================================================

      ! Hardwire output file for now
      INPUT_FILE = 'gctm.forcing.YYYYMMDD.hhmm'

      ! Initialize some variables
      NCOUNT(:)     = 0
      TRACER(:,:,:) = 0e0

      !=================================================================
      ! Open checkpoint file and read top-of-file header
      !=================================================================

      ! Copy input file name to a local variable
      FILENAME = TRIM( INPUT_FILE )

      ! Replace YYYY, MM, DD, HH tokens in FILENAME w/ actual values
      CALL EXPAND_DATE( FILENAME, YYYYMMDD, HHMMSS )

      ! Add ADJ_DIR prefix to name
      FILENAME = TRIM( OPTDATA_DIR ) // TRIM( FILENAME )

      ! Echo some input to the screen
      WRITE( 6, '(a)'   ) REPEAT( '=', 79 )
      WRITE( 6, '(a,/)' ) 'F O R C I N G   F I L E   I N P U T'

      WRITE( 6, 100 ) TRIM( FILENAME )
100   FORMAT( '     - READ_FORCE_U_FILE: Reading ', a )

      ! MAKE SURE FILE EXISTS

      INQUIRE( FILE = FILENAME, EXIST=EX)

      IF( .NOT. EX ) THEN
         CALL MAKE_FORCE_U_FILE(YYYYMMDD,HHMMSS)
         RETURN
      ENDIF

      ! Open the binary punch file for input
      CALL OPEN_BPCH2_FOR_READ( IU_RST, FILENAME )

 !     DO N = 1, N_TRACERS
      READ( IU_RST, IOSTAT=IOS )
     &        MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180

                                ! IOS < 0 is end-of-file, so exit
!      IF ( IOS < 0 ) EXIT

         ! IOS > 0 is a real I/O error -- print error message
      IF ( IOS > 0 ) CALL IOERROR( IOS,IU_RST,'read_checkpt_file:7' )

      READ( IU_RST, IOSTAT=IOS )
     & CATEGORY, NTRACER,  UNIT, ZTAU0,  ZTAU1,  RESERVED,
     & NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
     &         NSKIP

      IF ( IOS /= 0 ) CALL IOERROR(IOS,IU_RST,'read_checkpt_file:8' )

      READ( IU_RST, IOSTAT=IOS ) ( ( ( TRACER(I,J,L), I=1,NI )
     &, J=1,NJ ), L=1,NL )

      PRINT *,'WEAK_CONSTRAINT: TRACER ',MINVAL(TRACER),MAXVAL(TRACER)

      IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_checkpt_file:9')

         !==============================================================
         ! Assign data from the TRACER array to the FORCE_U_FULLGRID array.
         !==============================================================

         ! Only process checkpoint data (i.e. mixing ratio)
      !IF ( CATEGORY(1:8) == 'IJ-CHK-$' ) THEN

            ! Make sure array dimensions are of global size
            ! (NI=IIPAR; NJ=JJPAR, NL=LLPAR), or stop the run
      !      CALL CHECK_DIMENSIONS( NI, NJ, NL )


      DO L = 1, LLPAR
         DO J = 1, JJPAR
            DO I = 1, IIPAR
               FORCE_U_FULLGRID(I,J,L) = TRACER(I,J,L)
            ENDDO
         ENDDO
      ENDDO

         !ENDIF
!      ENDDO

      ! Close file

      CLOSE( IU_RST )

      END SUBROUTINE READ_FORCE_U_FILE

      END MODULE WEAK_CONSTRAINT_MOD
