      MODULE INV_HESSIAN_MOD
!
!*****************************************************************************
!  Module INV_HESSIAN_MOD contains all the subroutines that are used for
!  calculating the approximate inverse Hessian.  (dkh, 05/15/07, adj32_012)
!
!  Module Variables:
!  ============================================================================
!  (1 ) IIMAP          (INTEGER)   : 4D to 1D mapping array
!  (2 ) EMS_SF_OLD      (REAL*8)   : Scaling factors at previous iteration
!  (2 ) EMS_SF_ADJ_OLD  (REAL*8)   : Gradients at previous iteration
!
!  Module Routines
!  ============================================================================
!  (1 ) UPDATE_HESSIAN             : Updates inv Hessian estimate
!  (2 ) MAKE_HESS_FILE             : Saves inv Hessian to file
!  (3 ) INIT_INV_HESSIAN           : Allocates and intializes module variables
!  (3 ) CLEANUP_INV_HESSIAN        : Deallocates  module variables
!
!  NOTES:
!  (1 ) Now make the working arrays allocatable (nb, dkh, 08/02/12, adj33g)
!*****************************************************************************
!
      IMPLICIT NONE

#     include "define_adj.h"    ! obs operators

      !====================================================================
      ! MODULE VARIABLES  ( those that used to be program variables )
      !====================================================================
      INTEGER, ALLOCATABLE :: IIMAP(:,:,:,:)
      REAL*8,  ALLOCATABLE :: EMS_SF_OLD(:,:,:,:)
      REAL*8,  ALLOCATABLE :: EMS_SF_ADJ_OLD(:,:,:,:)
      INTEGER, ALLOCATABLE :: MAPI(:), MAPJ(:)
      INTEGER, ALLOCATABLE :: MAPM(:), MAPN(:)
      REAL*8,  ALLOCATABLE :: HINV(:,:)
      REAL*8,  ALLOCATABLE :: Y(:)
      REAL*8,  ALLOCATABLE :: S(:)
      REAL*8,  ALLOCATABLE :: SST(:,:)
      REAL*8,  ALLOCATABLE :: HINVY(:)
      REAL*8,  ALLOCATABLE :: YTHINV(:)
      REAL*8,  ALLOCATABLE :: HINVYYTHINV(:,:)
      REAL*8,  ALLOCATABLE :: FILTER(:,:)

      !====================================================================
      ! MODULE ROUTINES
      !====================================================================

      CONTAINS

!-----------------------------------------------------------------------------

      SUBROUTINE UPDATE_HESSIAN( )
!
!******************************************************************************
!  Subroutine UPDATE_HESSIAN constructs an approximation of the inverse
!  Hessian using the DFP formula (see Muller and Stavrakou, 2005, eqn 18).
!  (dkh, 05/15/07)
!
!  This routine is set up to be used offline so that the Hessian is
!  only approximated after the results from a completed optimization have
!  been obtained.

!  To implement, first do a normal optimization with LINVH = .FALSE., and
!  keep the cost function, scaling factor and gradient files in OptData.
!
!  Next, set the LINVH flag in input.gcadj to TRUE and rerun from X=1 to
!  XSTOP=z, where z is the final number of optimization steps previously
!  completed. Now execute the run script. The outputs are in diagadj
!  directory.
!
!  The routine will label the output files according to function evaluation
!  number (X), although it will only includ the accepted iterations in the
!  calculation, not the line search evaluations.
!
!  The initial estimate of HINV can be identiy matrix or an initial
!  estimate of uncertainty.  At the moment it is hardwired into the
!  intial definition of HINV in the code below.
!
!  WARNING: It is easy to max the dimension of the inverse Hessian so large
!  that your code will crash.  It may not even compile (error like
!  "relocation truncated to fit: R_X86_64_32S against...").
!
!  For example, the inverse Hessian will require y Mb of memory, where
!    y    = 3 * HMAX ^ 2 * 8 / 10^6
!    HMAX = IPAR * JJPAR * MMSCL * NNEMS
!
!  The 8 comes from 8 bits/byte (could be half this if used REAL*4), and the 3
!  comes from the fact that we have 3 arrays that are size(HMAX,HMAX).  Thus,
!  at 4x5 resolution with 33 emissions sectors (NNEMS) and 1 time group (MMSCL),
!  the memory requirements in double precision are nearly 300 Gb!  Or > 4 Gb for
!  NNEMS = 1 at 2x2.5.
!
!  Thus, if it takes too long, or too much memory, to consider all possible correlations,
!  one can apply a filter when developing the mapping array, and then set HMAX
!  to an appropriate value (a dry run may be necessary to see what value HMAX
!  should be).  But this is cheating with bad math, so you feel ashamed.
!
!  If you need to compile with arrays that are larger than the available memory,
!  utlizile swap space instead (warning: could get slow) with the -mcmodel compile
!  flag (ifort).
!
!  Module Variable as Input:
!  ============================================================================
!  (1 ) EMS_SF       : Emissions scaling factors at the current iteration
!  (2 ) EMS_SF_ADJ   : Emissions gradients at the current iteration
!  (3 ) N_CALC       : Current interation number
!
!  NOTES:
!  (1 ) Updated for adj32 (dkh, 01/11/12).
!  (2 ) Now initialize inv Hess to y^T s / y^T y (nb, dkh, 08/02/12, adj33g)
!******************************************************************************
!
      ! Reference to f90 modules
      USE ADJ_ARRAYS_MOD, ONLY : MMSCL, NNEMS
      USE ADJ_ARRAYS_MOD, ONLY : N_CALC
      USE ADJ_ARRAYS_MOD, ONLY : EMS_SF
      USE ADJ_ARRAYS_MOD, ONLY : EMS_SF_ADJ
      USE ERROR_MOD,      ONLY : ERROR_STOP

#     include "CMN_SIZE"

      ! Local variables
      ! Here I've hardwired this to IIPAR * JJPAR, as we filter for just NNEMS = 1
      ! below, and MMSCL = 1.
      !INTEGER, PARAMETER      :: HMAX = 72 * 46

      !Now HMAX is a variable (nab, 7/16/12)
      !HMAX = IIPAR * JJPAR * MMSCL * NNEMS

      INTEGER                 :: HMAX


      INTEGER                 :: I, J, M, N, II, JJ, NITR

      LOGICAL, SAVE           :: FIRST  = .TRUE.
      LOGICAL, SAVE           :: SECOND = .TRUE.


      REAL*8                  :: YTS_INV, YTHINVY,YTS,YTY
      REAL*8                  ::  YTHINVY_INV


      !=================================================================
      ! UPDATE_HESSIAN begins here!
      !=================================================================

      HMAX = IIPAR * JJPAR * MMSCL * NNEMS


      PRINT*, ' UPDATE HESSIAN AT ITERATE ', N_CALC

      IF ( FIRST ) THEN

         ! allocate and initialize arrays
         CALL INIT_INV_HESSIAN( HMAX )

         print*, ' FILTER sum = ', sum(FILTER)

         ! Initialize HINV to the identity matrix (or initial unc. est)
         HINV(:,:)   = 0d0

         DO JJ = 1, HMAX
         DO II = 1, HMAX

           ! for example, 30% uncertainty
           IF ( II == JJ ) HINV(II,II) = 1d0

         ENDDO
         ENDDO

         II = 0

         DO N = 1, NNEMS
         DO M = 1, MMSCL
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            !==============================================
            ! Apply filters
            !==============================================

            ! Only in places where emissions are nonzero
            !IF ( ABS(ADJ_EMS(I,J,M,N)) < 1d-4 ) CYCLE
            !IF ( ABS(EMS_SF_ADJ(I,J,M,N)) < 1d-5 ) CYCLE
            !IF ( FILTER(I,J) < 0.99 ) CYCLE

            ! Only correlation of the first emissions sector with itself
            !IF ( N /= 1 ) CYCLE

               ! Update vector index
               II = II + 1

               ! Save mapping arrays
               IIMAP(I,J,M,N) = II
               MAPI(II) = I
               MAPJ(II) = J
               MAPM(II) = M
               MAPN(II) = N

            !ENDIF

         ENDDO
         ENDDO
         ENDDO
         ENDDO


         EMS_SF_OLD(:,:,:,:)     = EMS_SF(:,:,:,:)
         EMS_SF_ADJ_OLD(:,:,:,:) = EMS_SF_ADJ(:,:,:,:)

         print*, ' UPDATE HESSIAN, pts founds = ', II

         CALL MAKE_HESS_FILE( HINV, HMAX, 1 )

         FIRST = .FALSE.

         RETURN
      ENDIF


      DO II = 1, HMAX

         I = MAPI(II)
         J = MAPJ(II)
         M = MAPM(II)
         N = MAPN(II)

         ! find s_k = f_{k+1} - f_{k}
         S(II) = EMS_SF(I,J,M,N) - EMS_SF_OLD(I,J,M,N)

         ! find y_k = grad_{k+1} - grad_{k}
         Y(II) = EMS_SF_ADJ(I,J,M,N) - EMS_SF_ADJ_OLD(I,J,M,N)

      ENDDO

      print*, ' UPDATE HESSIAN, pts founds = ', II

      ! Rotate
      EMS_SF_OLD(:,:,:,:)     = EMS_SF(:,:,:,:)
      EMS_SF_ADJ_OLD(:,:,:,:) = EMS_SF_ADJ(:,:,:,:)

      !----------------------------------------------------------
      ! Update inverse Hessian
      !----------------------------------------------------------

      ! y^T*s
      YTS = 0d0
      DO II = 1, HMAX

         YTS = YTS + Y(II) * S(II)

      ENDDO

      print*, ' YTS = ', YTS , N_CALC

      ! Initialize inv Hessian to y^T s / y^T y
      IF ( SECOND ) THEN

         print*, ' Initialize inv Hessian to y^T s / y^T y  = ',YTS/YTY

         ! y^T * y
         YTY = 0d0
         DO II = 1, HMAX
            YTY = YTY + Y(II) * Y(II)
         ENDDO

         IF ( YTY < 1D-38 ) THEN
            CALL ERROR_STOP('underflow','inv_hessian')
         ENDIF

         DO II = 1, HMAX
            HINV(II,II) = YTS / YTY
         ENDDO

         SECOND = .FALSE.

      ENDIF

      ! s * s^T / YTS
      DO II = 1, HMAX
      DO JJ = 1, HMAX

         SST(II,JJ) = S(II) * S(JJ)

      ENDDO
      ENDDO

      ! HINV * y
      DO II = 1, HMAX

         HINVY(II) = 0D0

         DO JJ = 1, HMAX

            HINVY(II) = HINVY(II) + HINV(II,JJ) * Y(JJ)

         ENDDO
      ENDDO

      ! y^T * HINV
      DO JJ = 1, HMAX

         YTHINV(JJ) = 0d0

         DO II = 1, HMAX

            YTHINV(JJ) = YTHINV(JJ) + Y(II) * HINV(II,JJ)

         ENDDO
      ENDDO


      ! HINVY * YTHINV
      DO JJ = 1, HMAX
      DO II = 1, HMAX

            HINVYYTHINV(II,JJ) = HINVY(II) * YTHINV(JJ)

      ENDDO
      ENDDO


      ! YT * HINVY
      YTHINVY = 0d0
      DO II = 1, HMAX
         YTHINVY = YTHINVY + Y(II) * HINVY(II)
      ENDDO
      print*, 'YTHINVY = ', YTHINVY

      ! HINV = HINV + SST * (1/YTS) - HINVYYTHINV * (1/YTHINVY)
      YTS_INV     = 1 / YTS
      YTHINVY_INV = 1 / YTHINVY
      DO JJ = 1, HMAX
      DO II = 1, HMAX

         HINV(II,JJ) = HINV(II,JJ)
     &               + SST(II,JJ)         * YTS_INV
     &               - HINVYYTHINV(II,JJ) * YTHINVY_INV

      ENDDO
      ENDDO

      print*, ' MAX HINV = ', MAXVAL(HINV(:,:))
      print*, ' MIN HINV = ', MINVAL(HINV(:,:))

      NITR = N_CALC

      CALL MAKE_HESS_FILE( HINV, HMAX, NITR )

      ! Return to calling program
      END SUBROUTINE UPDATE_HESSIAN
!------------------------------------------------------------------------------

      SUBROUTINE MAKE_HESS_FILE( HINV, HMAX, NITR )
!
!******************************************************************************
!  Subroutine MAKE_HESS_FILE creates a binary file of selected elements
!  of the approximate inverse hessian. (dkh, 05/15/07)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) HINV      : Current estimate of inverse hessian
!  (2 ) HMAX      : Dimension
!  (3 ) NITR      : Current iteration
!
!  Module Variable as Input:
!  ============================================================================
!  (1 ) IIMAP     : 3D to 1D mappying array
!
!  NOTES:
!  (1 ) Just like MAKE_GDT_FILE except
!        - pass NITR as an argument
!  (2 ) Updated for adj32 (dkh, 01/11/12)
!******************************************************************************
!
      ! References to F90 modules
      USE ADJ_ARRAYS_MOD,    ONLY : MMSCL, NNEMS
      USE ADJ_ARRAYS_MOD,    ONLY : EXPAND_NAME
      USE BPCH2_MOD
      USE DIRECTORY_ADJ_MOD, ONLY : DIAGADJ_DIR
      USE ERROR_MOD,         ONLY : DEBUG_MSG, ERROR_STOP
      USE FILE_MOD,          ONLY : IU_RST,      IOERROR
      USE GRID_MOD,          ONLY : GET_XOFFSET, GET_YOFFSET
      USE LOGICAL_MOD,       ONLY : LPRT
      USE LOGICAL_ADJ_MOD,   ONLY : LICS
      USE LOGICAL_ADJ_MOD,   ONLY : LADJ_EMS
      USE TIME_MOD,          ONLY : EXPAND_DATE, GET_TAU

#     include "CMN_SIZE"          ! Size parameters


      ! Arguments
      INTEGER                    :: HMAX
      REAL*8                     :: HINV(HMAX,HMAX)
      INTEGER                    :: NITR

      ! Local Variables
      INTEGER                    :: I,    I0, IOS, J
      INTEGER                    :: J0, L, M, N, II, JJ
      INTEGER                    :: YYYY, MM, DD,  HH, SS
      REAL*4                     :: TRACER(IIPAR,JJPAR,LLPAR)
      REAL*4                     :: EMS_3D(IIPAR,JJPAR,MMSCL)
      CHARACTER(LEN=255)         :: FILENAME

      ! For binary punch file, version 2.0
      REAL*4                     :: LONRES, LATRES
      INTEGER, PARAMETER         :: HALFPOLAR = 1
      INTEGER, PARAMETER         :: CENTER180 = 1

      CHARACTER(LEN=20)          :: OUTPUT_GDT_FILE
      CHARACTER(LEN=20)          :: MODELNAME
      CHARACTER(LEN=40)          :: CATEGORY
      CHARACTER(LEN=40)          :: UNIT
      CHARACTER(LEN=40)          :: RESERVED = ''
      CHARACTER(LEN=80)          :: TITLE

      !=================================================================
      ! MAKE_HESS_FILE begins here!
      !=================================================================

      ! Clear intermediate arrays
      EMS_3D(:,:,:) = 0d0

      ! Hardwire output file for now
      OUTPUT_GDT_FILE = 'gctm.invhess.NN'

      ! Define variables for BINARY PUNCH FILE OUTPUT
      TITLE    = 'GEOS-CHEM Adjoint File: ' //
     &           'Inverse hessian  '
      UNIT     = 'none'
      CATEGORY = 'IJ-GDE-$'
      LONRES   = DISIZE
      LATRES   = DJSIZE

      ! Call GET_MODELNAME to return the proper model name for
      ! the given met data being used (bmy, 6/22/00)
      MODELNAME = GET_MODELNAME()

      ! Get the nested-grid offsets
      I0 = GET_XOFFSET( GLOBAL=.TRUE. )
      J0 = GET_YOFFSET( GLOBAL=.TRUE. )

      !=================================================================
      ! Open the adjoint file for output -- binary punch format
      !=================================================================

      ! Copy the output observation file name into a local variable
      FILENAME = TRIM( OUTPUT_GDT_FILE )

      ! Append the iteration number suffix to the file name
      CALL EXPAND_NAME( FILENAME, NITR )

      ! Add the OPT_DATA_DIR prefix to the file name
      FILENAME = TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )

      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - MAKE_HESS_FILE: Writing ', a )

      ! Open checkpoint file for output
      CALL OPEN_BPCH2_FOR_WRITE( IU_RST, FILENAME, TITLE )

      IF ( LICS ) THEN

         CALL ERROR_STOP( 'inverse hessian not supported ',
     &                    ' MAKE_HESS_FILE, inverse_mod.f')

      ELSEIF ( LADJ_EMS ) THEN

         !=================================================================
         ! Write the standard error of each optimized scaling factor
         !=================================================================
         DO N = 1, NNEMS

            !Temporarily store quantities in the TRACER array
            EMS_3D(:,:,:) = 0d0

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, M, II )
            DO M = 1, MMSCL
            DO J = 1, JJPAR
            DO I = 1, IIPAR


               II = IIMAP(I,J,M,N)
               IF ( II == 0 ) CYCLE

                  IF ( HINV(II,II) > 0 )  THEN
                     EMS_3D(I,J,M) = REAL(SQRT(HINV(II,II)))
                  ELSE
                     print*, I, J, M, N, II
                     CALL ERROR_STOP('non positive hessian diagonal ',
     &                               'inverse_mod.f')
                  ENDIF

            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

            CALL BPCH2( IU_RST,    MODELNAME, LONRES, LATRES,
     &                  HALFPOLAR, CENTER180, CATEGORY,  N,
     &                  UNIT,      GET_TAU(), GET_TAU(), RESERVED,
     &                  IIPAR,     JJPAR,     MMSCL,     I0+1,
     &                  J0+1,      1,         EMS_3D )

         ENDDO

!         ! Reset CATEGORY as labeling in gamap is different
!         CATEGORY = 'IJ-COREL'
!
!         !=================================================================
!         ! Write correlation of optimized scale factors with a particular
          ! target cell, selected manually below.
!         !=================================================================
!         DO N = 1, NNEMS
!
!            ! target cell
!            JJ = IIMAP(13,33,1,IDADJEMS_ENH3_an)
!
!            !Temporarily store quantities in the TRACER array
!            EMS_3D(I,J,M) = 0d0
!
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( I, J, M, II )
!            DO M = 1, MMSCL
!            DO J = 1, JJPAR
!            DO I = 1, IIPAR
!
!
!               II = IIMAP(I,J,M,N)
!               !IF ( II == 0 ) CYCLE
!               IF ( II == 0 ) THEN
!                  EMS_3D(I,J,M) = 0d0
!               ELSE
!                  EMS_3D(I,J,M) = REAL(HINV(II,JJ)/(SQRT(HINV(II,II))
!     &                          * SQRT(HINV(JJ,JJ))))
!               ENDIF
!
!            ENDDO
!            ENDDO
!            ENDDO
!!$OMP END PARALLEL DO
!
!            CALL BPCH2( IU_RST,    MODELNAME, LONRES, LATRES,
!     &                  HALFPOLAR, CENTER180, CATEGORY,  N,
!     &                  UNIT,      GET_TAU(), GET_TAU(), RESERVED,
!     &                  IIPAR,     JJPAR,     MMSCL,     I0+1,
!     &                  J0+1,      1,         EMS_3D )
!
!         ENDDO
      ELSE
         CALL ERROR_STOP( 'simulation type not defined!',
     &                    'MAKE_HESS_FILE' )
      ENDIF

      ! Close file
      CLOSE( IU_RST )

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### MAKE_HESS_FILE: wrote file' )

      ! Return to calling program
      END SUBROUTINE MAKE_HESS_FILE

!-----------------------------------------------------------------------------

      SUBROUTINE INIT_INV_HESSIAN(HMAX)
!
!******************************************************************************
!  Subroutine INIT_INV_HESSIAN initializes and zeros all allocatable arrays
!
!  NOTES:
!******************************************************************************

      ! References to F90 modules
      USE ADJ_ARRAYS_MOD, ONLY : NNEMS, MMSCL
      USE ERROR_MOD, ONLY      : ALLOC_ERR

#     include "CMN_SIZE"       ! Size parameters

      ! Local variables
      INTEGER                 :: AS, I,HMAX

      !=================================================================
      ! INIT_INV_HESSIAN begins here!
      !=================================================================

      ALLOCATE( IIMAP(IIPAR,JJPAR,MMSCL,NNEMS), STAT = AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'IIMAP' )
      IIMAP = 0

      ALLOCATE( EMS_SF_OLD(IIPAR,JJPAR,MMSCL,NNEMS), STAT = AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EMS_SF_OLD' )
      EMS_SF_OLD = 0d0

      ALLOCATE( EMS_SF_ADJ_OLD(IIPAR,JJPAR,MMSCL,NNEMS), STAT = AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EMS_SF_ADJ_OLD' )
      EMS_SF_ADJ_OLD = 0d0

      ALLOCATE (MAPI(HMAX), STAT = AS)
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'MAPI' )
      MAPI = 0

      ALLOCATE (MAPJ(HMAX), STAT = AS)
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'MAPJ' )
      MAPJ = 0

      ALLOCATE (MAPM(HMAX), STAT = AS)
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'MAPM' )
      MAPM = 0

      ALLOCATE (MAPN(HMAX), STAT = AS)
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'MAPN' )
      MAPN = 0

      ALLOCATE (HINV(HMAX,HMAX), STAT = AS)
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'HINV' )
      HINV = 0

      ALLOCATE (HINVYYTHINV(HMAX,HMAX), STAT = AS)
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'HINVYYTHINV' )
      HINVYYTHINV = 0

      ALLOCATE (S(HMAX), STAT = AS)
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'S' )
      S = 0

      ALLOCATE (Y(HMAX), STAT = AS)
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'Y' )
      Y = 0

      ALLOCATE (SST(HMAX,HMAX), STAT = AS)
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SST' )
      SST = 0

      ALLOCATE (HINVY(HMAX), STAT = AS)
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'HINVY' )
      HINVY = 0

      ALLOCATE (YTHINV(HMAX), STAT = AS)
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'YTHINV' )
      YTHINV = 0

      ALLOCATE( FILTER(IIPAR,JJPAR), STAT = AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'FILTER' )
      FILTER = 0d0

       !CALL READ_FILTER()

      END SUBROUTINE INIT_INV_HESSIAN

!------------------------------------------------------------------------------

      ! Return to calling program
      SUBROUTINE CLEANUP_INV_HESSIAN
!
!******************************************************************************
!  Subroutine CLEANUP_INV_HESSIAN deallocates all previously allocated arrays
!  for inverse_mod -- call at the end of the program (dkh, 01/11/12)
!
!  NOTES:
!******************************************************************************
!
      !=================================================================
      ! CLEANUP_INV_HESSIAN begins here!
      !=================================================================
      IF ( ALLOCATED( IIMAP          ) ) DEALLOCATE( IIMAP          )
      IF ( ALLOCATED( EMS_SF_OLD     ) ) DEALLOCATE( EMS_SF_OLD     )
      IF ( ALLOCATED( EMS_SF_ADJ_OLD ) ) DEALLOCATE( EMS_SF_ADJ_OLD )

      ! Return to calling program
      END SUBROUTINE CLEANUP_INV_HESSIAN

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

      SUBROUTINE READ_FILTER( )
!
!******************************************************************************
!
      ! References to F90 modules
      USE RESTART_MOD,       ONLY : CHECK_DIMENSIONS
      USE BPCH2_MOD,         ONLY : OPEN_BPCH2_FOR_READ
      USE ERROR_MOD,         ONLY : DEBUG_MSG
      USE FILE_MOD,          ONLY : IU_RST, IOERROR
      USE LOGICAL_MOD,       ONLY : LPRT
      USE TIME_MOD,          ONLY : EXPAND_DATE
      USE TRACER_MOD,        ONLY : N_TRACERS
      USE ADJ_ARRAYS_MOD,    ONLY : NNEMS

#     include "CMN_SIZE"   ! Size parameters

      ! Local Variables
      INTEGER             :: I, IOS, J, L, N
      INTEGER             :: NCOUNT(NNPAR)
      REAL*4              :: TRACER(IIPAR,JJPAR,1)
      REAL*8              :: SUMTC
      CHARACTER(LEN=255)  :: FILENAME
      CHARACTER(LEN=255)  :: UNZIP_FILE_CMD
      CHARACTER(LEN=255)  :: ZIP_FILE_CMD

      ! For binary punch file, version 2.0
      INTEGER             :: NI,     NJ,     NL
      INTEGER             :: IFIRST, JFIRST, LFIRST
      INTEGER             :: NTRACER,   NSKIP
      INTEGER             :: HALFPOLAR, CENTER180
      REAL*4              :: LONRES,    LATRES
      REAL*8              :: ZTAU0,     ZTAU1
      CHARACTER(LEN=20)   :: MODELNAME
      CHARACTER(LEN=40)   :: CATEGORY
      CHARACTER(LEN=40)   :: UNIT
      CHARACTER(LEN=40)   :: RESERVED

      !=================================================================
      ! READ_OBS_FILE begins here!
      !=================================================================

      ! Hardwire output file for now
      FILENAME = TRIM('gctm.filter.3293')

      ! Initialize some variables
      NCOUNT(:)     = 0
      TRACER(:,:,:) = 0e0

      !=================================================================
      ! Open observation file and read top-of-file header
      !=================================================================

      ! Echo some input to the screen
      WRITE( 6, '(a)'   ) REPEAT( '=', 79 )
      WRITE( 6, '(a,/)' ) 'F I L T E R     F I L E   I N P U T'

      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - READ_FILTER: Reading ', a )

      ! Open the binary punch file for input
      CALL OPEN_BPCH2_FOR_READ( IU_RST, FILENAME )

      !=================================================================
      ! Read concentrations -- store in the TRACER array
      !=================================================================
      !DO N = 1, NOBS
      DO N = 1, NNEMS
         READ( IU_RST, IOSTAT=IOS )
     &     MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180

         ! IOS < 0 is end-of-file, so exit
         IF ( IOS < 0 ) EXIT

         ! IOS > 0 is a real I/O error -- print error message
         IF ( IOS > 0 ) CALL IOERROR( IOS,IU_RST,'read_obs_file:4' )

         READ( IU_RST, IOSTAT=IOS )
     &        CATEGORY, NTRACER,  UNIT, ZTAU0,  ZTAU1,  RESERVED,
     &        NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
     &        NSKIP

         IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_obs_file:5')

         READ( IU_RST, IOSTAT=IOS )
     &        ( ( ( TRACER(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )

         IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_obs_file:6')

         !==============================================================
         ! Assign data from the TRACER array to the STT array.
         !==============================================================

         ! Only process observation data (i.e. aerosol and precursors)
         IF ( N == 3 ) THEN

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
            DO J = 1, JJPAR
            DO I = 1, IIPAR
               FILTER(I,J) = TRACER(I,J,1)
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

         ENDIF
      ENDDO

      ! Close file
      CLOSE( IU_RST )

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### READ_FILTER: read file' )

      ! Return to calling program
      END SUBROUTINE READ_FILTER
!------------------------------------------
      END MODULE INV_HESSIAN_MOD


