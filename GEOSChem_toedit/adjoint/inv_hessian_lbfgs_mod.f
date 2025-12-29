!$ID$
      MODULE INV_HESSIAN_LBFGS_MOD
!
!*****************************************************************************
!  Module INV_HESSIAN_LBFGS_MOD contains all the subroutines that are used for
!  calculating the diagonal of the approximate L-BFGS inverse Hessian.
!  (nab, 03/23/12, adj32_027)
!
!  Module Variables:
!  ============================================================================
!  (1 ) IIMAP          (INTEGER)   : 4D to 1D mapping array
!  (2 ) EMS_SF_OLD      (REAL*8)   : Scaling factors at previous iteration
!  (3 ) EMS_SF_ADJ_OLD  (REAL*8)   : Gradients at previous iteration
!
!  Module Routines
!  ============================================================================
!  (1 ) LBFGS_INV_HESSIAN          : Updates inv Hessian estimate
!  (2 ) READ_GDT_FILE_AT           : Read gradient file
!  (3 ) READ_SF_FILE_AT            : Read scaling factor file
!  (4 ) MAKE_HESS_DIAG_FILE        : Saves diagonal of inv Hessian to file
!  (5 ) INIT_INV_HESSIAN           : Allocates and intializes module variables
!  (6 ) CLEANUP_INV_HESSIAN        : Deallocates  module variables
!
!  NOTES:
!  (1 ) Now calculate the inverse Hessian approximation for all emissions (nab)
!
!
!  For any question please contact:
!  Contact: Nicolas Bousserez (Nicolas.Bousserez@colorado.edu)
!******************************************************************************
!
      IMPLICIT NONE

#     include "define_adj.h"    ! obs operators

      !====================================================================
      ! MODULE VARIABLES  ( those that used to be program variables )
      !====================================================================
      INTEGER, ALLOCATABLE     :: IIMAP(:,:,:,:)
      REAL*8,  ALLOCATABLE     :: EMS_SF_OLD(:,:,:,:)
      REAL*8,  ALLOCATABLE     :: EMS_SF_ADJ_OLD(:,:,:,:)
      REAL*8,  ALLOCATABLE     :: ICS_SF_OLD(:,:,:,:)
      REAL*8,  ALLOCATABLE     :: ICS_SF_ADJ_OLD(:,:,:,:)

      INTEGER, ALLOCATABLE     :: MAPI(:), MAPJ(:), MAPL(:)
      INTEGER, ALLOCATABLE     :: MAPM(:), MAPN(:)
      REAL*8 , ALLOCATABLE     :: HINVD(:)

      !====================================================================
      ! MODULE ROUTINES
      !====================================================================

      CONTAINS

!-----------------------------------------------------------------------------
      SUBROUTINE LBFGS_INV_HESSIAN( MK )
!
!******************************************************************************
!  Compute the diagonal terms of the posterior error covariance
!  matrix using the L-BFGS formula:
!
! Byrd, R. H.; Lu, P.; Nocedal, J.; Zhu, C. (1995). "A Limited Memory Algorithm
! for Bound Constrained Optimization". SIAM Journal on Scientific Computing 16 (5): 1190
!
!  Variable as Input:
!  ============================================================================
!  (1 ) MK           : Number of last previous iteration used in the inverse
!                      Hessian L-BFGS approximation
!
!  Module variables as Output:
!  ============================================================================
!
!  NOTES:
!
!******************************************************************************
!
      ! Reference to f90 modules
      USE ADJ_ARRAYS_MOD,       ONLY : MMSCL, NNEMS
      USE ADJ_ARRAYS_MOD,       ONLY : N_CALC
      USE ADJ_ARRAYS_MOD,       ONLY : EMS_SF,ICS_SF
      USE ADJ_ARRAYS_MOD,       ONLY : EMS_SF_ADJ,ICS_SF_ADJ
      USE ERROR_MOD,            ONLY : ALLOC_ERR , ERROR_STOP
      USE LOGICAL_ADJ_MOD,      ONLY : LICS, LADJ_EMS
      USE MKL95_BLAS
      USE MKL95_LAPACK
      USE TRACER_MOD,           ONLY : N_TRACERS

#     include "CMN_SIZE"

      ! Arguments
      INTEGER, OPTIONAL, INTENT(IN) :: MK

      ! Local variables
      REAL(8)                       :: TMPDP1, TMPDP2
      REAL(8)                       :: THETA,THETA_DENO,THETA_NUM
      REAL(8), ALLOCATABLE          :: WS(:,:)
      REAL(8), ALLOCATABLE          :: WY(:,:)
      REAL(8), ALLOCATABLE          :: RK(:,:)
      REAL(8), ALLOCATABLE          :: DK(:,:)
      REAL(8), ALLOCATABLE          :: MM(:,:)
      CHARACTER(LEN=255)            :: MSG
      INTEGER                       :: KK, INFO
      INTEGER                       :: MM1
      INTEGER                       :: I, J,L, M, N, II, JJ, NITR, AS
      INTEGER                       :: HMAX

      !=================================================================
      ! LBFGS_INV_HESSIAN begins here!
      !=================================================================

      PRINT*,'********************************************'
      PRINT*,'STARTING L-BFGS INVERSE HESSIAN CALCULATION'
      PRINT*,'********************************************'

      IF ( LADJ_EMS ) THEN
        HMAX = IIPAR * JJPAR * MMSCL * NNEMS
      ELSEIF ( LICS ) THEN
        HMAX = IIPAR * JJPAR * LLPAR * N_TRACERS
      ENDIF

      ! allocate and initialize arrays
      CALL INIT_INV_HESSIAN( HMAX )

      KK = N_CALC

      MM1 = KK
      IF ( PRESENT(MK) ) MM1 = MK
      MM1      = MIN( KK - 1, MM1 )
      HINVD(:) = 0d0

      IF( KK == 0 ) THEN

        HINVD(:) = 1d0

        RETURN

      ENDIF

      ! ALLOCATIONS
      ALLOCATE( WS(MM1,HMAX), STAT = AS )
      IF ( AS /= 0) CALL ALLOC_ERR( 'WS', AS )
      WS = 0d0

      ALLOCATE(WY(HMAX,MM1),STAT=AS)
      IF ( AS /= 0) CALL ALLOC_ERR( 'WY', AS )
      WY = 0d0


      ALLOCATE( RK(MM1,MM1), STAT = AS )
      IF ( AS /= 0) CALL ALLOC_ERR( 'RK', AS )
      RK = 0d0

      ALLOCATE( DK(MM1,MM1), STAT = AS )
      IF ( AS /= 0) CALL ALLOC_ERR( 'DK', AS )
      DK = 0d0

      ALLOCATE( MM(MM1,HMAX), STAT = AS )
      IF ( AS /= 0) CALL ALLOC_ERR( 'MM', AS )
      MM = 0d0

      II = 0

      IF ( LADJ_EMS ) THEN

         DO N = 1, NNEMS
         DO M = 1, MMSCL
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            !==============================================
            ! Apply filters
            !==============================================

            ! Only in places where emissions are nonzero
            !IF ( ABS(ADJ_EMS(I,J,M,N)) < 1d-4 ) CYCLE

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


         !==================================================
         ! Start L-BFGS inverse Hessian diagonal elements extraction
         !==================================================

         DK = 0d0


         DO JJ = 1, MM1

            CALL READ_GDT_FILE_AT( JJ + ( KK - 1 - MM1 ) )
            EMS_SF_ADJ_OLD(:,:,:,:) = EMS_SF_ADJ(:,:,:,:)
            CALL READ_GDT_FILE_AT( JJ + ( KK - 1 - MM1 ) + 1 )

            CALL READ_SF_FILE_AT( JJ + ( KK - 1 - MM1 ) )
            EMS_SF_OLD(:,:,:,:) = EMS_SF(:,:,:,:)
            CALL READ_SF_FILE_AT( JJ + ( KK - 1 - MM1 ) + 1 )


            DO II = 1, HMAX

               I = MAPI(II)
               J = MAPJ(II)
               M = MAPM(II)
               N = MAPN(II)


               ! s_k = f_{k+1} - f_{k}
               WS(JJ,II) = EMS_SF(I,J,M,N) - EMS_SF_OLD(I,J,M,N)

               ! y_k = grad_{k+1} - grad_{k}
               WY(II,JJ) = EMS_SF_ADJ(I,J,M,N) - EMS_SF_ADJ_OLD(I,J,M,N)

            ENDDO
         ENDDO


      ELSEIF ( LICS ) THEN

         DO N = 1, N_TRACERS
         DO J = 1, JJPAR
         DO I = 1, IIPAR
         DO L = 1, LLPAR

            !==============================================
            ! Apply filters
            !==============================================

            ! Only in places where emissions are nonzero
            !IF ( ABS(ADJ_EMS(I,J,M,N)) < 1d-4 ) CYCLE

               ! Update vector index
               II = II + 1

               ! Save mapping arrays
               IIMAP(I,J,L,N) = II
               MAPI(II) = I
               MAPJ(II) = J
               MAPN(II) = N
               MAPL(II) = L

            !ENDIF

         ENDDO
         ENDDO
         ENDDO
         ENDDO

         !==================================================
         ! Start L-BFGS inverse Hessian diagonal elements extraction
         !==================================================

         DK = 0d0


         DO JJ = 1, MM1

            CALL READ_GDT_FILE_AT( JJ + ( KK - 1 - MM1 ) )
            ICS_SF_ADJ_OLD(:,:,:,:) = ICS_SF_ADJ(:,:,:,:)
            CALL READ_GDT_FILE_AT( JJ + ( KK - 1 - MM1 ) + 1 )

            CALL READ_SF_FILE_AT( JJ + ( KK - 1 - MM1 ) )
            ICS_SF_OLD(:,:,:,:) = ICS_SF(:,:,:,:)
            CALL READ_SF_FILE_AT( JJ + ( KK - 1 - MM1 ) + 1 )

            DO II = 1, HMAX

               I = MAPI(II)
               J = MAPJ(II)
               N = MAPN(II)
               L = MAPL(II)

               ! s_k = f_{k+1} - f_{k}
               WS(JJ,II) = ICS_SF(I,J,L,N) - ICS_SF_OLD(I,J,L,N)

               ! y_k = grad_{k+1} - grad_{k}
               WY(II,JJ) = ICS_SF_ADJ(I,J,L,N) -
     &                 ICS_SF_ADJ_OLD(I,J,L,N)

            ENDDO
         ENDDO

      ENDIF

      DO JJ = 1, MM1
         DK(JJ,JJ) = DOT( WS(JJ,:), WY(:,JJ) )
      ENDDO


      !theta = WY^T.WY/WY^T.WS
      !for inverse Hessian uses 1/theta
      THETA_DENO = DOT( WY(:,MM1), WY(:,MM1) )
      THETA_NUM  = DOT( WY(:,MM1), WS(MM1,:) )
      THETA      = THETA_NUM / THETA_DENO


!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( JJ , II )
      DO JJ = 1, MM1
      DO II = 1, MM1

         IF (  II <= JJ ) THEN
            RK(II,JJ) = DOT( WS(II,:), WY(:,JJ) )
         ENDIF

      ENDDO
      ENDDO
!$OMP END PARALLEL DO



!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( II )
      DO II=1, HMAX
          HINVD(II) = THETA
      ENDDO
!$OMP END PARALLEL DO

      ! assuming H_0=theta*I

      CALL GEMM( WY(:,1:MM1), WY(:,1:MM1), DK(1:MM1,1:MM1),
     &     TRANSA = 'T', ALPHA = THETA, BETA = 1d0 )

      ! solve for R^{-1}S' -> m x n
      CALL TRTRS( RK(1:MM1,1:MM1), WS(1:MM1,:), UPLO = 'U', INFO = INFO)

      IF( INFO < 0 ) THEN
          WRITE(MSG,'(a,i4,a)')'the ',i,
     &   '-th parameter had an illegal value.'
          CALL ERROR_STOP(MSG,'lbfgs_mod, get_pst_cov_diag')
      ENDIF

      ! (DK+Y'H_0Y)R^{-1}S'
      CALL GEMM( DK(1:MM1,1:MM1), WS(1:MM1,:), MM(1:MM1,:) )

      ! M-Y'H_0
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( JJ)
      DO JJ = 1, MM1
          CALL AXPY( WY(:,JJ), MM(JJ,:), A=-THETA )
      ENDDO
!$OMP END PARALLEL DO

       !H=H_0+SR^{-T}M
!      CALL gemm(WS(1:KK,:),MM(1:KK,:),HINVD,transa='T',beta=1d0)

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( II )
      DO II = 1, HMAX
         HINVD(II)=HINVD(II) + DOT( WS(1:MM1,II), MM(1:MM1,II) )
      ENDDO
!$OMP END PARALLEL DO

       !H=H-YH_0R^{-1}S'
!      CALL gemm(WY(:,1:KK),WS(1:KK,:),HINVD,alpha=-theta,beta=1d0)

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( II )
      DO II = 1, HMAX
         HINVD(II)=HINVD(II) - THETA * DOT( WY(II,1:MM1), WS(1:MM1,II) )
      ENDDO
!$OMP END PARALLEL DO

      PRINT*, ' MAX HINVD = ', MAXVAL(HINVD(:))
      PRINT*, ' MIN HINVD = ', MINVAL(HINVD(:))

      NITR = N_CALC

      PRINT*,'************************************'
      PRINT*,'NOW MAKING HESSIAN DIAGONAL BINARY FILE'
      PRINT*,'************************************'

      CALL MAKE_HESS_DIAG_FILE( HINVD,HMAX, NITR )


      IF(ALLOCATED(WS))DEALLOCATE(WS)
      IF(ALLOCATED(WY))DEALLOCATE(WY)
      IF(ALLOCATED(RK))DEALLOCATE(RK)
      IF(ALLOCATED(DK))DEALLOCATE(DK)
      IF(ALLOCATED(MM))DEALLOCATE(MM)


      CONTAINS

       FUNCTION SAFE_DIV( A, B ) RESULT( ANS )
       IMPLICIT NONE
       REAL*8, INTENT(IN) :: A, B
       REAL*8             :: ANS

       ANS = A * B / ( B * B + 1d-20 )

       RETURN

       END FUNCTION SAFE_DIV


      END SUBROUTINE LBFGS_INV_HESSIAN
!------------------------------------------------------------------------------

      SUBROUTINE READ_GDT_FILE_AT ( ITE )
!
!******************************************************************************
! Subroutine READ_GDT_FILE_AT reads the gctm.gdt file for a given iteration number
! into ADJ_xxx
!   (DKh, 9/17/04)

!
!  Module Variable as Input:
!  ============================================================================
!  (1 ) ITE    :  iteration number
!
!  Notes
!  (1 ) now CALLed GDT instead of ADJ
!  (2 ) Added ACTIVE_VARS == 'EMISSIONS' case.  (DKh, 11/27/04)
!  (3 ) Added ACTIVE_VARS == 'FDTEST' case. (DKh, 02/17/05)
!  (4 ) Now use CATEGORY = 'IJ-GDE-$' for EMISSIONS case. (DKh, 03/29/05)
!  (5 ) No longer pass COST_FUNC in the header; use cnf.* files. (DKh, 02/13/06)
!  (6 ) Now support strat fluxes LADJ_STRAT (hml, DKh, 02/20/12, adj32 _025)
!  (7)  Modified to read gctm.gdt files at a specific iteration number
!       (nab, 03/24/12, adj32_027) !

!******************************************************************************
!
      ! References to F90 modules
      USE ADJ_ARRAYS_MOD,      ONLY : ICS_SF_ADJ, EMS_SF_ADJ
      USE ADJ_ARRAYS_MOD,      ONLY : NNEMS, MMSCL, N_CALC
      USE ADJ_ARRAYS_MOD,      ONLY : EXPAND_NAME
      USE ADJ_ARRAYS_MOD,      ONLY : PROD_SF_ADJ,LOSS_SF_ADJ
      USE ADJ_ARRAYS_MOD,      ONLY : NSTPL
      USE BPCH2_MOD,           ONLY : OPEN_BPCH2_FOR_READ
      USE DIRECTORY_ADJ_MOD,   ONLY : OPTDATA_DIR
      USE ERROR_MOD,           ONLY : DEBUG_MSG, ERROR_STOP
      USE ERROR_MOD,           ONLY : DEBUG_MSG, ERROR_STOP
      USE FILE_MOD,            ONLY : IU_RST, IOERROR
      USE LOGICAL_ADJ_MOD,     ONLY : LICS, LADJ_EMS
      USE LOGICAL_ADJ_MOD,     ONLY : LADJ_STRAT
      USE LOGICAL_MOD,         ONLY : LPRT
      USE RESTART_MOD,         ONLY : CHECK_DIMENSIONS
      USE TIME_MOD,            ONLY : EXPAND_DATE
      USE TRACER_MOD,          ONLY : N_TRACERS

#     include "CMN_SIZE"   ! Size parameters

      ! Local Variables
      INTEGER  , INTENT(IN)  :: ITE
      INTEGER             :: I, IOS, J, L, M, N
      INTEGER             :: NCOUNT(NNPAR)
      REAL*4              :: TRACER(IIPAR,JJPAR,LLPAR)
      REAL*4              :: EMS_3D(IIPAR,JJPAR,MMSCL)
      REAL*4              :: PROD_3D(IIPAR,JJPAR,MMSCL)
      REAL*4              :: LOSS_3D(IIPAR,JJPAR,MMSCL)
      REAL*8              :: SUMTC
      CHARACTER(LEN=255)  :: FILENAME

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
      CHARACTER(LEN=20)   :: INPUT_GDT_FILE

      !=================================================================
      ! READ_GDT_FILE begins here!
      !=================================================================

      ! Hardwire output file for now
      INPUT_GDT_FILE = 'gctm.gdt.NN'

      ! Initialize some variables
      NCOUNT(:)     = 0
      TRACER(:,:,:) = 0e0

      !=================================================================
      ! Open gradient file and read top-of-file header
      !=================================================================

      ! Copy input file name to a local variable
      FILENAME = TRIM( INPUT_GDT_FILE )

      ! Replace NN tokens in FILENAME w/ actual values
      CALL EXPAND_NAME( FILENAME, ITE )

      ! Add OPTDATA_DIR prefix to FILENAME
      FILENAME = TRIM( OPTDATA_DIR ) // TRIM( FILENAME )

      ! Echo some input to the screen
      WRITE( 6, '(a)'   ) REPEAT( '=', 79 )
      WRITE( 6, '(a,/)' ) 'G D T   F I L E   I N P U T'
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - READ_GDT_FILE: Reading ', a )

      ! Open the binary punch file for input
      CALL OPEN_BPCH2_FOR_READ( IU_RST, FILENAME )


      IF ( LICS ) THEN
         !=================================================================
         ! Read adjoints -- store in the TRACER array
         !=================================================================
         DO N = 1, N_TRACERS
            READ( IU_RST, IOSTAT=IOS )
     &        MODELNAME, LONRES , LATRES, HALFPOLAR, CENTER180

            ! IOS < 0 is end-of-file, so exit
            IF ( IOS < 0 ) EXIT

            ! IOS > 0 is a REAL I/O error -- print error message
            IF ( IOS > 0 ) CALL IOERROR( IOS,IU_RST,'read_gdt_file:4' )

            READ( IU_RST, IOSTAT=IOS )
     &           CATEGORY, NTRACER,  UNIT, ZTAU0,  ZTAU1,  RESERVED,
     &           NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
     &           NSKIP

            IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_gdt_file:5')

            READ( IU_RST, IOSTAT=IOS )
     &           ( ( ( TRACER(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )

            IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_gdt_file:6')

            !==============================================================
            ! Assign data from the TRACER array to the ADJ_STT array.
            !==============================================================

            ! Only process observation data (i.e. aerosol and precursors)
            IF ( CATEGORY(1:8) == 'IJ-GDT-$' ) THEN

               ! Make sure array dimensions are of global size
               ! (NI=IIPAR; NJ=JJPAR, NL=LLPAR), or stop the run
               CALL CHECK_DIMENSIONS( NI, NJ, NL )

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
               DO L = 1, LLPAR
               DO J = 1, JJPAR
               DO I = 1, IIPAR
                  ICS_SF_ADJ(I,J,L,N) = TRACER(I,J,L)
               ENDDO
               ENDDO
               ENDDO
!$OMP END PARALLEL DO

            ENDIF
         ENDDO

      ENDIF
      IF ( LADJ_EMS ) THEN

         !=================================================================
         ! Read adjoints -- store in the TRACER array
         !=================================================================
         DO N = 1, NNEMS
            READ( IU_RST, IOSTAT=IOS )
     &        MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180

            ! IOS < 0 is end-of-file, so exit
            IF ( IOS < 0 ) EXIT

            ! IOS > 0 is a REAL I/O error -- print error message
            IF ( IOS > 0 ) CALL IOERROR( IOS,IU_RST,'read_gdt_file:4' )

            READ( IU_RST, IOSTAT=IOS )
     &           CATEGORY, NTRACER,  UNIT, ZTAU0,  ZTAU1,  RESERVED,
     &           NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
     &           NSKIP

            IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_gdt_file:5')
! debug (nab)
 !           PRINT*,'MMSCL: ',MMSCL
!!
            READ( IU_RST, IOSTAT=IOS )
     &           ( ( ( EMS_3D(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )

            IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_gdt_file:7')

            !==============================================================
            ! Assign data from the TRACER array to the ADJ_STT array.
            !==============================================================

            ! Only process observation data (i.e. aerosol and precursors)
            IF ( CATEGORY(1:8) == 'IJ-GDE-$' ) THEN
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, M )
               DO M = 1, MMSCL
               DO J = 1, JJPAR
               DO I = 1, IIPAR
                  EMS_SF_ADJ(I,J,M,N) = EMS_3D(I,J,M)
               ENDDO
               ENDDO
               ENDDO
!$OMP END PARALLEL DO

            ENDIF
         ENDDO

         IF ( LADJ_STRAT ) THEN
            !==============================================================
            ! Read adjoints -- store in the TRACER array
            !==============================================================
            DO N = 1, NSTPL
               READ( IU_RST, IOSTAT=IOS )
     &           MODELNAME, LONRES , LATRES, HALFPOLAR, CENTER180

               ! IOS < 0 is end-of-file, so exit
               IF ( IOS < 0 ) EXIT

               ! IOS > 0 is a REAL I/O error -- print error message
               IF ( IOS > 0 ) CALL IOERROR
     &                             ( IOS,IU_RST,'read_gdt_file:8' )

               READ( IU_RST, IOSTAT=IOS )
     &              CATEGORY, NTRACER,  UNIT, ZTAU0,  ZTAU1,  RESERVED,
     &              NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
     &              NSKIP

               IF ( IOS /= 0 ) CALL IOERROR
     &                              ( IOS,IU_RST,'read_gdt_file:9')

               READ( IU_RST, IOSTAT=IOS )
     &              ( ( ( PROD_3D(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )
               IF ( IOS /= 0 ) CALL IOERROR
     &                              ( IOS,IU_RST,'read_gdt_file:10')


               !===========================================================
               ! Assign data from the TRACER array to the ADJ_STT array.
               !===========================================================

               ! Only process observation data (i.e. aerosol and precursors)
               IF ( CATEGORY(1:8) == 'IJ-GDP-$' ) THEN

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, M )
                  DO M = 1, MMSCL
                  DO J = 1, JJPAR
                  DO I = 1, IIPAR
                     PROD_SF_ADJ(I,J,M,N) = PROD_3D(I,J,M)
                  ENDDO
                  ENDDO
                  ENDDO
!$OMP END PARALLEL DO
               ENDIF
            ENDDO

            DO N = 1, NSTPL
               READ( IU_RST, IOSTAT=IOS )
     &           MODELNAME, LONRES , LATRES, HALFPOLAR, CENTER180

               ! IOS < 0 is end-of-file, so exit
               IF ( IOS < 0 ) EXIT

               ! IOS > 0 is a REAL I/O error -- print error message
               IF ( IOS > 0 ) CALL IOERROR
     &                             ( IOS,IU_RST,'read_gdt_file:8b')

               READ( IU_RST, IOSTAT=IOS )
     &              CATEGORY, NTRACER,  UNIT, ZTAU0,  ZTAU1,  RESERVED,
     &              NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
     &              NSKIP

               IF ( IOS /= 0 ) CALL IOERROR
     &                              ( IOS,IU_RST,'read_gdt_file:9b')

               READ( IU_RST, IOSTAT=IOS )
     &              ( ( ( LOSS_3D(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )
               IF ( IOS /= 0 ) CALL IOERROR
     &                              (IOS,IU_RST,'read_gdt_file:10b')

               ! Only process observation data (i.e. aerosol and precursors)
               IF ( CATEGORY(1:8) == 'IJ-GDL-$' ) THEN
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, M )
                  DO M = 1, MMSCL
                  DO J = 1, JJPAR
                  DO I = 1, IIPAR
                     LOSS_SF_ADJ(I,J,M,N) = LOSS_3D(I,J,M)
                  ENDDO
                  ENDDO
                  ENDDO
!$OMP END PARALLEL DO
               ENDIF
            ENDDO
         ENDIF
      ENDIF

      ! Close file
      CLOSE( IU_RST )

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### READ_GDT_FILE: read file' )

      ! Return to CALLing program
      END SUBROUTINE READ_GDT_FILE_AT


!!------------------------------------------------------------------------------
!
      SUBROUTINE READ_SF_FILE_AT ( ITE )
!
!******************************************************************************
! Subroutine READ_SF_FILE_AT reads the gctm.sf file for a given iteration number
! into ADJ_xxx
!   (DKh, 9/17/04)
!
!  Module Variable as Input:
!  ============================================================================
!  (1 ) ITE    : iteration number
!
!  Notes
!  (1 ) Add support for ACTIVE_VARS == 'EMISSIONS' case (DKh, 11/27/04)
!  (2 ) Add support for ACTIVE_VARS == 'FDTEST' case (DKh, 02/17/05)
!  (3 ) Now use CATEGORY = 'IJ-EMS-$' for ACTIVE_VARS == 'EMISSIONS' case.
!       (DKh, 03/28/05)
!  (4 ) Change name from ICS to SF, replace CMN_ADJ (DKh, ks, mak, cs  06/08/09)
!  (5 ) Now support strat fluxes LADJ_STRAT (hml, DKh, 02/20/12, adj32_025)
!  (4)  Modified to read gctm.sf files at a specific iteration number
!       (nab, 03/24/12, adj32_027) !
!******************************************************************************
!
      ! References to F90 modules
      USE ADJ_ARRAYS_MOD,     ONLY : EMS_SF, ICS_SF
      USE ADJ_ARRAYS_MOD,     ONLY : NNEMS, MMSCL
      USE ADJ_ARRAYS_MOD,     ONLY : EXPAND_NAME
      USE ADJ_ARRAYS_MOD,     ONLY : PROD_SF, LOSS_SF
      USE ADJ_ARRAYS_MOD,     ONLY : NSTPL
      USE BPCH2_MOD,          ONLY : OPEN_BPCH2_FOR_READ
      USE DIRECTORY_ADJ_MOD,  ONLY : OPTDATA_DIR
      USE ERROR_MOD,          ONLY : DEBUG_MSG, ERROR_STOP
      USE FILE_MOD,           ONLY : IU_RST, IOERROR
      USE LOGICAL_ADJ_MOD,    ONLY : LICS, LADJ_EMS
      USE LOGICAL_ADJ_MOD,    ONLY : LADJ_STRAT
      USE LOGICAL_MOD,        ONLY : LPRT
      USE RESTART_MOD,        ONLY : CHECK_DIMENSIONS
      USE TIME_MOD,           ONLY : EXPAND_DATE
      USE TRACER_MOD,         ONLY : N_TRACERS

#     include "CMN_SIZE"   ! Size parameters
#     include "CMN"        ! LPRT

      ! Local Variables
      INTEGER, INTENT(IN) :: ITE
      INTEGER             :: I, IOS, J, L, M, N
      INTEGER             :: NCOUNT(NNPAR)
      REAL*4              :: TRACER(IIPAR,JJPAR,LLPAR)
      REAL*8              :: SUMTC
      CHARACTER(LEN=255)  :: FILENAME

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
      CHARACTER(LEN=20)   :: INPUT_SF_FILE

      !=================================================================
      ! READ_SF_FILE begins here!
      !=================================================================

      ! Hardwire output file for now
      INPUT_SF_FILE = 'gctm.sf.NN'

      ! Initialize some variables
      NCOUNT(:)     = 0
      TRACER(:,:,:) = 0e0

      !=================================================================
      ! Open SF file and read top-of-file header
      !=================================================================

      ! Copy input file name to a local variable
      FILENAME = TRIM( INPUT_SF_FILE )

      ! Replace NN token w/ actual value
      CALL EXPAND_NAME( FILENAME, ITE )

      ! Add OPTDATA_DIR prefix to FILENAME
      FILENAME = TRIM( OPTDATA_DIR ) // TRIM( FILENAME )

      ! can hardwire this to read a specific file from another run:
      !FILENAME = TRIM( 'opt_ics/ADJv27fi04r10/gctm.ics.16' )

      ! Echo some input to the screen
      WRITE( 6, '(a)'   ) REPEAT( '=', 79 )
      WRITE( 6, '(a,/)' ) 'S F   F I L E   I N P U T'
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( 'READ_SF_FILE: Reading ', a )

      ! Open the binary punch file for input
      CALL OPEN_BPCH2_FOR_READ( IU_RST, FILENAME )

      IF ( LICS ) THEN

         !=================================================================
         ! Read initial conditions -- store in the TRACER array
         !=================================================================
         DO N = 1, N_TRACERS
            READ( IU_RST, IOSTAT=IOS )
     &        MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180

            ! IOS < 0 is end-of-file, so exit
            IF ( IOS < 0 ) EXIT

            ! IOS > 0 is a REAL I/O error -- print error message
            IF ( IOS > 0 ) CALL IOERROR( IOS,IU_RST,'read_ics_file:4' )

            READ( IU_RST, IOSTAT=IOS )
     &           CATEGORY, NTRACER,  UNIT, ZTAU0,  ZTAU1,  RESERVED,
     &           NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
     &           NSKIP

            IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_ics_file:5')

            READ( IU_RST, IOSTAT=IOS )
     &           ( ( ( TRACER(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )

            IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_ics_file:6')

            !==============================================================
            ! Assign data from the TRACER array to the xxx_IC array.
            !==============================================================

            ! Only process observation data (i.e. aerosol and precursors)
            IF ( CATEGORY(1:8) == 'IJ-ICS-$' ) THEN

               ! Make sure array dimensions are of global size
               ! (NI=IIPAR; NJ=JJPAR, NL=LLPAR), or stop the run
               CALL CHECK_DIMENSIONS( NI, NJ, NL )

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
               DO L = 1, LLPAR
               DO J = 1, JJPAR
               DO I = 1, IIPAR
                  ICS_SF (I,J,L,N) = TRACER(I,J,L)
               ENDDO
               ENDDO
               ENDDO
!$OMP END PARALLEL DO

             ENDIF
         ENDDO

      ENDIF

      IF ( LADJ_EMS ) THEN

         !=================================================================
         ! Read emission scale factors -- store in the TRACER array
         !=================================================================
         DO N = 1, NNEMS
            READ( IU_RST, IOSTAT=IOS )
     &        MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180

            ! IOS < 0 is end-of-file, so exit
            IF ( IOS < 0 ) EXIT

            ! IOS > 0 is a REAL I/O error -- print error message
            IF ( IOS > 0 ) CALL IOERROR( IOS,IU_RST,'read_ics_file:4' )

            READ( IU_RST, IOSTAT=IOS )
     &           CATEGORY, NTRACER,  UNIT, ZTAU0,  ZTAU1,  RESERVED,
     &           NI,       NJ,       NL,   IFIRST, JFIRST, LFIRST,
     &           NSKIP

            IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_ics_file:5')

            READ( IU_RST, IOSTAT=IOS )
     &           ( ( ( TRACER(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )

            IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,'read_ics_file:6')

            !==============================================================
            ! Assign data from the TRACER array to the xxx_IC array.
            !==============================================================

            ! Only process observation data (i.e. aerosol and precursors)
            IF ( CATEGORY(1:8) == 'IJ-EMS-$' ) THEN

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, M )
               DO M = 1, MMSCL
               DO J = 1, JJPAR
               DO I = 1, IIPAR
                  EMS_SF(I,J,M,N) = TRACER(I,J,M)
               ENDDO
               ENDDO
               ENDDO
!$OMP END PARALLEL DO

             ENDIF
         ENDDO

         ! Strat prod and loss (hml)
         IF ( LADJ_STRAT ) THEN

            !=================================================================
            ! Read strat prod & loss scale factors -- store in the TRACER array
            !=================================================================
            DO N = 1, NSTPL
               READ( IU_RST, IOSTAT=IOS )
     &           MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180

               ! IOS < 0 is end-of-file, so exit
               IF ( IOS < 0 ) EXIT

               ! IOS > 0 is a REAL I/O error -- print error message
               IF ( IOS > 0 ) CALL IOERROR( IOS,IU_RST,
     &                                      'read_strat_file:4' )

               READ( IU_RST, IOSTAT=IOS )
     &              CATEGORY, NTRACER,  UNIT, ZTAU0,  ZTAU1,  RESERVED,
     &              NI,       NJ,       NL,    IFIRST, JFIRST, LFIRST,
     &              NSKIP

               IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,
     &                                       'read_strat_file:5')

               READ( IU_RST, IOSTAT=IOS )
     &              ( ( ( TRACER(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )

               IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,
     &                                       'read_strat_file:6')

               !==============================================================
               ! Assign data from the TRACER array to the xxx_STR array.
               !==============================================================

               ! Only process observation data (i.e. aerosol and precursors)

               IF ( CATEGORY(1:8) == 'IJ-STRP$' ) THEN

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, M )
                  DO M = 1, MMSCL
                  DO J = 1, JJPAR
                  DO I = 1, IIPAR
                     PROD_SF(I,J,M,N) = TRACER(I,J,M)
                  ENDDO
                  ENDDO
                  ENDDO
!$OMP END PARALLEL DO
               ENDIF
            ENDDO

            !=================================================================
            ! Read strat prod & loss scale factors -- store in the TRACER array
            !=================================================================
            DO N = 1, NSTPL
               READ( IU_RST, IOSTAT=IOS )
     &           MODELNAME, LONRES, LATRES, HALFPOLAR, CENTER180

               ! IOS < 0 is end-of-file, so exit
               IF ( IOS < 0 ) EXIT

               ! IOS > 0 is a REAL I/O error -- print error message
               IF ( IOS > 0 ) CALL IOERROR( IOS,IU_RST,
     &                                      'read_strat_file:4' )

               READ( IU_RST, IOSTAT=IOS )
     &              CATEGORY, NTRACER,  UNIT, ZTAU0,  ZTAU1,  RESERVED,
     &              NI,       NJ,       NL,    IFIRST, JFIRST, LFIRST,
     &              NSKIP

               IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,
     &                                       'read_strat_file:5')

               READ( IU_RST, IOSTAT=IOS )
     &              ( ( ( TRACER(I,J,L), I=1,NI ), J=1,NJ ), L=1,NL )

               IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_RST,
     &                                       'read_strat_file:6')

               IF ( CATEGORY(1:8) == 'IJ-STRL$' ) THEN

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, M )
                  DO M = 1, MMSCL
                  DO J = 1, JJPAR
                  DO I = 1, IIPAR
                     LOSS_SF(I,J,M,N) = TRACER(I,J,M)
                  ENDDO
                  ENDDO
                  ENDDO
!$OMP END PARALLEL DO
                ENDIF
             ENDDO
         ENDIF
      ENDIF

      ! Close file
      CLOSE( IU_RST )

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### READ_SF_FILE: read file' )

      ! Return to CALLing program
      END SUBROUTINE READ_SF_FILE_AT
!------------------------------------------------------------------------------

      SUBROUTINE MAKE_HESS_DIAG_FILE ( HINV,HMAX,NITR )
!
!******************************************************************************
!  Subroutine MAKE_HESS_DIAG_FILE creates a binary file of selected elements
!  of the approximate inverse hessian. (DKh, 05/15/07)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) HINV      : Current estimate of diagonal of inverse hessian
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
!  (2 ) Updated for adj32 (DKh, 01/11/12)
!  (3 ) Updated for adj32_027 as HINV has now dimension (HMAX) (nab, 25/03/12)
!  (4 ) Updated to include calculation of the DFS
!       and to filter out model pixels with no emissions (nab, 8/17/2012)
!******************************************************************************
!
      ! References to F90 modules
      USE ADJ_ARRAYS_MOD,    ONLY : MMSCL, NNEMS
      USE TRACER_MOD , 		 ONLY : N_TRACERS
      USE ADJ_ARRAYS_MOD,    ONLY : EXPAND_NAME
      USE ADJ_ARRAYS_MOD,    ONLY : EMS_ERROR, ICS_ERROR
      USE BPCH2_MOD
      USE DIRECTORY_ADJ_MOD, ONLY : DIAGADJ_DIR
      USE ERROR_MOD,         ONLY : DEBUG_MSG, ERROR_STOP
      USE FILE_MOD,          ONLY : IU_RST,IU_FILE,IOERROR
      USE GRID_MOD,          ONLY : GET_XOFFSET, GET_YOFFSET
      USE LOGICAL_MOD,       ONLY : LPRT
      USE LOGICAL_ADJ_MOD,   ONLY : LICS,LADJ_STRAT
      USE LOGICAL_ADJ_MOD,   ONLY : LADJ_EMS
      USE TIME_MOD,          ONLY : EXPAND_DATE, GET_TAU
      USE LOGICAL_ADJ_MOD,      ONLY : LICS, L4DVAR, LADJ_EMS
      USE ADJ_ARRAYS_MOD,       ONLY : EMS_SF_ADJ,ICS_SF_ADJ

#     include "CMN_SIZE"          ! Size parameters


      ! Arguments
      INTEGER                    :: HMAX
      REAL*8                     :: HINV(HMAX)
      INTEGER                    :: NITR
      REAL*8                     :: DFS, S2_INV,INFL

      ! Local Variables
      INTEGER                    :: I,    I0, IOS, J
      INTEGER                    :: J0, L, M, N, II, JJ
      INTEGER                    :: YYYY, MM, DD,  HH, SS
      REAL*4                     :: TRACER(IIPAR,JJPAR,LLPAR)
      REAL*4                     :: EMS_3D(IIPAR,JJPAR,MMSCL)
      REAL*4                     :: ICS_3D(IIPAR,JJPAR,LLPAR)

      CHARACTER(LEN=255)         :: FILENAME1,FILENAME2

      ! For binary punch file, version 2.0
      REAL*4                     :: LONRES, LATRES
      INTEGER, PARAMETER         :: HALFPOLAR = 1
      INTEGER, PARAMETER         :: CENTER180 = 1

      CHARACTER(LEN=20)          :: OUTPUT_GDT_FILE1,OUTPUT_GDT_FILE2
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
      ICS_3D(:,:,:) = 0d0

      ! Hardwire output file for now
      OUTPUT_GDT_FILE1 = 'gctm.posterror.NN'
      OUTPUT_GDT_FILE2 = 'gctm.diagmrm.NN'

      ! Define variables for BINARY PUNCH FILE OUTPUT
      TITLE    = 'GEOS-CHEM Adjoint File: ' //
     &           'Inverse hessian  '
      UNIT     = 'none'
      CATEGORY = 'IJ-GDE-$'
      LONRES   = DISIZE
      LATRES   = DJSIZE

      ! CALL GET_MODELNAME to return the proper model name for
      ! the given met data being used (bmy, 6/22/00)
      MODELNAME = GET_MODELNAME()

      ! Get the nested-grid offsets
      I0 = GET_XOFFSET( GLOBAL=.TRUE. )
      J0 = GET_YOFFSET( GLOBAL=.TRUE. )

      !=================================================================
      ! Open the adjoint file for output -- binary punch format
      !=================================================================

      ! Copy the output observation file name into a local variable
      FILENAME1 = TRIM( OUTPUT_GDT_FILE1 )
      FILENAME2 = TRIM( OUTPUT_GDT_FILE2 )

      ! Append the iteration number suffix to the file name
      CALL EXPAND_NAME( FILENAME1, NITR )
      CALL EXPAND_NAME( FILENAME2, NITR )

      ! Add the OPT_DATA_DIR prefix to the file name
      FILENAME1 = TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME1 )
      FILENAME2 = TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME2 )

      WRITE( 6, 100 ) TRIM( FILENAME1 )
      WRITE( 6, 100 ) TRIM( FILENAME2 )
 100  FORMAT( '     - MAKE_HESS_FILE: Writing ', a )

      ! Open checkpoint file for output
      CALL OPEN_BPCH2_FOR_WRITE( IU_RST, FILENAME1, TITLE )
      CALL OPEN_BPCH2_FOR_WRITE( IU_FILE, FILENAME2, TITLE )


      DFS = 0d0

      IF ( LADJ_STRAT ) THEN
         CALL ERROR_STOP( 'inverse hessian not supported ',
     &                    ' MAKE_HESS_FILE, inverse_mod.f')
      ELSEIF ( LICS ) THEN

           DO N = 1, N_TRACERS
             ICS_3D(:,:,:) = 0d0
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J,L, II )
            DO J = 1, JJPAR
            DO I = 1, IIPAR
            DO L = 1 , LLPAR

               II = IIMAP(I,J,L,N)
               IF ( II == 0 ) CYCLE
#if defined ( LOG_OPT )

                    S2_INV = 1d0 / (2 * LOG( ICS_ERROR(MAPN(II)) ) )
#else

                    S2_INV = 1d0 / (ICS_ERROR(MAPN(II)) )**2
#endif

               IF ( ABS( ICS_SF_ADJ(I,J,L,N)
     &               /MAXVAL(ICS_SF_ADJ(:,:,L,N)) ) < 0.01 ) THEN
                    HINV(II) = 1/S2_INV
               ENDIF

                  IF ( HINV(II) > 0 )  THEN
                     ICS_3D(I,J,L) = REAL(SQRT(HINV(II)))
                     DFS = DFS + HINV(II)*S2_INV
                  ELSE
                     print*, I, J, M, N, II
                     print*,'non positive hess diagonal:'
                     print*,HINV(II)
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
     &                  IIPAR,     JJPAR,     LLPAR,     I0+1,
     &                  J0+1,      1,         ICS_3D )


         ENDDO

      ELSEIF ( LADJ_EMS ) THEN

         !=================================================================
         ! WRITE the standard error of each optimized scaling factor
         !=================================================================
         DO N = 1, NNEMS

            !Temporarily store quantities in the TRACER array
            EMS_3D(:,:,:) = 0d0

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, M, II,S2_INV )
            DO M = 1, MMSCL
            DO J = 1, JJPAR
            DO I = 1, IIPAR

               II = IIMAP(I,J,M,N)
               IF ( II == 0 ) CYCLE
#if defined ( LOG_OPT )

                    S2_INV = 1d0 / (2 * LOG( EMS_ERROR(MAPN(II)) ) )
#else

                    S2_INV = 1d0 / ( EMS_ERROR(MAPN(II)) )**2
#endif

               IF ( ABS( EMS_SF_ADJ(I,J,M,N)
     &               /MAXVAL(EMS_SF_ADJ(:,:,M,N)) ) < 0.01 ) THEN
                    HINV(II) = 1/S2_INV
               ENDIF

                  IF ( HINV(II) > 0 )  THEN
                     EMS_3D(I,J,M) = REAL(SQRT(HINV(II)))
                     DFS = DFS + HINV(II)*S2_INV
                  ELSE
                     print*, I, J, M, N, II
                     print*,'non positive hess diagonal:'
                     print*,HINV(II)
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

!         ! Reset CATEGORY as labeling in gamap is dIFferent
!         CATEGORY = 'IJ-COREL'
!
!         !=================================================================
!         ! WRITE correlation of optimized scale factors with a particular
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
      CLOSE( IU_FILE )

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### MAKE_HESS_FILE: wrote file' )

             DFS = HMAX - DFS
             INFL = DFS/HMAX

      PRINT*,''
      PRINT*,'=============================================='
      PRINT*,'Degree of Freedom for Signal (DFS): ',DFS
      PRINT*,'Dimension (N): ',HMAX
      PRINT*,'Global average influence: DFS/N = ',INFL
      PRINT*,'=============================================='


      ! Return to CALLing program
      END SUBROUTINE MAKE_HESS_DIAG_FILE

!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------

      SUBROUTINE INIT_INV_HESSIAN(HMAX)
!
!******************************************************************************
!  Subroutine INIT_INV_HESSIAN initializes and zeros all ALLOCATABLE arrays
!
!  NOTES:
!******************************************************************************

      ! References to F90 modules
      USE ADJ_ARRAYS_MOD, ONLY : NNEMS,MMSCL
      USE TRACER_MOD , ONLY : N_TRACERS
      USE ERROR_MOD, ONLY      : ALLOC_ERR
      USE LOGICAL_ADJ_MOD,     ONLY : LICS,LADJ_EMS

#     include "CMN_SIZE"       ! Size parameters

      ! Local variables
      INTEGER                 :: AS, I, HMAX

      !=================================================================
      ! INIT_INV_HESSIAN begins here!
      !=================================================================

      IF ( LADJ_EMS) THEN
       ALLOCATE( IIMAP(IIPAR,JJPAR,MMSCL,NNEMS), STAT = AS )
       IF ( AS /= 0 ) CALL ALLOC_ERR( 'IIMAP' )
       IIMAP = 0
      ENDIF

      IF ( LICS ) THEN
       ALLOCATE( IIMAP(IIPAR,JJPAR,LLPAR,N_TRACERS), STAT = AS )
       IF ( AS /= 0 ) CALL ALLOC_ERR( 'IIMAP' )
       IIMAP = 0
      ENDIF

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

      ALLOCATE (MAPL(HMAX), STAT = AS)
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'MAPL' )
      MAPL = 0

      ALLOCATE (HINVD(HMAX), STAT = AS)
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'HINVD' )
      HINVD = 0d0

      ALLOCATE( EMS_SF_OLD(IIPAR,JJPAR,MMSCL,NNEMS), STAT = AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EMS_SF_OLD' )
      EMS_SF_OLD = 0d0

      ALLOCATE( EMS_SF_ADJ_OLD(IIPAR,JJPAR,MMSCL,NNEMS), STAT = AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'EMS_SF_ADJ_OLD' )
      EMS_SF_ADJ_OLD = 0d0

      ALLOCATE( ICS_SF_OLD(IIPAR,JJPAR,LLPAR,N_TRACERS), STAT = AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ICS_SF_OLD' )
      ICS_SF_OLD = 0d0

      ALLOCATE( ICS_SF_ADJ_OLD(IIPAR,JJPAR,LLPAR,N_TRACERS), STAT = AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ICS_SF_ADJ_OLD' )
      ICS_SF_ADJ_OLD = 0d0



      END SUBROUTINE INIT_INV_HESSIAN

!------------------------------------------------------------------------------

      ! Return to CALLing program
      SUBROUTINE CLEANUP_INV_HESSIAN
!
!******************************************************************************
!  Subroutine CLEANUP_INV_HESSIAN deALLOCATEs all previously ALLOCATEd arrays
!  for inverse_mod -- CALL at the end of the program (DKh, 01/11/12)
!
!  NOTES:
!******************************************************************************
!
      !=================================================================
      ! CLEANUP_INV_HESSIAN begins here!
      !=================================================================
      IF ( ALLOCATED( IIMAP          ) ) DEALLOCATE( IIMAP          )
      IF ( ALLOCATED( ICS_SF_OLD     ) ) DEALLOCATE( ICS_SF_OLD     )
      IF ( ALLOCATED( ICS_SF_ADJ_OLD ) ) DEALLOCATE( ICS_SF_ADJ_OLD )
      IF ( ALLOCATED( EMS_SF_OLD     ) ) DEALLOCATE( EMS_SF_OLD     )
      IF ( ALLOCATED( EMS_SF_ADJ_OLD ) ) DEALLOCATE( EMS_SF_ADJ_OLD )
      IF ( ALLOCATED( MAPI          ) ) DEALLOCATE( MAPI          )
      IF ( ALLOCATED( MAPJ          ) ) DEALLOCATE( MAPJ          )
      IF ( ALLOCATED( MAPM          ) ) DEALLOCATE( MAPM          )
      IF ( ALLOCATED( MAPL         ) ) DEALLOCATE( MAPL          )
      IF ( ALLOCATED( MAPN          ) ) DEALLOCATE( MAPN          )
      IF ( ALLOCATED( HINVD          ) ) DEALLOCATE( HINVD          )


      ! Return to CALLing program
      END SUBROUTINE CLEANUP_INV_HESSIAN

!------------------------------------------------------------------------------


      END MODULE INV_HESSIAN_LBFGS_MOD
