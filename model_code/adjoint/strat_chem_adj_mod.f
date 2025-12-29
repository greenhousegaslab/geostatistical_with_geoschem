!$ID$
!
!  Subroutine STRAT_CHEM_ADJ_MOD performs adjoint of strat chem.
!  Based on forward model routine STRAT_CHEM_MOD.
!
! !INTERFACE:
!
      MODULE STRAT_CHEM_ADJ_MOD
!
! !USES:
!
      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: DO_STRAT_CHEM_ADJ
!
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
      ! Scalars
      !REAL*8               :: DTCHEM

      ! Parameters
      !INTEGER, PARAMETER   :: NTR_GMI = 120 ! Number of species
                              ! 118 as output from GMI + NOx + Ox families

      !INTEGER, PARAMETER   :: MAX_FM  = 1 ! Max number of species in a fam
      ! Vestigial, as NOx and Ox families pre-processed, but may be useful
      ! for future uses, e.g., ClOx.

      ! Arrays
      !REAL*8,  ALLOCATABLE :: PROD(:,:,:,:)
      !REAL*8,  ALLOCATABLE :: LOSS(:,:,:,:)
      !INTEGER, ALLOCATABLE :: GMI_TO_GC(:,:)
      !INTEGER, SAVE        :: ncID_strat_rates

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement
      !=================================================================
      CONTAINS
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: DO_STRAT_CHEM_ADJ
!
! !DESCRIPTION: Function DO\_STRAT\_CHEM is the driver routine for computing
!     the simple linearized stratospheric chemistry scheme for a host of species
!     whose prod/loss rates were determined from the GMI combo model. Ozone is
!     treated using either Linoz or Synoz.
!
! !INTERFACE:
!
      SUBROUTINE DO_STRAT_CHEM_ADJ
!
! !USES:
!
      USE DAO_MOD,        ONLY : AD, CONVERT_UNITS
      USE ERROR_MOD,      ONLY : DEBUG_MSG
      USE LOGICAL_MOD,    ONLY : LLINOZ, LPRT
      USE TIME_MOD,       ONLY : GET_MONTH, TIMESTAMP_STRING
      USE TRACER_MOD,     ONLY : ITS_A_FULLCHEM_SIM, ITS_A_TAGOX_SIM
      USE TRACER_MOD,     ONLY : N_TRACERS, STT, TCVV, TRACER_MW_KG
      USE TRACERID_MOD,   ONLY : IDTOX
      USE TROPOPAUSE_MOD, ONLY : GET_MIN_TPAUSE_LEVEL, ITS_IN_THE_TROP
      ! adj_group (hml, 07/20/11)
      USE STRAT_CHEM_MOD, ONLY : PROD_0, LOSS_0
      USE STRAT_CHEM_MOD, ONLY : PROD, LOSS
      USE STRAT_CHEM_MOD, ONLY : DTCHEM
      USE STRAT_CHEM_MOD, ONLY : NSCHEM
      USE STRAT_CHEM_MOD, ONLY : Strat_TrID_GC
      USE STRAT_CHEM_MOD, ONLY : GET_RATES
      USE STRAT_CHEM_MOD, ONLY : GET_RATES_INTERP
      USE ADJ_ARRAYS_MOD, ONLY : STT_ADJ
      USE ADJ_ARRAYS_MOD, ONLY : PROD_SF,     LOSS_SF
      USE ADJ_ARRAYS_MOD, ONLY : PROD_SF_ADJ, LOSS_SF_ADJ
      USE ADJ_ARRAYS_MOD, ONLY : ID_LOSS
      USE ADJ_ARRAYS_MOD, ONLY : IFD, JFD, LFD
      USE ADJ_ARRAYS_MOD, ONLY : NSTPL
      USE TIME_MOD,       ONLY : ITS_A_NEW_MONTH
      USE BPCH2_MOD,      ONLY : GET_TAU0,        READ_BPCH2
      USE LINOZ_ADJ_MOD,  ONLY : DO_LINOZ_ADJ
      USE CHECKPOINT_MOD, ONLY : READ_BEFSTRAT_CHKFILE
      USE TIME_MOD,       ONLY : GET_NHMS
      USE TIME_MOD,       ONLY : GET_NYMD
      USE TRACER_MOD,     ONLY : STT_STRAT_TMP
      USE LOGICAL_ADJ_MOD,ONLY : LADJ_STRAT

#     include "define.h"
#     include "CMN_SIZE"
!
!EOP
!------------------------------------------------------------------------------
!
! !LOCAL VARIABLES:
!
      INTEGER, SAVE             :: LASTSEASON = -1
      INTEGER                   :: I, J, L, N, LMIN
      INTEGER                   :: IORD, JORD, KORD
      INTEGER                   :: NN, NS, NSL
      REAL*8                    :: dt, P, k, M0
      REAL*8                    :: P_ADJ, k_ADJ, M0_ADJ
      REAL*8                    :: LOSS_ADJ, PROD_ADJ
      CHARACTER(LEN=16)         :: STAMP
      INTEGER                   :: NHMS
      INTEGER                   :: NYMD


      !===============================
      ! DO_STRAT_CHEM_ADJ begins here!
      !===============================

      STAMP = TIMESTAMP_STRING()
      WRITE( 6, 100 ) STAMP
 100  FORMAT( '     - DO_STRAT_CHEM_ADJ: Strat chemistry at ', a )

      !================================================
      ! Determine the rates from disk; merge families
      !================================================

      ! Get the minimum level extent of the tropopause
      LMIN = GET_MIN_TPAUSE_LEVEL()

      ! Use ITS_A_NEW_MONTH instead, which works for forward and adjoint
      !IF ( GET_MONTH() /= LASTMONTH ) THEN
      IF ( ITS_A_NEW_MONTH() ) THEN

         WRITE(6,*) 'Getting new strat rates for month: ',GET_MONTH()

         IF ( LPRT ) CALL DEBUG_MSG( '### STRAT_CHEM_ADJ: at GET_RATES')

         ! Read rates for this month
         IF ( ITS_A_FULLCHEM_SIM() ) THEN
#if defined( GRID4x5 ) || defined( GRID2x25 )
               CALL GET_RATES( GET_MONTH() )
#else
               ! For resolutions finer than 2x2.5, nested,
               ! or otherwise exotic domains and resolutions
               CALL GET_RATES_INTERP( GET_MONTH() )
#endif

         ENDIF
      ENDIF

      IF ( LPRT )
     &   CALL DEBUG_MSG( '### STRAT_CHEM_ADJ: at DO_STRAT_CHEM_ADJ' )

      ! READING STT FROM CHECKPOINT FILE (hml, 07/31/11)
      NHMS     = GET_NHMS()
      NYMD     = GET_NYMD()
      CALL READ_BEFSTRAT_CHKFILE( NYMD, NHMS )

      WRITE(6,*) '-----------------------------------------------------'
      write(6,*) '    Doing strat chem ajdiont (STRAT_CHEM_ADJ_MOD)    '
      WRITE(6,*) '-----------------------------------------------------'

      !================================================================
      ! Full chemistry simulations
      !================================================================
      IF ( ITS_A_FULLCHEM_SIM() ) THEN

         !=============================================================
         ! Do chemical production and loss for non-ozone species for
         ! which we have explicit prod/loss rates from GMI
         !=============================================================

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N, k, P, dt, M0, NN, NS )
!$OMP+PRIVATE( k_ADJ,   P_ADJ,   M0_ADJ )
!$OMP+PRIVATE( LOSS_ADJ,     PROD_ADJ   )
!$OMP+SCHEDULE( DYNAMIC )
         DO J = 1,JJPAR
         DO I = 1,IIPAR

            DO L = LMIN,LLPAR

               IF ( ITS_IN_THE_TROP( I, J, L ) ) CYCLE

               DO N=1,NSCHEM ! Tracer index of active strat chem species
                  NN = Strat_TrID_GC(N) ! Tracer index in STT

                  ! Include something to expediate skipping past species
                  ! that we do not have strat chem for. Prob put tracer on
                  ! outermost loop.

                  ! Skip Ox; we'll always use either Linoz or Synoz
                  ! Now we will use GMI rate for Ox if LINOZ is off (hml, 10/31/11)
                  IF ( ITS_A_FULLCHEM_SIM() .and. (NN .eq. IDTOx) .and.
     &                 LLINOZ) CYCLE

                  ! adj_group:  make a version that applies scaling factors
                  ! and use this if the stratosphere adjoint ID #'s are active
                  IF ( LADJ_STRAT ) THEN
                     DO NS = 1, NSTPL

                        NSL = ID_LOSS(NS) ! same for ID_PROD(NS)

                        IF ( NN .EQ. NSL ) THEN

                        PROD(I,J,L,N) = PROD_0(I,J,L,N)
     &                                * PROD_SF(I,J,1,NS)
                        LOSS(I,J,L,N) = LOSS_0(I,J,L,N)
     &                                * LOSS_SF(I,J,1,NS)
                        ENDIF
                     ENDDO
                  ENDIF

                  ! recalculate forward values to use for adjoint code (hml)
                  dt = DTCHEM                              ! timestep [s]
                  k  = LOSS(I,J,L,N)                       ! loss freq [s-1]
                  P  = PROD(I,J,L,N) * AD(I,J,L) / TCVV(NN)! production term [kg s-1]
                  ! Use checkpointed value
                  M0 = STT_STRAT_TMP(I,J,L,NN)             ! initial mass [kg]

                  ! debug test
                  !IF ( I == IFD .and. J == JFD .and. L == LFD ) THEN
                  !   print*, ' IFD, JFD, LFD  = ', IFD, JFD, LFD
                  !   print*, NN,' STRAT TEST adj: k = ', k
                  !   print*, NN,' STRAT TEST adj: P = ', P
                  !   print*, NN,' STRAT TEST adj: M0= ', M0
                  !ENDIF

                  ! No prod or loss at all
                  IF ( k .eq. 0d0 .and. P .eq. 0d0 ) CYCLE

                  ! Simple analytic solution to dM/dt = P - kM over [0,t]
                  IF ( k .gt. 0d0 ) THEN
                     ! fwd code:
                     !STT(I,J,L,N) = M0 * exp(-k*t) + (P/k)*(1d0-exp(-k*t))
                     ! adj code:
                     M0_ADJ = STT_ADJ(I,J,L,NN) * exp(-k*dt)
                     P_ADJ  = STT_ADJ(I,J,L,NN) * (1d0 - exp(-k*dt))/k
                     k_ADJ  = STT_ADJ(I,J,L,NN)
     &                      * ( -p/(k**2) + p/(k**2)*exp(-k*dt)
     &                      + (p*dt/k)*exp(-k*dt) - dt*exp(-k*dt)*M0 )
                  ELSE
                     ! fwd code:
                     !STT(I,J,L,N) = M0 + P*t
                     ! adj code:
                     M0_ADJ = STT_ADJ(I,J,L,NN)
                     P_ADJ  = STT_ADJ(I,J,L,NN) * dt
                  ENDIF

                  ! fwd code:
                  !k = LOSS(I,J,L,N)                       ! loss freq [s-1]
                  !P = PROD(I,J,L,N) * AD(I,J,L) / TCVV(N) ! production term [kg s-1]
                  !M0 = STT(I,J,L,N)                       ! initial mass [kg]
                  ! adj code:
                  LOSS_ADJ          = K_ADJ
                  PROD_ADJ          = P_ADJ * AD(I,J,L) / TCVV(NN)
                  STT_ADJ (I,J,L,NN) = M0_ADJ

                  IF ( LADJ_STRAT ) THEN
                     DO NS = 1, NSTPL

                        NSL = ID_LOSS(NS) ! same for ID_PROD(NS)

                        IF ( NN .EQ. NSL ) THEN

                           ! fwd code:
                           !PROD(I,J,L,N) = PROD_0(I,J,L,N) * PROD_SF(I,J,1,N)
                           !LOSS(I,J,L,N) = LOSS_0(I,J,L,N) * LOSS_SF(I,J,1,N)
                           ! adj code:
                           PROD_SF_ADJ(I,J,1,NS) = PROD_SF_ADJ(I,J,1,NS)
     &                                           + PROD_0(I,J,L,N)
     &                                           * PROD_ADJ
                           LOSS_SF_ADJ(I,J,1,NS) = LOSS_SF_ADJ(I,J,1,NS)
     &                                           + LOSS_0(I,J,L,N)
     &                                           * LOSS_ADJ
                        ENDIF
                     ENDDO
                  ENDIF

                ENDDO ! N
             ENDDO ! L
          ENDDO ! I
          ENDDO ! J
!$OMP END PARALLEL DO


         !===================================
         ! Ozone
         !===================================

         ! fwd code: Put ozone in v/v
         !STT(:,:,:,IDTOX ) = STT(:,:,:,IDTOX) * TCVV( IDTOX ) / AD
         ! adj code: Put ozone back to kg
         STT_ADJ(:,:,:,IDTOX ) =
     &                         STT_ADJ(:,:,:,IDTOX) * AD / TCVV( IDTOX )

         IF ( LLINOZ ) THEN
            CALL DO_LINOZ_ADJ    ! Linoz
         ELSE
            ! must use Linoz or strat chem Ox fluxes for the adjoint
         ENDIF

         ! fwd code: Put ozone back to kg
         !STT(:,:,:,IDTOX) = STT(:,:,:,IDTOX) * AD / TCVV( IDTOX )
         ! adj code: Put ozone in v/v
         STT_ADJ(:,:,:,IDTOX) =
     &                       STT_ADJ(:,:,:,IDTOX)* TCVV( IDTOX ) / AD

      ELSE IF ( ITS_A_TAGOX_SIM() ) THEN

         ! fwd code:
         !CALL CONVERT_UNITS( 1, N_TRACERS, TCVV, AD, STT ) ! kg -> v/v
         ! adj code:
         CALL CONVERT_UNITS( 2, N_TRACERS, TCVV, AD, STT_ADJ ) ! v/v -> kg

         ! adjoint LINOZ  does not support tagged Ox simulation for now (hml, 10/05/11)
         IF ( LLINOZ ) THEN
            CALL DO_LINOZ_ADJ       ! Linoz
         ELSE
            ! must use Linoz or strat chem Ox fluxes for the adjoint
         ENDIF

         ! fwd code:
         !CALL CONVERT_UNITS( 2, N_TRACERS, TCVV, AD, STT ) ! v/v -> kg
         ! adj code:
         CALL CONVERT_UNITS( 1, N_TRACERS, TCVV, AD, STT_ADJ ) ! kg -> v/v

      ENDIF

      END SUBROUTINE DO_STRAT_CHEM_ADJ
!EOC
!------------------------------------------------------------------------------
      END MODULE STRAT_CHEM_ADJ_MOD
