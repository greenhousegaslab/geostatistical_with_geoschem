! $Id: fvdas_convect_adj_mod.f,v 1.4 2010/07/30 23:47:04 daven Exp $
      MODULE FVDAS_CONVECT_ADJ_MOD
!
!******************************************************************************
!  Module FVDAS_CONVECT_MOD contains routines (originally from NCAR) which
!  perform shallow and deep convection for the GEOS-4/fvDAS met fields.
!  These routines account for shallow and deep convection, plus updrafts
!  and downdrafts.  (pjr, dsa, bmy, 6/26/03, 1/21/04)
!
!  Module Variables:
!  ============================================================================
!  (1 ) RLXCLM   (LOGICAL) : Logical to relax column versus cloud triplet
!  (2 ) LIMCNV   (INTEGER) : Maximum CTM level for HACK convection
!  (3 ) CMFTAU   (REAL*8 ) : Characteristic adjustment time scale for HACK [s]
!  (4 ) EPS      (REAL*8 ) : A very small number [unitless]
!  (5 ) GRAV     (REAL*8 ) : Gravitational constant [m/s2]
!  (6 ) SMALLEST (REAL*8 ) : The smallest double-precision number
!  (7 ) TINYNUM  (REAL*8 ) : 2 times the machine epsilon for dble-precision
!  (8 ) TINYALT  (REAL*8 ) : arbitrary small num used in transport estimates
!
!  Module Routines:
!  ============================================================================
!  (1 ) INIT_FVDAS_CONVECT : Initializes fvDAS convection scheme
!  (2 ) FVDAS_CONVECT      : fvDAS convection routine, called from MAIN
!  (3 ) HACK_CONV          : HACK convection scheme routine
!  (4 ) ARCCONVTRAN        : Sets up fields for ZHANG/MCFARLANE convection
!  (5 ) CONVTRAN           : ZHANG/MCFARLANE convection scheme routine
!  (6 ) WHENFGT            : Test function -- not sure what this does?
!
!  GEOS-CHEM modules referenced by fvdas_convect_mod.f:
!  ============================================================================
!  (1 ) pressure_mod.f     : Module containing routines to compute P(I,J,L)
!
!  NOTES:
!  (1 ) Contains new updates for GEOS-4/fvDAS convection.  Also added OpenMP
!        parallel loop over latitudes in FVDAS_CONVECT. (swu, bmy, 1/21/04)
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables
      ! and routines from being seen outside "fvdas_convect_mod.f"
      !=================================================================

      ! Declare everything PRIVATE ...
      PRIVATE

      ! ... except routines INIT_FVDAS_CONVECT and FVDAS_CONVECT
      PUBLIC :: FVDAS_CONVECT_ADJ

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================

      ! Variables
      INTEGER            :: LIMCNV

      ! Constants
      LOGICAL, PARAMETER :: RLXCLM   = .TRUE.
      REAL*8,  PARAMETER :: CMFTAU   = 3600.d0
      REAL*8,  PARAMETER :: EPS      = 1.0d-13
      REAL*8,  PARAMETER :: GRAV     = 9.8d0
      REAL*8,  PARAMETER :: SMALLEST = TINY(1D0)
      REAL*8,  PARAMETER :: TINYALT  = 1.0d-36
      REAL*8,  PARAMETER :: TINYNUM  = 2*EPSILON(1D0)

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE ARCONVTRAN( DP,      MU,  MD,  EU,    MUG,
     &                       MDG,     DUG, EUG, EDG,   DPG,
     &                       DSUBCLD, JTG, JBG, IDEEP, LENGATH )
!
!******************************************************************************
!  Subroutine ARCONVTRAN sets up the convective transport using archived mass
!  fluxes from the ZHANG/MCFARLANE convection scheme.  The setup involves:
!    (1) Gather mass flux arrays.
!    (2) Calc the mass fluxes that are determined by mass balance.
!    (3) Determine top and bottom of convection.
!  (pjr, dsa, swu, bmy, 6/26/03, 1/21/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) DP      (REAL*8 ) : Delta pressure between interfaces [Pa]     Pa
!  (2 ) MU      (REAL*8 ) : Mass flux up                      [kg/m2/s]Pa/s
!  (3 ) MD      (REAL*8 ) : Mass flux down                    [kg/m2/s]Pa/s
!  (4 ) EU      (REAL*8 ) : Mass entraining from updraft      [1/s]    Pa/s
!
!  Arguments as Output:
!  ============================================================================
!  (5 ) MUG     (REAL*8 ) : Gathered mu (lon-alt array)
!  (6 ) MDG     (REAL*8 ) : Gathered md (lon-alt array)
!  (7 ) DUG     (REAL*8 ) : Mass detraining from updraft (lon-alt array)
!  (8 ) EUG     (REAL*8 ) : Gathered eu (lon-alt array)
!  (9 ) EDG     (REAL*8 ) : Mass entraining from downdraft (lon-alt array)
!  (10) DPG     (REAL*8 ) : Gathered .01*dp (lon-alt array)
!  (11) DSUBCLD (REAL*8 ) : Delta pressure from cloud base to sfc (lon-alt arr)
!  (12) JTG     (INTEGER) : Cloud top layer for columns undergoing conv.
!  (13) JBG     (INTEGER) : Cloud bottom layer for columns undergoing conv.
!  (14) IDEEP   (INTEGER) : Index of longitudes where deep conv. happens
!  (15) LENGATH (INTEGER) : Length of gathered arrays
!
!  NOTES:
!  (1 ) Removed NSTEP from arg list; it's not used.  Also zero arrays in order
!        to prevent them from being filled with compiler junk for latitudes
!        where no convection occurs at all. (bmy, 1/21/04)
!******************************************************************************
!
#     include "CMN_SIZE"   ! Size parameters

      ! Arguments
      INTEGER, INTENT(OUT) :: JTG(IIPAR)
      INTEGER, INTENT(OUT) :: JBG(IIPAR)
      INTEGER, INTENT(OUT) :: IDEEP(IIPAR)
      INTEGER, INTENT(OUT) :: LENGATH
      REAL*8,  INTENT(IN)  :: DP(IIPAR,LLPAR)
      REAL*8,  INTENT(IN)  :: MU(IIPAR,LLPAR)
      REAL*8,  INTENT(IN)  :: MD(IIPAR,LLPAR)
      REAL*8,  INTENT(IN)  :: EU(IIPAR,LLPAR)
      REAL*8,  INTENT(OUT) :: MUG(IIPAR,LLPAR)
      REAL*8,  INTENT(OUT) :: MDG(IIPAR,LLPAR)
      REAL*8,  INTENT(OUT) :: DUG(IIPAR,LLPAR)
      REAL*8,  INTENT(OUT) :: EUG(IIPAR,LLPAR)
      REAL*8,  INTENT(OUT) :: EDG(IIPAR,LLPAR)
      REAL*8,  INTENT(OUT) :: DPG(IIPAR,LLPAR)
      REAL*8,  INTENT(OUT) :: DSUBCLD(IIPAR)

      ! Local variables
      INTEGER              :: I, K, LENPOS
      INTEGER              :: INDEX(IIPAR)
      REAL*8               :: SUM(IIPAR)
      REAL*8               :: RDPG(IIPAR,LLPAR)

      !=================================================================
      ! ARCONVTRAN begins here!
      !=================================================================

      ! Initialize arrays
      DPG     = 0d0
      DSUBCLD = 0d0
      DUG     = 0d0
      EDG     = 0d0
      EUG     = 0d0
      JTG     = LLPAR
      JBG     = 1
      MDG     = 0d0
      MUG     = 0d0
      RDPG    = 0d0
      SUM     = 0d0

      !=================================================================
      ! First test if convection exists in the lon band I=1,IIPAR
      !=================================================================

      ! Sum all upward mass fluxes in the longitude band I=1,IIPAR
      DO K = 1, LLPAR
      DO I = 1, IIPAR
         SUM(I) = SUM(I) + MU(I,K)
      ENDDO
      ENDDO

      ! IDEEP is the index of longitudes where SUM( up mass flux ) > 0
      ! LENGATH is the # of values where SUM( up mass flux ) > 0
      CALL WHENFGT( IIPAR, SUM, 1, 0d0, IDEEP, LENGATH )

      ! Return if there is no convection the longitude band
      IF ( LENGATH == 0 ) RETURN

      !=================================================================
      ! Gather input mass fluxes in places where there is convection
      !=================================================================
      DO K = 1, LLPAR
      DO I = 1, LENGATH
         DPG(I,K)  = 0.01d0 * DP(IDEEP(I),K)              !convert Pa->hPa
         RDPG(I,K) = 1.d0 / DPG(I,K)
         MUG(I,K)  = MU(IDEEP(I),K) *  0.01d0             !convert Pa/s->hPa/s
         MDG(I,K)  = MD(IDEEP(I),K) *  0.01d0
         EUG(I,K)  = EU(IDEEP(I),K) *  0.01d0 * RDPG(I,K) !convert Pa/s->1/s
      ENDDO
      ENDDO

      !=================================================================
      ! Calc DU and ED in places where there is convection
      !=================================================================
      DO K = 1, LLPAR-1
      DO I = 1, LENGATH
         DUG(I,K) = EUG(I,K) - ( MUG(I,K) - MUG(I,K+1) ) * RDPG(I,K)
         EDG(I,K) = ( MDG(I,K) - MDG(I,K+1) ) * RDPG(I,K)
      ENDDO
      ENDDO

      DO I = 1, LENGATH
         DUG(I,LLPAR) = EUG(I,LLPAR) - MUG(I,LLPAR) * RDPG(I,LLPAR)
         EDG(I,LLPAR) = 0.0d0
      ENDDO

      DO K = 1, LLPAR
      DO I = 1, LENGATH
         IF ( DUG(I,K) < 1.d-7*EUG(I,K) ) DUG(I,K) = 0.0d0
      ENDDO
      ENDDO

      !=================================================================
      ! Find top and bottom layers with updrafts.
      !=================================================================
      DO I = 1, LENGATH
         JTG(I) = LLPAR
         JBG(I) = 1
      ENDDO

      ! Loop over altitudes
      DO K = 2, LLPAR

         ! Find places in the gathered array where MUG > 0
         CALL WHENFGT( LENGATH, MUG(:,K), 1, 0D0, INDEX, LENPOS )

         ! Compute top & bottom layers
         DO I = 1, LENPOS
            JTG(INDEX(I)) = MIN( K-1, JTG(INDEX(I)) )
            JBG(INDEX(I)) = MAX( K,   JBG(INDEX(I)) )
         ENDDO
      ENDDO

      !=================================================================
      ! Calc delta p between srfc and cloud base.
      !=================================================================
      DO I = 1, LENGATH
         DSUBCLD(I) = DPG(I,LLPAR)
      ENDDO

      DO K = LLPAR-1, 2, -1
      DO I = 1, LENGATH
         IF ( JBG(I) <= K ) THEN
            DSUBCLD(I) = DSUBCLD(I) + DPG(I,K)
         ENDIF
      ENDDO
      ENDDO

      ! Return to calling program
      END SUBROUTINE ARCONVTRAN

!------------------------------------------------------------------------------

      SUBROUTINE WHENFGT( N, ARRAY, INC, TARGET, INDEX, NVAL )
!
!******************************************************************************
!  Subroutine WHENFGT examines a 1-D vector and returns both an index array
!  of elements and the number of elements which are greater than a certain
!  target value.  This routine came with the fvDAS convection code, we just
!  cleaned it up and added comments. (swu, bmy, 1/21/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) N      (INTEGER) : Number of elements in ARRAY
!  (2 ) ARRAY  (REAL*8 ) : 1-D vector to be examined
!  (3 ) INC    (INTEGER) : Increment for stepping thru ARRAY
!  (4 ) TARGET (REAL*8 ) : Value that ARRAY will be tested against
!
!  Arguments as Output:
!  ============================================================================
!  (5 ) INDEX  (INTEGER) : Array  of places where ARRAY(I) > TARGET
!  (6 ) NVAL   (INTEGER) : Number of places where ARRAY(I) > TARGET
!
!  NOTES:
!  (1 ) Updated comments.  Now use F90 style declarations.  (bmy, 1/21/04)
!******************************************************************************
!
      ! Arguments
      INTEGER, INTENT(IN)  :: N, INC
      REAL*8,  INTENT(IN)  :: ARRAY(N), TARGET
      INTEGER, INTENT(OUT) :: INDEX(N), NVAL

      ! Local variables
      INTEGER              :: I, INA

      !=================================================================
      ! WHENFGT begins here!
      !=================================================================

      ! Initialize
      INA      = 1
      NVAL     = 0
      INDEX(:) = 0

      ! Loop thru the array
      DO I = 1, N

         ! If the current element of ARRAY is greater than TARGET,
         ! then increment NVAL and save the element # in INDEX
         IF ( ARRAY(INA) > TARGET ) THEN
	    NVAL        = NVAL + 1
	    INDEX(NVAL) = I
         ENDIF

         ! Skip ahead by INC elements
         INA = INA + INC
      ENDDO

      ! Return to calling program
      END SUBROUTINE WHENFGT

!------------------------------------------------------------------------------

      SUBROUTINE FVDAS_CONVECT_ADJ( TDT,   NTRACE, Q,   RPDEL, ETA,
     &                          BETA,  MU,     MD,    EU,    DP,
     &                          NSTEP, FRACIS, TCVV,  INDEXSOL, ADQ )
!
!******************************************************************************

#     include "CMN_SIZE"     ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN)    :: NSTEP, NTRACE
      INTEGER, INTENT(IN)    :: INDEXSOL(NTRACE)
      REAL*8,  INTENT(IN)    :: TDT
      !  CHK_STT_CONV  is only REAL*4 (dkh, 09/15/08)
      !REAL*8,  INTENT(INOUT) :: Q(IIPAR,JJPAR,LLPAR,NTRACE)
      REAL*4,  INTENT(INOUT) :: Q(IIPAR,JJPAR,LLPAR,NTRACE)
      REAL*8,  INTENT(INOUT) :: ADQ(IIPAR,JJPAR,LLPAR,NTRACE)
      REAL*8,  INTENT(IN)    :: RPDEL(IIPAR,JJPAR,LLPAR)
      REAL*8,  INTENT(IN)    :: ETA(IIPAR,JJPAR,LLPAR)
      REAL*8,  INTENT(IN)    :: BETA(IIPAR,JJPAR,LLPAR)
      REAL*8,  INTENT(IN)    :: MU(IIPAR,JJPAR,LLPAR)
      REAL*8,  INTENT(IN)    :: MD(IIPAR,JJPAR,LLPAR)
      REAL*8,  INTENT(IN)    :: EU(IIPAR,JJPAR,LLPAR)
      REAL*8,  INTENT(IN)    :: DP(IIPAR,JJPAR,LLPAR)
      REAL*8,  INTENT(IN)    :: FRACIS(IIPAR,JJPAR,LLPAR,NTRACE)
      REAL*8,  INTENT(IN)    :: TCVV(NTRACE)

      ! Local variables
      INTEGER                :: I, J, L, N
      INTEGER                :: LENGATH, ISTEP
      INTEGER                :: JT(IIPAR)
      INTEGER                :: MX(IIPAR)
      INTEGER                :: IDEEP(IIPAR)
      REAL*8                 :: DSUBCLD(IIPAR)
      REAL*8                 :: DPG(IIPAR,LLPAR)
      REAL*8                 :: DUG(IIPAR,LLPAR)
      REAL*8                 :: EDG(IIPAR,LLPAR)
      REAL*8                 :: EUG(IIPAR,LLPAR)
      REAL*8                 :: MDG(IIPAR,LLPAR)
      REAL*8                 :: MUG(IIPAR,LLPAR)
      REAL*8                 :: QTMP(IIPAR,LLPAR,NTRACE)
      REAL*8                 :: FTMP(IIPAR,LLPAR,NTRACE)

      integer ip1,ip2,ip3
      double precision adqtmp(iipar,llpar,ntrace)

c!==========================================================

! mak has reinstated OMP statements here, but still sees big
! differences between parallel and sequential. Need to investigate
! further (dkh, 01/06/10)
! BUG FIX:  also declase QTMP and FTMP thread private (dkh, 07/30/10).
! This fixes the discrepancies noted in comment above.
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N, ISTEP, ADQTMP, MUG, MDG )
!$OMP+PRIVATE( DUG, EUG, EDG, DPG, DSUBCLD, JT, MX, IDEEP, LENGATH )
!$OMP+PRIVATE( QTMP, FTMP )
!$OMP+SCHEDULE( DYNAMIC )

      ! Loop over latitudes
      DO J = JJPAR,1,-1

         ! Save lat slices of Q & FRACIS into QTMP & FTMP
         DO N = 1, NTRACE
         DO L = 1, LLPAR
         DO I = 1, IIPAR
            ADQTMP(I,L,N) = ADQ(I,J,L,N)
            QTMP(I,L,N)   = Q(I,J,L,N)
            FTMP(I,L,N)   = FRACIS(I,J,L,N)
         ENDDO
         ENDDO
         ENDDO

         ! Gather mass flux arrays, compute mass fluxes, and determine top
         ! and bottom of Z&M convection.  LENGATH = # of longitudes in the
         ! band I=1,IIPAR where deep convection happens at latitude J.
         CALL ARCONVTRAN( DP(:,J,:),  MU(:,J,:),  MD(:,J,:),
     &                    EU(:,J,:),  MUG,        MDG,
     &                    DUG,        EUG,        EDG,
     &                    DPG,        DSUBCLD,    JT,
     &                    MX,         IDEEP,      LENGATH )

         ! Loop over internal convection timestep

         DO ISTEP = NSTEP, 1, -1

            CALL HACK_CONV_ADJ( TDT, RPDEL(:,J,:), ETA(:,J,:),
     &           BETA(:,J,:), NTRACE, QTMP, ADQTMP )


            IF ( LENGATH > 0 ) THEN

               ! Only call CONVTRAN where convection happens
               ! (i.e. at latitudes where LENGATH > 0)
               CALL CONVTRAN_ADJ( NTRACE,  QTMP,    MUG,      MDG,
     &                            DUG,     EUG,     EDG,      DPG,
     &                            DSUBCLD, JT,      MX,       IDEEP,
     &                            1,       LENGATH, NSTEP,    0.5D0*TDT,
     &                            FTMP,    TCVV,    INDEXSOL, J,
     &                            ADQTMP )
            ENDIF

         ENDDO

         ! Save latitude slice QTMP back into global Q array
         DO N = 1, NTRACE
         DO L = 1, LLPAR
         DO I = 1, IIPAR
            ADQ(I,J,L,N) = ADQTMP(I,L,N)
         ENDDO
         ENDDO
         ENDDO

      ENDDO
!$OMP END PARALLEL DO

c!==========================================================

      ! Return to calling program
      END SUBROUTINE FVDAS_CONVECT_ADJ

!-----------------------------------------------------------------------
C
      subroutine hack_conv_adj( tdt, rpdel, eta, beta, ntrace, q, adq )
C***************************************************************
C***************************************************************
C** This routine was generated by the                         **
C** Tangent linear and Adjoint Model Compiler,  TAMC 5.3.2    **
C***************************************************************
C***************************************************************
C==============================================
C all entries are defined explicitly
C==============================================
      implicit none

#     include "CMN_SIZE"      ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN)    :: NTRACE
      REAL*8,  INTENT(IN)    :: TDT
      REAL*8,  INTENT(IN)    :: RPDEL(IIPAR,LLPAR)
      REAL*8,  INTENT(IN)    :: ETA(IIPAR,LLPAR)
      REAL*8,  INTENT(IN)    :: BETA(IIPAR,LLPAR)
      REAL*8,  INTENT(IN) :: Q(IIPAR,LLPAR,NTRACE)
      REAL*8,  INTENT(INOUT) :: ADQ(IIPAR,LLPAR,NTRACE)

C==============================================
C define local variables
C==============================================
      double precision tmp(iipar,llpar,ntrace)
      double precision adbotflx
      double precision adcmrc(iipar)
      double precision adcmrh(iipar,llpar+1)
      double precision addcmr1(iipar)
      double precision addcmr2(iipar)
      double precision addcmr3(iipar)
      double precision adefac1
      double precision adefac2
      double precision adefac3
      double precision adjfac
      double precision adt1
      double precision adtopflx
      double precision botflx
      double precision cmrc(iipar)
      double precision cmrh(iipar,llpar+1)
      double precision dcmr1(iipar)
      double precision dcmr2(iipar)
      double precision dcmr3(iipar)
      double precision efac1
      double precision efac2
      double precision efac3
      double precision etagdt(iipar)
      integer i
      integer ii
      integer ii2
      integer indx1(iipar)
      integer indx1h(iipar)
      integer ip1
      integer ip2
      integer k
      integer k2
      integer len1
      integer m
      double precision t1
      double precision temp
      double precision topflx

C----------------------------------------------
C RESET LOCAL ADJOINT VARIABLES
C----------------------------------------------
      adbotflx = 0.d0
      do ip1 = 1, iipar
        adcmrc(ip1) = 0.d0
      end do
      do ip2 = 1, llpar+1
        do ip1 = 1, iipar
          adcmrh(ip1,ip2) = 0.d0
        end do
      end do
      do ip1 = 1, iipar
        addcmr1(ip1) = 0.d0
      end do
      do ip1 = 1, iipar
        addcmr2(ip1) = 0.d0
      end do
      do ip1 = 1, iipar
        addcmr3(ip1) = 0.d0
      end do
      adefac1 = 0.d0
      adefac2 = 0.d0
      adefac3 = 0.d0
      adt1 = 0.d0
      adtopflx = 0.d0

C----------------------------------------------
C ROUTINE BODY
C----------------------------------------------
C----------------------------------------------
C FUNCTION AND TAPE COMPUTATIONS
C----------------------------------------------
      if (tdt .gt. cmftau) then
        temp = tdt
      else
        temp = cmftau
      endif
      if (rlxclm) then
        adjfac = tdt/temp
      else
        adjfac = 1.d0
      endif

C----------------------------------------------
C ADJOINT COMPUTATIONS
C----------------------------------------------
      do k = limcnv+1, llpar-1
        len1 = 0
        do i = 1, iipar
          if (eta(i,k) .ne. 0.) then
            etagdt(i) = eta(i,k)*grav*tdt*0.01d0
            len1 = len1+1
            indx1(len1) = i
         else
            etagdt(i) = 0.d0
          endif
        end do

        if (len1 .le. 0) then
        else
          do m = ntrace, 1, -1
            do ii = len1, 1, -1
              i = indx1(ii)
              if (q(i,k+1,m) .lt. 0. .or. q(i,k,m) .lt. 0. .or. q(i,k-1,
     $m) .lt. 0.) then
              else
                cmrh(i,k) = 0.5d0*(q(i,k-1,m)+q(i,k,m))
                cmrh(i,k+1) = 0.5d0*(q(i,k,m)+q(i,k+1,m))
                cmrc(i) = q(i,k+1,m)
                botflx = etagdt(i)*(cmrc(i)-cmrh(i,k+1))*adjfac
                topflx = beta(i,k)*etagdt(i)*(cmrc(i)-cmrh(i,k))*adjfac
                dcmr1(i) = -(botflx*rpdel(i,k+1))
                efac1 = 1.d0
                efac2 = 1.d0
                efac3 = 1.d0
                if (q(i,k+1,m)+dcmr1(i) .lt. 0.) then
                  t1 = q(i,k+1,m)/dcmr1(i)
                  if (t1 .lt. 0.) then
                    t1 = -t1
                  endif
                  t1 = t1-eps
                  if (tinyalt .gt. t1) then
                    efac1 = tinyalt
                  else
                    efac1 = t1
                  endif
                endif
                if (efac1 .eq. tinyalt .or. efac1 .gt. 1.) then
                  efac1 = 0.d0
                endif
                dcmr2(i) = (efac1*botflx-topflx)*rpdel(i,k)
                if (q(i,k,m)+dcmr2(i) .lt. 0.) then
                  t1 = q(i,k,m)/dcmr2(i)
                  if (t1 .lt. 0.) then
                    t1 = -t1
                  endif
                  t1 = t1-eps
                  if (tinyalt .gt. t1) then
                    efac2 = tinyalt
                  else
                    efac2 = t1
                  endif
                endif
                if (efac2 .eq. tinyalt .or. efac2 .gt. 1.) then
                  efac2 = 0.d0
                endif
                dcmr3(i) = efac2*topflx*rpdel(i,k-1)
                if (q(i,k-1,m)+dcmr3(i) .lt. 0.) then
                  t1 = q(i,k-1,m)/dcmr3(i)
                  if (t1 .lt. 0.) then
                    t1 = -t1
                  endif
                  t1 = t1-eps
                  if (tinyalt .gt. t1) then
                    efac3 = tinyalt
                  else
                    efac3 = t1
                  endif
                endif
                if (efac3 .eq. tinyalt .or. efac3 .gt. 1.) then
                  efac3 = 0.d0
                endif
                if (efac2 .gt. efac3) then
                  efac3 = efac2
                endif
                addcmr3(i) = addcmr3(i)+adq(i,k-1,m)
                addcmr2(i) = addcmr2(i)+adq(i,k,m)
                addcmr1(i) = addcmr1(i)+adq(i,k+1,m)
                adefac3 = adefac3+addcmr3(i)*topflx*rpdel(i,k-1)
                adtopflx = adtopflx+addcmr3(i)*efac3*rpdel(i,k-1)
                addcmr3(i) = 0.d0
                adbotflx = adbotflx+addcmr2(i)*efac1*rpdel(i,k)
                adefac1 = adefac1+addcmr2(i)*botflx*rpdel(i,k)
                adefac3 = adefac3-addcmr2(i)*topflx*rpdel(i,k)
                adtopflx = adtopflx-addcmr2(i)*efac3*rpdel(i,k)
                addcmr2(i) = 0.d0
                efac3 = 1.d0
                if (q(i,k-1,m)+dcmr3(i) .lt. 0.) then
                  t1 = q(i,k-1,m)/dcmr3(i)
                  if (t1 .lt. 0.) then
                    t1 = -t1
                  endif
                  t1 = t1-eps
                  if (tinyalt .gt. t1) then
                    efac3 = tinyalt
                  else
                    efac3 = t1
                  endif
                endif
                if (efac3 .eq. tinyalt .or. efac3 .gt. 1.) then
                  efac3 = 0.d0
                endif
                if (efac2 .gt. efac3) then
                  adefac2 = adefac2+adefac3
                  adefac3 = 0.d0
                endif
                efac3 = 1.d0
                if (q(i,k-1,m)+dcmr3(i) .lt. 0.) then
                  t1 = q(i,k-1,m)/dcmr3(i)
                  if (t1 .lt. 0.) then
                    t1 = -t1
                  endif
                  t1 = t1-eps
                  if (tinyalt .gt. t1) then
                    efac3 = tinyalt
                  else
                    efac3 = t1
                  endif
                endif
                if (efac3 .eq. tinyalt .or. efac3 .gt. 1.) then
                  adefac3 = 0.d0
                endif
                if (q(i,k-1,m)+dcmr3(i) .lt. 0.) then
                  if (tinyalt .gt. t1) then
                    adefac3 = 0.d0
                  else
                    adt1 = adt1+adefac3
                    adefac3 = 0.d0
                  endif
                  t1 = q(i,k-1,m)/dcmr3(i)
                  if (t1 .lt. 0.) then
                    adt1 = -adt1
                  endif
                  addcmr3(i) = addcmr3(i)-adt1*(q(i,k-1,m)/(dcmr3(i)*
     $dcmr3(i)))
                  adq(i,k-1,m) = adq(i,k-1,m)+adt1/dcmr3(i)
                  adt1 = 0.d0
                endif
                adefac2 = adefac2+addcmr3(i)*topflx*rpdel(i,k-1)
                adtopflx = adtopflx+addcmr3(i)*efac2*rpdel(i,k-1)
                addcmr3(i) = 0.d0
                adbotflx = adbotflx+addcmr2(i)*efac1*rpdel(i,k)
                adefac1 = adefac1+addcmr2(i)*botflx*rpdel(i,k)
                adefac2 = adefac2-addcmr2(i)*topflx*rpdel(i,k)
                adtopflx = adtopflx-addcmr2(i)*efac2*rpdel(i,k)
                addcmr2(i) = 0.d0
                efac2 = 1.d0
                if (q(i,k,m)+dcmr2(i) .lt. 0.) then
                  t1 = q(i,k,m)/dcmr2(i)
                  if (t1 .lt. 0.) then
                    t1 = -t1
                  endif
                  t1 = t1-eps
                  if (tinyalt .gt. t1) then
                    efac2 = tinyalt
                  else
                    efac2 = t1
                  endif
                endif
                if (efac2 .eq. tinyalt .or. efac2 .gt. 1.) then
                  adefac2 = 0.d0
                endif
                if (q(i,k,m)+dcmr2(i) .lt. 0.) then
                  if (tinyalt .gt. t1) then
                    adefac2 = 0.d0
                  else
                    adt1 = adt1+adefac2
                    adefac2 = 0.d0
                  endif
                  t1 = q(i,k,m)/dcmr2(i)
                  if (t1 .lt. 0.) then
                    adt1 = -adt1
                  endif
                  addcmr2(i) = addcmr2(i)-adt1*(q(i,k,m)/(dcmr2(i)*
     $dcmr2(i)))
                  adq(i,k,m) = adq(i,k,m)+adt1/dcmr2(i)
                  adt1 = 0.d0
                endif
                adbotflx = adbotflx+addcmr2(i)*efac1*rpdel(i,k)
                adefac1 = adefac1+addcmr2(i)*botflx*rpdel(i,k)
                adtopflx = adtopflx-addcmr2(i)*rpdel(i,k)
                addcmr2(i) = 0.d0
                adbotflx = adbotflx-addcmr1(i)*efac1*rpdel(i,k+1)
                adefac1 = adefac1-addcmr1(i)*botflx*rpdel(i,k+1)
                addcmr1(i) = 0.d0
                efac1 = 1.d0
                if (q(i,k+1,m)+dcmr1(i) .lt. 0.) then
                  t1 = q(i,k+1,m)/dcmr1(i)
                  if (t1 .lt. 0.) then
                    t1 = -t1
                  endif
                  t1 = t1-eps
                  if (tinyalt .gt. t1) then
                    efac1 = tinyalt
                  else
                    efac1 = t1
                  endif
                endif
                if (efac1 .eq. tinyalt .or. efac1 .gt. 1.) then
                  adefac1 = 0.d0
                endif
                if (q(i,k+1,m)+dcmr1(i) .lt. 0.) then
                  if (tinyalt .gt. t1) then
                    adefac1 = 0.d0
                  else
                    adt1 = adt1+adefac1
                    adefac1 = 0.d0
                  endif
                  t1 = q(i,k+1,m)/dcmr1(i)
                  if (t1 .lt. 0.) then
                    adt1 = -adt1
                  endif
                  addcmr1(i) = addcmr1(i)-adt1*(q(i,k+1,m)/(dcmr1(i)*
     $dcmr1(i)))
                  adq(i,k+1,m) = adq(i,k+1,m)+adt1/dcmr1(i)
                  adt1 = 0.d0
                endif
                adefac3 = 0.d0
                adefac2 = 0.d0
                adefac1 = 0.d0
                adbotflx = adbotflx-addcmr1(i)*rpdel(i,k+1)
                addcmr1(i) = 0.d0
                adcmrc(i) = adcmrc(i)+adtopflx*beta(i,k)*etagdt(i)*
     $adjfac
                adcmrh(i,k) = adcmrh(i,k)-adtopflx*beta(i,k)*etagdt(i)*
     $adjfac
                adtopflx = 0.d0
                adcmrc(i) = adcmrc(i)+adbotflx*etagdt(i)*adjfac
                adcmrh(i,k+1) = adcmrh(i,k+1)-adbotflx*etagdt(i)*adjfac
                adbotflx = 0.d0
                adq(i,k+1,m) = adq(i,k+1,m)+adcmrc(i)
                adcmrc(i) = 0.d0
                adq(i,k+1,m) = adq(i,k+1,m)+0.5d0*adcmrh(i,k+1)
                adq(i,k,m) = adq(i,k,m)+0.5d0*adcmrh(i,k+1)
                adcmrh(i,k+1) = 0.d0
                adq(i,k-1,m) = adq(i,k-1,m)+0.5d0*adcmrh(i,k)
                adq(i,k,m) = adq(i,k,m)+0.5d0*adcmrh(i,k)
                adcmrh(i,k) = 0.d0
              endif
            end do
          end do
        endif
      end do

      end subroutine hack_conv_adj

!------------------------------------------------------------------------------

      SUBROUTINE CONVTRAN_ADJ( NTRACE, Q,      MU,   MD,       DU,
     &     EU,     ED,     DP,   DSUBCLD,  JT,
     &     MX,     IDEEP,  IL1G, IL2G,     NSTEP,
     &     DELT,   FRACIS, TCVV, INDEXSOL, LATI_INDEX,
     &     ADQ)
!
!******************************************************************************
!
#     include "CMN_SIZE"     ! Size parameters
#     include "CMN_DIAG"     ! ND38, LD38

      ! Arguments
      INTEGER, INTENT(IN)    :: NTRACE
      INTEGER, INTENT(IN)    :: JT(IIPAR)
      INTEGER, INTENT(IN)    :: MX(IIPAR)
      INTEGER, INTENT(IN)    :: IDEEP(IIPAR)
      INTEGER, INTENT(IN)    :: IL1G
      INTEGER, INTENT(IN)    :: IL2G
      INTEGER, INTENT(IN)    :: NSTEP
      REAL*8,  INTENT(IN) :: Q(IIPAR,LLPAR,NTRACE)
      REAL*8,  INTENT(INOUT) :: ADQ(IIPAR,LLPAR,NTRACE)
      REAL*8,  INTENT(IN)    :: MU(IIPAR,LLPAR)
      REAL*8,  INTENT(IN)    :: MD(IIPAR,LLPAR)
      REAL*8,  INTENT(IN)    :: DU(IIPAR,LLPAR)
      REAL*8,  INTENT(IN)    :: EU(IIPAR,LLPAR)
      REAL*8,  INTENT(IN)    :: ED(IIPAR,LLPAR)
      REAL*8,  INTENT(IN)    :: DP(IIPAR,LLPAR)
      REAL*8,  INTENT(IN)    :: DSUBCLD(IIPAR)
      REAL*8,  INTENT(IN)    :: DELT
      REAL*8,  INTENT(IN)    :: FRACIS(IIPAR,LLPAR,NTRACE)
      REAL*8,  INTENT(IN)    :: TCVV(NTRACE)
      INTEGER, INTENT(IN)    :: INDEXSOL(NTRACE)
      INTEGER, INTENT(IN)    :: LATI_INDEX

      ! Local variables
      INTEGER                :: I,     K,      KBM,     KK,     KKP1
      INTEGER                :: KM1,   KP1,    KTM,     M,      ISTEP
      INTEGER                :: II,    JJ,     LL,      NN
      REAL*8                 :: CABV,  CBEL,   CDIFR,   CD2,    DENOM
      REAL*8                 :: SMALL, MBSTH,  MUPDUDP, MINC,   MAXC
      REAL*8                 :: QN,    FLUXIN, FLUXOUT, NETFLUX
      REAL*8                 :: CHAT(IIPAR,LLPAR)
      REAL*8                 :: COND(IIPAR,LLPAR)
      REAL*8                 :: CMIX(IIPAR,LLPAR)
      REAL*8                 :: FISG(IIPAR,LLPAR)
      REAL*8                 :: CONU(IIPAR,LLPAR)
      REAL*8                 :: DCONDT(IIPAR,LLPAR)

      integer ip1,ip2,m1
      real*8 adcabv
      real*8 adcbel
      real*8 adcdifr
      real*8 adchat(iipar,llpar)
      real*8 adcmix(iipar,llpar)
      real*8 adcond(iipar,llpar)
      real*8 adconu(iipar,llpar)
      real*8 addcondt(iipar,llpar)
      real*8 adfluxin
      real*8 adfluxout
      real*8 admaxc
      real*8 adminc
      real*8 adnetflux
      real*8 adqn
      real*8 adtmp
      real*8 tmp
      double precision adconuh
      double precision adcondh
      double precision qtmp(iipar,llpar,ntrace)

C----------------------------------------------
C RESET LOCAL ADJOINT VARIABLES
C----------------------------------------------
      adcabv = 0.
      adcbel = 0.
      adcdifr = 0.
      do ip2 = 1, llpar
        do ip1 = 1, iipar
          adchat(ip1,ip2) = 0.
        end do
      end do
      do ip2 = 1, llpar
        do ip1 = 1, iipar
          adcmix(ip1,ip2) = 0.
        end do
      end do
      do ip2 = 1, llpar
        do ip1 = 1, iipar
          adcond(ip1,ip2) = 0.
        end do
      end do
      do ip2 = 1, llpar
        do ip1 = 1, iipar
          adconu(ip1,ip2) = 0.
        end do
      end do
      do ip2 = 1, llpar
        do ip1 = 1, iipar
          addcondt(ip1,ip2) = 0.
        end do
      end do
      adfluxin = 0.
      adfluxout = 0.
      admaxc = 0.
      adminc = 0.
      adnetflux = 0.
      adqn = 0.
      adtmp = 0.

      !=================================================================
      ! CONVTRAN begins here!
      !=================================================================

      ! A small number
      SMALL = 1.d-36

      ! Threshold below which we treat the mass fluxes as zero (in mb/s)
      MBSTH = 1.d-15

      !=================================================================
      ! Find the highest level top and bottom levels of convection
      !=================================================================
      KTM = LLPAR
      KBM = LLPAR
      DO I = IL1G, IL2G
         KTM = MIN( KTM, JT(I) )
         KBM = MIN( KBM, MX(I) )
      ENDDO

      ! Loop ever each tracer
      DO M = NTRACE,1,-1
         do k = 1, llpar
            do i = il1g, il2g
               cmix(i,k) = q(ideep(i),k,m)
               if (cmix(i,k) .lt. 4.d0*smallest) then
                  cmix(i,k) = 0.d0
               endif
               fisg(i,k) = fracis(ideep(i),k,m)
            end do
         end do
         do k = 1, llpar
            KM1 = MAX( 1, K-1 )
            do i = il1g, il2g
               MINC = MIN( CMIX(I,KM1), CMIX(I,K) )
              MAXC = MAX( CMIX(I,KM1), CMIX(I,K) )

              IF ( MINC < 0d0 ) THEN
                 CDIFR = 0.d0
              ELSE
                 CDIFR = ABS( CMIX(I,K)-CMIX(I,KM1) ) / MAX(MAXC,SMALL)
              ENDIF

              IF ( CDIFR > 1.d-6 ) THEN

                  ! If the two layers differ significantly.
                  ! use a geometric averaging procedure
                 CABV = MAX( CMIX(I,KM1), MAXC*TINYNUM, SMALLEST )
                 CBEL = MAX( CMIX(I,K),   MAXC*TINYNUM, SMALLEST )

                 CHAT(I,K) = LOG( CABV / CBEL)
     &                /   ( CABV - CBEL)
     &                *     CABV * CBEL

              ELSE

                  ! Small diff, so just arithmetic mean
                 CHAT(I,K) = 0.5d0 * ( CMIX(I,K) + CMIX(I,KM1) )
              ENDIF

            conu(i,k) = chat(i,k)
            cond(i,k) = chat(i,k)
            dcondt(i,k) = 0.d0
          end do
        end do
        k = 2
        km1 = 1
        kk = llpar
        do i = il1g, il2g
          mupdudp = mu(i,kk)+du(i,kk)*dp(i,kk)
          if (mupdudp .gt. mbsth) then
            conu(i,kk) = eu(i,kk)*cmix(i,kk)*dp(i,kk)/mupdudp
          endif
          if (md(i,k) .lt. (-mbsth)) then
            cond(i,k) = (-(ed(i,km1)*cmix(i,km1)*dp(i,km1)))/md(i,k)
          endif
        end do
        do kk = llpar-1, 1, -1
          if (llpar .gt. kk+1) then
            kkp1 = kk+1
          else
            kkp1 = llpar
          endif
          do i = il1g, il2g
            mupdudp = mu(i,kk)+du(i,kk)*dp(i,kk)
            if (mupdudp .gt. mbsth) then
              conu(i,kk) = (mu(i,kkp1)*conu(i,kkp1)*fisg(i,kk)+eu(i,kk)*
     $cmix(i,kk)*dp(i,kk))/mupdudp
            endif
          end do
        end do
        do k = 3, llpar
           KM1 = MAX( 1, K-1 )
          do i = il1g, il2g
            if (md(i,k) .lt. (-mbsth)) then
              cond(i,k) = (md(i,km1)*cond(i,km1)-ed(i,km1)*cmix(i,km1)*
     $dp(i,km1))/md(i,k)
            endif
          end do
        end do
        do k = ktm, llpar
            KM1 = MAX( 1,     K-1 )
            KP1 = MIN( LLPAR, K+1 )
            do i = il1g, il2g
            fluxin = mu(i,kp1)*conu(i,kp1)*fisg(i,k)+(mu(i,k)+md(i,k))*
     $cmix(i,km1)-md(i,k)*cond(i,k)
            fluxout = mu(i,k)*conu(i,k)+(mu(i,kp1)+md(i,kp1))*cmix(i,k)-
     $md(i,kp1)*cond(i,kp1)
            netflux = fluxin-fluxout

            IF ( ABS(NETFLUX) < MAX(FLUXIN,FLUXOUT)*TINYNUM) THEN
               NETFLUX = 0.D0
            ENDIF

            dcondt(i,k) = netflux/dp(i,k)
          end do
        end do
        do k = kbm, llpar
          if (k-1 .gt. 1) then
            km1 = k-1
          else
            km1 = 1
          endif
          do i = il1g, il2g
            if (k .eq. mx(i)) then
              fluxin = (mu(i,k)+md(i,k))*cmix(i,km1)-md(i,k)*cond(i,k)
              fluxout = mu(i,k)*conu(i,k)
              netflux = fluxin-fluxout

              IF ( ABS(NETFLUX) < MAX(FLUXIN,FLUXOUT)*TINYNUM) THEN
                 NETFLUX = 0.D0
              ENDIF
              dcondt(i,k) = netflux/dp(i,k)
           else if (k .gt. mx(i)) then
              dcondt(i,k) = 0.d0
           endif
          end do
        end do
        do k = 1, llpar
          adqn = 0.
          do i = il1g, il2g
            adqn = 0.
            qn = cmix(i,k)+dcondt(i,k)*delt
            adqn = adqn+adq(ideep(i),k,m)
            adq(ideep(i),k,m) = 0.d0
            if (qn .lt. 0.d0) then
              adqn = 0.
            endif
            adcmix(i,k) = adcmix(i,k)+adqn
            addcondt(i,k) = addcondt(i,k)+adqn*delt
            adqn = 0.
          end do
        end do
        do k = llpar, kbm, -1
            KM1 = MAX( 1, K-1 )
          do i = il2g, il1g, -1
            if (k .eq. mx(i)) then
              fluxin = (mu(i,k)+md(i,k))*cmix(i,km1)-md(i,k)*cond(i,k)
              fluxout = mu(i,k)*conu(i,k)
              netflux = fluxin-fluxout
              adnetflux = adnetflux+addcondt(i,k)/dp(i,k)
              addcondt(i,k) = 0.
              if (fluxin .gt. fluxout) then
                if (netflux .lt. fluxin*tinynum) then
                  adnetflux = 0.
                endif
              else
                if ((-netflux) .lt. fluxout*tinynum) then
                  adnetflux = 0.
                endif
              endif
              adfluxin = adfluxin+adnetflux
              adfluxout = adfluxout-adnetflux
              adnetflux = 0.
              adconu(i,k) = adconu(i,k)+adfluxout*mu(i,k)
              adfluxout = 0.
              adcmix(i,km1) = adcmix(i,km1)+adfluxin*(mu(i,k)+md(i,k))
              adcond(i,k) = adcond(i,k)-adfluxin*md(i,k)
              adfluxin = 0.
            else if (k .gt. mx(i)) then
              addcondt(i,k) = 0.
            endif
          end do
        end do
        do k = ktm, llpar
          adfluxin = 0.
          adfluxout = 0.
          adnetflux = 0.
          KM1 = MAX( 1,     K-1 )
          KP1 = MIN( LLPAR, K+1 )
          do i = il1g, il2g
            adfluxin = 0.
            adfluxout = 0.
            adnetflux = 0.
            fluxin = mu(i,kp1)*conu(i,kp1)*fisg(i,k)+(mu(i,k)+md(i,k))*
     $cmix(i,km1)-md(i,k)*cond(i,k)
            fluxout = mu(i,k)*conu(i,k)+(mu(i,kp1)+md(i,kp1))*cmix(i,k)-
     $md(i,kp1)*cond(i,kp1)
            netflux = fluxin-fluxout
            adnetflux = adnetflux+addcondt(i,k)/dp(i,k)
            addcondt(i,k) = 0.
            if (fluxin .gt. fluxout) then
              if (netflux .lt. fluxin*tinynum) then
                adnetflux = 0.
              endif
            else
              if ((-netflux) .lt. fluxout*tinynum) then
                adnetflux = 0.
              endif
            endif
            adfluxin = adfluxin+adnetflux
            adfluxout = adfluxout-adnetflux
            adnetflux = 0.
            adcmix(i,k) = adcmix(i,k)+adfluxout*(mu(i,kp1)+md(i,kp1))
            adcond(i,kp1) = adcond(i,kp1)-adfluxout*md(i,kp1)
            adconu(i,k) = adconu(i,k)+adfluxout*mu(i,k)
            adfluxout = 0.
            adcmix(i,km1) = adcmix(i,km1)+adfluxin*(mu(i,k)+md(i,k))
            adcond(i,k) = adcond(i,k)-adfluxin*md(i,k)
            adconu(i,kp1) = adconu(i,kp1)+adfluxin*mu(i,kp1)*fisg(i,k)
            adfluxin = 0.
          end do
        end do
        do k = llpar, 3, -1
           KM1 = MAX( 1,     K-1 )
          do i = il1g, il2g
            if (md(i,k) .lt. (-mbsth)) then
              adcondh = adcond(i,k)
              adcond(i,k) = 0.
              adcmix(i,km1) = adcmix(i,km1)-adcondh*(ed(i,km1)*dp(i,km1)
     $/md(i,k))
              adcond(i,km1) = adcond(i,km1)+adcondh*(md(i,km1)/md(i,k))
            endif
          end do
        end do
        do kk = 1, llpar-1
          if (llpar .gt. kk+1) then
            kkp1 = kk+1
          else
            kkp1 = llpar
          endif
          do i = il1g, il2g
            mupdudp = mu(i,kk)+du(i,kk)*dp(i,kk)
            if (mupdudp .gt. mbsth) then
              adconuh = adconu(i,kk)
              adconu(i,kk) = 0.
              adcmix(i,kk) = adcmix(i,kk)+adconuh*(eu(i,kk)*dp(i,kk)/
     $mupdudp)
              adconu(i,kkp1) = adconu(i,kkp1)+adconuh*(mu(i,kkp1)*
     $fisg(i,kk)/mupdudp)
            endif
          end do
        end do
        k = 2
        km1 = 1
        kk = llpar
        do i = il1g, il2g
          mupdudp = mu(i,kk)+du(i,kk)*dp(i,kk)
          if (md(i,k) .lt. (-mbsth)) then
            adcmix(i,km1) = adcmix(i,km1)-adcond(i,k)*(ed(i,km1)*dp(i,
     $km1)/md(i,k))
            adcond(i,k) = 0.
          endif
          if (mupdudp .gt. mbsth) then
            adcmix(i,kk) = adcmix(i,kk)+adconu(i,kk)*(eu(i,kk)*dp(i,kk)/
     $mupdudp)
            adconu(i,kk) = 0.
          endif
        end do
        do k = llpar, 1, -1
          if (k-1 .gt. 1) then
            km1 = k-1
          else
            km1 = 1
          endif
          do i = il2g, il1g, -1
            if (cmix(i,km1) .gt. cmix(i,k)) then
              minc = cmix(i,k)
              maxc = cmix(i,km1)
            else
              minc = cmix(i,km1)
              maxc = cmix(i,k)
            endif
            if (minc .lt. 0.d0) then
              cdifr = 0.d0
            else
              if (maxc .gt. small) then
                tmp = maxc
              else
                tmp = small
              endif
              if (cmix(i,k) .gt. cmix(i,km1)) then
                cdifr = cmix(i,k)-cmix(i,km1)/tmp
              else
                cdifr = cmix(i,km1)-cmix(i,k)/tmp
              endif
            endif
            addcondt(i,k) = 0.
            adchat(i,k) = adchat(i,k)+adcond(i,k)
            adcond(i,k) = 0.
            adchat(i,k) = adchat(i,k)+adconu(i,k)
            adconu(i,k) = 0.
            if (cdifr .gt. 1.d-6) then
              if (maxc*tinynum .gt. smallest) then
                if (cmix(i,km1) .gt. maxc*tinynum) then
                  cabv = cmix(i,km1)
                else
                  cabv = maxc*tinynum
                endif
                if (cmix(i,k) .gt. maxc*tinynum) then
                  cbel = cmix(i,k)
                else
                  cbel = maxc*tinynum
                endif
              else
                if (cmix(i,km1) .gt. smallest) then
                  cabv = cmix(i,km1)
                else
                  cabv = smallest
                endif
                if (cmix(i,k) .gt. smallest) then
                  cbel = cmix(i,k)
                else
                  cbel = smallest
                endif
              endif
              adcabv = adcabv+adchat(i,k)*(log(cabv/cbel)/(cabv-cbel)+
     $(1./(cabv/cbel)/cbel/(cabv-cbel)-log(cabv/cbel)/((cabv-cbel)*
     $(cabv-cbel)))*cabv)*cbel
              adcbel = adcbel+adchat(i,k)*(log(cabv/cbel)/(cabv-cbel)*
     $cabv+((-(1./(cabv/cbel)*(cabv/(cbel*cbel))/(cabv-cbel)))+log(cabv/
     $cbel)/((cabv-cbel)*(cabv-cbel)))*cabv*cbel)
              adchat(i,k) = 0.
              if (maxc*tinynum .gt. smallest) then
                if (cmix(i,k) .gt. maxc*tinynum) then
                  adcmix(i,k) = adcmix(i,k)+adcbel
                  adcbel = 0.
                else
                  admaxc = admaxc+adcbel*tinynum
                  adcbel = 0.
                endif
                if (cmix(i,km1) .gt. maxc*tinynum) then
                  adcmix(i,km1) = adcmix(i,km1)+adcabv
                  adcabv = 0.
                else
                  admaxc = admaxc+adcabv*tinynum
                  adcabv = 0.
                endif
              else
                if (cmix(i,k) .gt. smallest) then
                  adcmix(i,k) = adcmix(i,k)+adcbel
                  adcbel = 0.
                else
                  adcbel = 0.
                endif
                if (cmix(i,km1) .gt. smallest) then
                  adcmix(i,km1) = adcmix(i,km1)+adcabv
                  adcabv = 0.
                else
                  adcabv = 0.
                endif
              endif
            else
              adcmix(i,k) = adcmix(i,k)+0.5d0*adchat(i,k)
              adcmix(i,km1) = adcmix(i,km1)+0.5d0*adchat(i,k)
              adchat(i,k) = 0.
            endif
            if (minc .lt. 0.d0) then
              adcdifr = 0.
            else
              if (cmix(i,k) .gt. cmix(i,km1)) then
                adcmix(i,k) = adcmix(i,k)+adcdifr
                adcmix(i,km1) = adcmix(i,km1)-adcdifr/tmp
                adtmp = adtmp+adcdifr*(cmix(i,km1)/(tmp*tmp))
                adcdifr = 0.
              else
                adcmix(i,k) = adcmix(i,k)-adcdifr/tmp
                adcmix(i,km1) = adcmix(i,km1)+adcdifr
                adtmp = adtmp+adcdifr*(cmix(i,k)/(tmp*tmp))
                adcdifr = 0.
              endif
              if (maxc .gt. small) then
                admaxc = admaxc+adtmp
                adtmp = 0.
              else
                adtmp = 0.
              endif
            endif
            if (cmix(i,km1) .gt. cmix(i,k)) then
              adcmix(i,km1) = adcmix(i,km1)+admaxc
              admaxc = 0.
              adcmix(i,k) = adcmix(i,k)+adminc
              adminc = 0.
            else
              adcmix(i,k) = adcmix(i,k)+admaxc
              admaxc = 0.
              adcmix(i,km1) = adcmix(i,km1)+adminc
              adminc = 0.
            endif
          end do
        end do
        do k = 1, llpar
          do i = il1g, il2g
            cmix(i,k) = q(ideep(i),k,m)
            if (cmix(i,k) .lt. 4.d0*smallest) then
              adcmix(i,k) = 0.
            endif
            adq(ideep(i),k,m) = adq(ideep(i),k,m)+adcmix(i,k)
            adcmix(i,k) = 0.
          end do
        end do

      ENDDO  !M ; End of tracer loop

      ! Return to calling program
      END SUBROUTINE CONVTRAN_ADJ

!-----------------------------------------------------------------------------

      END MODULE FVDAS_CONVECT_ADJ_MOD
