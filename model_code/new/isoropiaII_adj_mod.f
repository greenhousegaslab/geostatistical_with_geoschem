!------------------------------------------------------------------------------
!         Caltech Department of Chemical Engineering / Seinfeld Group         !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: isoropiaii_adj_mod
!
! !DESCRIPTION: Module ISOROPIAII_ADJ\_MOD contains the routines that provide
!  the interface between ISORROPIA II and GEOS-Chem.
!\\
!\\
!  The actual ANISORROPIA code which performs Na-SO4-NH3-NO3-Cl
!  aerosol thermodynamic equilibrium is in \texttt{isoropiaIIcode_adj.f}.
!\\
!\\
! !INTERFACE:
!
      MODULE ISOROPIAII_ADJ_MOD
!
! !USES:
!
      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: CLEANUP_ISOROPIAII
      PUBLIC  :: DO_ISOROPIAII
      PUBLIC  :: DO_ISOROPIAII_ADJ
      PUBLIC  :: GET_GNO3
      PUBLIC  :: GET_ISRINFO
!
! !PRIVATE MEMBER FUNCTIONS:
!
      PRIVATE :: GET_HNO3
      PRIVATE :: INIT_ISOROPIAII
      PRIVATE :: SAFELOG10
      PRIVATE :: SET_HNO3
!
! !REMARKS:
!  Original Author:
!  *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
!  *** GEORGIA INSTITUTE OF TECHNOLOGY
!  *** WRITTEN BY ATHANASIOS NENES
!  *** UPDATED BY CHRISTOS FOUNTOUKIS
!                                                                             .
!  Original v1.3 isoropia implementation into GEOS-Chem by
!  Becky Alexander and Bob Yantosca (bec, bmy, 4/12/05, 11/2/05)
!                                                                             .
!  For Ca,K,Mg = 0, ISOROPIA II performs exactly like ISOROPIAv1.7
!  Ca, K, Mg, Na from dust is not currently considered
!                                                                             .
!  To implement ISOROPIA II into GEOS-Chem:
!    * cleanup_isoropiaII needs to be called from cleanup.f
!    * DO_ISOROPIA needs to be replaced with DO_ISOROPIAII in chemistry_mod.f
!    * Change ISOROPIA to ISOROPIAII in sulfate_mod.f
!    * add isoropiaII_mod.f, isoropiaIIcode.f, and irspia.inc to Makefile
!                                                                             .
!  ISOROPIA II implementation notes by Havala O.T. Pye:
!  (1) The original isoropia code from T.Nenes is left as unmodified as
!       possible. Original isoropia code can be found in isoropiaIIcode.f
!       and common blocks can be found in isrpia.inc. For future upgrades
!       to isoropia, replace isrpia.inc and isoropiaIIcode.f with the new
!       version of isoropia and modify the call to ISOROPIA in this module.
!       Please let the original author know of any changes made to ISOROPIA.
!  (2) As of Nov 2007, routines using non-zero Ca, K, and Mg do not always
!       conserve mass. Ca, K, and Mg are set to zero.
!                                                                             .
!  NOTE: ISORROPIA is Greek for "equilibrium", in case you were wondering.
!
!  ANISORROPIA implementation in GEOS-Chem adjoint by
!  Shannon Capps (slc, 8/22/2011)
!  (1) As of Aug 2011, only Na-SO4-NH3-NO3-Cl routines have an adjoint.
!  (2) Adjoint calculations require online activity coefficient calculation
!       unlike the default configuration in GEOS-Chem that uses look up tables.
!      Reference: doi:10.5194/acp-12-527-2012
!
! !REVISION HISTORY:
!  06 Jul 2007 - H. O. T. Pye - Initial version
!  29 Jan 2010 - R. Yantosca  - Added ProTeX headers
!  21 Apr 2010 - R. Yantosca  - Bug fix in DO_ISOROPIAII for offline aerosol
!  22 Aug 2011 - S. Capps - ANISORROPIA implementation
!
!        *** VERY IMPORTANT PORTING WARNING (slc.1.2012) ***
!  ANISORROPIA code is optimized for adjoint frameworks and will not
!   perform commensurately with publicly released ISORROPIAII code.
!
!  Please visit http://nenes.eas.gatech.edu/ISORROPIA for current
!   releases of ISORROPIAII for forward modeling.
!
!EOP
!------------------------------------------------------------------------------
!BOC
      ! Array for offline HNO3 (for relaxation of M.M.)
      REAL*8,  ALLOCATABLE :: HNO3_sav(:,:,:)

      ! Array for offline use in sulfate_mod (SEASALT_CHEM)
      REAL*8,  ALLOCATABLE :: GAS_HNO3(:,:,:)

      ! AEROPH: Save information related to aerosol pH (hotp 8/11/09)
      REAL*8,  ALLOCATABLE :: PH_SAV(:,:,:)
      REAL*8,  ALLOCATABLE :: HPLUS_SAV(:,:,:)
      REAL*8,  ALLOCATABLE :: WATER_SAV(:,:,:)
      REAL*8,  ALLOCATABLE :: SULRAT_SAV(:,:,:)
      REAL*8,  ALLOCATABLE :: NARAT_SAV(:,:,:)
      REAL*8,  ALLOCATABLE :: ACIDPUR_SAV(:,:,:)

      CONTAINS
!EOC
!------------------------------------------------------------------------------
!         Caltech Department of Chemical Engineering / Seinfeld Group         !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_isoropiaii
!
! !DESCRIPTION: Subroutine DO\_ISOROPIAII is the interface between the
!  GEOS-Chem model and the aerosol thermodynamical equilibrium routine
!  ISORROPIA II.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE DO_ISOROPIAII
!
! !USES:
!
      USE CHECKPT_MOD,     ONLY : ANISO_IN
      USE DAO_MOD,         ONLY : AIRVOL, RH, T
      USE ERROR_MOD,       ONLY : DEBUG_MSG,       ERROR_STOP
      USE ERROR_MOD,       ONLY : SAFE_DIV
      USE GLOBAL_HNO3_MOD, ONLY : GET_GLOBAL_HNO3
      USE LOGICAL_MOD,     ONLY : LPRT
      USE TIME_MOD,        ONLY : GET_MONTH,       ITS_A_NEW_MONTH
      USE TRACER_MOD
      USE TRACERID_MOD,    ONLY : IDTHNO3, IDTNIT, IDTNH4, IDTNH3
      USE TRACERID_MOD,    ONLY : IDTSALA, IDTSO4
      USE TROPOPAUSE_MOD,  ONLY : ITS_IN_THE_STRAT
!
      USE ADJ_ARRAYS_MOD,  ONLY : IFD, JFD, LFD
      USE LOGICAL_ADJ_MOD, ONLY : LPRINTFD
      USE LOGICAL_ADJ_MOD, ONLY : LADJ


#     include "CMN_SIZE"        ! Size parameters
!
! !REMARKS:
!  Original isoropia v1.3 implmentation: (rjp, bec, bmy, 12/17/01, 8/22/05)
!
! !REVISION HISTORY:
!  24 Aug 2007 - H. O. T. Pye - Initial version, in ISORROPIA II
!  18 Dec 2009 - H. O. T. Pye - Added division checks
!  29 Jan 2010 - R. Yantosca  - Added ProTeX headers
!  21 Apr 2010 - E. Sofen     - Prevent out-of-bounds errors for offline
!                               aerosol simulations where HNO3 is undefined
!  23 Jul 2010 - R. Yantosca  - Bug fix: corrected typo in ND42 diag section
!  22 Aug 2011 - S. Capps     - ANISORROPIA implementation
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      ! Array dimensions
      INTEGER, PARAMETER       :: NOTHERA  =  9
      INTEGER, PARAMETER       :: NCTRLA   =  2
      INTEGER, PARAMETER       :: NCOMPA   =  8
      INTEGER, PARAMETER       :: NIONSA   = 10
      INTEGER, PARAMETER       :: NGASAQA  =  3
      INTEGER, PARAMETER       :: NSLDSA   = 19

      ! Concentration lower limit [mole/m3]
      REAL*8,  PARAMETER       :: CONMIN = 1.0d-30
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE            :: FIRST = .TRUE.
      INTEGER                  :: I,    J,    L,    N
      REAL*8                   :: ANO3, GNO3, RHI,  TEMPI
      REAL*8                   :: TCA,  TMG,  TK,   HNO3_DEN
      REAL*8                   :: TNA,  TCL,  TNH3, TNH4
      REAL*8                   :: TNIT, TNO3, TSO4, VOL
      REAL*8                   :: AERLIQ(NIONSA+NGASAQA+2)
      REAL*8                   :: AERSLD(NSLDSA)
      REAL*8                   :: GAS(NGASAQA)
      REAL*8                   :: OTHER(NOTHERA)
      REAL*8                   :: WI(NCOMPA)
      REAL*8                   :: WT(NCOMPA)
      REAL*8                   :: CNTRL(NCTRLA)
      CHARACTER(LEN=255)       :: X
      CHARACTER(LEN=15)        :: SCASI

      ! Flag and integer indicative of ANISORROPIA internal error system
      LOGICAL                  :: TRUSTISO

      !Temporary variables to check if division is safe
      REAL*8                   :: NUM_SAV, DEN_SAV

      ! AEROPH: Temporary variable for pH (hotp 8/11/09)
      REAL*8                   :: HPLUSTEMP

      ! debug variables
      INTEGER                  :: Itemp, Jtemp, Ltemp
      INTEGER                  :: ISOERRCOUNT,ISOCALLCOUNT
      INTEGER                  :: NERR, NERR22, NERR33, NERR44, NERR100
      INTEGER                  :: NERR101, NERR102, NERR103, NERR104
      INTEGER                  :: NERR50, NERROTHER, COTHER
      INTEGER                  :: CA, CB, CC, CD, CE, CF, CG, CH, CI, CJ
      LOGICAL, SAVE            :: FIRSTCHECK = .TRUE.

      !=================================================================
      ! DO_ISOROPIAII begins here!
      !=================================================================

      ! Location string
      X = 'DO_ISOROPIAII (isoropiaII_mod.f)'
      WRITE(6,*) X

      ! First-time initialization
      IF ( FIRST ) THEN

         ! Make sure certain tracers are defined
         IF ( IDTSO4  == 0 ) CALL ERROR_STOP( 'IDTSO4 is undefined!', X)
         IF ( IDTNH3  == 0 ) CALL ERROR_STOP( 'IDTNH3 is undefined!', X)
         IF ( IDTNH4  == 0 ) CALL ERROR_STOP( 'IDTNH4 is undefined!', X)
         IF ( IDTNIT  == 0 ) CALL ERROR_STOP( 'IDTNIT is undefined!', X)
         IF ( IDTSALA == 0 ) CALL ERROR_STOP( 'IDTSALA is undefined!',X)

         ! Initialize arrays
         CALL INIT_ISOROPIAII
         !WRITE(*,*) 'Successfully finished INIT_ISOROPIAII'

         ! Reset first-time flag
         FIRST = .FALSE.

         ! Reset error counting flag
         ISOERRCOUNT = 0

      ENDIF

      !=================================================================
      ! Check to see if we have to read in monthly mean HNO3
      !=================================================================
      IF ( IDTHNO3 == 0 ) THEN

         IF ( ITS_A_FULLCHEM_SIM() ) THEN

            ! Coupled simulation: stop w/ error since we need HNO3
            CALL ERROR_STOP( 'IDTHNO3 is not defined!', X )

         ELSE IF ( ITS_AN_AEROSOL_SIM() ) THEN

            ! Offline simulation: read monthly mean HNO3
            IF ( ITS_A_NEW_MONTH() ) THEN
               CALL GET_GLOBAL_HNO3( GET_MONTH() )
            ENDIF

            ! Initialize for each timestep (bec, bmy, 4/15/05)
            GAS_HNO3 = 0d0

         ELSE

            ! Otherwise stop w/ error
            CALL ERROR_STOP( 'Invalid simulation type!', X )

         ENDIF
      ENDIF

      ! AEROPH: Initialize arrays all the way up to LLPAR for
      ! aeroph. Arrays go up to LLPAR due to ND42 use (hotp 8/11/09)
      PH_SAV      = 0d0
      HPLUS_SAV   = 0d0
      WATER_SAV   = 0d0
      SULRAT_SAV  = 0d0
      NARAT_SAV   = 0d0
      ACIDPUR_SAV = 0d0

      ! Initialize the error distribution flags
      NERR22 = 0
      NERR33 = 0
      NERR44 = 0
      NERR100 = 0
      NERR101 = 0
      NERR102 = 0
      NERR103 = 0
      NERR104 = 0
      NERROTHER = 0

      ISOCALLCOUNT = 0
      ISOERRCOUNT = 0

      CA = 0
      CB = 0
      CC = 0
      CD = 0
      CE = 0
      CF = 0
      CG = 0
      CH = 0
      CI = 0
      CJ = 0
      COTHER = 0

      IF ( LADJ ) THEN                      ! adj_group
         ANISO_IN(:,:,:,1:14) = 0.d0
      ENDIF

      !WRITE(*,*) 'ANISO_IN: ',ANISO_IN(1,1,1,:)

      !=================================================================
      ! Loop over grid boxes and call ISOROPIA (see comments in the
      ! ISOROPIA routine ISOROPIAIICODE.f which describes
      ! the input/output args)
      !=================================================================

      ! AEROPH: add HPLUSTEMP as private (hotp 8/11/09)
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I,    J,      L,       N,      WI,   WT,  GAS,  TEMPI )
!$OMP+PRIVATE( RHI,  VOL,    TSO4,    TNH3,   TNA,  TCL, ANO3, GNO3  )
!$OMP+PRIVATE( TCA,  TMG,    TK,      CNTRL,  SCASI,     TRUSTISO    )
!$OMP+PRIVATE( TNO3, AERLIQ, AERSLD,  OTHER,  TNH4, TNIT,      NERR  )
!$OMP+PRIVATE( HPLUSTEMP,    NUM_SAV, DEN_SAV, HNO3_DEN              )
!$OMP+SCHEDULE( DYNAMIC )

      DO L = 1, LLTROP
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Skip strat boxes
         IF ( ITS_IN_THE_STRAT( I, J, L ) ) CYCLE

         ! Initialize WI, WT
         DO N = 1, NCOMPA
            WI(N) = 0d0
            WT(N) = 0d0
         ENDDO

         ! Initialize GAS
         DO N = 1, NGASAQA
            GAS(N) = 0d0
         ENDDO

         ! Temperature [K]
         TEMPI    = T(I,J,L)

         ! Relative humidity [unitless]
         RHI      = RH(I,J,L) * 1.d-2

         ! Force RH in the range 0.01 - 0.98
         RHI      = MAX( 0.01d0, RHI )
         RHI      = MIN( 0.98d0, RHI )

         ! Volume of grid box [m3]
         VOL      = AIRVOL(I,J,L)

         !---------------------------------
         ! Compute quantities for ISOROPIA
         !---------------------------------

         ! Total SO4 [mole/m3]
         !    Convert from kg to mole/m3 air
         TSO4     = STT(I,J,L,IDTSO4) * 1.d3 / ( 96.d0 * VOL )

         ! Total NH3 [mole/m3]
         !    Convert from kg to mole/m3 air
         TNH3     = STT(I,J,L,IDTNH4) * 1.d3 / ( 18.d0 * VOL )  +
     &              STT(I,J,L,IDTNH3) * 1.d3 / ( 17.d0 * VOL )

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%% NOTE: The error-trap statement above will halt execution if IDTSALA is
!%%% undefined.  Therefore this IF statement is superfluous.  Comment out
!%%% for clarity.  (hotp, bmy, 2/1/10)
!%%%
!%%%         IF ( IDTSALA > 0 ) THEN

            ! Total Na+ (30.61% by weight of seasalt) [mole/m3]
            TNA      = STT(I,J,L,IDTSALA) * 0.3061d0 * 1.d3 /
     &                                    ( 22.99d0  * VOL  )

            ! Total Cl- (55.04% by weight of seasalt) [mole/m3]
            TCL      = STT(I,J,L,IDTSALA) * 0.5504d0 * 1.d3 /
     &                                    ( 35.45d0  * VOL  )

!==============================================================================
!=== NOTE: As of 11/2007, ISORROPIAII does not conserve mass when Ca,K,Mg are
!=== non-zero. If you would like to consider Ca, K, Mg from seasalt and dust,
!=== isoropiaIIcode.f ISRP4F routines must be debugged.  (hotp, bmy, 2/1/10)
!===
!===            ! Total Ca2+ (1.16% by weight of seasalt) [mole/m3]
!===            TCA      = STT(I,J,L,IDTSALA) * 0.0116d0 * 1.d3 /
!===     &                                 ( 40.08d0  * VOL  )
!===
!===            ! Total K+   (1.1% by weight of seasalt)  [mole/m3]
!===            TK       = STT(I,J,L,IDTSALA) * 0.0110d0 * 1.d3 /
!===     &                                 ( 39.102d0 * VOL  )
!===
!===            ! Total Mg+  (3.69% by weight of seasalt) [mole/m3]
!===            TMG      = STT(I,J,L,IDTSALA) * 0.0369d0 * 1.d3 /
!===     &                                 ( 24.312d0 * VOL  )

            ! Set Ca, K, Mg to zero for time being (hotp, bmy, 2/1/10)
            TCA      = 0d0
            TK       = 0d0
            TMG      = 0d0
!==============================================================================
!%%%         ELSE
!%%%
!%%%            ! no seasalt, set to zero
!%%%            TNA = 0.d0
!%%%            TCL = 0.d0
!%%%            TCA = 0.d0
!%%%            TK  = 0.d0
!%%%            TMG = 0.d0
!%%%
!%%%         ENDIF
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         ! Compute gas-phase NO3
         IF ( IDTHNO3 > 0 ) THEN

            !---------------------
            ! COUPLED SIMULATION
            !---------------------

            ! Compute gas-phase HNO3 [mole/m3] from HNO3 tracer
            GNO3  = STT(I,J,L,IDTHNO3)
            GNO3  = MAX( GNO3 * 1.d3 / ( 63.d0 * VOL ), CONMIN )

            ! Aerosol-phase NO3 [mole/m3]
            ANO3     = STT(I,J,L,IDTNIT) * 1.d3 / ( 62.d0 * VOL )

            ! Total NO3 [mole/m3]
            TNO3    = GNO3 + ANO3

         ELSE

            !---------------------
            ! OFFLINE SIMULATION
            !---------------------

            ! Convert total inorganic NO3 from [ug/m3] to [mole/m3].
            ! GET_HNO3, lets HNO3 conc's evolve, but relaxes to
            ! monthly mean values every 3h.
            TNO3  = GET_HNO3( I,J,L ) * 1.d-6 / 63.d0

         ENDIF

         !---------------------------------
         ! Call ISOROPIAII
         !---------------------------------

         ! set type of ISOROPIA call
         ! Forward problem, do not change this value
         ! 0d0 represents forward problem
         CNTRL(1) = 0.0d0

         ! Metastable for now
         ! 1d0 represents metastable problem
         CNTRL(2) = 1.0d0

         ! Insert concentrations [mole/m3] into WI & prevent underflow
         WI(1)    = MAX( TNA,  CONMIN )
         WI(2)    = MAX( TSO4, CONMIN )
         WI(3)    = MAX( TNH3, CONMIN )
         WI(4)    = MAX( TNO3, CONMIN )
         WI(5)    = MAX( TCL,  CONMIN )
         WI(6)    = MAX( TCA,  CONMIN )
         WI(7)    = MAX( TK,   CONMIN )
         WI(8)    = MAX( TMG,  CONMIN )

         IF ( LPRINTFD
     &      .and. J == JFD .AND. L == LFD .AND. I == IFD) THEN
            WRITE(6,*) 'Forward, CELL(',I,',',J,',',L,')',WI(1:5)
            WRITE(6,*) 'Temp',TEMPI, '   RHI',RHI
         ENDIF

         ! Perform aerosol thermodynamic equilibrium
         ! ISOROPIAII can be found in isoropiaIIcode_adj.f
         ! inputs are WI, RHI, TEMPI, CNTRL

         ! adj_group: call special version for adjoint (slc.09.2011)
         IF ( .not. LADJ ) THEN
            CALL ISOROPIAII (WI,    RHI,  TEMPI,  CNTRL,
     &                       WT,    GAS,  AERLIQ, AERSLD,
     &                       SCASI, OTHER, TRUSTISO,NERR)

         ELSE
            ! Checkpoint ANISORROPIA input
            ANISO_IN(I,J,L,1:8) = WI(:)
            ANISO_IN(I,J,L,9)   = RHI
            ANISO_IN(I,J,L,10)  = TEMPI

            CALL ISOROPIAII (WI,    RHI,  TEMPI,  CNTRL,
     &                       WT,    GAS,  AERLIQ, AERSLD,
     &                       SCASI, OTHER, TRUSTISO,NERR)

            ! Debug ANISO checkpoint
            IF ( LPRINTFD
     &         .and. J == JFD .AND. L == LFD .AND. I == IFD) THEN
               WRITE(6,*) ' After ISOROPIAII ', I, J, L
               ! (slc.10.2011) debug - output of ISORROPIA
               WRITE(6,*) 'GAS ',GAS(:)
               WRITE(6,*) 'AERLIQ: 3,5,6,7 ',AERLIQ(3),
     &         AERLIQ(5),AERLIQ(6),AERLIQ(7)
               WRITE(6,*) 'TRUSTISO ',TRUSTISO
            ENDIF

         ENDIF        ! Checkpointing for adjoint

         !---------------------------------
         ! Save back into tracer array
         !---------------------------------

         ! Convert ISOROPIA output from [mole/m3] to [kg]
         TSO4 = MAX( 96.d-3 * VOL *   WT(2),            CONMIN )
         TNH3 = MAX( 17.d-3 * VOL *   GAS(1),           CONMIN )

         ! (slc.9.2011) - for adjoint to work without WT_ADJ
         ! TNH4 = MAX( 18.d-3 * VOL * ( WT(3) - GAS(1) ), CONMIN )
         ! TNIT = MAX( 62.d-3 * VOL * ( WT(4) - GAS(2) ), CONMIN )
         TNH4 = MAX( 18.d-3 * VOL *   AERLIQ(3),        CONMIN )
         TNIT = MAX( 62.d-3 * VOL *   AERLIQ(7),        CONMIN )

         !------------------------------------
         ! Check as to whether error occurred.
         !------------------------------------

         IF ( TRUSTISO )  THEN

            ! Save tracers back into STT array [kg]
            ! no longer save TSO4 back into STT. SO4 is all aerosol phase
            ! (hotp 11/7/07)
            ! STT(I,J,L,IDTSO4) = TSO4
            STT(I,J,L,IDTNH3) = TNH3
            STT(I,J,L,IDTNH4) = TNH4
            STT(I,J,L,IDTNIT) = TNIT

         ! slc.debug
            IF ( LADJ ) THEN                      ! adj_group
               IF ( 17.d-3 * VOL * GAS(1) < CONMIN ) THEN

               !WRITE(*,*) 'CONMIN > NH3', GAS(1)
               !WRITE(*,*) 'CELL:(',I,',',J,',',L,')'
                  ANISO_IN(I,J,L,11)   = 0.d0

               ELSE

                  ANISO_IN(I,J,L,11)   = 1.d0

               ENDIF

               IF ( 18.d-3 * VOL * AERLIQ(3) < CONMIN ) THEN

               !WRITE(*,*) 'CONMIN > NH4', AERLIQ(3)
               !WRITE(*,*) 'CELL:(',I,',',J,',',L,')'
                  ANISO_IN(I,J,L,12)   = 0.d0

               ELSE

                  ANISO_IN(I,J,L,12)   = 1.d0

               ENDIF

               IF ( 62.d-3 * VOL * AERLIQ(7) < CONMIN ) THEN

               !WRITE(*,*) 'CONMIN > NIT', AERLIQ(7)
               !WRITE(*,*) 'CELL:(',I,',',J,',',L,')'

                  ANISO_IN(I,J,L,13)   = 0.d0

               ELSE

                  ANISO_IN(I,J,L,13)   = 1.d0

               ENDIF

               IF ( 96.d-3 * VOL * WT(2) < CONMIN ) THEN

               !WRITE(*,*) 'CONMIN > SUL', WT(2)
               !WRITE(*,*) 'CELL:(',I,',',J,',',L,')'

                  ANISO_IN(I,J,L,15)   = 0.d0

               ELSE

                  ANISO_IN(I,J,L,15)   = 1.d0

               ENDIF

            ENDIF

         ELSE

            ! Echo location of NAN (probably leave this commented out
            ! unless you are getting lots of ADJ_NAN warnings
            !WRITE(6,*) 'Can't trust ANISO at I,J,L,N = ',I,J,L,N

            IF ( LADJ ) THEN                      ! adj_group

               ANISO_IN(I,J,L,11:15)   = 0.d0

            ENDIF

            !WRITE(*,*) 'ANISO_IN when TRUSTISO = .F.',
            !&                  ANISO_IN(I,J,L,11:14)

!!$OMP CRITICAL
            ! Show TRUSTISO flag so that a warning is echoed to screen
            TRUSTISO = .FALSE.
!!$OMP END CRITICAL

            ! Count number of errors and total calls
            ISOERRCOUNT = ISOERRCOUNT + 1


            SELECT CASE (NERR)
               CASE (22)
                  NERR22 = NERR22 + 1
               CASE (33)
                  NERR33 = NERR33 + 1
               CASE (50)
                  NERR50 = NERR50 + 1
               CASE (100)
                  NERR100 = NERR100 + 1
               CASE (101)
                  NERR101 = NERR101 + 1
               CASE (102)
                  NERR102 = NERR102 + 1
               CASE (103)
                  NERR103 = NERR103 + 1
               CASE (104)
                  NERR104 = NERR104 + 1
                CASE DEFAULT
                   NERROTHER = NERROTHER + 1
            END SELECT


            ! Do not replace original value
            !STT(I,J,L,IDTNH3) = STT(I,J,L,IDTNH3)
            !STT(I,J,L,IDTNH4) = STT(I,J,L,IDTNH4)

         ENDIF

         ! slc.debug

         ISOCALLCOUNT = ISOCALLCOUNT + 1

         SELECT CASE (SCASI)
            CASE("A2")
               CA = CA + 1
            CASE("B4")
              CB = CB + 1
            CASE("C2")
              CC = CC + 1
            CASE("D3")
              CD = CD + 1
            CASE("E4")
              CE = CE + 1
            CASE("F2")
              CF = CF + 1
            CASE("G5")
              CG = CG + 1
            CASE("H6")
              CH = CH + 1
            CASE("I6")
              CI = CI + 1
            CASE("J3")
              CJ = CJ + 1
            CASE DEFAULT
              COTHER = COTHER + 1
         END SELECT

         ! Special handling for HNO3 [kg]
         IF ( IDTHNO3 > 0 ) THEN

            !---------------------
            ! COUPLED SIMULATION
            !---------------------

            !------------------------------------
            ! Check as to whether error occurred.
            !------------------------------------

            IF ( TRUSTISO )  THEN

               ! HNO3 [mole/m3] is in GAS(2); convert & store in STT [kg]
               STT(I,J,L,IDTHNO3) = MAX( 63.d-3 * VOL * GAS(2), CONMIN )

            ! slc.debug
            IF ( LADJ ) THEN                      ! adj_group

               IF  ( 63.d-3 * VOL * GAS(2) < CONMIN ) THEN

               !WRITE(*,*) 'CONMIN > HNO3', STT(I,J,L,IDTHNO3)
                  ANISO_IN(I,J,L,14)   = 0.d0

               ELSE

                  ANISO_IN(I,J,L,14)   = 1.d0

               ENDIF

            ENDIF

            ! Save for use in DEN_SAV expression below (sofen, 4/21/10)
               HNO3_DEN           = STT(I,J,L,IDTHNO3)

            ENDIF

         ELSE

            !---------------------------------
            ! Check for trustworthiness.
            !---------------------------------

            IF ( TRUSTISO )  THEN

            !---------------------
            ! OFFLINE SIMULATION:
            !---------------------

            ! Convert total inorganic nitrate from [mole/m3] to [ug/m3]
            ! and save for next time
            ! WT(4) is in [mole/m3] -- unit conv is necessary!
            CALL SET_HNO3( I, J, L, 63.d6 * WT(4) )

            ! Save for use in sulfate_mod (SEASALT_CHEM) for offline
            ! aerosol simulations (bec, 4/15/05)
            GAS_HNO3(I,J,L) = GAS(2)

            ! Save for use in DEN_SAV expression below (sofen, 4/21/10)
            HNO3_DEN        = GAS(2) * VOL * 63d-3

                !---------------------------------
                ! Check for trustworthiness.
                !---------------------------------

                !IF ( .NOT. TRUSTISO )  THEN

                !   STT(I,J,L,IDTHNO3) = STT(I,J,L,IDTHNO3)
                !   STT(I,J,L,IDTNIT)  = STT(I,J,L,IDTNIT)

            ENDIF

         ENDIF

         !---------------------------------
         ! Check for trustworthiness.
         !---------------------------------

!         IF ( TRUSTISO )  THEN
!
!         !-------------------------
!         ! ND42 diagnostic arrays
!         !-------------------------
!
!         ! AEROPH: get pH related info to SAV arrays (hotp 8/11/09)
!         ! HPLUSTEMP is H+ in mol/L water, AERLIQ1 is H, AERLIQ8 is H2O
!         ! in mol/m3 air --> convert to mol/L water
!         IF ( AERLIQ(8) < 1d-32 ) THEN
!            ! Aerosol is dry so HPLUSTEMP and PH_SAV are undefined
!            ! We force HPLUSTEMP to 1d20 and PH_SAV to -999d0.
!            ! (hotp, ccc, 12/18/09)
!            HPLUSTEMP       = 1d20
!            !-------------------------------------------------------------
!            ! Prior to 7/23/10:
!            ! Bug fix: this should be PH_SAV(I,J,L) (sofen, bmy, 7/12/10)
!            !PH_SAV          = -999d0
!            !-------------------------------------------------------------
!            PH_SAV(I,J,L)   = -999d0
!         ELSE
!            HPLUSTEMP       = AERLIQ(1) / AERLIQ(8) * 1d3/18d0
!
!            ! Use SAFELOG10 to prevent NAN
!            PH_SAV(I,J,L)   = -1d0 * SAFELOG10( HPLUSTEMP )
!         ENDIF
!
!         ! Additional Info
!         HPLUS_SAV(I,J,L)   = AERLIQ(1)
!         WATER_SAV(I,J,L)   = AERLIQ(8)
!         SULRAT_SAV(I,J,L)  = OTHER(2)
!         NARAT_SAV(I,J,L)   = OTHER(4)
!
!         NUM_SAV            = ( STT(I,J,L,IDTNH3) /17d0         +
!     &                          STT(I,J,L,IDTNH4) /18d0         +
!     &                          STT(I,J,L,IDTSALA)*0.3061d0/23.0d0 )
!
!         DEN_SAV            = ( STT(I,J,L,IDTSO4)  / 96d0   * 2d0     +
!     &                          STT(I,J,L,IDTNIT)  / 62d0             +
!     &                          HNO3_DEN           / 63d0             +
!     &                          STT(I,J,L,IDTSALA) * 0.55d0 / 35.45d0 )
!
!         ! Value if DEN_SAV and NUM_SAV too small.
!         ACIDPUR_SAV(I,J,L) = SAFE_DIV(NUM_SAV, DEN_SAV,
!     &                                 0d0,
!     &                                 999d0)
!
!         ENDIF

      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !WRITE(*,*) 'Finished with OMP loop in ISOII'  !slc.debug

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### ISOROPIAII: a AERO_THERMO' )

!      WRITE(6,*)  'ISO calls: ',ISOCALLCOUNT
!      WRITE(6,*)  'ISO error occurrences: ',ISOERRCOUNT
!      WRITE(6,*)  'Specific error codes: '
!      WRITE(6,*)  'Error 22: ',NERR22
!      WRITE(6,*)  'Error 33: ',NERR33
!      WRITE(6,*)  'Error 44: ',NERR44
!      WRITE(6,*)  'Error 100: ',NERR100
!      WRITE(6,*)  'Error 101: ',NERR101
!      WRITE(6,*)  'Error 102: ',NERR102
!      WRITE(6,*)  'Error 103: ',NERR103
!      WRITE(6,*)  'Error 104: ',NERR104
!      WRITE(6,*)  'Error Other: ', NERROTHER
!
!      WRITE(6,*) '____________ Case Distribution ____________'
!      WRITE(6,*) 'A: ',CA
!      WRITE(6,*) 'B: ',CB
!      WRITE(6,*) 'C: ',CC
!      WRITE(6,*) 'D: ',CD
!      WRITE(6,*) 'E: ',CE
!      WRITE(6,*) 'F: ',CF
!      WRITE(6,*) 'G: ',CG
!      WRITE(6,*) 'H: ',CH
!      WRITE(6,*) 'I: ',CI
!      WRITE(6,*) 'J: ',CJ
!      WRITE(6,*) 'Other: ', COTHER
!
      ! Return to calling program
      END SUBROUTINE DO_ISOROPIAII
!EOC

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_isoropiaii_adj
!
! !DESCRIPTION: Subroutine DO\_ISOROPIAII_ADJ is the interface between the
!  GEOS-Chem model and the adjoint of the aerosol thermodynamical
!  equilibrium routine ISORROPIA II.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE DO_ISOROPIAII_ADJ
!
! !USES:
!
      USE ADJ_ARRAYS_MOD,  ONLY : IFD, JFD, LFD
      USE ADJ_ARRAYS_MOD,  ONLY : STT_ADJ
      USE CHECKPT_MOD,     ONLY : ANISO_IN
      USE DAO_MOD,         ONLY : AIRVOL, RH, T
      USE ERROR_MOD,       ONLY : DEBUG_MSG,       ERROR_STOP
      USE ERROR_MOD,       ONLY : SAFE_DIV,        IT_IS_NAN
      USE GLOBAL_HNO3_MOD, ONLY : GET_GLOBAL_HNO3
      USE LOGICAL_ADJ_MOD, ONLY : LPRINTFD
      USE LOGICAL_MOD,     ONLY : LPRT
      USE TIME_MOD,        ONLY : GET_MONTH,       ITS_A_NEW_MONTH
      USE TRACER_MOD
      USE TRACERID_MOD,    ONLY : IDTHNO3, IDTNIT, IDTNH4, IDTNH3
      USE TRACERID_MOD,    ONLY : IDTSALA, IDTSO4
      USE TROPOPAUSE_MOD,  ONLY : ITS_IN_THE_STRAT

#     include "CMN_SIZE"        ! Size parameters
!
! !REMARKS:
!  Original isoropia v1.3 implementation: (rjp, bec, bmy, 12/17/01, 8/22/05)
!
! !REVISION HISTORY:
!  30 Aug 2011 - S. Capps     - Interface ANISORROPIA with adjoint
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      ! Array dimensions
      INTEGER, PARAMETER       :: NOTHERA  =  9
      INTEGER, PARAMETER       :: NCTRLA   =  2
      INTEGER, PARAMETER       :: NCOMPA   =  8
      INTEGER, PARAMETER       :: NIONSA   = 10
      INTEGER, PARAMETER       :: NGASAQA  =  3
      INTEGER, PARAMETER       :: NSLDSA   = 19

      ! Concentration lower limit [mole/m3]
      REAL*8,  PARAMETER       :: CONMIN = 1.0d-30

      ! Adjoint parameters
      INTEGER, PARAMETER       ::  MAX_ALLOWED_NAN      = 10
      INTEGER, PARAMETER       ::  MAX_ALLOWED_EXPLD    = 10
      REAL*8,  PARAMETER       ::  MAX_ALLOWED_INCREASE = 10.0D10

!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE            :: FIRSTADJ = .TRUE.
      INTEGER                  :: I,    J,    L,    N
      REAL*8                   :: ANO3, GNO3, RHI,  TEMPI
      REAL*8                   :: TCA,  TMG,  TK,   HNO3_DEN
      REAL*8                   :: TNA,  TCL,  TNH3, TNH4
      REAL*8                   :: TNIT, TNO3, TSO4, VOL
      REAL*8                   :: AERLIQ(NIONSA+NGASAQA+2)
      REAL*8                   :: AERSLD(NSLDSA)
      REAL*8                   :: GAS(NGASAQA)
      REAL*8                   :: OTHER(NOTHERA)
      REAL*8                   :: WI(NCOMPA)
      REAL*8                   :: WT(NCOMPA)
      REAL*8                   :: CNTRL(NCTRLA)
      CHARACTER(LEN=255)       :: X
      CHARACTER(LEN=15)        :: SCASI
      LOGICAL                  :: TRUSTISO

      ! Adjoint variables
      LOGICAL                  :: ADJ_NAN = .FALSE.
      INTEGER                  :: ADJ_NAN_COUNT, ADJ_EXPLD_COUNT
      REAL*8                   :: WT_ADJ(NCOMPA)
      REAL*8                   :: WI_ADJ(NCOMPA)
      REAL*8                   :: GAS_ADJ(NGASAQA)
      REAL*8                   :: AERLIQ_ADJ(NIONSA+NGASAQA+2)
      REAL*8                   :: TNH3_ADJ, TNH4_ADJ, TNO3_ADJ
      REAL*8                   :: TSO4_ADJ, TNIT_ADJ, HNO3_ADJ
      REAL*8                   :: TCA_ADJ, TMG_ADJ, TK_ADJ
      REAL*8                   :: TNA_ADJ, TCL_ADJ
      REAL*8                   :: ANO3_ADJ, GNO3_ADJ
      REAL*8                   :: MAX_ADJ_TMP        ! Temp max value used for error checking

      !Temporary variables to check if division is safe
      REAL*8                   :: NUM_SAV, DEN_SAV

      ! AEROPH: Temporary variable for pH (hotp 8/11/09)
      REAL*8                   :: HPLUSTEMP

      ! debug variables
      INTEGER                  :: Itemp, Jtemp, Ltemp
      INTEGER                  :: ANISOERRCOUNT, ANISOCALLCOUNT
      INTEGER                  :: NERR, NERR22, NERR33, NERR44, NERR100
      INTEGER                  :: NERR101, NERR102, NERR103, NERR104
      INTEGER                  :: NERR50, NERROTHER, COTHER
      INTEGER                  :: CA, CB, CC, CD, CE, CF, CG, CH, CI, CJ
      LOGICAL, SAVE            :: FIRSTCHECK = .TRUE.

      !=================================================================
      ! DO_ISOROPIAII_ADJ begins here!
      !=================================================================

      WRITE(6,*) 'Inside DO_ISOROPIAII_ADJ'
      ! Location string
      X = 'DO_ISOROPIAII_ADJ (isoropiaII_adj_mod.f)'

      ! First-time initialization
      IF ( FIRSTADJ ) THEN

         ! Make sure certain tracers are defined
         IF ( IDTSO4  == 0 ) CALL ERROR_STOP( 'IDTSO4 is undefined!', X)
         IF ( IDTNH3  == 0 ) CALL ERROR_STOP( 'IDTNH3 is undefined!', X)
         IF ( IDTNH4  == 0 ) CALL ERROR_STOP( 'IDTNH4 is undefined!', X)
         IF ( IDTNIT  == 0 ) CALL ERROR_STOP( 'IDTNIT is undefined!', X)
         IF ( IDTSALA == 0 ) CALL ERROR_STOP( 'IDTSALA is undefined!',X)

         !  debug - slc.1.2012
         !! Initialize ADJ_NAN_COUNT
         !ADJ_NAN_COUNT   = 0
         !ADJ_EXPLD_COUNT = 0

         ! Reset first-time flag
         FIRSTADJ = .FALSE.

         ! Initialize error count flag
         ANISOERRCOUNT = 0
      ENDIF

      ! Save maximum adjoint for error checking later
      MAX_ADJ_TMP = MAXVAL( ABS(STT_ADJ) )

      !  debug - slc.1.2012
      !WRITE(*,*) 'Successfully initialized'

      !=================================================================
      ! Check to see if we have to read in monthly mean HNO3
      !=================================================================
      IF ( IDTHNO3 == 0 ) THEN

         IF ( ITS_A_FULLCHEM_SIM() ) THEN

            ! Coupled simulation: stop w/ error since we need HNO3
            CALL ERROR_STOP( 'IDTHNO3 is not defined!', X )

         ELSE IF ( ITS_AN_AEROSOL_SIM() ) THEN

            ! Offline simulation: read monthly mean HNO3
            IF ( ITS_A_NEW_MONTH() ) THEN
               CALL GET_GLOBAL_HNO3( GET_MONTH() )
            ENDIF

            ! Initialize for each timestep (bec, bmy, 4/15/05)
            GAS_HNO3 = 0d0

         ELSE

            ! Otherwise stop w/ error
            CALL ERROR_STOP( 'Invalid simulation type!', X )

         ENDIF
      ENDIF

      !  debug - slc.1.2012
      !WRITE(*,*) 'Successfully checked HNO3'

      ! AEROPH: Initialize arrays all the way up to LLPAR for
      ! aeroph. Arrays go up to LLPAR due to ND42 use (hotp 8/11/09)
      PH_SAV      = 0d0
      HPLUS_SAV   = 0d0
      WATER_SAV   = 0d0
      SULRAT_SAV  = 0d0
      NARAT_SAV   = 0d0
      ACIDPUR_SAV = 0d0

      ! Initialize the error distribution flags
      NERR22 = 0
      NERR33 = 0
      NERR44 = 0
      NERR100 = 0
      NERR101 = 0
      NERR102 = 0
      NERR103 = 0
      NERR104 = 0
      NERROTHER = 0

      ANISOCALLCOUNT = 0
      ANISOERRCOUNT = 0

      CA = 0
      CB = 0
      CC = 0
      CD = 0
      CE = 0
      CF = 0
      CG = 0
      CH = 0
      CI = 0
      CJ = 0
      COTHER = 0

      !=================================================================
      ! Loop over grid boxes and call ISOROPIA (see comments in the
      ! ISOROPIA routine ISOROPIAIICODE.f which describes
      ! the input/output args)
      !=================================================================

      ! AEROPH: add HPLUSTEMP as private (hotp 8/11/09)
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I,    J,      L,       N,      WI,   WT,  GAS,  TEMPI )
!$OMP+PRIVATE( RHI,  VOL,    TSO4,    TNH3,   TNA,  TCL, ANO3, GNO3  )
!$OMP+PRIVATE( TCA,  TMG,    TK,      CNTRL,  SCASI,     TRUSTISO    )
!$OMP+PRIVATE( TNO3, AERLIQ, AERSLD,  OTHER,  TNH4, TNIT,      NERR  )
!$OMP+PRIVATE( HPLUSTEMP,    NUM_SAV, DEN_SAV, HNO3_DEN              )
!$OMP+PRIVATE( WI_ADJ, WT_ADJ, GAS_ADJ, AERLIQ_ADJ, TSO4_ADJ         )
!$OMP+PRIVATE( TMG_ADJ, TK_ADJ, TCA_ADJ, TCL_ADJ, TNO3_ADJ, TNH3_ADJ )
!$OMP+PRIVATE( TNA_ADJ, GNO3_ADJ, ANO3_ADJ                           )

!$OMP+SCHEDULE( DYNAMIC )
      DO L = 1, LLTROP
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Skip strat boxes
         IF ( ITS_IN_THE_STRAT( I, J, L ) ) CYCLE


         ! BEGIN RECALCULATION OF FORWARD VALUES -->
         ! Initialize WI, WT
         WI(:) = 0.d0
         WT(:) = 0.d0

         ! Initialize adjoint variables WI_ADJ, WT_ADJ, GAS_ADJ, AERLIQ_ADJ
         WI_ADJ(:)     = 0d0
         WT_ADJ(:)     = 0d0
         GAS_ADJ(:)    = 0d0
         AERLIQ_ADJ(:) = 0d0

         ! Initialize GAS
         GAS(:) = 0.d0


         ! Volume of grid box [m3]
         VOL      = AIRVOL(I,J,L)


        ! ! Compute gas-phase NO3
        ! IF ( IDTHNO3 > 0 ) THEN
        !
        !    !---------------------
        !    ! COUPLED SIMULATION
        !    !---------------------

        !    ! Compute gas-phase HNO3 [mole/m3] from HNO3 tracer
        !    GNO3  = STT(I,J,L,IDTHNO3)
        !    GNO3  = MAX( GNO3 * 1.d3 / ( 63.d0 * VOL ), CONMIN )

        !    ! Aerosol-phase NO3 [mole/m3]
        !    ANO3     = STT(I,J,L,IDTNIT) * 1.d3 / ( 62.d0 * VOL )

        !    ! Total NO3 [mole/m3]
        !    TNO3    = GNO3 + ANO3

        ! ELSE

        !    !---------------------
        !    ! OFFLINE SIMULATION - no adjoint for this type of run
        !    !---------------------

        !    ! Convert total inorganic NO3 from [ug/m3] to [mole/m3].
        !    ! GET_HNO3, lets HNO3 conc's evolve, but relaxes to
        !    ! monthly mean values every 3h.
        !    TNO3  = GET_HNO3( I,J,L ) * 1.d-6 / 63.d0

        ! ENDIF

         !---------------------------------
         ! Call ANISORROPIA
         !---------------------------------

         ! set type of ANISORROPIA call
         ! Forward problem, do not change this value
         ! 0d0 represents forward problem
         CNTRL(1) = 0.0d0

         ! Metastable for now
         ! 1d0 represents metastable problem
         CNTRL(2) = 1.0d0

         ! From checkpointed files, gather input values (slc.09.27.2011)
         ! Load IN from ANISO_IN

         WI(:)  = ANISO_IN(I,J,L,1:8)

         ! Load parameters from ANISO_IN
         RHI    = ANISO_IN(I,J,L,9)
         TEMPI  = ANISO_IN(I,J,L,10)

         !WRITE(*,*) 'ISO_ADJ, ANISO_IN: ',ANISO_IN(I,J,L,:)
         !WRITE(*,*) 'STT_ADJ(IDTNIT): ', STT_ADJ(I,J,L,IDTNIT)
         !WRITE(*,*) 'STT_ADJ(IDTHNO3): ', STT_ADJ(I,J,L,IDTHNO3)
         !WRITE(*,*) 'STT_ADJ(IDTSO4): ', STT_ADJ(I,J,L,IDTSO4)
         !WRITE(*,*) 'STT_ADJ(IDTNH4): ', STT_ADJ(I,J,L,IDTNH4)
         !WRITE(*,*) 'STT_ADJ(IDTNH3): ', STT_ADJ(I,J,L,IDTNH3)

         !<--- END LOADING OF FORWARD VALUES

         !---> BEGIN ADJOINT CALCULATION


         ! adj code
         IF ( IDTHNO3 > 0 ) THEN
         !   IF ( TRUSTISO )  THEN       ! not defined !

               ! fwd code:
               ! STT(I,J,L,IDTHNO3) = MAX( 63.d-3 * VOL * GAS(2), CONMIN )
               ! adj code:
             IF ( ANISO_IN(I,J,L,14) .GT. 0.d0 ) THEN
                GAS_ADJ(2) = STT_ADJ(I,J,L,IDTHNO3) * 63.d-3 * VOL
             ELSE
                GAS_ADJ(2) = 0.d0
             ENDIF
         !   ENDIF
         ELSE
            CALL ERROR_STOP('adj not supported for offline', X )
         ENDIF

         !IF ( TRUSTISO )  THEN          ! not defined !

             ! fwd code:
             !STT(I,J,L,IDTSO4) = TSO4 - not in forward, but adding for
             !                        adjoint forcing only - slc.4.2013
             !STT(I,J,L,IDTNH3) = TNH3
             !STT(I,J,L,IDTNH4) = TNH4
             !STT(I,J,L,IDTNIT) = TNIT
             ! adj code:

             IF ( ANISO_IN(I,J,L,15) .GT. 0.d0 ) THEN
                TSO4_ADJ = STT_ADJ(I,J,L,IDTSO4)
             ELSE
                TSO4_ADJ = 0.d0
             ENDIF

             IF ( ANISO_IN(I,J,L,13) .GT. 0.d0 ) THEN
                TNIT_ADJ = STT_ADJ(I,J,L,IDTNIT)
             ELSE
                TNIT_ADJ = 0.d0
             ENDIF

             IF ( ANISO_IN(I,J,L,11) .GT. 0.d0 ) THEN
                TNH3_ADJ = STT_ADJ(I,J,L,IDTNH3)
             ELSE
                TNH3_ADJ = 0.d0
             ENDIF

             IF ( ANISO_IN(I,J,L,12) .GT. 0.d0 ) THEN
                TNH4_ADJ = STT_ADJ(I,J,L,IDTNH4)
             ELSE
                TNH4_ADJ = 0.d0
             ENDIF


            ! Debug ANISO checkpoint
            !IF ( LPRINTFD
       !&         .and. J == JFD .AND. L == LFD .AND. I == IFD) THEN
      !       IF ( ( ABS(STT_ADJ(I,J,L,IDTNH4)) .GT. 1d10) .OR.
      !&           ( ABS(STT_ADJ(I,J,L,IDTNH3)) .GT. 1d10) .OR.
      !&           ( ABS(STT_ADJ(I,J,L,IDTHNO3)) .GT. 1d10) .OR.
      !&           ( ABS(STT_ADJ(I,J,L,IDTNIT)) .GT. 1d10) .OR.
      !&           ( ABS(STT_ADJ(I,J,L,IDTSO4)) .GT. 1d10) ) THEN

      !         print*, ' Before ISOROPIAII_ADJ ', I, J, L
      !         print*, ' STT_ADJ(NIT) ', STT_ADJ(I,J,L,IDTNIT)
      !         print*, ' STT_ADJ(HNO3) ', STT_ADJ(I,J,L,IDTHNO3)
      !         print*, ' STT_ADJ(NH4) ', STT_ADJ(I,J,L,IDTNH4)
      !         print*, ' STT_ADJ(NH3) ', STT_ADJ(I,J,L,IDTNH3)
      !         print*, ' STT_ADJ(SO4) ', STT_ADJ(I,J,L,IDTSO4)
      !      ENDIF


         !ENDIF


         ! fwd code:
         !TSO4 = MAX( 96.d-3 * VOL *   WT(2),            CONMIN )
         !TNH3 = MAX( 17.d-3 * VOL *   GAS(1),           CONMIN )
       ! Changing for use of the adjoint without WT_ADJ
       !  !TNH4 = MAX( 18.d-3 * VOL * ( WT(3) - GAS(1) ), CONMIN )
       !  !TNIT = MAX( 62.d-3 * VOL * ( WT(4) - GAS(2) ), CONMIN )
         !TNH4 = MAX( 18.d-3 * VOL *   AERLIQ(3),        CONMIN )
         !TNIT = MAX( 62.d-3 * VOL *   AERLIQ(7),        CONMIN )
         ! adj code (note that we don't overwrite GAS_ADJ(2),
         ! which has already been assigned a value:
         AERLIQ_ADJ(5)  =   96.d-3 * VOL * TSO4_ADJ  ! SO4
         AERLIQ_ADJ(6)  =   97.d-3 * VOL * TSO4_ADJ  ! HSO4
         AERLIQ_ADJ(7)  =   62.d-3 * VOL * TNIT_ADJ
         AERLIQ_ADJ(3)  =   18.d-3 * VOL * TNH4_ADJ
         GAS_ADJ(1)     =   17.d-3 * VOL * TNH3_ADJ

         ! Changes implemented above (slc.4.2013)
         !!! Always zero because the TSO4_ADJ = nothing
         !!! WT_ADJ(2)  =   96.d-3 * VOL * TSO4_ADJ

         !IF ( LPRINTFD
         !&         .and. J == JFD .AND. L == LFD .AND. I == IFD) THEN
            !WRITE(*,*) 'Adjoint, CELL(',I,',',J,',',L,')',WI(1:5)
            !WRITE(*,*) 'Temp ',TEMPI ,' RH ', RHI
            !WRITE(*,*) '------------------------'
            !WRITE(*,*) 'adjoint forcing vectors '
            !WRITE(*,*) 'NH4 force', AERLIQ_ADJ(3)
            !WRITE(*,*) 'NIT force', AERLIQ_ADJ(7)
            !WRITE(*,*) 'NH3 force', GAS_ADJ(1)
            !WRITE(*,*) 'HNO3 force', GAS_ADJ(2)
         !ENDIF

         ! Perform aerosol thermodynamic equilibrium
         ! ISOROPIAII_ADJ can be found in ISOROPIAIICODE_ADJ.f
         ! inputs are WI, RHI, TEMPI, CNTRL, ADJ_GAS, ADJ_AERLIQ
         CALL ISOROPIAII_ADJ(WI, WI_ADJ, RHI, TEMPI,  CNTRL,
     &                       WT,  GAS, GAS_ADJ, AERLIQ, AERLIQ_ADJ,
     &                       AERSLD, SCASI, OTHER, TRUSTISO,NERR)

         IF ( TRUSTISO ) THEN              ! no ISOROPIAII_ADJ errors
            ! fwd code:
            !WI(1)    = MAX( TNA,  CONMIN )
            !WI(2)    = MAX( TSO4, CONMIN )
            !WI(3)    = MAX( TNH3, CONMIN )
            !WI(4)    = MAX( TNO3, CONMIN )
            !WI(5)    = MAX( TCL,  CONMIN )
            !WI(6)    = MAX( TCA,  CONMIN )
            !WI(7)    = MAX( TK,   CONMIN )
            !WI(8)    = MAX( TMG,  CONMIN )
            ! adj code (not sure if need all these, but include anyways to be complete):

            ! Modification for testing ANISO with no seasalt or dust
            ! adjoint - slc.4.2012

            TMG_ADJ  = 0.d0   ! WI_ADJ(8)
            TK_ADJ   = 0.d0   ! WI_ADJ(7)
            TCA_ADJ  = 0.d0   ! WI_ADJ(6)
            TCL_ADJ  = 0.d0   ! WI_ADJ(5)
            TNO3_ADJ = WI_ADJ(4)
            TNH3_ADJ = WI_ADJ(3)
            TSO4_ADJ = WI_ADJ(2)
            TNA_ADJ  = 0.d0   ! WI_ADJ(1)   - end of changes - slc.4.2012

            IF ( LPRINTFD
     &         .and. J == JFD .AND. L == LFD .AND. I == IFD) THEN
                WRITE(*,*) 'Adjoint execution, WI_ADJ'
                WRITE(*,*) 'SO4: ',WI_ADJ(2)
                WRITE(*,*) 'NH4: ',WI_ADJ(3)
                WRITE(*,*) 'NIT: ',WI_ADJ(4)
            ENDIF

            IF ( IDTHNO3 > 0 ) THEN

               ! fwd code:
               !TNO3    = GNO3 + ANO3
               ! adj code:
               GNO3_ADJ = TNO3_ADJ
               ANO3_ADJ = TNO3_ADJ

               ! fwd code:
               !ANO3     = STT(I,J,L,IDTNIT) * 1.d3 / ( 62.d0 * VOL )
               ! adj code:
               STT_ADJ(I,J,L,IDTNIT) = ANO3_ADJ * 1.d3 / ( 62.d0 * VOL )

               ! fwd code:
               !GNO3  = STT(I,J,L,IDTHNO3)
               !GNO3  = MAX( GNO3 * 1.d3 / ( 63.d0 * VOL ), CONMIN )
               GNO3_ADJ               = GNO3_ADJ * 1.d3 / ( 63.d0 * VOL)
               STT_ADJ(I,J,L,IDTHNO3) = GNO3_ADJ


             ELSE
               CALL ERROR_STOP('adj not supported for offline', X )
             ENDIF

             ! fwd code:
             !TCA      = 0d0
             !TK       = 0d0
             !TMG      = 0d0
             ! adj code:
             TCA_ADJ = 0d0
             TK_ADJ  = 0d0
             TMG_ADJ = 0d0

             ! Keep commented until seasalt adjoint is developed.
             ! (slc.1.2012)
             ! fwd code:
             !TCL      = STT(I,J,L,IDTSALA) * 0.5504d0 * 1.d3 /
             !                              ( 35.45d0  * VOL  )
             ! adj code (would add this once we have SALA ADJ:
             !STT_ADJ(I,J,L,IDTSALA) = TCL_ADJ * 0.5504d0 * 1.d3 /
             !                              ( 35.45d0  * VOL  )

             ! fwd code:
             !TNA      = STT(I,J,L,IDTSALA) * 0.3061d0 * 1.d3 /
             !                              ( 22.99d0  * VOL  )
             ! adj code
             !STT_ADJ(I,J,L,IDTSALA) = TNA_ADJ * 0.3061d0 * 1.d3 /
             !                              ( 22.99d0  * VOL  )
             TNA_ADJ = 0d0
             TCL_ADJ = 0d0

             !STT_ADJ(I,J,L,IDTDST1) = 0.d0 + STT_ADJ(I,J,L,IDTDST1) !TCA_ADJ
             !STT_ADJ(I,J,L,IDTDST2) = 0.d0 + STT_ADJ(I,J,L,IDTDST2) !TK_ADJ
             !STT_ADJ(I,J,L,IDTDST3) = 0.d0 + STT_ADJ(I,J,L,IDTDST3) !TMG_ADJ
             !STT_ADJ(I,J,L,IDTDST4) = 0.d0 + STT_ADJ(I,J,L,IDTDST4) !TNA_ADJ
             !STT_ADJ(I,J,L,IDTSALA) = 0.d0 + STT_ADJ(I,J,L,IDTSALA) !TCL_ADJ
             !STT_ADJ(I,J,L,IDTSALC) = 0.d0 + STT_ADJ(I,J,L,IDTSALC) !TCL_ADJ

             ! fwd code:

             !TNH3     = STT(I,J,L,IDTNH4) * 1.d3 / ( 18.d0 * VOL )  +
             !           STT(I,J,L,IDTNH3) * 1.d3 / ( 17.d0 * VOL )
             STT_ADJ(I,J,L,IDTNH4) =  TNH3_ADJ * 1.d3 / ( 18.d0 * VOL )

             STT_ADJ(I,J,L,IDTNH3) =  TNH3_ADJ * 1.d3 / ( 17.d0 * VOL )


             ! fwd code:
             !TSO4     = STT(I,J,L,IDTSO4) * 1.d3 / ( 96.d0 * VOL )
             ! adj code:
             STT_ADJ(I,J,L,IDTSO4) = TSO4_ADJ * 1.d3 / ( 96.d0 * VOL )

            IF ( LPRINTFD
     &         .and. J == JFD .AND. L == LFD .AND. I == IFD) THEN
                WRITE(*,*) 'After ANISORROPIA, STT_ADJ'
                WRITE(*,*) 'SO4 ',STT_ADJ(I,J,L,IDTSO4)
                WRITE(*,*) 'NH4 ',STT_ADJ(I,J,L,IDTNH4)
                WRITE(*,*) 'NH3 ',STT_ADJ(I,J,L,IDTNH3)
                WRITE(*,*) 'HNO3',STT_ADJ(I,J,L,IDTHNO3)
                WRITE(*,*) 'NIT ',STT_ADJ(I,J,L,IDTNIT)
                !PAUSE
            ENDIF

          ELSE

             ! Count the number of error flags & calls to reverse
             ANISOERRCOUNT = ANISOERRCOUNT + 1
             SELECT CASE (NERR)
                CASE (22)
                   NERR22 = NERR22 + 1
                CASE (33)
                   NERR33 = NERR33 + 1
                CASE (50)
                   NERR50 = NERR50 + 1
                CASE (100)
                   NERR100 = NERR100 + 1
                CASE (101)
                   NERR101 = NERR101 + 1
                CASE (102)
                   NERR102 = NERR102 + 1
                CASE (103)
                   NERR103 = NERR103 + 1
                CASE (104)
                   NERR104 = NERR104 + 1
                CASE DEFAULT
                   NERROTHER = NERROTHER + 1
             END SELECT

          ENDIF                     ! no ISOROPIAII_ADJ errors

             !  debug - slc.1.2012

             ANISOCALLCOUNT = ANISOCALLCOUNT + 1
             SELECT CASE (SCASI)
                CASE("A2")
                  CA = CA + 1
                CASE("B4")
                  CB = CB + 1
                CASE("C2")
                  CC = CC + 1
                CASE("D3")
                  CD = CD + 1
                CASE("E4")
                  CE = CE + 1
                CASE("F2")
                  CF = CF + 1
                CASE("G5")
                  CG = CG + 1
                CASE("H6")
                  CH = CH + 1
                CASE("I6")
                  CI = CI + 1
                CASE("J3")
                  CJ = CJ + 1
                CASE DEFAULT
                  COTHER = COTHER + 1
             END SELECT



          ! fwd code:
          !DO N = 1, NCOMPA
          !  WI(N) = 0d0
          !  WT(N) = 0d0
          !ENDDO
          ! adj code (reset values for safety)
          WI_ADJ(:)     = 0d0
          WT_ADJ(:)     = 0d0
          GAS_ADJ(:)    = 0d0
          AERLIQ_ADJ(:) = 0d0

         !<--- END ADJOINT CALCULATION

      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! More error checking: warn of exploding adjoit values, except
      ! the first jump up from zero (MAX_ADJ_TMP = 0 first few times)
      IF ( MAXVAL(ABS(STT_ADJ)) > (MAX_ADJ_TMP * MAX_ALLOWED_INCREASE)
     &   .AND. ( MAX_ADJ_TMP > 0d0 )  ) THEN

         WRITE(6,*)' *** - WARNING: EXPLODING adjoints in ADJ_AEROSOL'
         WRITE(6,*)' *** - MAX(ADJ_STT) before = ',MAX_ADJ_TMP
         WRITE(6,*)' *** - MAX(ADJ_STT) after  = ',MAXVAL(ABS(STT_ADJ))

         ADJ_EXPLD_COUNT = ADJ_EXPLD_COUNT + 1

         IF (ADJ_EXPLD_COUNT > MAX_ALLOWED_EXPLD )
     &       CALL ERROR_STOP('Too many exploding adjoints',
     &                       'ADJ_AEROSOL, adjoint_mod.f')

       ENDIF

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### ISOROPIAII_ADJ: AERO_THERMO_ADJ')
!      WRITE(6,*)  'ANISO calls: ',ANISOCALLCOUNT
!      WRITE(6,*)  'ANISO error occurrences: ',ANISOERRCOUNT
!      WRITE(6,*)  'Specific error codes: '
!      WRITE(6,*)  'Error 22: ',NERR22
!      WRITE(6,*)  'Error 33: ',NERR33
!      WRITE(6,*)  'Error 44: ',NERR44
!      WRITE(6,*)  'Error 100: ',NERR100
!      WRITE(6,*)  'Error 101: ',NERR101
!      WRITE(6,*)  'Error 102: ',NERR102
!      WRITE(6,*)  'Error 103: ',NERR103
!      WRITE(6,*)  'Error 104: ',NERR104
!      WRITE(6,*)  'Error Other: ', NERROTHER
!
!      WRITE(6,*) '____________ Case Distribution ____________'
!      WRITE(6,*) 'A: ',CA
!      WRITE(6,*) 'B: ',CB
!      WRITE(6,*) 'C: ',CC
!      WRITE(6,*) 'D: ',CD
!      WRITE(6,*) 'E: ',CE
!      WRITE(6,*) 'F: ',CF
!      WRITE(6,*) 'G: ',CG
!      WRITE(6,*) 'H: ',CH
!      WRITE(6,*) 'I: ',CI
!      WRITE(6,*) 'J: ',CJ
!      WRITE(6,*) 'Other: ', COTHER

      ! Return to calling program
      END SUBROUTINE DO_ISOROPIAII_ADJ
!EOC
!------------------------------------------------------------------------------
!         Caltech Department of Chemical Engineering / Seinfeld Group         !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: safelog10
!
! !DESCRIPTION: Calculates the LOG (base 10) of a number X.  Returns a minimum
!  value if X is too small, in order to avoid NaN or Infinity problems.
!\\
!\\
! !INTERFACE:
!
      FUNCTION SAFELOG10( X ) RESULT ( SAFLOG )
!
! !INPUT PARAMETERS:
!
      REAL*8, INTENT(IN) :: X        ! Argument for LOG10 function
!
! !RETURN VALUE:
!
      REAL*8             :: SAFLOG   ! LOG10 output --
!
! !REVISION HISTORY:
!  11 Aug 2009 - H. O. T. Pye - Initial version, in ISORROPIA II
!  29 Jan 2010 - R. Yantosca  - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC

      IF ( X <= 1d-20 ) THEN
          SAFLOG = -1d0*20d0   ! if X<0, make pH 20
      ELSE
          SAFLOG = LOG10(X)
      ENDIF

      END FUNCTION SAFELOG10
!EOC
!------------------------------------------------------------------------------
!         Caltech Department of Chemical Engineering / Seinfeld Group         !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_isrinfo
!
! !DESCRIPTION: Subroutine GET\_ISRINFO returns information related to
!  aerosol pH.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_ISRINFO( I, J, L, N ) RESULT ( RETURNVALUE )
!
! !INPUT PARAMETERS:
!
      INTEGER, INTENT(IN) :: I   ! GEOS-Chem longitude index
      INTEGER, INTENT(IN) :: J   ! GEOS-Chem latitude index
      INTEGER, INTENT(IN) :: L   ! GEOS-Chem level index
      INTEGER, INTENT(IN) :: N   ! Flag for which information is desired
!
! !RETURN VALUE:
!
      REAL*8              :: RETURNVALUE
!
! !REVISION HISTORY:
!  11 Aug 2009 - H. O. T. Pye - Initial version
!  29 Jan 2010 - R. Yantosca  - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      IF     ( N == 1 ) THEN
         RETURNVALUE = PH_SAV( I, J, L )
      ELSEIF ( N == 2 ) THEN
         RETURNVALUE = HPLUS_SAV( I, J, L )
      ELSEIF ( N == 3 ) THEN
         RETURNVALUE = WATER_SAV( I, J, L )
      ELSEIF ( N == 4 ) THEN
         RETURNVALUE = SULRAT_SAV( I, J, L )
      ELSEIF ( N == 5 ) THEN
         RETURNVALUE = NARAT_SAV( I, J, L )
      ELSEIF ( N == 6 ) THEN
         RETURNVALUE = ACIDPUR_SAV( I, J, L )
      ELSE
         ! return large value to indicate problem
         RETURNVALUE = 99999d0
         !FP_ISOP
         WRITE(*,*) 'VALUE NOT DEFINED IN GET_ISRINFO'
      ENDIF

      END FUNCTION GET_ISRINFO
!EOC
!------------------------------------------------------------------------------
!         Caltech Department of Chemical Engineering / Seinfeld Group         !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_hno3
!
! !DESCRIPTION: Subroutine GET\_HNO3 allows the HNO3 concentrations to evolve
!  with time, but relaxes back to the monthly mean concentrations every 3
!  hours.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_HNO3( I, J, L ) RESULT ( HNO3_UGM3 )
!
! !USES:
!
      USE GLOBAL_HNO3_MOD, ONLY : GET_HNO3_UGM3
      USE TIME_MOD,        ONLY : GET_ELAPSED_MIN
!
! !INPUT PARAMETERS:
!
      INTEGER, INTENT(IN) :: I  ! GEOS-Chem longitude index
      INTEGER, INTENT(IN) :: J  ! GEOS-Chem latitude index
      INTEGER, INTENT(IN) :: L  ! GEOS-Chem level index
!
! !REVISION HISTORY:
!  16 Dec 2002 - R. Yantosca  - Initial version, in ISORROPIA I
!  24 Mar 2003 - R. Yantosca  - Now use function GET_ELAPSED_MIN() from the
!                               new "time_mod.f" to get the elapsed minutes
!                               since the start of run.
!  06 Jul 2007 - H. O. T. Pye - Initial version, in ISORROPIA II
!  29 Jan 2010 - R. Yantosca  - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*8 :: HNO3_UGM3

      !=================================================================
      ! GET_HNO3 begins here!
      !=================================================================

      ! Relax to monthly mean HNO3 concentrations every 3 hours
      ! Otherwise just return the concentration in HNO3_sav
      IF ( MOD( GET_ELAPSED_MIN(), 180 ) == 0 ) THEN
         HNO3_UGM3 = GET_HNO3_UGM3( I, J, L )
      ELSE
         HNO3_UGM3 = HNO3_sav(I,J,L)
      ENDIF

      ! Return to calling program
      END FUNCTION GET_HNO3
!EOC
!------------------------------------------------------------------------------
!         Caltech Department of Chemical Engineering / Seinfeld Group         !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_hno3
!
! !DESCRIPTION: Subroutine SET\_HNO3 stores the modified HNO3 value back
!  into the HNO3\_sav array for the next timestep.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SET_HNO3( I, J, L, HNO3_UGM3 )
!
! !INPUT PARAMETERS:
!
      INTEGER, INTENT(IN) :: I           ! GEOS-Chem longitude index
      INTEGER, INTENT(IN) :: J           ! GEOS-Chem longitude index
      INTEGER, INTENT(IN) :: L           ! GEOS-Chem longitude index
      REAL*8,  INTENT(IN) :: HNO3_UGM3   ! HNO3 concentration [ug/m3]
!
! !REVISION HISTORY:
!  16 Dec 2002 - R. Yantosca  - Initial version, in ISORROPIA I
!  06 Jul 2007 - H. O. T. Pye - Initial version, in ISORROPIA II
!  29 Jan 2010 - R. Yantosca  - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      HNO3_sav(I,J,L) = HNO3_UGM3

      END SUBROUTINE SET_HNO3
!EOC
!------------------------------------------------------------------------------
!         Caltech Department of Chemical Engineering / Seinfeld Group         !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_gno3
!
! !DESCRIPTION: Function GET\_GNO3 returns the gas-phase HNO3 [v/v] for
!  calculation of sea-salt chemistry in sulfate\_mod (SEASALT\_CHEM).
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GET_GNO3( I, J, L, HNO3_kg )
!
! !USES:
!
      USE DAO_MOD, ONLY : AIRVOL, AD
!
! !INPUT PARAMETERS:
!
      INTEGER, INTENT(IN)  :: I       ! GEOS-Chem longitude index
      INTEGER, INTENT(IN)  :: J       ! GEOS-Chem latitude index
      INTEGER, INTENT(IN)  :: L       ! GEOS-Chem level index
!
! !OUTPUT PARAMETERS:
!
      REAL*8,  INTENT(OUT) :: HNO3_kg ! Gas-phase HNO3 [kg]
!
! !REVISION HISTORY:
!  15 Apr 2005 - B. Alexander - Initial version, in ISORROPIA I
!  06 Jul 2007 - H. O. T. Pye - Initial version, in ISORROPIA II
!  29 Jan 2010 - R. Yantosca  - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      ! Zero variables
      HNO3_kg  = 0.D0

      ! convert from [mole/m3] to [kg]
      HNO3_kg = GAS_HNO3(I,J,L) * 63.d-3 * AIRVOL(I,J,L)

      ! Return to calling program
      END SUBROUTINE GET_GNO3
!EOC
!------------------------------------------------------------------------------
!         Caltech Department of Chemical Engineering / Seinfeld Group         !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_isoropiaII
!
! !DESCRIPTION: Subroutine INIT\_ISOROPIAII initializes all module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_ISOROPIAII
!
! !USES:
!
      USE ERROR_MOD, ONLY : ALLOC_ERR

#     include "CMN_SIZE"  ! Size parameters
!
! !REVISION HISTORY:
!  06 Jul 2007 - H. O. T. Pye - Initial version
!  29 Jan 2010 - R. Yantosca  - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: AS

      !=================================================================
      ! INIT_ISOROPIAII begins here!
      !=================================================================

      WRITE(*,*) 'INIT_ISOROPIAII'

      ALLOCATE( HNO3_sav( IIPAR, JJPAR, LLTROP ) , STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'HNO3_sav' )
      HNO3_sav = 0d0

      ALLOCATE( GAS_HNO3( IIPAR, JJPAR, LLTROP ) , STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'GAS_HNO3' )
      GAS_HNO3 = 0d0

      ! AEROPH: diagnostic info (hotp 8/11/09)
      ! Allocate up to LLPAR, but zero above LLTROP
      ALLOCATE( PH_SAV( IIPAR, JJPAR, LLPAR ) , STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'PH_SAV' )
      PH_SAV = 0d0

      ALLOCATE( HPLUS_SAV( IIPAR, JJPAR, LLPAR ) , STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'HPLUS_SAV' )
      HPLUS_SAV = 0d0

      ALLOCATE( WATER_SAV( IIPAR, JJPAR, LLPAR ) , STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'WATER_SAV' )
      WATER_SAV = 0d0

      ALLOCATE( SULRAT_SAV( IIPAR, JJPAR, LLPAR ) , STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SULRAT_SAV' )
      SULRAT_SAV = 0d0

      ALLOCATE( NARAT_SAV( IIPAR, JJPAR, LLPAR ) , STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'NARAT_SAV' )
      NARAT_SAV = 0d0

      ALLOCATE( ACIDPUR_SAV( IIPAR, JJPAR, LLPAR ) , STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ACIDPUR_SAV' )
      ACIDPUR_SAV = 0d0

      ! Return to calling program
      END SUBROUTINE INIT_ISOROPIAII
!EOC
!------------------------------------------------------------------------------
!         Caltech Department of Chemical Engineering / Seinfeld Group         !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_isoropiaII
!
! !DESCRIPTION: Subroutine CLEANUP\_ISOROPIAII deallocates all module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_ISOROPIAII
!
! !REVISION HISTORY:
!  06 Jul 2007 - H. O. T. Pye - Initial version
!  29 Jan 2010 - R. Yantosca  - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC

      IF ( ALLOCATED( HNO3_sav    ) ) DEALLOCATE( HNO3_sav )
      IF ( ALLOCATED( GAS_HNO3    ) ) DEALLOCATE( GAS_HNO3 )
      ! AEROPH: Deallocate arrays for pH (hotp 8/11/09)
      IF ( ALLOCATED( PH_SAV      ) ) DEALLOCATE( PH_SAV     )
      IF ( ALLOCATED( HPLUS_SAV   ) ) DEALLOCATE( HPLUS_SAV  )
      IF ( ALLOCATED( WATER_SAV   ) ) DEALLOCATE( WATER_SAV  )
      IF ( ALLOCATED( SULRAT_SAV  ) ) DEALLOCATE( SULRAT_SAV )
      IF ( ALLOCATED( NARAT_SAV   ) ) DEALLOCATE( NARAT_SAV  )
      IF ( ALLOCATED( ACIDPUR_SAV ) ) DEALLOCATE( ACIDPUR_SAV)

      END SUBROUTINE CLEANUP_ISOROPIAII
!EOC
      END MODULE ISOROPIAII_ADJ_MOD
