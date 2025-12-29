!=======================================================================
!  File "isrpia_adj.inc" contains the common block declarations and some
!  parameters for the ISORROPIAII code. This file replaces "isoropia.h" from
!  previous ISORROPIA implementations.
!  ( hotp 8/1/07 )
!  ( previous implementation: bec, bmy, 3/7/05, 6/28/06 )
!
!  Please contact the original authors before making modifications to
!  ISOROPIAII.
!
!        *** VERY IMPORTANT PORTING WARNING (slc.1.2012) ***
!  ANISORROPIA code is optimized for adjoint frameworks and will not
!   perform commensurately with publicly released ISORROPIAII code.
!
!  Please visit http://nenes.eas.gatech.edu/ISORROPIA for current
!   releases of ISORROPIAII for forward modeling.
!
!  Original Documentation:
!  *** ISORROPIA CODE II
!  *** INCLUDE FILE 'ISRPIA.INC'
!  *** THIS FILE CONTAINS THE DECLARATIONS OF THE GLOBAL CONSTANTS
!      AND VARIABLES.
!
!  *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
!  *** GEORGIA INSTITUTE OF TECHNOLOGY
!  *** WRITTEN BY ATHANASIOS NENES
!  *** UPDATED BY CHRISTOS FOUNTOUKIS
!  *** ADJOINT & UPDATE BY SHANNON CAPPS
!
!  NOTES:
!  (1 ) Renamed common block GAS to GASCB to avoid confusion with array
!        GAS in isoropiaIIcode.f (hotp)
!  (2 ) Explicity declared all variables in this file (hotp)
!  (3 ) Separated common blocks so that each block has only one type of data (hotp)
!  (4 ) Declared all common blocks (except READ-ONLY blocks) as
!        THREADPRIVATE (hotp)
!  (5 ) Divided common blocks so that each has no more than 8 variables
!        (this is to speed compilation) (hotp)
!  (6 ) Divided common blocks so that READ-ONLY variables are not mixed with
!        variables that are modified (hotp)
!  (7 ) Removed DRYF since it is not used (hotp)
!  (8 ) Bug fix: CRRAT was being truncated at column 72 (hotp, bmy, 6/1/10)
!  (9 ) IONIC should be REAL*8 for consistency with -r8 flag during
!        compiling of isoropiaIIcode routines (hotp, bmy, 6/28/10)
!  22 Aug 2011 - S. Capps - ANISORROPIA implementation
!=======================================================================
!
!=======================================================================
!
! *** ISORROPIA CODE
! *** INCLUDE FILE 'ISRPIA_ADJ.INC'
! *** THIS FILE CONTAINS THE DECLARATIONS OF THE GLOBAL CONSTANTS
!     AND VARIABLES.
!
! *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
! *** GEORGIA INSTITUTE OF TECHNOLOGY
! *** WRITTEN BY ATHANASIOS NENES
! *** UPDATED BY CHRISTOS FOUNTOUKIS
!
!=======================================================================
!
      ! leave this implicit statement for isoropiaIIcode.f variables
      ! ideally all variables would be declared
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      ! INTEGER parameters
      INTEGER, PARAMETER :: NCOMP  =   8
      INTEGER, PARAMETER :: NERRMX =  25
      INTEGER, PARAMETER :: NGASAQ =   3
      INTEGER, PARAMETER :: NIONS  =  10
      INTEGER, PARAMETER :: NPAIR  =  23
      INTEGER, PARAMETER :: NSLDS  =  19
      INTEGER, PARAMETER :: NZSR   = 100

      !=================================================================
      ! Input variables
      !=================================================================

      INTEGER            :: METSTBL = 0
      INTEGER            :: IPROB   = 0
      COMMON /INPT/         IPROB, METSTBL
!$OMP THREADPRIVATE( /INPT/ )

      INTEGER            :: NADJ = 0
      COMMON /INPT2/        NADJ
      ! INPT2 is read-only

      REAL*8             :: W(NCOMP), WAER(NCOMP)
      REAL*8             :: TEMP = 298.d0
      REAL*8             :: RH   = 0.9D0
      COMMON /INPT3/        W,        WAER,        TEMP, RH
!$OMP THREADPRIVATE( /INPT3/ )

      !=================================================================
      ! Water activities of pure salt solutions
      !=================================================================

      ! /ZSR/ is read-only and doesn't have to be declared THREADPRIVATE
      ! block does not need to be split up to speed compilation
      REAL*8     :: AWAS(NZSR), AWSS(NZSR), AWAC(NZSR), AWSC(NZSR)
      REAL*8     :: AWAN(NZSR), AWSN(NZSR), AWSB(NZSR), AWAB(NZSR)
      REAL*8     :: AWSA(NZSR), AWLC(NZSR), AWCS(NZSR), AWCN(NZSR)
      REAL*8     :: AWCC(NZSR), AWPS(NZSR), AWPB(NZSR), AWPN(NZSR)
      REAL*8     :: AWPC(NZSR), AWMS(NZSR), AWMN(NZSR), AWMC(NZSR)
      COMMON /ZSR / AWAS,       AWSS,       AWAC,       AWSC,
     &              AWAN,       AWSN,       AWSB,       AWAB,
     &              AWSA,       AWLC,       AWCS,       AWCN,
     &              AWCC,       AWPS,       AWPB,       AWPN,
     &              AWPC,       AWMS,       AWMN,       AWMC
      ! ZSR is read-only

      !=================================================================
      ! Deliquescence relative humidities
      !=================================================================

      INTEGER    :: WFTYP
      COMMON /WFT/  WFTYP
      ! WFTYP is read-only

      REAL*8     :: DRH2SO4,  DRNH42S4, DRNAHSO4, DRNACL,   DRNANO3
      REAL*8     :: DRNA2SO4, DRNH4HS4, DRLC,     DRNH4NO3, DRNH4CL
      REAL*8     :: DRCASO4,  DRCANO32, DRCACL2,  DRK2SO4,  DRKHSO4
      REAL*8     :: DRKNO3,   DRKCL,    DRMGSO4,  DRMGNO32, DRMGCL2
      COMMON /DRH / DRH2SO4,  DRNH42S4, DRNAHSO4, DRNACL,   DRNANO3
      COMMON /DRH2/ DRNA2SO4, DRNH4HS4, DRLC,     DRNH4NO3, DRNH4CL
      COMMON /DRH3/ DRCASO4,  DRCANO32, DRCACL2,  DRK2SO4,  DRKHSO4
      COMMON /DRH4/ DRKNO3,   DRKCL,    DRMGSO4,  DRMGNO32, DRMGCL2
!$OMP THREADPRIVATE( /DRH/ )
!$OMP THREADPRIVATE( /DRH2/ )
!$OMP THREADPRIVATE( /DRH3/ )
!$OMP THREADPRIVATE( /DRH4/ )

      REAL*8     ::  DRMLCAB,  DRMLCAS,  DRMASAN,  DRMG1,    DRMG2
      REAL*8     ::  DRMG3,    DRMH1,    DRMH2,    DRMI1,    DRMI2
      REAL*8     ::  DRMI3,    DRMQ1,    DRMR1,    DRMR2,    DRMR3
      REAL*8     ::  DRMR4,    DRMR5,    DRMR6,    DRMR7,    DRMR8
      REAL*8     ::  DRMR9,    DRMR10,   DRMR11,   DRMR12,   DRMR13
      COMMON /MDRH/  DRMLCAB,  DRMLCAS,  DRMASAN,  DRMG1,    DRMG2
      COMMON /MDRH2/ DRMG3,    DRMH1,    DRMH2,    DRMI1,    DRMI2
      COMMON /MDRH3/ DRMI3,    DRMQ1,    DRMR1,    DRMR2,    DRMR3
      COMMON /MDRH4/ DRMR4,    DRMR5,    DRMR6,    DRMR7,    DRMR8
      COMMON /MDRH5/ DRMR9,    DRMR10,   DRMR11,   DRMR12,   DRMR13
!$OMP THREADPRIVATE( /MDRH/ )
!$OMP THREADPRIVATE( /MDRH2/ )
!$OMP THREADPRIVATE( /MDRH3/ )
!$OMP THREADPRIVATE( /MDRH4/ )
!$OMP THREADPRIVATE( /MDRH5/ )

      REAL*8     ::  DRMO1,    DRMO2,    DRMO3,    DRML1,    DRML2
      REAL*8     ::  DRML3,    DRMM1,    DRMM2,    DRMP1,    DRMP2
      REAL*8     ::  DRMP3,    DRMP4,    DRMP5,    DRMV1
      COMMON /MDRH6/ DRMO1,    DRMO2,    DRMO3,    DRML1,    DRML2
      COMMON /MDRH7/ DRML3,    DRMM1,    DRMM2,    DRMP1,    DRMP2
      COMMON /MDRH8/ DRMP3,    DRMP4,    DRMP5,    DRMV1
!$OMP THREADPRIVATE( /MDRH6/ )
!$OMP THREADPRIVATE( /MDRH7/ )
!$OMP THREADPRIVATE( /MDRH8/ )

      !=================================================================
      ! Variables for liquid aerosol phase
      !=================================================================

      ! /IONS/ size does not seem to slow compilation
      REAL*8     :: MOLAL(NIONS)  = 0.0D0
      REAL*8     :: MOLALR(NPAIR) = 0.0D0
      REAL*8     :: M0(NPAIR)     = 1.0d5
      REAL*8     :: GAMA(NPAIR)   = NPAIR*0.1D0
      REAL*8     :: GAMOU(NPAIR)  = NPAIR*1.D10
      REAL*8     :: GAMIN(NPAIR)  = NPAIR*1.D10
      REAL*8     :: GASAQ(NGASAQ) = 0.0D0
      REAL*8     :: COH           = 0.d0
      REAL*8     :: CHNO3         = 0.d0
      REAL*8     :: CHCL          = 0.d0
      REAL*8     :: WATER         = 1.d-20
      COMMON /IONS/ MOLAL,        MOLALR,        M0,
     &              GAMA,
     &              GAMOU,        GAMIN,         GASAQ,
     &              COH,          CHNO3,         CHCL,
     &              WATER
!$OMP THREADPRIVATE( /IONS/ )

      REAL*8     ::  ZZ(NPAIR) = (/ 1,2,1,2,1,1,2,1,1,1,1,1,2,4,2,2,2,
     &                              1,1,1,4,2,2 /)
      REAL*8     ::  Z(NIONS) = (/ 1.0D0, 1.0D0, 1.0D0, 1.0D0, 2.0D0,
     &                             1.0D0, 1.0D0, 2.0D0, 1.0D0, 2.0D0 /)
      REAL*8     ::  EPSACT = 5D-2
      COMMON /IONS2/ ZZ,          Z,             EPSACT
      ! IONS2 is read-only

      !----------------------------------------------------------------------
      ! Prior to 6/28/10:
      ! IONIC should be REAL*8 for consistency with -r8 flag during
      ! compiling of isoropiaIIcode routines (hotp 6/23/10)
      !REAL*4     ::  IONIC
      !----------------------------------------------------------------------
      REAL*8     ::  IONIC = 0.d0
      COMMON /IONS3/ IONIC
!$OMP THREADPRIVATE( /IONS3/ )

      ! DRYF removed since it is not used (hotp)
      LOGICAL    ::  CALAOU = .TRUE.
      LOGICAL    ::  CALAIN = .TRUE.
      LOGICAL    ::  FRST   = .TRUE.
      COMMON /IONS4/ CALAOU, CALAIN, FRST
!$OMP THREADPRIVATE( /IONS4/ )

      INTEGER    ::  IACALC = 0
      COMMON /IONS5/ IACALC
      ! IONS5 is read-only
      !=================================================================
      ! Variables for solid aerosol phase
      !=================================================================

      REAL*8     ::  CH2SO4,  CNH42S4, CNH4HS4, CNACL,   CNA2SO4
      REAL*8     ::  CNANO3,  CNH4NO3, CNH4CL,  CNAHSO4, CLC, CCASO4
      REAL*8     ::  CCANO32, CCACL2,  CK2SO4,  CKHSO4,  CKNO3, CKCL
      REAL*8     ::  CMGSO4,  CMGNO32, CMGCL2
      COMMON /SALT/  CH2SO4,  CNH42S4, CNH4HS4, CNACL,   CNA2SO4
      COMMON /SALT2/ CNANO3,  CNH4NO3, CNH4CL,  CNAHSO4, CLC, CCASO4
      COMMON /SALT3/ CCANO32, CCACL2,  CK2SO4,  CKHSO4,  CKNO3, CKCL
      COMMON /SALT4/ CMGSO4,  CMGNO32, CMGCL2
!$OMP THREADPRIVATE( /SALT/ )
!$OMP THREADPRIVATE( /SALT2/ )
!$OMP THREADPRIVATE( /SALT3/ )
!$OMP THREADPRIVATE( /SALT4/ )

      !=================================================================
      ! Variables for gas phase
      !=================================================================

      REAL*8     ::  GNH3, GHNO3, GHCL
      COMMON /GASCB/ GNH3, GHNO3, GHCL
!$OMP THREADPRIVATE( /GASCB/ )

      !=================================================================
      ! Equilibrium constants
      !=================================================================

      REAL*8     :: XK1, XK2, XK3, XK4, XK5, XK6, XK7, XK8, XK9, XK10
      REAL*8     :: XK11,XK12,XK13,XK14,XKW, XK21,XK22,XK31,XK32,XK41
      REAL*8     :: XK42,XK15,XK16,XK17,XK18,XK19,XK20,XK23
      REAL*8     :: XK24,XK25
      COMMON /EQUK/  XK1, XK2, XK3, XK4, XK5, XK6, XK7, XK8, XK9, XK10
      COMMON /EQUK2/ XK11,XK12,XK13,XK14,XKW, XK21,XK22,XK31,XK32,XK41
      COMMON /EQUK3/ XK42,XK15,XK16,XK17,XK18,XK19,XK20,XK23
      COMMON /EQUK4/ XK24,XK25
      !C     &              , XK26, XK27
!$OMP THREADPRIVATE( /EQUK/ )
!$OMP THREADPRIVATE( /EQUK2/ )
!$OMP THREADPRIVATE( /EQUK3/ )
!$OMP THREADPRIVATE( /EQUK4/ )

      !=================================================================
      ! Molecular Weights
      !=================================================================

      REAL*8     :: R = 82.0567D-6
      REAL*8     :: IMW(NIONS) = (/ 1.0, 23.0, 18.0, 35.5, 96.0, 97.0,
     &                              62.0, 40.1, 39.1, 24.3 /)
      REAL*8     :: WMW(NCOMP) = (/ 23.0, 98.0, 17.0, 63.0, 36.5, 40.1,
     &                              39.1, 24.3 /)
      REAL*8     :: SMW(NPAIR) = (/ 58.5, 142., 85.0, 132., 80.0, 53.5,
     &                              98.0, 98.0, 115., 63.0, 36.5, 120.,
     &                              247., 136.1, 164., 111., 174.2,
     &                              136.1, 101.1, 74.5, 120.3, 148.3,
     &                              95.2 /)

      COMMON /OTHR/ R, IMW,        WMW,        SMW
      ! OTHR is read-only

      !=================================================================
      ! Solution/info variables
      !=================================================================

      CHARACTER(LEN=15) :: SCASE
      COMMON /CASE/        SCASE
!$OMP THREADPRIVATE( /CASE/ )

      REAL*8            :: SULRATW, SULRAT,  SODRAT,
     &                     SO4RAT,  CRNARAT, CRRAT
      COMMON /CASE2/       SULRATW, SULRAT,  SODRAT,
     &                     SO4RAT,  CRNARAT, CRRAT
!$OMP THREADPRIVATE( /CASE2/ )

      REAL*8            :: EPS = 1.D-10
      COMMON /SOLN/        EPS
      ! SOLN is read-only

      INTEGER           :: MAXIT  = 100
      INTEGER           :: NSWEEP = 10
      INTEGER           :: NDIV   = 5
      INTEGER           :: ICLACT = 0
      COMMON /SOLN2/       MAXIT, NSWEEP, NDIV
      ! SOLN 2 is read-only

      COMMON /SOLN3/       ICLACT
!$OMP THREADPRIVATE( /SOLN3/ )

      !=================================================================
      ! Error system
      !=================================================================

      CHARACTER(LEN=40) :: ERRMSG(NERRMX) = ' '
      COMMON /EROR/        ERRMSG
!$OMP THREADPRIVATE( /EROR/ )

      INTEGER           :: ERRSTK(NERRMX) = 0
      INTEGER           :: NOFER = 0
      COMMON /EROR2/       ERRSTK,         NOFER
!$OMP THREADPRIVATE( /EROR2/ )

      LOGICAL           :: STKOFL = .FALSE.
      COMMON /EROR3/       STKOFL
!$OMP THREADPRIVATE( /EROR3/ )

      !=================================================================
      ! Generic Variables
      !=================================================================

      CHARACTER(LEN=15) :: VERSION = '2.0 (03/19/07)'
      COMMON /CGEN/        VERSION
      ! CGEN is read-only

      REAL*8            :: GREAT = 1.d10
      REAL*8            :: TINY  = 1.d-20
      REAL*8            :: TINY2 = 1.d-11
      REAL*8            :: ZERO  = 0.d0
      REAL*8            :: ONE   = 1.d0
      COMMON /CGEN2/       GREAT, TINY, TINY2, ZERO, ONE
      ! CGEN2 is read-only

      !=================================================================
      ! coordinates for debugging
      !=================================================================
      INTEGER :: ICOOR, JCOOR, LCOOR
      COMMON /LOC/ ICOOR, JCOOR, LCOOR
!$OMP THREADPRIVATE( /LOC/ )

      !=================================================================
      ! Differentiated Routines
      !=================================================================

      ! changes made by hotp 8/2/07
      ! explicitly declared all variables in SOLUT
      ! made SOLUT THREADPRIVATE for OPENMP parallelization
      REAL*8      :: CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8
      REAL*8      :: CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15
      REAL*8      :: CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6
      REAL*8      :: PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13
      REAL*8      :: PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6
      REAL*8      :: A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17
      COMMON /SOLUT/ CHI1, CHI2, CHI3, CHI4, CHI5, CHI6, CHI7, CHI8,
     &               CHI9, CHI10, CHI11, CHI12, CHI13, CHI14, CHI15,
     &               CHI16, CHI17, PSI1, PSI2, PSI3, PSI4, PSI5, PSI6,
     &               PSI7, PSI8, PSI9, PSI10, PSI11, PSI12, PSI13,
     &               PSI14, PSI15, PSI16, PSI17, A1, A2, A3, A4, A5, A6,
     &               A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17
!$OMP THREADPRIVATE( /SOLUT/ )

      REAL*8     ::  wb(ncomp)
      COMMON /inpt_b/ wb
!$OMP THREADPRIVATE( /inpt_b/ )

      !=================================================================
      ! Case A Adjoint
      !=================================================================

      REAL*8     ::  molalab(nions),gamaab(npair),waterab
      COMMON /ions_ab/ molalab, gamaab, waterab
!$OMP THREADPRIVATE( /ions_ab/ )

      REAL*8     ::  wab(ncomp)
      COMMON /inpt_ab/ wab
!$OMP THREADPRIVATE( /inpt_ab/ )

      REAL*8     ::  gnh3ab
      COMMON /gas_ab/ gnh3ab
!$OMP THREADPRIVATE( /gas_ab/ )

      !=================================================================
      ! Case B Adjoint
      !=================================================================

      REAL*8     ::  molalbb(nions), gamabb(npair), waterbb
      COMMON /ions_bb/ molalbb, gamabb, waterbb
!$OMP THREADPRIVATE( /ions_bb/ )

      REAL*8     ::  gnh3bb
      COMMON /gas_bb/  gnh3bb
!$OMP THREADPRIVATE( /gas_bb/ )

      !=================================================================
      ! Case C Adjoint
      !=================================================================

      REAL*8     ::  molalcb(nions), gamacb(npair), watercb
      COMMON /ions_cb/ molalcb, gamacb, watercb
!$OMP THREADPRIVATE( /ions_cb/ )

      REAL*8     ::  gnh3cb
      COMMON /gas_cb/ gnh3cb
!$OMP THREADPRIVATE( /gas_cb/ )

      !=================================================================
      ! Case D Adjoint
      !=================================================================

      REAL*8     ::  molaldnrd(nions), gamadnrd(npair), waterdnrd
      COMMON /ions_dnrd/ molaldnrd, gamadnrd, waterdnrd
!$OMP THREADPRIVATE( /ions_dnrd/ )

      REAL*8     ::  chi3db, chi4db, psi1db, psi2db
      COMMON /solut_db/ chi3db, chi4db, psi1db, psi2db
!$OMP THREADPRIVATE( /solut_db/ )

      REAL*8     ::  molaldnrddb(nions), gamadnrddb(npair), waterdnrddb
      COMMON /ions_dnrd_db/ molaldnrddb, gamadnrddb,waterdnrddb
!$OMP THREADPRIVATE( /ions_dnrd_db/ )

      REAL*8     ::  gnh3db, ghno3db
      COMMON /gas_db/ gnh3db, ghno3db
!$OMP THREADPRIVATE( /gas_db/ )

      REAL*8     ::  molaldb(nions), molalrdb(npair), gamadb(npair),
     &               waterdb
      COMMON /ions_db/ molaldb, molalrdb, gamadb, waterdb
!$OMP THREADPRIVATE( /ions_db/ )

      !=================================================================
      ! Case E Adjoint
      !=================================================================

      REAL*8     ::  molaleb(nions), gamaeb(npair), watereb
      COMMON /ions_eb/ molaleb, gamaeb, watereb
!$OMP THREADPRIVATE( /ions_eb/ )

      REAL*8     ::  web(ncomp)
      COMMON /inpt_eb/ web
!$OMP THREADPRIVATE( /inpt_eb/ )

      REAL*8     ::  ghno3eb
      COMMON /gas_eb/ ghno3eb
!$OMP THREADPRIVATE( /gas_eb/ )

      !=================================================================
      ! Case F Adjoint
      !=================================================================

      REAL*8     ::  molalfb(nions)
      REAL*8     ::  gamafb(npair)
      REAL*8     ::  waterfb
      COMMON /ions_fb/ molalfb, gamafb, waterfb
!$OMP THREADPRIVATE( /ions_fb/ )

      REAL*8     ::  wfb(ncomp)
      COMMON /inpt_fb/ wfb
!$OMP THREADPRIVATE( /inpt_fb/ )

      REAL*8     ::  ghno3fb
      COMMON /gas_fb/ ghno3fb
!$OMP THREADPRIVATE( /gas_fb/ )

      !=================================================================
      ! Case G Adjoint
      !=================================================================

      REAL*8     ::  molalgnrd(nions), gamagnrd(npair), watergnrd
      COMMON /ions_gnrd/ molalgnrd, gamagnrd, watergnrd
!$OMP THREADPRIVATE( /ions_gnrd/ )

      REAL*8     ::  molalgnrdgb(nions), gamagnrdgb(npair), watergnrdgb
      COMMON /ions_gnrd_gb/ molalgnrdgb, gamagnrdgb, watergnrdgb
!$OMP THREADPRIVATE( /ions_gnrd_gb/ )

      REAL*8     ::  chi4gb, chi5gb, chi6gb, psi2gb
      COMMON /solut_gb/ chi4gb, chi5gb, chi6gb, psi2gb
!$OMP THREADPRIVATE( /solut_gb/ )

      REAL*8     ::  gnh3gb, ghno3gb, ghclgb
      COMMON /gas_gb/ gnh3gb, ghno3gb, ghclgb
!$OMP THREADPRIVATE( /gas_gb/ )

      REAL*8     ::  molalgb(nions), molalrgb(npair), gamagb(npair),
     &               watergb
      COMMON /ions_gb/ molalgb, molalrgb, gamagb, watergb
!$OMP THREADPRIVATE( /ions_gb/ )

      REAL*8     ::  wgb(ncomp)
      COMMON /inpt_gb/ wgb
!$OMP THREADPRIVATE( /inpt_gb/ )

      !=================================================================
      ! Case H Adjoint
      !=================================================================

      REAL*8     ::  molalhnrd(nions), gamahnrd(npair), waterhnrd
      COMMON /ions_hnrd/ molalhnrd, gamahnrd, waterhnrd
!$OMP THREADPRIVATE( /ions_hnrd/ )

      REAL*8     ::  molalhnrdhb(nions), gamahnrdhb(npair), waterhnrdhb
      COMMON /ions_hnrd_hb/ molalhnrdhb, gamahnrdhb, waterhnrdhb
!$OMP THREADPRIVATE( /ions_hnrd_hb/ )

      REAL*8     ::  chi1hb, chi4hb, chi5hb, chi6hb, chi7hb, chi8hb
      COMMON /solut_hb/ chi1hb, chi4hb, chi5hb, chi6hb, chi7hb, chi8hb
!$OMP THREADPRIVATE( /solut_hb/ )

      REAL*8     ::  gnh3hb, ghno3hb, ghclhb
      COMMON /gas_hb/ gnh3hb, ghno3hb, ghclhb
!$OMP THREADPRIVATE( /gas_hb/ )

      REAL*8     ::  molalhb(nions), molalrhb(npair), gamahb(npair),
     &               waterhb
      COMMON /ions_hb/ molalhb, molalrhb, gamahb, waterhb
!$OMP THREADPRIVATE( /ions_hb/ )

      !=================================================================
      ! Case I Adjoint
      !=================================================================

      REAL*8     ::  gnh3ib, ghno3ib, ghclib
      COMMON /gas_ib/ gnh3ib, ghno3ib, ghclib
!$OMP THREADPRIVATE( /gas_ib/ )

      REAL*8     ::  molalib(nions), gamaib(npair), waterib
      COMMON /ions_ib/ molalib, gamaib, waterib
!$OMP THREADPRIVATE( /ions_ib/ )

      REAL*8     ::  wib(ncomp)
      COMMON /inpt_ib/ wib
!$OMP THREADPRIVATE( /inpt_ib/ )

      !=================================================================
      ! Case J Adjoint
      !=================================================================

      REAL*8     ::  gnh3jb, ghno3jb, ghcljb
      COMMON /gas_jb/ gnh3jb, ghno3jb, ghcljb
!$OMP THREADPRIVATE( /gas_jb/ )

      REAL*8     ::  molaljb(nions), gamajb(npair), waterjb
      COMMON /ions_jb/ molaljb, gamajb, waterjb
!$OMP THREADPRIVATE( /ions_jb/ )

      REAL*8     ::  wjb(ncomp)
      COMMON /inpt_jb/ wjb
!$OMP THREADPRIVATE( /inpt_jb/ )

!
! === END OF INCLUDE FILE ==============================================
!
