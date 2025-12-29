! $Id: global_ch4_mod.f,v 1.1 2012/03/01 22:00:27 daven Exp $
      MODULE GLOBAL_CH4_MOD
!
!******************************************************************************
!  Module GLOBAL_CH4_MOD contains variables and routines for simulating
!  CH4 chemistry in the troposphere (jsw, bnd, bmy, 1/17/01, 10/1/09)
!
!  Module Variables:
!  ===========================================================================
!  (1 ) N_CH4      (INTEGER) : Number of budget items in TCH4
!  (2 ) BAIRDENS   (REAL*8 ) : Array for air density [molec/cm3]
!  (3 ) BOH        (REAL*8 ) : Array for OH values [molec/cm3]
!  (4 ) COPROD     (REAL*8 ) : Array for zonal mean P(CO) [v/v/s]
!  (5 ) PAVG       (REAL*8 ) : Array for 24-h avg surface pressure [mb]
!  (6 ) TAVG       (REAL*8 ) : Array for 24-h avg temperature [K]
!  (7 ) TCH4       (REAL*8 ) : Array for CH4 budget (N_CH4 items)
!  (8 ) NCMSALTS   (INTEGER) : # of altitudes for CMS climatological OH
!  (9 ) NCMSLATS   (INTEGER) : # of latitudes for CMS climatological OH
!  (10) CMSALTS    (REAL*8 ) : Altitude values for CMS climatological OH
!  (11) CMSLATS    (REAL*8 ) : Latitude values for CMS climatological OH
!  (12) AVGOH      (REAL*8 ) : Array for CMS climatological OH [molec/cm3]
!  (13) FMOL_CH4   (REAL*8 ) : Molecular weight of CH4 [kg/mole]
!  (14) XNUMOL_CH4 (REAL*8 ) : Molecules CH4 / kg CH4
!  (15) CH4_EMIS   (REAL*8 ) : Array for CH4 Emissions
!
!  Module Routines:
!  ===========================================================================
!  (1 ) GET_GLOBAL_CH4        : Computes latitudinal, yearly CH4 gradient
!  (2 ) CH4_AVGTP             : Computes 24-h average pressure & temperature
!  (3 ) EMISSCH4              : Handles CH4 emissions
!  (4 ) CHEMCH4               : Handles CH4 chemistry (various sinks)
!  (7 ) CH4_DECAY             : Computes decay of CH4 w/ OH in the troposphere
!  (8 ) CH4_OHSAVE            : Saves OH conc. for CH3CCl3 diagnostic
!  (9 ) CH4_STRAT             : Computes loss of CH4 in the stratosphere
!  (10) CH4_BUDGET            : Computes global CH4 budgets, sources & sinks
!  (11) SUM_CH4               : Sums a sub-region of the TCH4 budget array
!  (12) INIT_CH4              : Allocates and zeroes module arrays
!  (13) CLEANUP_CH4           : Deallocates module arrays
!  (14) WETLAND_EMIS          : Computes CH4 emissions from Wetland
!  (16) CH4_DISTRIB           : Distributes chemical loss of CH4 to all tracers
!  (17) BIOBURN_EMIS          : Gets CH4 emissions from GFED2 biomass burning
!  (18) RICE_EMIS             : Gets and scales CH4 rice emissions
!  (19) BIOFUEL_EMIS          : Gets CH4 emissions from Yevich and Logan 2003
!  (20) ASEASONAL_ANTHRO_EMIS : Gets aseasonal anthropogenic CH4 emissions
!  (21) ASEASONAL_NATURAL_EMIS: Gets aseasonal natural CH4 emissions
!
!  GEOS-CHEM modules referenced by global_ch4_mod.f
!  ============================================================================
!  (1 ) bpch2_mod.f      : Module containing routines for binary punch file I/O
!  (2 ) diag_mod.f       : Module containing GEOS-CHEM diagnostic arrays
!  (3 ) dao_mod.f        : Module containing arrays for DAO met fields
!  (4 ) diag_oh_mod.f    : Module containing arrays for mean OH & CH3CCl3 life
!  (4 ) error_mod.f      : Module containing NaN and other error check routines
!  (5 ) grid_mod.f       : Module containing horizontal grid information
!  (6 ) pressure_mod.f   : Module containing routines to compute P(I,J,L)
!  (7 ) time_mod.f       : Module containing routines to compute date & time
!  (8 ) tracer_mod.f     : Module containing information on tracers and
!                          concentration array.
!  (9 ) logical_mod.f    : Module containing logical variables.
!  (10) directory_mod.f  : Module containing directory variables.
!  (11) file_mod.f       : Module containing file unit numbers.
!  (12) transfer_mod.f   : Module containing routines to change type of array
!                          data.
!  (13) regrid_1x1_mod.f : Module containing regridding routines
!  (14) diag_pl_mod.f    : Module containing variables and routines for prod
!                          and loss of chemical families.
!  (15) diag_oh_mod.f    : Module containing routines to archive OH mass.
!
!  NOTES:
!  (1 ) Merged routines from jsw's CH4 code  into "global_ch4_mod.f"
!        (bmy, 1/16/01)
!  (2 ) XNUMOL_CH4 and TCH4 have to be public - all other variables can
!        be made private, so as not to conflict with other common-block
!        definitions (bmy, 1/17/01)
!  (3 ) Minor fixes from jsw added (jsw, bmy, 2/17/01)
!  (4 ) Removed some F90 module references from EMISSCH4 (bmy, 3/20/01)
!  (5 ) Eliminate obsolete commented-out code (bmy, 4/20/01)
!  (6 ) Updated comments (bmy, 9/4/01)
!  (7 ) Fixes for binary punch file in READ_COPROD (bmy, 9/26/01)
!  (8 ) Removed obsolete code from READ_COPROD (bmy, 10/24/01)
!  (9 ) Minor bug fixes for compilation on ALPHA (bmy, 11/15/01)
!  (10) Eliminate obsolete code from 11/01 (bmy, 2/27/02)
!  (11) Now eliminate PS from the arg list to CH4_AVGTP (4/11/02)
!  (12) Now divide module header into MODULE PRIVATE, MODULE VARIABLES, and
!        MODULE ROUTINES sections.  Updated comments (bmy, 5/28/02)
!  (13) Replaced all instances of IM with IIPAR and JM with JJPAR, in order
!        to prevent namespace confusion for the new TPCORE (bmy, 6/25/02)
!  (14) Now reference "file_mod.f".  Also removed obsolete code. (bmy, 6/27/02)
!  (15) Now references "pressure_mod.f" (bmy, 8/21/02)
!  (16) Now reference AD and T from "dao_mod.f".  Now reference "error_mod.f".
!        Remove obsolete code from various routines.  Remove reference to
!        header file "comtrid.h" -- it's not used. (bmy, 11/6/02)
!  (17) Minor bug fix in FORMAT statements (bmy, 3/23/03)
!  (18) Now references "grid_mod.f" and "time_mod.f" (bmy, 3/27/03)
!  (19) Updates to GET_GLOBAL_CH4 (bmy, 7/1/03)
!  (20) Now references "directory_mod.f", "tracer_mod.f", and "diag_oh_mod.f"
!        (bmy, 7/20/04)
!  (21) Now can read data for both GEOS and GCAP grids (bmy, 8/16/05)
!  (22) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (23) Updated CH4 simulation (wecht, cph, ccarouge, 10/1/09)
!       completely replace this file (kjw, dkh, 02/12/12, adj32_023)
!******************************************************************************
!
      IMPLICIT NONE

      PUBLIC

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables
      ! and routines from being seen outside "global_ch4_mod.f"
      !=================================================================

      ! PRIVATE module variables
      PRIVATE :: BAIRDENS,  BOH
      PRIVATE :: COPROD,   PAVG,      TAVG
      PRIVATE :: NSEAS,    NCMSALTS,  NCMSLATS
      PRIVATE :: CMSALTS,  CMSLATS,   AVGOH
      PRIVATE :: FMOL_CH4, CH4_EMIS

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================

      ! Number of CH4 budget types
      INTEGER, PARAMETER   :: N_CH4 = 12

      ! Various arrays
      REAL*8,  ALLOCATABLE :: BAIRDENS(:,:,:)
      REAL*8,  ALLOCATABLE :: BOH(:,:,:,:)
      REAL*8,  ALLOCATABLE :: CH4LOSS(:,:,:,:)
      REAL*8,  ALLOCATABLE :: COPROD(:,:,:)
      REAL*8,  ALLOCATABLE :: PAVG(:,:,:)
      REAL*8,  ALLOCATABLE :: TAVG(:,:,:)
      REAL*8,  ALLOCATABLE :: TCH4(:,:,:,:)

      ! For Clarisa's Climatological OH
      INTEGER, PARAMETER   :: NSEAS    = 4
      INTEGER, PARAMETER   :: NCMSALTS = 7
      INTEGER, PARAMETER   :: NCMSLATS = 24

      REAL*8               :: CMSALTS(NCMSALTS) =
     &    (/ 1000d0, 900d0, 800d0, 700d0, 500d0, 300d0, 200d0 /)

      REAL*8               :: CMSLATS(NCMSLATS) =
     &    (/ 90d0,  84d0,  76d0,  68d0,  60d0,  52d0,  44d0,  36d0,
     &       28d0,  20d0,  12d0,   4d0,  -4d0, -12d0, -20d0, -28d0,
     &      -36d0, -44d0, -52d0, -60d0, -68d0, -76d0, -84d0, -90d0 /)

      REAL*8,  ALLOCATABLE  :: AVGOH(:,:,:)

      ! FMOL_CH4   =  kg CH4    / mole CH4
      ! XNUMOL_CH4 =  molec CH4 / kg CH4
      REAL*8, PARAMETER    :: FMOL_CH4   = 16d-3
      REAL*8, PARAMETER    :: XNUMOL_CH4 = 6.0221d23 / 16d-3

      REAL*8,  ALLOCATABLE  :: CH4_EMIS(:,:,:)

      REAL*8           :: TROPOCH4
      REAL*8           :: STRATOCH4

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE CH4_AVGTP
!
!******************************************************************************
!  Subroutine CH4_AVGTP gets the 24-h average surface pressure and temperature
!  needed for the CH4 simulation. (jsw, bnd, bmy, 1/16/01, 7/20/04)
!
!  NOTES:
!  (1 ) Created by Bryan Duncan (1/99).  Adapted for CH4 chemistry and
!        placed into module "global_ch4_mod.f" by Bob Yantosca. (bmy, 1/16/01)
!  (2 ) CH4_AVGTP is independent of "CMN_OH", "CMN_CO", and "CMN_CO_BUDGET".
!        (bmy, 1/16/01)
!  (3 ) Removed duplicate definition for NTDT, NMIN (bmy, 11/15/01)
!  (4 ) Removed PS from argument list.  Now use P(I,J)+PTOP instead of
!        PS, this ensures that we have consistency between P and AD.
!        (bmy, 4/11/02)
!  (5 ) Removed obsolete code (bmy, 6/27/02)
!  (6 ) Now uses GET_PCENTER from "pressure_mod.f" to return the pressure
!        at the midpoint of the box (I,J,L).  Also added parallel DO-loops.
!        Updated comments. (dsa, bdf, bmy, 8/21/02)
!  (7 ) Now reference T from "dao_mod.f".  Now reference GEOS_CHEM_STOP from
!        "error_mod.f" (bmy, 10/15/02)
!  (8 ) Removed NTDT, NMIN from the arg list.  Now uses functions GET_TS_DYN,
!        GET_TS_CHEM, and GET_ELAPSED_MIN from "time_mod.f" (bmy, 3/27/03)
!  (9 ) Remove reference to CMN, it's not needed (bmy, 7/20/04)
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,      ONLY : T
      USE ERROR_MOD,    ONLY : GEOS_CHEM_STOP
      USE PRESSURE_MOD, ONLY : GET_PCENTER
      USE TIME_MOD,     ONLY : GET_TS_DYN, GET_TS_CHEM, GET_ELAPSED_MIN
      USE TIME_MOD,     ONLY : GET_NYMD, GET_NHMS, GET_NYMDe

#     include "CMN_SIZE"  ! Size parameters

      ! Local variables
      INTEGER             :: NTDT, NMIN
      INTEGER             :: I, J, L, NTIMES, MNDT, K, M, N
      INTEGER, SAVE       :: NNCOUNT
      INTEGER, SAVE       :: NNEW
      REAL*8              :: Ptemp(IIPAR,JJPAR,LLPAR)

      !=================================================================
      ! CH4_AVGTP begins here!
      !=================================================================

      ! Get quantities from "time_mod.f"
      NTDT   = GET_TS_DYN() * 60
      NMIN   = GET_ELAPSED_MIN()
      MNDT   = NTDT  / 60
      NTIMES = GET_TS_CHEM() / MNDT

      ! NTIMES is the number of dynamic timesteps in a chem timestep
      ! If we're in the first day, NTIMES = NTIMES + 1

      IF ( NMIN .LE. GET_TS_CHEM() ) NTIMES = NTIMES + 1
      ! If we're in the final day, NTIMES = NTIMES - 1
      ! It doesn't really matter because CH4 chem is not done at the end of the final day.
      ! IF ( GET_NYMD()+1 .eq. GET_NYMDe() ) NTIMES = NTIMES - 1


      ! At the start of the run...
      IF ( NMIN == 0 ) THEN

         ! Initialize NNEW
	 NNEW = 0

         ! Error check -- need chem timestep (1440) to be divisible by
         ! dyn timestep
         IF ( mod( GET_TS_CHEM(), MNDT ) /= 0 ) THEN
            WRITE(*,*) ' '
            WRITE(*,*) 'CH4-OH parameterization option (i.e., NSRCX=5)!'
            WRITE(*,*) 'The chemistry time step (i.e., 24 hours) is'
            WRITE(*,*) 'not evenly divisible by the meteorological'
            WRITE(*,*) 'data read-in time step (i.e., 6 hours).  This'
            WRITE(*,*) 'will mess up SR avgtp which calculates a 24-'
            WRITE(*,*) 'hour average temperature and pressure to be'
            WRITE(*,*) 'used by SR getinfo.'
            WRITE(*,*) ' '
            CALL GEOS_CHEM_STOP
         ENDIF

         ! If NCHEM < NTDT then stop program.
         IF ( GET_TS_CHEM() < MNDT ) THEN
            WRITE(*,*) ' '
            WRITE(*,*) 'When using the CH4-OH parameterization'
            WRITE(*,*) 'option (i.e., NSRCX=5), take a 24-hour'
            WRITE(*,*) 'time step (i.e., NCHEM=1440 min.) because'
            WRITE(*,*) 'the OH parameterization produces a 24-hour'
            WRITE(*,*) 'average [OH]'
            WRITE(*,*) ' '
            CALL GEOS_CHEM_STOP
         ENDIF
      ENDIF

      !=================================================================
      ! If a new 24-hr period, set Pavg = 0, and reset NNEW, NCOUNT
      !=================================================================
      IF ( NNEW == 0 ) THEN
         Pavg(:,:,:) = 0d0
         Tavg(:,:,:) = 0d0
	 NNEW        = 1
	 NNCOUNT     = 0
      ENDIF

      !=================================================================
      ! Archive quantities
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, PTEMP )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Archive pressure
         Pavg(I,J,L) = Pavg(I,J,L) + GET_PCENTER(I,J,L)

         ! Archive temperature
         Tavg(I,J,L) = Tavg(I,J,L) + T(I,J,L)
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !================================================================
      ! Keep track to see if at end of NCHEM time step.
      ! If so, divide PAVG & TAVG by the number of times archived.
      !=================================================================

      NNCOUNT = NNCOUNT + 1

      IF ( NNCOUNT == NTIMES ) THEN
         Pavg(:,:,1:LLPAR) = Pavg(:,:,1:LLPAR) / DBLE( NTIMES )
         Tavg(:,:,1:LLPAR) = Tavg(:,:,1:LLPAR) / DBLE( NTIMES )
         NNEW              = 0
      ENDIF

      ! Return to calling program
      END SUBROUTINE CH4_AVGTP

!------------------------------------------------------------------------------

      SUBROUTINE EMISSCH4
!
!******************************************************************************
!  Subroutine EMISSCH4 places emissions of CH4 [kg] into the STT array.
!  (jsw, bnd, bey, bmy, 1/16/01, 10/3/05)
!
!  WARNING: Soil absorption has to be the 11th field in CH4_EMIS
!
!  NOTES:
!  (1 ) Created by Bryan Duncan (1/99).  Adapted for CH4 chemistry by
!        James Wang (7/00).  Inserted into module "global_ch4_mod.f"
!        by Bob Yantosca. (bmy, 1/16/01)
!  (2 ) EMISSCH4 is independent of "CMN_OH", "CMN_CO", and "CMN_CO_BUDGET".
!        (bmy, 1/16/01)
!  (3 ) GLOBASEAEMIS, GLOBSEAEMIS are diagnostics by jsw.
!  (4 ) Do not multiply CO emissions by 1.28 anymore (jsw, bmy, 2/12/01)
!  (5 ) Renamed input files to CH4_monthly.geos.{RES} and
!        CH4_aseasonal.geos.{RES}. (bmy, 2/12/01)
!  (6 ) Add reference to "CMN_SETUP" for the DATA_DIR variable (bmy, 2/13/01)
!  (7 ) Removed references to "biofuel_mod.f" and "biomass_mod.f"; these
!        weren't necessary (bmy, 3/20/01)
!  (8 ) Now reference IU_FILE and IOERROR from "file_mod.f".  Now use IU_FILE
!        instead of IUNIT as the file unit #. (bmy, 6/27/02)
!  (9 ) Now reference BXHEIGHT and SUNCOS from "dao_mod.f".  Remove reference
!        to header file "comtrid.h" -- it's not used.  Make FIRSTEMISS a local
!        SAVEd variable.  Also use MONTH from "CMN" instead of the variable
!        LMN. (bmy, 11/15/02)
!  (10) Now replace DXYP(JREF)*1d4 with routine GET_AREA_CM2 of "grid_mod.f".
!        Now use function GET_MONTH and GET_TS_EMIS from "time_mod.f".
!        Now use functions GET_XOFFSET and GET_YOFFSET from "grid_mod.f".
!        I0 and J0 are now local variables. (bmy, 3/27/03)
!  (11) Now reference STT from "tracer_mod.f".  Now reference DATA_DIR from
!        "directory_mod.f". (bmy, 7/20/04)
!  (12) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (13) Add non-local PBL capability (ccc, 8/31/09)
!  (14) Add adjoint scaling to all emis subroutines called (kjw, 2/22/10)
!******************************************************************************
!
      ! References to F90 modules
      USE TIME_MOD,      ONLY : ITS_A_NEW_MONTH, ITS_A_NEW_YEAR
      USE TIME_MOD,      ONLY : GET_MONTH,       GET_YEAR
      USE TIME_MOD,      ONLY : GET_TS_EMIS,     ITS_A_NEW_DAY
      USE GRID_MOD,      ONLY : GET_AREA_CM2,    GET_XOFFSET
      USE GRID_MOD,      ONLY : GET_YOFFSET
      USE TRACER_MOD,    ONLY : STT
      USE LOGICAL_MOD,   ONLY : LSPLIT
      USE LOGICAL_MOD,   ONLY : LGFED3BB, LDAYBB3
      USE DIAG_MOD,      ONLY : AD58
      USE ERROR_MOD,     ONLY : GEOS_CHEM_STOP,  IT_IS_NAN
      USE TRACER_MOD,    ONLY : N_TRACERS, ID_TRACER
!      USE LOGICAL_MOD,   ONLY : LWETL,           LBMCH4,       LRICE
!      USE LOGICAL_MOD,   ONLY : LBFCH4
      USE ADJ_ARRAYS_MOD,  ONLY : ADCH4EMS, EMS_SF, N_CALC
      USE LOGICAL_ADJ_MOD, ONLY : LADJ,     LADJ_EMS
      USE TIME_MOD,        ONLY : GET_NYMD, GET_NHMS
      USE LOGICAL_ADJ_MOD, ONLY : LFD_GLOB
      !kjw iterative optimization
      !USE LOGICAL_ADJ_MOD, ONLY : LCH4_ITERATE


!     Non-Local PBL mixing scheme is not available in adjoint (kjw, 2/20/2010)
!      USE VDIFF_PRE_MOD, ONLY : EMIS_SAVE ! (ccc, 08/31/09)
!      USE LOGICAL_MOD,   ONLY : LNLPBL    ! (ccc, 08/31/09)


#     include "CMN_SIZE"     ! Size parameters
#     include "CMN_DIAG"     ! Diagnostic switches

      ! Local Variables
      INTEGER                :: I, I0, IREF, J, J0, JREF, N, M
      REAL*8                 :: DTSRCE,       AREA_CM2

      ! Get nested-grid offsets
      I0  = GET_XOFFSET()
      J0  = GET_YOFFSET()


      !===================================================================
      ! Emissions are read or calculated at the first of every:
      !   1) Emission time step - Natural Wetlands (from J Kaplan)
      !   2) Month              - Biomass Burning and Rice
      !   3) Year               - All other sources
      !
      ! Emissions are stored at each time step in a 3D array:
      !   EMIS_CH4(IIPAR,JJPAR,N).
      !     Where N = 1:12
      !       1. Total Emissions (including soil absorption, counted neg.)
      !       2. Oil and Gas Processing
      !       3. Coal Mining
      !       4. Livestock
      !       5. Waste
      !       6. Biofuel
      !       7. Rice
      !       8. Other Anthropogenic
      !       9. Biomass Burning
      !       10. Wetlands
      !       11. Soil Absorption
      !       12. Other Natural
      !
      ! Emissions are added to STT array and AD58 (emission diagnostic)
      !   at every time step.
      !                                                      (kjw, 6/4/09)
      !===================================================================


      print*,'% --- ENTERING EMISSCH4!'


      ! ==================================================================
      ! 1)  Get Wetland Emissions
      !     NOTES:                                          (kjw, 5/28/09)
      !        Emissions calculated online every timestep in WETLAND_EMIS.
      !        WETLAND_EMIS adapted to GEOS-Chem by Jerome Drevet (3/06)
      !        from a wetland methane scheme provided by Jed O. Kaplan
      !        See subroutine WETLAND_EMIS for more information
      ! ==================================================================


      !4.1 Wetland emissions
      CALL WETLAND_EMIS( CH4_EMIS )

      ! ==================================================================
      ! 2)  Get Daily Varying CH4 Emissions
      !     NOTES:                                          (kjw, 9/24/12)
      !        Biomass burning emissions from GFED3
      !        Can be daily or monthly emissions
      ! ==================================================================
      IF ( ( ITS_A_NEW_DAY()   .AND. LDAYBB3  )   .OR.
     &     ( ITS_A_NEW_MONTH() .AND. LGFED3BB ) ) THEN

         !4.2 Biomass Burning emissions (CH4_BBN, #9)
         CALL BIOBURN_EMIS( CH4_EMIS )

      ENDIF

      ! ==================================================================
      ! 2)  Get Monthly Varying CH4 Emissions
      !     NOTES:                                          (kjw, 5/28/09)
      !        Biomass burning emissions from GFED2 or GFED3
      !        Biomass burning available from 1997-2007 (5/28/09)
      !        Rice emissions from EDGAR v4, modified by GEOS soil wetness
      ! ==================================================================
      IF ( ITS_A_NEW_MONTH() ) THEN

         !4.3 Rice emissions (CH4_RIC, #7)
         CALL RICE_EMIS( CH4_EMIS )

      ENDIF


      ! ==================================================================
      ! 3)  Get Aseasonal CH4 Emissions
      !     NOTES:                                      (kjw, 5/28/09)
      !        Anthropogenic emissions from EDGAR v4 except biofuel
      !        emissions which are from Yevich and Logan 2003.
      !        Soil absorption from Fung et. al. 1991.
      !        Other natural emissions include:
      !            termites from Fung et. al. 1991
      ! ==================================================================

      IF ( ITS_A_NEW_YEAR() ) THEN

        !4.4 Biofuel emissions (CH4_BFL, #6)
        !kjw replace with EDGARv4 biofuels in ASEASONAL_ANTHRO_EMIS
        !    (kjw, 11/17/11)
        !CALL BIOFUEL_EMIS( CH4_EMIS )

        !4.5 Aseasonal Anthropogenic emissions
        ! (CH4_OAG, #2; CH4_COL, #3; CH4_LIV, #4; CH4_WST, #5; CH4_OTA, #8)
        CALL ASEASONAL_ANTHRO_EMIS( CH4_EMIS )

        !4.6 Aseasonal Natural emissions (CH4_SAB, #11; CH4_OTN, #12)
        CALL ASEASONAL_NATURAL_EMIS( CH4_EMIS )

      ENDIF


      ! Total emission: sum of all emissions - (2*soil absorption)
      ! We have to substract soil absorption twice because it is added
      ! to other emissions in the SUM function. (ccc, 7/23/09)
      CH4_EMIS(:,:,1) = 0d0
      CH4_EMIS(:,:,1) = SUM(CH4_EMIS, 3) - (2 * CH4_EMIS(:,:,11))

      ! However, since we don't optimize for soil absorption, don't
      ! take it out of total emissions yet. Wait until after EMS_SF scaling
      ! to remove soil absorption. So, add soil absorption back to give
      ! total emissions (excluding soil absorption)
      CH4_EMIS(:,:,1) = CH4_EMIS(:,:,1) + CH4_EMIS(:,:,11)


      ! =================================================================
      ! Do Adjoint Scaling
      ! Remove LCH4_ITERATE for compatibility with v35c
      !IF ( LADJ .AND. (LFD_GLOB .OR. LADJ_EMS .OR. LCH4_ITERATE) ) THEN
      IF ( LADJ .AND. (LFD_GLOB .OR. LADJ_EMS) ) THEN

         ! Determine scale group (temporal)
         M = GET_SCALE_GROUP()
         WRITE(6,*) ' % - Rescale emissions: use SCALE_GROUP ',  M

         ! Rescale emissions
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            CH4_EMIS(I,J,1) = CH4_EMIS(I,J,1) * EMS_SF(I,J,M,ADCH4EMS)
         ENDDO
         ENDDO

      ENDIF
      ! =================================================================

      ! Now that we've done the emission factor scaling, add soil absorption
      ! back to the total emissions array
      CH4_EMIS(:,:,1) = CH4_EMIS(:,:,1) - CH4_EMIS(:,:,11)


      ! =================================================================
      ! Modify the STT with emissions rates.               (kjw, 5/29/09)
      ! There are 12 tracers in the multi-tracer run.
      ! One tracer for total CH4 and one for each emission category.
      !
      !  1. Total CH4
      !  2. Gas and Oil
      !  3. Coal
      !  4. Livestock
      !  5. Waste
      !  6. Biofuel
      !  7. Rice
      !  8. Other Anthropogenic
      !  9. Biomass Burning
      !  10. Wetlands
      !  11. Soil Absorption
      !  12. Other Natural
      ! =================================================================

      WRITE( 6, '(a)' ) '% EMISSCH4 --- Adding Emissions to STT array.'


      ! Emission timestep [s]
      DTSRCE = GET_TS_EMIS() * 60d0   !timestep in s.

      ! J0 and I0 are global variables, both set = 0.
      DO J = 1, JJPAR
         JREF = J + J0

         ! Grid box surface area [cm2]
         AREA_CM2 = GET_AREA_CM2( J )

      DO I = 1, IIPAR
         IREF = I + I0

         DO N = 1, N_TRACERS
               STT(IREF,JREF,1,N) = STT(IREF,JREF,1,N) +
     &            CH4_EMIS(IREF,JREF,ID_TRACER(N))
     &            / XNUMOL_CH4 * DTSRCE * AREA_CM2
         ENDDO


         IF ( ND58 > 0 ) THEN
            ! All emission sources except soil absorption
            AD58(IREF,JREF,1) = AD58(IREF,JREF,1) +
     &         ( CH4_EMIS(IREF,JREF,1) + CH4_EMIS(IREF,JREF,11) )
     &         / XNUMOL_CH4 * DTSRCE * AREA_CM2

            DO N = 2, PD58
               AD58(IREF,JREF,N) = AD58(IREF,JREF,N) +
     &            CH4_EMIS(IREF,JREF,N)
     &            / XNUMOL_CH4 * DTSRCE * AREA_CM2
            ENDDO

         ENDIF

      ENDDO
      ENDDO


      !===============================================================
      ! Sum up CH4 budgets
      !
      ! TCH4 - # molecules emitted from different sources
      !===============================================================


      DO J = 1,JJPAR
         JREF = J + J0

         ! Grid box surface area [cm2]
         AREA_CM2 = GET_AREA_CM2( J )

      DO I = 1,IIPAR
         IREF = I + I0
         ! Gas, oil, mine
         TCH4(I,J,1,5) = TCH4(I,J,1,5) +
     &        ( ( CH4_EMIS(I,J,2) + CH4_EMIS(I,J,3) ) *
     &                          AREA_CM2 * DTSRCE )

         ! agriculture (rice, animals, waste)
         TCH4(I,J,1,6) = TCH4(I,J,1,6) +
     &        ( ( CH4_EMIS(I,J,4) + CH4_EMIS(I,J,7) +
     &            CH4_EMIS(I,J,5) ) * AREA_CM2 * DTSRCE )

         ! Biomass burning (and biofuel?)
         TCH4(I,J,1,7) = TCH4(I,J,1,7)+
     &        ( ( CH4_EMIS(I,J,9) + CH4_EMIS(I,J,6) ) *
     &                          AREA_CM2 * DTSRCE )

         ! Termites
         TCH4(I,J,1,8) = TCH4(I,J,1,8)+
     &        ( CH4_EMIS(I,J,12) * AREA_CM2 * DTSRCE )

         ! Wetlands
         TCH4(I,J,1,9) = TCH4(I,J,1,9)+
     &        ( CH4_EMIS(I,J,10) * AREA_CM2 * DTSRCE )

         ! Soil Absorption
         TCH4(I,J,1,10) = TCH4(I,J,1,10)+
     &        ( CH4_EMIS(I,J,11) * AREA_CM2 * DTSRCE )


         TCH4(I,J,1,4) = TCH4(I,J,1,4) +
     &       ( CH4_EMIS(I,J,1)  + ( 2 * CH4_EMIS(IREF,JREF,11) ))
     &       * AREA_CM2 * DTSRCE

      ENDDO
      ENDDO

      ! Return to calling program
      END SUBROUTINE EMISSCH4

!------------------------------------------------------------------------------

      SUBROUTINE WETLAND_EMIS( EMIS3D )
!
!******************************************************************************
!  Subroutine WETLAND_CH4 calculates emissions of CH4 [kg] by Wetland.
!
!  NOTES:
!  (1 ) Adapted by Jérôme Drevet (3/06) from the BIOME-TG Wetland-Methane
!       scheme provided by Jed O. Kaplan.
!  (2 ) CH4 Emissions from Wetland depend on:
!		a - Soil Carbon content.
!		b - Vegetation type
!		c - Wetland area (%)
!		d - Soil moisture.
!       a, b, c are taken from the LPJ, a vegetation model. Data are provided
!	by J.O.Kaplan. Soil moisture is read from GEOS Met input files.
!
!  (3 ) Corrected order of DO loops (bmy, 10/1/09)
!  (4 ) Add adjoint scaling to emissions
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,       ONLY : GWETTOP,      LWI
      USE DAO_MOD,       ONLY : TSKIN,        TS
      USE DAO_MOD,       ONLY : FRLAND,       FRLAKE
      USE DAO_MOD,       ONLY : FROCEAN,      FRLANDIC
      USE BPCH2_MOD,     ONLY : GET_RES_EXT,  GET_MODELNAME
      USE BPCH2_MOD,     ONLY : GET_TAU0,     READ_BPCH2
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT, GET_NAME_EXT_2D
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE FILE_MOD,      ONLY : IU_FILE,      IOERROR
      USE GRID_MOD,      ONLY : GET_AREA_M2
      USE TIME_MOD,      ONLY : GET_MONTH,    GET_YEAR,    GET_TS_EMIS
      USE TIME_MOD,      ONLY : ITS_A_NEW_MONTH, ITS_A_NEW_YEAR
      USE TRANSFER_MOD,  ONLY : TRANSFER_2D
      USE DIAG_MOD,      ONLY : AD60, AD58
      USE TIME_MOD,      ONLY : GET_DIRECTION

#     include "CMN_SIZE"     ! Size parameters
#     include "CMN_DIAG"     ! Diagnostic switches
#     include "CMN"          ! PD58

      ! Arguments
      REAL*8, INTENT(INOUT)     :: EMIS3D(IIPAR,JJPAR,PD58)

      ! Local Variables
      INTEGER :: I, J, L
      INTEGER :: GM, YEAR
      REAL*4  :: ARRAY(IIPAR,JJPAR)
      REAL*8  :: WETFRAC(IIPAR,JJPAR)
      REAL*8  :: REALWET(IIPAR,JJPAR)
      REAL*8  :: EFF_GWET(IIPAR,JJPAR)
      REAL*8  :: SOIL_C(IIPAR,JJPAR)
      REAL*8  :: LITTER_C(IIPAR,JJPAR)
      REAL*8  :: litterfast
      REAL*8  :: litterslow
      REAL*8  :: soilfast
      REAL*8  :: soilslow
      REAL*8  :: HETEROR
      REAL*8  :: F_TEMP
      REAL*8  :: MEAN_T(IIPAR,JJPAR)
      REAL*8  :: METHANE_OUT(IIPAR,JJPAR)
      REAL*8  :: XTAU
      REAL*8  :: TROPICNESS
      REAL*8  :: EMIT_TROPIC
      REAL*8  :: EMIT_TEMPER
      REAL*8  :: MOIST_SCALE
      REAL*8  :: EMIT_FACT
      INTEGER :: MONTHDATES(12) = (/ 31, 28, 31, 30,
     &                               31, 30, 31, 31,
     &                               30, 31, 30, 31 /)
      CHARACTER(LEN=255) :: FILENAME
      CHARACTER(LEN=4)   :: CYEAR

      !=================================================================
      ! WETLAND_CH4 begins here!
      !=================================================================

      !4.10 Wetland emissions

      !===================================================================
      ! Get wetland fraction data
      !===================================================================

      WRITE( 6, '(a)' ) '% EMISSCH4 --- BEGIN WETLAND_EMIS'



#if   defined( GRID05x0666 ) && defined( NESTED_NA )
         FILENAME = TRIM( DATA_DIR )  // 'CH4_201203/wetlands/' //
     &              'Wetfrac.'        // GET_RES_EXT()          //
     &              '_NA.bpch'
#else
         FILENAME = TRIM( DATA_DIR )  // 'CH4_201203/wetlands/' //
     &              'Wetfrac.'        // GET_RES_EXT()          // 
     &              '.bpch'
#endif

      WRITE( 6, 91 ) TRIM ( FILENAME )
 91   FORMAT( '     - WL_CH4: Reading WET-FRAC: ', a )
      CALL FLUSH( 6 )

      XTAU = GET_TAU0( 1, 1, 2000 )

      CALL READ_BPCH2( FILENAME, 'WET-FRAC',     1,
     &                     XTAU,      IIPAR,     JJPAR,
     &                        1,      ARRAY,     QUIET=.TRUE.)

      CALL TRANSFER_2D( ARRAY(:,:), WETFRAC(:,:) )

      ! WETFRAC is maximum inundatable area in a box
      WETFRAC =  WETFRAC / 100d0


      !===================================================================
      ! Calculate inundated fraction
      !
      ! REALWET calculation is based on maximum inundatable area (WETFRAC)
      ! and top soil moisture information
      !
      ! NOTE: LWI (land/water/ice flag) definition has changed between
      !   GEOS4 and GEOS5.  This contributes to the variance between GEOS4
      !   and GEOS5 wetland emissions.  Below is Jerome Drevet's and Jed
      !   Kaplan's original calculation of REALWET using GEOS4 and a
      !   modified calculation using GEOS5.
      !                                                     (kjw, 6/10/09)
      !===================================================================

      ! REALWET is the actual inundated fraction of a box
      REALWET(:,:) = 0d0

      ! GEOS4 calculation
#if   defined( GEOS_4 )
      DO J = 1, JJPAR
      DO I = 1, IIPAR
	 ! We don't want emissions in frozen regions
	 IF (TSKIN(I,J) > 273) THEN
         ! We want emissions from land boxes only
         IF (LWI(I,J) == 1) THEN
	    ! If wetness>0.1, the wetland fraction is equal
            ! to the maximal potential wetland fraction
            IF (GWETTOP(I,J) > 0.1) THEN
		REALWET(I,J) = WETFRAC(I,J)
            ELSE
		REALWET(I,J) = 0.
            ENDIF
	 ENDIF
	 ENDIF
      ENDDO
      ENDDO

      ! GEOS5 Calculation
#elif defined( GEOS_5 ) || defined( GEOS_FP ) || defined( MERRA2 )

      DO J = 1, JJPAR
      DO I = 1, IIPAR
	 ! We don't want emissions in frozen regions
	 IF (TSKIN(I,J) > 273) THEN
         ! We want emissions from any box that contains some land
         ! FRLAND is fraction of grid box that is land
         IF (FRLAND(I,J) > 0) THEN
            ! Actual wetness of land /= GWETTOP because GWETTOP includes
            ! wetness in lakes, ocean, and ice.  Below is a scheme to
            ! calculate effective GWETTOP of the land fraction
            EFF_GWET(I,J) = ( GWETTOP(I,J) -
     &           ( FROCEAN(I,J) + FRLAKE(I,J) + FRLANDIC(I,J) ) )
     &                       / FRLAND(I,J)

            ! Catch for negative EFF_GWET
            IF ( EFF_GWET(I,J) < 0 ) THEN
               EFF_GWET(I,J) = 0d0
            ENDIF

	    ! If wetness>0.1, the wetland fraction is equal
            ! to the maximal potential wetland fraction
            IF (EFF_GWET(I,J) > 0.1) THEN
		REALWET(I,J) = WETFRAC(I,J)
            ELSE
		REALWET(I,J) = 0.
            ENDIF
	 ENDIF
	 ENDIF
      ENDDO
      ENDDO


#endif

         GM = GET_MONTH()
!kjw_adjoint
      ! Update Wetland Fraction Diagnostic if in the forward simulation
      IF ( GET_DIRECTION() .EQ. 1 ) THEN
         IF ( ND60 > 0 ) THEN
            AD60(:,:) = AD60(:,:) + REALWET(:,:)/(24d0*MONTHDATES(GM))
         ENDIF
      ENDIF
!kjw_adjoint

      !===================================================================
      ! Get litter carbon and soil carbon from LPJ DGVM (in gC/m2).
      !===================================================================

      ! Carbon litter and carbon soil files both have tau date of
      ! Jan. 1, 2000

      XTAU = GET_TAU0( 1, 1, 2000 )

#if defined( GRID05x0666 ) && defined( NESTED_NA )
      FILENAME = TRIM( DATA_DIR ) // 'CH4_201203/wetlands/' //
     &           'Carbon_litter.'     // GET_RES_EXT() //
     &           '_NA.bpch'
#else
      FILENAME = TRIM( DATA_DIR ) // 'CH4_201203/wetlands/' //
     &           'Carbon_litter.'     // GET_RES_EXT() //
     &           '.bpch'
#endif

      CALL READ_BPCH2( FILENAME, 'CO--SRCE',   1,
     &                        XTAU,      IIPAR,   JJPAR,
     &                           1,      ARRAY,   QUIET=.TRUE.)

      CALL TRANSFER_2D( ARRAY, LITTER_C )



#if defined( GRID05x0666 ) && defined( NESTED_NA )
      FILENAME = TRIM( DATA_DIR )  // 'CH4_201203/wetlands/' //
     &           'Carbon_soil.'    // GET_RES_EXT() //
     &           '_NA.bpch'
#else
      FILENAME = TRIM( DATA_DIR )  // 'CH4_201203/wetlands/' //
     &           'Carbon_soil.'    // GET_RES_EXT() //
     &           '.bpch'
#endif

      CALL READ_BPCH2( FILENAME, 'CO--SRCE',   1,
     &                 XTAU,      IIPAR,   JJPAR,
     &                 1,      ARRAY,   QUIET=.TRUE.)

      CALL TRANSFER_2D( ARRAY, SOIL_C )


      !===================================================================
      ! Get annual mean skin temperature.
      !===================================================================

      YEAR = GET_YEAR()

#if   defined( GEOS_5 )
         IF ( YEAR .LT. 2004 ) YEAR  = 2004
         IF ( YEAR .GT. 2010 ) YEAR  = 2010
         WRITE( CYEAR, '(i4)' ) YEAR

         FILENAME = TRIM( DATA_DIR ) // 'CH4_201305/wetlands/' //
     &              'TSKIN.'         // GET_NAME_EXT()         // 
     &              '.'              // GET_RES_EXT()          //
     &              '.'              // CYEAR                  // 
     &              '.bpch'

#elif defined( GEOS_4 )
         IF ( YEAR .LT. 2000 ) YEAR  = 2000
         IF ( YEAR .GT. 2006 ) YEAR  = 2006
         WRITE( CYEAR, '(i4)' ) YEAR

         FILENAME = TRIM( DATA_DIR )  // 'CH4_200911/wetlands/' //
     &              'TSKIN.'          // CYEAR                  // 
     &              '.'               // GET_NAME_EXT()         // 
     &              '.'               // GET_RES_EXT()
#endif

#if   defined( GRID05x0666 ) && defined( NESTED_NA )
         IF ( YEAR .LT. 2004 ) YEAR  = 2004
         IF ( YEAR .GT. 2004 ) YEAR  = 2009
         IF ( YEAR .GE. 2010 ) YEAR  = 2010
         WRITE( CYEAR, '(i4)' ) YEAR

         FILENAME = TRIM( DATA_DIR ) // 'CH4_201305/wetlands/' //
     &              'TSKIN.'         // GET_NAME_EXT()         //
     &              '.'              // GET_RES_EXT()          //
     &              '_NA.'           // CYEAR                  //
     &              '.bpch' 
#endif

      XTAU = GET_TAU0( 1, 1, YEAR )

      CALL READ_BPCH2( FILENAME,  'GMAO-2D',    2,
     &                     XTAU,      IIPAR,    JJPAR,
     &                        1,      ARRAY,    QUIET=.TRUE.)

      CALL TRANSFER_2D( ARRAY, MEAN_T )


      !===================================================================
      ! Calculate CH4 emissions!
      !===================================================================

      METHANE_OUT = 0d0

      DO J = 1, JJPAR
      DO I = 1, IIPAR

	 IF ( tskin(I,J) < 233. ) THEN
	    F_TEMP = 0
	 ELSE
	    F_TEMP = exp(308.56*(1.0/56.02-
     &           1.0/(tskin(I,J)-227.13))) !Lloyd & Taylor 1994
	 ENDIF

         ! Calculate Heterotrophic respiration
         litterfast = 0.985 * LITTER_C(i,j)
         litterslow = 0.015 * LITTER_C(i,j)
         soilfast =  0.985 * SOIL_C(i,j)
         soilslow =  0.015 * SOIL_C(i,j)

	 HETEROR = 1e3* F_TEMP *( litterfast*0.3
     &                          + litterslow*0.05
     &                          + soilfast*0.03
     &                          + soilslow*0.001 ) * 0.34 / 12.


         ! Calculate "tropicness" of each box
       	 TROPICNESS = exp((MEAN_T(I,J) - 303.15) / 8.)
	 IF ( TROPICNESS < 0 ) THEN
            TROPICNESS = 0
	 ENDIF
	 IF ( TROPICNESS > 1 ) THEN
            TROPICNESS = 1
	 ENDIF


         EMIT_TROPIC = 0.0
         EMIT_TEMPER = 0.0

         ! (moist_scale can be between 0.07 and 0.14)
         ! (emit_fact can be between 0.001 and 0.005)
         ! the lines above are comments by Jerome.  His paper publishes
         ! a value of 0.19 for MOIST_SCALE (kjw, 6/9/09)
         MOIST_SCALE = 0.205
         EMIT_FACT   = 0.018

         EMIT_TROPIC = HETEROR * MOIST_SCALE * REALWET(I,J)

         EMIT_TEMPER = HETEROR * EMIT_FACT   * REALWET(I,J)

         METHANE_OUT(I,J) = TROPICNESS      * EMIT_TROPIC +
     &                     (1 - TROPICNESS) * EMIT_TEMPER  !gCH4/m2/mth

         IF (METHANE_OUT(I,J) < 0 ) THEN
            METHANE_OUT(I,J)=0
         ENDIF

      ENDDO
      ENDDO


      ! METHANE_OUT:  g/m2/y --> molec/cm2/s
      METHANE_OUT = METHANE_OUT/16d0/1e4/(24d0*MONTHDATES(GM)*3.6e3)
     &              * 6.023e23

      EMIS3D(:,:,10) = METHANE_OUT


      ! Return to calling program
      END SUBROUTINE WETLAND_EMIS


!------------------------------------------------------------------------------

      SUBROUTINE BIOBURN_EMIS( EMIS3D )
!
!******************************************************************************
!  Subroutine BIOBURN_EMIS calculates CH4 emissions from GFED2 or GFED3 biomass
!  burning. (kjw, 6/03/09)
!
!  NOTES:
!
!******************************************************************************
!
      ! References to F90 modules
      USE BIOMASS_MOD,      ONLY : BIOMASS,       IDBCH4
      USE LOGICAL_MOD,      ONLY : LGFED2BB,      LGFED3BB

#     include "CMN_SIZE"     ! Size parameters
#     include "CMN_DIAG"     ! Diagnostic switches
#     include "CMN"          ! PD58

      ! Arguments
      REAL*8, INTENT(INOUT)     :: EMIS3D(IIPAR,JJPAR,PD58)

      ! Local Variables
      REAL*8                 :: E_CH4
      INTEGER                :: I, J

      WRITE( 6, '(a)' ) '% EMISSCH4 --- BEGIN BIOBURN_EMIS'

      !=================================================================
      ! BIOBURN_EMIS begins here!
      !=================================================================

      !4.9 Biomass Burning emissions. Calculate emissions from monthly
      !     GFED-2 or GFED-3.
      IF ( LGFED2BB .OR. LGFED3BB ) THEN
         DO I=1,IIPAR
         DO J=1,JJPAR
            ! Biomass burning emissions [molec/cm2/s]
            E_CH4 = BIOMASS( I,J,IDBCH4 )

            ! Place into CH4_EMIS array
            CH4_EMIS(I,J,9) = E_CH4

         ENDDO
         ENDDO
      ENDIF

      ! Return to calling program
      END SUBROUTINE BIOBURN_EMIS


!------------------------------------------------------------------------------

      SUBROUTINE RICE_EMIS( EMIS3D )
!
!******************************************************************************
!  Subroutine RICE_EMIS calculates CH4 emissions from rice and places CH4 [kg]
!  into the STT array. (kjw, 6/03/09)
!
!  Rice Emissions are scaled to GEOS soil wetness.  Scaling sceme developed
!     and implemented by Jerome Drevet.
!  Wetland emissions are modified by the presence of rice emissions.  Sceme
!     developed by Jerome Drevet.
!
!  NOTES:
!  (1 ) CH4 emissions from rice calculated with a routine created by Jerome
!       Drevet.  Adapted as its own subroutine by Kevin Wecht (6/03/09)
!  (2 ) Corrected ordering of DO loops (bmy, 10/1/09)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,     ONLY : GET_RES_EXT,  GET_MODELNAME
      USE BPCH2_MOD,     ONLY : GET_TAU0,     READ_BPCH2
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT, GET_NAME_EXT_2D
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE GRID_MOD,      ONLY : GET_AREA_CM2
      USE TIME_MOD,      ONLY : GET_MONTH,    GET_YEAR
      USE LOGICAL_MOD,   ONLY : LSPLIT
      USE TRANSFER_MOD,  ONLY : TRANSFER_2D

#     include "CMN_SIZE"     ! Size parameters
#     include "CMN_DIAG"     ! Diagnostic switches
#     include "CMN"          ! PD58

      ! Arguments
      REAL*8, INTENT(INOUT)     :: EMIS3D(IIPAR,JJPAR,PD58)

      ! Local Variables
      INTEGER                :: I,J,M,YEAR
      REAL*4                 :: scale
      REAL*4                 :: ARRAY(IIPAR,JJPAR)
      REAL*8                 :: DTSRCE,       AREA_CM2
      REAL*8                 :: MEAN_GWETTOP(IIPAR,JJPAR)
      REAL*8                 :: MONTH_GWETTOP(IIPAR,JJPAR)
      REAL*8                 :: wet_ratio
      REAL*8                 :: XTAU
      REAL*8                 :: RICE_OUT(IIPAR,JJPAR)
      REAL*8                 :: WETL_OUT(IIPAR,JJPAR)
      CHARACTER(LEN=255)     :: FILENAME
      CHARACTER(LEN=4)       :: CYEAR


      WRITE( 6, '(a)' ) '% EMISSCH4 --- BEGIN RICE_EMIS'

      !=================================================================
      ! RICE_EMIS begins here!
      !=================================================================

      !4.7 Rice Emissions
      ! For now, we only have emissions from 2004.

      WRITE( CYEAR, '(i4)' ) GET_YEAR()
      XTAU = GET_TAU0( 1, 1, GET_YEAR() )

      IF ( GET_YEAR() .LT. 2004 ) THEN
         CYEAR='2004'
         XTAU = GET_TAU0( 1, 1, 2004 )
      ENDIF
      IF ( GET_YEAR() .GT. 2008 ) THEN
         CYEAR='2008'
         XTAU = GET_TAU0( 1, 1, 2008 )
      ENDIF

      ! Read Rice Emissions
#if   defined( GRID05x0666 ) && defined( NESTED_NA )
      FILENAME = TRIM( DATA_DIR ) // 'CH4_201305/'     //
     &           'rice.'          // GET_RES_EXT()     //
     &           '_NA.'           // CYEAR             //
     &           '.bpch' 
#else
      FILENAME = TRIM( DATA_DIR ) // 'CH4_201305/'     //
     &           'rice.'          // GET_RES_EXT()     //
     &           '.'              // CYEAR             //
     &           '.bpch'
#endif

      CALL READ_BPCH2( FILENAME, 'CH4-EMIS', 1,
     &                 XTAU,      IIPAR,     JJPAR,
     &                 1,         ARRAY,     QUIET=.TRUE.)
      CALL TRANSFER_2D( ARRAY, RICE_OUT )

      ! Reset CYEAR to current year
      WRITE( CYEAR, '(i4)' ) GET_YEAR()

      ! Get annual and monthly mean soil wetness from GEOS
      ! One file contains both monthly and annual mean GWETTOP
      YEAR = GET_YEAR()
#if   defined( GEOS_5 )
      IF ( YEAR .LT. 2004 ) YEAR  = 2004
      IF ( YEAR .GT. 2010 ) YEAR  = 2010
      WRITE( CYEAR, '(i4)' ) YEAR

      FILENAME = TRIM( DATA_DIR )  // 'CH4_201305/wetlands/'  //
     &           'GWETTOP.'        // GET_NAME_EXT()         // 
     &           '.'               // GET_RES_EXT()          // 
     &           '.'               // CYEAR                  // 
     &           '.bpch'

#elif defined( GEOS_4 )
      IF ( YEAR .LT. 2000 ) YEAR  = 2000
      IF ( YEAR .GT. 2006 ) YEAR  = 2006
      WRITE( CYEAR, '(i4)' ) YEAR

      FILENAME = TRIM( DATA_DIR )  // 'CH4_200911/GWETTOP/'  //
     &           'GWETTOP.'        // CYEAR                  // 
     &           '.'               // GET_NAME_EXT()         // 
     &           '.'               // GET_RES_EXT()
#endif

#if   defined( GRID05x0666 ) && defined( NESTED_NA )
      IF ( YEAR .LT. 2004 ) YEAR  = 2004
      IF ( YEAR .GT. 2004 ) YEAR  = 2009
      IF ( YEAR .GE. 2010 ) YEAR  = 2010
      WRITE( CYEAR, '(i4)' ) YEAR

      FILENAME = TRIM( DATA_DIR )  // 'CH4_201305/wetlands/'  //
     &           'GWETTOP.'        // GET_NAME_EXT()         // 
     &           '.'               // GET_RES_EXT()          // 
     &           '_NA.'            // CYEAR                  //
     &           '.bpch' 
#endif

      ! Annual mean GWETTOP
      XTAU = GET_TAU0( 1, 1, YEAR )
      CALL READ_BPCH2( FILENAME, 'GMAO-2D',  2,
     &                 XTAU,      IIPAR,     JJPAR,
     &                    1,      ARRAY,     QUIET=.TRUE.)
      CALL TRANSFER_2D( ARRAY, MEAN_GWETTOP )

      ! Monthly mean GWETTOP
      XTAU = GET_TAU0( GET_MONTH(), 1, YEAR )
      CALL READ_BPCH2( FILENAME, 'GMAO-2D',  1,
     &                 XTAU,      IIPAR,     JJPAR,
     &                    1,      ARRAY,     QUIET=.TRUE.)
      CALL TRANSFER_2D( ARRAY, MONTH_GWETTOP )

      !scale rice emissions (by Jerome Drevet)
      DO J = 1, JJPAR
      DO I = 1, IIPAR
	 wet_ratio = MONTH_GWETTOP(I,J)/MEAN_GWETTOP(I,J)-1
         wet_ratio = wet_ratio * 2.
         wet_ratio = wet_ratio +1.
	 if (wet_ratio < 0) wet_ratio = 0
         RICE_OUT(I,J)=RICE_OUT(I,J)*wet_ratio
      ENDDO
      ENDDO

      ! Set Wetland Emissions during this time step to WETL_OUT array
      WETL_OUT = EMIS3D(:,:,10)

      ! Subtract rice contribution from wetland source since wetland maps
      ! used count rice paddies in their area (I think). (kjw, 2,20,2010)
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         if (RICE_OUT(I,J) > 0) THEN   ! If rice > 0
         if (WETL_OUT(I,J) > 0) THEN   ! If wtl  > 0
         if (WETL_OUT(I,J) > RICE_OUT(I,J)) THEN
             WETL_OUT(I,J) = WETL_OUT(I,J) - RICE_OUT(I,J)
         endif
         endif
         endif
      enddo
      enddo


      EMIS3D(:,:,7) = RICE_OUT


      ! Return to calling program
      END SUBROUTINE RICE_EMIS


!------------------------------------------------------------------------------
!
!      SUBROUTINE BIOFUEL_EMIS( EMIS3D )
!!
!!******************************************************************************
!!  Subroutine BIOFUEL_EMIS calculates CH4 emissions from anthropogenic
!!  biofuels in the Yevich and Logan 2003 inventory (kjw, 6/03/09)
!!
!!  CO Emissions are read from the inventory of Yevich and Logan 2003.
!!  CH4 Emissions are calculated from emission factors in
!!  Andreae and Merlet, 2001 (6/4/09).
!!
!!  NOTES:
!!
!!******************************************************************************
!!
!      ! References to F90 modules
!      USE BPCH2_MOD,      ONLY : READ_BPCH2,    GET_TAU0
!      USE BPCH2_MOD,      ONLY : GET_RES_EXT
!      USE DIRECTORY_MOD,  ONLY : DATA_DIR
!      USE GRID_MOD,       ONLY : GET_AREA_CM2
!      USE TIME_MOD,       ONLY : GET_MONTH,     GET_YEAR
!      USE TIME_MOD,       ONLY : EXPAND_DATE
!      USE TRANSFER_MOD,   ONLY : TRANSFER_2D
!
!
!#     include "CMN_SIZE"     ! Size parameters
!#     include "CMN_DIAG"     ! Diagnostic switches
!#     include "CMN"          ! PD58
!
!      ! Arguments
!      REAL*8, INTENT(INOUT)     :: EMIS3D(IIPAR,JJPAR,PD58)
!
!      ! Local Variables
!      INTEGER                :: I, J, YYYY
!      REAL*4                 :: ARRAY(IIPAR,JJPAR)
!      REAL*8                 :: BIOF_OUT(IIPAR,JJPAR)
!      REAL*8                 :: AREA_CM2
!      REAL*8                 :: TAU0, TAU1, SECONDS
!      CHARACTER(LEN=255)     :: FILENAME
!
!
!      WRITE( 6, '(a)' ) '% EMISSCH4 --- BEGIN BIOFUEL_EMIS'
!
!      !=================================================================
!      ! BIOFUEL_EMIS begins here!
!      !=================================================================
!
!      !4.6 Biofuel emissions
!
!      !=================================================================
!      ! Read monthly biofuel emissions [kgCO/box/year]
!      !=================================================================
!
!      !TAU value for biofuel bpch files in data directory
!      TAU0     = GET_TAU0( 1, 1, 1985 )
!
!      ! File name with GFED2 C emissions
!      FILENAME = TRIM( DATA_DIR ) // 'biofuel_200202/biofuel.geos.' //
!     &              GET_RES_EXT()
!
!      ! Read CO Biofuel emissions [kg/box/year]
!      CALL READ_BPCH2( FILENAME, 'BIOFSRCE',   4,
!     &                 TAU0,      IIPAR,       JJPAR,
!     &                 1,         ARRAY,       QUIET=.TRUE. )
!      CALL TRANSFER_2D( ARRAY, BIOF_OUT )
!
!      !=================================================================
!      ! Convert [kgCO/box/year] to [molecCH4/cm2/s]
!      !=================================================================
!
!      ! [kgCO/box/year] --> [kgCO/cm2/year]
!      DO J = 1, JJPAR
!         AREA_CM2 = GET_AREA_CM2( J )
!         BIOF_OUT(:,J) = BIOF_OUT(:,J) / AREA_CM2
!      ENDDO
!
!      ! [kgCO/cm2/year] --> [kgCO/cm2/s]
!      YYYY = GET_YEAR()
!      TAU0 = GET_TAU0(1, 1, YYYY)
!      YYYY = YYYY + 1
!      TAU1 = GET_TAU0(1, 1, YYYY)
!      SECONDS  = ( TAU1 - TAU0 ) * 3600d0  ! # seconds in the current year
!      BIOF_OUT = BIOF_OUT / SECONDS
!
!      ! [kgCO/cm2/s] --> [kgCH4/cm2/s]
!      BIOF_OUT = BIOF_OUT * 61d-1 / 78d0   ! Andreae and Merlet 2001
!
!      ! [kgCH4/cm2/s] --> [molecCH4/cm2/s]
!      BIOF_OUT = BIOF_OUT * XNUMOL_CH4
!
!
!
!      EMIS3D(:,:,6) = BIOF_OUT
!
!
!      ! Return to calling program
!      END SUBROUTINE BIOFUEL_EMIS
!
!
!------------------------------------------------------------------------------

      SUBROUTINE ASEASONAL_ANTHRO_EMIS( EMIS3D )
!
!******************************************************************************
!
!  Subroutine ASEASONAL_ANTHRO_EMIS reads CH4 emissions from anthropogenic
!  sources. (kjw, 6/03/09)
!
!  Aseasonal anthropogenic emissions currently include EDGAR v4 categories
!  that are not called in their own subroutines.  Current emission categories
!  read in this subroutine are: gas & oil, coal, livestock, waste, and other
!  anthropogenic sources.
!
!  NOTES:
!  (1 )
!
!******************************************************************************

      ! References to F90 modules
      USE BPCH2_MOD,     ONLY : GET_RES_EXT,  GET_MODELNAME
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT_2D
      USE BPCH2_MOD,     ONLY : GET_TAU0,     READ_BPCH2
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE TIME_MOD,      ONLY : GET_YEAR
      USE TRANSFER_MOD,  ONLY : TRANSFER_2D
      !USE LOGICAL_MOD,   ONLY : LGAO,         LCOL,         LLIV
      !USE LOGICAL_MOD,   ONLY : LWAST,        LOTANT


#     include "CMN_SIZE"     ! Size parameters
#     include "CMN_DIAG"     ! Diagnostic switches
#     include "CMN"          ! PD58

      ! Arguments
      REAL*8, INTENT(INOUT)     :: EMIS3D(IIPAR,JJPAR,PD58)

      ! Local Variables
      REAL*4                 :: ARRAY(IIPAR,JJPAR)
      REAL*8                 :: XTAU
      REAL*8                 :: COAL_OUT(IIPAR,JJPAR)
      REAL*8                 :: GAO_OUT(IIPAR,JJPAR)
      REAL*8                 :: WST_OUT(IIPAR,JJPAR)
      REAL*8                 :: LIV_OUT(IIPAR,JJPAR)
      REAL*8                 :: OTH_OUT(IIPAR,JJPAR)
      REAL*8                 :: BFL_OUT(IIPAR,JJPAR)
      CHARACTER(LEN=255)     :: FILENAME
      CHARACTER(LEN=4  )     :: CYEAR


      WRITE( 6, '(a)' ) '% EMISSCH4 --- BEGIN ASEASONAL_ANTHRO_EMIS'

      !=================================================================
      ! ASEASONAL_ANTHRO_EMIS begins here!
      !=================================================================


      ! For now, we only have emissions from 2004
      WRITE( CYEAR, '(i4)' ) GET_YEAR()
      XTAU = GET_TAU0( 1, 1, GET_YEAR() )

      IF ( GET_YEAR() .LT. 2004 ) THEN
         CYEAR='2004'
         XTAU = GET_TAU0( 1, 1, 2004 )
      ENDIF
      IF ( GET_YEAR() .GT. 2008 ) THEN
         CYEAR='2008'
         XTAU = GET_TAU0( 1, 1, 2008 )
      ENDIF

      !4.2 Gas and Oil emissions
      !IF ( LGAO ) THEN
#if   defined( GRID05x0666 ) && defined( NESTED_NA )
         FILENAME = TRIM( DATA_DIR ) // 'CH4_201305/'     //
     &              'oil_gas.'       // GET_RES_EXT()     //
     &              '_NA.'           // CYEAR             //
     &              '.bpch'
#else
         FILENAME = TRIM( DATA_DIR ) // 'CH4_201305/'     //
     &              'oil_gas.'       // GET_RES_EXT()    //
     &              '.'              // CYEAR            //
     &              '.bpch'
#endif

         CALL READ_BPCH2( FILENAME, 'CH4-EMIS',   1,
     &                       XTAU,      IIPAR,    JJPAR,
     &                       1,         ARRAY,    QUIET=.TRUE.)
         CALL TRANSFER_2D( ARRAY, GAO_OUT )
      !ENDIF

      !4.3 Coal Mine emissions
      !IF ( LCOL ) THEN
#if   defined( GRID05x0666 ) && defined( NESTED_NA )
         FILENAME = TRIM( DATA_DIR ) // 'CH4_201305/'     //
     &              'coal.'          // GET_RES_EXT()     //
     &              '_NA.'           // CYEAR             //
     &              '.bpch'
#else
         FILENAME = TRIM( DATA_DIR ) // 'CH4_201305/'     //
     &              'coal.'          // GET_RES_EXT()     //
     &              '.'              // CYEAR             //
     &              '.bpch'
#endif

         CALL READ_BPCH2( FILENAME, 'CH4-EMIS',   1,
     &                       XTAU,      IIPAR,    JJPAR,
     &                       1,         ARRAY,    QUIET=.TRUE.)
         CALL TRANSFER_2D( ARRAY, COAL_OUT )
      !ENDIF

      !4.4 Livestock emissions
      !IF ( LLIV ) THEN
#if   defined( GRID05x0666 ) && defined( NESTED_NA )
         FILENAME = TRIM( DATA_DIR ) // 'CH4_201305/'     //
     &              'livestock.'     // GET_RES_EXT()     //
     &              '_NA.'           // CYEAR             //
     &              '.bpch'
#else
         FILENAME = TRIM( DATA_DIR ) // 'CH4_201305/'     //
     &              'livestock.'     // GET_RES_EXT()     //
     &              '.'              // CYEAR             //
     &              '.bpch'
#endif

         CALL READ_BPCH2( FILENAME, 'CH4-EMIS',   1,
     &                       XTAU,      IIPAR,    JJPAR,
     &                       1,         ARRAY,    QUIET=.TRUE.)
         CALL TRANSFER_2D( ARRAY, LIV_OUT )
      !ENDIF


      !4.5 Waste emissions
      !IF ( LWAST ) THEN
#if   defined( GRID05x0666 ) && defined( NESTED_NA )
         FILENAME = TRIM( DATA_DIR ) // 'CH4_201305/'     //
     &              'waste.'         // GET_RES_EXT()     //
     &              '_NA.'           // CYEAR             //
     &              '.bpch'
#else
         FILENAME = TRIM( DATA_DIR ) // 'CH4_201305/'     //
     &              'waste.'         // GET_RES_EXT()     //
     &              '.'              // CYEAR             //
     &              '.bpch'
#endif

         CALL READ_BPCH2( FILENAME, 'CH4-EMIS',   1,
     &                       XTAU,      IIPAR,    JJPAR,
     &                       1,         ARRAY,    QUIET=.TRUE.)
         CALL TRANSFER_2D( ARRAY, WST_OUT )
      !ENDIF

      !4.6 Biofuel emissions
      !IF ( LBFCH4 ) THEN
#if   defined( GRID05x0666 ) && defined( NESTED_NA )
         FILENAME = TRIM( DATA_DIR ) // 'CH4_201305/'     //
     &              'resident.'      // GET_RES_EXT()     //
     &              '_NA.'           // CYEAR             //
     &              '.bpch'
#else
         FILENAME = TRIM( DATA_DIR ) // 'CH4_201305/'     //
     &              'resident.'      // GET_RES_EXT()     //
     &              '.'              // CYEAR             //
     &              '.bpch'
#endif

         CALL READ_BPCH2( TRIM(FILENAME), 'CH4-EMIS',   1,
     &        XTAU,      IIPAR,    JJPAR,
     &        1,         ARRAY,    QUIET=.TRUE.)

         CALL TRANSFER_2D( ARRAY, BFL_OUT )
      !ENDIF

      !4.8 Other Anthropogenic Emissions
      !IF ( LOTANT ) THEN
#if   defined( GRID05x0666 ) && defined( NESTED_NA )
         FILENAME = TRIM( DATA_DIR ) // 'CH4_201305/'     //
     &              'other.'         // GET_RES_EXT()     //
     &              '_NA.'           // CYEAR             //
     &              '.bpch'
#else
         FILENAME = TRIM( DATA_DIR ) // 'CH4_201305/'     //
     &              'other.'         // GET_RES_EXT()     //
     &              '.'              // CYEAR             //
     &              '.bpch'
#endif

         CALL READ_BPCH2( FILENAME, 'CH4-EMIS',   1,
     &                       XTAU,      IIPAR,    JJPAR,
     &                       1,         ARRAY,    QUIET=.TRUE.)
         CALL TRANSFER_2D( ARRAY, OTH_OUT )
      !ENDIF


      ! Add emissions to EMIS3D array
      EMIS3D(:,:,2) = GAO_OUT(:,:)
      EMIS3D(:,:,3) = COAL_OUT(:,:)
      EMIS3D(:,:,4) = LIV_OUT(:,:)
      EMIS3D(:,:,5) = WST_OUT(:,:)
      EMIS3D(:,:,6) = BFL_OUT(:,:)
      EMIS3D(:,:,8) = OTH_OUT(:,:)

      ! Return to calling program
      END SUBROUTINE ASEASONAL_ANTHRO_EMIS


!------------------------------------------------------------------------------

      SUBROUTINE ASEASONAL_NATURAL_EMIS( EMIS3D )
!
!******************************************************************************
!  Subroutine ASEASONAL_NATURAL_EMIS reads CH4 emissions from natural sources.
!  (kjw, 6/03/09)
!
!  Aseasonal natural emissions currently include termites (Fung et. al. 1991)
!  and soil absorption (Fung et. al. 1991).  Future additions may include
!  emissions from permafrost, clathrates, thermokarst lakes, or
!  geothermal vents.
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,     ONLY : GET_RES_EXT,  GET_MODELNAME
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT_2D
      USE BPCH2_MOD,     ONLY : GET_TAU0,     READ_BPCH2
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE TIME_MOD,      ONLY : GET_YEAR
      USE TRANSFER_MOD,  ONLY : TRANSFER_2D
      !USE LOGICAL_MOD,   ONLY : LSOABS,       LOTNAT

#     include "CMN_SIZE"     ! Size parameters
#     include "CMN_DIAG"     ! Diagnostic switches
#     include "CMN"          ! PD58

      ! Arguments
      REAL*8, INTENT(INOUT)     :: EMIS3D(IIPAR,JJPAR,PD58)

      ! Local Variables
      REAL*4                 :: ARRAY(IIPAR,JJPAR)
      REAL*8                 :: SOIL_OUT(IIPAR,JJPAR)
      REAL*8                 :: OTH_OUT(IIPAR,JJPAR)
      REAL*8                 :: XTAU
      CHARACTER(LEN=255)     :: FILENAME
      CHARACTER(LEN=4  )     :: CYEAR


      WRITE( 6, '(a)' ) '% EMISSCH4 --- BEGIN ASEASONAL_NATURAL_EMIS'

      !=================================================================
      ! ASEASONAL_NATURAL_EMIS begins here!
      !=================================================================


      ! We only have one year of soil absorption and other natural
      ! CH4 emissions.  These have the date Jan . 1, 1985
      XTAU = GET_TAU0( 1, 1, 1985 )

      !4.11 Soil Absorption
      !IF ( LSOABS ) THEN
#if   defined( GRID05x0666 ) && defined( NESTED_NA )
         FILENAME = TRIM( DATA_DIR ) // 'CH4_201203/'     //
     &              'soilabs.'       // GET_RES_EXT()     // 
     &              '_NA.bpch'
#else
         FILENAME = TRIM( DATA_DIR ) // 'CH4_200911/'     //
     &              'soilabs.'       // GET_NAME_EXT_2D() //
     &              '.'              // GET_RES_EXT() 
#endif

         CALL READ_BPCH2( FILENAME, 'CH4-EMIS', 1,
     &                       XTAU,      IIPAR,    JJPAR,
     &                       1,         ARRAY,    QUIET=.TRUE.)
         CALL TRANSFER_2D( ARRAY, SOIL_OUT )
      !ENDIF

      !4.12 Other Natural Emissions
      !IF ( LOTNAT ) THEN
#if   defined( GRID05x0666 ) && defined( NESTED_NA )
         FILENAME = TRIM( DATA_DIR ) // 'CH4_201203/'     //
     &              'termites.'      // GET_RES_EXT()     //
     &              '_NA.bpch'
#else
         FILENAME = TRIM( DATA_DIR ) // 'CH4_200911/'     //
     &              'termites.'      // GET_NAME_EXT_2D() //
     &              '.'              // GET_RES_EXT() 
#endif

         CALL READ_BPCH2( FILENAME, 'CH4-EMIS', 1,
     &                       XTAU,      IIPAR,    JJPAR,
     &                       1,         ARRAY,    QUIET=.TRUE.)
         CALL TRANSFER_2D( ARRAY, OTH_OUT )
      !ENDIF


      ! Add emissions to array
      EMIS3D(:,:,11) = SOIL_OUT
      EMIS3D(:,:,12) = OTH_OUT


      ! Return to calling program
      END SUBROUTINE ASEASONAL_NATURAL_EMIS


!------------------------------------------------------------------------------

      SUBROUTINE CHEMCH4
!
!******************************************************************************
!  Subroutine CHEMCH4 computes the chemical loss of CH4 (sources - sinks).
!  (jsw, bnd, bmy, 6/8/00, 10/3/05)
!
!  CH4 SOURCES
!  ============================================================================
!  (1 ) Oxidation of methane, isoprene and monoterpenes (SRCO_fromHCs).
!  (2 ) Direct emissions of CO from fossil fuel combustion, biomass
!        burning and wood (for fuel) burning (SR SETEMIS).
!  (3 ) Emissions.
!
!  CH4 SINKS:
!  ============================================================================
!  (1 ) Removal of CO by OH (SR OHparam & CO_decay).
!  (2 ) CO uptake by soils (neglected).
!  (3 ) Transport of CO to stratosphere from troposphere
!        (in dynamical subroutines).
!  (4 ) Removal by OH (Clarissa's OH--climatol_OH.f and CO_decay.f)
!  (5 ) Transport of CH4 between troposphere and stratosphere, and
!        destruction in strat (CH4_strat.f).
!
!  NOTES:
!  (1 ) Created by Bryan Duncan (1/99).  Adapted for CH4 chemistry by
!        James Wang (6/8/00).  Inserted into module "global_ch4_mod.f"
!        by Bob Yantosca. (bmy, 1/16/01)
!  (2 ) CHEMCH4 is independent of "CMN_OH", "CMN_CO", and "CMN_CO_BUDGET".
!        (bmy, 1/16/01)
!  (3 ) Updated comments (jsw, bmy, 2/12/01)
!  (4 ) LD43 is already declared in CMN_DIAG; don't redefine it (bmy, 11/15/01)
!  (5 ) Replaced all instances of IM with IIPAR and JM with JJPAR, in order
!        to prevent namespace confusion for the new TPCORE (bmy, 6/25/02)
!  (6 ) Now reference AD from "dao_mod.f".  Now reference GEOS_CHEM_STOP from
!        "error_mod.f"  Now make FIRSTCHEM a local SAVEd variable.  Now
!        reference ALBD from "dao_mod.f".  Now use MONTH and JDATE from "CMN"
!        instead of LMN and LDY. (bmy, 11/15/02)
!  (7 ) Remove NYMDb, NYMDe from the arg list.  Now use functions GET_MONTH,
!        GET_NYMDb, GET_NYMDe, GET_MONTH, GET_DAY from the new "time_mod.f"
!        (bmy, 3/27/03)
!  (8 ) Now reference DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!  (9 ) Remove reference to BPCH2_MOD, it's not needed (bmy, 10/3/05)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,     ONLY : READ_BPCH2, GET_TAU0
      USE DAO_MOD,       ONLY : AD
      USE DIAG_MOD,      ONLY : AD43
      USE DIAG_PL_MOD,   ONLY : AD65
      USE ERROR_MOD,     ONLY : GEOS_CHEM_STOP, IT_IS_NAN, IT_IS_FINITE
      USE TIME_MOD,      ONLY : GET_DAY, GET_MONTH, GET_NYMDb, GET_NYMDe
      USE TRACER_MOD,    ONLY : STT
      USE LOGICAL_MOD,   ONLY : LSPLIT, LCH4BUD
      USE GLOBAL_OH_MOD, ONLY : GET_GLOBAL_OH, OH
      USE TRANSFER_MOD,  ONLY : TRANSFER_3D

#     include "CMN_SIZE"     ! Size parameters
#     include "CMN"          ! LPAUSE
#     include "CMN_DIAG"     ! ND43, AD43

      ! Local variables
      LOGICAL                :: FIRSTCHEM = .TRUE.
      INTEGER                :: I, J, L, K, M, N
      INTEGER                :: IJ, JJ, NPART, III, JJJ
      INTEGER                :: NOHDO
      INTEGER, SAVE          :: NTALDT

      CHARACTER(LEN=255)     :: FILENAME
      CHARACTER(LEN=4)       :: CYEAR
      REAL*4                 :: ARRAY(IIPAR,JJPAR,LGLOB)
      INTEGER                :: TROPP
      REAL*8                 :: XTAU
      INTEGER                :: LMN
      REAL*8                 :: PREVCH4(IIPAR, JJPAR, LLPAR)


      ! Number of days per month
      INTEGER                :: NODAYS(12) = (/ 31, 28, 31, 30,
     &                                          31, 30, 31, 31,
     &                                          30, 31, 30, 31 /)

      ! External functions
      REAL*8 , EXTERNAL      :: BOXVL

      ! Weight of air (taken from "comode.h")
      REAL*8, PARAMETER      :: WTAIR = 28.966d0

      !=================================================================
      ! CHEMCH4 begins here!
      !=================================================================
      WRITE( 6, '(a)' ) '% --- ENTERING CHEMCH4! ---'

      !=================================================================
      ! (0) Calculate each box's air density [molec/cm3]
      !        do this for saving mean OH concentrations (kjw, 6/12/09)
      !=================================================================

      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         BAIRDENS(I,J,L) = AD(I,J,L) * 1000d0   / BOXVL(I,J,L) *
     &                                 6.023D23 / WTAIR
      ENDDO
      ENDDO
      ENDDO

      !=================================================================
      ! (1) If the first time step ...
      !=================================================================
      IF ( FIRSTCHEM ) THEN

         ! Counter for total number of timesteps per month for CO budget.
         NTALDT = 1

         ! Now read CH4 loss frequencies instead of CO production
         ! (kjw, 12/2/11)
         ! Zero CO Production array
         !COPROD(:,:,:) = 0d0
         !print*,'READ_COPROD'
         !! Read zonally-averaged CO production [v/v/s]
         !CALL READ_COPROD
         !print*,'READ_COPROD DONE'
         CH4LOSS(:,:,:,:) = 0d0
         CALL READ_CH4LOSS

         ! Added following line to increase strat. sink strength
         ! Hmm, the values I printed out above for COprod are very small,
         ! all less than 1e-15. (jsw)
c         COprod = COprod * 3d0
         ! Commented the above line because it was giving me negative CH4
         ! concentrations in the stratosphere.

         ! Initialize the CH4 burden TCH4
         ! (ccc, 7/23/09)
         TCH4(:,:,:,1) = STT(:,:,:,1) * XNUMOL_CH4
      ENDIF

      ! Initialize current month
      LMN = GET_MONTH()

      ! Increment counter of timesteps
      NTALDT = NTALDT + 1

      !=================================================================
      ! (2) Calculate the production and destruction of CO from
      !     gas-phase chemistry only.
      !
      ! Concerning O3, there are 3 options:   if (m lt 9) then MM_add = '0'
      !                                       else MM_add = ''
      !   A) The OH parameterization is calculated using GEOS monthly
      !      means (NCLIMATOLOGY=0) for the independent variable O3.
      !      The O3 column above independent variable is determined
      !      using jal's O3  climatologies for both the tropospheric
      !      and stratospheric portions of the O3 column
      !      (NCLIMATOLOGY2=1).
      !
      !   B) The O3 variable is determined from jal's O3 climatolgies
      !      (tropospheric portion) and the o3 column above variable
      !      is determined from jal's O3 climatolgies (NCLIMATOLOGY=1 &
      !      NCLIMATOLOGY2=1).
      !
      !=================================================================

      !================================================================
      ! (3) get parameterized OH fields or monthly mean fields.
      !
      ! Variables of note:
      ! ---------------------------------------------------------------
      ! (1) BOH = storage array for OH fields.
      !
      ! (2) NOHDO = switch
      !       ONLY USE CASE 1 as of 5/28/08 (kjw)
      !       = 1 : Get GEOS-Chem OH (v5-07-08) (kjw, 5/28/08)
      !
      ! (3) LPAUSE =  the vertical level of the tropopause.  Above this
      !     level, no [OH] is calculated.  The user can feed this
      !     SR a high value for LPAUSE which effectively turns this
      !     option off (i.e., LPAUSE > MVRTBX). If the [OH] = -999
      !     then the [OH] was not calculated.
      !================================================================

      ! Change value of NOHDO as listed above
      NOHDO = 1

      SELECT CASE ( NOHDO )

         ! NOHDO = 1: GEOS-Chem OH v5-07-08
         CASE ( 1 )

            ! If first of month, read monthly mean OH
            IF ( FIRSTCHEM ) THEN

               ! 3D OH Field
               BOH(:,:,:,:) = 0d0

               ! Loop over each month, reading OH
               DO M=1,12

                  ! Global OH
                  CALL GET_GLOBAL_OH( M )

                  ! Assign to module variable BOH
                  BOH(:,:,:,M) = OH(:,:,:)

               ENDDO

            ENDIF

         CASE DEFAULT
            WRITE( 6, '(a)' ) 'Invalid selection for NOHDO!'
            WRITE( 6, '(a)' ) 'Halting execution in CHEMCH4!'
            CALL GEOS_CHEM_STOP

      END SELECT

      !=================================================================
      ! (3.1) ND43 diagnostics...save [OH] in molecules/cm3
      !=================================================================

      IF ( ND43 > 0 ) THEN
         DO L = 1, LD43
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            IF ( L < LPAUSE(I,J) ) THEN
               AD43(I,J,L,1) = AD43(I,J,L,1) + BOH(I,J,L,LMN)
            ENDIF
         ENDDO
         ENDDO
         ENDDO
      ENDIF

      !=================================================================
      ! (4) Save OH concentrations for printing of global mean [OH] at
      !     end of simulation.
      !=================================================================
      CALL CH4_OHSAVE

      !=================================================================
      ! (5) If multi-CH4 tracers, we store the CH4 total conc. to
      !     distribute the sink after the chemistry. (ccc, 2/10/09)
      !=================================================================
      IF ( LSPLIT ) THEN
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            PREVCH4(I,J,L) = STT(I,J,L,1)
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF

      !=================================================================
      ! (6) calculate rate of decay of CH4 by OH oxidation.
      !=================================================================
      CALL CH4_DECAY

      !=================================================================
      ! (7) calculate CH4 chemistry in layers above tropopause.
      !=================================================================
      CALL CH4_STRAT

      !=================================================================
      ! (8) distribute the chemistry sink from total CH4 to other CH4
      !     tracers. (ccc, 2/10/09)
      !=================================================================
      IF ( LSPLIT ) THEN
         CALL CH4_DISTRIB(PREVCH4)
      ENDIF

      !=================================================================
      ! (9) write budget (i.e., monthly average fields).
      !
      ! Check to make sure the start and end times are on the
      ! first of a month.  If not the SR CO_budget will not
      ! work properly!
      !=================================================================
      NPART = GET_NYMDb() / 100

      IF ( LCH4BUD .and. ( GET_NYMDb() - NPART*100 ) /= 1 ) THEN
         print*,'Start date not equal to 1st of month!!!'
         print*,'  Therefore, SR CO_budget will not work!!!'
         CALL GEOS_CHEM_STOP
      ENDIF


      NPART = GET_NYMDe() /100

      IF ( LCH4BUD .and. ( GET_NYMDe() - NPART*100 ) /= 1 ) THEN
         print*,'End date not equal to 1st of month!!!'
         print*,'  Therefore, SR CO_budget will not work!!!'
         CALL GEOS_CHEM_STOP
      ENDIF

      ! Call CH4_BUDGET on the last day of the month
      IF ( LCH4BUD .and. GET_DAY() == NODAYS( GET_MONTH() ) ) THEN
         CALL FLUSH ( 6 )

         CALL CH4_BUDGET

	 NTALDT  = 0
         call flush(6)

      ENDIF

      call flush(6)

      ! Set FIRSTCHEM to FALSE
      FIRSTCHEM = .FALSE.

      ! Return to calling program
      END SUBROUTINE CHEMCH4

!------------------------------------------------------------------------------

      SUBROUTINE READ_COPROD
!
!*****************************************************************************
!  Subroutine READ_COPROD reads production and destruction rates for CO in
!  the stratosphere. (bnd, bmy, 1/17/01, 10/3/05)
!
!  Module Variables:
!  ===========================================================================
!  (1) COPROD (REAL*8) : Array containing P(CO) for all 12 months [v/v/s]
!
!  NOTES:
!  (1 ) Created by Bryan Duncan (1/99).  Adapted for CH4 chemistry by
!        James Wang (6/8/00).  Inserted into module "global_ch4_mod.f"
!        by Bob Yantosca. (bmy, 1/16/01)
!  (2 ) READ_COPROD is independent of "CMN_OH", "CMN_CO", and "CMN_CO_BUDGET".
!        (bmy, 1/16/01)
!  (3 ) ARRAY needs to be dimensioned (1,JJPAR,LGLOB) (bmy, 9/26/01)
!  (4 ) Remove obsolete code from 9/01 (bmy, 10/24/01)
!  (5 ) Now reference DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!  (6 ) Now reads data for both GEOS and GCAP grids (bmy, 8/16/05)
!  (7 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!*****************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT_2D, GET_RES_EXT
      USE BPCH2_MOD,     ONLY : GET_TAU0,        READ_BPCH2
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT,    GET_MODELNAME
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE TRANSFER_MOD,  ONLY : TRANSFER_ZONAL

      IMPLICIT NONE

#     include "CMN_SIZE"  ! Size parameters


      ! Local variables
      INTEGER            :: I, J, L, M
      REAL*4             :: ARRAY(1,JJPAR,LGLOB)
      REAL*4             :: DUMMY_IN(JJPAR,LGLOB)
      REAL*8             :: XTAU
      CHARACTER(LEN=255) :: FILENAME
      REAL*8             :: DUMMY_OUT(JJPAR,LGLOB)


      !=================================================================
      ! READ_COPROD begins here!
      !
      ! Read P(CO) for all 12 months
      !=================================================================
      DO M = 1, 12

         ! TAU value at the start of month M -- Use "generic" year 1985
         XTAU = GET_TAU0( M, 1, 1985 )

         ! Construct filename
#if defined( GRID05x0666 ) && defined( NESTED_NA )
         FILENAME = TRIM( DATA_DIR ) // 'pco_lco_200203/'  //
     &              'COprod.GEOS5.05x0666_NA.trim'
#else
         FILENAME = TRIM( DATA_DIR ) // 'pco_lco_200203/' //
     &              'COprod.'        // GET_NAME_EXT()    //
     7              '.'              // GET_RES_EXT()
#endif

         WRITE( 6, 93 ) TRIM ( FILENAME )
 93      FORMAT( '     - READ_COPROD: Reading COprod: ', a )
	 CALL FLUSH( 6 )

         CALL READ_BPCH2( TRIM(FILENAME), 'PORL-L=$', 9,
     &                    XTAU,      1,         JJPAR,
     &                    LLPAR,     ARRAY,     QUIET=.TRUE. )

         ! use 2D arrays for TRANSFER ZONAL
         DUMMY_IN(:,:) = ARRAY(1,:,:)

         ! Copy REAL*4 to REAL*8 data, and resize from (JJPAR,LGLOB)
         ! to (JJPAR,LLPAR) -- vertically regrid if necessary
         CALL TRANSFER_ZONAL( DUMMY_IN, DUMMY_OUT )

         COPROD(:,:,M) = DUMMY_OUT(:,:)

      ENDDO


      ! Return to calling program
      END SUBROUTINE READ_COPROD


!------------------------------------------------------------------------------

      SUBROUTINE READ_CH4LOSS
!
!*****************************************************************************
!  Subroutine READ_CH4LOSS reads CH4 loss frequencies in the stratosphere.
!  These values constitute a linearized stratospheric CH4 chemistry scheme.
!  Loss frequencies from 4x5 degree output from the GMI model. Thanks to Lee
!  Murray for the ch4 loss frequencies. (kjw, 11/19/2011)
!
!  Module Variables:
!  ===========================================================================
!  (1) CH4LOSS (REAL*8) : Array containing ch4 loss frequencies for all 12 months [1/s]
!
!  NOTES:
!  (1 ) Created by Bryan Duncan (1/99).  Adapted for CH4 chemistry by
!        James Wang (6/8/00).  Inserted into module "global_ch4_mod.f"
!        by Bob Yantosca. (bmy, 1/16/01)
!  (2 ) READ_CH4LOSS is independent of "F77_CMN_OH", "F77_CMN_CO", and "F77_CMN_CO_BUDGET".
!        (bmy, 1/16/01)
!  (3 ) ARRAY needs to be dimensioned (1,JJPAR,LGLOB) (bmy, 9/26/01)
!  (4 ) Remove obsolete code from 9/01 (bmy, 10/24/01)
!  (5 ) Now reference DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!  (6 ) Now reads data for both GEOS and GCAP grids (bmy, 8/16/05)
!  (7 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (8 ) Treat MERRA in the same way as for GEOS-5 (bmy, 8/13/10)
!*****************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT_2D, GET_RES_EXT
      USE BPCH2_MOD,     ONLY : GET_TAU0,        READ_BPCH2
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT,    GET_MODELNAME
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE TRANSFER_MOD,  ONLY : TRANSFER_3D


      IMPLICIT NONE
#     include "define.h"
#     include "CMN_SIZE"


      ! Local variables
      INTEGER            :: I, J, L, M
      REAL*4             :: ARRAY(IIPAR,JJPAR,LLPAR)
      REAL*8             :: XTAU
      CHARACTER(LEN=255) :: FILENAME

      !=================================================================
      ! READ_CH4LOSS begins here!
      !
      ! Read P(CO) for all 12 months
      !=================================================================
      ! Construct filename
#if   defined( GRID05x0666 ) && defined( NESTED_NA )
      FILENAME = TRIM( DATA_DIR )         // 'CH4_201203/' //
     &           'gmi.ch4loss.geos5_47L.' // GET_RES_EXT() //
     &           '_NA.bpch'
#else
      FILENAME = TRIM( DATA_DIR )         // 'CH4_201203/' //
     &           'gmi.ch4loss.geos5_47L.' // GET_RES_EXT() // 
     &           '.bpch'
#endif

      WRITE( 6, 93 ) TRIM ( FILENAME )
 93   FORMAT( '     - READ_CH4LOSS: Reading Ch4loss: ', a )
      CALL FLUSH( 6 )

      ! Read data for each month
      DO M = 1, 12

         ! TAU value at the start of month M -- Use "generic" year 1985
         XTAU = GET_TAU0( M, 1, 1985 )

         ! Read Loss frequencies in units of [1/s].  drevet.
         CALL READ_BPCH2( FILENAME, 'IJ-AVG-$', 1,
     &                    XTAU,      IIPAR,     JJPAR,
     &                    LLPAR,     ARRAY,     QUIET=.TRUE. )

         ! Place array into CH4LOSS module variable
         CH4LOSS(:,:,:,M) = ARRAY(:,:,:)

      ENDDO

      ! Return to calling program
      END SUBROUTINE READ_CH4LOSS

!------------------------------------------------------------------------------

      SUBROUTINE CH4_DECAY
!
!******************************************************************************
!  Subroutine CH4_DECAY calculates the decay rate of CH4 by OH.  OH is the
!  only sink for CH4 considered here. (jsw, bnd, bmy, 1/16/01, 7/20/04)
!
!  The annual mean tropopause is stored in the LPAUSE array
!  (from header file "CMN").  LPAUSE is defined such that:
!
!  Levels            1 <= L <= LPAUSE(I,J) - 1 are tropospheric
!          LPAUSE(I,J) <= L <= LLPAR           are stratospheric
!
!  We now use LPAUSE instead of NSKIPL to denote the strat/trop boundary.
!  (bmy, 4/18/00)
!
!  Monthly loss of CH4 is summed in TCH4(3)
!     TCH4(3)  = CH4 sink by OH
!
!  Module Variables:
!  ============================================================================
!  (1) BOH        (REAL*8) : Array holding global OH concentrations
!  (2) XNUMOL_CH4 (REAL*8) : Molec CH4 / kg CH4
!
!  NOTES:
!  (1 ) Created by Bryan Duncan (1/99).  Adapted for CH4 chemistry by
!        James Wang (7/00).  Inserted into module "global_ch4_mod.f"
!        by Bob Yantosca. (bmy, 1/16/01)
!  (2 ) CH4_DECAY is independent of "CMN_OH", "CMN_CO", and "CMN_CO_BUDGET".
!        (bmy, 1/16/01)
!  (3 ) Now use function GET_TS_CHEM from "time_mod.f" (bmy, 3/27/03)
!  (4 ) Now references STT from "tracer_mod.f" (bmy, 7/20/04)
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,    ONLY : AIRVOL, T
      USE TIME_MOD,   ONLY : GET_TS_CHEM, ITS_A_NEW_YEAR
      USE TRACER_MOD, ONLY : STT
cdrevet
      USE DIAG_MOD,   ONLY : AD19
cdrevet
      USE TIME_MOD,   ONLY : GET_NYMD, GET_NHMS, GET_MONTH


#     include "CMN_SIZE"       ! Size parameters
#     include "CMN"            ! STT, LPAUSE
#     include "CMN_DIAG"       ! ND19

      ! Local variables
      LOGICAL          :: FIRST_DECAY=.TRUE.
      INTEGER          :: I, J, L, M, N, LMN
      REAL*8           :: DT, GCH4, STT2GCH4, KRATE

      ! External variables
      REAL*8, EXTERNAL :: BOXVL


      !=================================================================
      ! CH4_DECAY begins here!
      !=================================================================

      ! Chemistry timestep in seconds
      DT = GET_TS_CHEM() * 60d0

      ! Initialize current month
      LMN = GET_MONTH()

      !=================================================================
      ! Compute decay of CH4 by OH in the troposphere
      !
      ! The decay for CH4 is calculated by:
      !    OH + CH4 -> CH3 + H2O
      !    k = 2.45E-12 exp(-1775/T)
      !
      ! This is from JPL '97.
      ! JPL '00 & '06 do not revise '97 value. (jsw, kjw)
      !=================================================================
      IF (ITS_A_NEW_YEAR()) THEN
	TROPOCH4=0d0
      ENDIF

      DO L = 1, MAXVAL( LPAUSE )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Only consider tropospheric boxes
         IF ( L < LPAUSE(I,J) ) THEN

            ! Use 24-hr avg temperature to calc. rate coeff.
            ! citation needed
            KRATE = 2.45d-12 * EXP( -1775d0 / T(I,J,L) )

            ! Conversion from [kg/box] --> [molec/cm3]
            ! [kg CH4/box] * [box/cm3] * XNUMOL_CH4 [molec CH4/kg CH4]
            STT2GCH4 = 1d0 / AIRVOL(I,J,L) / 1d6 * XNUMOL_CH4

            ! CH4 in [molec/cm3]
            GCH4 = STT(I,J,L,1) * STT2GCH4

            ! Sum loss in TCH4(3) (molecules/box)
            TCH4(I,J,L,3) = TCH4(I,J,L,3)+
     &           ( GCH4 * BOXVL(I,J,L)* KRATE * BOH(I,J,L,LMN) * DT)

            TROPOCH4=TROPOCH4+GCH4*KRATE*BOH(I,J,L,LMN)*DT/STT2GCH4

            ! Modify AD19 Diagnostic
            ! How much CH4 (kg) is lost by reaction with OH
	    IF ( ND19 > 0 ) THEN  ! --> [kg/box]
	    	AD19(I,J,12) = AD19(I,J,12) +
     &	            ( GCH4 * KRATE * BOH(I,J,L,LMN) * DT ) / STT2GCH4
	    ENDIF

            ! Calculate new CH4 value: [CH4]=[CH4](1-k[OH]*delt)
            GCH4 = GCH4 * ( 1d0 - KRATE * BOH(I,J,L,LMN) * DT )

            ! Convert back from [molec/cm3] --> [kg/box]
            STT(I,J,L,1) = GCH4 / STT2GCH4

         ENDIF
      ENDDO
      ENDDO
      ENDDO


      print*,'% --- CHEMCH4: CH4_DECAY: TROP DECAY (Tg): ',TROPOCH4/1e9



      ! Return to calling program
      END SUBROUTINE CH4_DECAY

!------------------------------------------------------------------------------

      SUBROUTINE CH4_OHSAVE
!
!*****************************************************************************
!  Subroutine CH4_OHSAVE archives the CH3CCl3 lifetime from the OH
!  used in the CH4 simulation. (bnd, jsw, bmy, 1/16/01, 7/20/04)
!
!  Subroutine CH4_OHSAVE now ONLY archives OH concentrations to be printed
!  as global mean OH by PRINT_DIAG_OH at the end of the simulation.  The
!  CH3CCl3 lifetime capability was disabled many years ago. (kjw, 6/12/09)
!
!  The annual mean tropopause is stored in the LPAUSE array
!  (from header file "CMN").  LPAUSE is defined such that:
!
!  Levels            1 <= L <= LPAUSE(I,J) - 1 are tropospheric
!          LPAUSE(I,J) <= L <= LLPAR           are stratospheric
!
!  Module Variables
!  ===========================================================================
!  (1) BOH (REAL*8) : Array containing global OH field
!
!  NOTES:
!  (1 ) Created by Bryan Duncan (1/99).  Adapted for CH4 chemistry by
!        James Wang (7/00).  Inserted into module "global_ch4_mod.f"
!        by Bob Yantosca. (bmy, 1/16/01)
!  (2 ) CH4_OHSAVE is independent of "CMN_OH", "CMN_CO", and "CMN_CO_BUDGET".
!        (bmy, 1/16/01)
!  (3 ) Now call DO_DIAG_OH_CH4 to pass OH diagnostic info to the
!        "diag_oh_mod.f" (bmy, 7/20/04)
!*****************************************************************************
!
      ! References to F90 modules
      USE DIAG_OH_MOD, ONLY : DO_DIAG_OH_CH4
      USE TIME_MOD,    ONLY : GET_MONTH
      USE GRID_MOD,    ONLY : GET_AREA_CM2
      USE DAO_MOD,     ONLY : T
      USE TRACER_MOD,  ONLY : STT


#     include "CMN_SIZE"  ! Size parameters
#     include "CMN"       ! LPAUSE

      ! Local variables
      INTEGER          :: I, J, L, LMN
      REAL*8           :: MASST,   AREA_CM2
      REAL*8           :: KCLO,    LOSS,       OHMASS  
      REAL*8           :: KCH4,    CH4LOSS,    CH4MASS
      REAL*8           :: CH4EMIS, CH4TROPMASS

      ! External functions
      REAL*8, EXTERNAL :: BOXVL

      !=================================================================
      ! CH4_OHSAVE begins here!
      !
      ! (1) Pass OH mass, total air mass, and  to "diag_oh_mod.f"
      ! (2) ND59: Diagnostic for CH3CCl3 calculation
      !=================================================================

      ! Initialize current month
      LMN = GET_MONTH()

      ! 1. Calculate OH mass and total air mass
      DO L = 1, MAXVAL( LPAUSE )
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         ! Only process tropospheric boxes (bmy, 4/17/00)
         IF ( L < LPAUSE(I,J) ) THEN

            ! Calculate OH mass [molec / box]
            OHMASS = BOH(I,J,L,LMN) * BAIRDENS(I,J,L) * BOXVL(I,J,L,LMN)

            ! Calculate total air mass [molec / box]
            MASST  = BAIRDENS(I,J,L) * BOXVL(I,J,L,LMN)

            ! Calculate CH3CCl3 + OH rate constant from JPL '06
            ! [cm3 / molec / s]
            KCLO   = 1.64d-12 * EXP( -1520.d0 / T(I,J,L) )

            ! Calculate Loss term [molec / box / s]
            LOSS   = KCLO            * BOH(I,J,L,LMN)  *
     &               BAIRDENS(I,J,L) * BOXVL(I,J,L)


            ! Calculate CH4 emissions [molec / box / s]
            !   Only for surface level
            ! Grid box surface area [cm2]
            IF ( L .GT. 1 ) THEN 
               CH4EMIS = 0d0
            ELSE
               AREA_CM2 = GET_AREA_CM2( J )
               CH4EMIS  = SUM(CH4_EMIS(I,J,2:10)) + CH4_EMIS(I,J,12)
               CH4EMIS  = CH4EMIS * AREA_CM2 ! [molec/cm2/s] --> [molec/box/s]
            ENDIF

         ELSE

            OHMASS      = 0d0
            MASST       = 0d0
            LOSS        = 0d0
            CH4LOSS     = 0d0
            CH4TROPMASS = 0d0
            CH4EMIS     = 0d0
            CH4MASS     = STT(I,J,L,1) * XNUMOL_CH4 

         ENDIF

            ! Pass OH mass, total mass, and loss to "diag_oh_mod.f",
            ! which calculates mass-weighted mean [OH] and CH3CCl3
            ! lifetime.
            CALL DO_DIAG_OH_CH4( I, J, L, OHMASS, MASST, LOSS,
     &                CH4LOSS, CH4TROPMASS, CH4EMIS, CH4MASS )

      ENDDO
      ENDDO
      ENDDO


      ! Return to calling program
      END SUBROUTINE CH4_OHSAVE

!------------------------------------------------------------------------------

      SUBROUTINE CH4_STRAT
!
!*****************************************************************************
!  Subroutine CH4_STRAT calculates uses production rates for CH4 to
!  calculate loss of CH4 in above the tropopause.
!  (jsw, bnd, bmy, 1/16/01, 7/20/04)
!
!  Production (mixing ratio/sec) rate provided by Dylan Jones.
!  Only production by CH4 + OH is considered.
!
!  The annual mean tropopause is stored in the LPAUSE array
!  (from header file "CMN").  LPAUSE is defined such that:
!
!  Levels           1 <= L <= LPAUSE(I,J) - 1 are tropospheric
!         LPAUSE(I,J) <= L <= LLPAR           are stratospheric (bmy, 4/18/00)
!
!  NOTES:
!  (1 ) Created by Bryan Duncan (1/99).  Adapted for CH4 chemistry by
!        James Wang (7/00).  Inserted into module "global_ch4_mod.f"
!        by Bob Yantosca. (bmy, 1/16/01)
!  (2 ) CH4_STRAT is independent of "CMN_OH", "CMN_CO", and "CMN_CO_BUDGET".
!        (bmy, 1/16/01)
!  (3 ) Removed LMN from the arg list and made it a local variable.  Now use
!        functions GET_MONTH and GET_TS_CHEM from "time_mod.f" (bmy, 3/27/03)
!  (4 ) Now references STT from "tracer_mod.f" (bmy, 7/20/04)
!*****************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,    ONLY : AIRVOL
      USE TIME_MOD,   ONLY : GET_MONTH, GET_TS_CHEM
      USE TRACER_MOD, ONLY : STT

#     include "CMN_SIZE"       ! Size parameters
#     include "CMN"            ! STT, LPAUSE

      ! Local variables
      INTEGER             :: I, J, L, LMN
      REAL*8              :: DT, GCH4, STT2GCH4, LRATE
      CHARACTER*20        :: STT_TEST
      CHARACTER*20        :: STT2GCH4_CHAR

      ! External functions
      REAL*8, EXTERNAL    :: BOXVL

      !=================================================================
      ! CH4_STRAT begins here!
      !=================================================================

      ! Chemistry timestep [s]
      DT  = GET_TS_CHEM() * 60d0

      ! Current month
      LMN = GET_MONTH()

      !=================================================================
      ! Loop over stratospheric boxes only
      !=================================================================
      DO L = MINVAL( LPAUSE ), LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         IF ( L >= LPAUSE(I,J) ) THEN

            ! Conversion factor [kg/box] --> [molec/cm3]
            ! [kg/box] / [AIRVOL * 1e6 cm3] * [XNUMOL_CH4 molec/mole]
            STT2GCH4 = 1d0 / AIRVOL(I,J,L) / 1d6 * XNUMOL_CH4

            ! CH4 in [molec/cm3]
            GCH4 = STT(I,J,L,1) * STT2GCH4

            ! Loss rate [molec/cm3/s]
            LRATE = GCH4 * CH4LOSS( I,J,L,LMN )

            ! CH4 in [molec/cm3]
            GCH4 = GCH4 - ( LRATE * DT )

            ! Convert back from [molec CH4/cm3] --> [kg/box]
            STT(I,J,L,1) = GCH4 / STT2GCH4

         ENDIF
      ENDDO
      ENDDO
      ENDDO

      ! Return to calling program
      END SUBROUTINE CH4_STRAT

!------------------------------------------------------------------------------

      SUBROUTINE CH4_BUDGET
!
!******************************************************************************
!  Subroutine CH4_BUDGET calculates the budget of CH4.  This SR only works
!  for monthly averages, so be sure to start on the first of the month
!  and run to another first of the month!!!  (jsw, bnd, bmy, 1/16/01, 10/3/05)
!
!  Modified for the run with new emissions (j drevet, 03/06)
!
!  Store the sources/sinks of CH4 in TCH4 in total molecules
!           ( 1) = Initial burden
!           ( 2) = Final burden
!  SINKS
!           ( 3) = Tropospheric CH4 sink by OH
!  SOURCES
!           ( 4) = Total Sources
!           ( 5) = Industrial (Gas+Oil+Mine)
!           ( 6) = Agriculture (Enteric fermentation+Manure+Rice+Waste+Waste water)
!           ( 7) = Biomass burning
!           ( 8) = Termites
!           ( 9) = Wetland
!           (10) = Soil absorption
!           (11) = Interhemispheric Exchange (+ = northward)
!
!  NOTES:
!  (1 ) Created by Bryan Duncan (1/99).  Adapted for CH4 chemistry by
!        James Wang (7/00).  Inserted into module "global_ch4_mod.f"
!        by Bob Yantosca. (bmy, 1/16/01)
!  (2 ) CH4_BUDGET is independent of "CMN_OH", "CMN_CO", and "CMN_CO_BUDGET".
!        (bmy, 1/16/01)
!  (3 ) Updated comments (jsw, bmy, 2/13/01)
!  (4 ) Renamed XLABEL to LABEL so as not to conflict w/ "CMN"
!  (5 ) Now use functions GET_MONTH, GET_YEAR, GET_DIAGb, and GET_CT_DYN from
!        "time_mod.f".  Removed LMN from the arg list and made it a local
!        variable.  Use functions GET_XOFFSET and GET_YOFFSET from
!        "grid_mod.f".  (bmy, 3/27/03)
!  (6 ) Now references STT from "tracer_mod.f" (bmy, 7/20/04)
!  (7 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,  ONLY : BPCH2,       BPCH2_HDR,   GET_MODELNAME
      USE GRID_MOD,   ONLY : GET_XOFFSET, GET_YOFFSET
      USE TIME_MOD,   ONLY : GET_MONTH,   GET_YEAR
      USE TIME_MOD,   ONLY : GET_DIAGb,   GET_CT_DYN
      USE TRACER_MOD, ONLY : STT

#     include "CMN_SIZE"       ! Size parameters
#     include "CMN"            ! STT, LPAUSE

      ! Local variables
      INTEGER                :: I, J, K, L, M, NERROR, UD, LMN

      REAL*8                 :: STTCONV, TGS, SCALEDYN
      REAL*8                 :: NTP, NTQ, NTP2, NTQ2
      REAL*8                 :: SOURCES, SINKS

      CHARACTER(LEN=17)      :: MERGE
      CHARACTER(LEN=13)      :: MERGE2

      ! For binary punch file, v. 2.0
      REAL*4                 :: ARRAY(IIPAR, JJPAR, LLPAR)
      REAL*4                 :: LONRES, LATRES

      INTEGER                :: IFIRST, JFIRST, LFIRST
      INTEGER, PARAMETER     :: HALFPOLAR = 1
      INTEGER, PARAMETER     :: CENTER180 = 1

      CHARACTER (LEN=20)     :: MODELNAME
      CHARACTER (LEN=40)     :: UNIT
      CHARACTER (LEN=40)     :: RESERVED = ''
      CHARACTER (LEN=40)     :: CATEGORY
      CHARACTER (LEN=80)     :: LABEL

      ! External functions
      REAL*8, EXTERNAL       :: BOXVL

      !=================================================================
      ! CH4_BUDGET begins here!
      !
      ! Initialize quantities
      !=================================================================
      IFIRST    = GET_XOFFSET() + 1
      JFIRST    = GET_YOFFSET() + 1
      LFIRST    = 1
      LONRES    = DISIZE
      LATRES    = DJSIZE

      ! Current month
      LMN       = GET_MONTH()

      ! Make up a category name for GAMAP (use 8 characters)
      CATEGORY  = 'CH4BUDT'

      ! Get the proper model name for the binary punch file
      MODELNAME = GET_MODELNAME()

      ! Descriptor string
      LABEL    = 'GEOS-CHEM -- CH4 Budget output (jsw, bmy, 1/16/01)'

      ! Unit of quantity being saved
      UNIT      = 'Tg'  !(NOTE: check w/ bnd to get the right units!!!)

      ! Scale factor for dynamic time steps
      SCALEDYN  = FLOAT( GET_CT_DYN() ) + 1D-20

      !=================================================================
      ! Store the final burden of CH4 in TCH4(2)
      ! Convert kg CH4/box to molecules/box.
      !=================================================================
      TCH4(:,:,:,2) = 0d0
      TCH4(:,:,:,2) = STT(:,:,:,1) * XNUMOL_CH4

      !=================================================================
      ! Write GLOBAL AVERAGES for all layers to ASCII file
      !=================================================================
      WRITE( MERGE, 2 ) GET_MONTH(), GET_YEAR()
 2    FORMAT( 'CH4budget.', I2.2, '.',I4 )

      OPEN( 189, FILE=MERGE, STATUS='UNKNOWN' )
      REWIND( 189 )

      TGS     = 1.D-9
      STTCONV = XNUMOL_CH4/TGS
      SOURCES = 0.D0
      SINKS   = 0.D0
      NERROR  = 0

      WRITE(189,18)
      WRITE(189,1801)
 1801 FORMAT('*************************')
      WRITE(189,1800)
 1800 FORMAT('LAYERS 1 - 20')
      WRITE(189,1801)
      WRITE(189,18)

      WRITE(189,18)
      WRITE(189,38)
      WRITE(189,18)
      WRITE(189,19)
      WRITE(189,1990)
 1990 FORMAT('Tropospheric Burden')

      NTP=SUM_CH4(1,IIPAR,1,JJPAR,1,LLPAR,1,1,1)
      WRITE(189,20)NTP,NTP/STTCONV

      NTP2=SUM_CH4(1,IIPAR,1,JJPAR,1,LLPAR,2,2,1)
      WRITE(189,21)NTP2,NTP2/STTCONV

      WRITE(189,18)
      WRITE(189,1991)
 1991 FORMAT('Stratospheric Burden')

      NTP=SUM_CH4(1,IIPAR,1,JJPAR,1,LLPAR,1,1,0)
      WRITE(189,20) NTP,NTP/STTCONV

      NTP=SUM_CH4(1,IIPAR,1,JJPAR,1,LLPAR,2,2,0)
      WRITE(189,21) NTP,NTP/STTCONV

      WRITE(189,18)
      WRITE(189,31)

c Sinks   jsw has checked correctness of code for sinks.
      NTP=SUM_CH4(1,IIPAR,1,JJPAR,1,LLPAR,3,3,1)
      NTQ=SUM_CH4(1,IIPAR,1,JJPAR,1,LLPAR,3,3,0)
      SINKS=NTP+NTQ

      WRITE(189,22) NTP,NTP/SINKS*100.D0,NTP/STTCONV
      WRITE(189,220) NTQ,NTQ/SINKS*100.D0,NTQ/STTCONV
      WRITE(189,29)
      WRITE(189,34) SINKS,SINKS/STTCONV  !Just OH sink
      WRITE(189,18)
      WRITE(189,30)

C Sources
      NTP=SUM_CH4(1,IIPAR,1,JJPAR,1,1,4,4,1)
      SOURCES=NTP

      NTP=SUM_CH4(1,IIPAR,1,JJPAR,1,1,5,5,1)
      WRITE(189,24) NTP,NTP/SOURCES*100.D0,NTP/STTCONV

      NTP=SUM_CH4(1,IIPAR,1,JJPAR,1,1,9,9,1)
      WRITE(189,27) NTP,NTP/SOURCES*100.D0,NTP/STTCONV

      NTP=SUM_CH4(1,IIPAR,1,JJPAR,1,1,6,6,1)
      WRITE(189,39) NTP,NTP/SOURCES*100.D0,NTP/STTCONV

      NTP=SUM_CH4(1,IIPAR,1,JJPAR,1,1,7,7,1)
      WRITE(189,25) NTP,NTP/SOURCES*100.D0,NTP/STTCONV

      NTP=SUM_CH4(1,IIPAR,1,JJPAR,1,1,8,8,1)
      WRITE(189,26) NTP,NTP/SOURCES*100.D0,NTP/STTCONV

Cjsw Following lines added by jsw.
      NTP=SUM_CH4(1,IIPAR,1,JJPAR,1,1,10,10,1)
      WRITE(189,35) NTP,NTP/SOURCES*100.D0,NTP/STTCONV

      SINKS=SINKS-NTP  !Minus sign because soil absorption is negative.

      WRITE(189,29)
      WRITE(189,28) SOURCES,SOURCES/STTCONV
      WRITE(189,18)

      NTP=SUM_CH4(1,IIPAR,1,JJPAR,1,LLPAR,1,1,1)
      NTP2=SUM_CH4(1,IIPAR,1,JJPAR,1,LLPAR,2,2,1)
      NTQ=SUM_CH4(1,IIPAR,1,JJPAR,1,LLPAR,1,1,0)
      NTQ2=SUM_CH4(1,IIPAR,1,JJPAR,1,LLPAR,2,2,0)
      WRITE(189,18)
      WRITE(189,288) NTP-NTP2+NTQ-NTQ2+SOURCES-SINKS,
     *     (NTP-NTP2+NTQ-NTQ2+SOURCES-SINKS)/STTCONV
      WRITE(189,18)
      WRITE(189,289) -(NTP-NTP2+NTQ-NTQ2),
     *     -(NTP-NTP2+NTQ-NTQ2)/STTCONV

      !=================================================================
      ! Write SOUTHERN HEMISPHERE averages to ASCII file
      ! jsw:  I have not modified the remaining code for CH4.
      !=================================================================

      SOURCES = 0.D0
      SINKS   = 0.D0

      WRITE(189,18)
      WRITE(189,18)
      WRITE(189,36)
      WRITE(189,18)
      WRITE(189,19)
      WRITE(189,1990)

      NTP=SUM_CH4(1,IIPAR,1,JJPAR/2,1,LLPAR,1,1,1)
      WRITE(189,20) NTP,NTP/STTCONV

      NTP=SUM_CH4(1,IIPAR,1,JJPAR/2,1,LLPAR,2,2,1)
      WRITE(189,21) NTP,NTP/STTCONV

      WRITE(189,18)
      WRITE(189,1991)
      NTP=SUM_CH4(1,IIPAR,1,JJPAR/2,1,LLPAR,1,1,0)
      WRITE(189,20) NTP,NTP/STTCONV

      NTP=SUM_CH4(1,IIPAR,1,JJPAR/2,1,LLPAR,2,2,0)
      WRITE(189,21) NTP,NTP/STTCONV
      WRITE(189,18)
      WRITE(189,31)

      ! Sinks
      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR/2+1,1,LLPAR,11,11,1)
      IF( NTP > 0d0) SINKS = SINKS + NTP

      NTP2=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR/2+1,1,LLPAR,11,11,0)
      IF( NTP2 > 0d0 ) SINKS = SINKS + NTP2

      NTQ=SUM_CH4(1,IIPAR,1,JJPAR/2,1,LLPAR,3,3,1)
      NTQ2=SUM_CH4(1,IIPAR,1,JJPAR/2,1,LLPAR,3,3,0)
      SINKS=SINKS+NTQ+NTQ2
      WRITE(189,22) NTQ,NTQ/SINKS*100.D0,NTQ/STTCONV
      WRITE(189,220) NTQ2,NTQ2/SINKS*100.D0,NTQ2/STTCONV

      IF(NTP.GE.0.D0) THEN
         WRITE(189,270) NTP,NTP/SINKS*100.D0,NTP/STTCONV
      ENDIF

      IF(NTP.GE.0.D0) THEN
         WRITE(189,2700) NTP,NTP/SINKS*100.D0,NTP/STTCONV
      ENDIF

      WRITE(189,29)
      WRITE(189,34) SINKS,SINKS/STTCONV
      WRITE(189,18)
      WRITE(189,30)

      ! Sources
      NTQ=SUM_CH4(1,IIPAR,1,JJPAR/2,1,LLPAR,5,9,1)
      NTQ2=SUM_CH4(1,IIPAR,1,JJPAR/2,1,LLPAR,4,4,0)
      SOURCES=NTQ+NTQ2
      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR/2+1,1,LLPAR,11,11,1)
      IF(NTP.LT.0.D0) THEN
         SOURCES=SOURCES-NTP
      ENDIF

      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR/2+1,1,LLPAR,11,11,0)
      IF(NTP.LT.0.D0) THEN
         SOURCES=SOURCES-NTP
      ENDIF

      NTP=SUM_CH4(1,IIPAR,1,JJPAR/2,1,LLPAR,4,4,1)
      NTP2=SUM_CH4(1,IIPAR,1,JJPAR/2,1,LLPAR,4,4,0)
      NTP=SUM_CH4(1,IIPAR,1,JJPAR/2,1,1,5,5,1)
      WRITE(189,24) NTP,NTP/SOURCES*100.D0,NTP/STTCONV

      NTP=SUM_CH4(1,IIPAR,1,JJPAR/2,1,1,6,6,1)
      WRITE(189,39) NTP,NTP/SOURCES*100.D0,NTP/STTCONV

      NTP=SUM_CH4(1,IIPAR,1,JJPAR/2,1,1,7,7,1)
      WRITE(189,25) NTP,NTP/SOURCES*100.D0,NTP/STTCONV

      NTP=SUM_CH4(1,IIPAR,1,JJPAR/2,1,1,8,8,1)
      WRITE(189,26) NTP,NTP/SOURCES*100.D0,NTP/STTCONV

      NTP=SUM_CH4(1,IIPAR,1,JJPAR/2,1,1,9,9,1)
      WRITE(189,27) NTP,NTP/SOURCES*100.D0,NTP/STTCONV

      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR/2+1,1,LLPAR,11,11,1)
      IF(NTP.LT.0.D0) THEN
         WRITE(189,270) -NTP,-NTP/SOURCES*100.D0,-NTP/STTCONV
      ENDIF

      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR/2+1,1,LLPAR,11,11,0)
      NTP2=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR/2+1,1,LLPAR,11,11,0)
      IF( NTP2 > 0d0 ) SINKS = SINKS + NTP2

      NTQ=SUM_CH4(1,IIPAR,1,JJPAR/2,1,LLPAR,3,3,1)
      NTQ2=SUM_CH4(1,IIPAR,1,JJPAR/2,1,LLPAR,3,3,0)
      SINKS=SINKS+NTQ+NTQ2
      WRITE(189,22) NTQ,NTQ/SINKS*100.D0,NTQ/STTCONV

      IF(NTP.LT.0.D0) THEN
         WRITE(189,2700) -NTP,-NTP/SOURCES*100.D0,-NTP/STTCONV
      ENDIF

      WRITE(189,29)
      WRITE(189,28) SOURCES,SOURCES/STTCONV
      WRITE(189,18)

      NTP=SUM_CH4(1,IIPAR,1,JJPAR/2,1,LLPAR,1,1,1)
      NTP2=SUM_CH4(1,IIPAR,1,JJPAR/2,1,LLPAR,2,2,1)
      NTQ=SUM_CH4(1,IIPAR,1,JJPAR/2,1,LLPAR,1,1,0)
      NTQ2=SUM_CH4(1,IIPAR,1,JJPAR/2,1,LLPAR,2,2,0)
      WRITE(189,18)
      WRITE(189,288) (NTP-NTP2+NTQ-NTQ2+SOURCES-SINKS),
     *     (NTP-NTP2+NTQ-NTQ2+SOURCES-SINKS)/STTCONV
      WRITE(189,18)
      WRITE(189,289) -(NTP-NTP2+NTQ-NTQ2),
     *     -(NTP-NTP2+NTQ-NTQ2)/STTCONV

      !=================================================================
      ! Write NORTHERN HEMISPHERE averages to ASCII file
      ! jsw:  I have not modified the remaining code for CH4.
      !=================================================================

      SOURCES = 0.D0
      SINKS   = 0.D0

      WRITE(189,18)
      WRITE(189,18)
      WRITE(189,37)
      WRITE(189,18)
      WRITE(189,19)
      WRITE(189,1990)

      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR,1,LLPAR,1,1,1)
      WRITE(189,20) NTP,NTP/STTCONV
      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR,1,LLPAR,2,2,1)
      WRITE(189,21) NTP,NTP/STTCONV

      WRITE(189,18)
      WRITE(189,1991)

      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR,1,LLPAR,1,1,0)
      WRITE(189,20) NTP,NTP/STTCONV

      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR,1,LLPAR,2,2,0)
      WRITE(189,21) NTP,NTP/STTCONV

      WRITE(189,18)
      WRITE(189,31)
c Sinks
      NTQ=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR,1,LLPAR,3,3,1)
      NTQ2=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR,1,LLPAR,3,3,0)
      SINKS=NTQ+NTQ2

      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR/2+1,1,LLPAR,11,11,1)
      IF(NTP.LT.0.D0) THEN
         SINKS=SINKS-NTP
      ENDIF
      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR/2+1,1,LLPAR,11,11,0)
      IF(NTP.LT.0.D0) THEN
         SINKS=SINKS-NTP
      ENDIF

      WRITE(189,22) NTQ,NTQ/SINKS*100.D0,NTQ/STTCONV
      WRITE(189,220) NTQ2,NTQ2/SINKS*100.D0,NTQ2/STTCONV

      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR/2+1,1,LLPAR,11,11,1)
      IF(NTP.LT.0.D0) THEN
         WRITE(189,270) -NTP,-NTP/SINKS*100.D0,-NTP/STTCONV
      ENDIF

      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR/2+1,1,LLPAR,11,11,0)
      IF(NTP.LT.0.D0) THEN
         WRITE(189,2700) -NTP,-NTP/SINKS*100.D0,-NTP/STTCONV
      ENDIF

      WRITE(189,29)
      WRITE(189,34)SINKS,SINKS/STTCONV
      WRITE(189,18)
      WRITE(189,30)
C Sources
      NTQ=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR,1,LLPAR,5,9,1)
      NTQ2=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR,1,LLPAR,4,4,0)
      SOURCES=NTQ+NTQ2

      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR/2+1,1,LLPAR,11,11,1)
      IF(NTP.GE.0.D0) THEN
         SOURCES=SOURCES+NTP
      ENDIF
      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR/2+1,1,LLPAR,11,11,0)
      IF(NTP.GE.0.D0) THEN
         SOURCES=SOURCES+NTP
      ENDIF

      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR,1,LLPAR,4,4,1)

      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR,1,LLPAR,4,4,0)

      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR,1,1,5,5,1)
      WRITE(189,24) NTP,NTP/SOURCES*100.D0,NTP/STTCONV

      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR,1,1,6,6,1)
      WRITE(189,39) NTP,NTP/SOURCES*100.D0,NTP/STTCONV

      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR,1,1,7,7,1)
      WRITE(189,25) NTP,NTP/SOURCES*100.D0,NTP/STTCONV

      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR,1,1,8,8,1)
      WRITE(189,26) NTP,NTP/SOURCES*100.D0,NTP/STTCONV

      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR,1,1,9,9,1)
      WRITE(189,27) NTP,NTP/SOURCES*100.D0,NTP/STTCONV

      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR/2+1,1,LLPAR,11,11,1)
      IF(NTP.GE.0.D0) THEN
         WRITE(189,270) NTP,NTP/SOURCES*100.D0,NTP/STTCONV
      ENDIF

      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR/2+1,1,LLPAR,11,11,0)
      IF(NTP.GE.0.D0) THEN
         WRITE(189,2700) NTP,NTP/SOURCES*100.D0,NTP/STTCONV
      ENDIF

      WRITE(189,29)
      WRITE(189,28) SOURCES,SOURCES/STTCONV
      WRITE(189,18)

      NTP=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR,1,LLPAR,1,1,1)
      NTP2=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR,1,LLPAR,2,2,1)
      NTQ=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR,1,LLPAR,1,1,0)
      NTQ2=SUM_CH4(1,IIPAR,JJPAR/2+1,JJPAR,1,LLPAR,2,2,0)
      WRITE(189,18)
      WRITE(189,288) (NTP-NTP2+NTQ-NTQ2+SOURCES-SINKS),
     *     (NTP-NTP2+NTQ-NTQ2+SOURCES-SINKS)/STTCONV
      WRITE(189,18)
      WRITE(189,289) -(NTP-NTP2+NTQ-NTQ2),
     *     -(NTP-NTP2+NTQ-NTQ2)/STTCONV

 18   FORMAT()
 19   FORMAT('                    #Molecules               TG')
 20   FORMAT('  Start of Month  :',E10.3,10x,F10.3)
 21   FORMAT('  End of Month    :',E10.3,10x,F10.3)
 22   FORMAT('  CH4 decay-trop   :',E10.3,2x,F6.1,2x,F10.3)
 220  FORMAT('  CH4 decay-strat  :',E10.3,2x,F6.1,2x,F10.3)
 24   FORMAT('  Industrial      :',E10.3,2x,F6.1,2x,F10.3)
 25   FORMAT('  Biomass Burning :',E10.3,2x,F6.1,2x,F10.3)
 26   FORMAT('  Termites        :',E10.3,2x,F6.1,2x,F10.3)
 27   FORMAT('  Wetland         :',E10.3,2x,F6.1,2x,F10.3)
 270  FORMAT('  N-S Ex.-trop    :',E10.3,2x,F6.1,2x,F10.3)
 2700 FORMAT('  N-S Ex.-strat   :',E10.3,2x,F6.1,2x,F10.3)
 28   FORMAT('Total Sources     :',E10.3,10x,F10.3)
 288  FORMAT('Initial-Final+Sources-Sinks=',E10.3,2x,F10.3)
 289  FORMAT('Net Gain          : ',E10.3,10x,F10.3)
 29   FORMAT('                     ---------')
 30   FORMAT('SOURCES                          %Source')
 31   FORMAT('SINKS                            %Sink')
 34   FORMAT('Total Sinks       :',E10.3,10x,F10.3)
 35   FORMAT('  Soil absorption :',E10.3,2x,F6.1,2x,F10.3)
 39   FORMAT('  Agriculture     :',E10.3,2x,F6.1,2x,F10.3)

 36   FORMAT('*****  Southern Hemisphere  *****')
 37   FORMAT('*****  Northern Hemisphere  *****')
 38   FORMAT('*****  Global  *****')

      CLOSE(189)

!     !=================================================================
!     ! Also save to binary punch file. Don't save the bpunch file
!     ! anymore, because it's not used. Keep the code for reference.
!     ! The code creates the bpunch file fort.190. Should use a diag.
!     ! instead. (ccc, 8/14/09)
!     !=================================================================
!     CALL BPCH2_HDR( 190, LABEL )
!
!     DO K = 1, N_CH4
!
!        ! Cast REAL*8 into REAL*4, convert from molec to Tg
!        ARRAY(:,:,:) = TCH4(:,:,:,K) / STTCONV
!
!        ! Save the data block
!        CALL BPCH2( 190,       MODELNAME,   LONRES,      LATRES,
!    &               HALFPOLAR, CENTER180,   CATEGORY,    K,
!    &               UNIT,      GET_DIAGB(), GET_DIAGb(), RESERVED,
!    &               IIPAR,     JJPAR,       LLPAR,       IFIRST,
!    &               JFIRST,    LFIRST,      ARRAY )
!     ENDDO
!
!     CLOSE(190)

      !=================================================================
      ! Final burden at last of month equals initial burden
      ! of next month.  Also set TCH4 = 0 for next month.
      !=================================================================
      TCH4(:,:,:,1      ) = TCH4(:,:,:,2)
      TCH4(:,:,:,2:N_CH4) = 0d0

      ! Return to calling program
      END SUBROUTINE CH4_BUDGET

!------------------------------------------------------------------------------

      REAL*8 FUNCTION SUM_CH4( I1, I2, J1, J2, L1, L2, K1, K2, UPDOWN )
!
!******************************************************************************
!  Function SUM_CH4 sums a section of the TCH4 array bounded by the input
!  variables I1, I2, J1, J2, L1, L2, K1, K2.  SUM_CH4 is called by
!  module subroutine CH4_BUDGET. (jsw, bnd, bmy, 1/16/01)
!
!  Store the sources/sinks of CH4 in TCH4 in total molecules
!           ( 1) = Initial burden
!           ( 2) = Final burden
!  SINKS
!           ( 3) = Tropospheric CH4 sink by OH
!  SOURCES
!           ( 4) = Total Source
!           ( 5) = Industral
!           ( 6) = Agriculture
!           ( 7) = Biomass Burning
!           ( 8) = Termites
!           ( 9) = Wetland
!           (10) = Soil absorption
!           (11) = Interhemispheric Exchange (+ = northward)
!           (12) = ...
!
!  Levels            1 <= L <= LPAUSE(I,J) - 1 are tropospheric
!          LPAUSE(I,J) <= L <= LLPAR           are stratospheric (bmy, 4/17/00)
!
!  Arguments as Input:
!  ============================================================================
!  (1-2) I1, I2 (INTEGER) : Min and max longitude indices of TCH4 to sum
!  (3-4) J1, J2 (INTEGER) : Min and max latitude  indices of TCH4 to sum
!  (5-6) L1, L2 (INTEGER) : Min and max altitude  indices of TCH4 to sum
!  (7-8) K1, K2 (INTEGER) : Min and max tracer    indices of TCH4 to sum
!  (9  ) UPDOWN (INTEGER) : Sum in troposphere (=1) or in stratosphere (=0)
!
!  NOTES:
!  (1 ) Created by Bryan Duncan (1/99).  Adapted for CH4 chemistry by
!        James Wang (7/00).  Inserted into module "global_ch4_mod.f"
!        by Bob Yantosca. (bmy, 1/16/01)
!  (2 ) CH4_BUDGET is independent of "CMN_OH", "CMN_CO", and "CMN_CO_BUDGET".
!        (bmy, 1/16/01)
!  (3 ) Updated comments (jsw, bmy, 2/12/01)
!******************************************************************************
!
#     include "CMN_SIZE"       ! Size parameters
#     include "CMN"            ! LPAUSE

      ! Arguments
      INTEGER, INTENT(IN) :: I1, I2, J1, J2, L1, L2
      INTEGER, INTENT(IN) :: K1, K2, UPDOWN

      ! Local variables
      INTEGER             :: I, J, K, L, LPAUSE_MIN, LPAUSE_MAX

      !=================================================================
      ! SUM_CH4 begins here!
      !=================================================================

      ! Compute the minimum value of LPAUSE once for use in
      ! the DO-loops below (bmy, 4/18/00)
      LPAUSE_MIN = MINVAL( LPAUSE )
      LPAUSE_MAX = MAXVAL( LPAUSE )

      ! Initialize SUM_CH4
      SUM_CH4 = 0d0

      ! Test on UPDOWN
      IF ( UPDOWN == 1 ) THEN

         !=============================================================
         ! UPDOWN = 1: Sum up from the surface to the tropopause
         !=============================================================
         DO K = K1, K2
         DO L = L1, LPAUSE_MAX
         DO J = J1, J2
         DO I = I1, I2
            IF ( L < LPAUSE(I,J) ) THEN
               SUM_CH4 = SUM_CH4 + TCH4(I,J,L,K)
            ENDIF
         ENDDO
         ENDDO
         ENDDO
         ENDDO

      ELSE

         !=============================================================
         ! UPDOWN = 0: Sum up from the tropopause to the atm top
         !=============================================================
         DO K = K1,         K2
         DO L = LPAUSE_MIN, L2
         DO J = J1,         J2
         DO I = I1,         I2
            IF ( L >= LPAUSE(I,J) ) THEN
               SUM_CH4 = SUM_CH4 + TCH4(I,J,L,K)
            ENDIF
         ENDDO
         ENDDO
         ENDDO
         ENDDO
      ENDIF

      ! Return to calling program
      END FUNCTION SUM_CH4

!------------------------------------------------------------------------------

      SUBROUTINE CH4_DISTRIB(PREVCH4)
!
!******************************************************************************
!  Subroutine CH4_DISTRIB allocates the chemistry sink to different
!  emission tracers.
!  (ccc, 10/2/09)
!
!  Arguments as Input:
!  ============================================================================
!  PREVCH4(IIPAR, JJPAR, LLPAR) (REAL*8) : Store CH4 concentration before
!                                          chemistry
!
!******************************************************************************
!
      USE TRACER_MOD,    ONLY : STT, N_TRACERS
      USE ERROR_MOD,     ONLY : SAFE_DIV

      IMPLICIT NONE

#     include "CMN_SIZE"     ! Size parameters

      !Arguments
      REAL*8                 :: PREVCH4(IIPAR, JJPAR, LLPAR)

      !Local variables
      INTEGER                :: N, I, J, L

      !========================================================================
      ! CH4_DISTRIB begins here
      !========================================================================

      DO N=2,N_TRACERS

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            STT(I,J,L,N) = SAFE_DIV(STT(I,J,L,N),PREVCH4(I,J,L),0.d0)
     &                     * STT(I,J,L,1)
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO

      ENDDO

      ! Return to calling program
      END SUBROUTINE CH4_DISTRIB


!------------------------------------------------------------------------------

      FUNCTION GET_SCALE_GROUP( ) RESULT( CURRENT_GROUP )
!
!********************************************************************************
! Subroutine GET_SCALE_GROUP determines which predifined scaling index corresponds
! to the current time and location  (dkh, 12/02/04)
!
! NOTES
! (1 ) CURRENT_GROUP is currently only a function of TAU
! (2 ) Get rid of I,J as argument. (dkh, 03/28/05)
!
!********************************************************************************

      ! Reference to f90 modules
      USE TIME_MOD,       ONLY : GET_TAU, GET_TAUe, GET_TAUb, GET_MONTH
      USE ADJ_ARRAYS_MOD, ONLY : MMSCL

#     include "CMN_SIZE" ! Size stuff

      ! Arguments
      INTEGER      :: I, J

      ! Local Variables
      REAL*8       :: TOTAL_HR, CURRENT_HR, GROUP_LENGTH
      REAL*8       :: TAU, TAUe, TAUb

      ! Function variable
      INTEGER      :: CURRENT_GROUP
      LOGICAL, SAVE :: MONTHLY = .TRUE.
      INTEGER, SAVE :: MONTH_SAVE
      INTEGER, SAVE :: GROUP_SAVE
      LOGICAL, SAVE :: FIRST = .TRUE.

      !============================================================
      ! GET_SCALE_GROUP begins here!
      !============================================================

      ! Currently there is no spatial grouping

      ! Determine temporal grouping
      IF ( MMSCL == 1 ) THEN
         CURRENT_GROUP = 1
         RETURN
      ENDIF

      IF ( MONTHLY ) THEN
         IF (FIRST) THEN
            MONTH_SAVE = GET_MONTH()
            CURRENT_GROUP = MMSCL
            GROUP_SAVE = MMSCL
            FIRST = .FALSE.
         ENDIF
         IF ( MONTH_SAVE /= GET_MONTH() ) THEN
            MONTH_SAVE = GET_MONTH()
            GROUP_SAVE = GROUP_SAVE - 1
            CURRENT_GROUP = GROUP_SAVE
         ELSE
               CURRENT_GROUP = GROUP_SAVE
         ENDIF

      ELSE
         ! Retrieve time parameters
         TAUe       = GET_TAUe()
         TAUb       = GET_TAUb()
         TAU        = GET_TAU()
         TOTAL_HR   = TAUe - TAUb
         CURRENT_HR = TAU  - TAUb


         ! The last time step always belongs to the last group
         IF ( TAU == TAUe ) THEN
            CURRENT_GROUP = MMSCL
            RETURN
         ELSE

            ! Determine the length of each group
            GROUP_LENGTH = REAL( TOTAL_HR / MMSCL )

            ! Index is the current time divided by the group length, plus one
            CURRENT_GROUP = SNGL( CURRENT_HR / GROUP_LENGTH ) + 1

         ENDIF

      ENDIF

      END FUNCTION GET_SCALE_GROUP


!------------------------------------------------------------------------------

      SUBROUTINE INIT_GLOBAL_CH4
!
!******************************************************************************
!  Subroutine INIT_GLOBAL_CH4 allocates and zeroes module arrays.
!  (bmy, 1/16/01, 10/15/02)
!
!  NOTES:
!  (1 ) Now references ALLOC_ERR from "error_mod.f" (bmy, 10/15/02)
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD, ONLY : ALLOC_ERR

#     include "CMN_SIZE"
#     include "CMN_DIAG"

      ! Local variables
      INTEGER :: AS

      ALLOCATE( AVGOH( NSEAS, NCMSLATS, NCMSALTS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AVGOH' )
      AVGOH = 0d0

      ALLOCATE( BAIRDENS( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BAIRDENS' )
      BAIRDENS = 0d0

      ALLOCATE( BOH( IIPAR, JJPAR, LLPAR, 12 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BOH' )
      BOH = 0d0

      ALLOCATE( COPROD( JJPAR, LLPAR, 12 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'COPROD' )
      COPROD = 0d0

      ALLOCATE( PAVG( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'PAVG' )
      PAVG = 0d0

      ALLOCATE( TAVG( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'TAVG' )
      TAVG = 0d0

      ALLOCATE( CH4LOSS( IIPAR, JJPAR, LLPAR, 12 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'CH4LOSS' )
      CH4LOSS = 0d0

      ALLOCATE( TCH4( IIPAR, JJPAR, LLPAR, N_CH4 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'TCH4' )
      TCH4 = 0d0

      ALLOCATE( CH4_EMIS( IIPAR, JJPAR, PD58), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'CH4_EMIS' )
      CH4_EMIS = 0d0

      ! Return to calling program
      END SUBROUTINE INIT_GLOBAL_CH4

!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_GLOBAL_CH4
!
!******************************************************************************
!  Subroutine CLEANUP_GLOBAL_CH4 deallocates module arrays. (bmy, 1/16/01)
!******************************************************************************
!
      IF ( ALLOCATED( BAIRDENS  ) ) DEALLOCATE( BAIRDENS  )
      IF ( ALLOCATED( BOH       ) ) DEALLOCATE( BOH       )
      IF ( ALLOCATED( CH4LOSS   ) ) DEALLOCATE( CH4LOSS   )
      IF ( ALLOCATED( COPROD    ) ) DEALLOCATE( COPROD    )
      IF ( ALLOCATED( TCH4      ) ) DEALLOCATE( TCH4      )
      IF ( ALLOCATED( TAVG      ) ) DEALLOCATE( TAVG      )
      IF ( ALLOCATED( CH4_EMIS  ) ) DEALLOCATE( CH4_EMIS  )

      END SUBROUTINE CLEANUP_GLOBAL_CH4

!------------------------------------------------------------------------------

      END MODULE GLOBAL_CH4_MOD
