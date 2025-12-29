!$Id: strat_chem_mod.f,v 1.2 2012/07/13 20:09:14 nicolas Exp $
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: strat_chem_mod
!
! !DESCRIPTION: Module STRAT\_CHEM\_MOD contains variables and routines for
!  performing a simple linearized chemistry scheme for more realistic
!  upper boundary conditions. Archived 3D monthly climatological production
!  rates and loss frequencies are applied from the GMI combo model.
!
!  In the original schem code (schem.f), only the following species
!  were destroyed by photolysis in the stratosphere:
!    PAN, H2O2, ACET, MEK, ALD2, RCHO, MVK, MACR, R4N2, CH2O, N2O5, HNO4, MP
!  and by reaction with OH for:
!    ALK4, ISOP, H2O2, ACET, MEK, ALD2, RCHO, MVK, MACR, PMN, R4N2,
!    PRPE, C3H8, CH2O, C2H6, HNO4, MP
!
!  The updated code includes at least all of these, and many more. The code
!  is flexible enough to automatically apply the rate to any new tracers
!  for future simulations that share the name in tracer_mod with the
!  GMI name.  (See Documentation).
!
!\\
!\\
! !INTERFACE:
!
      MODULE STRAT_CHEM_MOD
!
! !USES:
!
      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: DO_STRAT_CHEM
      PUBLIC  :: CLEANUP_STRAT_CHEM

      ! hml
      PUBLIC  :: PROD_0
      PUBLIC  :: LOSS_0
      PUBLIC  :: PROD
      PUBLIC  :: LOSS
      PUBLIC  :: DTCHEM
      PUBLIC  :: NSCHEM
      PUBLIC  :: GET_RATES
      PUBLIC  :: GET_RATES_INTERP
      PUBLIC  :: Strat_TrID_GC

!
! !PRIVATE MEMBER FUNCTIONS:
!
      PRIVATE :: INIT_STRAT_CHEM

!
! !ADJOINT GROUP:
!
!
! !PUBLIC DATA MEMBERS:
!
! !REMARKS:
!
!  References:
!  ============================================================================
!  (1 )
! !REVISION HISTORY:
!  1 Feb 2011 - L. Murray - Initial version
!  22 Oct 2011 - H.-M. Lee - Modified to implement in adjoint.
!                Now we can calculte strat prod and loss sensitivity. adj32_025
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!

      ! Scalars
      REAL*8               :: dTchem          ! chemistry time step [s]
      INTEGER, PARAMETER   :: NTR_GMI  = 120  ! Number of species from GMI model
      INTEGER              :: NSCHEM          ! Number of species upon which to
                                              ! apply P's & k's in GEOS-Chem

      ! Arrays
      REAL*8,  ALLOCATABLE :: PROD(:,:,:,:)
      REAL*8,  ALLOCATABLE :: LOSS(:,:,:,:)
      REAL*8,  ALLOCATABLE :: STRAT_OH(:,:,:) ! Monthly mean OH [v/v]
      INTEGER, SAVE        :: ncID_strat_rates

      ! For adjoint
      REAL*8,  ALLOCATABLE :: PROD_0(:,:,:,:)
      REAL*8,  ALLOCATABLE :: LOSS_0(:,:,:,:)

      CHARACTER(LEN=16)    :: GMI_TrName(NTR_GMI)     ! Tracer names in GMI
      INTEGER              :: Strat_TrID_GC(NTR_GMI)  ! Maps 1:NSCHEM to STT index
      INTEGER              :: Strat_TrID_GMI(NTR_GMI) ! Maps 1:NSCHEM to GMI index
                         ! (At most NTR_GMI species could overlap between G-C & GMI)


      ! Variables used to calculate the strat-trop exchange flux
      REAL*8               :: TauInit             ! Initial time
      INTEGER              :: NymdInit, NhmsInit  ! Initial date
      REAL*8               :: TpauseL_Cnt         ! Tropopause counter
      REAL*8, ALLOCATABLE  :: TpauseL(:,:)        ! Tropopause level aggregator
      REAL*8, ALLOCATABLE  :: MInit(:,:,:,:)      ! Init. atm. state for STE period
      REAL*8, ALLOCATABLE  :: Before(:,:,:)       ! Init. atm. state each chem. dt
      REAL*8, ALLOCATABLE  :: SChem_Tend(:,:,:,:) ! Stratospheric chemical tendency
                                                  !   (total P - L) [kg period-1]

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement
      !=================================================================
      CONTAINS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: DO_STRAT_CHEM
!
! !DESCRIPTION: Function DO\_STRAT\_CHEM is the driver routine for computing
!     the simple linearized stratospheric chemistry scheme for a host of species
!     whose prod/loss rates were determined from the GMI combo model. Ozone is
!     treated using either Linoz or Synoz.
!
! !INTERFACE:
!
      SUBROUTINE DO_STRAT_CHEM
!
! !USES:
!
      USE DAO_MOD,        ONLY : AD, CONVERT_UNITS
      USE ERROR_MOD,      ONLY : DEBUG_MSG, GEOS_CHEM_STOP
      USE LOGICAL_MOD,    ONLY : LLINOZ, LPRT
      USE LINOZ_MOD,      ONLY : DO_LINOZ
      USE TIME_MOD,       ONLY : GET_MONTH, TIMESTAMP_STRING
      USE TRACER_MOD,     ONLY : ITS_A_FULLCHEM_SIM, ITS_A_TAGOX_SIM
      USE TRACER_MOD,     ONLY : N_TRACERS, STT, TCVV, TRACER_MW_KG
      USE TRACERID_MOD,   ONLY : IDTOX
      USE TROPOPAUSE_MOD, ONLY : GET_MIN_TPAUSE_LEVEL, GET_TPAUSE_LEVEL
      USE TROPOPAUSE_MOD, ONLY : ITS_IN_THE_TROP
      USE NETCDF_UTIL_MOD

      ! adj_group (hml, 07/25/11)
      USE ADJ_ARRAYS_MOD, ONLY : PROD_SF, LOSS_SF
      USE ADJ_ARRAYS_MOD, ONLY : ID_LOSS
      USE ADJ_ARRAYS_MOD, ONLY : IFD, JFD, LFD, DO_CHK_FILE
      USE ADJ_ARRAYS_MOD, ONLY : NSTPL
      USE TRACER_MOD,     ONLY : STT_STRAT_TMP
      USE LOGICAL_ADJ_MOD,ONLY : LADJ
      USE LOGICAL_ADJ_MOD,ONLY : LADJ_STRAT
      USE CHECKPOINT_MOD, ONLY : MAKE_BEFSTRAT_CHKFILE
      USE TIME_MOD,       ONLY : GET_NHMS
      USE TIME_MOD,       ONLY : GET_NYMD
      USE TIME_MOD,       ONLY : GET_TAU
      USE UPBDFLX_MOD,    ONLY : UPBDFLX_O3, INIT_UPBDFLX


#     include "CMN_SIZE"
!
! !REMARKS:
!
! !REVISION HISTORY:
!  1 Feb 2011 - L. Murray - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE             :: FIRST = .TRUE.
      INTEGER, SAVE             :: LASTMONTH = -999
      INTEGER, SAVE             :: LASTSEASON = -1
      INTEGER                   :: I, J, L, N, LMIN
      INTEGER                   :: IORD, JORD, KORD
      INTEGER                   :: NHMS
      INTEGER                   :: NYMD
      INTEGER                   :: NN, NS, NSL
      REAL*8                    :: TAU
      REAL*8                    :: dt, P, k, M0
      REAL*8                    :: STT0(IIPAR,JJPAR,LLPAR,N_TRACERS)
      CHARACTER(LEN=16)         :: STAMP


      !===============================
      ! DO_STRAT_CHEM begins here!
      !===============================

      STAMP = TIMESTAMP_STRING()
      WRITE( 6, 10 ) STAMP
 10   FORMAT( '    - DO_STRAT_CHEM: Linearized strat chemistry at ', a )

      IF ( FIRST ) THEN

         ! Allocate all module arrays
         CALL INIT_STRAT_CHEM

#if    defined( GEOS_3 )
         ! Initialize some Synoz variables
         IF ( .NOT. ( LLINOZ ) ) THEN
            CALL GET_ORD( IORD, JORD, KORD )
            CALL INIT_UPBDFLX( IORD, JORD, KORD )
         ENDIF
#endif

      ENDIF

      ! Get the minimum level extent of the tropopause
      LMIN = GET_MIN_TPAUSE_LEVEL()

      IF ( GET_MONTH() /= LASTMONTH ) THEN

         IF ( LPRT ) CALL DEBUG_MSG( '### STRAT_CHEM: at GET_RATES' )

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

         ! Save month for next iteration
         LASTMONTH = GET_MONTH()
      ENDIF

      ! Set first-time flag to false
      FIRST = .FALSE.

      IF ( LPRT ) CALL DEBUG_MSG( '### STRAT_CHEM: at DO_STRAT_CHEM' )

      WRITE(6,*) '-----------------------------------------------------'
      write(6,*) '    Doing stratospheric chemistry (STRAT_CHEM_MOD)   '
      WRITE(6,*) '-----------------------------------------------------'


      !================================================================
      ! Full chemistry simulations
      !================================================================
      IF ( ITS_A_FULLCHEM_SIM() ) THEN

         !! Advance counter for number of times we've sampled the tropopause level
         !TpauseL_CNT = TpauseL_CNT + 1d0

         !=============================================================
         ! Do chemical production and loss for non-ozone species for
         ! which we have explicit prod/loss rates from GMI
         !=============================================================

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N, k, P, dt, M0, NS, NN )
!$OMP+SCHEDULE( DYNAMIC )
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! Add to tropopause level aggregator for later determining STE flux
            TpauseL(I,J) = TpauseL(I,J) + GET_TPAUSE_LEVEL(I,J)

            DO L = LMIN, LLPAR

               IF ( ITS_IN_THE_TROP( I, J, L ) ) CYCLE

               DO N = 1, NSCHEM ! Tracer index of active strat chem species
                  NN = Strat_TrID_GC(N) ! Tracer index in STT

                  ! Skip Ox; we'll always use either Linoz or Synoz
                  ! Adj use GMI rate for Ox if LINOZ is off (hml, 10/31/11)
                  !IF ( ITS_A_FULLCHEM_SIM() .and. NN .eq. IDTOx ) CYCLE
                  IF ( ITS_A_FULLCHEM_SIM() .and. (NN .eq. IDTOx) .and.
     &                 LLINOZ ) CYCLE

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

                  ! Check point values of STT
                  STT_STRAT_TMP(I,J,L,NN) = STT(I,J,L,NN)

                  dt = DTCHEM                              ! timestep [s]
                  k  = LOSS(I,J,L,N)                       ! loss freq [s-1]
                  P  = PROD(I,J,L,N) * AD(I,J,L) / TCVV(NN)! production term [kg s-1]
                  M0 = STT(I,J,L,NN)                       ! initial mass [kg]

                  ! debug test
                  !IF ( I == IFD .and. J == JFD .and. L == LFD ) THEN
                  !   print*, NN,' STRAT TEST fwd: k = ', k
                  !   print*, NN,' STRAT TEST fwd: P = ', P
                  !   print*, NN,' STRAT TEST fwd: M0= ', M0
                  !ENDIF

                  ! No prod or loss at all
                  IF ( k .eq. 0d0 .and. P .eq. 0d0 ) CYCLE

                  ! Simple analytic solution to dM/dt = P - kM over [0,dt]
                  IF ( k .gt. 0d0 ) THEN
                     STT(I,J,L,NN) = M0 * exp(-k*dt) + (P/k)
     &                             * (1d0-exp(-k*dt))
                  ELSE
                     STT(I,J,L,NN) = M0 + P*dt
                  ENDIF

                   ! Aggregate chemical tendency [kg box-1]
                   SCHEM_TEND(I,J,L,NN) = SCHEM_TEND(I,J,L,NN)
     &                                  + ( STT(I,J,L,NN) - M0 )

               ENDDO ! N
            ENDDO ! L
         ENDDO ! I
         ENDDO ! J
!$OMP END PARALLEL DO

         ! Make check point file
         IF ( DO_CHK_FILE() ) THEN
            NHMS     = GET_NHMS()
            NYMD     = GET_NYMD()
            TAU      = GET_TAU()
            CALL MAKE_BEFSTRAT_CHKFILE( NYMD, NHMS, TAU )
         ENDIF

         !===================================
         ! Ozone
         !===================================

         ! Make note of inital state for determining tendency later
         BEFORE = STT(:,:,:,IDTOX )
         ! Put ozone in v/v
         STT(:,:,:,IDTOX ) = STT(:,:,:,IDTOX) * TCVV( IDTOX ) / AD

         IF ( LLINOZ ) THEN
            CALL DO_LINOZ       ! Linoz
         ELSE IF ( .not. LADJ ) THEN
            ! Must use Linoz or strat chem Ox fluxes for the adjoint
            CALL UPBDFLX_O3     ! Synoz
         ENDIF

         ! Now move unit conversion into LINOZ (hml, 11/06/11)
         ! Put ozone back to kg
         STT(:,:,:,IDTOX) = STT(:,:,:,IDTOX) * AD / TCVV( IDTOX )

         ! Put tendency into diagnostic array [kg box-1]
         SCHEM_TEND(:,:,:,IDTOX) = SCHEM_TEND(:,:,:,IDTOX)
     &                           + ( STT(:,:,:,IDTOX) - BEFORE )


      !======================================================================
      ! Tagged Ox simulation
      !======================================================================

      ELSE IF ( ITS_A_TAGOX_SIM() ) THEN

         ! Intial conditions
         STT0(:,:,:,:) = STT(:,:,:,:)

         CALL CONVERT_UNITS( 1, N_TRACERS, TCVV, AD, STT ) ! kg -> v/v

         IF ( LLINOZ ) THEN
            CALL DO_LINOZ       ! Linoz
         ELSE IF ( .not. LADJ ) THEN
            ! must use Linoz or strat chem Ox fluxes for the adjoint
            CALL UPBDFLX_O3     ! Synoz
         ENDIF

         CALL CONVERT_UNITS( 2, N_TRACERS, TCVV, AD, STT ) ! v/v -> kg

         ! Add to tropopause level aggregator for later determining STE flux
         TpauseL_CNT = TpauseL_CNT + 1d0
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
         DO I=1,IIPAR
            DO J=1,JJPAR
               TpauseL(I,J) = TpauseL(I,J) + GET_TPAUSE_LEVEL(I,J)
            ENDDO
         ENDDO
!$OMP END PARALLEL DO

         ! Aggregate chemical tendency [kg box-1]
         DO N=1,NSCHEM
            NN = Strat_TrID_GC(N)
               SCHEM_TEND(:,:,:,NN) = SCHEM_TEND(:,:,:,NN)
     &                             + ( STT(:,:,:,NN) - STT0(:,:,:,NN) )
         ENDDO

      ELSE

         ! The code will need to be modified for other tagged simulations
         ! (e.g., CO). Simulations like CH4, CO2 with standard tracer names
         ! should probably just work as is with the full chemistry code above,
         ! but would need to be tested.
         WRITE( 6, * ) 'Strat chemistry needs to be activated for ' //
     &         'your simulation type.'
         !WRITE( 6, * ) 'Please see GeosCore/strat_chem_mod.F90' // &
         WRITE( 6, * ) 'Please see new/strat_chem_mod.f' //
     &         'or disable in input.geos'
         CALL GEOS_CHEM_STOP

      ENDIF

      END SUBROUTINE DO_STRAT_CHEM
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GET_RATES
!
! !DESCRIPTION: Function GET\_RATES reads from disk the chemical production
!  and loss rates for the species of interest
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GET_RATES( THISMONTH )
!
! !USES:
!
      USE BPCH2_MOD,       ONLY : GET_NAME_EXT, GET_RES_EXT
      USE BPCH2_MOD,       ONLY : GET_TAU0, READ_BPCH2
      USE DIRECTORY_MOD,   ONLY : DATA_DIR
      USE TIME_MOD,        ONLY : GET_MONTH
      USE LOGICAL_MOD,     ONLY : LLINOZ
      USE TRACER_MOD,      ONLY : N_TRACERS, TRACER_NAME, TRACER_COEFF
      USE TRANSFER_MOD,    ONLY : TRANSFER_3D
      USE NETCDF_UTIL_MOD, ONLY : NCDF_GET_VAR
      USE NETCDF_UTIL_MOD, ONLY : NCDF_GET_VARID
      USE NETCDF_UTIL_MOD, ONLY : NCDF_OPEN_FOR_READ
      USE NETCDF_UTIL_MOD, ONLY : NCDF_CLOSE


      IMPLICIT NONE

#     include "CMN_SIZE"

!
! !INPUT PARAMETERS:
!
      ! Arguments
      INTEGER,INTENT(IN) :: THISMONTH
!
! !REVISION HISTORY:
!  1 Feb 2011 - L. Murray - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      CHARACTER(LEN=255) :: FILENAME
      CHARACTER(LEN=6)   :: SPNAME( NTR_GMI )
      REAL*4             :: ARRAY( IIPAR, JJPAR, LGLOB ) ! Full vertical res
      REAL*8             :: ARRAY2( IIPAR, JJPAR, LLPAR )! Actual vertical res
      REAL*8             :: XTAU
      INTEGER            :: N, M, S, F, NN, fileID
      INTEGER            :: prodID, lossID
      INTEGER            :: ohID


      !=================================================================
      ! GET_RATES begins here
      !=================================================================

      ! Initialize arrays
      LOSS = 0d0
      PROD = 0d0

      WRITE(6, 11  )
     &       '       - Getting new strat prod/loss rates for month: ',
     &       THISMONTH
11     FORMAT( a, I2.2 )

      M = THISMONTH

      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ! Get stratospheric OH mixing ratio [v/v]
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      FILENAME = 'strat_chem_201206/gmi.clim.OH.' //
     &                 GET_NAME_EXT() // '.' // GET_RES_EXT() // '.nc'
      FILENAME = TRIM( DATA_DIR ) // TRIM( FILENAME )
      WRITE(6,'(a)') '       => Reading from file: ' // trim(filename)
      call ncdf_open_for_read( fileID, TRIM( FILENAME ) )

      ohID = ncdf_get_varid( fileID, 'species' )

      call ncdf_get_var( fileID, ohID, array,
     &                        (/     1,     1,     1,  m /),   ! Start
     &                        (/ iipar, jjpar, lglob,  1 /)  ) ! Count
      call ncdf_close( fileID )

      ! Cast from REAL*4 to REAL*8 and resize to 1:LLPAR
      call transfer_3D( array, array2 )

      STRAT_OH(:,:,:) = ARRAY2

      DO N=1,NSCHEM
         NN = Strat_TrID_GMI(N)

         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         ! Open individual species file
         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         FILENAME = 'strat_chem_201206/gmi.clim.' //
     &      TRIM( GMI_TrName(NN) ) // '.' //
     &      GET_NAME_EXT() // '.' // GET_RES_EXT() // '.nc'
         FILENAME = TRIM( DATA_DIR ) // TRIM( FILENAME )
         WRITE(6,'(a)') '      => Reading from file: ' // trim(filename)
         call ncdf_open_for_read( fileID, TRIM( FILENAME ) )

         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         ! Read production rate [v/v/s]
         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         ! Get the variable IDs for the species, prod and loss rates
         prodID = ncdf_get_varid( fileID, 'prod' )

         call ncdf_get_var( fileID, prodID, array,
     &                              (/     1,     1,     1,  m  /),   ! Start
     &                              (/ iipar, jjpar, lglob,  1  /)  ) ! Count

         ! Cast from REAL*4 to REAL*8 and resize to 1:LLPAR
         call transfer_3D( array, array2 )

         PROD(:,:,:,N) = ARRAY2

         ! Save rates from file to respective arrays (hml, 09/15/11)
         PROD_0(:,:,:,N) = PROD(:,:,:,N)

         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         ! Read loss frequencies [s^-1]
         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         lossID = ncdf_get_varid( fileID, 'loss' )

         call ncdf_get_var( fileID, lossID, array,
     &         (/     1,     1,     1,  m  /),   ! Start
     &         (/ iipar, jjpar, lglob,  1  /)  ) ! Count

         ! Cast from REAL*4 to REAL*8 and resize to 1:LLPAR
         call transfer_3D( array, array2 )

         LOSS(:,:,:,N) = ARRAY2

         ! Save rates from file to respective arrays (hml, 09/15/11)
         LOSS_0(:,:,:,N) = LOSS(:,:,:,N)

         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         ! Close species file
         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         call ncdf_close( fileID )

      ENDDO

      END SUBROUTINE GET_RATES
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GET_RATES_INTERP
!
! !DESCRIPTION: Function GET\_RATES\_INTERP reads from disk the chemical
! production and loss rates for the species of interest to resolutions finer
! than 2 x 2.5 (e.g., nested simluations) via simple nearest-neighbor mapping.
!\\
!\\
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GET_RATES_INTERP( THISMONTH )
!
! !USES:
!
      USE BPCH2_MOD,       ONLY : GET_NAME_EXT, GET_RES_EXT
      USE DIRECTORY_MOD,   ONLY : DATA_DIR
      USE GRID_MOD,        ONLY : GET_YMID, GET_XMID
      USE LOGICAL_MOD,     ONLY : LLINOZ
      USE TRACER_MOD,      ONLY : N_TRACERS, TRACER_NAME, TRACER_COEFF
      USE TRANSFER_MOD,    ONLY : TRANSFER_3D
      USE ADJ_ARRAYS_MOD,  ONLY : NSTPL
      USE NETCDF_UTIL_MOD, ONLY : NCDF_GET_VAR
      USE NETCDF_UTIL_MOD, ONLY : NCDF_GET_VARID
      USE NETCDF_UTIL_MOD, ONLY : NCDF_OPEN_FOR_READ
      USE NETCDF_UTIL_MOD, ONLY : NCDF_CLOSE

#     include "define.h"
#     include "CMN_SIZE"
!
! !INPUT PARAMETERS:
!
      ! Arguments
      INTEGER,INTENT(IN) :: THISMONTH
!
! !REVISION HISTORY:
!  1 Feb 2011 - L. Murray - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      CHARACTER(LEN=255) :: FILENAME
      CHARACTER(LEN=6)   :: SPNAME( NTR_GMI )
      REAL*4             :: ARRAY( IIPAR, JJPAR, LGLOB )
      REAL*8             :: ARRAY2( IIPAR, JJPAR, LLPAR )
      INTEGER            :: N, M, S, F
      INTEGER            :: NN
      INTEGER            :: ohID
      INTEGER            :: prodID, lossID
      INTEGER            :: lat_varID, lon_varID

      REAL*4             :: XMID_COARSE(144), YMID_COARSE(91)
      INTEGER            :: I_f2c(IIPAR), J_f2c(JJPAR) ! f2c = fine to coar map'ng
      INTEGER            :: I, J, fileID
      INTEGER            :: II(1), JJ(1)
      REAL*4             :: COLUMN( LGLOB )



      !=================================================================
      ! GET_RATES_INTERP begins here
      !=================================================================

      ! In the original schem code, the following species were destroyed
      ! by photolysis in the stratosphere:
      !  PAN, H2O2, ACET, MEK, ALD2, RCHO, MVK, MACR, R4N2, CH2O,
      !  N2O5, HNO4, MP
      ! And by reaction with OH for:
      !  ALK4, ISOP, H2O2, ACET, MEK, ALD2, RCHO, MVK, MACR, PMN, R4N2,
      !  PRPE, C3H8, CH2O, C2H6, HNO4, MP
      ! The updated code includes at least all of these, and several more.

      ! Initialize arrays
      LOSS = 0d0
      PROD = 0d0

      ! first read in the OH file so that we can get the lat and long
      ! values to populate XMID_COARSE and YMID_COARSE
      ! Path to input data, use 2 x 2.5 file
#if defined( GEOS_FP )
      FILENAME = 'strat_chem_201206/gmi.clim.OH.geos5' //
     &           '.2x25.nc'
#else
      FILENAME = 'strat_chem_201206/gmi.clim.OH.' // GET_NAME_EXT() //
     &           '.2x25.nc'
#endif
      !FILENAME = TRIM( DATA_DIR_1x1 ) // TRIM( FILENAME )
      !FILENAME = TRIM( DATA_DIR ) // TRIM('../GEOS_2x25/') //
      FILENAME = TRIM( DATA_DIR ) // TRIM('../GEOS_2x2.5/') //
     &           TRIM( FILENAME )

      WRITE(6, 11  )
     &       '       - Getting new strat prod/loss rates for month: ',
     &       THISMONTH
11     FORMAT( a, I2.2 )

      ! Open the netCDF file containing the rates
      WRITE(6,'(a)')
     &     '         => Interpolate to resolution from file: '
     &    // trim(filename)
      call ncdf_open_for_read( fileID, TRIM( filename ) )

      ! Get the lat and lon centers of the 2x2.5 GMI climatology
      ! WARNING MAKE 2x25 after testing
      !call NcRd( XMID_COARSE, fileID, 'longitude', (/1/),  (/144/) )
      !call NcRd( YMID_COARSE, fileID, 'latitude',  (/1/),  (/91/) )
      lat_varID = ncdf_get_varid( fileID, 'latitude' )
      lon_varID = ncdf_get_varid( fileID, 'longitude' )

      !call NcRd( XMID_COARSE, fileID, 'longitude', (/1/),  (/144/) )
      !call NcRd( YMID_COARSE, fileID, 'latitude',  (/1/),  (/91/) )
      call ncdf_get_var( fileID, lon_varID, XMID_COARSE, (/1/), (/144/))
      call ncdf_get_var( fileID, lat_varID, YMID_COARSE, (/1/), (/91/) )

      ! For each fine grid index, determine the closest coarse (2x2.5) index
      ! Note: This doesn't do anything special for the date line, and may
      ! therefore not pick the exact closest if it is on the other side.
      ! Note: CMN_SIZE_MOD claims in its comments that IIPAR < IGLOB, but
      ! in actuality, IIPAR = IGLOB and JJPAR = JGLOB, the dimensions of the nested
      ! region.
      DO I=1,IGLOB
         II = MINLOC( ABS( GET_XMID(I) - XMID_COARSE ) )
         I_f2c(I) = II(1)
         !print*,'I:',I,'->',II(1)
      ENDDO
      DO J=1,JGLOB
         JJ = MINLOC( ABS( GET_YMID(J) - YMID_COARSE ) )
         J_f2c(J) = JJ(1)
         !print*,'J:',J,'->',JJ(1)
      ENDDO

!      call ncdf_close( fileID )

      M = THISMONTH

      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      !Get Stratospheric OH concentrations [v/v]
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      ohID = ncdf_get_varid( fileID, 'species' )

      DO I=1,IGLOB
         DO J=1,JGLOB

            call ncdf_get_var( fileID, ohID, column,
     &                 (/ I_f2c(I), J_f2c(J),     1, m /),   ! Start
     &                 (/        1,        1, lglob, 1 /)  ) ! Count
            array( I, J, : ) = column

         ENDDO
      ENDDO
      call ncdf_close( fileID )
      call transfer_3D( array, array2 )
      STRAT_OH(:,:,:) = ARRAY2


      DO N=1,NSCHEM
         NN = Strat_TrID_GMI(N)

         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         ! Open individual species file
         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         ! Path to input data, use 2 x 2.5 file
         ! (lzh,02/01/2015) add geosfp - use geos-5
#if defined( GEOS_FP )
         FILENAME = 'strat_chem_201206/gmi.clim.' //
     &         TRIM( GMI_TrName(NN) ) // '.geos5' //
     &         '.2x25.nc'
#else
         FILENAME = 'strat_chem_201206/gmi.clim.' //
     &         TRIM( GMI_TrName(NN) ) // '.' // GET_NAME_EXT() //
     &         '.2x25.nc'
#endif
         !FILENAME = TRIM( DATA_DIR ) // TRIM('../GEOS_2x25/') //
         FILENAME = TRIM( DATA_DIR ) // TRIM('../GEOS_2x2.5/') //
     &              TRIM( FILENAME )

         WRITE(6,'(a)')
     &         '         => Interpolate to resolution from file: ' //
     &         trim(filename)

         call ncdf_open_for_read( fileID, TRIM( FILENAME ) )

         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         ! Read production rate [v/v/s]
         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         array = 0.0

         prodID = ncdf_get_varid( fileID, 'prod' )

         DO I=1,IGLOB
            DO J=1,JGLOB

               call ncdf_get_var( fileID, prodID, column,
     &                          (/ I_f2c(I), J_f2c(J),     1,  m /),   ! Start
     &                          (/        1,        1, lglob,  1 /)  ) ! Count
               array( I, J, : ) = column

            ENDDO
         ENDDO

         ! Cast from REAL*4 to REAL*8 and resize to 1:LLPAR if necessary
         call transfer_3D( array, array2 )

         PROD(:,:,:,N) = ARRAY2

         ! Save rates from file to respective arrays (hml, 09/15/11)
         PROD_0(:,:,:,N) = PROD(:,:,:,N)

         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         ! Read loss frequencies [s^-1]
         !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         array = 0.0

         lossID = ncdf_get_varid( fileID, 'loss' )

         DO I=1,IGLOB
            DO J=1,JGLOB

               call ncdf_get_var( fileID, lossID, column,
     &                          (/ I_f2c(I), J_f2c(J),     1,  m /),   ! Start
     &                          (/        1,        1, lglob,  1 /)  ) ! Count
               array( I, J, : ) = column

            ENDDO
         ENDDO

         ! Cast from REAL*4 to REAL*8 and resize to 1:LLPAR if necessary
         call transfer_3D( array, array2 )

         LOSS(:,:,:,N) = ARRAY2

         ! Save rates from file to respective arrays (hml, 09/15/11)
         LOSS_0(:,:,:,N) = LOSS(:,:,:,N)

         call ncdf_close( fileID )

      ENDDO

      !call ncdf_close( ncID_strat_rates )

      END SUBROUTINE GET_RATES_INTERP
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: calc_ste
!
! !DESCRIPTION: Subroutine CALC\_STE estimates what the stratosphere-to-
!               troposphere exchange flux must have been since the last time
!               it was reset
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CALC_STE
!
! !USES:
!
      USE TRACER_MOD, ONLY : STT, TRACER_MW_KG, N_TRACERS, TRACER_NAME
      USE TIME_MOD,   ONLY : GET_TAU, GET_NYMD, GET_NHMS, EXPAND_DATE

      IMPLICIT NONE

#include "define.h"
#include "CMN_SIZE"

!
! !REVISION HISTORY:
!  28 Apr 2012 - L. Murray - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC

      REAL*8             :: M1(IIPAR,JJPAR,LLPAR), M2(IIPAR,JJPAR,LLPAR)
      REAL*8             :: dStrat, STE, Tend, tauEnd, dt
      INTEGER            :: N, I, J, L, NN,  LTP(IIPAR,JJPAR)
      CHARACTER(LEN=256) :: dateStart, dateEnd

      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ! By simple mass balance, dStrat/dt = P - L - STE,
      ! where STE is the net stratosphere-to-troposphere mass exchange.
      !
      ! Therefore, we estimate STE as
      !   STE = (P-L) - dStrat/dt
      !
      ! As the tropopause is dynamic, we use the mean tropopause level during
      ! the period for determining initial and end stratospheric masses.
      ! (ltm, 04/28/2012)
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if defined( NESTED_NA ) || defined( NESTED_CH ) || defined( NESTED_EU )
      ! This method only works for a global domain.
      ! It could be modified for nested domains if the total mass flux across the
      ! boundaries during the period is taken into account.
      RETURN
#endif

      ! Determine mean tropopause level for the period
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I,  J  )
      DO I = 1,IIPAR
         DO J = 1,JJPAR
            LTP(I,J) = NINT( TPauseL(I,J) / TPauseL_Cnt )
         ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Period over which STE is being determined [a]
      tauEnd = GET_TAU() ! [h]
      dt = ( tauEnd - tauInit ) / 24d0 / 365.25d0

      dateStart = 'YYYY-MM-DD hh:mm'
      CALL EXPAND_DATE(dateStart,NymdInit,NhmsInit)
      dateEnd = 'YYYY-MM-DD hh:mm'
      CALL EXPAND_DATE(dateEnd,GET_NYMD(),GET_NHMS())

      ! Print to output
      WRITE( 6, * ) ''
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )
      WRITE( 6, '(a)' ) '  Strat-Trop Exchange'
      WRITE( 6, '(a)' ) REPEAT( '-', 79 )
      WRITE( 6, '(a)' )
     &     '  Global stratosphere-to-troposphere fluxes estimated over'
      WRITE( 6, 100 ) TRIM(dateStart), TRIM(dateEnd)
      WRITE( 6, * ) ''
      WRITE( 6, 110 ) 'Species','[moles a-1]','* [g/mol]','= [Tg a-1]'

  100 FORMAT( 2x,a16,' to ',a16 )
  110 FORMAT( 2x,a8,':',4x,a11  ,4x,a9  ,4x,  a11 )

      ! Loop through each species
      DO N=1,N_TRACERS

         ! Populate before (M1) and after (M2) state for the species [kg]
         M1 = MInit(:,:,:,N)
         M2 =   STT(:,:,:,N)

         ! Zero out tropopshere and determine total change in the stratospheric
         ! burden of species N (dStrat) [kg]
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I,  J  )
         DO I=1,IIPAR
            DO J=1,JJPAR
               M2(I,J,1:LTP(I,J)) = 0d0
               M1(I,J,1:LTP(I,J)) = 0d0
            ENDDO
         ENDDO
!$OMP END PARALLEL DO
         dStrat   = SUM(M2)-SUM(M1)

         ! The total chemical tendency (P-L) over the period for species N [kg]
         Tend   = SUM(Schem_tend(:,:,:,N))

         ! Calculate flux as STE = (P-L) - dStrat/dt
         STE = (Tend-dStrat)/dt ! [kg a-1]

         ! Print to standard output
         WRITE(6,120) TRIM(TRACER_NAME(N)),
     &         STE/TRACER_MW_KG(N),            ! mol/a-1
     &         TRACER_MW_KG(N)*1d3,            ! g/mol
     &         STE*1d-9                        ! Tg a-1

      ENDDO

  120 FORMAT( 2x,a8,':',4x,e11.3,4x,f9.1,4x,f11.4 )

      WRITE( 6, * ) ''
      WRITE( 6, '(a)'   ) REPEAT( '=', 79 )
      WRITE( 6, * ) ''

      ! Reset variables for next STE period
      NymdInit             = GET_NYMD()
      NhmsInit             = GET_NHMS()
      TauInit              = GET_TAU()
      TPauseL_Cnt          = 0d0
      TPauseL(:,:)         = 0d0
      SChem_tend(:,:,:,:)  = 0d0
      MInit(:,:,:,:)       = STT(:,:,:,:)

      END SUBROUTINE CALC_STE
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_strat_chem
!
! !DESCRIPTION: Subroutine INIT\_STRAT\_CHEM allocates all module arrays.
!  It also opens the necessary rate files.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_STRAT_CHEM
!
! !USES:
!
      USE ERROR_MOD,       ONLY : ALLOC_ERR
      USE LOGICAL_MOD,     ONLY : LLINOZ
      USE TRACER_MOD,      ONLY : ITS_A_FULLCHEM_SIM, ITS_A_TAGOX_SIM
      USE TRACER_MOD,      ONLY : N_TRACERS, TRACER_NAME, TRACER_COEFF
      USE TRACER_MOD,      ONLY : TRACER_COEFF, STT
      USE BPCH2_MOD,       ONLY : GET_NAME_EXT, GET_RES_EXT
      USE DIRECTORY_MOD,   ONLY : DATA_DIR
      USE TIME_MOD,        ONLY : EXPAND_DATE
      USE TIME_MOD,        ONLY : GET_TAU, GET_NYMD, GET_NHMS
      USE TIME_MOD,        ONLY : GET_TS_CHEM
      USE LOGICAL_ADJ_MOD, ONLY : LADJ
      USE UPBDFLX_MOD,     ONLY : UPBDFLX_O3, INIT_UPBDFLX

#     include "CMN_SIZE"
!
! !REVISION HISTORY:
!  1 Feb 2011 - L. Murray - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!

      CHARACTER(LEN=16)     :: sname
      INTEGER               :: AS, N, NN

      CHARACTER(LEN=255)    :: FILENAME, FILENAMEOUT
      CHARACTER(LEN=6)      :: SPNAME( NTR_GMI )
      INTEGER               :: spname_varID

      !=================================================================
      ! INIT_STRAT_CHEM begins here!
      !=================================================================

      ! Initialize counters, initial times, mapping arrays
      TpauseL_Cnt       = 0.
      NSCHEM            = 0
      TauInit           = GET_TAU()
      NymdInit          = GET_NYMD()
      NhmsInit          = GET_NHMS()
      strat_trID_GC(:)  = 0
      strat_trID_GMI(:) = 0

      ! Initialize timestep for chemistry [s]
      dTchem = GET_TS_CHEM() * 60d0

      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ! Determine the mapping for the GMI to the GC variables based on
      ! tracer name, which only needs to be done once per model run.
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      ! List of available tracers with archived monthly climatological
      ! production rates, loss frequencies, and mixing ratios from the
      ! GMI Combo model (tracer names here are as used in GMI).
      GMI_TrName = (/
     &     'A3O2',     'ACET',   'ACTA',   'ALD2',    'ALK4',  'ATO2',
     &     'B3O2',       'Br',   'BrCl',    'BrO',  'BrONO2',  'C2H6',
     &     'C3H8',     'CCl4', 'CF2Br2', 'CF2Cl2', 'CF2ClBr', 'CF3Br',
     &   'CFC113',   'CFC114', 'CFC115',  'CFCl3',    'CH2O', 'CH3Br',
     &  'CH3CCl3',    'CH3Cl',    'CH4',     'CO',      'Cl',   'Cl2',
     &    'Cl2O2',      'ClO', 'ClONO2',    'EOH',    'ETO2',   'ETP',
     &     'GCO3',     'GLYC',   'GLYX',     'GP',    'GPAN',     'H',
     &       'H2',    'H2402',    'H2O',   'H2O2',     'HAC',   'HBr',
     & 'HCFC141b', 'HCFC142b', 'HCFC22',  'HCOOH',     'HCl',  'HNO2',
     &     'HNO3',     'HNO4',    'HO2',   'HOBr',    'HOCl',  'IALD',
     &     'IAO2',      'IAP',   'INO2',   'INPN',    'ISN1',  'ISNP',
     &     'ISOP',      'KO2',   'MACR',   'MAN2',    'MAO3',  'MAOP',
     &      'MAP',     'MCO3',    'MEK',   'MGLY',     'MO2',   'MOH',
     &       'MP',     'MRO2',    'MRP',    'MVK',    'MVN2',     'N',
     &      'N2O',     'N2O5',     'NO',    'NO2',     'NO3',   'NOx',
     &        'O',      'O1D',     'O3',   'OClO',      'OH',    'Ox',
     &      'PAN',      'PMN',    'PO2',     'PP',     'PPN',  'PRN1',
     &     'PRPE',     'PRPN',   'R4N1',   'R4N2',    'R4O2',   'R4P',
     &     'RA3P',     'RB3P',   'RCHO',   'RCO3',   'RCOOH',  'RIO1',
     &     'RIO2',      'RIP',    'ROH',     'RP',    'VRO2',   'VRP' /)


      !===========================!
      ! Full chemistry simulation !
      !===========================!
      IF ( ITS_A_FULLCHEM_SIM() ) THEN

         DO NN = 1, NTR_GMI

            sname = TRIM(GMI_TrName(NN))

            DO N = 1, N_TRACERS

                IF ( TRIM(TRACER_NAME(N)) .eq. TRIM(sname) ) THEN

                   IF ( LLINOZ .and.
     &             TRIM(TRACER_NAME(N)) .eq. 'Ox' ) THEN
                      WRITE(6,*), TRIM(TRACER_NAME(N)) // ' (via Linoz)'

                   ! For adjoint
                   !ELSE IF ( TRIM(TRACER_NAME(N)) .eq. 'Ox' ) THEN
                   !   WRITE(6,*) TRIM(TRACER_NAME(N)) // ' (via Synoz)'
                   ELSEIF ( TRIM(TRACER_NAME(N)) .eq. 'Ox' ) THEN
                      WRITE(6,*),TRIM(TRACER_NAME(N)) // ' (Ox via GMI)'

                   ELSE
                      WRITE(6,*), TRIM(TRACER_NAME(N)) //
     &                           ' (via GMI rates)'
                   ENDIF

                   NSCHEM                 = NSCHEM + 1
                   Strat_TrID_GC(NSCHEM)  = N  ! Maps 1:NSCHEM to STT index
                   Strat_TrID_GMI(NSCHEM) = NN ! Maps 1:NSCHEM to GMI_TrName index

                ENDIF

             ENDDO
          ENDDO

          ! Allocate array to hold monthly mean OH mixing ratio
          ALLOCATE( STRAT_OH( IIPAR, JJPAR, LLPAR ), STAT=AS )
          IF ( AS /=0 ) CALL ALLOC_ERR( 'STRAT_OH' )
          STRAT_OH = 0d0

      !===========!
      ! Tagged Ox !
      !===========!
      ELSE IF ( ITS_A_TAGOX_SIM() ) THEN

         IF ( LLINOZ ) THEN
            WRITE(6,*) 'Linoz ozone performed on: '
         ELSE IF ( .not. LADJ ) THEN
            ! must use Linoz or strat chem Ox fluxes for the adjoint
            CALL UPBDFLX_O3     ! Synoz
            WRITE(6,*) 'Synoz ozone performed on: '
         ENDIF

         DO N = 1, N_TRACERS
            IF ( TRIM(TRACER_NAME(N)) .eq. 'Ox' .or.
     &           TRIM(TRACER_NAME(N)) .eq. 'OxStrt' ) THEN
               NSCHEM = NSCHEM + 1
               Strat_TrID_GC(NSCHEM) = N
               WRITE(6,*) TRIM(TRACER_NAME(N))
            ENDIF
         ENDDO

      ENDIF

      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
      ! Allocate and initialize prod & loss arrays         !
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

      ! Allocate PROD -- array for clim. production rates [v/v/s]
      ALLOCATE( PROD( IIPAR, JJPAR, LLPAR, NSCHEM ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'PROD' )
      PROD = 0d0

      ! Allocate LOSS -- array for clim. loss freq [s-1]
      ALLOCATE( LOSS( IIPAR, JJPAR, LLPAR, NSCHEM ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'LOSS' )
      LOSS = 0d0

      ! For adjoint
      !ALLOCATE( PROD_0( IIPAR, JJPAR, LLPAR, N_TRACERS ), STAT=AS )
      ALLOCATE( PROD_0( IIPAR, JJPAR, LLPAR, NSCHEM ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'PROD_0' )
      PROD_0 = 0d0

      !ALLOCATE( LOSS_0( IIPAR, JJPAR, LLPAR, N_TRACERS ), STAT=AS )
      ALLOCATE( LOSS_0( IIPAR, JJPAR, LLPAR, NSCHEM ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'LOSS_0' )
      LOSS_0 = 0d0


      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
      ! Allocate and initialize arrays for STE calculation !
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

      ! Array to hold initial state of atmosphere at the beginning
      ! of the period over which to estimate STE. Populate with
      ! initial atm. conditions from restart file [kg].
      ALLOCATE( MInit( IIPAR, JJPAR, LLPAR, N_TRACERS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'MInit' )
      MInit = STT

      ! Array to determine the mean tropopause level over the period
      ! for which STE is being estimated.
      ALLOCATE( TPAUSEL( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'TPAUSEL' )
      TPAUSEL = 0d0

      ! Array to save chemical state before each chemistry time step [kg]
      ALLOCATE( BEFORE( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BEFORE' )
      BEFORE = 0d0

      ! Array to aggregate the stratospheric chemical tendency [kg period-1]
      ALLOCATE( SCHEM_TEND(IIPAR,JJPAR,LLPAR,N_TRACERS), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'SCHEM_TEND' )
      SCHEM_TEND = 0d0


      END SUBROUTINE INIT_STRAT_CHEM
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_strat_chem
!
! !DESCRIPTION: Subroutine CLEANUP\_STRAT\_CHEM deallocates all module
!  arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_STRAT_CHEM
!
! !USES:
!      USE NETCDF_UTIL_MOD
! !REVISION HISTORY:
!  1 Feb 2011 - L. Murray - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC

      IF ( ALLOCATED( PROD       ) ) DEALLOCATE( PROD       )
      IF ( ALLOCATED( LOSS       ) ) DEALLOCATE( LOSS       )
      IF ( ALLOCATED( PROD_0     ) ) DEALLOCATE( PROD_0     )
      IF ( ALLOCATED( LOSS_0     ) ) DEALLOCATE( LOSS_0     )
      IF ( ALLOCATED( STRAT_OH   ) ) DEALLOCATE( STRAT_OH   )
      IF ( ALLOCATED( MInit      ) ) DEALLOCATE( MInit      )
      IF ( ALLOCATED( TPAUSEL    ) ) DEALLOCATE( TPAUSEL    )
      IF ( ALLOCATED( BEFORE     ) ) DEALLOCATE( BEFORE     )
      IF ( ALLOCATED( SCHEM_TEND ) ) DEALLOCATE( SCHEM_TEND )



      END SUBROUTINE CLEANUP_STRAT_CHEM
!EOC
      END MODULE STRAT_CHEM_MOD
