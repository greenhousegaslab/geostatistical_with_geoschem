!------------------------------------------------------------------------------
!          Matthew Cooper - Dalhousie University                              !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: regrid_a2a_mod.F90
!
! !DESCRIPTION: Module REGRID\_A2A\_MOD uses an algorithm adapted from MAP\_A2A
!   code to regrid from one horizonatal grid to another.
!\\
!\\
! !INTERFACE:
!
MODULE REGRID_A2A_MOD
!
! !USES:
!
  IMPLICIT NONE
  PRIVATE
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: XMAP
  PRIVATE :: YMAP
  PRIVATE :: READ_INPUT_GRID
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: DO_REGRID_A2A
  PUBLIC  :: DO_REGRID_DKH
  PUBLIC  :: MAP_A2A
!
! !REVISION HISTORY:
!  13 Mar 2012 - M. Cooper   - Initial version
!  03 Apr 2012 - M. Payer    - Now use functions GET_AREA_CM2(I,J,L),
!                              GET_YEDGE(I,J,L) and GET_YSIN(I,J,L) from the
!                              new grid_mod.F90
!  22 May 2012 - L. Murray   - Implemented several bug fixes
!  23 Aug 2012 - R. Yantosca - Add capability for starting from hi-res grids
!                              (generic 0.5x0.5, generic 0.25x0.25, etc.)
!  23 Aug 2012 - R. Yantosca - Add subroutine READ_INPUT_GRID, which reads the
!                              grid parameters (lon & lat edges) w/ netCDF
!  27 Aug 2012 - R. Yantosca - Now parallelize key DO loops
!EOP
!------------------------------------------------------------------------------
!BOC
  CONTAINS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_regrid_a2a
!
! !DESCRIPTION: Subroutine DO\_REGRID\_A2A regrids 2-D data in the
!  horizontal direction.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DO_REGRID_A2A( FILENAME, IM, JM, INGRID, OUTGRID, IS_MASS, &
                            netCDF )
!
! !USES:
!
    USE GRID_MOD,   ONLY : GET_XEDGE
    USE GRID_MOD,   ONLY : GET_YSIN
    USE GRID_MOD,   ONLY : GET_AREA_CM2
    USE GRID_MOD,   ONLY : GET_IJ
    USE GRID_MOD,   ONLY : GET_IJ_GLOBAL
    USE FILE_MOD,   ONLY : IOERROR
    USE inquireMod, ONLY : findFreeLUN

#     include "CMN_SIZE"  ! Size parameters
#     include "CMN_GCTM"  ! Size parameters
!
! !INPUT PARAMETERS:
!
    ! Name of file with lon and lat edge information on the INPUT GRID
    CHARACTER(LEN=*), INTENT(IN)    :: FILENAME

    ! Number of lon centers and lat centers on the INPUT GRID
    INTEGER,          INTENT(IN)    :: IM
    INTEGER,          INTENT(IN)    :: JM

    ! Data array on the input grid
    REAL*8,           INTENT(IN)    :: INGRID(IM,JM)

    ! IS_MASS=0 if data is units of concentration (molec/cm2/s, unitless, etc.)
    ! IS_MASS=1 if data is units of mass (kg/yr, etc.); we will need to convert
    !           INGRID to per unit area
    INTEGER,          INTENT(IN)    :: IS_MASS

    ! Read from netCDF file?  (needed for debugging, will disappear later)
    LOGICAL, OPTIONAL,INTENT(IN)    :: netCDF
!
! !OUTPUT PARAMETERS:
!
    ! Data array on the OUTPUT GRID
    REAL*8,           INTENT(OUT)   :: OUTGRID(IIPAR,JJPAR)
!
! !REVISION HISTORY:

!  13 Mar 2012 - M. Cooper   - Initial version
!  22 May 2012 - L. Murray   - Bug fix: INSIN should be allocated w/ JM+1.
!  22 May 2012 - R. Yantosca - Updated comments, cosmetic changes
!  25 May 2012 - R. Yantosca - Bug fix: declare the INGRID argument as
!                              INTENT(IN) to preserve the values of INGRID
!                              in the calling routine
!  06 Aug 2012 - R. Yantosca - Now make IU_REGRID a local variable
!  06 Aug 2012 - R. Yantosca - Move calls to findFreeLUN out of DEVEL block
!  23 Aug 2012 - R. Yantosca - Now use f10.4 format for hi-res grids
!  23 Aug 2012 - R. Yantosca - Now can read grid info from netCDF files
!  27 Aug 2012 - R. Yantosca - Add parallel DO loops
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER           :: I,        J
    INTEGER           :: IOS,      M
    INTEGER           :: IU_REGRID
    REAL*8            :: INAREA,   RLAT
    CHARACTER(LEN=15) :: HEADER1
    CHARACTER(LEN=20) :: FMT_LAT,  FMT_LON, FMT_LEN
    LOGICAL           :: USE_NETCDF

    ! Arrays
    REAL*8            :: INLON  (IM   +1)  ! Lon edges        on INPUT GRID
    REAL*8            :: INSIN  (JM   +1)  ! SIN( lat edges ) on INPUT GRID
    REAL*8            :: OUTLON (IIPAR+1)  ! Lon edges        on OUTPUT GRID
    REAL*8            :: OUTSIN (JJPAR+1)  ! SIN( lat edges ) on OUTPUT GRID
    REAL*8            :: IN_GRID(IM,JM  )  ! Shadow variable for INGRID

    ! dkh debug
    REAL*8            :: total
    REAL*4            :: LON
    REAL*4            :: LAT
    INTEGER           :: IIJJ(2)

    !======================================================================
    ! Initialization
    !
    ! NOTE: In the near future ASCII input will be replaced by netCDF!
    !======================================================================

    ! Save value of netCDF to shadow variable
    IF ( PRESENT( netCDF ) ) THEN
       USE_netCDF = netCDF
    ELSE
       USE_netCDF = .FALSE.
    ENDIF

    ! Longitude edges on the OUTPUT GRID
    ! NOTE: May have to make OUTLON a 2-D array later for the GI model
    DO I = 1, IIPAR+1
       OUTLON(I) = GET_XEDGE( I )
    ENDDO

    ! SIN( lat edges ) on the OUTPUT GRID
    ! NOTE: May have to make OUTSIN a 2-D array later for the GI model
    DO J = 1, JJPAR+1
       OUTSIN(J) = GET_YSIN( 1, J, 1 )
    ENDDO

    ! Read the input grid specifications
    IF ( USE_netCDF ) THEN

       !------------------------------------------
       ! %%% FROM NETCDF FILE %%%
       !------------------------------------------

       ! Read the grid specifications from a netCDF file
       CALL READ_INPUT_GRID( IM, JM, FILENAME, INLON, INSIN )

    ELSE

       !------------------------------------------
       ! %%% FROM ASCII FILE %%%
       !
       ! NOTE: Deprecated, will be removed later.
       !------------------------------------------

       ! Find a free file LUN
       IU_REGRID = findFreeLUN()

       ! Open file containing lon & lat edges on the INPUT GRID
       OPEN( IU_REGRID, FILE=TRIM( FILENAME ), STATUS='OLD', IOSTAT=IOS )
       IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_REGRID, 'latlonread' )

       ! Create the approprate FORMAT strings
       WRITE(FMT_LEN,*) IM+1

       ! NOTE: If the resolution of the grid is high enough, we have
       ! to allow for an extra digit in the input file.  This will
       ! become obsolete once we migrate to netCDF format (bmy, 8/23/12)
       IF ( IM > 1000 ) THEN
          FMT_LON='(' // TRIM ( FMT_LEN ) // 'F10.4)'   ! For hi-res grids
       ELSE
          FMT_LON='(' // TRIM ( FMT_LEN ) // 'F9.3)'    ! For all other grids
       ENDIF

       WRITE(FMT_LEN,*) JM
       FMT_LAT='(' // TRIM ( FMT_LEN ) // 'F15.10)'

       ! Read lon edges & SIN( lat edges ) on the INPUT GRID
       READ( IU_REGRID, '(A15)',IOSTAT=IOS ) HEADER1
       READ( IU_REGRID,FMT_LON,IOSTAT=IOS  ) ( INLON(M), M=1,IM+1 )
       READ( IU_REGRID,FMT_LAT,IOSTAT=IOS  ) ( INSIN(M), M=1,JM+1 )

       ! Close file
       CLOSE( IU_REGRID )

    ENDIF



    !======================================================================
    ! Regridding
    !======================================================================

    ! Copy the input argument INGRID to a local shadow variable,
    ! so that we can preserve the value of INGRID in the calling routine
    IN_GRID = INGRID

    ! Convert input to per area units if necessary
    IF ( IS_MASS == 1 ) THEN

       !$OMP PARALLEL DO                   &
       !$OMP DEFAULT( SHARED             ) &
       !$OMP PRIVATE( I, J, RLAT, INAREA )
       DO J = 1, JM
          RLAT   = INSIN(J+1) - INSIN(J)
          INAREA = ( 2d0 * PI * Re * RLAT * 1d4 * Re ) / DBLE( IM )
          DO I = 1, IM
             IN_GRID(I,J) = IN_GRID(I,J) / INAREA
          ENDDO
       ENDDO
       !$OMP END PARALLEL DO

    ENDIF

    ! dkh debug
    ! display total ingrid
    IF ( IS_MASS == 3 ) THEN
       total = 0d0
       DO J = 1, JM
          RLAT   = INSIN(J+1) - INSIN(J)
          INAREA = ( 2d0 * PI * Re * RLAT * 1d4 * Re ) / DBLE( IM )
          DO I = 1, IM
             total = total + IN_GRID(I,J) * INAREA * 1d-4 ! 1d-4 because INAREA is cm2, but HTAP is per m2
          ENDDO
       ENDDO
       ! and convert from kg per s to kg per month (July)
       print*, ' HTAP CO sum in = ', total * 60d0 * 60d0 * 24d0 * 31d0
    ENDIF

!       DO J = 1, JM
!          RLAT   = INSIN(J+1) - INSIN(J)
!          INAREA = ( 2d0 * PI * Re * RLAT * 1d4 * Re ) / DBLE( IM )
!          DO I = 1, IM
!             LAT = -179.95d0 + (JM - 1 ) * 0.1d0
!             LON = - 89.95d0 + (IM - 1 ) * 0.1d0
!             IIJJ = GET_IJ(LON, LAT)
!             OUTGRID(IIJJ(1),IIJJ(2)) = OUTGRID(IIJJ(1),IIJJ(2))  &
!                                      + IN_GRID(I,J) * INAREA * 1d-4 ! 1d-4 because INAREA is cm2, but HTAP is per m2
!          ENDDO
!       ENDDO
!       DO J = 1, JJPAR
!       DO I = 1, IIPAR
!          OUTGRID(I,J) = OUTGRID(I,J) / GET_AREA_CM2( J ) * 1d4 ! cm2 to m2
!       ENDDO
!       ENDDO

    ! Call MAP_A2A to do the regridding
    CALL MAP_A2A( IM,    JM,    INLON,  INSIN,  IN_GRID,        &
                  IIPAR, JJPAR, OUTLON, OUTSIN, OUTGRID, 0, 0 )

    ! Convert back from "per area" if necessary
    IF ( IS_MASS == 1 ) THEN

       !$OMP PARALLEL DO       &
       !$OMP DEFAULT( SHARED ) &
       !$OMP PRIVATE( I, J   )
       DO J = 1, JJPAR
       DO I = 1, IIPAR
          OUTGRID(I,J) = OUTGRID(I,J) * GET_AREA_CM2( J )
       ENDDO
       ENDDO
       !$OMP END PARALLEL DO

    ENDIF

    ! dkh debug
    ! display total ingrid
    IF ( IS_MASS == 3 ) THEN
       total = 0d0
       DO J = 1, JJPAR
       DO I = 1, IIPAR
          total = total + OUTGRID(I,J) * GET_AREA_CM2( J ) * 1d-4 ! cm2 to m2
       ENDDO
       ENDDO
       print*, ' HTAP CO sum in = ', total * 60d0 * 60d0 * 24d0 * 31d0
    ENDIF
  END SUBROUTINE DO_REGRID_A2A
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_regrid_dkh
!
! !DESCRIPTION: Subroutine DO\_REGRID\_DKH regrids 2-D data in the
!  horizontal direction.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DO_REGRID_DKH( FILENAME, IM, JM, INGRID, OUTGRID, IS_MASS, &
                            netCDF )
!
! !USES:
!
    USE GRID_MOD,   ONLY : GET_XEDGE
    USE GRID_MOD,   ONLY : GET_YSIN
    USE GRID_MOD,   ONLY : GET_AREA_CM2
    USE GRID_MOD,   ONLY : GET_IJ
    USE GRID_MOD,   ONLY : GET_IJ_GLOBAL
    USE FILE_MOD,   ONLY : IOERROR
    USE inquireMod, ONLY : findFreeLUN
    USE GRID_MOD,   ONLY : GET_XOFFSET, GET_YOFFSET

#     include "CMN_SIZE"  ! Size parameters
#     include "CMN_GCTM"  ! Size parameters
!
! !INPUT PARAMETERS:
!
    ! Name of file with lon and lat edge information on the INPUT GRID
    CHARACTER(LEN=*), INTENT(IN)    :: FILENAME

    ! Number of lon centers and lat centers on the INPUT GRID
    INTEGER,          INTENT(IN)    :: IM
    INTEGER,          INTENT(IN)    :: JM

    ! Data array on the input grid
    REAL*8,           INTENT(IN)    :: INGRID(IM,JM)

    ! =1 if we need to convert INGRID to per unit area
    INTEGER,          INTENT(IN)    :: IS_MASS

    ! Read from netCDF file?  (needed for debugging, will disappear later)
    LOGICAL, OPTIONAL,INTENT(IN)    :: netCDF
!
! !OUTPUT PARAMETERS:
!
    ! Data array on the OUTPUT GRID
    REAL*8,           INTENT(OUT)   :: OUTGRID(IIPAR,JJPAR)
!
! !REVISION HISTORY:

!  17 Nov 2013 - D. Henze    - Initial version, based on A2A
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER           :: I,        J
    INTEGER           :: IOS,      M
    INTEGER           :: IU_REGRID
    REAL*8            :: INAREA,   RLAT
    CHARACTER(LEN=15) :: HEADER1
    CHARACTER(LEN=20) :: FMT_LAT,  FMT_LON, FMT_LEN
    LOGICAL           :: USE_NETCDF

    ! Arrays
    REAL*8            :: INLON  (IM   +1)  ! Lon edges        on INPUT GRID
    REAL*8            :: INSIN  (JM   +1)  ! SIN( lat edges ) on INPUT GRID
    REAL*8            :: OUTLON (IIPAR+1)  ! Lon edges        on OUTPUT GRID
    REAL*8            :: OUTSIN (JJPAR+1)  ! SIN( lat edges ) on OUTPUT GRID
    REAL*8            :: IN_GRID(IM,JM  )  ! Shadow variable for INGRID

    ! (yd 2015/6/19 Global offset for nested)
    INTEGER           :: P,        Q
    REAL*8            :: P0,       Q0
    REAL*8            :: OUTGRID_GLOB(IIPAR_L,JJPAR_L)
    !(quzhen 2016/2/1)
    INTEGER           :: LON_COUNT
    INTEGER           :: IIJJ_PREV !lon index of previous longitude
    REAL*8            :: RATIO

    ! dkh debug
    REAL*8            :: total
    REAL*4            :: LON
    REAL*4            :: LAT
    INTEGER           :: IIJJ(2)

    !======================================================================
    ! Initialization
    !
    ! NOTE: In the near future ASCII input will be replaced by netCDF!
    !======================================================================

    ! Save value of netCDF to shadow variable
    IF ( PRESENT( netCDF ) ) THEN
       USE_netCDF = netCDF
    ELSE
       USE_netCDF = .FALSE.
    ENDIF

    ! Longitude edges on the OUTPUT GRID
    ! NOTE: May have to make OUTLON a 2-D array later for the GI model
    DO I = 1, IIPAR+1
       OUTLON(I) = GET_XEDGE( I )
    ENDDO

    ! SIN( lat edges ) on the OUTPUT GRID
    ! NOTE: May have to make OUTSIN a 2-D array later for the GI model
    DO J = 1, JJPAR+1
       OUTSIN(J) = GET_YSIN( 1, J, 1 )
    ENDDO

    ! Read the input grid specifications
    IF ( USE_netCDF ) THEN

       !------------------------------------------
       ! %%% FROM NETCDF FILE %%%
       !------------------------------------------

       ! Read the grid specifications from a netCDF file
       CALL READ_INPUT_GRID( IM, JM, FILENAME, INLON, INSIN )

    ELSE

       !------------------------------------------
       ! %%% FROM ASCII FILE %%%
       !
       ! NOTE: Deprecated, will be removed later.
       !------------------------------------------

       ! Find a free file LUN
       IU_REGRID = findFreeLUN()

       ! Open file containing lon & lat edges on the INPUT GRID
       OPEN( IU_REGRID, FILE=TRIM( FILENAME ), STATUS='OLD', IOSTAT=IOS )
       IF ( IOS /= 0 ) CALL IOERROR( IOS, IU_REGRID, 'latlonread' )

       ! Create the approprate FORMAT strings
       WRITE(FMT_LEN,*) IM+1

       ! NOTE: If the resolution of the grid is high enough, we have
       ! to allow for an extra digit in the input file.  This will
       ! become obsolete once we migrate to netCDF format (bmy, 8/23/12)
       IF ( IM > 1000 ) THEN
          FMT_LON='(' // TRIM ( FMT_LEN ) // 'F10.4)'   ! For hi-res grids
       ELSE
          FMT_LON='(' // TRIM ( FMT_LEN ) // 'F9.3)'    ! For all other grids
       ENDIF

       WRITE(FMT_LEN,*) JM
       FMT_LAT='(' // TRIM ( FMT_LEN ) // 'F15.10)'

       ! Read lon edges & SIN( lat edges ) on the INPUT GRID
       READ( IU_REGRID, '(A15)',IOSTAT=IOS ) HEADER1
       READ( IU_REGRID,FMT_LON,IOSTAT=IOS  ) ( INLON(M), M=1,IM+1 )
       READ( IU_REGRID,FMT_LAT,IOSTAT=IOS  ) ( INSIN(M), M=1,JM+1 )

       ! Close file
       CLOSE( IU_REGRID )

    ENDIF



    !======================================================================
    ! Regridding
    !======================================================================

    ! Copy the input argument INGRID to a local shadow variable,
    ! so that we can preserve the value of INGRID in the calling routine
    IN_GRID = INGRID

    ! Convert input to per area units if necessary
    IF ( IS_MASS == 1 ) THEN

       !$OMP PARALLEL DO                   &
       !$OMP DEFAULT( SHARED             ) &
       !$OMP PRIVATE( I, J, RLAT, INAREA )
       DO J = 1, JM
          RLAT   = INSIN(J+1) - INSIN(J)
          INAREA = ( 2d0 * PI * Re * RLAT * 1d4 * Re ) / DBLE( IM )
          DO I = 1, IM
             IN_GRID(I,J) = IN_GRID(I,J) / INAREA
          ENDDO
       ENDDO
       !$OMP END PARALLEL DO

    ENDIF

       OUTGRID_GLOB = 0d0

       DO J = 1, JM
          RLAT   = INSIN(J+1) - INSIN(J)
          INAREA = ( 2d0 * PI * Re * RLAT * 1d4 * Re ) / DBLE( IM )
          DO I = 1, IM
             LON =    0.05d0 + (I - 1 ) * 0.1d0
             IF ( LON > 180d0 ) LON = LON - 360d0
             LAT = - 89.95d0 + (J - 1 ) * 0.1d0
             !(quzhen 2014/12/31)
             IIJJ = GET_IJ_GLOBAL(LON, LAT)

#if defined(NESTED_CH)

             !(quzhen 2016/2/1 ONLY DO THIS IN NESTED_CH DOMAIN)
             RATIO = 1d0
             IF (LON .GT. 69.65 .AND. LON .LT. 150.35)  THEN  !in domain
               IF (MOD(IIJJ(1), 3) .EQ. 1) THEN
                 ! lower boundary
                 IF (IIJJ(1) .NE. IIJJ_PREV) THEN ! a new grid
                   RATIO = (0.75 - 2./3.)/0.1
                   LON_COUNT = 0
                 ENDIF
                 LON_COUNT = LON_COUNT + 1
                 IF(LON_COUNT .EQ. 6) THEN
                   ! upper boundary
                   RATIO = (1./3. - 0.25)/0.1
                 ENDIF
               ENDIF

               IF (MOD(IIJJ(1), 3) .EQ. 2) THEN ! 2nd pattern
                 ! lower boundary
                 IF (IIJJ(1) .NE. IIJJ_PREV) THEN
                   RATIO = (0.35 - 1./3.)/0.1
                   LON_COUNT = 0
                 ENDIF
                 LON_COUNT = LON_COUNT + 1
                 IF(LON_COUNT .EQ. 7) THEN
                   ! upper boundary
                   RATIO = 0.05/0.1
                 ENDIF
               ENDIF

               IF (MOD(IIJJ(1), 3) .EQ. 0) THEN ! 3rd pattern
                 ! lower boundary
                 IF (IIJJ(1) .NE. IIJJ_PREV) THEN
                   !INAREA = INAREA * 0.05/0.1
                   RATIO = 0.05/0.1
                   LON_COUNT = 0
                 ENDIF
                 LON_COUNT = LON_COUNT + 1
                 !print*,'quzhen 3rd pattern,lon_count',lon_count, ratio
                 IF(LON_COUNT .EQ. 7) THEN
                   ! upper boundary
                   !INAREA = INAREA * (2./3. - 0.65)/0.1
                   RATIO = (2./3. - 0.65)/0.1
                 ENDIF
               ENDIF

             ENDIF
             IIJJ_PREV = IIJJ(1)

#endif
             !(quzhen 2015/1/23 regrid globally first)
             OUTGRID_GLOB(IIJJ(1),IIJJ(2)) = OUTGRID_GLOB(IIJJ(1),IIJJ(2))  &
                                      + IN_GRID(I,J) * INAREA * 1d-4 * RATIO ! 1d-4 because INAREA is cm2, but HTAP is per m2
          ENDDO
       ENDDO

#if (defined (NESTED_NA) || defined (NESTED_CH) || defined (NESTED_SD) )
       !(quzhen 2015/1/23 select data in the window area)
       P0 = GET_XOFFSET( GLOBAL=.TRUE. )
       Q0 = GET_YOFFSET( GLOBAL=.TRUE. )

       DO P = 1, IIPAR
       DO Q = 1, JJPAR
          OUTGRID( P, Q ) = OUTGRID_GLOB ( P+P0, Q+Q0)
       ENDDO
       ENDDO
#else
       OUTGRID = OUTGRID_GLOB
#endif

       DO J = 1, JJPAR
       DO I = 1, IIPAR
          OUTGRID(I,J) = OUTGRID(I,J) / GET_AREA_CM2( J ) * 1d4 ! cm2 to m2
       ENDDO
       ENDDO

    ! Convert back from "per area" if necessary
    IF ( IS_MASS == 1 ) THEN

       !$OMP PARALLEL DO       &
       !$OMP DEFAULT( SHARED ) &
       !$OMP PRIVATE( I, J   )
       DO J = 1, JJPAR
       DO I = 1, IIPAR
          OUTGRID(I,J) = OUTGRID(I,J) * GET_AREA_CM2( J )
       ENDDO
       ENDDO
       !$OMP END PARALLEL DO

    ENDIF

  END SUBROUTINE DO_REGRID_DKH
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: map_a2a
!
! !DESCRIPTION: Subroutine MAP\_A2A is a horizontal arbitrary grid to arbitrary
!  grid conservative high-order mapping regridding routine by S-J Lin.
!\\
!\\
! !INTERFACE:
!
!  (1 ) INLON   (REAL*8   ) : Longitude edges of input grid
!  (2 ) INSIN   (REAL*8   ) : Sine of input grid latitude edges
!  (3 ) INGRID  (REAL*8   ) : Data array to be regridded

  SUBROUTINE map_a2a( im, jm, lon1, sin1, q1, &
                      in, jn, lon2, sin2, q2, ig, iv)
!
! !INPUT PARAMETERS:
!
    ! Longitude and Latitude dimensions of INPUT grid
    INTEGER, INTENT(IN)  :: im, jm

    ! Longitude and Latitude dimensions of OUTPUT grid
    INTEGER, INTENT(IN)  :: in, jn

    ! IG=0: pole to pole;
    ! IG=1 J=1 is half-dy north of south pole
    INTEGER, INTENT(IN)  :: ig

    ! IV=0: Regrid scalar quantity
    ! IV=1: Regrid vector quantity
    INTEGER, INTENT(IN)  :: iv

    ! Longitude edges (degrees) of INPUT and OUTPUT grids
    REAL*8,  INTENT(IN)  :: lon1(im+1), lon2(in+1)

    ! Sine of Latitude Edges (radians) of INPUT and OUTPUT grids
    REAL*8,  INTENT(IN)  :: sin1(jm+1), sin2(jn+1)

    ! Quantity on INPUT grid
    REAL*8,  INTENT(IN)  :: q1(im,jm)
!
! !OUTPUT PARAMETERS:
!
    ! Regridded quantity on OUTPUT grid
    REAL*8,  INTENT(OUT) :: q2(in,jn)
!
!  !REVISION HISTORY:
!  (1) Original subroutine by S-J Lin.  Converted to F90 freeform format
!      and inserted into "Geos3RegridModule" by Bob Yantosca (9/21/00)
!  (2) Added F90 type declarations to be consistent w/ TypeModule.f90.
!      Also updated comments. (bmy, 9/21/00)
!  21 Sep 2000 - R. Yantosca - Initial version
!  27 Aug 2012 - R. Yantosca - Add parallel DO loops
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: i,j,k
    REAL*8  :: qtmp(in,jm)

    !===================================================================
    ! E-W regridding
    !===================================================================
    IF ( im .eq. in ) THEN

       ! Don't call XMAP if both grids have the same # of longitudes
       ! but save the input data in the QTMP array
       !$OMP PARALLEL DO       &
       !$OMP DEFAULT( SHARED ) &
       !$OMP PRIVATE( I, J )
       DO j=1,jm-ig
       DO i=1,im
          qtmp(i,j+ig) = q1(i,j+ig)
       ENDDO
       ENDDO
       !$OMP END PARALLEL DO

    ELSE

       ! Otherwise, call XMAP to regrid in the E-W direction
       CALL xmap(im, jm-ig, lon1, q1(1,1+ig),in, lon2, qtmp(1,1+ig) )

    ENDIF

    !===================================================================
    ! N-S regridding
    !===================================================================
    IF ( jm .eq. jn ) THEN

       ! Don't call XMAP if both grids have the same # of longitudes,
       ! but assign the value of QTMP to the output Q2 array
       !$OMP PARALLEL DO       &
       !$OMP DEFAULT( SHARED ) &
       !$OMP PRIVATE( I, J )
       DO j=1,jm-ig
       DO i=1,in
          q2(i,j+ig) = qtmp(i,j+ig)
       ENDDO
       ENDDO
       !$OMP END PARALLEL DO

    ELSE

       ! Otherwise, call YMAP to regrid in the N-S direction
       CALL ymap(in, jm, sin1, qtmp(1,1+ig), jn, sin2, q2(1,1+ig), ig, iv)

    ENDIF

  END SUBROUTINE map_a2a
!EOC
!------------------------------------------------------------------------------
!                   Prasad Kasibhatla - Duke University                       !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ymap
!
! !DESCRIPTION: Routine to perform area preserving mapping in N-S from an
!  arbitrary resolution to another.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ymap(im, jm, sin1, q1, jn, sin2, q2, ig, iv)
!
! !INPUT PARAMETERS:
!

    ! original E-W dimension
    INTEGER, INTENT(IN)  :: im

    ! original N-S dimension
    INTEGER, INTENT(IN)  :: jm

    ! Target N-S dimension
    INTEGER, INTENT(IN)  :: jn

    ! IG=0: scalars from SP to NP (D-grid v-wind is also IG=0)
    ! IG=1: D-grid u-wind
    INTEGER, INTENT(IN)  :: ig

    ! IV=0: scalar;
    ! IV=1: vector
    INTEGER, INTENT(IN)  :: iv

    ! Original southern edge of the cell sin(lat1)
    REAL*8,  INTENT(IN)  :: sin1(jm+1-ig)

    ! Original data at center of the cell
    REAL*8,  INTENT(IN)  :: q1(im,jm)

    ! Target cell's southern edge sin(lat2)
    REAL*8,  INTENT(IN)  :: sin2(jn+1-ig)
!
! !OUTPUT PARAMETERS:
!
    ! Mapped data at the target resolution
    REAL*8,  INTENT(OUT) :: q2(im,jn)
!
! !REMARKS:
!
!   sin1 (1) = -1 must be south pole; sin1(jm+1)=1 must be N pole.
!
!   sin1(1) < sin1(2) < sin1(3) < ... < sin1(jm) < sin1(jm+1)
!   sin2(1) < sin2(2) < sin2(3) < ... < sin2(jn) < sin2(jn+1)!
!
! !AUTHOR:
!   Developer: Prasad Kasibhatla
!   March 6, 2012
!
! !REVISION HISTORY
!  06 Mar 2012 - P. Kasibhatla - Initial version
!  27 Aug 2012 - R. Yantosca   - Added parallel DO loops
!  27 Aug 2012 - R. Yantosca   - Change REAL*4 variables to REAL*8 to better
!                                ensure numerical stability
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER              :: i, j0, m, mm, j
    REAL*8               :: dy1(jm)
    REAL*8               :: dy
!------------------------------------------------------------------------------
! Prior to 8/27/12:
! Change REAL*4 to REAL*8, to eliminate numerical noise (bmy, 8/27/12)
!    REAL*4               :: qsum, sum
!------------------------------------------------------------------------------
    REAL*8               :: qsum, sum

    ! YMAP begins here!
    do j=1,jm-ig
       dy1(j) = sin1(j+1) - sin1(j)
    enddo

    !===============================================================
    ! Area preserving mapping
    !===============================================================

    !$OMP PARALLEL DO                          &
    !$OMP DEFAULT( SHARED                    ) &
    !$OMP PRIVATE( I, J0, J, M, QSUM, MM, DY )
    do 1000 i=1,im
       j0 = 1
       do 555 j=1,jn-ig
       do 100 m=j0,jm-ig

          !=========================================================
          ! locate the southern edge: sin2(i)
          !=========================================================
          if(sin2(j) .ge. sin1(m) .and. sin2(j) .le. sin1(m+1)) then

             if(sin2(j+1) .le. sin1(m+1)) then

                ! entire new cell is within the original cell
                q2(i,j)=q1(i,m)
                j0 = m
                goto 555
             else

                ! South most fractional area
                qsum=(sin1(m+1)-sin2(j))*q1(i,m)

                do mm=m+1,jm-ig

                   ! locate the northern edge: sin2(j+1)
                   if(sin2(j+1) .gt. sin1(mm+1) ) then

                      ! Whole layer
                      qsum = qsum + dy1(mm)*q1(i,mm)
                   else

                      ! North most fractional area
                      dy = sin2(j+1)-sin1(mm)
                      qsum=qsum+dy*q1(i,mm)
                      j0 = mm
                      goto 123
                   endif
                enddo
                goto 123
             endif
          endif
100    continue
123    q2(i,j) = qsum / ( sin2(j+1) - sin2(j) )
555    continue
1000 continue
     !$OMP END PARALLEL DO

     !===================================================================
     ! Final processing for poles
     !===================================================================
     if ( ig .eq. 0 .and. iv .eq. 0 ) then

!------------------------------------------------------------------------------
! Prior to 8/27/12:
! Change REAL*4 to REAL*8, to eliminate numerical noise (bmy, 8/27/12)
!        ! South pole
!        sum = 0.
!        do i=1,im
!           sum = sum + q2(i,1)
!        enddo
!
!        sum = sum / float(im)
!        do i=1,im
!           q2(i,1) = sum
!        enddo
!
!        ! North pole:
!        sum = 0.
!        do i=1,im
!           sum = sum + q2(i,jn)
!        enddo
!
!        sum = sum / float(im)
!        do i=1,im
!           q2(i,jn) = sum
!        enddo
!------------------------------------------------------------------------------
        ! South pole
        sum = 0.d0
        do i=1,im
           sum = sum + q2(i,1)
        enddo

        sum = sum / DBLE( im )
        do i=1,im
           q2(i,1) = sum
        enddo

        ! North pole:
        sum = 0.d0
        do i=1,im
           sum = sum + q2(i,jn)
        enddo

        sum = sum / DBLE( im )
        do i=1,im
           q2(i,jn) = sum
        enddo

     endif

   END SUBROUTINE YMAP
!EOC
!------------------------------------------------------------------------------
!                   Prasad Kasibhatla - Duke University                       !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: xmap
!
! !DESCRIPTION: Routine to perform area preserving mapping in E-W from an
!  arbitrary resolution to another.
!  Periodic domain will be assumed, i.e., the eastern wall bounding cell
!  im is lon1(im+1) = lon1(1); Note the equal sign is true geographysically.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE xmap(im, jm, lon1, q1, in, lon2, q2)
!
! !INPUT PARAMETERS:
!
    ! Original E-W dimension
    INTEGER, INTENT(IN)  :: im

    ! Target E-W dimension
    INTEGER, INTENT(IN)  :: in

    ! Original N-S dimension
    INTEGER, INTENT(IN)  :: jm

    ! Original western edge of the cell
    REAL*8,  INTENT(IN)  :: lon1(im+1)

    ! Original data at center of the cell
    REAL*8,  INTENT(IN)  :: q1(im,jm)

    ! Target cell's western edge
    REAL*8,  INTENT(IN)  :: lon2(in+1)
!
! !OUTPUT PARAMETERS:
!
    ! Mapped data at the target resolution
    REAL*8,  INTENT(OUT) :: q2(in,jm)
!
! !REMARKS:
!   lon1(1) < lon1(2) < lon1(3) < ... < lon1(im) < lon1(im+1)
!   lon2(1) < lon2(2) < lon2(3) < ... < lon2(in) < lon2(in+1)
!
! !AUTHOR:
!   Developer: Prasad Kasibhatla
!   March 6, 2012
!
! !REVISION HISTORY
!  06 Mar 2012 - P. Kasibhatla - Initial version
!  27 Aug 2012 - R. Yantosca   - Added parallel DO loops
!  27 Aug 2012 - R. Yantosca   - Change REAL*4 variables to REAL*8 to better
!                                ensure numerical stability
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER              :: i1, i2, i, i0, m, mm, j
    REAL*8               :: qtmp(-im:im+im)
    REAL*8               :: x1(-im:im+im+1)
    REAL*8               :: dx1(-im:im+im)
    REAL*8               :: dx
!------------------------------------------------------------------------------
! Prior to 8/27/12:
! Change REAL*4 to REAL*8, to eliminate numerical noise (bmy, 8/27/12)
!    REAL*4               :: qsum
!------------------------------------------------------------------------------
    REAL*8               :: qsum
    LOGICAL              :: found

    ! XMAP begins here!
    do i=1,im+1
       x1(i) = lon1(i)
    enddo

    do i=1,im
       dx1(i) = x1(i+1) - x1(i)
    enddo

    !===================================================================
    ! check to see if ghosting is necessary
    ! Western edge:
    !===================================================================
    found = .false.
    i1 = 1
    do while ( .not. found )
       if( lon2(1) .ge. x1(i1) ) then
          found = .true.
       else
          i1 = i1 - 1
          if (i1 .lt. -im) then
             write(6,*) 'failed in xmap western edge '
             stop
          else
             x1(i1) = x1(i1+1) - dx1(im+i1)
             dx1(i1) = dx1(im+i1)
          endif
       endif
    enddo

    !===================================================================
    ! Eastern edge:
    !===================================================================
    found = .false.
    i2 = im+1
    do while ( .not. found )
       if( lon2(in+1) .le. x1(i2) ) then
          found = .true.
       else
          i2 = i2 + 1
          if (i2 .gt. 2*im) then
             write(6,*) 'failed in xmap eastern edge'
             stop
          else
             dx1(i2-1) = dx1(i2-1-im)
             x1(i2) = x1(i2-1) + dx1(i2-1)
          endif
       endif
    enddo

    !$OMP PARALLEL DO                                &
    !$OMP DEFAULT( SHARED                          ) &
    !$OMP PRIVATE( J, QTMP, I, I0, M, QSUM, MM, DX )
    do 1000 j=1,jm

       !=================================================================
       ! Area preserving mapping
       !================================================================

       qtmp(0)=q1(im,j)
       do i=1,im
          qtmp(i)=q1(i,j)
       enddo
       qtmp(im+1)=q1(1,j)

       ! check to see if ghosting is necessary
       ! Western edge
       if ( i1 .le. 0 ) then
          do i=i1,0
             qtmp(i) = qtmp(im+i)
          enddo
       endif

       ! Eastern edge:
       if ( i2 .gt. im+1 ) then
          do i=im+1,i2-1
             qtmp(i) = qtmp(i-im)
          enddo
       endif

       i0 = i1

       do 555 i=1,in
       do 100 m=i0,i2-1

          !=============================================================
          ! locate the western edge: lon2(i)
          !=============================================================
          if(lon2(i) .ge. x1(m) .and. lon2(i) .le. x1(m+1)) then

             if(lon2(i+1) .le. x1(m+1)) then

                ! entire new grid is within the original grid
                q2(i,j)=qtmp(m)
                i0 = m
                goto 555
             else

                ! Left most fractional area
                qsum=(x1(m+1)-lon2(i))*qtmp(m)
                do mm=m+1,i2-1

                   ! locate the eastern edge: lon2(i+1)
                   if(lon2(i+1) .gt. x1(mm+1) ) then

                      ! Whole layer
                      qsum = qsum + dx1(mm)*qtmp(mm)

                   else
                      ! Right most fractional area
                      dx = lon2(i+1)-x1(mm)
                      qsum=qsum+dx*qtmp(mm)
                      i0 = mm
                      goto 123
                   endif
                enddo
                goto 123
             endif
          endif
100    continue
123    q2(i,j) = qsum / ( lon2(i+1) - lon2(i) )
555    continue
1000 continue
     !$OMP END PARALLEL DO

  END SUBROUTINE xmap
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_input_grid
!
! !DESCRIPTION: Routine to read variables and attributes from a netCDF
!  file.  This routine was automatically generated by the Perl script
!  NcdfUtilities/perl/ncCodeRead.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE READ_INPUT_GRID( IM, JM, fileName, lon_edges, lat_sines )
!
! !USES:
!
    ! Modules for netCDF read
    USE m_netcdf_io_open
    USE m_netcdf_io_get_dimlen
    USE m_netcdf_io_read
    USE m_netcdf_io_readattr
    USE m_netcdf_io_close

    IMPLICIT NONE

#   include "netcdf.inc"
!
! !INPUT PARAMETERS:
!
    INTEGER,          INTENT(IN)  :: IM                ! # of longitudes
    INTEGER,          INTENT(IN)  :: JM                ! # of latitudes
    CHARACTER(LEN=*), INTENT(IN)  :: fileName          ! File w/ grid info
!
! !OUTPUT PARAMETERS:
!
    REAL*8,           INTENT(OUT) :: lon_edges(IM+1)   ! Lon edges [degrees]
    REAL*8,           INTENT(OUT) :: lat_sines(JM+1)   ! SIN( latitude edges )
!
! !REMARKS:
!  Created with the ncCodeRead script of the NcdfUtilities package,
!  with subsequent hand-editing.
!
! !REVISION HISTORY:
!  23 Aug 2012 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: fId                          ! netCDF file ID

    ! Arrays
    INTEGER            :: st1d(1), ct1d(1)             ! netCDF start & count

    !======================================================================
    ! Read data from file
    !======================================================================

    ! Open file for reading
    CALL Ncop_Rd( fId, TRIM( fileName ) )

    ! Read lon_edges from file
    st1d = (/ 1    /)
    ct1d = (/ IM+1 /)
    CALL NcRd( lon_edges, fId,  "lon_edges", st1d, ct1d )

    ! Read lat_sines from file
    st1d = (/ 1    /)
    ct1d = (/ JM+1 /)
    CALL NcRd( lat_sines, fId,  "lat_sines", st1d, ct1d )

    ! Close netCDF file
    CALL NcCl( fId )

  END SUBROUTINE READ_INPUT_GRID
!EOC
END MODULE REGRID_A2A_MOD
