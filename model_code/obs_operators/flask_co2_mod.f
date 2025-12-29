!$Id: flask_co2_mod.f,v 1.2 2011/02/23 00:08:48 daven Exp $
      MODULE FLASK_CO2_MOD

      !
      ! Observation operator for surface flask observations
      !
  
      IMPLICIT NONE 

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================   
 
      ! Parameters
      INTEGER, PARAMETER           :: MAXLEV = 20
      INTEGER, PARAMETER           :: nlev=1
      INTEGER, PARAMETER           :: MAXGOS = 40000


      ! Record to store data from each GOS obs
      TYPE FLASK_CO2_OBS 
         REAL*8                       :: QF(1)
         INTEGER                       :: DAY(1)
         REAL*8                       :: LAT(1)
         REAL*8                       :: LON(1)
         REAL*8                       :: TIME(1)
         REAL*8                       :: error(1)
         REAL*8                       :: CO2(1)
         REAL*8                       :: height(1)
         character(100)               :: id(1)

      ENDTYPE FLASK_CO2_OBS  

      TYPE(FLASK_CO2_OBS)                          :: FLASK(MAXGOS)

      ! IDTCO2 isn't defined in tracerid_mod because people just assume 
      ! it is one. Define it here for now as a temporary patch. 
      INTEGER, PARAMETER :: IDTCO2   = 1 
      ! Same thing for TCVV(IDTCO2) 
      REAL*8,  PARAMETER :: TCVV_CO2 = 28.97d0  / 44d0 

      CONTAINS
!------------------------------------------------------------------------------

      SUBROUTINE READ_FLASK_CO2_OBS( YYYYMMDD, NFLASK )
!
!******************************************************************************
!  Subroutine READ_GOS_CO2_OBS reads the file and passes back info contained
!  therein. (dkh, 10/12/10) 
! 
!  Based on READ_TES_NH3 OBS (dkh, 04/26/10) 
! 
!  Arguments as Input:
!  ============================================================================
!  (1 ) YYYYMMDD    INTEGER : Current year-month-day
!
!  Arguments as Output: 
!  ============================================================================
!  (1 ) NGOS      (INTEGER) : Number of GOS retrievals for current day 
!
!  Module variable as Output: 
!  ============================================================================
!  (1 ) GOS    (GOS_CO2_OBS) : CO2 retrieval for current day 
!     
!  NOTES:
!  (1 ) Add calculation of S_OER_INV, though eventually we probably want to 
!        do this offline. (dkh, 05/04/10) 
!******************************************************************************
!
      ! Reference to f90 modules
      USE DIRECTORY_MOD,          ONLY : DATA_DIR, OBS_DIR_FLASK
      USE DIRECTORY_MOD,          ONLY : FFILE_PREFIX
      USE NETCDF 
      USE TIME_MOD              
      IMPLICIT NONE



      ! Arguments
      INTEGER,            INTENT(IN)  :: YYYYMMDD
    
      ! local variables 
      INTEGER                         :: FID
      INTEGER                         :: NFLASK 
      INTEGER                         :: START0(1), COUNT0(1)
      INTEGER                         :: START1(2), COUNT1(2)
      INTEGER                         :: START2(3), COUNT2(3)
      INTEGER                         :: N, J
      INTEGER                         :: NT_ID
      INTEGER                         :: CO2_ID
      INTEGER                         :: OI_ID, ID_ID, AS_ID
      INTEGER :: LA_ID, LO_ID, TM_ID, LV_ID, QF_ID, height_ID
      CHARACTER(LEN=255)              :: READ_FILENAME

      REAL*8, PARAMETER               :: FILL = -999.0D0
      REAL*8, PARAMETER               :: TOL  = 1d-04
      INTEGER                         :: I, II, III, status

      real*4              ::temp1(1),temp2(MAXLEV),temp4(6)
	integer :: temp3(6)
      CHARACTER(LEN=5)                :: TMP
      CHARACTER(46)                  :: obsdir
      logical :: readFile, it_exists
	integer :: istat,iobs
	character(3):: stn
	integer :: year, mon, day, hr, sec
	real :: co2, lat, lon, height
	character(255):: tempdir

      !Observation filepath and filenames
      CHARACTER(LEN=255)             :: FILE_SUFFIX
      FILE_SUFFIX='.nc'

      !=================================================================
      ! READ_GOS_CO2_OBS begins here!
      !=================================================================
      MON=GET_MONTH()
      YEAR=GET_YEAR()
      day=get_day()
	write(6,*)TRIM(OBS_DIR_FLASK),TRIM(FFILE_PREFIX)
	write(READ_FILENAME,102)TRIM(OBS_DIR_FLASK),TRIM(FFILE_PREFIX)
     &,year,mon,day,TRIM(FILE_SUFFIX)
102    format(a,a,i4.4,i2.2,i2.2,a)

      WRITE(6,*) '    - READ_FLASK_CO2_OBS: reading file: ', 
     &   READ_FILENAME

      INQUIRE( FILE=TRIM( READ_FILENAME ), EXIST=IT_EXISTS )
      if (it_exists == .false.) then
                print *,'No observation data files for this date'
                NFLASK=0
                return
      endif
      ! Open file and assign file id (FID)
      CALL CHECK( NF90_OPEN( READ_FILENAME, NF90_NOWRITE, FID ), 0 )


      ! Get data record IDs
      !--------------------------------
      CALL CHECK( NF90_INQ_DIMID( FID, "obs",        NT_ID),  102 )
      CALL CHECK( NF90_INQ_VARID(FID,"value",CO2_ID), 107 )
      CALL CHECK( NF90_INQ_VARID( FID, "altitude", Height_ID ), 105 )
      CALL CHECK( NF90_INQ_VARID( FID, "latitude",        LA_ID ), 108 )
      CALL CHECK( NF90_INQ_VARID( FID, "longitude",       LO_ID ), 109 )
      CALL CHECK( NF90_INQ_VARID( FID, "time_components",TM_ID ), 111 )
      CALL CHECK(NF90_INQ_VARID(FID,"CT_assim",QF_ID),112)
      CALL CHECK( NF90_INQ_VARID( FID, "CT_MDM", OI_ID),115)
      CALL CHECK( NF90_INQ_VARID( FID, "obspack_id",ID_ID),115)


      ! READ number of retrievals, NGOS
      CALL CHECK( NF90_INQUIRE_DIMENSION( FID, NT_ID, len=NFLASK),202 )

      ! define record size
      START0 = (/1/)
      COUNT0 = (/1/)
      ! loop over records
      DO N = 1, NFLASK

         ! Update starting index
         START0(1) = N
         CALL CHECK( NF90_GET_VAR  ( FID, QF_ID,
     &      FLASk(N)%QF,          START0, COUNT0 ), 301 )

         CALL CHECK( NF90_GET_VAR  ( FID, LA_ID,
     &      FLASk(N)%LAT,          START0, COUNT0 ), 301 )

         CALL CHECK( NF90_GET_VAR  ( FID, LO_ID,
     &      FLASk(N)%LON,          START0, COUNT0 ), 302 )

         CALL CHECK( NF90_GET_VAR  ( FID, height_ID,
     &      FLASk(N)%height,          START0, COUNT0 ), 303 )
        CALL CHECK( NF90_GET_VAR  ( FID, ID_ID,
     &      FLASk(N)%ID,start=(/1,n/), count=(/100,1/)),304 )


         ! READ the error statistics of the apriori state
         CALL CHECK( NF90_GET_VAR  ( FID, OI_ID,
     &      temp1,START0, COUNT0 ), 307 )

            ! SCOT: some unassimilated observations have NA error
            ! values,and that causes GEOS-Chem to crash.
            ! Makes it difficult to calculate co-samples.
            ! I'm going to temporarily set the error to 1.
            !FLASk(N)%error=1
            FLASk(N)%error=temp1

         ! READ the error statistics of the xCO2 state
         CALL CHECK( NF90_GET_VAR  ( FID, CO2_ID,
     &      temp1, START0, COUNT0 ), 308 )
            FLASk(N)%CO2=real(temp1(1),kind=8)
         ! Assume that hour is in the 4th position of the time vector
         ! Assume that minute is in the 5th position of the time vector
         ! Assume that seconds are in the 6th position
         CALL CHECK( NF90_GET_VAR( FID, TM_ID,
     & temp4,    start=(/1,n/), count=(/6,1/)),308)
          FLASk(N)%TIME(1)=temp4(4)+temp4(5)/60.0+temp4(6)/(60*60)

      ENDDO



      END SUBROUTINE READ_FLASK_CO2_OBS
!------------------------------------------------------------------------------

      SUBROUTINE CHECK( STATUS, LOCATION )
!
!******************************************************************************
!  Subroutine CHECK checks the status of calls to netCDF libraries routines
!  (dkh, 02/15/09) 
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) STATUS    (INTEGER) : Completion status of netCDF library call    
!  (2 ) LOCATION  (INTEGER) : Location at which netCDF library call was made   
!     
!     
!  NOTES:
!
!******************************************************************************
!
      ! Reference to f90 modules 
      USE ERROR_MOD,    ONLY  : ERROR_STOP
      USE NETCDF 
 
      ! Arguments
      INTEGER, INTENT(IN)    :: STATUS 
      INTEGER, INTENT(IN)    :: LOCATION
    
      !=================================================================
      ! CHECK begins here!
      !=================================================================

      IF ( STATUS /= NF90_NOERR ) THEN 
        WRITE(6,*) TRIM( NF90_STRERROR( STATUS ) )
        WRITE(6,*) 'At location = ', LOCATION 
        CALL ERROR_STOP('netCDF error', 'tes_nh3_mod')
      ENDIF 

      ! Return to calling program
      END SUBROUTINE CHECK

!------------------------------------------------------------------------------

      SUBROUTINE CALC_FLASK_CO2_FORCE( COST_FUNC )
!
!******************************************************************************
!  Subroutine CALC_GOS_CO2_FORCE calculates the adjoint forcing from the GOSAT
!  CO2 observations and updates the cost function. (dkh, 10/12/10) 
! 
!
!  Arguments as Input/Output:
!  ============================================================================
!  (1 ) COST_FUNC (REAL*8) : Cost function                        [unitless]
!     
!     
!  NOTES:
!  
!******************************************************************************
!
      ! Reference to f90 modules
      USE ADJ_ARRAYS_MOD,     ONLY : STT_ADJ
      USE ADJ_ARRAYS_MOD,     ONLY : N_CALC
      USE ADJ_ARRAYS_MOD,     ONLY : EXPAND_NAME
      USE ADJ_ARRAYS_MOD,     ONLY : FLASKQFLAG
      USE CHECKPT_MOD,        ONLY : CHK_STT
      USE COMODE_MOD,         ONLY : CSPEC, JLOP
      USE DAO_MOD,            ONLY : AD
      USE DAO_MOD,            ONLY : AIRDEN
      USE DAO_MOD,            ONLY : BXHEIGHT
      USE DIRECTORY_ADJ_MOD,  ONLY : DIAGADJ_DIR
      USE DIRECTORY_ADJ_MOD,  ONLY : TOPOFN
      USE GRID_MOD,           ONLY : GET_IJ, XMID, YMID
      USE PRESSURE_MOD,       ONLY : GET_PCENTER, GET_PEDGE
      USE TIME_MOD,           ONLY : GET_NYMD,    GET_NHMS, GET_DAY
      USE TIME_MOD,           ONLY : GET_TS_CHEM, YMD_EXTRACT
      USE TRACER_MOD,         ONLY : XNUMOLAIR, STT
      USE TROPOPAUSE_MOD,     ONLY : ITS_IN_THE_TROP
      USE interpolate
      use netcdf


#     include      "CMN_SIZE"      ! Size params

      ! Arguments
      REAL*8, INTENT(INOUT)       :: COST_FUNC
   
      ! Local variables 
      INTEGER                     :: NTSTART, NTSTOP, NT 
      INTEGER                     :: IIJJ(2), I,      J
      INTEGER                     :: L,       LL,     LGOS
      INTEGER                     :: JLOOP
      REAL*8                      :: GC_PRES(LLPAR)
      REAL*8                      :: GC_CO2_NATIVE(LLPAR)
      REAL*8                      :: GC_CO2(MAXLEV)
      REAL*8                      :: GC_PSURF
      REAL*8                      :: MAP(LLPAR,MAXLEV)
      REAL*8                      :: CO2_HAT(MAXLEV)
      REAL*8                      :: CO2_PERT(MAXLEV)
      REAL*8                      :: FORCE(MAXLEV)
      REAL*8                      :: DIFF(MAXLEV)
      REAL*8                      :: NEW_COST(MAXGOS)
      REAL*8                      :: OLD_COST
      REAL*8, SAVE                :: TIME_FRAC(MAXGOS)
      INTEGER,SAVE                :: NFLASK
      REAL*8                      :: CO2CORRECTED

      REAL*8                      :: GC_CO2_NATIVE_ADJ(LLPAR)
      REAL*8                      :: CO2_HAT_ADJ(MAXLEV)
      REAL*8                      :: CO2_PERT_ADJ(MAXLEV)
      REAL*8                      :: GC_CO2_ADJ(MAXLEV)
      REAL*8                      :: DIFF_ADJ(MAXLEV)
      REAL*8                :: temp1(1,MAXLEV)
      REAL*8                 :: temp(MAXLEV)

	real :: tindex	
	integer :: itindex1, itindex2, ii, jj
	real :: temp2(IIPAR, LLPAR),ptemp2(IIPAR, LLPAR)
	real :: atemp2(IIPAR, LLPAR)
	real :: temp3(LLPAR), ptemp3(LLPAR), atemp3(LLPAR)
	real :: temp_lat, LONF(IIPAR), LATG(JJPAR)
	real :: ttemp4, ttemp3(IIPAR)
	real :: map_X(IIPAR), map_Y(JJPAR)
	  REAL:: TEMP_ADJ(JJPAR, LLPAR)
     
   
      LOGICAL, SAVE               :: FIRST = .TRUE. 
      INTEGER                     :: IOS
      CHARACTER(LEN=255)          :: FILENAME
      LOGICAL                     :: THERE
      !
      ! Junjie (random number)
      !
      REAL*8                      :: rand(MAXGOS), GC_HH_FRAC
      REAL*8                      :: xco2_before(1)
      REAL*4                      :: CO2_pseudo(MAXLEV)
      INTEGER                 :: YYYYMMDD, HH, MM, SS, DD
      character(35)                 :: filenameobs=
     &  'residuals_insitu.csv'
      character(35)                 :: filenameobs2=
     &  'residuals_situ_identifier.csv'
      REAL*8                      :: outvec(7)
      real :: PHIS(IIPAR, JJPAR)
        integer :: status, ncid, varid, kk



      !=================================================================
      ! CALC_GOS_CO2_FORCE begins here!
      !=================================================================

      print*, '     - CALC_FLASK_CO2_FORCE '
    
      ! Reset 
      NEW_COST = 0D0 

      ! Create files for the residuals
         OPEN(unit = 34, access = "sequential", action = "write",
     &       position='append',file=TRIM(filenameobs),
     &       form="formatted")

         OPEN(unit = 35, access = "sequential", action = "write",
     &       position='append',file=TRIM(filenameobs2),
     &       form="formatted")

      ! Save a value of the cost function first
      OLD_COST = COST_FUNC

      ! Check if it is the last hour of a day 
      IF ( GET_NHMS() == 236000 - GET_TS_CHEM() * 100 ) THEN 
      !IF ( GET_NHMS() == 0 ) THEN 
 
         ! Read the GOS CO2 file for this day 
	write(*,*)'GET_NYMD()=', GET_NYMD()
         CALL READ_FLASK_CO2_OBS( GET_NYMD(), NFLASK )
 
         ! TIME is YYYYMMDD.frac-of-day.  Subtract date and save just time fraction
         ! Don't need to adjust this for GOSAT CO2, for which TIME is already 
         ! just the time fraction. 
         TIME_FRAC(1:NFLASK) = FLASK(1:NFLASK)%TIME(1)/24.0

      ENDIF 

      CALL YMD_EXTRACT(GET_NHMS(), HH, MM, SS)
      GC_HH_FRAC=real(HH,kind=8)/24.0
	DD=GET_DAY()

      print*, ' GC_HH_FRAC = ', GC_HH_FRAC

      !
      ! Read in surface geopotential height and then convert it to height by dividing it by 9.8
      !
      INQUIRE( FILE=TRIM( TOPOFN ), EXIST=THERE)
      if (THERE == .false.) then
                print *,'ERROR: CANNOT FIND MERRA2.20150101.CN.4x5.NC'
                print *,'CHECK flask_co2_mod.f'
                stop
      endif
	status=nf90_open(trim(TOPOFN), nf90_nowrite, ncid)
	status=nf90_inq_varid(ncid, "PHIS", varid)
	status=nf90_get_var(ncid, varid, phis)
	status=nf90_close(ncid)
	phis=phis/9.8

      DO NT=1, NFLASK	

        ! For the OCO-2 MIP, we'll only use observations with a quality
        ! flag of 1. (data with a quality flag of 2 are withheld for
        ! model evaluation). A qualify flag of 0 means that the
        ! observations are not suitable for assimilation in an inverse
        ! model. 

	IF((GC_HH_FRAC .ge. (TIME_FRAC(NT)-0.5/24.0)) 
     & .and. (GC_HH_FRAC .lt. (TIME_FRAC(NT)+0.5/24.0)))then

        ! Input.geos contains a variable called FLASKQFLAG. If the value
        ! is 1, then only use observations with a QFlag of 1. If the
        ! value is 0, then use all observations.
        IF((FLASKQFLAG.eq.1 .AND. FLASK(NT)%QF(1).eq.1) .OR.
     &      FLASKQFLAG.eq.0)then

!          IF(FLASK(NT)%QF(1).eq.1)then
!	  IF(FLASK(NT)%QF(1).ge.1)then
!         IF(FLASK(NT)%QF(1).ge.0)then
         !print*, '     - CALC_GOS_CO2_FORCE: analyzing record ', NT 
         !print*, 'TIME_FRAC(NT)',TIME_FRAC(NT)
         !print*, TRIM(FLASK(NT)%ID(1))

	LGOS=1
         ! For safety, initialize these up to LLGOS
         GC_CO2(:)       = 0d0 
         MAP(:,:)        = 0d0 
         CO2_HAT_ADJ(:)  = 0d0 
         FORCE(:)        = 0d0 


         ! INTERpolate in x and y direction
         !IIJJ  = GET_IJ( REAL(FLASK(NT)%LON(1),4),
         !& REAL(FLASK(NT)%LAT(1),4))
             !I     = IIJJ(1)
             !J     = IIJJ(2)

	 temp_lat=flask(NT)%LAT(1)

	 LATG=real(YMID, kind=4)
	 LONF=real(XMID, kind=4)
	MAP_Y(1:JJPAR)=GET_INTMAP_Y(JJPAR, LATG, temp_lat)

	 temp_lat=flask(NT)%LON(1)
	MAP_X(1:IIPAR)=GET_INTMAP_X(IIPAR, LONF, temp_lat)
	
	temp2=0.0;ptemp2=0.0;atemp2=0.0
	DO L=1, LLPAR
	DO JJ=1,JJPAR
		DO ii=1, IIPAR
	  temp2(ii,L)=CHK_STT(II,jj, L, IDTCO2)*map_Y(jj)+temp2(ii,L)	
	  ptemp2(ii,L)=BXHEIGHT(II,jj, L)*map_Y(jj)+ptemp2(ii,L)	
	  Atemp2(ii,L)=AD(II,jj, L)*map_Y(jj)+atemp2(ii,L)	
		ENDDO
	ENDDO
	ENDDO
	temp3=0.0;ptemp3=0.0;Atemp3=0.0
	DO L=1, LLPAR
		DO ii=1, IIPAR
	  temp3(L)=temp2(ii,L)*map_X(ii)+temp3(L)	
	  ptemp3(L)=ptemp2(II, L)*map_X(ii)+ptemp3(L)	
	  Atemp3(L)=atemp2(II, L)*map_X(ii)+atemp3(L)	
		ENDDO
	ENDDO
	ttemp3=0.0;ttemp4=0.0
	DO JJ=1, JJPAR
		DO ii=1, IIPAR
	  ttemp3(ii)=PHIS(ii,JJ)*map_Y(JJ)+ttemp3(ii)	
		ENDDO
	ENDDO
	DO II=1, IIPAR
	    ttemp4=ttemp4+ttemp3(II)*map_X(ii)
	ENDDO

         ! Get GC surface pressure (mbar)
         GC_PSURF = ttemp4

         !
         ! Calculate the GEOS-Chem height at I, J
         !
         GC_PRES(1)=ttemp4
         DO L=2, LLPAR
             GC_PRES(L)=ptemp3(L)+GC_PRES(L-1)
         ENDDO
 
         ! Calculate the interpolation weight matrix 
         MAP(1:LLPAR,1:1) 
     &      = GET_INTMAP_HEIGHT( LLPAR, GC_PRES(:),   GC_PSURF, 
     &                    1,  REAL(FLASK(NT)%HEIGHT(1:1),kind=8)
     & , GC_PSURF  )

	!write(*,*)'map=', map(:,1)
	!write(*,*)'obs=', FLASK(NT)%lon(1), 
       ! & flask(nt)%lat, flask(nt)%height, GC_PSURF
        !    	write(*,*)'model height=', GC_PRES(:)

         ! Get CO2 values at native model resolution
         !GC_CO2_NATIVE(:) = STT(I,J,:,IDTCO2)
         GC_CO2_NATIVE(:) = temp3(:)

 
         ! Convert from kg/box to v/v
         !GC_CO2_NATIVE(:) = GC_CO2_NATIVE(:) * TCVV_CO2
         !&                    / (AD(I,J,:) )
         GC_CO2_NATIVE(:) = GC_CO2_NATIVE(:) * TCVV_CO2
     &                    / (Atemp3(:) )

	!write(*,*)'model CO2=', GC_CO2_NATIVE(:)

         ! Interpolate GC CO2 column to TES grid 
            GC_CO2(1) = 0d0 
            DO L = 1, LLPAR 
               GC_CO2(1) = GC_CO2(1) 
     &                    + MAP(L,1) * GC_CO2_NATIVE(L) 
            ENDDO


         !--------------------------------------------------------------
         ! Apply GOS observation operator
         !
         !   x_hat = x_a + A_k ( x_m - x_a ) 
         !  
         !  where  
         !    x_hat = GC modeled column as seen by TES [vmr]
         !    x_a   = GOS apriori column               [vmr]
         !    x_m   = GC modeled column                [vmr]
         !    A_k   = GOS averaging kernel 
         !--------------------------------------------------------------

         ! x_m - x_a
         DO L = 1,  1
           CO2_PERT(L) = GC_CO2(L) - FLASK(NT)%CO2(L)
	!   write(*,*)'model=', GC_CO2(L), 'flask=', FLASK(NT)%CO2(L)
         ENDDO
     
          
        DO L = 1,  nlev
            DIFF(L)=CO2_PERT(L)
         ENDDO
      ! endif

	
        ! From Scot: Originally, this code would only include obs
        ! with residuals less than 6 ppb. That created a problem
        ! For the inverse model -- LBFGS would make crazy flux estimtes
        ! with the purpose of exploring the uncertainty space, and this
        ! script would then dump all in situ observations from the 
        ! optimization. LBFGS gets thrown off if the observations
        ! that get included or excluded is changing with every
        ! iteration.
	IF(abs(DIFF(1))*1e6.lt.1e10)THEN

          outvec = (/real(GET_NYMD(),8),FLASK(NT)%LAT(1),
     &                    FLASK(NT)%LON(1),
     &                    GC_CO2(1)*1e6,FLASK(NT)%CO2(1)*1e6,
     &                    FLASK(NT)%error(1)*1e6,FLASK(NT)%QF(1)/)

          ! Formatting for CSV
 121  FORMAT(1x, *(g0, ", "))

          ! Write residuals to file
          WRITE(34,121) outvec
          WRITE(35,121) TRIM(FLASK(NT)%ID(1))

         ! Calculate 1/2 * DIFF^T * S_{obs}^{-1} * DIFF
         !DO L = 1, LGOS
         DO L = 1, nlev
               FORCE(L)     = 0d0
               !FORCE(L)  = FORCE(L) + 1.0/(1e-6*1e-6) * DIFF(L)
               FORCE(L) =FORCE(L) + 
     &1.0/(FLASK(NT)%error(1)*flask(NT)%error(1)) * DIFF(L)
               ! FORCE(L)  = FORCE(L) +  DIFF(L)
               NEW_COST(NT) = NEW_COST(NT) + 0.5d0 * DIFF(L) * FORCE(L)
         ENDDO

         !--------------------------------------------------------------
         ! Begin adjoint calculations
         !--------------------------------------------------------------

         ! The adjoint forcing  is S_{obs}^{-1} * DIFF = FORCE
         DIFF_ADJ(:) = FORCE(:)

         ! Adjoint of difference
         !DO L = 1, LGOS
         DO L = 1,  nlev
            CO2_HAT_ADJ(L) =  DIFF_ADJ(L)
         ENDDO
         DO L = 1,  nlev
            CO2_PERT_ADJ(L) =  DIFF_ADJ(L)
         ENDDO

         ! Adjoint of x_m - x_a
         DO L = 1, 1
           ! fwd code:
           !CO2_PERT(L) = GC_CO2(L) - GOS(NT)%PRIOR(L)
           ! adj code:
           GC_CO2_ADJ(L) = CO2_PERT_ADJ(L)
         ENDDO
 
         ! adjoint of interpolation
         DO L  = 1, LLPAR
            GC_CO2_NATIVE_ADJ(L) = 0d0
            DO LL = 1, LGOS
               GC_CO2_NATIVE_ADJ(L) = GC_CO2_NATIVE_ADJ(L)
     &                              + MAP(L,LL) * GC_CO2_ADJ(LL)
            ENDDO
         ENDDO


         ! Adjoint of unit conversion
         GC_CO2_NATIVE_ADJ(:) = GC_CO2_NATIVE_ADJ(:) * TCVV_CO2
      !&                        / (AD(I,J,:) )
     &                        / (atemp3(:) )



        ! Pass adjoint back to adjoint tracer array
         !STT_ADJ(I,J,:,IDTCO2)  = STT_ADJ(I,J,:,IDTCO2)
         !&                          + GC_CO2_NATIVE_ADJ(:)
	TEMP_ADJ=0.0
	DO LL=1, LLPAR
	    DO JJ=1, JJPAR
         TEMP_ADJ(JJ,LL)  = temp_ADJ(JJ,LL)
     &                          + GC_CO2_NATIVE_ADJ(LL)*MAP_Y(JJ)
	    ENDDO
	ENDDO
	DO LL=1, LLPAR
	    DO JJ=1, JJPAR
		DO II=1, IIPAR
         STT_ADJ(II,JJ,LL,IDTCO2)  = STT_ADJ(II,JJ,LL,IDTCO2)
     &                          + TEMP_ADJ(JJ,LL)*MAP_X(II)
		ENDDO
	    ENDDO
	ENDDO

 103  FORMAT(1X,d14.6,1X,d14.6,1X,d14.6)
 104  FORMAT(1X,d14.6)

!         WRITE(103,110) ( 1d6 * FLASK(NT)%CO2(LL),   LL=LGOS,1,-1)
!         WRITE(105,110) ( DIFF(LL),                LL=LGOS,1,-1)
!         WRITE(106,112) ( FORCE(LL),               LL=LGOS,1,-1)
!         WRITE(107,111) NT, LGOS
!         WRITE(108,112) ( CO2_PERT_ADJ(LL),        LL=LGOS,1,-1)
!         WRITE(109,112) ( GC_CO2_ADJ(LL),          LL=LGOS,1,-1)
!         WRITE(110,110) ( 1d6 * CO2_HAT(LL),       LL=LGOS,1,-1)
!         WRITE(111,110) ( GC_PRES(L),              L=LLPAR,1,-1)
!         WRITE(112,110) ( 1d6 * GC_CO2_NATIVE(L),  L=LLPAR,1,-1)
!         WRITE(113,110) ( 1d6 * GC_CO2(LL),        LL=LGOS,1,-1)
! 110     FORMAT(F18.6,1X)
! 111     FORMAT(i4,1X,i4,1x)
! 112     FORMAT(D14.6,1X)

      ! Update cost function 
         COST_FUNC = COST_FUNC + NEW_COST(NT)
	ENDIF !IF (quality control)
	ENDIF ! END quality control

	ENDIF


      ENDDO  ! NT
!!$OMP END PARALLEL DO

      ! Update cost function 
      !COST_FUNC = COST_FUNC + SUM(NEW_COST(NTSTOP:NTSTART))

      IF ( FIRST ) FIRST = .FALSE. 

      ! Close the model-data residual csv file
      CLOSE(34)
      CLOSE(35)

      print*, ' Updated value of COST_FUNC = ', COST_FUNC 
      print*, 'GOSAT contribution          = ', OLD_COST
      print*, 'IN SITU contribution        = ', COST_FUNC - OLD_COST  

      ! Return to calling program
      END SUBROUTINE CALC_FLASK_CO2_FORCE


      SUBROUTINE GET_NT_RANGE( NTES, HHMMSS, TIME_FRAC, NTSTART, NTSTOP)
!
!******************************************************************************
!  Subroutine GET_NT_RANGE retuns the range of retrieval records for the 
!  current model hour 
! 
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) NTES   (INTEGER) : Number of TES retrievals in this day 
!  (2 ) HHMMSS (INTEGER) : Current model time 
!  (3 ) TIME_FRAC (REAL) : Vector of times (frac-of-day) for the TES retrievals
!     
!  Arguments as Output:
!  ============================================================================
!  (1 ) NTSTART (INTEGER) : TES record number at which to start
!  (1 ) NTSTOP  (INTEGER) : TES record number at which to stop
!     
!  NOTES:
!
!******************************************************************************
!
      ! Reference to f90 modules
      USE ERROR_MOD,    ONLY : ERROR_STOP
      USE TIME_MOD,     ONLY : YMD_EXTRACT

      ! Arguments
      INTEGER, INTENT(IN)   :: NTES
      INTEGER, INTENT(IN)   :: HHMMSS
      REAL*8,  INTENT(IN)   :: TIME_FRAC(NTES)
      INTEGER, INTENT(OUT)  :: NTSTART
      INTEGER, INTENT(OUT)  :: NTSTOP
    
      ! Local variables 
      INTEGER, SAVE         :: NTSAVE
      LOGICAL               :: FOUND_ALL_RECORDS 
      INTEGER               :: NTEST
      INTEGER               :: HH, MM, SS
      REAL*8                :: GC_HH_FRAC
      REAL*8                :: H1_FRAC

      !=================================================================
      ! GET_NT_RANGE begins here!
      !=================================================================


      ! Initialize 
      FOUND_ALL_RECORDS  = .FALSE. 
      NTSTART            = 0
      NTSTOP             = 0

      ! set NTSAVE to NTES every time we start with a new file
      IF ( HHMMSS == 230000 ) NTSAVE = NTES - 100
   
      print*, ' co2 hack : skip lat 100 records, where out of order' 
      print*, ' co2 hack : skip lat 100 records, where out of order' 
      print*, ' co2 hack : skip lat 100 records, where out of order' 
      print*, ' co2 hack : skip lat 100 records, where out of order' 
      print*, ' co2 hack : skip lat 100 records, where out of order' 
  


      print*, ' GET_NT_RANGE for ', HHMMSS
      print*, ' NTSAVE ', NTSAVE
      print*, ' NTES   ', NTES
   
      CALL YMD_EXTRACT( HHMMSS, HH, MM, SS )


      ! Convert HH from hour to fraction of day 
      GC_HH_FRAC = REAL(HH,8) / 24d0 
 
      ! one hour as a fraction of day 
      H1_FRAC    = 1d0 / 24d0 


      ! dkh debug
      print*, ' co2 time frac = ', TIME_FRAC

      ! All records have been read already 
      IF ( NTSAVE == 0 ) THEN 

         print*, 'All records have been read already '
         RETURN 

      ! No records reached yet
      ELSEIF ( TIME_FRAC(NTSAVE) + H1_FRAC < GC_HH_FRAC ) THEN 
           
      
         print*, 'No records reached yet'
         RETURN

      !
      ELSEIF ( TIME_FRAC(NTSAVE) + H1_FRAC >=  GC_HH_FRAC ) THEN 
      
         ! Starting record found
         NTSTART = NTSAVE   

         print*, ' Starting : TIME_FRAC(NTSTART) ', 
     &               TIME_FRAC(NTSTART), NTSTART
 
         ! Now search forward to find stopping record
         NTEST = NTSTART

         DO WHILE ( FOUND_ALL_RECORDS == .FALSE. ) 
              
            ! Advance to the next record
            NTEST = NTEST - 1  
           
            ! Stop if we reach the earliest available record 
            IF ( NTEST == 0 ) THEN 
           
               NTSTOP            = NTEST + 1
               FOUND_ALL_RECORDS = .TRUE.

               print*, ' Records found '
               print*, ' NTSTART, NTSTOP = ', NTSTART, NTSTOP

               ! Reset NTSAVE 
               NTSAVE = NTEST

            ! When the combined test date rounded up to the nearest
            ! half hour is smaller than the current model date, the 
            ! stopping record has been passed. 
            ELSEIF (  TIME_FRAC(NTEST) + H1_FRAC <  GC_HH_FRAC ) THEN
          
               print*, ' Testing : TIME_FRAC ', 
     &                  TIME_FRAC(NTEST), NTEST
 
               NTSTOP            = NTEST + 1 
               FOUND_ALL_RECORDS = .TRUE. 

               print*, ' Records found '
               print*, ' NTSTART, NTSTOP = ', NTSTART, NTSTOP
                  
               ! Reset NTSAVE 
               NTSAVE = NTEST 

            ELSE 
               print*, ' still looking ', NTEST 
                  
            ENDIF 
                 
         ENDDO 
 
      ELSE

         CALL ERROR_STOP('problem', 'GET_NT_RANGE' ) 

      ENDIF 

      ! Return to calling program
      END SUBROUTINE GET_NT_RANGE

!------------------------------------------------------------------------------

      FUNCTION GET_INTMAP( LGC_TOP, GC_PRESC, GC_SURFP,
     &                     LTM_TOP, TM_PRESC, TM_SURFP  )
     *         RESULT      ( HINTERPZ )
!
!******************************************************************************
!  Function GET_INTMAP linearly interpolates column quatities
!   based upon the centered (average) pressue levels. 
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) LGC_TOP (TYPE) : Description                          [unit]
!  (2 ) GC_PRES (TYPE) : Description                          [unit]
!  (3 ) GC_SURFP(TYPE) : Description                          [unit]
!  (4 ) LTM_TOP (TYPE) : Description                          [unit]
!  (5 ) TM_PRES (TYPE) : Description                          [unit]
!  (6 ) TM_SURFP(TYPE) : Description                          [unit]
!     
!  Arguments as Output:
!  ============================================================================
!  (1 ) HINTERPZ (TYPE) : Description                          [unit]
!     
!  NOTES:
!  (1 ) Based on the GET_HINTERPZ_2 routine I wrote for read_sciano2_mod. 
!
!******************************************************************************
!
      ! Reference to f90 modules
      USE ERROR_MOD,     ONLY : ERROR_STOP
      USE PRESSURE_MOD,  ONLY : GET_BP

      ! Arguments
      INTEGER            :: LGC_TOP, LTM_TOP
      REAL*8             :: GC_PRESC(LGC_TOP)
      REAL*8             :: TM_PRESC(LTM_TOP) 
      REAL*8             :: GC_SURFP
      REAL*8             :: TM_SURFP
 
      ! Return value 
      REAL*8             :: HINTERPZ(LGC_TOP, LTM_TOP)

      ! Local variables 
      INTEGER  :: LGC, LTM
      REAL*8   :: DIFF, DELTA_SURFP
      REAL*8   :: LOW, HI

      !=================================================================
      ! GET_HINTERPZ_2 begins here!
      !=================================================================

      HINTERPZ(:,:) = 0D0 
  
!      ! Rescale GC grid according to TM surface pressure
!!         p1_A =     (a1 + b1 (ps_A - PTOP))
!!         p2_A =     (a2 + b2 (ps_A - PTOP))
!!         p1_B =     (a + b (ps_B - PTOP))
!!         p2_B =    *(a + b (ps_B - PTOP))
!!         pc_A = 0.5(a1+a2 +(b1+b2)*(ps_A - PTOP))
!!         pc_B = 0.5(a1+a2 +(b1+b2)*(ps_B - PTOP))
!!         pc_B - pc_A = 0.5(b1_b2)(ps_B-ps_A)
!!         pc_B = 0.5(b1_b2)(ps_B-ps_A) + pc_A
!      DELTA_SURFP   = 0.5d0 * ( TM_SURFP -GC_SURFP )
!
!      DO LGC = 1, LGC_TOP
!         GC_PRESC(LGC) = ( GET_BP(LGC) + GET_BP(LGC+1))
!     &               * DELTA_SURFP + GC_PRESC(LGC)
!         IF (GC_PRESC(LGC) < 0) THEN 
!            CALL ERROR_STOP( 'highly unlikey', 
!     &                       'read_sciano2_mod.f')
!         ENDIF 
!
!      ENDDO 
      

      ! Loop over each pressure level of TM grid
      DO LTM = 1, LTM_TOP
 
         ! Find the levels from GC that bracket level LTM
         DO LGC = 1, LGC_TOP - 1

            LOW = GC_PRESC(LGC+1)
            HI  = GC_PRESC(LGC)
            IF (LGC == 0) HI = TM_SURFP

            ! Linearly interpolate value on the LTM grid 
            IF ( TM_PRESC(LTM) <= HI .and. 
     &           TM_PRESC(LTM)  > LOW) THEN 

               DIFF                = HI - LOW  
               HINTERPZ(LGC+1,LTM) = ( HI - TM_PRESC(LTM)  ) / DIFF
               HINTERPZ(LGC  ,LTM) = ( TM_PRESC(LTM) - LOW ) / DIFF


            ENDIF 
 
            ! dkh debug
            !print*, 'LGC,LTM,HINT', LGC, LTM, HINTERPZ(LGC,LTM)

          ENDDO
       ENDDO

       ! Bug fix:  a more general version allows for multiples TES pressure
       ! levels to exist below the lowest GC pressure.  (dm, dkh, 09/30/10) 
       ! OLD code:
       !IF ( TM_PRESC(1) > GC_PRESC(1) ) THEN
       !   HINTERPZ(1,1)         = 1D0 
       !   HINTERPZ(2:LGC_TOP,1) = 0D0 
       !ENDIF
       ! New code:
       ! Loop over each pressure level of TM grid
       DO LTM = 1, LTM_TOP
          IF ( TM_PRESC(LTM) > GC_PRESC(1) ) THEN
             HINTERPZ(1,LTM)         = 1D0
             HINTERPZ(2:LGC_TOP,LTM) = 0D0
          ENDIF
       ENDDO

      ! Return to calling program
      END FUNCTION GET_INTMAP

!------------------------------------------------------------------------------

      FUNCTION GET_INTMAP_HEIGHT( LGC_TOP, GC_PRESC, GC_SURFP,
     &                     LTM_TOP, TM_PRESC, TM_SURFP  )
     *         RESULT      ( HINTERPZ )
!
!******************************************************************************
!  Function GET_INTMAP linearly interpolates column quatities
!   based upon the centered (average) pressue levels. 
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) LGC_TOP (TYPE) : Description                          [unit]
!  (2 ) GC_PRES (TYPE) : Description                          [unit]
!  (3 ) GC_SURFP(TYPE) : Description                          [unit]
!  (4 ) LTM_TOP (TYPE) : Description                          [unit]
!  (5 ) TM_PRES (TYPE) : Description                          [unit]
!  (6 ) TM_SURFP(TYPE) : Description                          [unit]
!     
!  Arguments as Output:
!  ============================================================================
!  (1 ) HINTERPZ (TYPE) : Description                          [unit]
!     
!  NOTES:
!  (1 ) Based on the GET_HINTERPZ_2 routine I wrote for read_sciano2_mod. 
!
!******************************************************************************
!
      ! Reference to f90 modules
      USE ERROR_MOD,     ONLY : ERROR_STOP
      USE PRESSURE_MOD,  ONLY : GET_BP

      ! Arguments
      INTEGER            :: LGC_TOP, LTM_TOP
      REAL*8             :: GC_PRESC(LGC_TOP)
      REAL*8             :: TM_PRESC(LTM_TOP) 
      REAL*8             :: GC_SURFP
      REAL*8             :: TM_SURFP
 
      ! Return value 
      REAL*8             :: HINTERPZ(LGC_TOP, LTM_TOP)

      ! Local variables 
      INTEGER  :: LGC, LTM
      REAL*8   :: DIFF, DELTA_SURFP
      REAL*8   :: LOW, HI

      !=================================================================
      ! GET_HINTERPZ_2 begins here!
      !=================================================================

      HINTERPZ(:,:) = 0D0 
  
!      ! Rescale GC grid according to TM surface pressure
!!         p1_A =     (a1 + b1 (ps_A - PTOP))
!!         p2_A =     (a2 + b2 (ps_A - PTOP))
!!         p1_B =     (a + b (ps_B - PTOP))
!!         p2_B =    *(a + b (ps_B - PTOP))
!!         pc_A = 0.5(a1+a2 +(b1+b2)*(ps_A - PTOP))
!!         pc_B = 0.5(a1+a2 +(b1+b2)*(ps_B - PTOP))
!!         pc_B - pc_A = 0.5(b1_b2)(ps_B-ps_A)
!!         pc_B = 0.5(b1_b2)(ps_B-ps_A) + pc_A
!      DELTA_SURFP   = 0.5d0 * ( TM_SURFP -GC_SURFP )
!
!      DO LGC = 1, LGC_TOP
!         GC_PRESC(LGC) = ( GET_BP(LGC) + GET_BP(LGC+1))
!     &               * DELTA_SURFP + GC_PRESC(LGC)
!         IF (GC_PRESC(LGC) < 0) THEN 
!            CALL ERROR_STOP( 'highly unlikey', 
!     &                       'read_sciano2_mod.f')
!         ENDIF 
!
!      ENDDO 
      

      ! Loop over each pressure level of TM grid
      DO LTM = 1, LTM_TOP
 
         ! Find the levels from GC that bracket level LTM
         DO LGC = 1, LGC_TOP - 1

            LOW = GC_PRESC(LGC+1)
            HI  = GC_PRESC(LGC)
            IF (LGC == 0) HI = TM_SURFP

            ! Linearly interpolate value on the LTM grid 
            IF ( TM_PRESC(LTM) >= HI .and. 
     &           TM_PRESC(LTM)  < LOW) THEN 

               DIFF                = HI - LOW  
               HINTERPZ(LGC+1,LTM) = ( HI - TM_PRESC(LTM)  ) / DIFF
               HINTERPZ(LGC  ,LTM) = ( TM_PRESC(LTM) - LOW ) / DIFF


            ENDIF 
 
            ! dkh debug
            !print*, 'LGC,LTM,HINT', LGC, LTM, HINTERPZ(LGC,LTM)

          ENDDO
       ENDDO
       DO LTM = 1, LTM_TOP
          IF ( TM_PRESC(LTM) < GC_PRESC(1) ) THEN
             HINTERPZ(1,LTM)         = 1D0
             HINTERPZ(2:LGC_TOP,LTM) = 0D0
          ENDIF
       ENDDO

      ! Return to calling program
      END FUNCTION GET_INTMAP_HEIGHT

!!------------------------------------------------------------------------------
      FUNCTION GET_INTMAP_Y(npoint, LATG,lat)
     *         RESULT      ( map )
!
!******************************************************************************
! Interpolate in x,or y direction
!
      USE interpolate

      ! Arguments
	integer :: npoint
	real :: latg(npoint)
	real :: lat

      ! Return value
      REAL*8             :: map(npoint)

	real :: tindex
	integer:: itindex1, itindex2


	map=0
         tindex=index_interp(lat
     &, LATG, npoint)


         if(tindex.eq.0.and.(lat.le.-89))then
                  itindex1=1
                  itindex2=1
         else if(lat.gt.89)then
                  itindex1=46
                  itindex2=46
	 else
                  itindex1=floor(tindex)
                  itindex2=ceiling(tindex)
         endif

                if(itindex1.ne.itindex2)then
			map(itindex1)=(itindex2-tindex)
			map(itindex2)=(tindex-itindex1)
		
                else
			map(itindex1)=1.0
                endif
      ! Return to calling program
      END FUNCTION GET_INTMAP_Y

!!=---------------------------------------
      FUNCTION GET_INTMAP_X(npoint, LONF,lon)
     *         RESULT      ( map )
!
!******************************************************************************
! Interpolate in x,or y direction
!
      USE interpolate

      ! Arguments
        integer :: npoint
        real :: lonf(npoint)
        real :: lon

      ! Return value
      REAL*8             :: map(npoint)

        real :: tindex
        integer:: itindex1, itindex2


        map=0
         tindex=index_interp(lon
     &, LONF, npoint)

           if(lon.gt.175.0)then
		map(npoint)=(180.0-lon)/5.0
		map(1)=(lon-175.0)/5.0
           else
               tindex=index_interp(LON
     &, LONF, npoint)
                   itindex1=floor(tindex)
                   itindex2=ceiling(tindex)
                if(itindex1.ne.itindex2)then
                        map(itindex1)=(itindex2-tindex)
                        map(itindex2)=(tindex-itindex1)
                else
                        map(itindex1)=1.0
                endif
           endif


      ! Return to calling program
      END FUNCTION GET_INTMAP_X




!      SUBROUTINE MAKE_O3_FILE(  )
!!

!!------------------------------------------------------------------------------
!      SUBROUTINE MAKE_O3_FILE(  )
!!
!!******************************************************************************
!!  Subroutine MAKE_O3_FILE saves O3 profiles that correspond to time and
!!  place of TES O3 obs. (dkh, 03/01/09) 
!!
!!  Module variables as Input:
!!  ============================================================================
!!  (1 ) O3_SAVE (REAL*8) : O3 profiles                             [ppmv]
!!     
!!  NOTES:
!!
!!******************************************************************************
!!
!      ! Reference to f90 modules
!      USE BPCH2_MOD
!      USE DIRECTORY_ADJ_MOD, ONLY : ADJTMP_DIR
!      USE ERROR_MOD,        ONLY : ERROR_STOP
!      USE GRID_MOD,    ONLY : GET_XOFFSET, GET_YOFFSET
!      USE TIME_MOD,    ONLY : EXPAND_DATE
!
!#     include "CMN_SIZE"    ! Size params
!      
!      ! Local variables    
!      INTEGER              :: I, J, I0, J0, L, NT
!      CHARACTER(LEN=120)   :: FILENAME
!      REAL*4               :: DAT(1,LLPAR,MAXTES)
!      INTEGER, PARAMETER   :: IUN = 88 
!      
!      ! For binary punch file, version 2.0
!      CHARACTER(LEN=20)    :: MODELNAME
!      CHARACTER(LEN=40)    :: CATEGORY
!      CHARACTER(LEN=40)    :: UNIT     
!      CHARACTER(LEN=40)    :: RESERVED = ''
!      CHARACTER(LEN=80)    :: TITLE 
!      REAL*4               :: LONRES, LATRES
!      INTEGER, PARAMETER   :: HALFPOLAR = 1
!      INTEGER, PARAMETER   :: CENTER180 = 1
!
!      !=================================================================
!      ! MAKE_O3_FILE begins here!
!      !=================================================================
!      
!      FILENAME = TRIM( 'nh3.bpch' )
!      
!      ! Append data directory prefix
!      FILENAME = TRIM( ADJTMP_DIR ) // TRIM( FILENAME )
!      
!      ! Define variables for BINARY PUNCH FILE OUTPUT
!      TITLE    = 'O3 profile '
!      CATEGORY = 'IJ-AVE-$'
!      LONRES   = DISIZE
!      LATRES   = DJSIZE
!      UNIT     = 'ppmv'
!
!      ! Call GET_MODELNAME to return the proper model name for
!      ! the given met data being used (bmy, 6/22/00)
!      MODELNAME = GET_MODELNAME()
!
!      ! Get the nested-grid offsets
!      I0 = GET_XOFFSET( GLOBAL=.TRUE. )
!      J0 = GET_YOFFSET( GLOBAL=.TRUE. )
!
!      !=================================================================
!      ! Open the checkpoint file for output -- binary punch format
!      !=================================================================
!
!
!      WRITE( 6, 100 ) TRIM( FILENAME )
! 100  FORMAT( '     - MAKE_O3_FILE: Writing ', a )
!
!      ! Open checkpoint file for output
!      CALL OPEN_BPCH2_FOR_WRITE( IUN, FILENAME, TITLE )
!
!      ! Temporarily store data in DAT as REAL4
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( SHARED )
!!$OMP+PRIVATE( NT ) 
!      DO NT = 1, MAXTES
!
!         DAT(1,:,NT) = REAL(O3_SAVE(:,NT))
!
!      ENDDO
!!$OMP END PARALLEL DO
!
!      CALL BPCH2( IUN,       MODELNAME, LONRES,    LATRES,
!     &            HALFPOLAR, CENTER180, CATEGORY,  1,
!     &            UNIT,      1d0,       1d0,       RESERVED,
!     &            1,         LLPAR,     MAXTES,     I0+1,
!     &            J0+1,      1,         DAT )
!
!      ! Close file
!      CLOSE( IUN )        
!
!      print*, ' O3_SAVE sum write = ', SUM(O3_SAVE(:,:))
!
!      ! Return to calling program
!      END SUBROUTINE MAKE_O3_FILE
!
!!------------------------------------------------------------------------------
!      SUBROUTINE READ_O3_FILE(  )
!!
!!******************************************************************************
!!  Subroutine READ_O3_FILE reads the GC modeled O3 profiles that correspond
!!  to the TES O3 times and locations. (dkh, 03/01/09) 
!!
!!  NOTES:
!!
!!******************************************************************************
!!
!      ! Reference to F90 modules
!      USE BPCH2_MOD,         ONLY : READ_BPCH2
!      USE DIRECTORY_ADJ_MOD, ONLY : ADJTMP_DIR
!
!
!#     include "CMN_SIZE"          ! Size parameters
!
!      ! Local variables
!      REAL*4                     :: DAT(1,LLPAR,MAXTES)
!      CHARACTER(LEN=255)         :: FILENAME
!
!      !=================================================================
!      ! READ_USA_MASK begins here!
!      !=================================================================
!
!      ! File name
!      FILENAME = TRIM( ADJTMP_DIR )           //
!     &           'nh3.bpch'
!
!      ! Echo info
!      WRITE( 6, 100 ) TRIM( FILENAME )
! 100  FORMAT( '     - READ_O3_FILE: Reading ', a )
!      
!      
!      ! USA mask is stored in the bpch file as #2
!      CALL READ_BPCH2( FILENAME, 'IJ-AVE-$', 1,
!     &                 1d0,            1,        LLPAR, 
!     &                 MAXTES,    DAT,      QUIET=.TRUE. )
!      
!      ! Cast to REAL*8
!      O3_SAVE(:,:) = DAT(1,:,:)
!      
!      print*, ' O3_SAVE sum read = ', SUM(O3_SAVE(:,:))
!
!      ! Return to calling program
!      END SUBROUTINE READ_O3_FILE
!
!!-----------------------------------------------------------------------------
!      FUNCTION GET_DOUBLED_O3( NYMD, NHMS, LON, LAT ) RESULT( O3_DBL )
!!
!!******************************************************************************
!!  Subroutine GET_DOUBLED_O3 reads and returns the nh3 profiles from 
!!  model run with doubled emissions. (dkh, 11/08/09) 
!!
!!  NOTES:
!!
!!******************************************************************************
!!
!      ! Reference to F90 modules
!      USE BPCH2_MOD,         ONLY : READ_BPCH2
!      USE DIRECTORY_MOD,     ONLY : DATA_DIR
!      USE TIME_MOD,          ONLY : EXPAND_DATE
!      USE TIME_MOD,          ONLY : GET_TAU
!
!
!#     include "CMN_SIZE"          ! Size parameters
!
!      ! Arguments   
!      INTEGER                    :: NYMD, NHMS
!      REAL*4                     :: LON,  LAT
!      
!      ! Function arg 
!      REAL*8                     :: O3_DBL(LLPAR)
!
!      ! Local variables
!      REAL*4                     :: DAT(144,91,20)
!      CHARACTER(LEN=255)         :: FILENAME
!      INTEGER                    :: IIJJ(2)
!
!      !=================================================================
!      ! GET_DOUBLED_O3 begins here!
!      !=================================================================
!
!      ! filename
!      FILENAME = 'nh3.YYYYMMDD.hhmm'
!
!      ! Expand filename
!      CALL EXPAND_DATE( FILENAME, NYMD, NHMS )
!
!      ! Full path to file
!      FILENAME = TRIM( DATA_DIR )           //
!     &           'doubled_nh3/'             // 
!     &           TRIM( FILENAME )           //
!     &           TRIM( '00'     )           
!
!      ! Echo info
!      WRITE( 6, 100 ) TRIM( FILENAME )
! 100  FORMAT( '     - GET_DOUBLED_O3: Reading ', a )
!      
!      ! dkh debug
!      print*, ' GET_TAU() = ', GET_TAU()
!      
!      ! Get data
!      CALL READ_BPCH2( FILENAME, 'IJ-AVG-$', 29, 
!     &                 GET_TAU(), 144,      91, 
!     &                 20,        DAT,      QUIET=.FALSE. )
!      
!      IIJJ = GET_IJ_2x25( LON, LAT )
! 
!      print*, ' found doubled in I/J = ', IIJJ
!
!      ! just the column for the present location, and convert ppb to ppm
!      O3_DBL(1:20)     = REAL(DAT(IIJJ(1),IIJJ(2),:),8) / 1000d0 
!      O3_DBL(21:LLPAR) = 0d0 
!     
!      print*, ' O3_DBL = ', O3_DBL
! 
!      ! Return to calling program
!      END FUNCTION GET_DOUBLED_O3
!
!!------------------------------------------------------------------------------
      FUNCTION GET_IJ_2x25( LON, LAT ) RESULT ( IIJJ )

!
!******************************************************************************
!  Subroutine GET_IJ_2x25 returns I and J index from the 2 x 2.5 grid for a 
!  LON, LAT coord. (dkh, 11/08/09) 
! 
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) LON (REAL*8) : Longitude                          [degrees]
!  (2 ) LAT (REAL*8) : Latitude                           [degrees]
!     
!  Function result
!  ============================================================================
!  (1 ) IIJJ(1) (INTEGER) : Long index                    [none]
!  (2 ) IIJJ(2) (INTEGER) : Lati index                    [none]
!     
!  NOTES:
!
!******************************************************************************
!     
      ! Reference to f90 modules
      USE ERROR_MOD,    ONLY : ERROR_STOP

      ! Arguments
      REAL*4    :: LAT, LON
      
      ! Return
      INTEGER :: I, J, IIJJ(2)
      
      ! Local variables 
      REAL*8              :: TLON, TLAT, DLON, DLAT
      REAL*8,  PARAMETER  :: DISIZE = 2.5d0
      REAL*8,  PARAMETER  :: DJSIZE = 2.0d0
      INTEGER, PARAMETER  :: IIMAX  = 144
      INTEGER, PARAMETER  :: JJMAX  = 91
      
      
      !=================================================================
      ! GET_IJ_2x25 begins here!
      !=================================================================

      TLON = 180d0 + LON + DISIZE
      TLAT =  90d0 + LAT + DJSIZE
      
      I = TLON / DISIZE
      J = TLAT / DJSIZE

      
      IF ( TLON / DISIZE - REAL(I)  >= 0.5d0 ) THEN
         I = I + 1
      ENDIF
      
      IF ( TLAT / DJSIZE - REAL(J)  >= 0.5d0 ) THEN
         J = J + 1
      ENDIF

      
      ! Longitude wraps around
      !IF ( I == 73 ) I = 1 
      IF ( I == ( IIMAX + 1 ) ) I = 1
      
      ! Check for impossible values 
      IF ( I > IIMAX .or. J > JJMAX .or. 
     &     I < 1     .or. J < 1          ) THEN
         CALL ERROR_STOP('Error finding grid box', 'GET_IJ_2x25')
      ENDIF
      
      IIJJ(1) = I
      IIJJ(2) = J
      
      ! Return to calling program
      END FUNCTION GET_IJ_2x25

!------------------------------------------------------------------------------

      SUBROUTINE SVD(A,N,U,S,VT)
!
!******************************************************************************
!  Subroutine SVD is a driver for the LAPACK SVD routine DGESVD. (dkh, 05/04/10) 
! 
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) A   (REAL*8) :  N x N matrix to decompose
!  (2 ) N  (INTEGER) :  N is dimension of A
!     
!  Arguments as Output:
!  ============================================================================
!  (1 ) U   (REAL*8) :  Array of left singular vectors
!  (2 ) S   (REAL*8) :  Vector of singular values
!  (3 ) VT  (REAL*8) :  Array of right singular vectors, TRANSPOSED 
!      
!     
!  NOTES:
!
*  Copyright (C) 2009-2010 Intel Corporation. All Rights Reserved.
*  The information and material ("Material") provided below is owned by Intel
*  Corporation or its suppliers or licensors, and title to such Material remains
*  with Intel Corporation or its suppliers or licensors. The Material contains
*  proprietary information of Intel or its suppliers and licensors. The Material
*  is protected by worldwide copyright laws and treaty provisions. No part of
*  the Material may be copied, reproduced, published, uploaded, posted,
*  transmitted, or distributed in any way without Intel's prior express written
*  permission. No license under any patent, copyright or other intellectual
*  property rights in the Material is granted to or conferred upon you, either
*  expressly, by implication, inducement, estoppel or otherwise. Any license
*  under such intellectual property rights must be express and approved by Intel
*  in writing.
*  =============================================================================
*
*  DGESVD Example.
*  ==============
*
*  Program computes the singular value decomposition of a general
*  rectangular matrix A:
*
*    8.79   9.93   9.83   5.45   3.16
*    6.11   6.91   5.04  -0.27   7.98
*   -9.15  -7.93   4.86   4.85   3.01
*    9.57   1.64   8.83   0.74   5.80
*   -3.49   4.02   9.80  10.00   4.27
*    9.84   0.15  -8.99  -6.02  -5.31
*
*  Description.
*  ============
*
*  The routine computes the singular value decomposition (SVD) of a real
*  m-by-n matrix A, optionally computing the left and/or right singular
*  vectors. The SVD is written as
*
*  A = U*SIGMA*VT
*
*  where SIGMA is an m-by-n matrix which is zero except for its min(m,n)
*  diagonal elements, U is an m-by-m orthogonal matrix and VT (V transposed)
*  is an n-by-n orthogonal matrix. The diagonal elements of SIGMA
*  are the singular values of A; they are real and non-negative, and are
*  returned in descending order. The first min(m, n) columns of U and V are
*  the left and right singular vectors of A.
*
*  Note that the routine returns VT, not V.
*
*  Example Program Results.
*  ========================
*
* DGESVD Example Program Results
*
* Singular values
*  27.47  22.64   8.56   5.99   2.01
*
* Left singular vectors (stored columnwise)
*  -0.59   0.26   0.36   0.31   0.23
*  -0.40   0.24  -0.22  -0.75  -0.36
*  -0.03  -0.60  -0.45   0.23  -0.31
*  -0.43   0.24  -0.69   0.33   0.16
*  -0.47  -0.35   0.39   0.16  -0.52
*   0.29   0.58  -0.02   0.38  -0.65
*
* Right singular vectors (stored rowwise)
*  -0.25  -0.40  -0.69  -0.37  -0.41
*   0.81   0.36  -0.25  -0.37  -0.10
*  -0.26   0.70  -0.22   0.39  -0.49
*   0.40  -0.45   0.25   0.43  -0.62
*  -0.22   0.14   0.59  -0.63  -0.44
*  =============================================================================
!******************************************************************************
!
      ! Arguements 
      INTEGER,INTENT(IN)     :: N
      REAL*8, INTENT(IN)     :: A(N,N)
      REAL*8, INTENT(OUT)    :: U(N,N)
      REAL*8, INTENT(OUT)    :: S(N)
      REAL*8, INTENT(OUT)    :: VT(N,N)

      ! Local variables 
      INTEGER, PARAMETER     :: LWMAX = MAXLEV * 35 
      INTEGER                :: INFO, LWORK
      DOUBLE PRECISION       :: WORK( LWMAX )

*     .. External Subroutines ..
      EXTERNAL               :: DGESVD

*     .. Intrinsic Functions ..
      INTRINSIC              :: INT, MIN

      !=================================================================
      ! SVD begins here!
      !=================================================================

*     .. Executable Statements ..
      !WRITE(*,*)'DGESVD Example Program Results'
*
*     Query the optimal workspace.
*
      print*, ' here 1 '
      LWORK = -1
      CALL DGESVD( 'All', 'All', N, N, A, N, S, U, N, VT, N,
     $             WORK, LWORK, INFO )
      LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
*
*     Compute SVD.
*
      print*, ' here 2 '
      CALL DGESVD( 'All', 'All', N, N, A, N, S, U, N, VT, N,
     $             WORK, LWORK, INFO )
*
*     Check for convergence.
*
      print*, ' here 3 '
      IF( INFO.GT.0 ) THEN
         WRITE(*,*)'The algorithm computing SVD failed to converge.'
         STOP
      END IF

!  Uncomment the following to print out singlular values, vectors (dkh, 05/04/10) 
!
!     Print singular values.
!
      CALL PRINT_MATRIX( 'Singular values', 1, N, S, 1 )
!
!     Print left singular vectors.
!
      CALL PRINT_MATRIX( 'Left singular vectors (stored columnwise)',
     $                   N, N, U, N   )
!
!     Print right singular vectors.
!
      CALL PRINT_MATRIX( 'Right singular vectors (stored rowwise)',
     $                   N, N, VT, N    )

      ! Return to calling program
      END SUBROUTINE SVD
!------------------------------------------------------------------------------
      SUBROUTINE DGESVD_EXAMPLE

*     .. Parameters ..
      INTEGER          M, N
      PARAMETER        ( M = 6, N = 5 )
      INTEGER          LDA, LDU, LDVT
      PARAMETER        ( LDA = M, LDU = M, LDVT = N )
      INTEGER          LWMAX
      PARAMETER        ( LWMAX = 1000 )
*
*     .. Local Scalars ..
      INTEGER          INFO, LWORK
*
*     .. Local Arrays ..
      DOUBLE PRECISION A( LDA, N ), U( LDU, M ), VT( LDVT, N ), S( N ),
     $                 WORK( LWMAX )
      DATA             A/
     $  8.79, 6.11,-9.15, 9.57,-3.49, 9.84,
     $  9.93, 6.91,-7.93, 1.64, 4.02, 0.15,
     $  9.83, 5.04, 4.86, 8.83, 9.80,-8.99,
     $  5.45,-0.27, 4.85, 0.74,10.00,-6.02,
     $  3.16, 7.98, 3.01, 5.80, 4.27,-5.31
     $                  /
*
*     .. External Subroutines ..
      EXTERNAL         DGESVD
      !EXTERNAL         PRINT_MATRIX
*
*     .. Intrinsic Functions ..
      INTRINSIC        INT, MIN
*
*     .. Executable Statements ..
      WRITE(*,*)'DGESVD Example Program Results'
*
*     Query the optimal workspace.
*
      LWORK = -1
      CALL DGESVD( 'All', 'All', M, N, A, LDA, S, U, LDU, VT, LDVT,
     $             WORK, LWORK, INFO )
      LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
*
*     Compute SVD.
*
      CALL DGESVD( 'All', 'All', M, N, A, LDA, S, U, LDU, VT, LDVT,
     $             WORK, LWORK, INFO )
*
*     Check for convergence.
*
      IF( INFO.GT.0 ) THEN
         WRITE(*,*)'The algorithm computing SVD failed to converge.'
         STOP
      END IF
*
*     Print singular values.
*
      CALL PRINT_MATRIX( 'Singular values', 1, N, S, 1 )
*
*     Print left singular vectors.
*
      CALL PRINT_MATRIX( 'Left singular vectors (stored columnwise)',
     $                   M, N, U, LDU )
*
*     Print right singular vectors.
*
      CALL PRINT_MATRIX( 'Right singular vectors (stored rowwise)',
     $                   N, N, VT, LDVT )

*
*     End of DGESVD Example.
      END SUBROUTINE DGESVD_EXAMPLE
!------------------------------------------------------------------------------
*
*     Auxiliary routine: printing a matrix.
*
      SUBROUTINE PRINT_MATRIX( DESC, M, N, A, LDA )
      CHARACTER*(*)    DESC
      INTEGER          M, N, LDA
      DOUBLE PRECISION A( LDA, * )
*
      INTEGER          I, J
*
      WRITE(*,*)
      WRITE(*,*) DESC
      DO I = 1, M
         WRITE(*,9998) ( A( I, J ), J = 1, N )
      END DO
*
! Change format of output (dkh, 05/04/10) 
! 9998 FORMAT( 11(:,1X,F6.2) )
 9998 FORMAT( 11(:,1X,E14.8) )
      RETURN

      END SUBROUTINE PRINT_MATRIX 
!------------------------------------------------------------------------------


      END MODULE FLASK_CO2_MOD
