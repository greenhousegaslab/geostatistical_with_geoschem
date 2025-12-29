!$Id: population_mod.f,v 1.1 2012/03/01 22:00:27 daven Exp $
       MODULE POPULATION_MOD
!
!******************************************************************************
!  Module POPULATION_MOD contains code for incorporating  population weighting 
!  into cost functions / exposure metrics. Population data taken from:
!
!    Center for International Earth Science Information Network (CIESIN), 
!    Columbia University; and Centro Internacional de Agricultura Tropical 
!    (CIAT). 2005. Gridded Population of the World, Version 3 (GPWv3): 
!    Population Count Grid. Palisades, NY: Socioeconomic Data and Applications 
!    Center (SEDAC), Columbia University. 
!    Available at http://sedac.ciesin.columbia.edu/gpw. 2/11/2012.
!
!  Steven Vogel, jk, dkh, 02/04/2012, adj32_024
!     
!  Module Variables:
!  ============================================================================
!  (1 ) POP_REDUCED   (REAL*8) : Array of census population
!
!  Module Routines:
!  ===========================================================================
!  (1 ) POP_WEIGHT_COST        : Computes population weighted cost function
!  (2 ) READ_IN_POPULATION     : Reads in gridded population data file
!  (3 ) INIT_POPULATOIN_MOD    : Allocates & zeroes module arrays
!  (4 ) CLEANUP_POPULATION_MOD : Deallocates module arrays
!
!  NOTES:
!
!*****************************************************************************
      IMPLICIT NONE

      PUBLIC

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================
      REAL*8,  ALLOCATABLE :: POP_REDUCED(:,:)

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE POP_WEIGHT_COST
!
!******************************************************************************
! This subroutine based on CALC_ADJ_FORCE_FOR_SENSE in geos_chem_adj_mod.f
! Calculates population weighted cost function when called in 
! geos_chem_adj_mod.f
!
!******************************************************************************
!
      ! References to F90 modules 
      USE ADJ_ARRAYS_MOD,       ONLY : N_CALC, COST_FUNC
      USE ADJ_ARRAYS_MOD,       ONLY : STT_ADJ
      USE ADJ_ARRAYS_MOD,       ONLY : GET_CF_REGION
      USE ADJ_ARRAYS_MOD,       ONLY : NSPAN
      USE ADJ_ARRAYS_MOD,       ONLY : OBS_THIS_TRACER
      USE ADJ_ARRAYS_MOD,       ONLY : NOBS
      USE ADJ_ARRAYS_MOD,       ONLY : TRACER_IND
      USE CHECKPT_MOD,          ONLY : CHK_STT
      USE DAO_MOD,              ONLY : AIRVOL, AD
      USE LOGICAL_MOD,          ONLY : LPRT
      USE TRACER_MOD,           ONLY : N_TRACERS
      USE GRID_MOD,             ONLY : GET_XMID, GET_YMID

      ! Header files
#     include "CMN_SIZE"             ! Size parameters

      ! Local variables
      REAL*8              :: ADJ_FORCE(IIPAR,JJPAR,LLPAR,N_TRACERS)
      INTEGER             :: I, J, N, NN
      LOGICAL, SAVE       :: FIRST = .TRUE.

      REAL*8,  DIMENSION(IIPAR,JJPAR,NOBS) :: COST_NUMM
      REAL*8,  DIMENSION(IIPAR,JJPAR)      :: DENOMM_POP
      REAL*8,  DIMENSION(IIPAR,JJPAR)      :: DENOMM_VOL
      REAL*4,  DIMENSION(IIPAR,JJPAR)      :: N_INCLUDE
      REAL*8              :: NEW_COST_SCALAR
      REAL*8              :: POP_TOT
      REAL*8              :: VOL_TOT
      REAL*8              :: N_TOT
      REAL*8              :: FACTORR

      !=================================================================
      ! POP_WEIGHT_COST begins here!
      !=================================================================

      ! Get population data     
      IF ( FIRST ) THEN

         CALL INIT_POPULATION_MOD 
         CALL READ_IN_POPULATION

         ! replace NOBS2STT with TRACER_IND

         FIRST = .FALSE.

      ENDIF

      IF ( LPRT ) THEN  
         print*, 'SEV DEBUG = ', maxval(POP_REDUCED)
         print*, 'SEV DEBUG = ', minval(POP_REDUCED)
      ENDIF 

      ! Initialze cost fnc variables
      NEW_COST_SCALAR = 0d0
      COST_NUMM       = 0d0
      DENOMM_POP      = 0d0
      DENOMM_VOL      = 0d0
      N_INCLUDE       = 0e0

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I,    J,   N,   NN )
      DO N = 1, NOBS
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         NN = TRACER_IND(N)

         ! Determine the contribution to the cost function in each grid cell
         ! from each species
         COST_NUMM(I,J,N)  = GET_CF_REGION(I,J,1)
     &                     * CHK_STT(I,J,1,NN)
     &                     * POP_REDUCED(I,J)

         ! Set denominator population and volume 
         IF ( N == 1 .and. GET_CF_REGION(I,J,1) > 0d0 
     &               .and. POP_REDUCED(I,J)     > 0d0  ) THEN

            DENOMM_POP(I,J)  =
     &                        POP_REDUCED(I,J) * GET_CF_REGION(I,J,1)

            DENOMM_VOL(I,J)  =
     &                        AIRVOL(I,J,1) * GET_CF_REGION(I,J,1)

            N_INCLUDE(I,J) = 1e0 

! For debugging:
!       WRITE(55,100) I, J, GET_XMID(I), GET_YMID(J), POP_REDUCED(I,J), 
!     &       AIRVOL(I,J,1), CHK_STT(I,J,1,27), 
!     &      CHK_STT(I,J,1,31), CHK_STT(I,J,1,32), CHK_STT(I,J,1,34:37)
! 100   FORMAT(1X,I6,1X,I6,1X,F16.8,1X,F16.8,1X,E16.8,1X,E16.8,1X,F16.8,1X,
!     &   F16.8,1X, F16.8,1X, F16.8,1X,F16.8,1X,F16.8,1X,F16.8)

         ENDIF

      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      POP_TOT = SUM(DENOMM_POP)
      VOL_TOT = SUM(DENOMM_VOL)
      N_TOT   = SUM(N_INCLUDE)

         print*, 'SEV DEBUG TOTAL VOLUME = ', VOL_TOT
         print*, 'SEV DEBUG TOTAL POP = ', POP_TOT
         print*, 'SEV DEBUG TOTAL N   = ', N_TOT
#if   defined ( GRID2x25 )
         print*, 'SEV DEBUG DENOM POP= ', DENOMM_POP(117,64)
#endif 

      FACTORR  = 1d9 / ( POP_TOT * VOL_TOT * NSPAN ) * N_TOT

      NEW_COST_SCALAR = SUM(COST_NUMM(:,:,:)) * FACTORR

      ! Update cost function
      COST_FUNC = COST_FUNC + NEW_COST_SCALAR

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I,    J,    N,   NN )
      DO N = 1, NOBS
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         NN = TRACER_IND(N)

         ! Force the adjoint variables x with dJ/dx=1
         ADJ_FORCE(I,J,1,NN)  = GET_CF_REGION(I,J,1)
     &                        * POP_REDUCED(I,J)
     &                        * FACTORR

         STT_ADJ(I,J,1,NN)   = STT_ADJ(I,J,1,NN) + ADJ_FORCE(I,J,1,NN)

      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Return to calling program
      END SUBROUTINE POP_WEIGHT_COST

!------------------------------------------------------------------------------
 
       SUBROUTINE READ_IN_POPULATION
!        
!******************************************************************************
!  Subroutine READ_IN_POPULATION reads in gridded population data. 
!  by Steven Vogel, based on code from Jamin Koo (dkh, 02/13/12, adj32_024) 
!
!  NOTES:
!
!******************************************************************************
!

      ! References to F90 modules 
      USE BPCH2_MOD,     ONLY : GET_RES_EXT
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE LOGICAL_MOD,   ONLY : LPRT

#     include "CMN_SIZE"             ! Size parameters

      ! Local variables
      CHARACTER(LEN=255)     :: FNAME 
      INTEGER                :: IOS, IOS2

      !=================================================================
      ! READ_IN_POPULATION begins here!
      !=================================================================

      ! Generate population data filename 
      FNAME = TRIM( DATA_DIR ) // 'population_201202/' // 
     &             'world_population.' // GET_RES_EXT()

      ! Read the population from ascii file.
      WRITE( 6, '(a)' ) '  Reading in population from ', FNAME

      OPEN( UNIT=11, FILE=FNAME, STATUS='OLD', IOSTAT=IOS)

      IF ( IOS /= 0 ) THEN
         PRINT *, 'ERROR opening weight'
      ELSE
         READ( UNIT=11, FMT=*, IOSTAT=IOS2 ) POP_REDUCED
         IF ( IOS2 < 0 ) THEN 
            WRITE( 6, '(a)' ) '  Unexpected End of File encountered  '
         ELSE IF ( IOS > 0 ) THEN
            WRITE( 6, '(a)' ) '  Error occurred reading pop data!  '
         ENDIF
      ENDIF

      CLOSE( UNIT=10 )

      IF ( LPRT ) THEN 
         PRINT *, 'sum of population',         sum(POP_REDUCED)
         PRINT *, 'Population Grid Test Max',  maxval(POP_REDUCED)
         PRINT *, 'Population Grid Test Min',  minval(POP_REDUCED)
         PRINT *, 'Population Grid Test Size', size(POP_REDUCED)
      ENDIF 

      !CALL MAKE_POP_FILE()

      ! Return to calling program
      END SUBROUTINE READ_IN_POPULATION      

!------------------------------------------------------------------------------
      SUBROUTINE MAKE_POP_FILE( )
!
!******************************************************************************
!  Subroutine MAKE_POP_FILE creates binary world population file. 
!  (dkh, 9/01/04)
!
!******************************************************************************
      ! References to F90 modules
      USE BPCH2_MOD
      USE ERROR_MOD,           ONLY : DEBUG_MSG
      USE FILE_MOD,            ONLY : IU_RST,      IOERROR
      USE GRID_MOD,            ONLY : GET_XOFFSET, GET_YOFFSET
      USE LOGICAL_MOD,         ONLY : LPRT
      USE TIME_MOD,            ONLY : GET_TAU

#     include "CMN_SIZE"   ! Size parameters
      
      ! Local Variables     
      INTEGER              :: I,    I0, IOS, J,  J0, L, N
      INTEGER              :: YYYY, MM, DD,  HH, SS
      CHARACTER(LEN=255)   :: FILENAME

      ! For binary punch file, version 2.0
      REAL*4               :: LONRES, LATRES
      INTEGER, PARAMETER   :: HALFPOLAR = 1
      INTEGER, PARAMETER   :: CENTER180 = 1

      CHARACTER(LEN=20)    :: MODELNAME
      CHARACTER(LEN=40)    :: CATEGORY
      CHARACTER(LEN=40)    :: UNIT    
      CHARACTER(LEN=40)    :: RESERVED = ''
      CHARACTER(LEN=80)    :: TITLE

      
      !=================================================================
      ! MAKE_POP_FILE begins here!
      !=================================================================
      
      ! Define variables for BINARY PUNCH FILE OUTPUT
      TITLE    = 'GEOS-CHEM OBS File: ' //
     &           'Observation Concentrations (kg/box)'
      UNIT     = 'people'
      CATEGORY = 'IJ-POP-$'
      LONRES   = DISIZE
      LATRES   = DJSIZE
      
      ! Call GET_MODELNAME to return the proper model name for
      ! the given met data being used (bmy, 6/22/00)
      MODELNAME = GET_MODELNAME()
      
      ! Get the nested-grid offsets
      I0 = GET_XOFFSET( GLOBAL=.TRUE. )
      J0 = GET_YOFFSET( GLOBAL=.TRUE. )
      
      !=================================================================
      ! Open the observation file for output -- binary punch format
      !=================================================================
      
      ! Copy the output observation file name into a local variable
      FILENAME = TRIM( 'pop.bpch' )
      
      WRITE( 6, 100 ) TRIM( FILENAME )
 100  FORMAT( '     - MAKE_POP_FILE: Writing ', a )
 
      ! Open checkpoint file for output
      CALL OPEN_BPCH2_FOR_WRITE( IU_RST, FILENAME, TITLE )
      
      !=================================================================
      ! Write pop
      !=================================================================
      
      CALL BPCH2( IU_RST,    MODELNAME, LONRES,    LATRES,
     &            HALFPOLAR, CENTER180, CATEGORY,  1,
     &            UNIT,      GET_TAU(), GET_TAU(), RESERVED,
     &            IIPAR,     JJPAR,     1,         I0+1,
     &            J0+1,      1,         REAL(POP_REDUCED,4) )
      
      
      ! Close file
      CLOSE( IU_RST )
      
      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### MAKE_POP_FILE: wrote file' )
      
      ! Return to calling program
      END SUBROUTINE MAKE_POP_FILE
!------------------------------------------------------------------------------

      SUBROUTINE INIT_POPULATION_MOD
!
!******************************************************************************
!  Subroutine INIT_POPULATION_MOD initializes and zeros all allocatable arrays
!  declared in "population_mod.f"
!     
!  NOTES:
!     
!******************************************************************************
!
      ! Reference to f90 modules
      USE ADJ_ARRAYS_MOD, ONLY : NOBS
      USE ERROR_MOD,      ONLY : ALLOC_ERR

#     include "CMN_SIZE"       ! Size parameters

      ! local variables 
      INTEGER                :: AS
      
      !=================================================================
      ! INIT_POPULATION_MOD
      !=================================================================

      ALLOCATE( POP_REDUCED(IIPAR,JJPAR) , STAT = AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'POP_REDUCED' ) 
      POP_REDUCED = 0d0
               
      ! Return to calling program
      END SUBROUTINE INIT_POPULATION_MOD
                  
!-----------------------------------------------------------------------------
      SUBROUTINE CLEANUP_POPULATION_MOD
!              
!******************************************************************************
!  Subroutine CLEANUP_POPULATION_MOD deallocates all previously allocated arrays 
!              
!  NOTES:         
!
!******************************************************************************
!              
      !=================================================================
      ! CLEANUP_POPULATION_MOD begins here!
      !=================================================================
      IF ( ALLOCATED( POP_REDUCED ) ) DEALLOCATE( POP_REDUCED )
               
      ! Return to calling program
      END SUBROUTINE CLEANUP_POPULATION_MOD
                  
!------------------------------------------------------------------------------
      
      END MODULE POPULATION_MOD
