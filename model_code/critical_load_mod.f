!
      MODULE CRITICAL_LOAD_MOD

!module critical_load contains variable  and routines to read acidity/N deposition critical loads over the US (SD) and generate map of exceedence


      IMPLICIT NONE

      ! Make everything private
      PRIVATE


      !except
      PUBLIC :: GET_CL_EXCEEDENCE
      PUBLIC :: CL_FILENAME
      PUBLIC :: GC_FILENAME

      CHARACTER*255 :: CL_FILENAME
      CHARACTER*255 :: GC_FILENAME

      !very large value used to generate critical load files
      !to make sure no exceedence is found in regions where CL is not defined
      INTEGER, PARAMETER :: MISSING_VALUE = 9D6

      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE READ_CRITICAL_LOAD(CRIT_L)
! Subroutine read_critical_load is used to read annual critical/load acidification files.
! Files are in netcdf format (spatial resolution is reduced nested domain over the US)

      USE TRANSFER_MOD,    ONLY : TRANSFER_2D
      USE NETCDF_UTIL_MOD

#     include "CMN_SIZE"

      REAL*8, INTENT(OUT) :: CRIT_L(IIPAR,JJPAR)
      REAL*4              :: ARRAY(IIPAR,JJPAR)
      INTEGER             :: I,J
      INTEGER             :: fileID, varID
      CHARACTER*255       :: FILENAME

      WRITE(6,*) '=========================='
      WRITE(6,*) '=== Read Critical load ==='
      WRITE(6,*) '=========================='


      call ncdf_open_for_read( fileID, TRIM(CL_FILENAME) )
      varID = ncdf_get_varid( fileID, 'ecoreg' )

      call ncdf_get_var( fileID, varID, array,
     &     start=(/     1,     1  /),
     &     count=(/ iipar, jjpar  /)  )

      ! real*4->real*8
      CALL TRANSFER_2D( ARRAY, CRIT_L )

      WRITE(*,*) 'Min Critical Load: ',MINVAL( CRIT_L )
      WRITE(*,*) 'Max Critical Load: ',MAXVAL( CRIT_L,
     &     MASK = CRIT_L < MISSING_VALUE )

      END SUBROUTINE READ_CRITICAL_LOAD

!------------------------------------------------------------------------------

      SUBROUTINE READ_GC_LOAD(GC_L)
! Subroutine read_annual deposition. Read GC calculated deposition from 3yr nested simulation (lzh)

      USE TRANSFER_MOD,    ONLY : TRANSFER_2D
      USE NETCDF_UTIL_MOD

      USE TRACERID_MOD,    ONLY : IDTNH3, IDTNH4, IDTNIT, IDTNITs
      USE TRACERID_MOD,    ONLY : IDTHNO3, IDTR4N2, IDTPMN, IDTPPN
      USE TRACERID_MOD,    ONLY : IDTPAN, IDTNOX, IDTN2O5
      USE TRACERID_MOD,    ONLY : IDTSO2, IDTSO4, IDTSO4s
      USE TRACER_MOD,      ONLY : TRACER_MW_KG, TRACER_NAME
      USE LOGICAL_ADJ_MOD, ONLY : LADJ_CL_ACID, LADJ_CL_NDEP
      USE GRID_MOD,        ONLY : GET_AREA_M2
      USE LOGICAL_MOD,     ONLY : LPRT

#     include "CMN_SIZE"

      REAL*8, INTENT(OUT) :: GC_L(IIPAR,JJPAR)
      REAL*4              :: ARRAY(IIPAR,JJPAR)
      REAL*8              :: ARRAY8(IIPAR,JJPAR)
      REAL*8              :: SWET(IIPAR,JJPAR)
      INTEGER             :: N, NMAX
      INTEGER             :: fileID, varID
      INTEGER             :: I,J
      REAL*8              :: TRACER(18)
      REAL*8              :: DRY(IIPAR,JJPAR)
      REAL*8              :: TEMP(IIPAR,JJPAR)


      WRITE(6,*) '===================='
      WRITE(6,*) '=== Read GC load ==='
      WRITE(6,*) '===================='

      TRACER( 1 )  = IDTHNO3
      TRACER( 2 )  = IDTR4N2
      TRACER( 3 )  = IDTPMN
      TRACER( 4 )  = IDTPPN
      TRACER( 5 )  = IDTPAN
      TRACER( 6 )  = IDTNOX
      TRACER( 7 )  = IDTN2O5
      TRACER( 8 )  = IDTN2O5 !repeat N2O5 as it contains two N
      TRACER( 9 )  = IDTNH3
      TRACER(10 )  = IDTNH4
      TRACER(11 )  = IDTNIT
      TRACER(12 )  = IDTNITs

      !repeat S tracers since they contain two equivalents
      TRACER(13) = IDTSO2
      TRACER(14) = IDTSO4
      TRACER(15) = IDTSO4s
      TRACER(16) = IDTSO2
      TRACER(17) = IDTSO4
      TRACER(18) = IDTSO4s

      GC_L(:,:) = 0d0
      SWET(:,:) = 0d0
      DRY(:,:)  = 0d0


      IF ( LADJ_CL_ACID ) THEN
         NMAX = 18
      ELSEIF (LADJ_CL_NDEP) THEN
         NMAX = 12
      ENDIF

!     open netcdf file
      call ncdf_open_for_read( fileID, TRIM(GC_FILENAME) )

      DO N = 1,NMAX

         ! read dry deposition
         IF (LPRT) THEN
            WRITE(6,100) TRIM(TRACER_NAME(TRACER(N)))
         ENDIF

         IF ( TRACER(N) .EQ. IDTNOx) THEN

            varID = ncdf_get_varid( fileID,
     &           'DRY_NO2' )

         ELSE

            varID = ncdf_get_varid( fileID,
     &           'DRY_' // TRIM(TRACER_NAME(TRACER(N))) )

         ENDIF

         call ncdf_get_var( fileID, varID, array,
     &        start=(/     1,     1  /),
     &        count=(/ iipar, jjpar  /)  )

         ! real*4->real*8
         CALL TRANSFER_2D( ARRAY, ARRAY8 )

         IF ( LADJ_CL_NDEP ) THEN

            ! convert from molec/cm2/s to kgN/ha/yr
            GC_L(:,:) = GC_L(:,:) +
     &           14D-3 / 6.022D23 *
     &           1D8 *
     &           86400D0 * 365D0 *
     &           ARRAY8(:,:)

            !for diagnostics (kgN/cm2/yr)
            DRY(:,:) = DRY(:,:) +
     &           14D-3 / 6.022D23 *
     &           86400D0 * 365D0 *
     &           ARRAY8(:,:)

            DO J = 1, JJPAR

               TEMP(:,J) =  14D-3 / 6.022D23 *
     &           86400D0 * 365D0 *
     &           ARRAY8(:,J) * GET_AREA_M2(J) *
     &           1D4

            ENDDO

         ELSEIF ( LADJ_CL_ACID ) THEN

            ! convert from molec/cm2/s to equiv/ha/yr
            GC_L(:,:) = GC_L(:,:) +
     &           1D0 / 6.022D23 *
     &           1D8 *
     &           86400D0 * 365D0 *
     &           ARRAY8(:,:)

            ! for diagnostics equiv/cm2/yr
            DRY(:,:) = DRY(:,:) +
     &           1D0 / 6.022D23 *
     &           86400D0 * 365D0 *
     &           ARRAY8(:,:)

         ENDIF

         ! read wet deposition
         IF ( LPRT ) THEN
            WRITE(6,101) TRIM( TRACER_NAME(TRACER(N)) )
         ENDIF

         ! no wet deposition for NOx so no need to read NOx
         IF ( IDTNOx .NE. TRACER(N) ) THEN

            varID = ncdf_get_varid( fileID,
     &           'WET_' // TRIM( TRACER_NAME( TRACER(N))) )

            call ncdf_get_var( fileID, varID, array,
     &           start=(/     1,     1  /),
     &           count=(/ iipar, jjpar  /)  )

            ! real*4->real*8
            CALL TRANSFER_2D( ARRAY, ARRAY8 )

            IF ( LADJ_CL_NDEP ) THEN

               ! convert from kg/s to kgN/yr
               SWET(:,:) = SWET(:,:) +
     &              14D-3 / TRACER_MW_KG( TRACER(N) ) *
     &              86400D0 * 365D0 *
     &              ARRAY8(:,:)

            ELSEIF ( LADJ_CL_ACID ) THEN

               ! convert from kg/s to equiv/yr
               SWET(:,:) = SWET(:,:) +
     &              1D0 / TRACER_MW_KG( TRACER(N) ) *
     &              86400D0 * 365D0 *
     &              ARRAY8(:,:)

            ENDIF

         ENDIF

      ENDDO

      !for diagnostics
      DO J = 1,JJPAR

         DRY(:,J) = DRY(:,J) *
     &        GET_AREA_M2(J) * 1D4

      ENDDO


      WRITE(*,*) 'GC deposition (wet): ',SUM( SWET )
      WRITE(*,*) 'GC deposition (dry): ',SUM( DRY )


      !convert SWET from equiv/yr to equiv/yr/ha
      !add to critical load

      DO J = 1, JJPAR

         GC_L(:,J) = GC_L(:,J) +
     &        SWET(:,J) / GET_AREA_M2(J) * 1D4

      ENDDO

      WRITE(*,*) 'Min GC Load: ', MINVAL( GC_L(:,:) )
      WRITE(*,*) 'Max GC Load: ', MAXVAL( GC_L(:,:) )


 100  FORMAT('Dry deposition :',a)
 101  FORMAT('Wet deposition :',a)

      END SUBROUTINE READ_GC_LOAD

!------------------------------------------------------------------------------
!
      SUBROUTINE GET_CL_EXCEEDENCE( EXCEEDENCE )
! generate maks for criticial load exceedence
!
      USE BPCH2_MOD
      USE LOGICAL_MOD, ONLY : LPRT
      USE FILE_MOD,    ONLY : IU_DEBUG

#     include "CMN_SIZE"

      REAL*8, INTENT(OUT) :: EXCEEDENCE(IIPAR,JJPAR)
      REAL*8              :: CRIT_L(IIPAR,JJPAR)
      REAL*8              :: GC_L(IIPAR,JJPAR)
      INTEGER             :: I,J,NE,DEFINED

      EXCEEDENCE(:,:) = 0d0
      NE              = 0
      DEFINED         = 0

      ! read GC load from deposition
      CALL READ_GC_LOAD(GC_L)

      ! read critical load (deposition)
      CALL READ_CRITICAL_LOAD(CRIT_L)

      DO I = 1, IIPAR
      DO J = 1, JJPAR

         IF ( GC_L(I,J) .GT. CRIT_L(I,J) ) THEN
            EXCEEDENCE(I,J) = 1D0
            NE = NE + 1
         ENDIF

         ! count the number of grid cells where a critical load is defined
         IF ( CRIT_L(I,J) .LT. MISSING_VALUE ) THEN
            DEFINED = DEFINED + 1
         ENDIF

      ENDDO
      ENDDO

      WRITE(6,*) 'Number of exceedences: ', NE
      WRITE(6,*) 'Fraction of grid cells with exceedence',
     &     REAL(NE)/REAL(DEFINED)


      ! write out exceedence to bpch file
      OPEN( IU_DEBUG, FILE='critical_load.log', STATUS='UNKNOWN' )

      DO I = 1, IIPAR
         DO J = 1, JJPAR
            WRITE( IU_DEBUG, '(F10.3,X)', advance='no' ) EXCEEDENCE(I,J)
         ENDDO
         WRITE( IU_DEBUG,* )
      END DO

      CLOSE( IU_DEBUG )


      END SUBROUTINE GET_CL_EXCEEDENCE

!------------------------------------------------------------------------------

      END MODULE CRITICAL_LOAD_MOD
