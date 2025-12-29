
      MODULE PARANOX_ADJ_MOD

      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC interpolate_lut2_adj
!
! !MODULE VARIABLES
!
      ! fracnox                 = look up table for fraction of NOx remaining
      !                           for ship emissions (gvinken, 6/6/10)
      ! intope                  = look up table for integrated Ozone Production
      !                           Efficiency for ship emiss (gvinken, 6/6/10)
      REAL*4 :: fracnox(4, 4, 4, 12, 12, 4, 5)
      REAL*4 :: intope(4, 4, 4, 12, 12, 4, 5)

      CONTAINS

!  Differentiation of interpolate_lut2 in reverse (adjoint) mode (with options i4 dr8 r8):
!   gradient     of useful results: int_ope fraction_nox
!   with respect to varying inputs: no2 no int_ope fraction_nox
!                o3
!   RW status of diff variables: no2:out no:out int_ope:in-zero
!                fraction_nox:in-zero o3:out
!
!
      SUBROUTINE INTERPOLATE_LUT2_ADJ(i, j, o3, o3b, no, nob, no2, no2b,
     &                                dens,jo1d, jno2, fraction_nox,
     &                                fraction_noxb, int_ope, int_opeb)

!
! !USES:
!
      USE ERROR_MOD,          ONLY : ERROR_STOP
      USE DAO_MOD,            ONLY : TS, SUNCOS, SUNCOS_5hr
      USE PARANOX_MOD,        ONLY : FRACNOX, INTOPE

#     include "CMN_SIZE"  ! Size parameters

!
! !INPUT PARAMETERS:
!
      INTEGER,    INTENT(IN)   :: I, J
      REAL*8,     INTENT(IN)   :: o3, no, no2, dens, jno2, jo1d
      REAL*4,     INTENT(IN)   :: fraction_nox, int_ope
      REAL*4,     INTENT(IN)   :: fraction_noxb, int_opeb

!
! !OUTPUT PARAMETERS
!
      REAL*8,     INTENT(OUT)  :: o3b, nob, no2b


!
! !REVISION HISTORY:
!     Aug 2013 - Yanko Davila - Initial version based on forward subroutine.
!                               Automatic code generated with TAPENADE 3.7
!                               adBufer.f and adStack.c from ISOROPIAII
!
! !LOCAL VARIABLES:
!
!     !=======================================================================
!     ! temp   : model temperature
!     ! jno2   : J(NO2) value
!     ! cao3   : concentration O3 in ambient air
!     ! alfa0  : solar zenith angle 5 hours ago
!     ! alfa5  : solar zenith angle at this time
!     ! jo1d   : ratio J(O1D)/J(NO2)
!     ! canox  : concentration NOx in ambient air
!     !
!     ! o3     : incoming o3 concentration
!     ! no     : incoming no
!     ! no2    : incoming no2
!     ! dens   : incoming air density
!     !=======================================================================

      INTEGER                    :: IJLOOP
      INTEGER, PARAMETER         :: ntemp  = 4
      INTEGER, PARAMETER         :: njno2  = 4
      INTEGER, PARAMETER         :: ncao3  = 4
      INTEGER, PARAMETER         :: nalfa0 = 12
      INTEGER, PARAMETER         :: nalfa5 = 12
      INTEGER, PARAMETER         :: njo1d  = 4
      INTEGER, PARAMETER         :: ncanox = 5


      REAL*4,  DIMENSION(ntemp)  :: templev
      REAL*4,  DIMENSION(njno2)  :: jno2lev
      REAL*4,  DIMENSION(ncao3)  :: cao3lev
      REAL*4,  DIMENSION(nalfa0) :: alfa0lev
      REAL*4,  DIMENSION(nalfa5) :: alfa5lev
      REAL*4,  DIMENSION(njo1d)  :: jo1dlev
      REAL*4,  DIMENSION(ncanox) :: canoxlev

      ! Temporary variable storage
      REAL*4                     :: temp_tmp,  jno2_tmp,  cao3_tmp
      REAL*4                     :: alfa0_tmp, alfa5_tmp, jo1d_tmp
      REAL*4                     :: canox_tmp

      ! ADJOINT of Temporary variable storage
      REAL*4                    :: temp_tmpb,  jno2_tmpb,  cao3_tmpb
      REAL*4                    :: alfa0_tmpb, alfa5_tmpb, jo1d_tmpb
      REAL*4                    :: canox_tmpb

      ! Interpolation parameters
      REAL*4,  DIMENSION(2)      :: xtemp,     xjno2,     xcao3, xalfa0
      REAL*4,  DIMENSION(2)      :: xalfa5,    xjo1d,     xcanox

      ! ADJOINT of Interpolation parameters
      REAL*4, DIMENSION(2)       :: xtempb,  xjno2b, xcao3b, xalfa0b
      REAL*4, DIMENSION(2)       :: xalfa5b, xjo1db, xcanoxb

      ! For loops
      INTEGER                    :: itemp,     ijno2,     icao3, ialfa0
      INTEGER                    :: ialfa5,    ijo1d,     icanox
      INTEGER                    :: i0,i1,i2,i3,i4,i5,i6,i7

      ! array contain temp, jno2, cao3, alfa_0, alfa_5, jo1d, canox
      REAL*4, DIMENSION(7)       :: var_array

      ! ADJOINT of array contain temp, jno2, cao3, alfa_0, alfa_5, jo1d, canox
      REAL*4, DIMENSION(7)       :: var_arrayb

      CHARACTER(len=255)         :: MSG

      ! TAPENADE generated variables
      INTEGER :: branch
      REAL*4 :: temp3
      REAL*4 :: temp2
      REAL*4 :: temp1
      REAL*4 :: temp0
      REAL*4 :: tempb2
      REAL*4 :: tempb1
      REAL*4 :: tempb0
      REAL*4 :: temp2b1
      REAL*4 :: temp2b0
      REAL*8 :: tempb
      REAL*4 :: temp2b
      REAL*4 :: temp
      REAL*4 :: temp4

      !=================================================================
      ! INTERPOLATE_LUT2_ADJ begins here!
      !=================================================================

      ! Set the levels that were chosen in the look up table
      templev  = (/ 275.,  280.,   285.,   310.         /)
      jno2lev  = (/ 5.e-4, 0.0025, 0.0050, 0.012        /)
      cao3lev  = (/   5.,   20.,    35.,    75.         /)
      alfa0lev = (/ -90.,  -60.,   -45.,   -30.,
     &              -15.,    0.,    15.,    30.,
     &               45.,   60.,    75.,    90.         /)
      alfa5lev = (/ -90.,  -60.,   -45.,   -30.,
     &              -15.,    0.,    15.,    30.,
     &               45.,   60.,    75.,    90.         /)
      jo1dlev  = (/ 5.e-4, 0.0015, 0.0025, 0.0055       /)
      canoxlev = (/  10.,  200.,   1000.,  2000., 6000. /)

!      PRINT*,"Temperature levels are: ",templev
!      PRINT*,"This is grid cell: ",I,J

      ! Temperature
!      PRINT*,"Temperature here is: ",TS(I,J)
!      PRINT*,"USA: ",TS(32,64)

      ! Tracer concentrations in v/v
!      PRINT*,"[O3] is: ",STT(I,J,1,IDTO3)/ State_Met%AD(I,J,1) * TCVV(IDTO3)
!      PRINT*,"[CO] is: ",STT(I,J,1,IDTCO)/ State_Met%AD(I,J,1) * TCVV(IDTCO)
!      PRINT*,"IDTO3 is: ", IDTO3
!      PRINT*,"IDO3 is: ", IDO3
!      PRINT*,"In USA: ",STT(32,64,1,IDTO3)/State_Met%AD(32,64,1) * TCVV(IDTO3)

      ! SOLAR ZENITH ANGLES IN DEGREES
!      IJLOOP = ( (J-1) * IIPAR ) + I
!      PRINT*,"Local Time: ",GET_LOCALTIME(I)
!      PRINT*,"Solar Zenith Angle at this location: ",
!     $            ASIND(SUNCOS(IJLOOP))
!      IJLOOP = ( (64-1) * IIPAR ) + 32
!      PRINT*,"Local USA time: ", GET_LOCALTIME(32)
!      PRINT*,"Solar Zenith Angle at USA: ",
!     &  ASIND(SUNCOS(I,J))
!      PRINT*,"Solar Zenith Angle at USA - 5: ",
!     &  ASIND(SUNCOS_5hr(IJLOOP))

      ! Set the variables
      IJLOOP       = ( (J-1) * IIPAR ) + I
      var_array(1) = TS(I,J)                          ! Temperature
      var_array(2) = JNO2                             ! J(NO2), 1/s
      var_array(3) = o3 / dens * 1.E9                 ! [O3] in ppbv
      var_array(4) = ASIND(SUNCOS(IJLOOP))            ! alfa0
      var_array(5) = ASIND(SUNCOS_5hr(IJLOOP))        ! alfa5
      var_array(6) = JO1D / JNO2                      ! J(O1D)/J(NO2)
      var_array(7) = (no + no2) / dens * 1.E12        ! [NOx] in pptv

      ! prevent NaN when jvalues are 0.
      IF (jno2 .EQ. 0.) THEN
         var_array(6) = 0.
         CALL PUSHCONTROL1B(0)
      ELSE
         CALL PUSHCONTROL1B(1)
      ENDIF

      ! First some error checking
      ! ########### MAYBE CHECK HERE FOR NEGATIVE VALUES?##########

      !
      ! Determine reference index ( itemp,  ijno2, icao3, ialfa0,
      !                             ialfa5, ialfa, ijo1d, icaco )
      !
      !========================================================================
      ! Find smallest temperature reference level (i) for which actual
      ! temperature is smaller, then do
      !
      ! x(1) = ( temperature_level(i+1) - actual temperature   )
      !        -------------------------------------------------
      !        ( temperature_level(i+1) - temperature_level(i) )
      !
      ! then x(2) = 1.0 - x(1)
      !
      !========================================================================


      !---------------------
      ! Temperature:
      !---------------------
      temp_tmp = var_array(1)

      ! If temperature larger than largest in LUT, assign largest temp
      IF (var_array(1) .GT. templev(ntemp)) THEN
         temp_tmp = templev(ntemp)
         CALL PUSHCONTROL1B(0)
      ELSE
         CALL PUSHCONTROL1B(1)
      ENDIF

      ! If temp smaller, assign smallest temp level
      IF (var_array(1) .LT. templev(1)) THEN
         temp_tmp = templev(1)
         CALL PUSHCONTROL1B(1)
      ELSE
         CALL PUSHCONTROL1B(0)
      ENDIF

      DO i0=1,ntemp-1
         itemp = i0
         IF( templev( itemp+1 ) > temp_tmp ) EXIT
      END DO

      xtemp(1) = ( templev( itemp+1 ) - temp_tmp       ) /
     $           ( templev( itemp+1 ) - templev( itemp ) )
      xtemp(2) = 1.0 - xtemp(1)

      !---------------------
      ! J(NO2):
      !---------------------
      jno2_tmp = var_array(2)

      ! If larger than largest in LUT, assign largest level values
      IF (var_array(2) .GT. jno2lev(njno2)) THEN
         jno2_tmp = jno2lev(njno2)
         CALL PUSHCONTROL1B(0)
      ELSE
         CALL PUSHCONTROL1B(1)
      ENDIF

      ! If smaller, assign smallest level value
      IF (var_array(2) .LT. jno2lev(1)) THEN
         jno2_tmp = jno2lev(1)
         CALL PUSHCONTROL1B(1)
      ELSE
         CALL PUSHCONTROL1B(0)
      ENDIF

      DO i0=1,njno2-1
         ijno2 = i0
         IF( jno2lev( ijno2+1 ) > jno2_tmp ) EXIT
      END DO

      xjno2(1) = ( jno2lev( ijno2+1 ) - jno2_tmp       ) /
     $           ( jno2lev( ijno2+1 ) - jno2lev( ijno2 ) )
      xjno2(2) = 1.0 - xjno2(1)

      !---------------------
      ! [O3]:
      !---------------------
      cao3_tmp = var_array(3)

      ! If larger than largest in LUT, assign largest level values
      IF (var_array(3) .GT. cao3lev(ncao3)) THEN
         cao3_tmp = cao3lev(ncao3)
         CALL PUSHCONTROL1B(0)
      ELSE
         CALL PUSHCONTROL1B(1)
      ENDIF

      ! If smaller, assign smallest level value
      IF (var_array(3) .LT. cao3lev(1)) THEN
         cao3_tmp = cao3lev(1)
         CALL PUSHCONTROL1B(1)
      ELSE
         CALL PUSHCONTROL1B(0)
      ENDIF

      DO i0=1,ncao3-1
         icao3 = i0
         IF( cao3lev( icao3+1 ) > cao3_tmp ) EXIT
      END DO

      xcao3(1) = ( cao3lev( icao3+1 ) - cao3_tmp       ) /
     $           ( cao3lev( icao3+1 ) - cao3lev( icao3 ) )
      xcao3(2) = 1.0 - xcao3(1)

      !---------------------
      ! alfa0:
      !---------------------
      alfa0_tmp = var_array(4)

      ! If larger than largest in LUT, assign largest level values
      IF (var_array(4) .GT. alfa0lev(nalfa0)) THEN
         alfa0_tmp = alfa0lev(nalfa0)
         CALL PUSHCONTROL1B(0)
      ELSE
         CALL PUSHCONTROL1B(1)
      ENDIF

      ! If smaller, assign smallest level value
      IF (var_array(4) .LT. alfa0lev(1)) THEN
         alfa0_tmp = alfa0lev(1)
         CALL PUSHCONTROL1B(1)
      ELSE
         CALL PUSHCONTROL1B(0)
      ENDIF

      DO i0=1,nalfa0-1
         ialfa0 = i0
         IF( alfa0lev( ialfa0+1 ) > alfa0_tmp ) EXIT
      END DO

      xalfa0(1) = ( alfa0lev( ialfa0+1 ) - alfa0_tmp        ) /
     $            ( alfa0lev( ialfa0+1 ) - alfa0lev( ialfa0 ) )
      xalfa0(2) = 1.0 - xalfa0(1)

      !---------------------
      ! alfa5:
      !---------------------
      alfa5_tmp = var_array(5)

      ! If larger than largest in LUT, assign largest level values
      IF (var_array(5) .GT. alfa5lev(nalfa5)) THEN
         alfa5_tmp = alfa5lev(nalfa5)
         CALL PUSHCONTROL1B(0)
      ELSE
         CALL PUSHCONTROL1B(1)
      END IF

      ! If smaller, assign smallest level value
      IF (var_array(5) .LT. alfa5lev(1)) THEN
         alfa5_tmp = alfa5lev(1)
         CALL PUSHCONTROL1B(1)
      ELSE
         CALL PUSHCONTROL1B(0)
      END IF

      DO i0=1,nalfa5-1
         ialfa5 = i0
         IF( alfa5lev( ialfa5+1 ) > alfa5_tmp ) EXIT
      END DO

      xalfa5(1) = ( alfa5lev( ialfa5+1 ) - alfa5_tmp        ) /
     $            ( alfa5lev( ialfa5+1 ) - alfa5lev( ialfa5 ) )
      xalfa5(2) = 1.0 - xalfa5(1)

      !---------------------
      ! jo1d:
      !---------------------
      jo1d_tmp = var_array(6)

      ! If larger than largest in LUT, assign largest level values
      IF (var_array(6) .GT. jo1dlev(njo1d)) THEN
         jo1d_tmp = jo1dlev(njo1d)
         CALL PUSHCONTROL1B(0)
      ELSE
         CALL PUSHCONTROL1B(1)
      ENDIF

      ! If smaller, assign smallest level value
      IF (var_array(6) .LT. jo1dlev(1)) THEN
         jo1d_tmp = jo1dlev(1)
         CALL PUSHCONTROL1B(1)
      ELSE
         CALL PUSHCONTROL1B(0)
      ENDIF

      DO i0=1,njo1d-1
         ijo1d = i0
         IF( jo1dlev( ijo1d+1 ) > jo1d_tmp ) EXIT
      END DO

      xjo1d(1) = ( jo1dlev( ijo1d+1 ) - jo1d_tmp       ) /
     $           ( jo1dlev( ijo1d+1 ) - jo1dlev( ijo1d ) )
      xjo1d(2) = 1.0 - xjo1d(1)

      !---------------------
      ! [NOx]:
      !---------------------
      canox_tmp = var_array(7)

      ! If larger than largest in LUT, assign largest level values
      IF (var_array(7) .GT. canoxlev(ncanox)) THEN
         canox_tmp = canoxlev(ncanox)
         CALL PUSHCONTROL1B(0)
      ELSE
         CALL PUSHCONTROL1B(1)
      ENDIF

      ! If smaller, assign smallest level value
      IF (var_array(7) .LT. canoxlev(1)) THEN
         canox_tmp = canoxlev(1)
         CALL PUSHCONTROL1B(1)
      ELSE
         CALL PUSHCONTROL1B(0)
      ENDIF

      DO i0=1,ncanox-1
         icanox = i0
         IF( canoxlev( icanox+1 ) > canox_tmp ) EXIT
      END DO

      xcanox(1) = ( canoxlev( icanox+1 ) - canox_tmp        ) /
     $            ( canoxlev( icanox+1 ) - canoxlev( icanox ) )
      xcanox(2) = 1.0 - xcanox(1)

!      PRINT*,"The i-values are:",    itemp, ijno2, icao3, ialfa0,
!     $                               ialfa5, ijo1d, icanox
!      PRINT*,"Variables are: ",      var_array
!      PRINT*,"For testing, xtemp: ", xtemp

      !======================
      ! Linear interpolation
      !======================


      xcanoxb = 0.0_4
      xtempb = 0.0_4
      xcao3b = 0.0_4
      xjo1db = 0.0_4
      xalfa0b = 0.0_4
      xjno2b = 0.0_4
      xalfa5b = 0.0_4

      DO i1=2,1,-1
      DO i2=2,1,-1
      DO i3=2,1,-1
      DO i4=2,1,-1
      DO i5=2,1,-1
      DO i6=2,1,-1
      DO i7=2,1,-1

         !IF ENCOUNTER -999 IN THE LUT PRINT ERROR!!
         IF ( ( fracnox( itemp+i1-1,  ijno2+i2-1,  icao3+i3-1,
     $                   ialfa0+i4-1, ialfa5+i5-1, ijo1d+i6-1,
     $                   icanox+i7-1  ) < 0. )     .or.
     $        ( fracnox( itemp+i1-1,  ijno2+i2-1,  icao3+i3-1,
     $                   ialfa0+i4-1, ialfa5+i5-1, ijo1d+i6-1,
     $                   icanox+i7-1  ) > 1. ) )   THEN

            PRINT*, 'INTERPOLATE_LUT2_ADJ: fracnox = ',
     $          fracnox( itemp+i1-1,  ijno2+i2-1,  icao3+i3-1,
     $                   ialfa0+i4-1, ialfa5+i5-1, ijo1d+i6-1,
     $                   icanox+i7-1  )

            MSG = 'LUT error: Fracnox should be between 0 and 1!'
            CALL ERROR_STOP( MSG,
     $          'INTERPOLATE_LUT2_ADJ ("paranox_adj_mod.F")' )
         ENDIF

                  temp1  =  xalfa5(i5)  *  xjo1d(i6)
                  temp0  =  xcao3(i3)   *  xalfa0(i4)
                  temp   =  xtemp(i1)   *  xjno2(i2)

                  tempb2 = fracnox( itemp+i1-1,  ijno2+i2-1, icao3+i3-1,
     &                              ialfa0+i4-1, ialfa5+i5-1,ijo1d+i6-1,
     &                              icanox+i7-1 ) * fraction_noxb

                  tempb0 =  temp0  *  temp1      *  tempb2
                  tempb1 =  temp   *  xcanox(i7) *  tempb2
                  temp4  =  xalfa5(i5) *  xjo1d(i6)
                  temp3  =  xcao3(i3)  *  xalfa0(i4)
                  temp2  =  xtemp(i1)  *  xjno2(i2)

                  temp2b = intope( itemp+i1-1,   ijno2+i2-1, icao3+i3-1,
     &                             ialfa0+i4-1,  ialfa5+i5-1,ijo1d+i6-1,
     &                             icanox+i7-1 ) * int_opeb

                  temp2b0 = temp3  *  temp4      *  temp2b
                  temp2b1 = temp2  *  xcanox(i7) *  temp2b

                  xtempb(i1)  = xtempb(i1)
     &                        + xcanox(i7) * xjno2(i2) * tempb0
     &                        + xcanox(i7) * xjno2(i2) * temp2b0

                  xjno2b(i2)  = xjno2b(i2)
     &                        + xcanox(i7) * xtemp(i1) * tempb0
     &                        + xcanox(i7) * xtemp(i1) * temp2b0

                  xcanoxb(i7) = xcanoxb(i7)
     &                        + temp  * tempb0
     &                        + temp2 * temp2b0

                  xcao3b(i3)  = xcao3b(i3)
     &                        + temp1 * xalfa0(i4) * tempb1
     &                        + temp4 * xalfa0(i4) * temp2b1

                  xalfa0b(i4) = xalfa0b(i4)
     &                        + temp1 * xcao3(i3) * tempb1
     &                        + temp4 * xcao3(i3) * temp2b1

                  xalfa5b(i5) = xalfa5b(i5)
     &                        + temp0 * xjo1d(i6) * tempb1
     &                        + temp3 * xjo1d(i6) * temp2b1

                  xjo1db(i6)  = xjo1db(i6)
     &                        + temp0 * xalfa5(i5) * tempb1
     &                        + temp3 * xalfa5(i5) * temp2b1

      END DO
      END DO
      END DO
      END DO
      END DO
      END DO
      END DO

      !---------------------
      ! [NOx]: ADJOINT
      !---------------------
      xcanoxb(1) = xcanoxb(1) - xcanoxb(2)
      xcanoxb(2) = 0.0_4
      canox_tmpb = -( xcanoxb(1)/
     &              ( canoxlev(icanox+1) - canoxlev(icanox) ) )

      CALL POPCONTROL1B(branch)
      IF (branch .NE. 0) canox_tmpb = 0.0_4

      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) canox_tmpb = 0.0_4

      var_arrayb = 0.0_4
      var_arrayb(7) = var_arrayb(7) + canox_tmpb

      !---------------------
      ! jo1d: ADJOINT
      !---------------------
      xjo1db(1) = xjo1db(1) - xjo1db(2)
      xjo1db(2) = 0.0_4
      jo1d_tmpb = -( xjo1db(1)/
     &             ( jo1dlev(ijo1d+1) - jo1dlev(ijo1d) ) )

      CALL POPCONTROL1B(branch)
      IF (branch .NE. 0) jo1d_tmpb = 0.0_4

      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) jo1d_tmpb = 0.0_4

      var_arrayb(6) = var_arrayb(6) + jo1d_tmpb

      !---------------------
      ! alfa5: ADJOINT
      !---------------------
      xalfa5b(1) = xalfa5b(1) - xalfa5b(2)
      xalfa5b(2) = 0.0_4

      alfa5_tmpb = -( xalfa5b(1)/
     &              ( alfa5lev(ialfa5+1) - alfa5lev(ialfa5) ) )

      CALL POPCONTROL1B(branch)
      IF (branch .NE. 0) alfa5_tmpb = 0.0_4

      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) alfa5_tmpb = 0.0_4

      var_arrayb(5) = var_arrayb(5) + alfa5_tmpb

      !---------------------
      ! alfa0: ADJOINT
      !---------------------
      xalfa0b(1) = xalfa0b(1) - xalfa0b(2)
      xalfa0b(2) = 0.0_4

      alfa0_tmpb = -( xalfa0b(1)/
     &              ( alfa0lev(ialfa0+1) - alfa0lev(ialfa0) ) )

      CALL POPCONTROL1B(branch)
      IF (branch .NE. 0) alfa0_tmpb = 0.0_4

      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) alfa0_tmpb = 0.0_4

      var_arrayb(4) = var_arrayb(4) + alfa0_tmpb

      !---------------------
      ! [O3]: ADJOINT
      !---------------------
      xcao3b(1) = xcao3b(1) - xcao3b(2)
      xcao3b(2) = 0.0_4

      cao3_tmpb = -( xcao3b(1)/
     &             ( cao3lev(icao3+1) - cao3lev(icao3) ) )

      CALL POPCONTROL1B(branch)
      IF (branch .NE. 0) cao3_tmpb = 0.0_4

      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) cao3_tmpb = 0.0_4

      var_arrayb(3) = var_arrayb(3) + cao3_tmpb

      !---------------------
      ! J(NO2): ADJOINT
      !---------------------
      xjno2b(1) = xjno2b(1) - xjno2b(2)
      xjno2b(2) = 0.0_4

      jno2_tmpb = -( xjno2b(1)/
     &             ( jno2lev(ijno2+1) - jno2lev(ijno2) ) )

      CALL POPCONTROL1B(branch)
      IF (branch .NE. 0) jno2_tmpb = 0.0_4

      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) jno2_tmpb = 0.0_4
      var_arrayb(2) = var_arrayb(2) + jno2_tmpb

      !---------------------
      ! Temperature: ADJOINT
      !---------------------
      xtempb(1) = xtempb(1) - xtempb(2)
      xtempb(2) = 0.0_4

      temp_tmpb = -( xtempb(1)/
     &             ( templev(itemp+1) - templev(itemp) ) )

      CALL POPCONTROL1B(branch)
      IF (branch .NE. 0) temp_tmpb = 0.0_4

      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) temp_tmpb = 0.0_4

      var_arrayb(1) = var_arrayb(1) + temp_tmpb

      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) var_arrayb(6) = 0.0_4


      tempb = 1.e12*var_arrayb(7)/dens
      nob = tempb
      no2b = tempb

      var_arrayb(7) = 0.0_4
      var_arrayb(6) = 0.0_4
      var_arrayb(5) = 0.0_4
      var_arrayb(4) = 0.0_4

      o3b = 1.e9*var_arrayb(3)/dens

      !int_opeb = 0.0_4
      !fraction_noxb = 0.0_4

      END SUBROUTINE INTERPOLATE_LUT2_ADJ

      END MODULE PARANOX_ADJ_MOD


