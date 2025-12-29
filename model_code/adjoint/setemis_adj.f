! $Id: setemis_adj.f,v 1.1 2010/04/01 07:09:43 daven Exp $
      SUBROUTINE SETEMIS_ADJ( )
!
!******************************************************************************
!  Subroutine SETEMIS_ADJ passes adjoints from SMVGEAR array adjoint REMIS_ADJ
!  back to GEOS-Chem emission adjoint arrays, e.g., BIOFUEL_ADJ
!
!  Based on forward code SETEMIS (lwh, jyl, gmg, djj, bdf, bmy, 6/8/98, 6/11/08)
!
!  Variables taken from F90 Modules:
!  ============================================================================
!  (1 ) BIOFUEL  (REAL*8 ) : Biofuel burning emissions  [molec (C)/cm3/s    ]
!  (2 ) BFTRACE  (INTEGER) : Index array for biofuels   [CTM tracer #       ]
!  (3 ) NBFTRACE (INTEGER) : Number of biofuel species  [unitless           ]
!  (4 ) BURNEMIS (REAL*8 ) : Biomass burning emissions  [molec (C)/cm3/s    ]
!  (5 ) BIOTRCE  (INTEGER) : Index array for bioburn    [CTM tracer #       ]
!  (6 ) NBIOTRCE (INTEGER) : Number of bioburn species  [unitless           ]
!  (7 ) JLOP     (INTEGER) : SMVGEAR grid box index     [unitless           ]
!  (8 ) REMIS    (REAL*8 ) : SMVGEAR emissions array    [molec species/cm3/s]
!  (9 ) VOLUME   (REAL*8 ) : SMVGEAR volume array       [cm3                ]
!
!  NOTES:
!  (1 ) Original code from Harvard Tropospheric Chemistry Module for 3-D
!        applications by Larry Horowitz, Jinyou Liang, Gerry Gardner,
!        Prof. Daniel Jacob of Harvard University (Release V2.0)
!  (2 ) New version 3.0 by Bob Yantosca to place NOx emissions into boxes
!        above the surface. (bmy, 6/8/98)
!  (3 ) Also now do chemistry up to the location of the annual mean
!         tropopause (bmy, 12/9/99)
!  (4 ) BURNEMIS is now dynamically allocatable and is contained in F90
!        module "biomass_mod.f".  BIOTRCE and NBIOTRCE are also contained
!        in "biomass_mod.f".  (bmy, 9/12/00)
!  (5 ) BIOFUEL is now dynamically allocatable and is contained in F90
!        module "biofuel_mod.f".  BFTRACE and NBFTRACE are also contained
!        in "biofuel_mod.f" (bmy, 9/12/00, 4/17/01)
!  (6 ) BURNEMIS and BIOFUEL are now treated as true global arrays,
!        and need to be referenced by the global offset variables
!        IREF = I + I0 and JREF = J + J0 (bmy, 9/12/00)
!  (7 ) Now reference JLOP, REMIS, VOLUME from F90 module "comode_mod.f",
!        in order to save memory (bmy, 10/19/00)
!  (8 ) Now add in up to NBFTRACE biofuel species (bmy, 4/17/01)
!  (9 ) Add new subroutine header, updated comments, cosmetic changes.
!        (bmy, 4/17/01)
!  (10) Updated comments -- GEMISNOX is [molec/cm3/s]. (bdf, bmy, 6/7/01)
!  (11) For GEOS-3, we now distributes surface emissions throughout the
!        boundary layer.  This is necessary since the first couple of GEOS-3
!        surface layers are very thin.  Piling up of emissions into a small
!        layer will cause SMVGEAR to choke.  (bdf, bmy, 6/15/01)
!  (12) Also now reference BFTRACE and NBFTRACE from "biofuel_mod.f",
!        and reference AD12 from "diag_mod.f". (bdf, bmy, 6/15/01)
!  (13) For GEOS-1, GEOS-STRAT, emit into the surface layer, as we did
!        in prior versions. (bmy, 6/26/01)
!  (14) Bug fix: corrected a typo for the biofuel emissions (bmy, 7/10/01)
!  (15) Bug fix: make sure BIOMASS and BIOFUEL, and SOIL NOx emissions have
!        units of [molec/box/s] before distributing thru the boundary layer.
!        This involves multiplication by VOLUME(JLOOP1) and division by
!        VOLUME(JLOOP). (bmy, 7/16/01)
!  (16) XTRA2(IREF,JREF,5) is now XTRA2(I,J).  BIOFUEL(:,IREF,JREF) is now
!        BIOFUEL(:,I,J). BURNEMIS(:,IREF,JREF) is now BURNEMIS(:,I,J).
!        Replace PW(I,J) with P(I,J). (bmy, 9/28/01)
!  (17) Removed obsolete code from 9/01 (bmy, 10/24/01)
!  (18) Now references GET_PEDGE from "pressure_mod.f", to compute P at
!        the bottom edge of grid box (I,J,L).  (dsa, bdf, bmy, 8/21/02)
!  (19) Now reference IDTNOX, IDENOX, etc from "tracerid_mod.f" (bmy, 11/6/02)
!  (20) Remove references to IREF, JREF (bmy, 2/11/03)
!  (21) NEMIS is now NEMIS(NCS) for SMVGEAR II (gcc, bdf, bmy, 4/1/03)
!  (22) Added parallel loop over N.  Also directly substituted JLOP(I,J,1)
!        for all instances of JLOOP1.  Updated comments. (hamid, bmy, 3/19/04)
!  (23) Bug fix for COMPAQ compiler...do not use EXIT from w/in parallel loop.
!        (auvray, bmy, 11/29/04)
!  (24) Now replace XTRA2 with GET_PBL_TOP_L in "pbl_mix_mod.f".  Now remove
!        reference to CMN, it's obsolete.  Now references GET_TPAUSE_LEVEL
!        from "tropopause_mod.f" (bmy, 8/22/05)
!  (25) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (26) Now updated for new "biomass_mod.f" (bmy, 4/5/06)
!  (27) Now account for the different definition of tropopause in case
!        of variable tropopause.   The BIOMASS array from "biomass_mod.f" is
!        now in units of [molec CO/cm2/s].  Adjust unit conversion accordingly.
!        Also replace NBIOMAX with NBIOMAX_GAS, since aerosol biomass is
!        handled elsewhere.  (bdf, phs, bmy, 9/27/06)
!  (28) Now replace GEMISNOX array (from CMN_NOX) with module arrays
!        EMIS_LI_NOx and EMIS_AC_NOx (ltm, bmy, 10/3/07)
!  (29) Bug fix: resize EMISRR to be consistent w/ CMN_O3 (bmy, jaf, 6/11/08)
!******************************************************************************
!
      ! References to F90 modules
      USE AIRCRAFT_NOX_MOD,  ONLY : EMIS_AC_NOx
      USE BIOFUEL_MOD,       ONLY : BIOFUEL,   BFTRACE, NBFTRACE
      USE BIOMASS_MOD,       ONLY : BIOMASS,   BIOTRCE, NBIOMAX_GAS
      USE COMODE_MOD,        ONLY : JLOP,      REMIS,   VOLUME
      USE COMODE_MOD,        ONLY : IYSAVE
      USE DIAG_MOD,          ONLY : AD12
      USE GRID_MOD,          ONLY : GET_AREA_CM2
      USE LOGICAL_MOD,       ONLY : LVARTROP
      USE LIGHTNING_NOX_MOD, ONLY : EMIS_LI_NOx
      USE PBL_MIX_MOD,       ONLY : GET_PBL_TOP_L
      USE PRESSURE_MOD,      ONLY : GET_PEDGE
      USE TRACERID_MOD,      ONLY : CTRMB,     IDEMIS,  IDENOX
      USE TROPOPAUSE_MOD,    ONLY : GET_TPAUSE_LEVEL

      ! Ship plume emissions moved from calcrate.F
      USE EDGAR_MOD,         ONLY : GET_EDGAR_NOx
      USE ICOADS_SHIP_MOD,   ONLY : GET_ICOADS_SHIP
      USE EMEP_MOD,          ONLY : GET_EMEP_ANTHRO, GET_EUROPE_MASK
      USE EMISSIONS_MOD,     ONLY : NOx_SCALING
      USE COMODE_MOD,        ONLY : CSPEC, AIRDENS, CSPEC_ADJ
      USE COMODE_MOD,        ONLY : CSPEC_ORIG, CSPEC_FOR_KPP
      USE LOGICAL_MOD,       ONLY : LDRYD,    LPRT
      USE LOGICAL_MOD,       ONLY : LICOADSSHIP, LEDGARSHIP, LEMEPSHIP
      USE DRYDEP_MOD,        ONLY : SHIPO3DEP
      USE ADJ_ARRAYS_MOD,    ONLY : SHIPO3DEP_ADJ
      USE TRACERID_MOD,      ONLY : IDTNO,  IDTNOX, IDTHNO3
      USE TRACERID_MOD,      ONLY : IDO3,   IDHNO3, IDNO2, IDNO
      USE TRACERID_MOD,      ONLY : IDEHNO3, IDEOX
      USE DAO_MOD,           ONLY : BXHEIGHT, SUNCOS, SUNCOS_5hr
      USE PARANOX_MOD,       ONLY : INTERPOLATE_LUT2
      USE PARANOX_ADJ_MOD,   ONLY : INTERPOLATE_LUT2_ADJ
      USE COMODE_MOD,        ONLY : CHK_CSPEC

      USE ADJ_ARRAYS_MOD,    ONLY : NADJ_EANTHRO
      USE ADJ_ARRAYS_MOD,    ONLY : NADJ_EBIOMASS
      USE ADJ_ARRAYS_MOD,    ONLY : NADJ_EBIOFUEL
      USE ADJ_ARRAYS_MOD,    ONLY : IDADJ_ENOX_so
      USE ADJ_ARRAYS_MOD,    ONLY : IDADJ_ENOX_li
      USE ADJ_ARRAYS_MOD,    ONLY : IDADJ_ENOX_ac
      USE ADJ_ARRAYS_MOD,    ONLY : EMS_SF_ADJ
      USE ADJ_ARRAYS_MOD,    ONLY : REMIS_ADJ

      IMPLICIT NONE

#     include "CMN_SIZE"  ! Size parameters
#     include "CMN_NOX"   ! GEMISNOX2
#     include "CMN_O3"    ! EMISRR, EMISRRN
#     include "comode.h"  ! IDEMS, NEMIS

      ! Local variables
      LOGICAL             :: IS_LI_NOx, IS_AC_NOx
      INTEGER             :: I, J,  JLOOP, JLOOP1, LTROP
      INTEGER             :: L, LL, N, NN,  NBB, NBF, TOP
      REAL*8              :: COEF1,   TOTPRES, DELTPRES
      REAL*8              :: EMIS_BL, NOXTOT,  TOTAL, A_CM2


      REAL*8              :: SHIP, TOTO3, JNO2, JO1D
      REAL*8              :: DTSRCE, AREA_CM2
      REAL*4              :: FRACTION_NOX, INT_OPE
      INTEGER             :: NK
      CHARACTER*8         :: SPECNAME
      REAL*8, EXTERNAL    :: FJFUNC

      REAL*8              :: EMIS_BL_ADJ
      REAL*4              :: FRACTION_NOX_ADJ
      REAL*4              :: INT_OPE_ADJ
      REAL*8              :: SHIP_ADJ
      REAL*8              :: TOTO3_ADJ
      REAL*8              :: O3_ADJ
      REAL*8              :: NO_ADJ
      REAL*8              :: NO2_ADJ

      INTEGER             :: M

      !=================================================================
      ! SETEMIS_ADJ begins here!
      !=================================================================

      ! Test if the EMIS_LI_NOx and EMIS_AC_NOx arrays are allocated
      IS_LI_NOx = ALLOCATED( EMIS_LI_NOx )
      IS_AC_NOX = ALLOCATED( EMIS_AC_NOX )

      M = 1

      NCS = NCSURBAN

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( N,       NN,     NBB,      NBF,          I       )
!$OMP+PRIVATE( J,       L,      JLOOP,    COEF1,        TOP     )
!$OMP+PRIVATE( TOTPRES, NOXTOT, DELTPRES, EMIS_BL,      A_CM2   )
!$OMP+PRIVATE( TOTO3,   SHIP,   AREA_CM2, FRACTION_NOX, INT_OPE )
!$OMP+PRIVATE( NK,      JNO2,   JO1D,     SPECNAME              )
!$OMP+PRIVATE( EMIS_BL_ADJ,     TOTO3_ADJ,INT_OPE_ADJ           )
!$OMP+PRIVATE( FRACTION_NOX_ADJ,SHIP_ADJ                        )
!$OMP+PRIVATE( O3_ADJ,  NO_ADJ, NO2_ADJ                         )
!$OMP+SCHEDULE( DYNAMIC )

      ! Loop over emission species
      DO N = 1, NEMIS(NCS)

         ! Get CTM tracer number NN corresponding to emission species N
         NN = IDEMS(N)
         IF ( NN == 0 ) CYCLE

         ! We have to search for the biomass burning species in
         ! BIOTRCE with the same CTM tracer number NN as in IDEMS
         NBB = 0
         IF ( ALLOCATED( BIOMASS ) ) THEN
            DO I = 1, NBIOMAX_GAS
               IF ( BIOTRCE(I) == NN ) THEN
                  NBB = I
#if   defined( COMPAQ )
                  ! COMPAQ has an issue with EXIT from w/in parallel loop
                  ! (auvray, bmy, 11/29/04)
#else
                  EXIT
#endif
               ENDIF
            ENDDO
         ENDIF

         ! We have to search for the biofuel burning species in
         ! BFTRACE with the same CTM tracer number NN as in IDEMS
         NBF = 0
         IF ( ALLOCATED( BIOFUEL ) ) THEN
            DO I = 1, NBFTRACE
               IF ( BFTRACE(I) == NN ) THEN
                  NBF = I
#if   defined( COMPAQ )
                  ! COMPAQ has an issue with EXIT from w/in parallel loop
                  ! (auvray, bmy, 11/29/04)
#else
                  EXIT
#endif
              ENDIF
            ENDDO
         ENDIF

         ! COEF1 = molecules of emission species / molecules of tracer
         COEF1 = 1.0 + CTRMB(NN, IDEMIS(NN))

         ! Loop over Lat and Long boxes
         DO J = 1, NLAT
         DO I = 1, NLONG

            !===========================================================
            ! For GEOS-3: distribute surface emissions throughout the
            ! entire boundary layer.  Define some variables here.
            ! (bdf, 6/15/01)
            !===========================================================

            ! Top level of the boundary layer
            ! guard for b.l. being in first level.
            TOP = FLOOR( GET_PBL_TOP_L( I, J ) )
            IF ( TOP == 0 ) TOP = 1

            ! Pressure thickness of entire boundary layer [hPa]
            TOTPRES = GET_PEDGE(I,J,1) - GET_PEDGE(I,J,TOP+1)

            !===========================================================
            ! Adjoint of biofuel burning source [molec/cm3/s]
            !===========================================================
            IF ( NBF /= 0 ) THEN
               DO L = 1, TOP
                  JLOOP   = JLOP(I,J,L)

                  IF ( JLOOP /= 0 ) THEN

                     ! fwd code:
                     !REMIS(JLOOP,N) = REMIS(JLOOP,N) +
                     !                 ( EMIS_BL / VOLUME(JLOOP) )
                     ! adj code:
                     EMIS_BL_ADJ = REMIS_ADJ(JLOOP,N) / VOLUME(JLOOP)

                     ! recalc DELTPRES
                     ! Thickness of level L [mb]
                     DELTPRES = GET_PEDGE(I,J,L) - GET_PEDGE(I,J,L+1)

                     ! fwd code:
                     ! Store in EMIS_BL.
                     !EMIS_BL  = ( BIOFUEL(NBF,I,J) *
                     !             VOLUME( JLOP(I,J,1) ) / COEF1 ) *
                     !             ( DELTPRES / TOTPRES )
                     !         *   BIOFUEL_ICS(I,J,NBF)
                     IF ( NADJ_EBIOFUEL(NN) > 0 ) THEN
                        EMS_SF_ADJ(I,J,M,NADJ_EBIOFUEL(NN))
     &                            = EMS_SF_ADJ(I,J,M,NADJ_EBIOFUEL(NN))
     &                            + ( BIOFUEL(NBF,I,J)
     &                            * VOLUME( JLOP(I,J,1) ) / COEF1 )
     &                            * ( DELTPRES / TOTPRES )
     &                            * EMIS_BL_ADJ
                     ENDIF

                  ENDIF
                  EMIS_BL_ADJ = 0d0
               ENDDO
            ENDIF

            !===========================================================
            ! Adjoint of biomass burning source [molec/cm3/s]
            !===========================================================
            IF ( NBB /= 0 ) THEN
               DO L = 1, TOP
                  JLOOP   = JLOP(I,J,L)

                  IF ( JLOOP /= 0 ) THEN

                     ! fwd code:
                     !REMIS(JLOOP,N) = REMIS(JLOOP,N) +
     &               !                 ( EMIS_BL / VOLUME(JLOOP) )
                     ! adj code:
                     EMIS_BL_ADJ = REMIS_ADJ(JLOOP,N) / VOLUME(JLOOP)

                     ! Thickness of level L [mb]
                     DELTPRES = GET_PEDGE(I,J,L) - GET_PEDGE(I,J,L+1)

                     ! Grid box area [cm2]
                     A_CM2    = GET_AREA_CM2( IYSAVE(JLOOP) )

                     ! fwd code:
                     !EMIS_BL  = ( BIOMASS(I,J,NBB) * A_CM2 / COEF1   ) *
                     !           ( DELTPRES                 / TOTPRES )
                     ! adj code:
                     IF ( NADJ_EBIOMASS(NN) > 0 ) THEN
                        EMS_SF_ADJ(I,J,M,NADJ_EBIOMASS(NN))
     &                        =   EMS_SF_ADJ(I,J,M,NADJ_EBIOMASS(NN))
     &                        +   BIOMASS(I,J,NBB) * A_CM2 / COEF1
     &                        *   ( DELTPRES / TOTPRES )
     &                        *   EMIS_BL_ADJ
                     ENDIF


                  ENDIF
                  EMIS_BL_ADJ = 0d0
               ENDDO
            ENDIF

            !===========================================================
            ! Adjoints of non-NOx sources
            !===========================================================
            IF ( N /= IDENOX ) THEN

               !========================================================
               ! Anthropogenic tracers other than NOx [molec/box/s]
               ! Distribute emissions thru the entire boundary layer
               !========================================================
               DO L = 1, TOP
                  JLOOP   = JLOP(I,J,L)

                  IF ( JLOOP /= 0 ) THEN

                     ! Thickness of level L [mb]
                     DELTPRES = GET_PEDGE(I,J,L) - GET_PEDGE(I,J,L+1)

                     ! fwd code:
                     !REMIS(JLOOP,N) = EMIS_BL / VOLUME(JLOOP)
                     ! adj code:
                     EMIS_BL_ADJ = REMIS_ADJ(JLOOP,N) / VOLUME(JLOOP)

                     ! fwd code:
                     !EMIS_BL        = ( EMISRR(I,J,N) / COEF1   ) *
                     !                 ( DELTPRES      / TOTPRES )
                     ! adj code:
                     IF ( NADJ_EANTHRO(NN) > 0 ) THEN
                        EMS_SF_ADJ(I,J,M,NADJ_EANTHRO(NN))
     &                              = EMS_SF_ADJ(I,J,M,NADJ_EANTHRO(NN))
     &                              + ( EMISRR(I,J,N) / COEF1   )
     &                              * ( DELTPRES      / TOTPRES )
     &                              * EMIS_BL_ADJ
                     ENDIF
                  ENDIF
               ENDDO

               ! fwd code:
               !EMIS_BL = 0d0
               ! adj code:
               EMIS_BL_ADJ = 0d0

            ! For NOx only....
            ELSEIF( N == IDENOX ) THEN

               !========================================================
               ! Ship NOx (emitted as NO, O3, HNO3 after plume evolution)
               ! [molec/cm3/s]
               ! Distribute emissions thru the entire boundary layer
               !========================================================

               ! Update only if ship emissions are turned on
               IF ( LICOADSSHIP .or. LEDGARSHIP .or. LEMEPSHIP ) THEN

                  ! Recompute EMIS_BL, SHIP, FRACTION_NOX, INT_OPE

                  ! Surface area of grid box
                  AREA_CM2 = GET_AREA_CM2( J )

                  ! Reset
                  SHIP = 0D0

                  ! handle global inventory first
                  IF ( LEDGARSHIP ) THEN

                      ! Get SHIP EDGAR emissions for NOx [molec/cm2/s]
                     SHIP = GET_EDGAR_NOx( I, J,
     &                      MOLEC_CM2_S=.TRUE., SHIP=.TRUE.)

                  ! ICOADS ship emissions (cklee,7/09/09)
                  ELSE IF ( LICOADSSHIP ) THEN

                     ! Get ICOADS  emissions for NOx [molec/cm2/s]
                     SHIP = GET_ICOADS_SHIP( I, J, IDTNOX,
     &                      MOLEC_CM2_S=.TRUE. )

                  ENDIF

                  ! Overwrite Europe
                  IF ( LEMEPSHIP ) THEN

                     ! Prevent overwriting ICOADS data with a 0 value from EMEP
                     IF (GET_EUROPE_MASK( I, J ) > 0d0) THEN

                        ! Get SHIP EMEP emissions for NOx [molec/cm2/s]
                        SHIP = GET_EMEP_ANTHRO( I, J, IDTNOX,
     &                         SHIP=.TRUE. )
                     ENDIF

                  ENDIF

                  ! Add possible scaling of NOx emissions
                  SHIP = SHIP * NOx_SCALING

                  !----------------------------------
                  ! Get J-Values for J(NO2) and J(O3)
                  !----------------------------------

                  ! Check if sun is up
                  !Need to replace with SUNCOS.  Careful -- SUNCOS has JLOOP 1D index.
                  !Need to convert from I,J
                  JLOOP = JLOP(I,J,1)
                  IF (SUNCOS(JLOOP) > 0d0 ) THEN

                     ! Loop over photolysis reactions to find NO2, O3
                     DO L = 1, JPHOTRAT(NCS)

                        ! Reaction number
                        NK  = NRATES(NCS) + L

                        ! Name of species being photolyzed
                        SPECNAME = NAMEGAS(IRM(1,NK,NCS))

                        ! Check if this is NO2 or O3, store values, 1/s
                        SELECT CASE ( TRIM( SPECNAME ) )
                        CASE ( 'NO2' )
                           JNO2 = FJFUNC(I,J,1,L,1,SPECNAME)
                        CASE ( 'O3' )
                           JO1D = FJFUNC(I,J,1,L,1,SPECNAME)
                        CASE DEFAULT
                        END SELECT

                     ENDDO

                  ELSE

                     ! J-values are zero when sun is down
                     JNO2 = 0d0
                     JO1D = 0d0

                  ENDIF

                  ! Determine fraction of NOx remaining and integrated Ozone
                  ! Production Efficiency for ship emiss (gvinken,mpayer,2/7/12)
                  ! Uses surface-layer concentrations of O3, NO, NO2 (molec/cm3)
                  ! from CSPEC and air density (molec/cm3)
                  ! These values of CSPEC need to be forward model values from
                  ! before KPP but after PARTITION
                  CALL INTERPOLATE_LUT2( I, J,
     &                 CSPEC_ORIG(JLOOP,IDO3),  CSPEC_ORIG(JLOOP,IDNO),
     &                 CSPEC_ORIG(JLOOP,IDNO2), AIRDENS(JLOOP),
     &                 JO1D,               JNO2,
     &                 FRACTION_NOx,       INT_OPE)

                  !-----------------
                  ! Ship O3
                  !-----------------

                  EMIS_BL = SHIP * (1D0 - FRACTION_NOX) * INT_OPE

                  !----------------------------
                  ! adjoint code begins here
                  !----------------------------
                  EMIS_BL_ADJ      = 0d0
                  TOTO3_ADJ        = 0d0
                  INT_OPE_ADJ      = 0d0
                  FRACTION_NOX_ADJ = 0d0
                  SHIP_ADJ         = 0d0

                  ! Ship plume chemistry can create or destroy net O3
                  ! Treat positive O3 production as emissions
                  ! Treat O3 *destruction* in plume as dry deposition
                  ! (cdh, 3/21/2013)
                  IF (EMIS_BL >= 0d0) THEN

                     ! Loop over the boundary layer
                     DO L = 1, TOP
                        JLOOP   = JLOP(I,J,L)

                        IF ( JLOOP == 0 ) CYCLE

                        ! Thickness of level L [mb]
                        DELTPRES = GET_PEDGE(I,J,L) - GET_PEDGE(I,J,L+1)


                        ! fwd code:
                        !REMIS(JLOOP,IDEOX) = REMIS(JLOOP,IDEOX)
                        !                   + EMIS_BL
                        !                   * ( DELTPRES / TOTPRES )
                        !                   * AREA_CM2 / VOLUME(JLOOP)
                        ! adj code:
                        EMIS_BL_ADJ =        EMIS_BL_ADJ
     &                                     + REMIS_ADJ(JLOOP,IDEOX)
     &                                     * ( DELTPRES / TOTPRES )
     &                                     * AREA_CM2 / VOLUME(JLOOP)
                     ENDDO

                     ! fwd code:
                     !IF (LDRYD) SHIPO3DEP(I,J) = 0d0
                     ! adj code:
                     SHIPO3DEP_ADJ(I,J) = 0d0


                  ELSE

                     ! No change in REMIS

                     ! Recalculate TOTO3

                     ! Initialize
                     TOTO3 = 0d0

                     ! Loop over the boundary layer
                     ! If LNLPBL, then TOP=1
                     DO L = 1, TOP
                        JLOOP   = JLOP(I,J,L)

                        IF ( JLOOP == 0 ) CYCLE

                        ! Cumulative O3 in the PBL, molec(O3)/cm2
                        TOTO3 = TOTO3 + ( BXHEIGHT(I,J,L) *
     &                          1d2 ) * CSPEC_ORIG(JLOOP,IDO3)

                     ENDDO


                     ! fwd code:
                     !IF (LDRYD) SHIPO3DEP(I,J) = ABS( EMIS_BL / TOTO3 )
                     ! adj code:
!$OMP CRITICAL
                     IF ( LDRYD .and. ( TOTO3 .ge. 1d0 ) ) THEN

                        IF ( EMIS_BL / TOTO3 > 0 ) THEN
                           EMIS_BL_ADJ = SHIPO3DEP_ADJ(I,J) / TOTO3
                           TOTO3_ADJ   =  - EMIS_BL * SHIPO3DEP_ADJ(I,J)
     &                                 / ( TOTO3**2 )
                        ELSE
                           EMIS_BL_ADJ =  - SHIPO3DEP_ADJ(I,J) / TOTO3
                           TOTO3_ADJ   =  EMIS_BL * SHIPO3DEP_ADJ(I,J)
     &                                 / ( TOTO3**2 )
                        ENDIF

                        ! ensure that TOTO3_ADJ isn't inf or NaN
                        !IF ( TOTO3 < 1d-150 ) TOTO3_ADJ = 0d0

                        SHIPO3DEP_ADJ(I,J) = 0d0

                     ELSEIF (LDRYD) THEN

                        EMIS_BL_ADJ = 0d0
                        TOTO3_ADJ   = 0d0
                        SHIPO3DEP_ADJ(I,J) = 0d0

                     ENDIF
!$OMP END CRITICAL

                     DO L = 1, TOP
                        JLOOP   = JLOP(I,J,L)

                        IF ( JLOOP == 0 ) CYCLE

                        ! fwd:
                        !TOTO3 = TOTO3 + ( BXHEIGHT(I,J,L) *
                        !        1d2 ) * CSPEC(JLOOP,IDO3)
                        ! adj:
                        CSPEC_ADJ(JLOOP,IDO3) = CSPEC_ADJ(JLOOP,IDO3)
     &                                        + TOTO3_ADJ
     &                                        * BXHEIGHT(I,J,L) * 1d2

                     ENDDO

                     ! fwd:
                     !TOTO3 = 0d0
                     ! adj:
                     TOTO3_ADJ = 0d0

                  ENDIF

                  ! fwd code:
                  !EMIS_BL = SHIP * (1D0 - FRACTION_NOX) * INT_OPE
                  ! adj code: (SHIP_ADJ not used yet, but go ahead and include it)
                  SHIP_ADJ = EMIS_BL_ADJ
     &                     * (1D0 - FRACTION_NOX) * INT_OPE
                  INT_OPE_ADJ  = EMIS_BL_ADJ
     &                     * SHIP * (1D0 - FRACTION_NOX)
                  FRACTION_NOX_ADJ = - EMIS_BL_ADJ
     &                     * SHIP * INT_OPE
                  EMIS_BL_ADJ = 0d0

                  !-----------------
                  ! Ship HNO3
                  !-----------------

                  ! Loop over the boundary layer
                  DO L = 1, TOP
                     JLOOP   = JLOP(I,J,L)

                     IF ( JLOOP == 0 ) CYCLE

                     ! Thickness of level L [mb]
                     DELTPRES = GET_PEDGE(I,J,L) - GET_PEDGE(I,J,L+1)


                     ! fwd:
                     !REMIS(JLOOP,IDEHNO3) = REMIS(JLOOP,IDEHNO3)
                     !                     + ( EMIS_BL * AREA_CM2 )
                     !                     / VOLUME(JLOOP)
                     ! adj:
                     EMIS_BL_ADJ = REMIS_ADJ(JLOOP,IDEHNO3)
     &                           * AREA_CM2 / VOLUME(JLOOP)

                     ! fwd:
                     !EMIS_BL = SHIP * ( 1D0      - FRACTION_NOX ) *
                     !                 ( DELTPRES / TOTPRES      )
                     ! adj:
                     SHIP_ADJ = SHIP_ADJ + EMIS_BL_ADJ
     &                        * ( 1D0 - FRACTION_NOX )
     &                        * ( DELTPRES / TOTPRES )
                     FRACTION_NOX_ADJ = FRACTION_NOX_ADJ
     &                        - SHIP * EMIS_BL_ADJ
     &                        * ( DELTPRES / TOTPRES )

                  ENDDO


                  !-----------------
                  ! Ship NO
                  !-----------------
                  ! Loop over the boundary layer
                  DO L = 1, TOP
                     JLOOP   = JLOP(I,J,L)
                     EMIS_BL = 0d0

                     IF ( JLOOP == 0 ) CYCLE

                     ! Thickness of level L [mb]
                     DELTPRES = GET_PEDGE(I,J,L) - GET_PEDGE(I,J,L+1)


                     ! fwd:
                     !REMIS(JLOOP,IDENOX) = REMIS(JLOOP,IDENOX)
                     !                    + ( EMIS_BL * AREA_CM2 )
                     !                    / VOLUME(JLOOP)
                     ! adj:
                     EMIS_BL_ADJ = REMIS_ADJ(JLOOP,IDENOX)
     &                           * AREA_CM2 / VOLUME(JLOOP)
                     ! fwd:
                     !EMIS_BL = ( SHIP     * FRACTION_NOX )
                     !        * ( DELTPRES / TOTPRES      )
                     ! adj:
                     SHIP_ADJ = SHIP_ADJ + EMIS_BL_ADJ
     &                        * ( DELTPRES / TOTPRES ) * FRACTION_NOX
                     FRACTION_NOX_ADJ = FRACTION_NOX_ADJ + EMIS_BL_ADJ
     &                        * ( DELTPRES / TOTPRES ) * SHIP

                  ENDDO

                  JLOOP   = JLOP(I,J,1)

                  ! fwd code:
                  !CALL INTERPOLATE_LUT2( I, J,
                  !     CSPEC(JLOOP,IDO3),  CSPEC(JLOOP,IDNO),
                  !     CSPEC(JLOOP,IDNO2), AIRDENS(JLOOP),
                  !     JO1D,               JNO2,
                  !     FRACTION_NOx,       INT_OPE)
                  ! adj code:
                  CALL INTERPOLATE_LUT2_ADJ(  I,       J,
     &                 CSPEC_ORIG(JLOOP,IDO3),  O3_ADJ,
     &                 CSPEC_ORIG(JLOOP,IDNO),  NO_ADJ,
     &                 CSPEC_ORIG(JLOOP,IDNO2), NO2_ADJ, AIRDENS(JLOOP),
     &                    JO1D,               JNO2,    FRACTION_NOX,
     &                    FRACTION_NOX_ADJ,   INT_OPE, INT_OPE_ADJ )


                  CSPEC_ADJ(JLOOP,IDO3) = CSPEC_ADJ(JLOOP,IDO3) + O3_ADJ
                  CSPEC_ADJ(JLOOP,IDNO) = CSPEC_ADJ(JLOOP,IDNO) + NO_ADJ
                  CSPEC_ADJ(JLOOP,IDNO2)= CSPEC_ADJ(JLOOP,IDNO2)+NO2_ADJ

               ENDIF

               !========================================================
               ! Adjoint of Aircraft and Lightning NOx [molec/cm3/s]
               !========================================================

               ! bdf - variable tropopause is a tropospheric box
               IF ( LVARTROP ) THEN
                  LTROP = GET_TPAUSE_LEVEL( I, J )
               ELSE
                  LTROP = GET_TPAUSE_LEVEL( I, J ) - 1
               ENDIF


               DO L = 1, LTROP
                  JLOOP   = JLOP(I,J,L)

                  IF ( JLOOP /= 0 ) THEN


                     !-----------------
                     ! Aircraft NOx
                     !-----------------
                     IF ( IS_AC_NOx ) THEN

                        ! fwd:
                        !REMIS(JLOOP,N) = REMIS(JLOOP,N) + EMIS_BL
                        ! adj:
                        EMIS_BL_ADJ = REMIS_ADJ(JLOOP,N)

                        ! fwd:
                        !EMIS_BL        = EMIS_AC_NOx(I,J,L) / COEF1
                        ! adj:
                        IF ( IDADJ_ENOX_ac > 0 ) THEN
                           EMS_SF_ADJ(I,J,M,IDADJ_ENOX_ac)
     &                        = EMS_SF_ADJ(I,J,M,IDADJ_ENOX_ac)
     &                        + EMIS_AC_NOx(I,J,L) / COEF1
     &                        * EMIS_BL_ADJ
                        ENDIF

                     ENDIF

                     !-----------------
                     ! Lightning NOx
                     !-----------------
                     IF ( IS_LI_NOx ) THEN

                        ! fwd code:
                        !REMIS(JLOOP,N) = REMIS(JLOOP,N) + EMIS_BL
                        ! adj code:
                        EMIS_BL_ADJ = REMIS_ADJ(JLOOP,N)

                        ! fwd code:
                        !EMIS_BL        = EMIS_LI_NOx(I,J,L) / COEF1
                        ! adj code:
                        IF ( IDADJ_ENOX_li > 0 ) THEN
                           EMS_SF_ADJ(I,J,M,IDADJ_ENOX_li)
     &                        = EMS_SF_ADJ(I,J,M,IDADJ_ENOX_li)
     &                        + EMIS_LI_NOx(I,J,L) / COEF1
     &                        * EMIS_BL_ADJ
                        ENDIF

                     ENDIF

                  ENDIF

                  ! fwd code:
                  !EMIS_BL = 0d0
                  ! adj code:
                  EMIS_BL_ADJ = 0d0

               ENDDO


               !========================================================
               ! Soil Nox emissions [molec/cm3/s]
               ! Distribute emissions thru the entire boundary layer
               !========================================================
               DO L = 1, TOP
                  JLOOP   = JLOP(I,J,L)

                  IF ( JLOOP /= 0 ) THEN

                     ! Thickness of level L [mb]
                     DELTPRES = GET_PEDGE(I,J,L) - GET_PEDGE(I,J,L+1)

                     ! fwd code:
                     !REMIS(JLOOP,N) = REMIS(JLOOP,N) +
                     !                 ( EMIS_BL / VOLUME(JLOOP) )
                     ! adj code:
                     EMIS_BL_ADJ = REMIS_ADJ(JLOOP,N) / VOLUME(JLOOP)


                     ! fwd code:
                     !EMIS_BL    = ( GEMISNOX2(I,J)                    *
                     !               VOLUME( JLOP(I,J,1) ) / COEF1   ) *
                     !               ( DELTPRES            / TOTPRES )
                     ! adj code:
                     IF ( IDADJ_ENOX_so > 0 ) THEN
                        EMS_SF_ADJ(I,J,M,IDADJ_ENOX_so)
     &                                 = EMS_SF_ADJ(I,J,M,IDADJ_ENOX_so)
     &                                 + GEMISNOX2(I,J)
     &                                 * VOLUME( JLOP(I,J,1) ) / COEF1
     &                                 * DELTPRES / TOTPRES
     &                                 * EMIS_BL_ADJ
                     ENDIF

                  ENDIF

                  ! fwd code:
                  !EMIS_BL = 0d0
                  ! adj code:
                  EMIS_BL_ADJ = 0d0

               ENDDO

               !========================================================
               ! Adjoint of Anthropogenic NOx emissions [molec/box/s]
               !========================================================

               NOXTOT = 0d0
               DO L = 1, NOXEXTENT
                  NOXTOT = NOXTOT + EMISRRN(I,J,L)
               ENDDO

               ! Loop over the boundary layer
               DO L = 1, TOP
                  JLOOP   = JLOP(I,J,L)

                  IF ( JLOOP /= 0 ) THEN

                     ! Thickness of level L [mb]
                     DELTPRES = GET_PEDGE(I,J,L) - GET_PEDGE(I,J,L+1)

                     ! fwd code:
                     !REMIS(JLOOP,N) = EMIS_BL / VOLUME(JLOOP)
                     ! adj code:
                     EMIS_BL_ADJ = REMIS_ADJ(JLOOP,N) / VOLUME(JLOOP)

                     ! fwd code:
                     !EMIS_BL        = ( NOXTOT   / COEF1   ) *
                     !                 ( DELTPRES / TOTPRES )
                     ! adj code:
                     IF ( NADJ_EANTHRO(NN) > 0 ) THEN
                        EMS_SF_ADJ(I,J,M,NADJ_EANTHRO(NN))
     &                             = EMS_SF_ADJ(I,J,M,NADJ_EANTHRO(NN))
     &                             + ( NOXTOT   / COEF1   ) *
     &                               ( DELTPRES / TOTPRES ) *
     &                               EMIS_BL_ADJ
                     ENDIF
                  ENDIF

                  ! fwd code:
                  !EMIS_BL = 0d0
                  ! adj code:
                  EMIS_BL_ADJ = 0d0

               ENDDO

               ! fwd code:
               !NOXTOT = 0d0
               !DO L = 1, NOXEXTENT
               !   NOXTOT = NOXTOT + EMISRRN(I,J,L)
               !ENDDO
               ! adj code: could use this to distinguish between surface and stack emissions
               !DO L = NOXEXTENT, 1, -1
               !   EMISRRN_ADJ(I,J,L) = NOXTOT_ADJ
               !ENDDO
               !! Reset adjoint
               !NOXTOT_ADJ = 0d0

            ENDIF

         ENDDO  ! I
         ENDDO  ! J

!         ! fwd code:
!         !DO JLOOP = 1, NTTLOOP
!         !   REMIS(JLOOP,N) = 0d0
!         !ENDDO
!         ! adj code: take out of N loop
!         DO JLOOP = 1, NTTLOOP
!            REMIS_ADJ(JLOOP,N) = 0d0
!         ENDDO

      ENDDO     ! N
!$OMP END PARALLEL DO

      REMIS_ADJ(:,:) = 0d0

      ! Return to calling program
      END SUBROUTINE SETEMIS_ADJ
