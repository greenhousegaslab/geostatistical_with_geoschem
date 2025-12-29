!$Id: input_adj_mod.f,v 1.21 2012/08/10 22:08:22 nicolas Exp $
      MODULE INPUT_ADJ_MOD
!
!******************************************************************************
!  Module INPUT_ADJ_MOD reads the GEOS-Chem ADJOINT input file (input.gcadj)
!  at the start of the inverse run and passes the information to several other
!  GEOS-Chem F90 modules. It complements input.geos with adjoint specific flags
!  and settings. Most of the code follows the convention from input_mod.f
!  (adj_group, 6/6/09)
!
!  Module Variables:
!  ============================================================================
!  (1 ) VERBOSE   (LOGICAL )  : Turns on echo-back of lines read from disk.
!  (2 ) FIRSTCOL  (INTEGER )  : First column of the input file (default=26)
!  (3 ) MAXDIM    (INTEGER )  : Maximum number of substrings to read in
!  (9 ) FILENAME  (CHAR*255)  : GEOS-CHEM adjoint input file name
!  (10) TOPTITLE  (CHAR*255)  : Top line of input file
!
!  Module Routines:
!  ============================================================================
!  (1 ) READ_INPUT_ADJ_FILE   : Driver routine for reading GEOS-CHEM input file
!  (2 ) READ_ONE_LINE         : Reads one line at a time
!  (3 ) SPLIT_ONE_LINE        : Splits one line into substrings (by spaces)
!  (4 ) READ_ADJ_SIMULATION_MENU  : Reads the GEOS-Chem adjoint simulation menu
!  (5 ) READ_FWD_MODEL_MENU   : Reads forward model options
!  (6 ) READ_ADJ_OPTIONS_MENU : Reads adjoint model options
!  (7 ) READ_ADJ_DIRECTORIES_MENU  : Reads the GEOS-Chem adj. directories
!  (8 ) READ_CONTROL_VARS_MENU: Reads what are control variables
!  (9 ) READ_OBSERVATION_MENU : Reads vars related to observations
!  (10) READ_FD_ MENU         : Reads finite difference test variables
!  (11) READ_ADJ_DIAGNOSTICS_MENU : Reads the GEOS-Chem adj. diagnostic menu
!  (12) VALIDATE_DIRECTORIES  : Makes sure all given directories are valid
!  (13) ARE_FLAGS_VALID       : Makes sure all flags are valid/not conflicting
!  (14) CHECK_DIRECTORY       : Checks a single directory for errors
!  (15) CLEAN_FILE_DIRS       : Clean out directories
!  (16) INIT_DEP_MAPS         : Make mapping arrays for dep adjoint forcing
!  (17) INIT_INPUT_ADJ        : Initializes directory & logical variables
!
!  GEOS-CHEM modules referenced by "input_adj_mod.f"
!  ============================================================================
!  (1 ) directory_adj_mod.f   : Module w/ GC adjoint directories
!  (2 ) error_mod.f           : Module w/ I/O error and NaN check routines
!  (3 ) file_mod.f            : Module w/ file unit numbers and error checks
!  (4 ) grid_mod.f            : Module w/ horizontal grid information
!  (5 ) logical_adj_mod.f     : Module w/ GC adjoint logical switches
!  (6 ) adj_arrays_mod.f      : Module w/ adj. arrays.
!  NOTES:
!  (1 ) Add LPOP_UGM3 (sev, dkh, 02/13/12, adj32_024)
!  (2 ) Add LINVH_BFGS (nab, 25/03/12 )
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables
      ! and routines from being seen outside "input_adj_mod.f"
      !=================================================================

      ! Make everything PRIVATE ...
      PRIVATE

      ! ... except these routines
      PUBLIC :: READ_INPUT_ADJ_FILE
      PUBLIC :: INIT_DEP_MAPS

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================
      LOGICAL            :: VERBOSE  = .FALSE.
      INTEGER, PARAMETER :: FIRSTCOL = 33
      INTEGER, PARAMETER :: MAXDIM   = 255
      INTEGER            :: CT1, CT2, CT3
      CHARACTER(LEN=255) :: FILENAME = 'input.gcadj'
      CHARACTER(LEN=255) :: TOPTITLE

      ! For RRATE list
      LOGICAL            :: READ_STR_ID = .FALSE.
      LOGICAL            :: READ_RXN_ID = .FALSE.

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE READ_INPUT_ADJ_FILE
!
!******************************************************************************
!  Subroutine READ_INPUT_ADJ_FILE is the driver program for reading the
!  GEOS_CHEM adjoint input file "input.gcadj" from disk. (adj_group, 6/07/09)
!
!  NOTES:
!  (1 ) Now call DO_GAMAP (dkh, 02/09/10)
!  (2 ) Now call INIT_TRACERID_ADJ (dkh, 03/30/10)
!******************************************************************************
!
      ! References to F90 modules
      USE ADJ_ARRAYS_MOD, ONLY : INIT_TRACERID_ADJ
      USE CHARPAK_MOD,    ONLY : STRREPL
      USE FILE_MOD,       ONLY : IU_GEOS, IOERROR
      USE INPUT_MOD,      ONLY : TRACERINFO, DIAGINFO
      USE GAMAP_MOD,      ONLY : DO_GAMAP

      ! Local variables
      LOGICAL                :: EOF
      INTEGER                :: IOS
      CHARACTER(LEN=1)       :: TAB   = ACHAR(9)
      CHARACTER(LEN=1)       :: SPACE = ' '
      CHARACTER(LEN=255)     :: LINE

      !=================================================================
      ! READ_INPUT_ADJ_FILE begins here!
      !=================================================================

      ! Echo output
      WRITE( 6, '(a  )' ) REPEAT( '=', 79 )
      WRITE( 6, '(a,/)' )'G E O S - C H E M   A D J O I N T   I N P U T'
      WRITE( 6, 100   ) TRIM( FILENAME )
 100  FORMAT( 'READ_INPUT_ADJ_FILE: Reading ', a )

      ! Initialize directory & logical variables
      CALL INIT_INPUT_ADJ

      ! Initialize adjoint tracer ID's to zero
      CALL INIT_TRACERID_ADJ

      ! Open file
      OPEN( IU_GEOS, FILE=TRIM( FILENAME ), STATUS='OLD', IOSTAT=IOS )
      IF ( IOS /= 0 ) CALL IOERROR( IOS,IU_GEOS,'read_input_adj_file:1')

      ! Read TOPTITLE for binary punch file
      TOPTITLE = READ_ONE_LINE( EOF  )
      IF ( EOF ) RETURN

      ! Loop until EOF
      DO

         ! Read a line from the file, exit if EOF
         LINE = READ_ONE_LINE( EOF )
         IF ( EOF ) EXIT

         ! Replace tab characters in LINE (if any) w/ spaces
         CALL STRREPL( LINE, TAB, SPACE )

         !=============================================================
         ! Call individual subroutines to read sections of the file
         !=============================================================
         IF      ( INDEX( LINE, 'ADJOINT SIMULATION MENU'  ) > 0 ) THEN
            CALL READ_ADJ_SIMULATION_MENU

         ELSE IF ( INDEX( LINE, 'FORWARD MODEL OPTIONS'    ) > 0 ) THEN
            CALL READ_FWD_MODEL_MENU

         ELSE IF ( INDEX( LINE, 'ADJOINT MODEL OPTIONS'    ) > 0 ) THEN
            CALL READ_ADJ_OPTIONS_MENU

         ELSE IF ( INDEX( LINE, 'DIRECTORIES'              ) > 0 ) THEN
            CALL READ_ADJ_DIRECTORIES_MENU

         ELSE IF ( INDEX( LINE, 'CONTROL VARIABLE MENU'   ) > 0 ) THEN
            CALL READ_CONTROL_VARS_MENU

         !mkeller: weak constraint menu
         ELSE IF ( INDEX( LINE, 'WEAK CONSTRAINT MENU'     ) > 0 ) THEN
            CALL READ_WEAK_CONSTRAINT_MENU

         ELSE IF ( INDEX( LINE, 'OBSERVATION MENU'         ) > 0 ) THEN
            CALL READ_OBSERVATION_MENU

         ELSE IF ( INDEX( LINE, 'FINITE DIFFERENCE MENU'   ) > 0 ) THEN
            CALL READ_FD_MENU

         ELSE IF ( INDEX( LINE, 'DIAGNOSTICS MENU'   ) > 0 ) THEN
            CALL READ_ADJ_DIAGNOSTICS_MENU

         ELSE IF ( INDEX( LINE, 'CRITICAL LOAD MENU' ) > 0 ) THEN
            CALL READ_ADJ_CRITICAL_LOAD_MENU

         ELSE IF ( INDEX( LINE, 'END OF FILE'      ) > 0 ) THEN
           EXIT

         ENDIF
      ENDDO

      ! Close input file
      CLOSE( IU_GEOS )

      !=================================================================
      ! Further error-checking and initialization
      !=================================================================

      ! Make sure all directories are valid
      CALL VALIDATE_DIRECTORIES

      ! Clean out file directories (rm *.chk.* , *.adj.* , *.sf.* and
      !  *.gdt.* files )
      CALL CLEAN_FILE_DIRS

      ! Are all the flags a valid combination?
      CALL ARE_FLAGS_VALID

      ! Echo output
      WRITE( 6, '(a)' ) REPEAT( '=', 79 )

      ! Now call this routine here so that adjoint names have been
      ! defined. (dkh, 02/09/10)
      CALL DO_GAMAP( DIAGINFO, TRACERINFO )

      ! Return to calling program
      END SUBROUTINE READ_INPUT_ADJ_FILE

!------------------------------------------------------------------------------

      FUNCTION READ_ONE_LINE( EOF, LOCATION ) RESULT( LINE )
!
!******************************************************************************
!  Subroutine READ_ONE_LINE reads a line from the input file.  If the global
!  variable VERBOSE is set, the line will be printed to stdout.  READ_ONE_LINE
!  can trap an unexpected EOF if LOCATION is passed.  Otherwise, it will pass
!  a logical flag back to the calling routine, where the error trapping will
!  be done. (bmy, 7/20/04)
!
!  Arguments as Output:
!  ===========================================================================
!  (1 ) EOF      (CHARACTER) : Logical flag denoting EOF condition
!  (2 ) LOCATION (CHARACTER) : Name of calling routine; traps premature EOF
!
!  Function value:
!  ===========================================================================
!  (1 ) LINE     (CHARACTER) : A line of text as read from the file
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE FILE_MOD, ONLY : IU_GEOS, IOERROR
      USE FILE_MOD, ONLY : IU_RXN !(hml, 04/03/13)
      USE FILE_MOD, ONLY : IU_STR !(hml, 05/22/13)

      ! Arguments
      LOGICAL,          INTENT(OUT)          :: EOF
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: LOCATION

      ! Local variables
      INTEGER                                :: IOS
      CHARACTER(LEN=255)                     :: LINE, MSG

      !=================================================================
      ! READ_ONE_LINE begins here!
      !=================================================================

      ! Initialize
      EOF = .FALSE.

      ! Read a line from the file (hml, 05/22/13)
      IF ( READ_STR_ID ) READ ( IU_STR, '(a)', IOSTAT=IOS ) LINE
      IF ( READ_RXN_ID ) READ ( IU_RXN, '(a)', IOSTAT=IOS ) LINE
      IF ( .NOT. READ_STR_ID .AND. .NOT. READ_RXN_ID ) THEN
                         READ ( IU_GEOS,'(a)', IOSTAT=IOS ) LINE
      ENDIF

      ! IO Status < 0: EOF condition
      IF ( IOS < 0 ) THEN
         EOF = .TRUE.

         ! Trap unexpected EOF -- stop w/ error msg if LOCATION is passed
         ! Otherwise, return EOF to the calling program
         IF ( PRESENT( LOCATION ) ) THEN
            MSG = 'READ_ONE_LINE: error at: ' // TRIM( LOCATION )
            WRITE( 6, '(a)' ) MSG
            WRITE( 6, '(a)' ) 'Unexpected end of file encountered!'
            WRITE( 6, '(a)' ) 'STOP in READ_ONE_LINE (input_mod.f)'
            WRITE( 6, '(a)' ) REPEAT( '=', 79 )
            STOP
         ELSE
            RETURN
         ENDIF
      ENDIF

      ! IO Status > 0: true I/O error condition (hml, 05/22/13)
      IF ( READ_STR_ID ) THEN
         IF ( IOS > 0 ) CALL IOERROR( IOS, IU_STR, 'read_one_line:1-a' )
      ENDIF
      IF ( READ_RXN_ID ) THEN
         IF ( IOS > 0 ) CALL IOERROR( IOS, IU_RXN, 'read_one_line:1-b' )
      ENDIF
      IF ( .NOT. READ_STR_ID .AND. .NOT. READ_RXN_ID ) THEN
         IF ( IOS > 0 ) CALL IOERROR( IOS, IU_GEOS,'read_one_line:1-c' )
      ENDIF

      ! Print the line (if necessary)
      IF ( VERBOSE ) WRITE( 6, '(a)' ) TRIM( LINE )

      ! Return to calling program
      END FUNCTION READ_ONE_LINE

!------------------------------------------------------------------------------

      SUBROUTINE SPLIT_ONE_LINE( SUBSTRS, N_SUBSTRS, N_EXP, LOCATION )
!
!******************************************************************************
!  Subroutine SPLIT_ONE_LINE reads a line from the input file (via routine
!  READ_ONE_LINE), and separates it into substrings. (bmy, 7/20/04)
!
!  SPLIT_ONE_LINE also checks to see if the number of substrings found is
!  equal to the number of substrings that we expected to find.  However, if
!  you don't know a-priori how many substrings to expect a-priori,
!  you can skip the error check.
!
!  Arguments as Input:
!  ===========================================================================
!  (3 ) N_EXP     (INTEGER  ) : Number of substrings we expect to find
!                               (N_EXP < 0 will skip the error check!)
!  (4 ) LOCATION  (CHARACTER) : Name of routine that called SPLIT_ONE_LINE
!
!  Arguments as Output:
!  ===========================================================================
!  (1 ) SUBSTRS   (CHARACTER) : Array of substrings (separated by " ")
!  (2 ) N_SUBSTRS (INTEGER  ) : Number of substrings actually found
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE CHARPAK_MOD, ONLY: STRSPLIT

      ! Arguments
      CHARACTER(LEN=255), INTENT(OUT) :: SUBSTRS(MAXDIM)
      INTEGER,            INTENT(OUT) :: N_SUBSTRS
      INTEGER,            INTENT(IN)  :: N_EXP
      CHARACTER(LEN=*),   INTENT(IN)  :: LOCATION

      ! Local varaibles
      LOGICAL                         :: EOF
      CHARACTER(LEN=255)              :: LINE, MSG

      !=================================================================
      ! SPLIT_ONE_LINE begins here!
      !=================================================================

      ! Create error msg
      MSG = 'SPLIT_ONE_LINE: error at ' // TRIM( LOCATION )

      !=================================================================
      ! Read a line from disk
      !=================================================================
      LINE = READ_ONE_LINE( EOF )

      ! STOP on End-of-File w/ error msg
      IF ( EOF ) THEN
         WRITE( 6, '(a)' ) TRIM( MSG )
         WRITE( 6, '(a)' ) 'End of file encountered!'
         WRITE( 6, '(a)' ) 'STOP in SPLIT_ONE_LINE (input_mod.f)!'
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         STOP
      ENDIF

      !=================================================================
      ! Split the lines between spaces -- start at column FIRSTCOL
      !=================================================================
      CALL STRSPLIT( LINE(FIRSTCOL:), ' ', SUBSTRS, N_SUBSTRS )

      ! Sometimes we don't know how many substrings to expect,
      ! if N_EXP is greater than MAXDIM, then skip the error check
      IF ( N_EXP < 0 ) RETURN

      ! Stop if we found the wrong
      IF ( N_EXP /= N_SUBSTRS ) THEN
         WRITE( 6, '(a)' ) TRIM( MSG )
         WRITE( 6, 100   ) N_EXP, N_SUBSTRS
         WRITE( 6, '(a)' ) 'STOP in SPLIT_ONE_LINE (input_mod.f)!'
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         STOP
 100     FORMAT( 'Expected ',i2, ' substrs but found ',i3 )
      ENDIF

      ! Return to calling program
      END SUBROUTINE SPLIT_ONE_LINE

!------------------------------------------------------------------------------

!MK-WEAK_CONSTRAINT:
      SUBROUTINE READ_WEAK_CONSTRAINT_MENU
!
!******************************************************************************
!  Subroutine READ_WEAK_CONSTRAINT_MENU reads the WEAK CONSTRAINT MENU section of
!  the GEOS-CHEM adjoint input file (adj_group, 6/07/09)
!
!  NOTES:
!  (1 ) first attempt at subroutine (mkeller)
!******************************************************************************
!
      ! References to F90 modules
      USE WEAK_CONSTRAINT_MOD, ONLY : DO_WEAK_CONSTRAINT
      USE WEAK_CONSTRAINT_MOD, ONLY : MIN_LON_U
      USE WEAK_CONSTRAINT_MOD, ONLY : MAX_LON_U
      USE WEAK_CONSTRAINT_MOD, ONLY : MIN_LAT_U
      USE WEAK_CONSTRAINT_MOD, ONLY : MAX_LAT_U
      USE WEAK_CONSTRAINT_MOD, ONLY : MIN_LEV_U_INDEX
      USE WEAK_CONSTRAINT_MOD, ONLY : MAX_LEV_U_INDEX
      USE WEAK_CONSTRAINT_MOD, ONLY : LEN_SUBWINDOW_U
      USE WEAK_CONSTRAINT_MOD, ONLY : N_TRACER_U
      USE TRACER_MOD,          ONLY : N_TRACERS
      USE ERROR_MOD,           ONLY : ERROR_STOP
      ! Local variables
      INTEGER            :: N
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)

      !=================================================================
      ! READ_WEAK_CONSTRAINT_MENU begins here!
      !=================================================================

      ! Check if we are running the weak constraint module at all
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_weak_constraint_menu:1')
      READ( SUBSTRS(1:N), * ) DO_WEAK_CONSTRAINT
      IF (.NOT. DO_WEAK_CONSTRAINT) THEN
         PRINT*, 'NOT RUNNING WEAK CONSTRAINT MODEL!'
         RETURN
      ENDIF

      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_weak_constraint_menu:2')
      READ( SUBSTRS(1:N), * ) N_TRACER_U

      IF ( N_TRACER_U > N_TRACERS ) THEN
        CALL ERROR_STOP( 'WC-Index bigger than total number of tracers',
     & 'read_weak_constraint_menu, input_adj_mod.f')
      ENDIF

      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_weak_constraint_menu:3')
      READ( SUBSTRS(1:N), * ) MIN_LON_U

      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_weak_constraint_menu:4')
      READ( SUBSTRS(1:N), * ) MAX_LON_U

      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_weak_constraint_menu:5')
      READ( SUBSTRS(1:N), * ) MIN_LAT_U

      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_weak_constraint_menu:6')
      READ( SUBSTRS(1:N), * ) MAX_LAT_U

      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_weak_constraint_menu:7')
      READ( SUBSTRS(1:N), * ) MIN_LEV_U_INDEX

      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_weak_constraint_menu:8')
      READ( SUBSTRS(1:N), * ) MAX_LEV_U_INDEX

      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_weak_constraint_menu:9')
      READ( SUBSTRS(1:N), * ) LEN_SUBWINDOW_U

      END SUBROUTINE READ_WEAK_CONSTRAINT_MENU

!------------------------------------------------------------------------------

      SUBROUTINE READ_ADJ_SIMULATION_MENU

!
!******************************************************************************
!  Subroutine READ_ADJ_SIMULATION_MENU reads the SIMULATION MENU section of
!  the GEOS-CHEM adjoint input file (adj_group, 6/07/09)
!
!  NOTES:
!  (1 ) Reordering and updates (dkh, 02/09/11)
!******************************************************************************
!
      ! References to F90 modules
      USE LOGICAL_MOD,       ONLY : LTRAN
      USE LOGICAL_ADJ_MOD,   ONLY : LADJ
      USE LOGICAL_ADJ_MOD,   ONLY : LFD_SPOT
      USE LOGICAL_ADJ_MOD,   ONLY : LFD_GLOB
      USE LOGICAL_ADJ_MOD,   ONLY : LFDTEST
      USE LOGICAL_ADJ_MOD,   ONLY : LSENS
      USE LOGICAL_ADJ_MOD,   ONLY : L4DVAR
      USE LOGICAL_ADJ_MOD,   ONLY : L3DVAR

      ! Local variables
      INTEGER            :: N
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)

      !=================================================================
      ! READ_ADJ_SIMULATION_MENU begins here!
      !=================================================================

      ! Check if we are running the adjoint at all
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_adj_sim_menu:1' )
      READ( SUBSTRS(1:N), * ) LADJ
      IF (.NOT. LADJ) THEN
         PRINT*, 'NOT RUNNING THE ADJOINT!'
         RETURN
      ENDIF

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_adj_sim_menu:2' )

      !! Doing transport adjoint
      !CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_adj_sim_menu:2' )
      !READ( SUBSTRS(1:N), * ) LADJ_TRAN

      ! Doing 4DVAR
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_adj_sim_menu:3' )
      READ( SUBSTRS(1:N), * ) L4DVAR

       ! Doing 3DVAR
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_adj_sim_menu:4' )
      READ( SUBSTRS(1:N), * ) L3DVAR

      ! Doing sensitivity run (no differences in cost function, just tracer conc)
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_adj_sim_menu:5' )
      READ( SUBSTRS(1:N), * ) LSENS

      ! Move to FORWARD MODEL menu (dkh, 02/09/11)
      !! Doing chemistry
      !CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_adj_sim_menu:3' )
      !READ( SUBSTRS(1:N), * ) LADJ_CHEM
      !
      !! Doing aerosol thermodynamics
      !CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_adj_sim_menu:4' )
      !READ( SUBSTRS(1:N), * ) LAERO_THERM

      IF ( .NOT. ( LSENS .OR. L4DVAR .OR. L3DVAR ) ) THEN
        PRINT*, '******************************************'
        PRINT*, 'HAVE TO PICK A SIMULATION, READ THE MANUAL!'
        PRINT*, '******************************************'
        RETURN
      ENDIF

      ! Check to see if its a finite difference calculation
      IF ( LSENS ) THEN

         ! Doing finite difference test in 1 gridbox
         CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_adj_menu:6' )
         READ( SUBSTRS(1:N), * ) LFD_SPOT

         ! Doing finite difference test in all grid boxes, turn transport off
         CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_adj_menu:7' )
         READ( SUBSTRS(1:N), * ) LFD_GLOB

         ! turn of transport for global FD test
         IF ( LFD_GLOB ) LTRAN = .FALSE.

         ! define a more generic LFDTEST flag if either method is true
         IF ( LFD_GLOB .OR. LFD_SPOT ) LFDTEST = .TRUE.

      ENDIF

      ! Move these to other menus (dkh, 02/09/11)
      !!=================================================================
      !! Include a priori term of the cost function (the one without the data)
      !! aka source term
      !! aka background term
      !! aka penalty term
      !CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_adj_sim_menu:10' )
      !READ( SUBSTRS(1:N), * ) LAPSRC
      !
      !! Compute background error covariance
      !CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_adj_sim_menu:11' )
      !READ( SUBSTRS(1:N), * ) LBKCOV
      !
      !! Compute approximation of inverse Hessian
      !CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_adj_sim_menu:12' )
      !READ( SUBSTRS(1:N), * ) LINVH
      !
      !! include LINOZ
      !! NOTE: This flag controls both forward and adjoint execution
      !CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_adj_sim_menu:13' )
      !READ( SUBSTRS(1:N), * ) LLINOZ
      !
      !! Check if we are running the adjoint at all
      !CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_adj_sim _menu:14' )
      !READ( SUBSTRS(1:N), * ) LRXNR

      !=================================================================
      ! Print to screen
      !=================================================================
      WRITE( 6, '(/,a)' ) 'ADJOINT SIMULATION MENU'
      WRITE( 6, '(  a)' ) '---------------'
      WRITE( 6, 100     ) 'Doing adjoint run             : ', LADJ
      !WRITE( 6, 100     ) 'Doing adjoint transport       : ', LADJ_TRAN
      !WRITE( 6, 100     ) 'Doing adjoint chemistry       : ', LADJ_CHEM
      !WRITE( 6, 100     ) 'Doing aerosol thermodynamics  : ',LAERO_THERM
      WRITE( 6, 100     ) 'Doing 4DVAR (inversion)       : ', L4DVAR
      WRITE( 6, 100     ) 'Doing 3DVAR                   : ', L3DVAR
      WRITE( 6, 100     ) 'Doing sensitivity run         : ', LSENS
      !WRITE( 6, 100     ) 'Include source term in J      : ', LAPSRC
      !WRITE( 6, 100     ) 'Compute background error cov  : ', LBKCOV
      !WRITE( 6, 100     ) 'Compute inverse Hessian       : ', LINVH
      !WRITE( 6, 100     ) 'Use LINOZ (fwd and adj)       : ', LLINOZ
      !WRITE( 6, 100     ) 'Include reaction rates LRXNR  : ', LRXNR
      WRITE( 6, 100     ) 'Doing finite diff check (1box): ', LFD_SPOT
      WRITE( 6, 100     ) 'Doing finite diff check (glob): ', LFD_GLOB


      ! Format statements
 100  FORMAT( A, L5             )
 110  FORMAT( A, I5             )

      !=================================================================
      ! Call setup routines from other GEOS-CHEM modules
      !=================================================================

      ! Set counter
      CT1 = CT1 + 1

      ! Return to calling program
      END SUBROUTINE READ_ADJ_SIMULATION_MENU

!------------------------------------------------------------------------------

      SUBROUTINE READ_FWD_MODEL_MENU
!
!******************************************************************************
!  Subroutine READ_FWD_MODEL_MENU reads the FORWARD MODEL OPTIONS section of
!  the GEOS-CHEM adjoint input file (dkh, 02/09/11)
!
!  NOTES:
!
!******************************************************************************
!
      ! References to F90 modules
      USE LOGICAL_MOD,       ONLY : LTRAN
      USE LOGICAL_ADJ_MOD,   ONLY : LADJ_CHEM
      !USE LOGICAL_ADJ_MOD,   ONLY : LLINOZ
      USE LOGICAL_ADJ_MOD,   ONLY : LAERO_THERM
      USE LOGICAL_ADJ_MOD,   ONLY : LISO

      ! Local variables
      INTEGER            :: N
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)

      !=================================================================
      ! READ_FWD_MODEL_MENU begins here!
      !=================================================================

      ! Doing chemistry
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_fwd_menu:1' )
      READ( SUBSTRS(1:N), * ) LADJ_CHEM

      ! Doing aerosol thermodynamics
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_fwd_menu:2' )
      READ( SUBSTRS(1:N), * ) LAERO_THERM

      ! Use ISORROPIAII
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_fwd_menu:3' )
      READ( SUBSTRS(1:N), * ) LISO

      ! Now use new strat_chem_mod (hml, dkh, 02/14/12, adj32_025)
      !! include LINOZ
      !! NOTE: This flag controls both forward and adjoint execution
      !CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_fwd_menu:3' )
      !READ( SUBSTRS(1:N), * ) LLINOZ

      !=================================================================
      ! Print to screen
      !=================================================================
      WRITE( 6, '(/,a)' ) 'FORWARD MODEL OPTIONS'
      WRITE( 6, '(  a)' ) '---------------'
      WRITE( 6, 100     ) 'Doing adjoint chemistry       : ', LADJ_CHEM
      WRITE( 6, 100     ) 'Doing aerosol thermodynamics  : ',LAERO_THERM
      WRITE( 6, 100     ) ' => ISORROPIAII               : ', LISO
      !WRITE( 6, 100     ) 'Use LINOZ (fwd and adj)       : ', LLINOZ

      ! Format statements
 100  FORMAT( A, L5             )
 110  FORMAT( A, I5             )

      !=================================================================
      ! Call setup routines from other GEOS-CHEM modules
      !=================================================================

      ! Set counter
      CT1 = CT1 + 1

      ! Return to calling program
      END SUBROUTINE READ_FWD_MODEL_MENU

!------------------------------------------------------------------------------

      SUBROUTINE READ_ADJ_OPTIONS_MENU
!
!******************************************************************************
!  Subroutine READ_ADJ_OPTIONS_MENU reads the ADJOINT MODEL OPTIONS section of
!  the GEOS-CHEM adjoint input file (adj_group, 6/07/09)
!
!  NOTES:
!  (1 ) Reordering and updates (dkh, 02/09/11)
!  (2 ) Add LINVH_BFGS (nab, 25/03/12 )
!******************************************************************************
!
      ! References to F90 modules
      USE LOGICAL_MOD,       ONLY : LTRAN
      USE LOGICAL_ADJ_MOD,   ONLY : LAPSRC
      USE LOGICAL_ADJ_MOD,   ONLY : LBKCOV
      USE LOGICAL_ADJ_MOD,   ONLY : LINVH, LINVH_BFGS
      USE LOGICAL_ADJ_MOD,   ONLY : LRXNR
      USE LOGICAL_ADJ_MOD,   ONLY : LDEL_CHKPT
      USE LOGICAL_ADJ_MOD,   ONLY : LFILL_ADJ

      ! Local variables
      INTEGER            :: N
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)

      !=================================================================
      ! READ_ADJ_OPTIONS_MENU begins here!
      !=================================================================

      ! Move these to other menus (dkh, 02/09/11)
      !=================================================================
      ! Include a priori term of the cost function (the one without the data)
      ! aka source term
      ! aka background term
      ! aka penalty term
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_adj_opt_menu:10' )
      READ( SUBSTRS(1:N), * ) LAPSRC

      ! Compute background error covariance
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_adj_opt_menu:11' )
      READ( SUBSTRS(1:N), * ) LBKCOV

      ! Compute approximation of inverse Hessian
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_adj_opt_menu:12' )
      READ( SUBSTRS(1:N), * ) LINVH

      ! Compute approximation of L-BFGS inverse Hessian
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_adj_opt_menu:13' )
      READ( SUBSTRS(1:N), * ) LINVH_BFGS

      ! Compute reaction rate constant sensitivities
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_adj_opt_menu:14' )
      READ( SUBSTRS(1:N), * ) LRXNR

      ! Delete checkpt files
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_adj_opt_menu:15' )
      READ( SUBSTRS(1:N), * ) LDEL_CHKPT

      ! Scale up and FILL adj transport
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_adj_opt_menu:15' )
      READ( SUBSTRS(1:N), * ) LFILL_ADJ

      !=================================================================
      ! Print to screen
      !=================================================================
      WRITE( 6, '(/,a)' ) 'ADJOINT MODEL OPTIONS'
      WRITE( 6, '(  a)' ) '---------------'
      WRITE( 6, 100     ) 'Include source term in J       : ', LAPSRC
      WRITE( 6, 100     ) 'Compute background error cov   : ', LBKCOV
      WRITE( 6, 100     ) 'Compute inverse Hessian        : ', LINVH
      WRITE( 6, 100     ) 'Compute L-BFGS inverse Hessian : '
     &                    , LINVH_BFGS
      WRITE( 6, 100     ) 'Include reaction rates LRXNR   : ', LRXNR
      WRITE( 6, 100     ) 'Delete chkpt files LDEL_CHKPT  : ',
     &    LDEL_CHKPT
      WRITE( 6, 100     ) 'Scale up and FILL adj transport: ', LFILL_ADJ


      ! Format statements
 100  FORMAT( A, L5             )
 110  FORMAT( A, I5             )

      !=================================================================
      ! Call setup routines from other GEOS-CHEM modules
      !=================================================================

      ! Set counter
      CT1 = CT1 + 1

      ! Return to calling program
      END SUBROUTINE READ_ADJ_OPTIONS_MENU

!------------------------------------------------------------------------------

      SUBROUTINE READ_ADJ_DIRECTORIES_MENU
!
!******************************************************************************
!  Subroutine READ_ADJ_DIRECTORIES_MENU reads the DIRECTORIES MENU section of
!  the GEOS-CHEM adjoint input file (adj_group, 6/07/09)
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE DIRECTORY_ADJ_MOD, ONLY : OPTDATA_DIR, ADJTMP_DIR, DIAGADJ_DIR

      ! Local variables
      INTEGER            :: N
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)

      !=================================================================
      ! READ_ADJ_DIRECTORIES_MENU begins here!
      !=================================================================

      ! Optimization output data dir
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_directories_menu:1' )
      READ( SUBSTRS(1:N), '(a)' ) OPTDATA_DIR

      ! Optimization temporary directory
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_directories_menu:2' )
      READ( SUBSTRS(1:N), '(a)' ) ADJTMP_DIR

      ! Optimization diagnostic file directory
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_directories_menu:2' )
      READ( SUBSTRS(1:N), '(a)' ) DIAGADJ_DIR

      WRITE( 6, '(/,a)' ) 'DIRECTORIES  MENU'
      WRITE( 6, '(  a)' ) '---------------'
      WRITE( 6, 110     ) 'Optimization output directory          : ',
     &     TRIM( OPTDATA_DIR )
      WRITE( 6, 110     ) 'Temporary adjoint directory            : ',
     &     TRIM( ADJTMP_DIR )
      WRITE( 6, 110     ) 'Diagnostic adjoint directory           : ',
     &     TRIM( DIAGADJ_DIR )

 110  FORMAT( A, A              )

      ! Set counter
      CT1 = CT1 + 1


      END SUBROUTINE READ_ADJ_DIRECTORIES_MENU
!---------------------------------------------------------------------------------------
!
!      SUBROUTINE READ_CONTROL_PARAMS_MENU
!!
!!******************************************************************************
!!  Subroutine READ_CONTROL_PARAMS_MENU reads the CONTROL PARAMETERS MENU section of
!!  the GEOS-CHEM adjoint input file (adj_group, 6/07/09)
!!
!!  NOTES:
!!  (1 ) Add ICS_SF_tmp, EMS_SF_tmp (mak, dkh, 10/01/09)
!!  (2 ) Merge this with CONTROL_VARS_MENU
!!******************************************************************************
!!
!      ! References to F90 modules
!      USE LOGICAL_ADJ_MOD,   ONLY : LADJ_EMS
!      USE LOGICAL_ADJ_MOD,   ONLY : LICS
!      USE ADJ_ARRAYS_MOD,    ONLY : ICS_SF_tmp
!      USE ADJ_ARRAYS_MOD,    ONLY : EMS_SF_tmp
!
!
!      ! Local variables
!      INTEGER                    :: N
!      CHARACTER(LEN=255)         :: SUBSTRS(MAXDIM)
!
!      !=================================================================
!      ! READ_ADJ_SIMULATION_MENU begins here!
!      !=================================================================
!
!      ! Optimizing emissions
!      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_control_params_menu:1' )
!      READ( SUBSTRS(1:N), * ) LADJ_EMS
!
!      ! Optimizing initial conditions
!      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_control_params_menu:3' )
!      READ( SUBSTRS(1:N), * ) LICS
!
!      ! Optimizing initial conditions
!      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_control_params_menu:4' )
!      READ( SUBSTRS(1:N), * ) ICS_SF_tmp
!
!      ! Optimizing initial conditions
!      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_control_params_menu:5' )
!      READ( SUBSTRS(1:N), * ) EMS_SF_tmp
!
!
!      !=================================================================
!      ! Print to screen
!      !=================================================================
!      WRITE( 6, '(/,a)' ) 'CONTROL PARAMETERS MENU'
!      WRITE( 6, '(  a)' ) '---------------'
!      WRITE( 6, 100     ) 'Optimizing emissions          : ', LADJ_EMS
!      WRITE( 6, 100     ) 'Optimizing initial conditions : ', LICS
!      WRITE( 6, 110     ) 'First guess for ICS_SF is     : ', ICS_SF_tmp
!      WRITE( 6, 110     ) 'First guess for EMS_SF is     : ', EMS_SF_tmp
!
!
!      ! Format statements
! 100  FORMAT( A, L5             )
! 110  FORMAT( A, f7.2           )
!
!
!      !=================================================================
!      ! Call setup routines from other GEOS-CHEM modules
!      !=================================================================
!
!      ! Set counter
!      CT1 = CT1 + 1
!
!      ! Return to calling program
!      END SUBROUTINE READ_CONTROL_PARAMS_MENU
!
!!------------------------------------------------------------------------------
      SUBROUTINE READ_CONTROL_VARS_MENU
!
!******************************************************************************
!  Subroutine READ_CONTROL_VARS_MENU reads the CONTROL VARIABLES MENU section of
!  the GEOS-CHEM adjoint input file (adj_group, 6/07/09)
!
!  NOTES:
!  (1 ) Reorder and update (dkh, 02/09/11)
!  (2 ) Add support for strat fluxes LADJ_STRAT (hml, dkh, 02/14/12, adj32_025)
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD,         ONLY : ALLOC_ERR
      USE ERROR_MOD,         ONLY : ERROR_STOP
      USE LOGICAL_ADJ_MOD,   ONLY : LADJ_EMS
      USE LOGICAL_ADJ_MOD,   ONLY : LICS
      USE ADJ_ARRAYS_MOD,    ONLY : NNEMS
      USE ADJ_ARRAYS_MOD,    ONLY : ID_ADEMS
      USE ADJ_ARRAYS_MOD,    ONLY : ADEMS_NAME
      USE ADJ_ARRAYS_MOD,    ONLY : TRACERID_ADJ
      USE ADJ_ARRAYS_MOD,    ONLY : INIT_ADJ_EMS
      USE ADJ_ARRAYS_MOD,    ONLY : OPT_THIS_EMS
      USE ADJ_ARRAYS_MOD,    ONLY : OPT_THIS_TRACER
      USE ADJ_ARRAYS_MOD,    ONLY : REG_PARAM_ICS
      USE ADJ_ARRAYS_MOD,    ONLY : REG_PARAM_EMS
      USE ADJ_ARRAYS_MOD,    ONLY : MMSCL
      USE ADJ_ARRAYS_MOD,    ONLY : EMS_ERROR
      USE ADJ_ARRAYS_MOD,    ONLY : COV_ERROR_LX, COV_ERROR_LY
      USE ADJ_ARRAYS_MOD,    ONLY : ICS_ERROR
      USE ADJ_ARRAYS_MOD,    ONLY : ICS_SF_DEFAULT
      USE ADJ_ARRAYS_MOD,    ONLY : EMS_SF_DEFAULT
      USE TRACER_MOD,        ONLY : N_TRACERS
      USE TRACER_MOD,        ONLY : TRACER_NAME

      ! for strat prod and loss SF (hml, 08/14/11)
      USE LOGICAL_ADJ_MOD,   ONLY : LADJ_STRAT
      USE ADJ_ARRAYS_MOD,    ONLY : NSTPL
      USE ADJ_ARRAYS_MOD,    ONLY : ID_PROD
      USE ADJ_ARRAYS_MOD,    ONLY : ID_LOSS
      USE ADJ_ARRAYS_MOD,    ONLY : PROD_NAME
      USE ADJ_ARRAYS_MOD,    ONLY : LOSS_NAME
      USE ADJ_ARRAYS_MOD,    ONLY : STRPID_ADJ
      USE ADJ_ARRAYS_MOD,    ONLY : STRLID_ADJ
      USE ADJ_ARRAYS_MOD,    ONLY : INIT_ADJ_STRAT
      USE ADJ_ARRAYS_MOD,    ONLY : OPT_THIS_PROD
      USE ADJ_ARRAYS_MOD,    ONLY : OPT_THIS_LOSS
      USE ADJ_ARRAYS_MOD,    ONLY : REG_PARAM_PROD
      USE ADJ_ARRAYS_MOD,    ONLY : REG_PARAM_LOSS
      USE ADJ_ARRAYS_MOD,    ONLY : PROD_ERROR
      USE ADJ_ARRAYS_MOD,    ONLY : LOSS_ERROR
      USE ADJ_ARRAYS_MOD,    ONLY : PROD_SF_DEFAULT
      USE ADJ_ARRAYS_MOD,    ONLY : LOSS_SF_DEFAULT

      ! for reaction rates (tww, 05/08/12)
      USE LOGICAL_ADJ_MOD,   ONLY : LADJ_RRATE
      USE ADJ_ARRAYS_MOD,    ONLY : NRRATES
      USE ADJ_ARRAYS_MOD,    ONLY : ID_RRATES
      USE ADJ_ARRAYS_MOD,    ONLY : RRATES_NAME
      USE ADJ_ARRAYS_MOD,    ONLY : INIT_ADJ_RRATES
      USE ADJ_ARRAYS_MOD,    ONLY : OPT_THIS_RATE
      USE ADJ_ARRAYS_MOD,    ONLY : REG_PARAM_RATE
      USE ADJ_ARRAYS_MOD,    ONLY : RATE_ERROR
      USE ADJ_ARRAYS_MOD,    ONLY : RATE_SF_DEFAULT
      USE GCKPP_ADJ_GLOBAL,  ONLY : NCOEFF_RATE
      ! (hml, 05/22/13)
      USE LOGICAL_ADJ_MOD,   ONLY : FI_STRID
      USE LOGICAL_ADJ_MOD,   ONLY : FI_RXNID
      USE FILE_MOD,          ONLY : IOERROR
      USE FILE_MOD,          ONLY : IU_STR
      USE FILE_MOD,          ONLY : IU_RXN

#     include "define_adj.h"

      ! Local variables
      INTEGER                    :: N, T, NSOPT, TMP, AS
      CHARACTER(LEN=255)         :: SUBSTRS(MAXDIM)
      ! For RRATE list (hml, 04/03/13)
      INTEGER                    :: IOS

      !=================================================================
      ! READ_CONTROL_VARS_MENU begins here!
      !=================================================================

      !=================================================================
      ! Allocate arrays
      !=================================================================
      ! First allocate OPT_THIS_TRACER  to be max species
      ALLOCATE( OPT_THIS_TRACER( N_TRACERS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'OPT_THIS_TRACER' )
      OPT_THIS_TRACER = .FALSE.

      ALLOCATE( REG_PARAM_ICS( N_TRACERS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'REG_PARAM_ICS' )
      REG_PARAM_ICS = 1d0

      ALLOCATE( ICS_ERROR( N_TRACERS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ICS_ERROR' )
      ICS_ERROR = 1d0
#if   defined ( LOG_OPT )
      ICS_ERROR = EXP(1d0)
#endif

      ALLOCATE( ICS_SF_DEFAULT( N_TRACERS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ICS_SF_DEFAULT' )
      ICS_SF_DEFAULT = 1d0

      !=================================================================
      ! Read menu
      !=================================================================

      ! Optimizing initial conditions
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_control_vars_menu:1' )
      READ( SUBSTRS(1:N), * ) LICS

      ! Optimizing emissions
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_control_vars_menu:2' )
      READ( SUBSTRS(1:N), * ) LADJ_EMS

      ! Optimizing strat prod & loss (hml, 08/11/11, adj32_025)
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_control_vars_menu:3a' )
      READ( SUBSTRS(1:N), * ) LADJ_STRAT

      ! Specifying reaction rates (tww, 05/08/12)
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_control_vars_menu:3d' )
      READ( SUBSTRS(1:N), * ) LADJ_RRATE

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_control_vars_menu:3b' )

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_control_vars_menu:3c')

      ! Number of species to optimize
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_control_vars_menu:4' )
      READ( SUBSTRS(1:N), * ) NSOPT

      IF ( LICS .AND. NSOPT .EQ. 0) THEN
         CALL ERROR_STOP( ' LICS is T but NSOPT is 0 ',
     &        ' READ_CONTROL_VARS_MENU, geos_chem_mod.f ' )
      ENDIF

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_control_vars_menu:5' )

      DO T = 1, NSOPT
         ! Split line into substrings
         CALL SPLIT_ONE_LINE( SUBSTRS, N, -1,'read_control_vars_menu:6')

         ! set OPT_THIS_TRACER to true for species we're optimizing
         READ( SUBSTRS(1), * ) TMP
         OPT_THIS_TRACER(TMP) = .TRUE.

         ! now move this to observation menu (dkh, 02/09/11)
         !! observe this species?
         !READ( SUBSTRS(3), *) OBS_THIS_SPECIES(TMP)

         ! Defualt scaling factor for this initial condition
          READ( SUBSTRS(3), *) ICS_SF_DEFAULT(TMP)

         ! REG_PARAM for this species
         READ( SUBSTRS(4), *) REG_PARAM_ICS(TMP)

         ! ICS_ERROR for this emission
         READ( SUBSTRS(5), *) ICS_ERROR(TMP)

      ENDDO

      ! Obsolete -- now we only list tracer that are observed
      ! compute number of observed species
      !NOBS = 0
      !DO T = 1, N_TRACERS
      !   IF ( OBS_THIS_SPECIES(T) ) THEN
      !      NOBS = NOBS + 1
      !   ENDIF
      !ENDDO

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_control_vars_menu:7' )

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_control_vars_menu:7b')

      ! Optimizing emissions
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_control_vars_menu:8' )
      READ( SUBSTRS(1:N), * ) NNEMS

      IF ( .NOT. LADJ_EMS ) NNEMS = 0

      ! If we're optimizing initial conditions, number of tracers is
      !N_TRACERS

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_control_vars_menu:9' )

      IF ( LADJ_EMS .AND. NNEMS .GT. 0) THEN

         CALL INIT_ADJ_EMS

      ELSEIF ( LADJ_EMS .AND. NNEMS .EQ. 0) THEN
         CALL ERROR_STOP( ' LADJ_EMS is T but NNEMS is 0 ',
     &        ' READ_CONTROL_VARS_MENU, geos_chem_mod.f ' )
      ENDIF

      !=================================================================
      ! Read emission ID
      !=================================================================
      IF ( LADJ_EMS ) THEN
         DO T = 1, NNEMS

            ! Split line into substrings
            CALL SPLIT_ONE_LINE( SUBSTRS, N, -1,
     &           'read_control_vars_menu:10')

            ! Save tracer number
            READ( SUBSTRS(1), * ) ID_ADEMS(T)

            ! Save tracer name
            ADEMS_NAME(T)  = TRIM( SUBSTRS(2) )

            ! optimize this emission?Q
            READ( SUBSTRS(3), *) OPT_THIS_EMS(T)

            ! Defualt scaling factor for this emission
            READ( SUBSTRS(4), *) EMS_SF_DEFAULT(T)

            ! REG_PARAM for this emission
            READ( SUBSTRS(5), *) REG_PARAM_EMS(T)

            ! EMS_ERROR for this emission
            READ( SUBSTRS(6), *) EMS_ERROR(T)

            ! CORR_LX for this emission
            READ( SUBSTRS(7), *) COV_ERROR_LX(T)

            ! CORR_LY for this emission
            READ( SUBSTRS(8), *) COV_ERROR_LY(T)

         ENDDO

         ! Number of temporal groups of the control vector,
         ! e.g. monthly optimization in a year-long simulation would have
         ! 12. If in doubt, set to 1
!         CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,
!     &      'read_control_vars_menu:11' )
!         READ( SUBSTRS(1:N), * ) MMSCL

         ! Strat prod and loss (hml, adj32_025)
         ! Separator line
         CALL SPLIT_ONE_LINE( SUBSTRS, N, -1,
     &                      'read_control_vars_menu:12b')

         ! Separator line
         CALL SPLIT_ONE_LINE( SUBSTRS, N, -1,
     &                      'read_control_vars_menu:12c')

         ! Optimizing strat prod & loss
         CALL SPLIT_ONE_LINE( SUBSTRS, N,  1,
     &                      'read_control_vars_menu:13' )
         READ( SUBSTRS(1:N), * ) NSTPL
         IF ( .NOT. LADJ_STRAT ) NSTPL = 0

         ! Read the list from file? (hml, 05/21/13)
         CALL SPLIT_ONE_LINE( SUBSTRS, N,  1,
     &                      'read_control_vars_menu:14' )
         READ( SUBSTRS(1:N), * ) FI_STRID

         ! Separator line
         CALL SPLIT_ONE_LINE( SUBSTRS, N, -1,
     &                      'read_control_vars_menu:15' )

         IF ( LADJ_STRAT .AND. NSTPL .GT. 0) THEN

            CALL INIT_ADJ_STRAT

         ELSEIF ( LADJ_STRAT .AND. NSTPL .EQ. 0) THEN
            CALL ERROR_STOP( ' LADJ_STRAT is T but NSTPL is 0 ',
     &              ' READ_CONTROL_VARS_MENU, geos_chem_mod.f ' )
         ENDIF

         !PRINT *, ' NSTPL = ' , NSTPL

         !=================================================================
         ! Read Stratospheric Tracers ID
         !=================================================================
         IF ( LADJ_STRAT .AND. .NOT. FI_STRID ) THEN ! (hml, 05/21/13)

            DO T = 1, NSTPL

               ! Split line into substrings
               CALL SPLIT_ONE_LINE( SUBSTRS, N, -1,
     &              'read_control_vars_menu:16')

               ! Save tracer number
               READ( SUBSTRS(1), * ) ID_PROD(T)

               ! Save tracer name
               PROD_NAME(T)  = TRIM( SUBSTRS(2) )

               ! optimize this strat prod & loss?
               READ( SUBSTRS(3), *) OPT_THIS_PROD(T)

               ! Defualt prod scaling factor for this strat tracer
               READ( SUBSTRS(4), *) PROD_SF_DEFAULT(T)

               ! REG_PARAM for this strat tracer
               READ( SUBSTRS(5), *) REG_PARAM_PROD(T)

               ! STR_ERROR for this strat tracer
               READ( SUBSTRS(6), *) PROD_ERROR(T)

            ENDDO

            DO T = 1, NSTPL

               ! Split line into substrings
               CALL SPLIT_ONE_LINE( SUBSTRS, N, -1,
     &              'read_control_vars_menu:16-b')

               ! Save tracer number
               READ( SUBSTRS(1), * ) ID_LOSS(T)

               ! Save tracer name
               LOSS_NAME(T)  = TRIM( SUBSTRS(2) )

               ! optimize this strat prod & loss?
               READ( SUBSTRS(3), *) OPT_THIS_LOSS(T)

               ! Defualt loss scaling factor for this strat tracer
               READ( SUBSTRS(4), *) LOSS_SF_DEFAULT(T)

               ! REG_PARAM for this strat tracer
               READ( SUBSTRS(5), *) REG_PARAM_LOSS(T)

               ! STR_ERROR for this strat tracer
               READ( SUBSTRS(6), *) LOSS_ERROR(T)

            ENDDO

         ELSE IF ( LADJ_STRAT .AND. FI_STRID ) THEN

            CALL READ_STRID_FILE

         ENDIF

         ! Separator line
         CALL SPLIT_ONE_LINE( SUBSTRS, N, -1,
     &                        'read_control_vars_menu:15' )

         ! Separator line
         CALL SPLIT_ONE_LINE( SUBSTRS, N, -1,
     &                        'read_control_vars_menu:16' )

         ! Specifying reaction rates
         CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,
     &                        'read_control_vars_menu:16b' )
         READ( SUBSTRS(1:N), * ) NRRATES

         IF ( .NOT. LADJ_RRATE ) NRRATES = 0

         ! Read the list from file?
         CALL SPLIT_ONE_LINE( SUBSTRS, N,  1,
     &                      'read_control_vars_menu:16c' )
         READ( SUBSTRS(1:N), * ) FI_RXNID

         ! Separator line
         CALL SPLIT_ONE_LINE( SUBSTRS, N, -1,
     &                        'read_control_vars_menu:16d' )

         IF ( LADJ_RRATE .AND. NRRATES .GT. 0) THEN

            CALL INIT_ADJ_RRATES

         ELSEIF ( LADJ_RRATE .AND. NRRATES .EQ. 0) THEN

            CALL ERROR_STOP( ' LADJ_RRATE is T but NRRATES is 0 ',
     &           ' READ_CONTROL_VARS_MENU, geos_chem_mod.f ' )

         ENDIF

         IF ( LADJ_RRATE .AND. ( NRRATES .NE. NCOEFF_RATE ) ) THEN

            print*, 'NRRATES =', NRRATES
            print*, 'NCOEFF_RATE =', NCOEFF_RATE
            CALL ERROR_STOP( 'NRRATES not equal NCOEFF_RATE
     &           Check gckpp_adj_Global.f90 ',
     &           ' READ_CONTROL_VARS_MENU, geos_chem_mod.f ' )

         ENDIF



         !=================================================================
         ! Read Reaction Rates ID
         !=================================================================
         IF ( LADJ_RRATE .AND. .NOT. FI_RXNID ) THEN

            ! Added block to read reaction rate entries (tww, 05/08/12)

            DO T = 1, NRRATES

               ! Split line into substrings
               CALL SPLIT_ONE_LINE( SUBSTRS, N, -1,
     &              'read_control_vars_menu:17')

               ! Save tracer number
               READ( SUBSTRS(1), * ) ID_RRATES(T)

               ! Save tracer name
               RRATES_NAME(T)  = TRIM( SUBSTRS(2) )

               ! Optimize this rate?
               READ( SUBSTRS(3), *) OPT_THIS_RATE(T)

               ! Default scaling factor for this rate
               READ( SUBSTRS(4), *) RATE_SF_DEFAULT(T)

               ! REG_PARAM for this rate
               READ( SUBSTRS(5), *) REG_PARAM_RATE(T)

               ! RATE_ERROR for this rate
               READ( SUBSTRS(6), *) RATE_ERROR(T)

            ENDDO

         ELSEIF ( LADJ_RRATE .AND. FI_RXNID ) THEN

            CALL READ_RXNID_FILE

         ENDIF

      ENDIF

      !=================================================================
      ! Print to screen
      !=================================================================
      WRITE( 6, '(/,a)' )
     &       'CONTROL VARIABLE MENU'
      WRITE( 6, '(  a)' ) '---------------'
      WRITE( 6, 130     ) 'Optimizing initial conditions : ', LICS
      WRITE( 6, 130     ) 'Optimizing emissions          : ', LADJ_EMS
      WRITE( 6, 130     ) 'Optimizing strat prod & loss  : ', LADJ_STRAT
      WRITE( 6, 130     ) 'Optimizing rxn rates          : ', LADJ_RRATE
      WRITE( 6, '(  a)' ) '>------------------------------<'

      IF ( LICS ) THEN
         WRITE( 6, '(  a)' )
     &    '  Tracers optimizing SF_DEFAULT REG_PARAM ERROR'
         ! Print info about each tracer
         DO T = 1, N_TRACERS

            IF( OPT_THIS_TRACER(T) ) THEN
               ! Write tracer number, name and it's default scaling factor
               WRITE( 6, 140 )  T, TRACER_NAME(T), ICS_SF_DEFAULT(T),
     &                          REG_PARAM_ICS(T), ICS_ERROR(T)
            ENDIF

         ENDDO

      ELSEIF ( LADJ_EMS ) THEN

         WRITE( 6, '(  a)' )
     &     '  # Emission               Opt  SF   REG   ERR'

         ! Print info about each tracer
         DO T = 1, NNEMS

            ! Write tracer number, name, optimize, default SF, reg param
            ! and error
            WRITE( 6, 120 ) ID_ADEMS(T), ADEMS_NAME(T), OPT_THIS_EMS(T),
     &         EMS_SF_DEFAULT(T), REG_PARAM_EMS(T), EMS_ERROR(T)

         ENDDO

         ! WRITE( 6, 110     ) 'Number of time contrl groups  : ', MMSCL

         ! Strat prod and loss (hml)
         IF ( LADJ_STRAT ) THEN

            WRITE( 6, '(  a)' )
     &        '  # Strat trc              Opt  SF    REG   ERR'

            ! Print info about each prod tracer
            DO T = 1, NSTPL

               ! Write tracer number, name, default SF of prod,
               ! reg param, and error
               WRITE( 6, 150 ) ID_PROD(T), PROD_NAME(T),
     &            OPT_THIS_PROD(T),        PROD_SF_DEFAULT(T),
     &            REG_PARAM_PROD(T),       PROD_ERROR(T)

            ENDDO

            CALL STRPID_ADJ

            ! Print info about each tracer loss
            DO T = 1, NSTPL

               ! Write tracer number, name, default SF of loss,
               ! reg param, and error
               WRITE( 6, 150 ) ID_LOSS(T), LOSS_NAME(T),
     &          OPT_THIS_LOSS(T),          LOSS_SF_DEFAULT(T),
     &          REG_PARAM_LOSS(T),         LOSS_ERROR(T)

            ENDDO

            CALL STRLID_ADJ

         ENDIF

         ! Print info about rxn rates (tww, 05/15/12)
         IF ( LADJ_RRATE ) THEN

            WRITE( 6, '(a)' )
     &             ' RXN ID    NAME       OPT  DEF_SF  REG  ERR'

            ! Print info about each tracer
            DO T = 1, NRRATES

               ! Write tracer number, name
               WRITE( 6, 150 ) ID_RRATES(T), RRATES_NAME(T),
     &                         OPT_THIS_RATE(T), RATE_SF_DEFAULT(T),
     &                         REG_PARAM_RATE(T), RATE_ERROR(T)

            ENDDO

         ENDIF

         !=================================================================
         ! Call setup routines from other F90 modules
         !=================================================================

         CALL TRACERID_ADJ

      ENDIF

      ! Set counter
      CT1 = CT1 + 1

      ! Format statements
 100  FORMAT( I3, 1x, A10, 6x, f5.2, 6x, f5.2 )
 110  FORMAT( A, I5             )
! 120  FORMAT( I3, 1x, A14, 6x, L5, 1x, f5.2, 1x, f5.2, 1x, f5.2 )
 120  FORMAT( I3, 1x, A14, 6x, L5, 1x, f5.2, 1x, f6.2, 1x, f5.2 )
 130  FORMAT( A, L5             )
 140  FORMAT( I3, 1x, A10, 6x, f5.2, 6x, f5.2, 6x f5.2 )
 150  FORMAT( I3, 1x, A14, 5x, L5, 1x, f5.2, 1x, f5.2, 1x, f5.2 )

      ! Return to calling program
      END SUBROUTINE READ_CONTROL_VARS_MENU

!------------------------------------------------------------------------------

      SUBROUTINE READ_OBSERVATION_MENU
!
!******************************************************************************
!  Subroutine READ_OBSERVATION_MENU reads the OBSERVATION OPTIONS MENU section of
!  the GEOS-CHEM adjoint input file (adj_group, 6/07/09)
!
!  NOTES:
!   a) calculate NSPAN using NYMDf, NHMSf (fp)
!******************************************************************************
!
      ! References to F90 modules
      USE ADJ_ARRAYS_MOD,    ONLY : OBS_FREQ
      USE ADJ_ARRAYS_MOD,    ONLY : OBS_THIS_SPECIES
      USE ADJ_ARRAYS_MOD,    ONLY : OBS_THIS_TRACER
      USE ADJ_ARRAYS_MOD,    ONLY : NSPAN
      USE ADJ_ARRAYS_MOD,    ONLY : NOBS
      USE ADJ_ARRAYS_MOD,    ONLY : NOBS_CSPEC
      USE ADJ_ARRAYS_MOD,    ONLY : IDCSPEC_ADJ
      USE ADJ_ARRAYS_MOD,    ONLY : ID2C
      USE ADJ_ARRAYS_MOD,    ONLY : GET_SPEC
      USE ADJ_ARRAYS_MOD,    ONLY : CNAME
      USE ADJ_ARRAYS_MOD,    ONLY : FORCING_MASK_FILE
      USE ADJ_ARRAYS_MOD,    ONLY : FORCING_MASK_VARIABLE
      USE ADJ_ARRAYS_MOD,    ONLY : FORCING_MASK_FILE_NC
      USE ADJ_ARRAYS_MOD,    ONLY : NB_MASK_VAR
      USE ADJ_ARRAYS_MOD,    ONLY : DEP_UNIT
      USE ADJ_ARRAYS_MOD,    ONLY : TRACER_IND
      USE ADJ_ARRAYS_MOD,    ONLY : NFD
      USE ERROR_MOD,         ONLY : ALLOC_ERR
      USE ERROR_MOD,         ONLY : ERROR_STOP
      USE LOGICAL_ADJ_MOD,   ONLY : LMAX_OBS
      USE LOGICAL_ADJ_MOD,   ONLY : LKGBOX
      USE LOGICAL_ADJ_MOD,   ONLY : LUGM3
      USE LOGICAL_ADJ_MOD,   ONLY : LPOP_UGM3
      USE LOGICAL_ADJ_MOD,   ONLY : LSTT_PPB
      USE LOGICAL_ADJ_MOD,   ONLY : LSTT_TROP_PPM
      USE LOGICAL_ADJ_MOD,   ONLY : LCSPEC_PPB
      USE LOGICAL_ADJ_MOD,   ONLY : LCSPEC_OBS
      USE LOGICAL_ADJ_MOD,   ONLY : LSENS
      USE LOGICAL_ADJ_MOD,   ONLY : LFDTEST
      USE LOGICAL_ADJ_MOD,   ONLY : LFD_GLOB
      USE LOGICAL_ADJ_MOD,   ONLY : LADJ_FDEP
      USE LOGICAL_ADJ_MOD,   ONLY : LADJ_DDEP_TRACER
      USE LOGICAL_ADJ_MOD,   ONLY : LADJ_DDEP_CSPEC
      USE LOGICAL_ADJ_MOD,   ONLY : LADJ_WDEP_CV
      USE LOGICAL_ADJ_MOD,   ONLY : LADJ_WDEP_LS
      USE LOGICAL_ADJ_MOD,   ONLY : LKGNHAYR
      USE LOGICAL_ADJ_MOD,   ONLY : LEQHAYR
      USE LOGICAL_ADJ_MOD,   ONLY : LMOLECCM2S
      USE LOGICAL_ADJ_MOD,   ONLY : LKGS
      USE LOGICAL_ADJ_MOD,   ONLY : LFORCE_MASK
      USE LOGICAL_ADJ_MOD,   ONLY : LFORCE_MASK_NC, LFORCE_MASK_BPCH
      USE LOGICAL_ADJ_MOD,   ONLY : LFD_GLOB
      USE TRACER_MOD,        ONLY : N_TRACERS, ITS_A_FULLCHEM_SIM
      USE TIME_MOD,          ONLY : GET_JD, NYMDe, NHMSe, NYMDb, NHMSb
      ! for flux based cost function (hml,06/13/12)
      USE LOGICAL_ADJ_MOD,   ONLY : LFLX_UGM2

#     include "CMN_SIZE"
#     include "comode.h"          ! IGAS, NAMEGAS

      ! Local variables
      INTEGER                    :: N,T,J
      INTEGER                    :: TMP
      INTEGER                    :: NUNIT_COUNT
      INTEGER                    :: AS
      CHARACTER(LEN=255)         :: SUBSTRS(MAXDIM)
      LOGICAL                    :: EOF
      CHARACTER(LEN=255)         :: LINE
      CHARACTER(LEN=15)          :: TNAME(N_TRACERS)
      REAL*8                     :: MASK_AREA
      INTEGER                    :: NHMSf, NYMDf              !fp
      REAL*8                     :: JDF, JDE, JDB             !fp

      !=================================================================
      ! READ_OBSERVATION_MENU begins here!
      !=================================================================

      ! First allocate OBS_THIS_TRACER to be max species
      ALLOCATE( OBS_THIS_TRACER( N_TRACERS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'OBS_THIS_TRACER' )
      OBS_THIS_TRACER = 0

      ! Also allocation the mapping between observed and all tracers
      ALLOCATE( TRACER_IND( N_TRACERS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'TRACER_IND' )
      TRACER_IND = 0

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_obs_menu:1'  )

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_obs_menu:2'  )

      ! Optimization output data dir
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_obs_menu:3'  )
      READ( SUBSTRS(1:N), * ) OBS_FREQ

      ! Maximum number of obs?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_obs_menu:4'  )
      READ( SUBSTRS(1:N), * ) LMAX_OBS

      ! Number of obs evaluations
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 2, 'read_observation_menu:5' )
      READ( SUBSTRS(1:N), * ) NYMDf, NHMSf

      !calculate nspan (fp)
      JDE = GET_JD(NYMDe, NHMSe)
      JDB = GET_JD(NYMDb, NHMSb)

! Need to review this with fp, not sure it gives desired behavior
!      !fp to avoid error with new definition of nspan when LMAX_OBS is false but LFD_GLOB is true
!      IF ( LFD_GLOB ) THEN
!         IF ( .not. LMAX_OBS ) THEN
!            LMAX_OBS = .TRUE.
!            NYMDF = NYMDB
!            NHMSF = NHMSB
!         ENDIF
!      ENDIF

      JDF = GET_JD(NYMDf, NHMSf)

      ! add error catch (fp)
      IF ( JDB .GT. JDF .and. LMAX_OBS) THEN
         CALL ERROR_STOP(
     &        ' You cannot force adjoint beyond simulation start time ',
     &        ' input_adj_mod.f ')
      ENDIF

      ! add error catch (yd)
      IF ( JDF .GT. JDE .and. LMAX_OBS ) THEN
         CALL ERROR_STOP(
     &        ' The forcing time period is outside of the run period',
     &        ' input_adj_mod.f ')
      ENDIF

      NSPAN = NINT( ( JDE - JDF ) * 24D0 * 60D0 /  OBS_FREQ  )

      ! Want to only evalute CF once for FD_GLOB test
      IF ( LFD_GLOB ) THEN
          LMAX_OBS = .TRUE.
          NSPAN    = 1
      ENDIF

      !=================================================================
      ! Cost function options
      !=================================================================
      NUNIT_COUNT = 0

      ! Separator line: COST FUNCTION options for LSENS:---
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_observation_menu:6' )

      ! Cost function STT in kg / box
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_observation_menu:7' )
      READ( SUBSTRS(1:N), * ) LKGBOX
      IF ( LKGBOX ) NUNIT_COUNT = NUNIT_COUNT + 1

      ! Cost function STT in ug / m3
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_observation_menu:8' )
      READ( SUBSTRS(1:N), * ) LUGM3
      IF ( LUGM3  ) NUNIT_COUNT = NUNIT_COUNT + 1

      ! Cost function STT in ppb
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_observation_menu:9' )
      READ( SUBSTRS(1:N), * ) LSTT_PPB
      IF ( LSTT_PPB ) NUNIT_COUNT = NUNIT_COUNT + 1

      ! Cost function STT in free trop in ppm
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_observation_menu:10' )
      READ( SUBSTRS(1:N), * ) LSTT_TROP_PPM
      IF ( LSTT_TROP_PPM ) NUNIT_COUNT = NUNIT_COUNT + 1

      ! Cost function CSPEC in ppb
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_observation_menu:11' )
      READ( SUBSTRS(1:N), * ) LCSPEC_PPB
      IF ( LCSPEC_PPB ) NUNIT_COUNT = NUNIT_COUNT + 1

      ! Cost function STT in population weighted ug / m3 (adj32_024)
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_observation_menu:12' )
      READ( SUBSTRS(1:N), * ) LPOP_UGM3
      IF ( LPOP_UGM3  ) NUNIT_COUNT = NUNIT_COUNT + 1

      ! Cost function STT in flux ug / m2 / hr (hml,06/13/12)
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_observation_menu:13' )
      READ( SUBSTRS(1:N), * ) LFLX_UGM2
      IF ( LFLX_UGM2  ) NUNIT_COUNT = NUNIT_COUNT + 1

      ! deposition based cost function?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_observation_menu:13.a' )
      READ( SUBSTRS(1:N), * ) LADJ_FDEP

      ! tracer dry deposition
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_observation_menu:13.b' )
      READ( SUBSTRS(1:N), * ) LADJ_DDEP_TRACER

      ! species dry deposition
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_observation_menu:13.b2')
      READ( SUBSTRS(1:N), * ) LADJ_DDEP_CSPEC

      ! wet deposition LS
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_observation_menu:13.c' )
      READ( SUBSTRS(1:N), * ) LADJ_WDEP_LS

      ! wet deposition CV
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_observation_menu:13.d' )
      READ( SUBSTRS(1:N), * ) LADJ_WDEP_CV

      ! wet deposition units
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_observation_menu:13.e1')
      READ( SUBSTRS(1:N), * ) LMOLECCM2S
      IF ( LMOLECCM2S .and. LADJ_FDEP ) NUNIT_COUNT = NUNIT_COUNT + 1

      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_observation_menu:13.e2')
      READ( SUBSTRS(1:N), * ) LKGNHAYR
      IF ( LKGNHAYR .and. LADJ_FDEP ) NUNIT_COUNT = NUNIT_COUNT + 1

      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_observation_menu:13.e3')
      READ( SUBSTRS(1:N), * ) LEQHAYR
      IF ( LEQHAYR .and. LADJ_FDEP ) NUNIT_COUNT = NUNIT_COUNT + 1

      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_observation_menu:13.e4')
      READ( SUBSTRS(1:N), * ) LKGS
      IF ( LKGS .and. LADJ_FDEP ) NUNIT_COUNT = NUNIT_COUNT + 1

      ! now we can define units of the deposition cost function
      ! fp 3/10/2013
      IF ( LADJ_FDEP ) THEN
         IF ( LKGNHAYR ) THEN
            DEP_UNIT = TRIM( 'kgN/ha/yr' )
         ELSEIF ( LEQHAYR ) THEN
            DEP_UNIT = TRIM( 'eq/ha/yr' )
         ELSEIF ( LMOLECCM2S ) THEN
            DEP_UNIT = TRIM( 'molec/cm2/s' )
         ELSEIF ( LKGS ) THEN
            DEP_UNIT = TRIM( 'kg/s' )
         ELSE
            CALL ERROR_STOP(' No unit selected for deposition ',
     &                   ' input_adj_mod.f ')
         ENDIF
      ELSE
         !set all deposition switches to false to avoid unwanted behavior with fd tests (fp)
         LKGS        = .FALSE.
         LEQHAYR     = .FALSE.
         LKGNHAYR    = .FALSE.
         LMOLECCM2S  = .FALSE.

         LADJ_DDEP_TRACER = .FALSE.
         LADJ_DDEP_CSPEC  = .FALSE.
         LADJ_WDEP_CV     = .FALSE.
         LADJ_WDEP_LS     = .FALSE.

      ENDIF

      ! Make sure that we haven't defined too many
      IF ( NUNIT_COUNT > 1 ) THEN
         CALL ERROR_STOP(' More than one choice for cost function ',
     &                   ' input_adj_mod.f ')


      ! Make sure that we have picked at least one.  For
      ! FD tests, the default is forced to be kg/box.
      ELSEIF ( NUNIT_COUNT == 0 .and. LSENS .and. ( .not. LFDTEST ) )
     &   THEN
         CALL ERROR_STOP(' Need to choose one option for units    ',
     &                   ' input_adj_mod.f ')
      ENDIF
      ! Make sure that if deposition is selected that at least one option
      ! is turned on.
      IF ( LADJ_FDEP .and. ( .not. LADJ_FDEP    )
     &               .and. ( .not. LADJ_WDEP_LS )
     &               .and. ( .not. LADJ_WDEP_CV ) ) THEN
         CALL ERROR_STOP(' No deposition option selected ',
     &                   ' input_adj_mod.f ')
      ENDIF

      ! Regional mask?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_observation_menu:13.f' )
      READ( SUBSTRS(1:N), * ) LFORCE_MASK

      IF ( LFORCE_MASK ) THEN

         !fp add option for nc file

         CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,
     &        'read_observation_menu:13.b1')

         READ( SUBSTRS(1:N), * ) LFORCE_MASK_BPCH

         IF ( LFORCE_MASK_BPCH ) THEN

            CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,
     &           'read_observation_menu:13.b2' )

            READ( SUBSTRS(1:N), '(a)' ) FORCING_MASK_FILE
            CALL CHECK_FILE( FORCING_MASK_FILE )

         ELSE

            CALL SPLIT_ONE_LINE( SUBSTRS, N, -1,
     &           'read_observation_menu:13.b2' )

         ENDIF

         CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,
     &        'read_observation_menu:13.b3')

         READ( SUBSTRS(1:N), * ) LFORCE_MASK_NC

         IF ( LFORCE_MASK_NC ) THEN

            CALL SPLIT_ONE_LINE( SUBSTRS, N, 1,
     &           'read_observation_menu:13.b4' )

            READ( SUBSTRS(1:N), '(a)' ) FORCING_MASK_FILE_NC

            CALL CHECK_FILE( FORCING_MASK_FILE_NC )

            CALL SPLIT_ONE_LINE( SUBSTRS, NB_MASK_VAR, -1,
     &           'read_observation_menu:13.b5' )

            ALLOCATE( FORCING_MASK_VARIABLE( NB_MASK_VAR ), STAT = AS )

            DO N = 1, NB_MASK_VAR
               FORCING_MASK_VARIABLE( N ) = TRIM( SUBSTRS(N) )
            ENDDO

         ELSE

            CALL SPLIT_ONE_LINE( SUBSTRS, N, -1,
     &           'read_observation_menu:13.b5' )

            CALL SPLIT_ONE_LINE( SUBSTRS, N, -1,
     &           'read_observation_menu:13.b5' )

         ENDIF

      ELSE

         ! skip lines
         ! bpch switch
         CALL SPLIT_ONE_LINE( SUBSTRS, N, -1,
     &      'read_control_vars_menu:13.a1')

         ! bpch file
         CALL SPLIT_ONE_LINE( SUBSTRS, N, -1,
     &      'read_control_vars_menu:13.a2')

         ! nc switch
         CALL SPLIT_ONE_LINE( SUBSTRS, N, -1,
     &      'read_control_vars_menu:13.a3')

         ! nc file
         CALL SPLIT_ONE_LINE( SUBSTRS, N, -1,
     &      'read_control_vars_menu:13.a4')

         ! nc variable
         CALL SPLIT_ONE_LINE( SUBSTRS, N, -1,
     &      'read_control_vars_menu:13.a5')

      ENDIF

      IF ( LFORCE_MASK ) THEN
         IF ( LFORCE_MASK_BPCH .and. LFORCE_MASK_NC ) THEN
            CALL ERROR_STOP(' Two mask files are defined',
     &           ' input_adj_mod.f ')
         ENDIF

         IF ( .not. LFORCE_MASK_BPCH .and. .not. LFORCE_MASK_NC ) THEN
            CALL ERROR_STOP(' No mask file is defined',
     &           ' input_adj_mod.f ')
         ENDIF

      ENDIF

      !=================================================================
      ! Tracer observations
      !=================================================================

      ! Separator line: >------------------------------<
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_control_vars_menu:11b')


      ! Number of tracers to observe
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_observation_menu:12' )
      READ( SUBSTRS(1:N), * ) NOBS

      ! Separator line:   => obs these tracers------>  : TRC# tracer_name
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_control_vars_menu:13' )

      IF ( NOBS > 0 ) THEN
         DO T = 1, NOBS

            ! Split line into substrings
            CALL SPLIT_ONE_LINE( SUBSTRS, N, -1,
     &      'read_control_vars_menu:14')

            ! tracer id number
            READ( SUBSTRS(1), *) TMP

            ! tracer name
            READ( SUBSTRS(2), *) TNAME(TMP)

            ! observe this species?
            OBS_THIS_TRACER(TMP) = .TRUE.

            ! track tracer index
            TRACER_IND(T) = TMP

         ENDDO

         ! Separator line
         CALL SPLIT_ONE_LINE( SUBSTRS, N, -1,
     &      'read_control_vars_menu:14b')

      ELSE

         ! Loop until at the next section
         DO

            ! Read a line from the file
            LINE = READ_ONE_LINE( EOF )

            ! Stop reading lines when we've passed the Tracer section
            IF ( .not. (INDEX( LINE, 'Tracer' ) > 0 ) ) EXIT

         ENDDO

      ENDIF

      !=================================================================
      ! Species observations
      !=================================================================

      ! Number of species to observe
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_observation_menu:15' )
      READ( SUBSTRS(1:N), * ) NOBS_CSPEC

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_control_vars_menu:15b')

      IF ( NOBS_CSPEC > 0 ) LCSPEC_OBS = .TRUE.

      IF ( ITS_A_FULLCHEM_SIM() .and. LCSPEC_OBS ) THEN

         ! First allocate OBS_THIS_SPECIES to be max species
         ALLOCATE( OBS_THIS_SPECIES( NOBS_CSPEC ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'OBS_THIS_SPECIES' )
         OBS_THIS_SPECIES = 0

         !
         ALLOCATE( CNAME( NOBS_CSPEC ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'CNAME' )
         CNAME = ''


         DO T = 1, NOBS_CSPEC

            ! Split line into substrings
            CALL SPLIT_ONE_LINE( SUBSTRS, N, -1,
     &         'read_observation_menu:17')

            ! Save species name
            CNAME(T)  = TRIM( SUBSTRS(1) )

            ! observe this species?
            OBS_THIS_SPECIES(T) = .TRUE.

         ENDDO

      ENDIF

      !=================================================================
      ! Print to screen
      !=================================================================
      WRITE( 6, '(/,a)' ) 'OBSERVATION MENU'
      WRITE( 6, '(  a)' ) '---------------'
      WRITE( 6, 110     ) 'Observation frequency         : ', OBS_FREQ
      IF ( LFD_GLOB ) THEN
         !print*,' *** FD_GLOB: enforce values on LMAX_OBS and NSPAN ***'
         print*,' *** FD_GLOB: enforce values on LMAX_OBS ***'
      ENDIF
      WRITE( 6, 100     ) 'Limit number of observations  : ', LMAX_OBS
      WRITE( 6, 150     ) 'Forcing time till             : ',
     &     NYMDf, NHMSf
      WRITE( 6, 110     ) '                     NSPAN    => ', NSPAN
      WRITE( 6, '(  a)' ) 'Cost function options         :--- '
      WRITE( 6, 100     ) ' tracer kg/box                : ', LKGBOX
      WRITE( 6, 100     ) ' tracer ug/m3                 : ', LUGM3
      WRITE( 6, 100     ) ' tracer ppb                   : ', LSTT_PPB
      WRITE( 6, 100     ) ' tracer ppm free trop         : ',
     &   LSTT_TROP_PPM
      WRITE( 6, 100     ) ' species ppb w/averaging      : ', LCSPEC_PPB
      WRITE( 6, 100     ) ' tracer ug/m3 pop weighted    : ', LPOP_UGM3
      WRITE( 6, 100     ) ' deposition based?            : ', LADJ_FDEP
      WRITE( 6, 100     ) '   => tracer  dry dep         : ',
     &   LADJ_DDEP_TRACER
      WRITE( 6, 100     ) '   => species dry dep         : ',
     &   LADJ_DDEP_CSPEC
      WRITE( 6, 100     ) '   => wet LS deposition       : ',
     &   LADJ_WDEP_LS
      WRITE( 6, 100     ) '   => wet CV deposition       : ',
     &   LADJ_WDEP_CV
      IF (LADJ_FDEP) THEN
        WRITE( 6, 140     ) '   Deposition                 : ', DEP_UNIT
      ELSE
        WRITE( 6, 140     ) '   Deposition                 :     NONE '
      ENDIF
      WRITE( 6, 100     ) ' Regional forcing mask?       : ',
     &   LFORCE_MASK
      IF ( LFORCE_MASK ) THEN
         IF ( LFORCE_MASK_BPCH ) THEN
            WRITE( 6, 140     ) '   => mask name               : ',
     &           TRIM(FORCING_MASK_FILE)
         ELSEIF ( LFORCE_MASK_NC ) THEN
            WRITE( 6, 140     ) '   => mask name               : ',
     &           TRIM(FORCING_MASK_FILE_NC)
            DO N = 1,NB_MASK_VAR
               WRITE( 6, 140     ) '   => varname                 : ',
     &              TRIM(FORCING_MASK_VARIABLE(N))
            ENDDO
         ENDIF
      ELSE
         WRITE( 6, 140     ) '   => mask name               :     ',
     &      'NOT USED'
      ENDIF
      WRITE( 6, '(  a)' ) '>------------------------------<'
      WRITE( 6, 110     ) 'Number of tracers to observe  : ', NOBS

      IF ( NOBS > 0 ) THEN
         WRITE( 6, '(  a)' ) '  Tracers to observe          '

         ! Print info about each tracer
         DO T = 1, N_TRACERS

            IF( OBS_THIS_TRACER(T) ) THEN
             ! Write tracer number, name and if it's observed
             WRITE( 6, 130 )  T, TNAME(T)
          ENDIF

         ENDDO

      ENDIF

      IF ( LCSPEC_OBS ) THEN
         WRITE( 6, '(  a)' ) REPEAT( '-', 48 )
         WRITE( 6, 110     ) 'Number of species to observe  : ',
     &      NOBS_CSPEC
         WRITE( 6, '(  a)' ) '  Species to observe      '

         ! Print info about each tracer
         DO T = 1, NOBS_CSPEC

             ! Write tracer number, name and if it's observed
             WRITE( 6, 120 )  T, CNAME(T)

         ENDDO

      ENDIF



 100  FORMAT( A, L5             )
 110  FORMAT( A, I5             )
 120  FORMAT( I3, 1x, A10       )
 130  FORMAT( I3, 1x, A10, 6x, I5 )
 140  FORMAT( A, A              )
 150  FORMAT( A, I8, 1x, I6     )

      ! Set counter
      CT1 = CT1 + 1


      END SUBROUTINE READ_OBSERVATION_MENU

!---------------------------------------------------------------------------------------

      SUBROUTINE READ_FD_MENU
!
!******************************************************************************
!  Subroutine READ_FD_MENU reads the FINITE DIFFERENCE MENU section of
!  the GEOS-CHEM adj input file (adj_group, 6/08/09)
!
!  NOTES:
!  (1 ) Add support for strat fluxes LADJ_STRAT (hml, dkh, 02/14/12, adj32_025)
!******************************************************************************
!
      ! References to F90 modules
      USE ADJ_ARRAYS_MOD,    ONLY : FD_DIFF
      USE ADJ_ARRAYS_MOD,    ONLY : LONFD
      USE ADJ_ARRAYS_MOD,    ONLY : LATFD
      USE ADJ_ARRAYS_MOD,    ONLY : IFD
      USE ADJ_ARRAYS_MOD,    ONLY : JFD
      USE ADJ_ARRAYS_MOD,    ONLY : LFD
      USE ADJ_ARRAYS_MOD,    ONLY : NFD
      USE ADJ_ARRAYS_MOD,    ONLY : ICSFD
      USE ADJ_ARRAYS_MOD,    ONLY : EMSFD
      USE ADJ_ARRAYS_MOD,    ONLY : MFD
      USE ADJ_ARRAYS_MOD,    ONLY : STRFD
      USE DRYDEP_MOD,        ONLY : NTRAIND
      USE DRYDEP_MOD,        ONLY : NUMDEP
      USE ERROR_MOD,         ONLY : ERROR_STOP
      USE GRID_MOD,          ONLY : GET_BOUNDING_BOX
      USE GRID_MOD,          ONLY : GET_XMID
      USE GRID_MOD,          ONLY : GET_YMID
      USE LOGICAL_ADJ_MOD,   ONLY : LPRINTFD
      USE LOGICAL_ADJ_MOD,   ONLY : LFD_SPOT
      USE LOGICAL_ADJ_MOD,   ONLY : LFD_GLOB
      USE LOGICAL_ADJ_MOD,   ONLY : LFDTEST
      USE LOGICAL_ADJ_MOD,   ONLY : LADJ_FDEP
      USE LOGICAL_ADJ_MOD,   ONLY : LADJ_DDEP_TRACER
      USE LOGICAL_ADJ_MOD,   ONLY : LADJ_DDEP_CSPEC
      USE LOGICAL_ADJ_MOD,   ONLY : LADJ_WDEP_CV
      USE LOGICAL_ADJ_MOD,   ONLY : LADJ_WDEP_LS
      USE LOGICAL_MOD,       ONLY : LTRAN
      USE WETSCAV_MOD,       ONLY : GET_WETDEP_IDWETD
      USE WETSCAV_MOD,       ONLY : NSOL
      USE ADJ_ARRAYS_MOD,    ONLY : RATFD

      ! Local variables
      INTEGER            :: N
      REAL*8             :: tmpbox(4)
      INTEGER            :: tmpbox1(4)
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)
      LOGICAL            :: USEINDEX = .FALSE.
      INTEGER            :: IFDTMP, JFDTMP

      !=================================================================
      ! READ_FD_MENU begins here!
      !=================================================================

      ! FD difference size
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_fd_menu:1' )
      READ( SUBSTRS(1:N), * ) FD_DIFF

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_fd_menu:1.5' )

      ! if we're doing global check, then exit
      ! But it is still nice to define IFD, JFD, LFD, etc, if LPRINTFD
      ! is on. Returning here makes these ind undefined,
      ! which lead to seg faults  (dkh, 06/11/09)
      !IF ( LFD_GLOB ) THEN
      !   PRINT*, 'All gridboxes are used in the global FD test'
      !   RETURN
      !ENDIF

      ! longitude of the FD gridbox
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_fd_menu:5' )
      READ( SUBSTRS(1:N), * ) LONFD

      ! latitude of the FD gridbox
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_fd_menu:6' )
      READ( SUBSTRS(1:N), * ) LATFD

      ! check if we're specifying indecies (as opposed to lat/lon)
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_fd_menu:2' )
      READ( SUBSTRS(1:N), * ) USEINDEX

      ! IFD gridbox
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_fd_menu:3' )
      READ( SUBSTRS(1:N), * ) IFDTMP

      ! JFD gridbox
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_fd_menu:4' )
      READ( SUBSTRS(1:N), * ) JFDTMP

      ! get corresponding box indecies for the LONFD and LATFD
      tmpbox(1) = LONFD
      tmpbox(2) = LATFD
      tmpbox(3) = LONFD
      tmpbox(4) = LATFD

      ! Move this below, as it doesn't work with nested domain  (dkh, 01/19/12, adj32_015 )
      !CALL GET_BOUNDING_BOX(tmpbox,tmpbox1)

      IF ( USEINDEX ) THEN
         IFD = IFDTMP
         JFD = JFDTMP

         ! now also adjust LONFD and LATFD (dkh, 02/11/11)
         LONFD = GET_XMID( IFD )
         LATFD = GET_YMID( JFD )

      ELSE

         ! Moved here (dkh, 01/19/12, adj32_015)
         CALL GET_BOUNDING_BOX(tmpbox,tmpbox1)

          IFD = tmpbox1(1)
          JFD = tmpbox1(2)
      ENDIF

      ! FD perturbation box level
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_fd_menu:7' )
      READ( SUBSTRS(1:N), * ) LFD

      ! FD perturbation species
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_fd_menu:8' )
      READ( SUBSTRS(1:N), * ) NFD

      ! Separator line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_fd_menu:8.5' )


      ! FD perturbation box temporal element
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_fd_menu:9' )
      READ( SUBSTRS(1:N), * ) MFD

      ! FD perturbation species
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_fd_menu:10' )
      READ( SUBSTRS(1:N), * ) EMSFD

      ! FD perturbation species
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_fd_menu:11' )
      READ( SUBSTRS(1:N), * ) ICSFD

      ! FD perturbation species (hml, 08/11/11, adj32_025)
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_fd_menu:12' )
      READ( SUBSTRS(1:N), * ) STRFD

      ! FD perturbation rate (tww, 05/15/12)
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_fd_menu:13' )
      READ( SUBSTRS(1:N), * ) RATFD

      ! Move these to adjoint menu       (dkh, 02/09/11)
      !! Doing finite difference test in 1 gridbox
      !CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_fd_menu:12' )
      !READ( SUBSTRS(1:N), * ) LFD_SPOT
      !
      !! Doing finite difference test in all grid boxes, turn transport off
      !CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_fd_menu:13' )
      !READ( SUBSTRS(1:N), * ) LFD_GLOB
      !
      !! turn of transport for global FD test
      !IF ( LFD_GLOB ) LTRAN = .FALSE.
      !
      !! define a more generic LFDTEST flag if either method is true
      !IF ( LFD_GLOB .OR. LFD_SPOT ) LFDTEST = .TRUE.

      !=================================================================
      ! Print to screen
      !=================================================================
      WRITE( 6, '(/,a)' ) 'FINITE DIFFERENCE MENU'
      WRITE( 6, '(  a)' ) '---------------'
      WRITE( 6, 100     ) 'Finite diff. increment FD_DIFF: ', FD_DIFF
      WRITE( 6, 120     ) 'Finite diff longitude LONFD   : ', LONFD
      WRITE( 6, 110     ) 'Finite diff long. index IFD   : ', IFD
      WRITE( 6, 120     ) 'Finite diff latitude LATFD    : ', LATFD
      WRITE( 6, 110     ) 'Finite diff lat. index JFD    : ', JFD
      WRITE( 6, 110     ) 'Finite diff vert index LFD    : ', LFD
      WRITE( 6, 110     ) 'FD species NFD                : ', NFD
      WRITE( 6, 110     ) 'FD time.group index MFD       : ', MFD
      WRITE( 6, 110     ) 'FD emiss EMSFD                : ', EMSFD
      WRITE( 6, 110     ) 'FD initial cond ICSFD         : ', ICSFD
      WRITE( 6, 110     ) 'FD strat prod & loss STRFD    : ', STRFD
      WRITE( 6, 110     ) 'FD reaction rate RATFD        : ', RATFD
      !WRITE( 6, 130     ) 'Doing finite diff check (1box): ', LFD_SPOT
      !WRITE( 6, 130     ) 'Doing finite diff check (glob): ', LFD_GLOB

      ! Format statements
 100  FORMAT( A, f11.6          )
 110  FORMAT( A, I4             )
 120  FORMAT( A, f7.2          )
 130  FORMAT( A, L5             )

      !=================================================================
      ! Call setup routines from other GEOS-CHEM modules
      !=================================================================

      ! Set counter
      CT1 = CT1 + 1

      ! Return to calling program
      END SUBROUTINE READ_FD_MENU

!------------------------------------------------------------------------------

      SUBROUTINE READ_ADJ_DIAGNOSTICS_MENU
!
!******************************************************************************
!  Subroutine READ_ADJ_DIAGNOSTICS_MENU reads the DIAGNOSTICS MENU section of
!  the GEOS-CHEM adjoint input file (adj_group, 6/07/09)
!
!  NOTES:
!  (1 ) Add LITR (zhe, dkh, 02/04/11)
!  (2 ) Add LTRAJ_SCALE (dkh, 02/09/11)
!  (3 ) Add LEMS_ABS, LTES_BLVMR (dkh, 02/17/11)
!  (4 ) Add support for strat fluxes LADJ_STRAT (hml, dkh, 02/14/12, adj32_025)
!******************************************************************************
!
      ! References to F90 modules
      USE CHARPAK_MOD,     ONLY : STRREPL
      USE LOGICAL_ADJ_MOD, ONLY : LADJDIAG
      USE LOGICAL_ADJ_MOD, ONLY : LJSAVE
      USE LOGICAL_ADJ_MOD, ONLY : LADJ_TRAJ
      USE LOGICAL_ADJ_MOD, ONLY : LDCOSAT
      USE LOGICAL_ADJ_MOD, ONLY : LHMOD
      USE LOGICAL_ADJ_MOD, ONLY : LhOBS
      USE LOGICAL_ADJ_MOD, ONLY : LHMODIFF
      USE LOGICAL_ADJ_MOD, ONLY : LADJ_FORCE
      USE LOGICAL_ADJ_MOD, ONLY : LMODBIAS
      USE LOGICAL_ADJ_MOD, ONLY : LOBS_COUNT
      USE LOGICAL_ADJ_MOD, ONLY : LDOFS
      USE LOGICAL_ADJ_MOD, ONLY : LPRINTFD
      USE LOGICAL_ADJ_MOD, ONLY : LDEL_CHKPT
      USE LOGICAL_ADJ_MOD, ONLY : LITR
      USE LOGICAL_ADJ_MOD, ONLY : LTRAJ_SCALE
      USE LOGICAL_ADJ_MOD, ONLY : LTES_BLVMR
      USE LOGICAL_ADJ_MOD, ONLY : LEMS_ABS
      USE LOGICAL_ADJ_MOD, ONLY : LSAT_HDF_L2, LSAT_HDF_L3

      ! Local variables
      INTEGER                  :: N
      CHARACTER(LEN=255)       :: SUBSTRS(MAXDIM)

      ! (dkh, 01/09/12, adj32_010)
      LOGICAL                  :: EOF
      INTEGER                  :: IOS
      CHARACTER(LEN=1)         :: TAB   = ACHAR(9)
      CHARACTER(LEN=1)         :: SPACE = ' '
      CHARACTER(LEN=255)       :: LINE

      !=================================================================
      ! READ_ADJ_SIMULATION_MENU begins here!
      !=================================================================

      LJSAVE     = .FALSE.
      LADJ_TRAJ  = .FALSE.
      LHMOD      = .FALSE.
      LhOBS      = .FALSE.
      LHMODIFF   = .FALSE.
      LADJ_FORCE = .FALSE.
      LMODBIAS   = .FALSE.
      LOBS_COUNT = .FALSE.
      LDOFS      = .FALSE.
      LITR       = .FALSE.
      LTRAJ_SCALE= .FALSE.
      LTES_BLVMR = .FALSE.
      LEMS_ABS   = .FALSE.
      LSAT_HDF_L2= .FALSE.
      LSAT_HDF_L3= .FALSE.


      ! Save any diagnostics? If not, exit subroutine with all flags FALSE
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_adj_diag_menu:1' )
      READ( SUBSTRS(1:N), * ) LADJDIAG

      IF ( .NOT. LADJDIAG ) THEN
         WRITE( 6, '(/,a)' ) 'SKIPPING DIAGNOSTICS MENU'
         RETURN
      ENDIF

      ! PRINT debug messages in FD cell files
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_adj_diag_menu:2' )
      READ( SUBSTRS(1:N), * ) LPRINTFD

      ! Move to other menu (dkh, 02/09/11)
      !! Delete checkpt files
      !CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_adj_diag_menu:2.1' )
      !READ( SUBSTRS(1:N), * ) LDEL_CHKPT

      ! SAVE .save and .sav2 files
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_adj_diag_menu:3' )
      READ( SUBSTRS(1:N), * ) LJSAVE

      ! Save adjoint trajectory files
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_adj_diag_menu:4' )
      READ( SUBSTRS(1:N), * ) LADJ_TRAJ

      ! save STT adjoints as scaling factor sensitivities?
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_adj_diag_menu:4.0' )
      READ( SUBSTRS(1:N), * ) LTRAJ_SCALE

      ! Save iteration information
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_adj_diag_menu:4.1' )
      READ( SUBSTRS(1:N), * ) LITR

      ! Save sense w.r.t absolute emis
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_adj_diag_menu:4.2' )
      READ( SUBSTRS(1:N), * ) LEMS_ABS

      ! CO satellite diagnostics? if not, don't read the next 7 lines
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_adj_diag_menu:5' )
      READ( SUBSTRS(1:N), * ) LDCOSAT

      IF ( LDCOSAT ) THEN

         ! Save H(model), model *ak
         CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_adj_diag_menu:6' )
         READ( SUBSTRS(1:N), * ) LHMOD

         ! Save h(obs), gridded and filtered observations
         CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_adj_diag_menu:7' )
         READ( SUBSTRS(1:N), * ) LhOBS

         ! Save H(mod) - h(obs)
         CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_adj_diag_menu:8' )
         READ( SUBSTRS(1:N), * ) LHMODIFF

         ! Save adjoint forcing
         CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_adj_diag_menu:9' )
         READ( SUBSTRS(1:N), * ) LADJ_FORCE

         ! Save model bias (H(model)-h(obs))/h(obs)
         CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_adj_diag_menu:10' )
         READ( SUBSTRS(1:N), * ) LMODBIAS

         ! Save observation count (array with count/box)
         CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_adj_diag_menu:11' )
         READ( SUBSTRS(1:N), * ) LOBS_COUNT

         ! Save gridded DOFs
         CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_adj_diag_menu:12' )
         READ( SUBSTRS(1:N), * ) LDOFs

      !----------------------------------------------------------------
      ! BUG FIX: Allow for proper reading of menu below the CO sub menu
      !  (dkh, 01/08/12, adj32_010)
      ! OLD CODE:
      !ENDIF
      !
      !! Separator line: TES NH3 diagnostics
      !CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_fd_menu:8.5' )
      ! NEW CODE:
      ELSE

         DO WHILE ( INDEX( LINE, 'TES NH3 diagnostics'  ) .le. 0 )

            ! still need to advance through the file
            LINE = READ_ONE_LINE( EOF )
            IF ( EOF ) EXIT

            ! Replace tab characters in LINE (if any) w/ spaces
            CALL STRREPL( LINE, TAB, SPACE )

         ENDDO

      ENDIF
      !----------------------------------------------------------------

      ! Save BLVMR
      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_adj_diag_menu:12' )
      READ( SUBSTRS(1:N), * ) LTES_BLVMR

      ! Separator line: >------------------------------<
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_control_vars_menu:13')

      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_adj_diag_menu:14' )
      READ( SUBSTRS(1:N), * ) LSAT_HDF_L2

      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_adj_diag_menu:15' )
      READ( SUBSTRS(1:N), * ) LSAT_HDF_L3

      !=================================================================
      ! Print to screen
      !=================================================================
      WRITE( 6, '(/,a)' ) 'DIAGNOSTICS MENU'
      WRITE( 6, '(  a)' ) '---------------'
      WRITE( 6, 100     ) 'Print adj debug LPRINTFD      : ', LPRINTFD
      !WRITE( 6, 100     ) 'Delete chkpt files LDEL_CHKPT : ', LDEL_CHKPT
      WRITE( 6, 100     ) 'Save .jsave and .jsave2 files : ', LJSAVE
      WRITE( 6, 100     ) 'Adjoint trajectory files      : ', LADJ_TRAJ
      WRITE( 6, 100     ) '   w.r.t. scaling factors     : ',
     &   LTRAJ_SCALE
      WRITE( 6, 100     ) 'Save iteration diagnostics    : ', LITR
      WRITE( 6, 100     ) 'Save sense w.r.t absolute emis: ', LEMS_ABS
      IF ( LEMS_ABS ) PRINT*, ' ### WARNING: LEMS_ABS only for SO2, BC'
      WRITE( 6, 100     ) 'Save CO sat. diagnostics      : ', LDCOSAT

      IF ( LDCOSAT) THEN
      WRITE( 6, 100     ) 'Save H(model)                 : ', LHMOD
      WRITE( 6, 100     ) 'Save h(obs)                   : ', LhOBS
      WRITE( 6, 100     ) 'Save H(model)-h(obs)          : ', LHMODIFF
      WRITE( 6, 100     ) 'Save adjoint forcing          : ', LADJ_FORCE
      WRITE( 6, 100     ) 'Save model bias               : ', LMODBIAS
      WRITE( 6, 100     ) 'Save number of obs/gridbox    : ', LOBS_COUNT
      WRITE( 6, 100     ) 'Save gridded DOFs             : ', LDOFS
      ENDIF

      WRITE( 6, 100     ) 'TES NH3 BLVMR                 : ', LTES_BLVMR
      WRITE( 6, 100     ) 'HDF Level 2                   : ',LSAT_HDF_L2
      WRITE( 6, 100     ) 'HDF Level 3                   : ',LSAT_HDF_L3

      ! Format statements
 100  FORMAT( A, L5             )

      !=================================================================
      ! Call setup routines from other GEOS-CHEM modules
      !=================================================================

      ! Set counter
      CT1 = CT1 + 1

      ! Return to calling program
      END SUBROUTINE READ_ADJ_DIAGNOSTICS_MENU

!------------------------------------------------------------------------------

      SUBROUTINE READ_ADJ_CRITICAL_LOAD_MENU
!fp
!add new menu to streamline the inputs for critical load sensitivity simulations

      USE CRITICAL_LOAD_MOD, ONLY : CL_FILENAME
      USE CRITICAL_LOAD_MOD, ONLY : GC_FILENAME
      USE LOGICAL_ADJ_MOD,   ONLY : LADJ_CL
      USE LOGICAL_ADJ_MOD,   ONLY : LADJ_CL_NDEP
      USE LOGICAL_ADJ_MOD,   ONLY : LADJ_CL_ACID

#     include "CMN_SIZE"

      INTEGER            :: N
      CHARACTER(LEN=255) :: SUBSTRS(MAXDIM)
      REAL*8              :: MASK(IIPAR,JJPAR)


      !=================================================================
      ! READ_ADJ_CRITICAL_LOAD begins here!
      !=================================================================

      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_critical_load_menu:1' )
      READ( SUBSTRS(1:N), * ) LADJ_CL

      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_critical_load_menu:2' )
      READ( SUBSTRS(1:N), * ) LADJ_CL_NDEP

      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_critical_load_menu:3' )
      READ( SUBSTRS(1:N), * ) LADJ_CL_ACID

      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_critical_load_menu:2' )
      READ( SUBSTRS(1:N), '(a)' ) CL_FILENAME

      CALL SPLIT_ONE_LINE( SUBSTRS, N, 1, 'read_critical_load_menu:3' )
      READ( SUBSTRS(1:N), '(a)' ) GC_FILENAME


      WRITE( 6, '(/,a)' ) 'CRITICAL LOAD MENU'
      WRITE( 6, '(  a)' ) '---------------'
      WRITE( 6, 100     ) 'Doing critical load run       : ',
     &     LADJ_CL
      WRITE( 6, 100     ) ' => based on N deposition     : ',
     &     LADJ_CL_NDEP
      WRITE( 6, 100     ) ' => based on acid deposition  : ',
     &     LADJ_CL_ACID
      WRITE( 6, '( a)'  ) ' Critical Load base file      : ',
     &     TRIM(CL_FILENAME)
      WRITE( 6, '( a)'  ) ' GC Load file                 : ',
     &     TRIM(GC_FILENAME)

 100  FORMAT( A, L5             )

      ! Return to calling program
      END SUBROUTINE READ_ADJ_CRITICAL_LOAD_MENU

!------------------------------------------------------------------------------

      SUBROUTINE ARE_FLAGS_VALID(  )
!
!******************************************************************************
!  Subroutine ARE_FLAGS_VALID checks to make sure that flags for the forward
!  calculation (set in input.geos) do not confict with flags for the adjoint
!  calculation (set in input.gcadj ). (dkh, 11/02/05, adj_group 6/07/09)
!
!  NOTES:
!  (1 ) Add support for strat fluxes LADJ_STRAT (hml, dkh, 02/14/12, adj32_025)
!  (2 ) Add LINVH_BFGS (nab, 25/03/12)
!
!******************************************************************************
!
      ! Reference to f90 modules
      USE ADJ_ARRAYS_MOD,  ONLY : OBS_FREQ
      USE ADJ_ARRAYS_MOD,  ONLY : IFD, JFD,  LFD, NFD
      USE ADJ_ARRAYS_MOD,  ONLY : ICSFD
      USE ADJ_ARRAYS_MOD,  ONLY : NNEMS,     EMSFD
      ! USE ADJ_ARRAYS_MOD,  ONLY : MMSCL
      USE ADJ_ARRAYS_MOD,  ONLY : N_CALC,    N_CALC_STOP
      USE ADJ_ARRAYS_MOD,  ONLY : IDADJ_EBCPI_an, IDADJ_EBCPO_an
      USE ADJ_ARRAYS_MOD,  ONLY : IDADJ_EOCPI_an, IDADJ_EOCPO_an
      USE ADJ_ARRAYS_MOD,  ONLY : IDADJ_EBCPI_bb, IDADJ_EBCPO_bb
      USE ADJ_ARRAYS_MOD,  ONLY : IDADJ_EOCPI_bb, IDADJ_EOCPO_bb
      USE ADJ_ARRAYS_MOD,  ONLY : IDADJ_EBCPI_bf, IDADJ_EBCPO_bf
      USE ADJ_ARRAYS_MOD,  ONLY : IDADJ_EOCPI_bf, IDADJ_EOCPO_bf
      USE ADJ_ARRAYS_MOD,  ONLY : IDADJ_ENH3_an
      USE ADJ_ARRAYS_MOD,  ONLY : IDADJ_ENH3_na
      USE ADJ_ARRAYS_MOD,  ONLY : IDADJ_ENH3_bb
      USE ADJ_ARRAYS_MOD,  ONLY : IDADJ_ENH3_bf
      USE ADJ_ARRAYS_MOD,  ONLY : IDADJ_ESO2_bf
      USE ADJ_ARRAYS_MOD,  ONLY : IDADJ_ESO2_bb
      USE ADJ_ARRAYS_MOD,  ONLY : IDADJ_ESO2_sh
      USE ADJ_ARRAYS_MOD,  ONLY : IDADJ_ESO2_an1
      USE ADJ_ARRAYS_MOD,  ONLY : IDADJ_ESO2_an2
      USE ADJ_ARRAYS_MOD,  ONLY : IDADJ_EDST1,    IDADJ_EDST2
      USE ADJ_ARRAYS_MOD,  ONLY : IDADJ_EDST3,    IDADJ_EDST4
      USE ADJ_ARRAYS_MOD,  ONLY : N_CARB_EMS_ADJ
      USE ADJ_ARRAYS_MOD,  ONLY : N_SULF_EMS_ADJ
      USE ADJ_ARRAYS_MOD,  ONLY : N_DUST_EMS_ADJ
      USE ADJ_ARRAYS_MOD,  ONLY : IS_CARB_EMS_ADJ
      USE ADJ_ARRAYS_MOD,  ONLY : IS_SULF_EMS_ADJ
      USE ADJ_ARRAYS_MOD,  ONLY : IS_DUST_EMS_ADJ
      USE ADJ_ARRAYS_MOD,  ONLY : NOBS_CSPEC
      USE ADJ_ARRAYS_MOD,  ONLY : IDCSPEC_ADJ
      USE ADJ_ARRAYS_MOD,  ONLY : OBS_THIS_SPECIES
      USE ADJ_ARRAYS_MOD,  ONLY : OBS_THIS_TRACER
      USE ADJ_ARRAYS_MOD,  ONLY : OPT_THIS_EMS
      USE ADJ_ARRAYS_MOD,  ONLY : OPT_THIS_TRACER
      USE ADJ_ARRAYS_MOD,  ONLY : EMS_ERROR
      USE ADJ_ARRAYS_MOD,  ONLY : ICS_ERROR
      USE ADJ_ARRAYS_MOD,  ONLY : NOBS
      USE ADJ_ARRAYS_MOD,  ONLY : NSTPL,     STRFD
      USE ADJ_ARRAYS_MOD,  ONLY : OPT_THIS_PROD
      USE ADJ_ARRAYS_MOD,  ONLY : OPT_THIS_LOSS
      USE ADJ_ARRAYS_MOD,  ONLY : PROD_ERROR
      USE ADJ_ARRAYS_MOD,  ONLY : LOSS_ERROR
      USE ADJ_ARRAYS_MOD,  ONLY : NRRATES, RATFD
      USE ADJ_ARRAYS_MOD,  ONLY : OBS_FREQ
      USE ADJ_ARRAYS_MOD,  ONLY : CNAME
      USE ADJ_ARRAYS_MOD,  ONLY : NSPAN
      USE ERROR_MOD,       ONLY : ERROR_STOP
      USE LOGICAL_MOD,     ONLY : LDRYD,     LCHEM,       LTURB,
     &                            LCHEM,     LWETD,       LTRAN,
     &                            LCONV,     LSOILNOX,    LSCHEM
      USE LOGICAL_ADJ_MOD, ONLY : LADJ_CHEM, LAERO_THERM, LADJ_TRAN,
     &                            LSENS,     LFDTEST,     L4DVAR,
     &                            LICS,      LADJ_EMS,    LFD_GLOB,
     &                            LBKCOV,    LADJ,        LLINOZ,
     &                            L3DVAR,    LCSPEC_PPB,  LCSPEC_OBS,
     &                            LEMS_ABS,  LAPSRC,      LINVH,
     &                            LINVH_BFGS,
     &                            LADJ_STRAT, LADJ_RRATE,
     &                            LADJ_FDEP,
     &                            LADJ_DDEP_TRACER,
     &                            LADJ_DDEP_CSPEC,
     &                            LADJ_WDEP_LS,
     &                            LADJ_WDEP_CV,
     &                            LMAX_OBS
      USE ADJ_ARRAYS_MOD,  ONLY : OBS_FREQ
      USE LOGICAL_ADJ_MOD, ONLY : LKGNHAYR
      USE LOGICAL_ADJ_MOD, ONLY : LEQHAYR
      USE LOGICAL_ADJ_MOD, ONLY : LMOLECCM2S
      USE LOGICAL_ADJ_MOD, ONLY : LKGS
      USE TIME_MOD,        ONLY : GET_TIME_AHEAD
      USE TIME_MOD,        ONLY : GET_NYMDe
      USE TIME_MOD,        ONLY : GET_NHMSe
      USE TRACER_MOD,      ONLY : N_TRACERS, ITS_A_FULLCHEM_SIM
      USE TRACER_MOD,      ONLY : ITS_A_TAGOX_SIM
      USE TRACER_MOD,      ONLY : SIM_TYPE
      USE TRACERID_MOD,    ONLY : IDTSO4,    IDTDST1,     IDTSOA1
      USE TRACERID_MOD,    ONLY : IDTSALA
      USE TRACERID_MOD,    ONLY : IDTNIT,  IDTNH4,  IDTNH3
      USE TRACERID_MOD,    ONLY : IDTHNO3, IDTSO2
      USE LOGICAL_ADJ_MOD, ONLY : LSAT_HDF_L2, LSAT_HDF_L3
      USE LOGICAL_ADJ_MOD, ONLY : LADJ_CL, LADJ_CL_NDEP, LADJ_CL_ACID
      USE TRACER_MOD,      ONLY : TRACER_NAME
      USE GCKPP_ADJ_GLOBAL, ONLY: NCOEFF_RATE

#     include "CMN_SIZE"        ! Size params
#     include "comode.h"        ! NAMEGAS, SMAL2
#     include "define_adj.h"
#     include "CMN_DIAG"        ! ND44



      ! local variables
      INTEGER                  :: N,T
      CHARACTER(LEN=255)       :: MSG
      INTEGER                  :: COUNT_ON

      INTEGER, PARAMETER       :: N_NDEP = 4 !number of tracers for N deposition
      INTEGER, PARAMETER       :: N_ACID = 2 !number of tracers for acid deposition (on top of N_NDEP tracers)
      INTEGER, PARAMETER       :: N_NDEP_CSPEC = 7
      INTEGER                  :: NDEP(N_NDEP), ACID(N_ACID)
      INTEGER                  :: DATE(2)
      CHARACTER*255            :: NDEP_CSPEC(N_NDEP_CSPEC)
      LOGICAL                  :: FOUND

      !NITS and SO4S are not supported at the moment for wet/dry deposition (fp 1/5/2013)

      NDEP(1) = IDTHNO3
      NDEP(2) = IDTNIT
      NDEP(3) = IDTNH3
      NDEP(4) = IDTNH4
      ACID(1) = IDTSO2
      ACID(2) = IDTSO4

      NDEP_CSPEC(1) = 'DRYHNO3'
      NDEP_CSPEC(2) = 'DRYNO2'
      NDEP_CSPEC(3) = 'DRYPAN'
      NDEP_CSPEC(4) = 'DRYPPN'
      NDEP_CSPEC(5) = 'DRYPMN'
      NDEP_CSPEC(6) = 'DRYN2O5'
      NDEP_CSPEC(7) = 'DRYR4N2'

      !=================================================================
      ! ARE_FLAGS_VALID begins here!
      !=================================================================

      ! check if we are even doing an adjoint run
      IF ( .not. LADJ ) RETURN

      !=================================================================
      ! Check forward model options
      !=================================================================
      ! first check if "input.geos" is set to a supported simulation:
      IF ( SIM_TYPE .NE. 7 .AND.   ! FULL CHEM
     &     SIM_TYPE .NE. 3 .AND.   ! TAGGED CO
     &     SIM_TYPE .NE. 9 .AND.   ! CH4 (kjw, adj32_023)
     &     SIM_TYPE .NE. 6 .and.   ! TAGGED OX  (lzh, 12/12/2009)
     &     SIM_TYPE .NE.10 .and.   ! Offline aerosol (adj32_013)
     &     SIM_TYPE .NE. 12) THEN  ! TAGGED CO2
         CALL ERROR_STOP( ' This simulation is not supported ',
     &        ' ARE_FLAGS_VALID, input_adj_mod.f ' )
      ENDIF

      ! So far only BC and OC will work with the offline aerosol sim,
      ! but not the other aerosols (well, dust might, but untested).
      ! (yhmao, dkh, 01/13/12, adj32_013)
      IF ( SIM_TYPE == 10 ) THEN
         IF ( IDTSO4 .or. IDTSALA .or. IDTSOA1 ) THEN
            CALL ERROR_STOP('offline aero adj only for dust and BC/OC',
     &                      ' ARE_FLAGS_VALID, input_adj_mod.f' )
         ENDIF
      ENDIF

      !=================================================================
      ! Check forward and adjoint process options
      !=================================================================
      ! Much of the relevant aerosol chemistry is DRYDEP, and adjoint
      ! of sulfate chemistry will get called if LADJ_CHEM is true,
      ! so we shouldn't have DRYDEP = FALSE and LADJ_CHEM = TRUE.
      ! Should this depend on LSULF at all?
!      IF ( ( LADJ_CHEM .AND. ( .NOT. LDRYD     ) ) .OR.
!     &     ( LDRYD     .AND. ( .NOT. LADJ_CHEM ) )     )  THEN
      ! I think we can have DRYD w/o chem
      IF (  ITS_A_FULLCHEM_SIM() .AND.
     &      LADJ_CHEM .AND. ( .NOT. LDRYD     ) ) THEN
         CALL ERROR_STOP( ' LADJ_CHEM and LDRYD inconsistent ',
     &                    ' ARE_FLAGS_VALID, input_adj_mod.f ' )
      ENDIF

! Not sure about this, leave it out for now (dkh, 06/24/09)
!      ! Don't know why, but if WETD, and CHEM are only fwd true
!      ! and LADJ_CHEM is only adj true, get error.  Have to turn
!      ! on LTRAN.  ( though just WETD, no adj_chem, no TRAN, seems ok).
!      ! something to do with RH?  I think there may be others that
!      ! require LTRAN.... The error pops up as "Invalid EXTRA", caused
!      ! because TS_DYN is 60.
!      IF ( LCHEM .AND. ( .NOT. LTRAN ) ) THEN
!         CALL ERROR_STOP( ' LCHEM and LTRAN inconsistent ',
!     &                    ' ARE_FLAGS_VALID, geos_chem_mod.f ' )
!      ENDIF

      ! LCHEM controls chemistry in the fwd calc, so need this on
      ! if want aerosol thermo or the rest of chemistry.
      IF ( ( LAERO_THERM .OR. LADJ_CHEM )
     &     .AND. ( .NOT. LCHEM ) ) THEN
         CALL ERROR_STOP( ' LCHEM, LADJ_CHEM, LAERO_THERM inconsistent',
     &                    ' ARE_FLAGS_VALID, input_adj_mod.f ' )
      ENDIF
      ! ... and the opposite...
      IF ( .not. ( LADJ_CHEM )
     &     .and.  ( LCHEM ) ) THEN
         CALL ERROR_STOP( ' LADJ_CHEM off but LCHEM is on! ',
     &                    ' ARE_FLAGS_VALID, input_adj_mod.f ' )
      ENDIF

      ! If you have LTURB and LTRAN, but nothing else, adjoints explode.
      ! (dkh, 11/22/05)
      IF ( LTURB .AND. LTRAN .AND. LTRAN .AND. ( .NOT. LCONV )
     &      .AND. ( .NOT. LWETD ) .AND. ( .NOT. LCHEM ) ) THEN
         CALL ERROR_STOP( ' LTURB and LTRAN lead to errors in adj? ',
     &                    ' ARE_FLAGS_VALID, input_adj_mod.f ' )
      ENDIF


      ! Now use new strat_chem_mod (hml, dkh, 02/14/12, adj32_025)
      !! Make sure that if strat fluxes are on, LINOZE adj is on (dkh, 04/25/10)
      !IF ( LUPBD /= LLINOZ ) THEN
      !   CALL ERROR_STOP( ' LUPBD and LLINOZ not consistent ',
      !                    ' ARE_FLAGS_VALID, input_adj_mod.f ' )
      !ENDIF


      ! Only include adjoint w.r.t strat fluxes if strat chem is turned on
      IF ( LADJ_STRAT .and. ( .not. LSCHEM )  ) THEN
         CALL ERROR_STOP( ' LADJ_STRAT needs LSCHEM on ',
     &                    ' ARE_FLAGS_VALID, input_adj_mod.f ' )
      ENDIF

      !=================================================================
      ! Check adjoint simulation type
      !
      ! Overall simulation type must be one and only one of:
      !  - 3DVAR
      !  - 4DVAR
      !  - SENS
      !=================================================================
      ! check at least one:
      IF (  (.not. LSENS ) .and. ( .not. L3DVAR )
     &                     .and. ( .not. L4DVAR ) ) THEN
         MSG = 'Invalid adj run options: no simulation type defined!'
         CALL ERROR_STOP( MSG, 'ARE_FLAGS_VALID ("input_adj_mod.f")' )
      ! check not more than one:
      ENDIF
      IF ( ( LSENS  .AND. L4DVAR ) .or.
     &     ( LSENS  .AND. L3DVAR ) .or.
     &     ( L4DVAR .AND. L3DVAR )      ) THEN
         CALL ERROR_STOP( 'Either sensitivity or a var, pick only one!',
     &        ' ARE_FLAGS_VALID, input_adj_mod.f ' )
      ENDIF


      !=================================================================
      ! Check adjoint simulation subtypes
      !=================================================================
#if defined ( PM_ATTAINMENT ) || defined ( SOMO35_ATTAINMENT )
      IF ( OBS_FREQ /= 60 ) THEN
         CALL ERROR_STOP( ' OBS_FREQ should be 60 for attainment ',
     &                    ' ARE_FLAGS_VALID, input_adj_mod.f ' )
      ENDIF
#endif

      !If none of the datasets are selected or  PSEUDO_OBS FLAG, then it should be
      ! 3DVAR and 4DVAR
      ! add IMPROVE_BC_OC_OBS (adj32_013), MODIS_AOD_OBS (adj32_011)
      ! add OMI_SO2_OBS ()
      ! add MOPITT_V5_CO_OBS (adj32_016)
      ! add CH4 (kjw, dkh, 02/12/12, adj32_023)
      IF ( L3DVAR .or. L4DVAR  ) THEN
#if !defined(MOPITT_V5_CO_OBS) && !defined(MOPITT_V6_CO_OBS) && !defined(AIRS_CO_OBS) && !defined(SCIA_BRE_CO_OBS) && !defined(TES_NH3_OBS)&& !defined(SCIA_DAL_SO2_OBS) && !defined(PM_ATTAINMENT) && !defined(IMPROVE_SO4_NIT_OBS) && !defined(CASTNET_NH4_OBS) && !defined(SOMO35_ATTAINMENT) && !defined(TES_O3_OBS)&& !defined(SCIA_KNMI_NO2_OBS) && !defined(SCIA_DAL_NO2_OBS) && !defined(PSEUDO_OBS) && !defined(GOSAT_CO2_OBS) & !defined(INSITU_CO2_OBS) & !defined(MODIS_AOD_OBS) && !defined(IMPROVE_BC_OC_OBS) && !defined(TES_CH4_OBS) && !defined(SCIA_CH4_OBS) && !defined(MEM_CH4_OBS) && !defined(LEO_CH4_OBS) && !defined(GEOCAPE_CH4_OBS)  && !defined( TES_O3_IRK ) && !defined( OMI_SO2_OBS ) && !defined( OMI_NO2_OBS )
         MSG = 'Invalid adj run options: need to define obs for xDVAR'
         CALL ERROR_STOP( MSG, 'ARE_FLAGS_VALID ("input_adj_mod.f")' )
#endif
      ENDIF

      ! Conversly, if any of the obs operators are defined, then make sure it is
      ! a 3DVAR or 4DVAR simulation
      ! add IMPROVE_BC_OC_OBS (adj32_013), MODIS_AOD_OBS (adj32_011)
      ! add OMI_SO2_OBS ()
      ! add MOPITT_V5_CO_OBS (adj32_016)
      ! add CH4 (kjw, dkh, 02/12/12, adj32_023)
#if defined(MOPITT_V5_CO_OBS) || defined(MOPITT_V6_CO_OBS) || defined(AIRS_CO_OBS) || defined(SCIA_BRE_CO_OBS) || defined(TES_NH3_OBS) || defined(SCIA_DAL_SO2_OBS) || defined(PM_ATTAINMENT) || defined(IMPROVE_SO4_NIT_OBS) || defined(CASTNET_NH4_OBS) || defined(SOMO35_ATTAINMENT) || defined(TES_O3_OBS) || defined(SCIA_KNMI_NO2_OBS) || defined(SCIA_DAL_NO2_OBS) || defined(GOSAT_CO2_OBS) || defined(INSITU_CO2_OBS) || defined(MODIS_AOD_OBS) || defined(IMPROVE_BC_OC_OBS) || defined(TES_CH4_OBS) || defined(SCIA_CH4_OBS) || defined(MEM_CH4_OBS) || defined(LEO_CH4_OBS) || defined(GEOCAPE_CH4_OBS)  || defined(TES_O3_IRK) || defined(OMI_SO2_OBS) || defined( OMI_NO2_OBS ) || defined(PSEUDO_OBS)

         IF ( .not. ( L3DVAR .or. L4DVAR ) ) THEN
            MSG = 'Invalid adj run options: need to define VAR for obs'
            CALL ERROR_STOP( MSG, 'ARE_FLAGS_VALID ("input_adj_mod.f")')
         ENDIF
#endif

      IF ( LSAT_HDF_L2 .or. LSAT_HDF_L3  ) THEN
#if !defined(MOPITT_V5_CO_OBS) && !defined(MOPITT_V6_CO_OBS) && !defined( OMI_SO2_OBS ) && !defined( OMI_NO2_OBS )
         MSG = 'Invalid adj run options: Satellite HDF diagnostics are
     & only supported by OMI, TES and MOPITT obs operator for xDVAR'
         CALL ERROR_STOP( MSG, 'ARE_FLAGS_VALID ("input_adj_mod.f")' )
#endif
      ENDIF

      ! If we are using real observations, make sure pseudo obs are commented (mak, dkh, 10/01/09)
      ! add IMPROVE_BC_OC_OBS (adj32_013), MODIS_AOD_OBS (adj32_011)
      ! add OMI_SO2_OBS ()
      ! add MOPITT_V5_CO_OBS (adj32_016)
      ! add CH4 (kjw, dkh, 02/12/12, adj32_023)
#if  defined(MOPITT_V5_CO_OBS) || defined(MOPITT_V6_CO_OBS) || defined(AIRS_CO_OBS) || defined(SCIA_BRE_CO_OBS) || defined(TES_NH3_OBS) || defined(SCIA_DAL_SO2_OBS) || defined(PM_ATTAINMENT) || defined(IMPROVE_SO4_NIT_OBS) || defined(CASTNET_NH4_OBS) || defined(SOMO35_ATTAINMENT) || defined(TES_O3_OBS) || defined(SCIA_KNMI_NO2_OBS) || defined(SCIA_DAL_NO2_OBS) || defined(GOSAT_CO2_OBS) || defined(INSITU_CO2_OBS) || defined(MODIS_AOD_OBS) || defined(IMPROVE_BC_OC_OBS) || defined(TES_CH4_OBS) || defined(SCIA_CH4_OBS) || defined(MEM_CH4_OBS) || defined(LEO_CH4_OBS) || defined(GEOCAPE_CH4_OBS) || defined(TES_O3_IRK) || defined(OMI_SO2_OBS) || defined( OMI_NO2_OBS )


#if defined(PSEUDO_OBS)
         IF ( L4DVAR  ) THEN
            MSG = 'Invalid adj options: define real or pseudo obs'
            CALL ERROR_STOP( MSG, 'ARE_FLAGS_VALID ("input_adj_mod.f")')
         ENDIF
#endif

#endif
      ! ( LFDTEST .AND. .NOT. LSENS ) LSENS = .TRUE.
      IF ( LFDTEST .AND. (.not. LSENS ) ) THEN
          CALL ERROR_STOP( 'FD tests are a subtpye of SENS',
     &        ' ARE_FLAGS_VALID, input_adj_mod.f ' )
       ENDIF


      IF ( LFDTEST .AND. LICS .AND. LADJ_EMS ) THEN
          CALL ERROR_STOP( 'FD test for ems AND ics not supported',
     &        ' ARE_FLAGS_VALID, input_adj_mod.f ' )
       ENDIF

      IF ( LFDTEST .and.
     &     ( ( N_CALC_STOP > 3 ) .or.
     &       ( N_CALC_STOP < 1 )      ) ) THEN
          CALL ERROR_STOP( 'FD tests need to have 1 < N_CALC_STOP < 3',
     &        ' ARE_FLAGS_VALID, input_adj_mod.f ' )
      ENDIF

      IF ( LFDTEST .AND. LFD_GLOB .AND. LTRAN ) THEN
          CALL ERROR_STOP( 'FD_GLOB should be done with transport off',
     &        ' ARE_FLAGS_VALID, input_adj_mod.f ' )
       ENDIF

      ! Estimating inv Hessian only supported for 4DVar (dkh, 01/12/12, adj32_012)
      IF ( ( LINVH .or. LINVH_BFGS ) .and. ( .not. L4DVAR ) ) THEN
         CALL ERROR_STOP( 'LINVH and LINVH_BFGS only with 4DVAR ',
     &        ' ARE_FLAGS_VALID, input_adj_mod.f ' )
      ENDIF

      ! RTFD should equal the "Rate #" listed in input.gcadj (or RF_IDX)
      ! corresponding to the listed rate we wish to test
      IF ( LFDTEST .and. LADJ_RRATE ) THEN
         IF ( RATFD > NRRATES ) THEN
            CALL ERROR_STOP('Invalid RTFD',
     &        ' ARE_FLAGS_VALID, input_adj_mod.f ' )
         ENDIF
      ENDIF

      IF ( .NOT. LADJ_RRATE .AND. ( NCOEFF_RATE .NE. 0 ) )
     &   CALL ERROR_STOP('Invalid NCOEFF_RATE', 'ARE_FLAGS_VALID,
     &    input_adj_mod.f ')

      !=================================================================
      ! Check adjoint control parameters
      !=================================================================
      IF ( (.not. LICS ) .AND. ( .not. LADJ_EMS ) ) THEN
          CALL ERROR_STOP( 'Must select either ICS or EMS ',
     &        ' ARE_FLAGS_VALID, input_adj_mod.f ' )
       ENDIF

      ! LADJ_STRAT is a sub-type of LADJ_EMS (dkh, 02/23/12, adj32_025)
      IF ( ( LADJ_STRAT ) .AND. ( .not. LADJ_EMS ) ) THEN
          CALL ERROR_STOP( 'LADJ_STRAT is a sub-type of LADJ_EMS',
     &        ' ARE_FLAGS_VALID, input_adj_mod.f ' )
       ENDIF

      ! check settings for tagged Ox sim
      IF ( ITS_A_TAGOX_SIM() ) THEN
          IF (  LICS ) THEN
             CALL ERROR_STOP( 'Tagged OX adjoint only LADJ_EMS ',
     &           ' ARE_FLAGS_VALID, input_adj_mod.f ' )
          ENDIF
!          IF (  MMSCL .ne. LLPAR ) THEN
!             CALL ERROR_STOP( 'Need MMSCL = LLPAR for tag ox adj ' ,
!     &           ' ARE_FLAGS_VALID, input_adj_mod.f ' )
!          ENDIF
      ENDIF

      IF ( IFD .GT. IIPAR ) THEN
         CALL ERROR_STOP( ' IFD has to be less than IIPAR !',
     &        ' ARE_FLAGS_VALID, input_adj_mod.f ' )
      ENDIF

      IF ( JFD .GT. JJPAR ) THEN
         CALL ERROR_STOP( ' JFD has to be less than JJPAR !',
     &        ' ARE_FLAGS_VALID, input_adj_mod.f ' )
      ENDIF

      IF ( LFD .GT. LLPAR ) THEN
         CALL ERROR_STOP( ' LFD has to be less than LLPAR !',
     &        ' ARE_FLAGS_VALID, input_adj_mod.f ' )
      ENDIF

      IF ( NFD .GT. N_TRACERS ) THEN
         CALL ERROR_STOP( ' NFD has to be less than number of tracers!',
     &        ' ARE_FLAGS_VALID, input_adj_mod.f ' )
      ENDIF

      IF ( ICSFD .GT. N_TRACERS ) THEN
         CALL ERROR_STOP( ' ICSFD has to be < number of tracers!',
     &        ' ARE_FLAGS_VALID, input_adj_mod.f ' )
      ENDIF

      ! (dkh, 11/11/09)
      IF ( LADJ_EMS ) THEN
         IF ( EMSFD .GT. NNEMS ) THEN
            CALL ERROR_STOP(
     &           ' EMSFD has to be < number of active adj emissons!',
     &           ' ARE_FLAGS_VALID, input_adj_mod.f ' )
         ENDIF
      ENDIF

      ! (dkh, 01/12/12, adj32_012)
!      IF ( LINVH .and. ( .not. LADJ_EMS .or. LICS ) ) THEN
!         CALL ERROR_STOP( ' LINVH only supported for LADJ_EMS ',
!     &                    ' ARE_FLAGS_VALID, input_adj_mod.f ' )
!      ENDIF
!      IF ( ( LINVH .or. LINVH_BFGS ) .and.
!     &          ( .not. LADJ_EMS .or. LICS ) ) THEN
!         CALL ERROR_STOP( ' LINVH only supported for LADJ_EMS ',
!     &                    ' ARE_FLAGS_VALID, input_adj_mod.f ' )
!      ENDIF

      ! Check to make sure error specifications are usable for LAPSRC (dkh, 02/22/11)
      IF ( LAPSRC ) THEN

         ! Check emissions
         IF ( LADJ_EMS ) THEN

            DO N = 1, NNEMS

               ! Skip emissions that are not included in optimization
               IF ( .not. OPT_THIS_EMS(N) ) CYCLE

#if defined ( LOG_OPT )
               IF ( EMS_ERROR(N) < ( 1d0 + SMAL2 ) ) THEN
                  print*, ' EMS_ERROR stop at N = ', N
                  CALL ERROR_STOP( ' EMS_ERROR is too small ',
     &                             ' input_adj_mod.f '         )
               ENDIF
#else
               IF ( EMS_ERROR(N) < ( SMAL2 ) ) THEN
                  print*, ' EMS_ERROR stop at N = ', N
                  CALL ERROR_STOP( ' EMS_ERROR is too small ',
     &                             ' input_adj_mod.f '         )
               ENDIF
#endif
            ENDDO

            ! Check strat prod and loss tracers (hml, adj32_025)
            IF  ( LADJ_STRAT ) THEN

               DO N = 1, NSTPL

                  ! Skip tracers that are not included in optimization
                  IF (.not. OPT_THIS_PROD(N) .AND.
     &                .not. OPT_THIS_LOSS(N)) CYCLE

                  IF ( PROD_ERROR(N) < ( SMAL2 ) ) THEN
                     print*, ' PROD_ERROR stop at N = ', N
                     CALL ERROR_STOP( ' PROD_ERROR is too small ',
     &                                ' input_adj_mod.f '         )
                  ENDIF

                  IF ( LOSS_ERROR(N) < ( SMAL2 ) ) THEN
                     print*, ' LOSS_ERROR stop at N = ', N
                     CALL ERROR_STOP( ' LOSS_ERROR is too small ',
     &                                ' input_adj_mod.f '         )
                  ENDIF
               ENDDO
            ENDIF

         ! Check tracers
         ELSEIF  ( LICS ) THEN

            DO N = 1, N_TRACERS

               ! Skip tracers that are not included in optimization
               IF ( .not. OPT_THIS_TRACER(N) ) CYCLE

#if defined ( LOG_OPT )
               IF ( ICS_ERROR(N) < ( 1d0 + SMAL2 ) ) THEN
                  print*, ' ICS_ERROR stop at N = ', N
                  CALL ERROR_STOP( ' ICS_ERROR is too small ',
     &                             ' input_adj_mod.f '         )
               ENDIF
#else
               IF ( ICS_ERROR(N) < ( SMAL2 ) ) THEN
                  print*, ' ICS_ERROR stop at N = ', N
                  CALL ERROR_STOP( ' ICS_ERROR is too small ',
     &                             ' input_adj_mod.f '         )
               ENDIF
#endif
            ENDDO
         ENDIF
      ENDIF

      !=================================================================
      ! Check observation settings
      !=================================================================
#if   defined ( SCIA_KNMI_NO2_OBS ) || defined ( SCIA_DAL_NO2_OBS )
      ! Since the NO2 obs operators will pass adjoints back
      ! to CSPEC via CSPEC_AFTER_CHEM_ADJ, we need to make sure that
      ! these species are listed as observed species
      FOUND = .FALSE.
      DO N = 1, NOBS_CSPEC

         IF ( TRIM( NAMEGAS( IDCSPEC_ADJ(N) ) ) == 'NO2' ) THEN
            FOUND = .TRUE.
         ENDIF

      ENDDO
      IF ( .not. FOUND ) THEN

         CALL ERROR_STOP( ' Need to list NO2 as observed species',
     &                    ' input_adj_mod ' )
      ENDIF

! BUG FIX: move this to INIT_CSPEC_ADJ, by which point the necessary
! CSPEC variables have been initialized (nb, dkh, 01/06/12, adj32_002)
!--------------------------------------------------------------------
!#elif defined ( TES_O3_OBS )
!      ! Since the O3 obs operators will pass adjoints back
!      ! to CSPEC via CSPEC_AFTER_CHEM_ADJ, we need to make sure that
!      ! these species are listed as observed species
!      FOUND = .FALSE.
!      DO N = 1, NOBS_CSPEC
!
!         IF ( TRIM( NAMEGAS( IDCSPEC_ADJ(N) ) ) == 'O3' ) THEN
!            FOUND = .TRUE.
!         ENDIF
!
!      ENDDO
!      IF ( .not. FOUND ) THEN
!
!         CALL ERROR_STOP( ' Need to list O3 as observed species',
!     &                    ' input_adj_mod ' )
!      ENDIF
!--------------------------------------------------------------------
#endif

      ! We only observe species in CSPEC for full chemistry runs
      IF ( .not. ITS_A_FULLCHEM_SIM() .and.
     &           NOBS_CSPEC /= 0            ) THEN
         CALL ERROR_STOP( ' NOBS_CSPEC needs to be zero',
     &                    ' input_adj_mod ' )
      ENDIF

      ! If we are using CSPEC for the cost function, then
      ! at least one species needs to be listed in the obsevation
      ! menu.
      IF ( LCSPEC_PPB .and. ( .not. LCSPEC_OBS ) ) THEN
         CALL ERROR_STOP(
     &      ' Need to observe a cspec species for LCSPEC_PPB',
     &                    ' input_adj_mod ' )
      ENDIF

      ! If we are doing a sensitivty calculation w.r.t. cspec
      ! observations, then make sure we have the cspec-based
      ! option selected.
      IF ( LSENS .and. LCSPEC_OBS .and. ( .not. LCSPEC_PPB )
     &     .and. ( .not. LADJ_DDEP_CSPEC ) ) THEN
         CALL ERROR_STOP(
     &      ' Need to select a cost function option that uses CSPEC',
     &                    ' input_adj_mod ' )
      ENDIF
#if   defined ( PSEUDO_OBS )
      IF ( LCSPEC_OBS ) THEN
         CALL ERROR_STOP(
     &      ' PSEUDO_OBS only implemented for tracer obs',
     &                    ' input_adj_mod ' )
      ENDIF
#endif

!      ! The deposition forcings are cummulative, and the coding
!      ! of the timing of the forcing assumes LMAX_OBS
!      IF ( LADJ_FDEP .and. ( .not. LMAX_OBS ) ) THEN
!         CALL ERROR_STOP (' Need LMAX_OBS = T and NSPAN for LADJ_FDEP',
!     &                    ' input_adj_mod ' )
!      ENDIF

      ! Deposition forcing FD tests use forward model diagnostics
      ! for evaluation of depo fluxes
      IF ( LFD_GLOB .and.  ( LADJ_DDEP_TRACER .or. LADJ_DDEP_CSPEC )
     &   .and. ( ND44 == 0 ) ) THEN
         CALL ERROR_STOP (' Turn on ND44 for dry dep forcing FD test ',
     &                    ' input_adj_mod ' )
      ENDIF

      ! Only allow one dep forcing option at a time for FD tests
      IF ( LFDTEST .and. LADJ_FDEP ) THEN
         COUNT_ON = 0
         IF ( LADJ_WDEP_LS     ) COUNT_ON = COUNT_ON + 1
         IF ( LADJ_WDEP_CV     ) COUNT_ON = COUNT_ON + 1
         IF ( LADJ_DDEP_TRACER ) COUNT_ON = COUNT_ON + 1
         IF ( LADJ_DDEP_CSPEC  ) COUNT_ON = COUNT_ON + 1
         IF ( COUNT_ON > 1 ) THEN
            CALL ERROR_STOP (' Only one dep forcing for FD test ',
     &                       ' input_adj_mod ' )
         ENDIF
         IF ( COUNT_ON == 0 ) THEN
            CALL ERROR_STOP (' Which dep forcing option do you want?',
     &                       ' input_adj_mod ' )
         ENDIF
      ENDIF

      ! Deposition forcing uses forward model diagnostics
      ! for evaluation of depo fluxes
      IF ( LADJ_WDEP_LS .and. ( ND39 == 0 ) ) THEN
         CALL ERROR_STOP (' Turn on ND39 for wet LS forcing ',
     &                    ' input_adj_mod ' )
      ENDIF

      ! Deposition forcing uses forward model diagnostics
      ! for evaluation of depo fluxes
      IF ( LADJ_WDEP_CV .and. ( ND38 == 0 ) ) THEN
         CALL ERROR_STOP (' Turn on ND38 for wet CV forcing ',
     &                    ' input_adj_mod ' )
      ENDIF

      ! Deposition forcing uses forward model diagnostics
      ! for evaluation of depo fluxes
      IF ( ( LADJ_DDEP_TRACER .or. LADJ_DDEP_CSPEC ) .and.
     &   ( ND44 == 0 ) ) THEN
         CALL ERROR_STOP (' Turn on ND44 for DDEP forcing ',
     &                    ' input_adj_mod ' )
      ENDIF

      ! FD test of the dry deposition adjoint only supported for molec/cm2/s
      IF ( LFD_GLOB .and.  ( LADJ_DDEP_TRACER .or. LADJ_DDEP_CSPEC )
     &               .and. ( .not. LMOLECCM2S ) ) THEN
         CALL ERROR_STOP (' Set units to molec/cm2/s for ddep FD test',
     &                    ' input_adj_mod ' )
      ENDIF

      ! FD test of the wet deposition adjoint only supported for kg/s
      IF ( LFD_GLOB .and.  ( LADJ_WDEP_LS .or. LADJ_WDEP_CV )
     &               .and. ( .not. LKGS ) ) THEN
         CALL ERROR_STOP (' Set units to kg/s for wdep FD test',
     &                    ' input_adj_mod ' )
      ENDIF

      ! Make sure that NFD matches the observed tracer or species
      IF ( LFDTEST ) THEN

         ! check species
         IF ( NOBS_CSPEC > 0 ) THEN

            IF ( NFD /= 1 .or. NOBS_CSPEC > 1 ) THEN
               CALL ERROR_STOP(
     &         ' For species FD, list only one species and set NFD = 1',
     &                         ' input_adj_mod' )
            ENDIF

         ! check tracers
         ELSE

            IF ( .not. OBS_THIS_TRACER(NFD) ) THEN
               CALL ERROR_STOP(' Observed tracer and NFD must match',
     &                         ' input_adj_mod' )
            ENDIF

            IF ( NOBS > 1 ) THEN
               CALL ERROR_STOP(' Only observe tracer NFD for FD test',
     &                         ' input_adj_mod' )
            ENDIF
         ENDIF

      ENDIF

      ! Check to make sure that our observation time range fits in the
      ! simulation time range
      IF ( LMAX_OBS ) THEN
         DATE =  GET_TIME_AHEAD( NSPAN * OBS_FREQ )
         print*, ' DDD DATE = ', DATE
         print*, ' DDD NYMDe= ', GET_NYMDe()
         print*, ' DDD NHMSe= ', GET_NHMSe()
         IF ( ( DATE(1) > GET_NYMDe() ) .or.
     &        ( DATE(1) == GET_NYMDe().and.
     &          DATE(2) > GET_NHMSe()      ) ) THEN
             CALL ERROR_STOP(' NSPAN too long! ',
     &                       ' input_adj_mod' )
         ENDIF
      ENDIF

      !=================================================================
      ! Check diagnostics
      !=================================================================
      IF ( LEMS_ABS .and. ( .not. LADJ_EMS ) ) THEN
         CALL ERROR_STOP (' LEMS_ABS only for active vars = emissions',
     &                    ' input_adj_mod ' )
      ENDIF

      !=================================================================
      ! Check if all emissions adjoint ID #'s are defined for particular
      ! sets of emissions species.
      !=================================================================

      ! Primary carbonaceous aerosol emissions
      IF ( IDADJ_EBCPI_an > 0 .and. IDADJ_EBCPO_an > 0 .and.
     &     IDADJ_EOCPI_an > 0 .and. IDADJ_EOCPO_an > 0 .and.
     &     IDADJ_EBCPI_bb > 0 .and. IDADJ_EBCPO_bb > 0 .and.
     &     IDADJ_EOCPI_bb > 0 .and. IDADJ_EOCPO_bb > 0 .and.
     &     IDADJ_EBCPI_bf > 0 .and. IDADJ_EBCPO_bf > 0 .and.
     &     IDADJ_EOCPI_bf > 0 .and. IDADJ_EOCPO_bf > 0       ) THEN
         IS_CARB_EMS_ADJ = .TRUE.
      ENDIF
      IF ( N_CARB_EMS_ADJ > 0 .and. ( .not. IS_CARB_EMS_ADJ ) ) THEN
         CALL ERROR_STOP( 'Not enough carbon emissions adjoint IDs ',
     &                    'ARE_FLAGS_VALID')
      ENDIF
      IF ( N_CARB_EMS_ADJ > 12 ) THEN
         CALL ERROR_STOP( 'Too many carbon emissions adjoint IDs ',
     &                    'ARE_FLAGS_VALID')
      ENDIF

      ! Sulfate aerosol (and precursor) emissions
      IF ( IDADJ_ENH3_bb > 0 .and. IDADJ_ENH3_bf  > 0 .and.
     &     IDADJ_ENH3_na > 0 .and. IDADJ_ENH3_an  > 0 .and.
     &     IDADJ_ESO2_bb > 0 .and. IDADJ_ESO2_an1 > 0 .and.
     &     IDADJ_ESO2_bf > 0 .and. IDADJ_ESO2_an2 > 0 .and.
     &     IDADJ_ESO2_sh > 0                                ) THEN
         IS_SULF_EMS_ADJ = .TRUE.
      ENDIF
      IF ( N_SULF_EMS_ADJ > 0 .and. ( .not. IS_SULF_EMS_ADJ ) ) THEN
         CALL ERROR_STOP(
     &            'Not enough sulfate aerosol emissions adjoint IDs ',
     &            'ARE_FLAGS_VALID')
      ENDIF
      IF ( N_SULF_EMS_ADJ > 9 ) THEN
         CALL ERROR_STOP( 'Too many sulfate emissions adjoint IDs ',
     &                    'ARE_FLAGS_VALID')
      ENDIF

      ! Dust aerosol emissions ( xxu, 11/01/10) (dkh, 01/09/12, adj32_011)
      IF ( IDADJ_EDST1 > 0 .and. IDADJ_EDST2  > 0 .and.
     &     IDADJ_EDST3 > 0 .and. IDADJ_EDST4  > 0  ) THEN
         IS_DUST_EMS_ADJ = .TRUE.
      ENDIF
      IF ( N_DUST_EMS_ADJ > 0 .and. ( .not. IS_DUST_EMS_ADJ ) ) THEN
         CALL ERROR_STOP(
     &            'Not enough Dust aerosol emissions adjoint IDs ',
     &            'ARE_FLAGS_VALID')
      ENDIF
      IF ( N_DUST_EMS_ADJ > 4 ) THEN
         CALL ERROR_STOP( 'Too many dust emissions adjoint IDs ',
     &                    'ARE_FLAGS_VALID')
      ENDIF

      !=================================================================
      ! Check consistency between input.gcadj and define_adj.h options
      !=================================================================

      IF ( LBKCOV ) THEN

      IF ( LICS ) CALL ERROR_STOP( 'Off-diagonal calculation only
     & works with LADJ_EMDS', 'ARE_FLAGS_VALID' )

#if ! defined ( LBKCOV_ERR )

      CALL ERROR_STOP( 'Off-diagonal calculation requires LBKCOV_ERR
     & to be set in define_adj.h ', 'ARE_FLAGS_VALID' )

#endif
      ENDIF

#if  defined ( LBKCOV_ERR )
      IF ( .not. LBKCOV ) THEN

      CALL ERROR_STOP( 'Off-diagonal calculation requires LBKCOV
     & to be set in input.gcadj ', 'ARE_FLAGS_VALID' )

      ENDIF
#endif


      IF ( LINVH_BFGS ) THEN
#if ! defined ( LBFGS_INV )

      CALL ERROR_STOP( 'L-BFGS calculation requires LBFGS_INV
     & to be set in define_adj.h ', 'ARE_FLAGS_VALID' )

#endif
      ENDIF

#if  defined ( LBFGS_INV )
      IF ( .not. LINVH_BFGS ) THEN

      CALL ERROR_STOP( 'L-BFGS calculation requires an option
     & to be set in input.gcadj ', 'ARE_FLAGS_VALID' )

      ENDIF
#endif

      ! fp check for wetdep sensitivities: these units only make sense
      ! if we observe one tracer or species at a time.
      IF (   ( LKGS                .OR.      LMOLECCM2S ) .AND.
     &     ( ( NOBS_CSPEC .GT. 1 ) .OR. ( NOBS .GT. 1 ) ) ) THEN

         CALL ERROR_STOP(' Does not seem to make sense ?',
     &                   'ARE_FLAGS_VALID')
      ENDIF

      !throw an error if kks and nobs_cspec>1
      IF ( LKGS .AND. ( NOBS_CSPEC .GT. 1) ) THEN
         CALL ERROR_STOP(' not implemented',
     &                   'ARE_FLAGS_VALID')
      ENDIF


      IF ( LEQHAYR .OR. LKGNHAYR ) THEN

         DO T = 1, NOBS_CSPEC

            FOUND = .FALSE.

            DO N = 1, N_NDEP_CSPEC

               IF ( TRIM(CNAME(T))
     &              .NE. TRIM(NDEP_CSPEC(N)) )
     &              FOUND = .TRUE.

            ENDDO

            IF ( .not. FOUND ) THEN
               WRITE(6,'( a )') CNAME(T)
               CALL ERROR_STOP(' Does not seem to make sense ?',
     &              'ARE_FLAGS_VALID')
            ENDIF

         ENDDO

         DO T = 1, N_TRACERS

            IF ( OBS_THIS_TRACER(T) ) THEN

               FOUND = .FALSE.

               DO N = 1,N_NDEP

                  IF ( T .EQ. NDEP( N ) ) FOUND = .TRUE.

               ENDDO

               IF ( LEQHAYR ) THEN

                  DO N = 1,N_ACID

                     IF ( T .EQ. ACID(N) ) FOUND = .TRUE.

                  ENDDO

               ENDIF

               IF ( .not. FOUND ) THEN
                  WRITE(6,*) 'TRACER: ',TRACER_NAME(T)
                  CALL ERROR_STOP(' Does not seem to make sense ?',
     &                 'ARE_FLAGS_VALID')
               ENDIF

            ENDIF
         ENDDO

      ENDIF

      IF ( LADJ_CL ) THEN

         IF ( .not. LADJ_CL_NDEP .and.
     &        .not. LADJ_CL_ACID ) THEN

            CALL ERROR_STOP(' Need to select N/Ac for Cl cost function',
     &           'ARE_FLAGS_VALID')

         ENDIF

         IF ( LADJ_CL_NDEP .and. .not. LKGNHAYR )
     &        CALL ERROR_STOP(' Units are inconsistent',
     &           'ARE_FLAGS_VALID')


         IF ( LADJ_CL_ACID .and. .not. LEQHAYR )
     &        CALL ERROR_STOP(' Units are inconsistent',
     &           'ARE_FLAGS_VALID')


         IF ( .not. LADJ_DDEP_TRACER .OR.
     &        .not. LADJ_DDEP_CSPEC  .OR.
     &        .not. LADJ_WDEP_LS     .OR.
     &        .not. LADJ_WDEP_CV     ) THEN

            CALL ERROR_STOP(
     &           ' All deposition flags need to be turned on',
     &           'ARE_FLAGS_VALID')

         ENDIF

         DO T = 1, N_NDEP

            IF ( .not. OBS_THIS_TRACER( NDEP(T) ) ) THEN

               WRITE(*,*) 'TRACER: ',TRACER_NAME(NDEP(T))
               CALL ERROR_STOP(
     &              'All N/Acid species need to be listed',
     &              'ARE_FLAGS_VALID')

            ENDIF

         ENDDO

         IF ( LADJ_CL_ACID ) THEN

            DO T = 1, N_ACID

               IF ( .not. OBS_THIS_TRACER( ACID(T) ) )  THEN

                  WRITE(*,*) 'TRACER: ',TRACER_NAME(ACID(T))

                  CALL ERROR_STOP(
     &                 ' All N/Acid species need to be listed',
     &                 'ARE_FLAGS_VALID')
               ENDIF

            ENDDO

         ENDIF

         DO N = 1, N_NDEP_CSPEC

            FOUND = .FALSE.

            DO T = 1, NOBS_CSPEC

               IF ( TRIM(CNAME(T))
     &              .NE. TRIM(NDEP_CSPEC(N)) )
     &              FOUND = .TRUE.

            ENDDO

            IF ( .not. FOUND ) THEN

               WRITE(*,*) 'CSPEC: ',TRIM(NDEP_CSPEC(N))
               CALL ERROR_STOP(
     &              ' All N/Acid species need to be listed',
     &              'ARE_FLAGS_VALID')
            ENDIF

         ENDDO


      ENDIF


      ! Return to calling program
      END SUBROUTINE ARE_FLAGS_VALID

!------------------------------------------------------------------------------

      SUBROUTINE VALIDATE_DIRECTORIES
!
!******************************************************************************
!  Subroutine VALIDATE_DIRECTORIES makes sure that each of the directories
!  that we have read from the GEOS-CHEM input file are valid.  Also, trailing
!  separator characters will be added. (bmy, 7/20/04, 8/4/06)
!
!  NOTES:
!  (1 ) Now make sure all USE statements are USE, ONLY.  Now also validate
!        GCAP and GEOS-5 directories. (bmy, 10/3/05)
!  (2 ) Now references DATA_DIR_1x1 from directory_mod.f (bmy, 10/24/05)
!  (3 ) Remove support for GEOS-1 and GEOS-STRAT met fields (bmy, 8/4/06)
!******************************************************************************
!
      ! References to F90 modules
      ! References to F90 modules
      USE DIRECTORY_ADJ_MOD, ONLY : OPTDATA_DIR, ADJTMP_DIR, DIAGADJ_DIR

      ! Local variables
      CHARACTER(LEN=255)     :: DIR

      !=================================================================
      ! VALIDATE_DIRECTORIES begins here!
      !=================================================================

      ! Check directories
      CALL CHECK_DIRECTORY( OPTDATA_DIR     )
      CALL CHECK_DIRECTORY( ADJTMP_DIR      )
      CALL CHECK_DIRECTORY( DIAGADJ_DIR     )

      ! Return to calling program
      END SUBROUTINE VALIDATE_DIRECTORIES

!------------------------------------------------------------------------------

      SUBROUTINE CHECK_DIRECTORY( DIR )
!
!******************************************************************************
!  Subroutine CHECK_DIRECTORY makes sure that the given directory
!  is valid.  Also a trailing slash character will be added if necessary.
!  (bmy, 3/20/03, 3/23/05)
!
!  Arguments as Input/Output:
!  ============================================================================
!  (1 ) DIR (CHARACTER) : Directory to be checked
!
!  NOTES:
!  (1 ) Now references FILE_EXISTS from "file_mod.f" (bmy, 3/23/05)
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD,     ONLY : ERROR_STOP
      USE FILE_MOD,      ONLY : FILE_EXISTS
      USE UNIX_CMDS_MOD, ONLY : SEPARATOR

      ! Arguments
      CHARACTER(LEN=*), INTENT(INOUT) :: DIR

      ! Local variables
      INTEGER                         :: C
      CHARACTER(LEN=255)              :: MSG

      !=================================================================
      ! CHECK_DIRECTORY begins here!
      !=================================================================

      ! Locate the last non-white-space character of NEWDIR
      C = LEN_TRIM( DIR )

      ! Add the trailing directory separator if it is not present
      IF ( DIR(C:C) /= TRIM( SEPARATOR ) ) THEN
         DIR(C+1:C+1) = TRIM( SEPARATOR )
      ENDIF

      !=================================================================
      ! Test if the directory actually exists
      !=================================================================

      ! If the directory does not exist then stop w/ an error message
      IF ( .not. FILE_EXISTS( DIR ) ) THEN
         MSG = 'Invalid directory: ' // TRIM( DIR )
         CALL ERROR_STOP( MSG, 'CHECK_DIRECTORY ("input_adj_mod.f")' )
      ENDIF

      ! Return to calling program
      END SUBROUTINE CHECK_DIRECTORY

!------------------------------------------------------------------------------

      SUBROUTINE CHECK_FILE( FILE )
!
!******************************************************************************
!  Subroutine CHECK_FILE makes sure that the given file exists. (dkh, 03/10/13)
!  Based on CHECK_DIR (bmy, 3/20/03, 3/23/05)
!
!  Arguments as Input/Output:
!  ============================================================================
!  (1 ) FILE (CHARACTER) : File to be checked
!
!  NOTES:
!
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD,     ONLY : ERROR_STOP
      USE FILE_MOD,      ONLY : FILE_EXISTS

      ! Arguments
      CHARACTER(LEN=*), INTENT(IN) :: FILE

      ! Local variables
      CHARACTER(LEN=255)              :: MSG

      !=================================================================
      ! CHECK_FILE begins here!
      !=================================================================

      ! If the directory does not exist then stop w/ an error message
      IF ( .not. FILE_EXISTS( FILE ) ) THEN
         MSG = 'Invalid file: ' // TRIM( FILE )
         CALL ERROR_STOP( MSG, 'CHECK_FILE ("input_adj_mod.f")' )
      ENDIF

      ! Return to calling program
      END SUBROUTINE CHECK_FILE

!------------------------------------------------------------------------------

      SUBROUTINE CLEAN_FILE_DIRS()
!
!******************************************************************************
! Subroutine CLEAN_FILE_DIRS gets rid of files in ADJTMP_DIR and in OptData that
! are left over from previous runs. (10/28/04)
!
!
!  NOTES:
!  (1 ) If the last run to be computed completed cleanly, there will not be
!       any *.chk.* files, and SYSTEM will complain a bit about this.  It's OK
!       (dkh, 10/03/04)
!  (2 ) Add caviot that if L_MAKE_CHK is false, don't delete old *chk* files
!  (3 ) Add feature to clean out OPTDATA_DIR (dkh, 10/28/04)
!  (4 ) Delete *.ics.* and *.gdt.* files during observation run. (dkh, 11/11/04)
!  (5 ) Delete cfn.* files during observation run. (dkh, 02/13/06)
!  (6 ) Move from inverse_mod.f to input_adj_mod.f (dkh, 07/28/09)
!  (7 ) Now clean out old ems.adj.* and gctm.iteration files (dkh, 02/17/11)
!  (8 ) Now keep files for offline inv hessian  (dkh, 01/12/12, adj32_012)
!  (9 ) Add support for strat fluxes LADJ_STRAT (hml, dkh, 02/14/12, adj32_025)
!  (10 ) Add LINVH_BFGS (nab, 25/03/12 )
!******************************************************************************
!
      ! Reference to F90 modules
      USE ADJ_ARRAYS_MOD,    ONLY : N_CALC
      USE ADJ_ARRAYS_MOD,    ONLY : N_CALC_STOP
      USE DIRECTORY_ADJ_MOD, ONLY : OPTDATA_DIR
      USE DIRECTORY_ADJ_MOD, ONLY : ADJTMP_DIR
      USE DIRECTORY_ADJ_MOD, ONLY : DIAGADJ_DIR
      USE ERROR_MOD,         ONLY : ERROR_STOP
      USE LOGICAL_ADJ_MOD,   ONLY : LEMS_ABS
      USE LOGICAL_ADJ_MOD,   ONLY : LITR
      USE LOGICAL_ADJ_MOD,   ONLY : LINVH,LINVH_BFGS

#     include "CMN_SIZE"     ! Size params

      ! Local variables
      CHARACTER(LEN=255)    :: REMOVE_OBS_FILE_CMD
      CHARACTER(LEN=255)    :: REMOVE_CHK_FILE_CMD
      CHARACTER(LEN=255)    :: REMOVE_ADJ_FILE_CMD
      CHARACTER(LEN=255)    :: REMOVE_OPT_FILE_CMD
      CHARACTER(LEN=255)    :: REMOVE_FD_FILE_CMD

      !============================================================
      ! CLEAN_FILE_DIRS starts here!
      !============================================================

      WRITE( 6, '(a)'   ) REPEAT( '=', 79 )
      WRITE( 6, '(a,/)'   ) 'C L E A N   O U T   O L D   F I L E S'

      IF ( N_CALC_STOP == 0 ) THEN

         ! Clear any old .obs. files
         REMOVE_OBS_FILE_CMD = 'rm ' //
     &                         TRIM( ADJTMP_DIR ) // '*.obs.*'

         WRITE( 6, 102 ) TRIM( REMOVE_OBS_FILE_CMD )
 102     FORMAT( '     - INVERSE: Executing: ',a )

         CALL SYSTEM( TRIM ( REMOVE_OBS_FILE_CMD ) )

         ! Clean out old *.gdt.*, *.ics.* and cnf.* files
         REMOVE_OPT_FILE_CMD = 'rm ' //
     &                       TRIM (OPTDATA_DIR) // '*.gdt.*'
     &                       // ' ' //
     &                       TRIM (OPTDATA_DIR) // '*.sf.*'
     &                       // ' ' //
     &                       TRIM (OPTDATA_DIR) // 'cfn.*'

         WRITE( 6, 102 ) TRIM( REMOVE_OPT_FILE_CMD )

         CALL SYSTEM ( TRIM( REMOVE_OPT_FILE_CMD ) )


      ELSE


         ! Clean out old .chk. files
         REMOVE_CHK_FILE_CMD = 'rm ' //
     &                         TRIM (ADJTMP_DIR) // '*.chk.*'

         WRITE( 6, 102 ) TRIM( REMOVE_CHK_FILE_CMD )

         CALL SYSTEM ( TRIM( REMOVE_CHK_FILE_CMD ) )


         ! Clean out old .adj. files
         ! BUG FIX: the *.adj.* files are in DAIGADJ_DIR (jk, dkh, 04/25/10)
         ! Update: be more specific here so that we don't delete ems.adj.NN
         !  (dkh, 02/18/11)
         ! Now keep these if doing inv Hessian update (dkh, 01/12/12, adj32_012)
         IF ( .not. ( LINVH .or. LINVH_BFGS ) ) THEN

            REMOVE_ADJ_FILE_CMD = 'rm ' //
     &                         TRIM (DIAGADJ_DIR) // 'gctm.adj.*'

            WRITE( 6, 102 ) TRIM( REMOVE_ADJ_FILE_CMD )

            CALL SYSTEM ( TRIM( REMOVE_ADJ_FILE_CMD ) )

            ! Remove optimization files now, as would have been done normally
            ! for the "REFERENCE" run at N_CALC_STOP = 0, as the JACOBIAN test
            ! run begins with N_CALC_STOP = 1.
            IF ( N_CALC_STOP == 1 ) THEN

               ! Clean out old *.gdt.*, *.ics.* and cnf.* files
               REMOVE_OPT_FILE_CMD = 'rm ' //
     &                             TRIM (OPTDATA_DIR) // '*.gdt.*'
     &                             // ' ' //
     &                             TRIM (OPTDATA_DIR) // '*.sf.*'
     &                             // ' ' //
     &                             TRIM (OPTDATA_DIR) // 'cfn.*'

               WRITE( 6, 102 ) TRIM( REMOVE_OPT_FILE_CMD )

               CALL SYSTEM ( TRIM( REMOVE_OPT_FILE_CMD ) )

               ! Clean out old *.fd.* files  (dkh, 06/24/09)
               REMOVE_FD_FILE_CMD = 'rm ' //
     &                             TRIM (DIAGADJ_DIR) // '*.fd.*'
     &                             // ' ' //
     &                             TRIM (DIAGADJ_DIR) // '*.fdglob.*'

               WRITE( 6, 102 ) TRIM( REMOVE_FD_FILE_CMD )

               CALL SYSTEM ( TRIM( REMOVE_FD_FILE_CMD ) )

               ! Clean out old ems.adj.* files (dkh, 02/17/11)
               IF ( LEMS_ABS ) THEN
                  REMOVE_ADJ_FILE_CMD = 'rm ' //
     &                                TRIM (DIAGADJ_DIR) // 'ems.adj.*'

                  WRITE( 6, 102 ) TRIM( REMOVE_ADJ_FILE_CMD )

                  CALL SYSTEM ( TRIM( REMOVE_ADJ_FILE_CMD ) )
               ENDIF

               ! Clean out old gctm.iteration file (dkh, 02/17/11)
               IF ( LITR ) THEN
                  REMOVE_ADJ_FILE_CMD = 'rm ' //

     &                          TRIM (DIAGADJ_DIR) // 'gctm.iteration'

                  WRITE( 6, 102 ) TRIM( REMOVE_ADJ_FILE_CMD )

                  CALL SYSTEM ( TRIM( REMOVE_ADJ_FILE_CMD ) )
               ENDIF

            ENDIF

         ENDIF ! LINV

      ENDIF

      END SUBROUTINE CLEAN_FILE_DIRS

!-----------------------------------------------------------------------------------------

      SUBROUTINE INIT_DEP_MAPS
!
!******************************************************************************
!  Subroutine INIT_DEP_MAPS creates mapping arrays for going from tracer and
!  species concentrations to deposition index. (dkh, 05/30/13)
!
!  NOTES:
!
!******************************************************************************
!
      ! References to F90 modules
      USE ADJ_ARRAYS_MOD,       ONLY : NOBS2NDEP
      USE ADJ_ARRAYS_MOD,       ONLY : NOBSCSPEC2NDEP
      USE ADJ_ARRAYS_MOD,       ONLY : NOBS2NWDEP
      USE ADJ_ARRAYS_MOD,       ONLY : NTR2NOBS
      USE ADJ_ARRAYS_MOD,       ONLY : TRACER_IND
      USE ADJ_ARRAYS_MOD,       ONLY : NOBS
      USE ADJ_ARRAYS_MOD,       ONLY : NOBS_CSPEC
      USE ADJ_ARRAYS_MOD,       ONLY : IDCSPEC_ADJ
      USE DRYDEP_MOD,           ONLY : NTRAIND
      USE DRYDEP_MOD,           ONLY : NUMDEP
      USE DRYDEP_MOD,           ONLY : DEPNAME
      USE ERROR_MOD,            ONLY : ERROR_STOP
      USE ERROR_MOD,            ONLY : ALLOC_ERR
      USE TRACER_MOD,           ONLY : N_TRACERS
      USE WETSCAV_MOD,          ONLY : GET_WETDEP_IDWETD
      USE WETSCAV_MOD,          ONLY : NSOL

#     include "CMN_SIZE"
#     include "comode.h"             ! IRM

      ! Local variables
      INTEGER                       :: N
      INTEGER                       :: NN
      INTEGER                       :: AS
      INTEGER                       :: JJ
      INTEGER                       :: NK
      LOGICAL                       :: FOUND

      !=================================================================
      ! INIT_DEP_MAPS begins here!
      !=================================================================

      ALLOCATE( NOBS2NDEP( NOBS ) , STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'NOBS2NDEP' )
      NOBS2NDEP = 0

      ALLOCATE( NOBSCSPEC2NDEP( NOBS_CSPEC ) , STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'NOBSCSPEC2NDEP' )
      NOBSCSPEC2NDEP = 0

      ALLOCATE( NOBS2NWDEP( NOBS ) , STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'NOBS2NWDEP' )
      NOBS2NWDEP = 0

      ALLOCATE( NTR2NOBS( N_TRACERS ) , STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'NTR2NOBS' )
      NTR2NOBS = 0

      ! NOBS2NDEP: Map from NOBS to N_DEP (drydep id)
      DO N = 1, NOBS

         DO NN = 1, NUMDEP
            IF ( NTRAIND(NN) == TRACER_IND(N) ) NOBS2NDEP(N) = NN
         ENDDO

      ENDDO

      ! NOBSCSPEC2NDEP: Map from NOBS_CSPEC to N_DEP (drydep id)
      DO N = 1, NOBS_CSPEC

         !DO NN = 1, NUMDEP
        !    this may not work since CNAME would by DRYNO2 and DEPNAME NO2? or NOx?
        !    IF ( DEPNAME(NN) == CNAME(N) ) NOBSCSPEC2NDEP(N) = NN
        ! ENDDO

         ! Determine drydep ID that corresponds to this species
         NOBSCSPEC2NDEP(N) = -999
         NCS     = NCSURBAN

         DO NN = 1, NUMDEP

            NK = NTDEP(NN)
            IF ( NK <= 0 ) CYCLE
            JJ = IRM(NPRODLO+1,NK,NCS)
            IF ( JJ == IDCSPEC_ADJ(N) ) NOBSCSPEC2NDEP(N) = NN

         ENDDO

         IF ( NOBSCSPEC2NDEP(N) < 0 ) THEN
            CALL ERROR_STOP('Species not in ND44','INIT_CSPEC_ADJ.f')
         ENDIF

      ENDDO

      ! NOBS2NWDEP: Map from NOBS to N_WDEP (wetdep id)
      DO N = 1, NOBS

         ! Get wetdep ID number for this observed tracer
         DO NN = 1, NSOL
            IF ( GET_WETDEP_IDWETD(NN) == TRACER_IND(N) )
     &           NOBS2NWDEP(N) = NN
         ENDDO

      ENDDO

      ! NTR2NOBS: Map from tracer index to observed tracer index
      DO NN = 1, N_TRACERS

         DO N = 1, NOBS

            IF ( TRACER_IND(N) == NN ) NTR2NOBS(NN) = N

         ENDDO

      ENDDO

!      ! NTR2NOBSCSPEC: Map from tracer index to observed species  index
!      DO NN = 1, NTRACER
!
!         DO N = 1, NOBS_CSPEC
!
!            IF ( TRACER_IND(N) == NN ) NTR2NOBSCSPEC(NN) = N
!
!         ENDDO
!
!      ENDDO

      ! Return to calling program
      END SUBROUTINE INIT_DEP_MAPS

!-----------------------------------------------------------------------------

      SUBROUTINE INIT_INPUT_ADJ
!
!******************************************************************************
!  Subroutine INIT_INPUT_ADJ initializes all variables from
!  "directory_adj_mod.f" and "logical_adj_mod.f" for safety's sake.
!  (adj_group, 6/07/09)
!
!  NOTES:
!  (1 ) Add LTES_PSO  (kjw, dkh, 02/12/12, adj32_023)
!  (2 ) Add support for strat fluxes LADJ_STRAT (hml, dkh, 02/14/12, adj32_025)
!  (3 ) Add LINVH_BFGS (nab, 25/03/12 )
!******************************************************************************
!
      ! References to F90 modules
      USE DIRECTORY_ADJ_MOD, ONLY : OPTDATA_DIR, ADJTMP_DIR, DIAGADJ_DIR
      USE LOGICAL_ADJ_MOD,   ONLY : LADJ,LADJ_TRAN, LADJ_CHEM,
     &                              LAERO_THERM, LFD_SPOT, LFD_GLOB,
     &                              LSENS, L4DVAR, L3DVAR, LAPSRC,
     &                              LBKCOV,LINVH,LINVH_BFGS,
     &                              LLINOZ, LFDTEST, LISO,
     &                              LICS, LRXNR, LADJDIAG, LJSAVE,
     &                              LDCOSAT, LHMOD, LHOBS,
     &                              LHMODIFF, LADJ_FORCE, LMODBIAS,
     &                              LOBS_COUNT, LDOFS, LADJ_EMS,
     &                              LDEL_CHKPT, LADJ_TRAJ, LITR,
     &                              LDEVOC, LTES_PSO, LADJ_STRAT,
     &                              LADJ_RRATE,
     &                              LADJ_CL, LADJ_CL_NDEP, LADJ_CL_ACID


      !=================================================================
      ! INIT_INPUT_ADJ begins here!
      !=================================================================

      ! Initialize directories
      OPTDATA_DIR  = ''
      DIAGADJ_DIR  = ''
      ADJTMP_DIR   = ''

      ! Initialize logicals
      LADJ        = .FALSE.
      LADJ_TRAN   = .FALSE.
      LADJ_CHEM   = .FALSE.
      LAERO_THERM = .FALSE.
      LFD_SPOT    = .FALSE.
      LFD_GLOB    = .FALSE.
      LSENS       = .FALSE.
      L4DVAR      = .FALSE.
      L3DVAR      = .FALSE.
      LAPSRC      = .FALSE.
      LBKCOV      = .FALSE.
      LINVH       = .FALSE.
      LINVH_BFGS  = .FALSE.
      LISO        = .FALSE.
      !LLINOZ      = .FALSE.
      LFDTEST     = .FALSE.
      LADJ_EMS    = .FALSE.
      LICS        = .FALSE.
      LRXNR       = .FALSE.
      LADJDIAG    = .FALSE.
      LJSAVE      = .FALSE.
      LADJ_TRAJ   = .FALSE.
      LDCOSAT     = .FALSE.
      LHMOD       = .FALSE.
      LHOBS       = .FALSE.
      LHMODIFF    = .FALSE.
      LADJ_FORCE  = .FALSE.
      LMODBIAS    = .FALSE.
      LOBS_COUNT  = .FALSE.
      LDOFS       = .FALSE.
      LDEL_CHKPT  = .FALSE.
      LITR        = .FALSE.
      LDEVOC      = .TRUE.
      LTES_PSO    = .FALSE.
      LADJ_STRAT  = .FALSE.
      LADJ_RRATE  = .FALSE.
      LADJ_CL     = .FALSE.
      LADJ_CL_NDEP= .FALSE.
      LADJ_CL_ACID= .FALSE.

      ! Initialize counters
      CT1          = 0
      CT2          = 0
      CT3          = 0

      ! Return to calling program
      END SUBROUTINE INIT_INPUT_ADJ

!------------------------------------------------------------------------------

      SUBROUTINE READ_STRID_FILE
!
!******************************************************************************
!  Subroutine READ_STRID_FILE reads the list of stratospheric production
!  and loss rates from STR_ID file in run directory.
!  (hml, 05/22/13)
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE ADJ_ARRAYS_MOD,    ONLY : NSTPL
      USE ADJ_ARRAYS_MOD,    ONLY : ID_PROD
      USE ADJ_ARRAYS_MOD,    ONLY : ID_LOSS
      USE ADJ_ARRAYS_MOD,    ONLY : OPT_THIS_PROD
      USE ADJ_ARRAYS_MOD,    ONLY : OPT_THIS_LOSS
      USE ADJ_ARRAYS_MOD,    ONLY : PROD_SF_DEFAULT
      USE ADJ_ARRAYS_MOD,    ONLY : LOSS_SF_DEFAULT
      USE ADJ_ARRAYS_MOD,    ONLY : PROD_NAME
      USE ADJ_ARRAYS_MOD,    ONLY : LOSS_NAME
      USE ADJ_ARRAYS_MOD,    ONLY : REG_PARAM_PROD
      USE ADJ_ARRAYS_MOD,    ONLY : REG_PARAM_LOSS
      USE ADJ_ARRAYS_MOD,    ONLY : PROD_ERROR
      USE ADJ_ARRAYS_MOD,    ONLY : LOSS_ERROR
      USE FILE_MOD,          ONLY : IOERROR, IU_STR


      ! local variables
      INTEGER                    :: IOS
      INTEGER                    :: T
      INTEGER                    :: N
      CHARACTER(LEN=255)         :: SUBSTRS(MAXDIM)


      !=================================================================
      ! READ_STRID_FILE begins here!
      !=================================================================

      ! Open STR_ID file containing list of 24x2 strat prod & loss
      OPEN( IU_STR, FILE='STR_ID', STATUS='OLD', IOSTAT=IOS )

      IF ( IOS /= 0 ) CALL IOERROR(IOS,IU_STR,'read_input_adj_strat:1')

      READ_STR_ID = .TRUE.

      ! Read a header line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_STR_ID:1')

      DO T = 1, NSTPL
         ! Split line into substrings
         CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_strat_p')

         ! Save tracer number
         READ( SUBSTRS(1), * ) ID_PROD(T)

         ! Save tracer name
         PROD_NAME(T)  = TRIM( SUBSTRS(2) )

         ! optimize this strat prod & loss?
         READ( SUBSTRS(3), *) OPT_THIS_PROD(T)

         ! Defualt prod scaling factor for this strat tracer
         READ( SUBSTRS(4), *) PROD_SF_DEFAULT(T)

         ! REG_PARAM for this strat tracer
         READ( SUBSTRS(5), *) REG_PARAM_PROD(T)

         ! STR_ERROR for this strat tracer
         READ( SUBSTRS(6), *) PROD_ERROR(T)

      ENDDO

      DO T = 1, NSTPL

         ! Split line into substrings
         CALL SPLIT_ONE_LINE( SUBSTRS, N, -1,'read_strat_l')

         ! Save tracer number
         READ( SUBSTRS(1), * ) ID_LOSS(T)

         ! Save tracer name
         LOSS_NAME(T)  = TRIM( SUBSTRS(2) )

         ! optimize this strat prod & loss?
         READ( SUBSTRS(3), *) OPT_THIS_LOSS(T)

         ! Defualt loss scaling factor for this strat tracer
         READ( SUBSTRS(4), *) LOSS_SF_DEFAULT(T)

         ! REG_PARAM for this strat tracer
         READ( SUBSTRS(5), *) REG_PARAM_LOSS(T)

         ! STR_ERROR for this strat tracer
         READ( SUBSTRS(6), *) LOSS_ERROR(T)

      ENDDO

      CLOSE(IU_STR)

      READ_STR_ID = .FALSE.

      ! Return to calling program
      END SUBROUTINE READ_STRID_FILE

!------------------------------------------------------------------------------

      SUBROUTINE READ_RXNID_FILE
!
!******************************************************************************
!  Subroutine READ_RXNID_FILE reads the list of kpp reactions
!  from RXN_ID file in run directory.
!  (hml, 05/22/13)
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE FILE_MOD,          ONLY : IOERROR, IU_RXN
      USE ADJ_ARRAYS_MOD,    ONLY : NRRATES
      USE ADJ_ARRAYS_MOD,    ONLY : ID_RRATES
      USE ADJ_ARRAYS_MOD,    ONLY : OPT_THIS_RATE
      USE ADJ_ARRAYS_MOD,    ONLY : RATE_SF_DEFAULT
      USE ADJ_ARRAYS_MOD,    ONLY : RRATES_NAME
      USE ADJ_ARRAYS_MOD,    ONLY : REG_PARAM_RATE
      USE ADJ_ARRAYS_MOD,    ONLY : RATE_ERROR

      ! local variables
      INTEGER                    :: IOS
      INTEGER                    :: N
      INTEGER                    :: T
      CHARACTER(LEN=255)         :: SUBSTRS(MAXDIM)


      !=================================================================
      ! READ_RXNID_FILE begins here!
      !=================================================================

      ! Open RXN_ID file containing list of 297 reactions
      OPEN( IU_RXN, FILE='RXN_ID', STATUS='OLD', IOSTAT=IOS )

      IF ( IOS /= 0 ) CALL IOERROR(
     &                     IOS,IU_RXN,'read_input_adj_rrate:1')

      READ_RXN_ID = .TRUE.

      ! Read a header line
      CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_RXN_ID:1')

      DO T = 1, NRRATES
         ! Split line into substrings
         CALL SPLIT_ONE_LINE( SUBSTRS, N, -1, 'read_RXN_ID:2')

         ! Save tracer number
         READ( SUBSTRS(1), * ) ID_RRATES(T)

         ! Save tracer name
         RRATES_NAME(T) = TRIM( SUBSTRS(2) )

         ! Optimize this rate?
         READ( SUBSTRS(3), * ) OPT_THIS_RATE(T)

         ! Default scaling factor for this rate
         READ( SUBSTRS(4), * ) RATE_SF_DEFAULT(T)

         ! REG_PARAM for this rate
         READ( SUBSTRS(5), * ) REG_PARAM_RATE(T)

         ! RATE_ERROR for this rate
         READ( SUBSTRS(6), * ) RATE_ERROR(T)

      ENDDO

      CLOSE(IU_RXN)

      READ_RXN_ID = .FALSE.


      ! Return to calling program
      END SUBROUTINE READ_RXNID_FILE

!------------------------------------------------------------------------------

      ! End of module
      END MODULE INPUT_ADJ_MOD
