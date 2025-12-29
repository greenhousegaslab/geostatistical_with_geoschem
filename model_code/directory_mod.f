! $Id: directory_mod.f,v 1.1 2009/06/09 21:51:53 daven Exp $
      MODULE DIRECTORY_MOD
!
!******************************************************************************
!  Module DIRECTORY_MOD contains the directory path variables used by
!  GEOS-CHEM. (bmy, 7/20/04, 10/24/05)
!
!  Module Variables:
!  ============================================================================
!  (1 ) DATA_DIR     (CHAR*255) : Main DAO met field directory
!  (2 ) DATA_DIR_1x1 (CHAR*255) : Root data dir for 1x1 emission fields
!  (2 ) GCAP_DIR     (CHAR*255) : Subdir where GCAP       met data are stored
!  (3 ) GEOS_1_DIR   (CHAR*255) : Subdir where GEOS-1     met data are stored
!  (4 ) GEOS_S_DIR   (CHAR*255) : Subdir where GEOS-STRAT met data are stored
!  (5 ) GEOS_3_DIR   (CHAR*255) : Subdir where GEOS-3     met data are stored
!  (6 ) GEOS_4_DIR   (CHAR*255) : Subdir where GEOS-4     met data are stored
!  (7 ) GEOS_5_DIR   (CHAR*255) : Subdir where GEOS-5     met data are stored
!  (8 ) TEMP_DIR     (CHAR*255) : Temporary directory for unzipping met dat
!  (9 ) RUN_DIR      (CHAR*255) : Run directory for GEOS-CHEM
!  (10) OH_DIR       (CHAR*255) : Dir containing OH files are stored
!  (11) O3PL_DIR     (CHAR*255) : Dir containing archived O3 P/L rate files
!  (12) TPBC_DIR     (CHAR*255) : TPCORE boundary conditions dir (nested grid)
!
!  NOTES:
!  (1 ) Added variables GCAP_DIR and GEOS_5_DIR (swu, bmy, 5/25/05)
!  (2 ) Added DATA_DIR_1x1 (bmy, 10/24/05)
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================
      CHARACTER(LEN=255) :: DATA_DIR
      CHARACTER(LEN=255) :: DATA_DIR_1x1
      CHARACTER(LEN=255) :: GCAP_DIR
      CHARACTER(LEN=255) :: GEOS_1_DIR
      CHARACTER(LEN=255) :: GEOS_S_DIR
      CHARACTER(LEN=255) :: GEOS_3_DIR
      CHARACTER(LEN=255) :: GEOS_4_DIR
      CHARACTER(LEN=255) :: GEOS_5_DIR
      CHARACTER(LEN=255) :: GEOS_FP_DIR     !! (lzh,04/09/2014)
      CHARACTER(LEN=255) :: MERRA2_DIR      !! (mps,10/23/2015)
      CHARACTER(LEN=255) :: TEMP_DIR
      CHARACTER(LEN=255) :: RUN_DIR
      CHARACTER(LEN=255) :: OH_DIR
      CHARACTER(LEN=255) :: O3PL_DIR
      CHARACTER(LEN=255) :: TPBC_DIR
      CHARACTER(LEN=255) :: TPBC_DIR_NA     ! (lzh, 02//01/2015, update nested runs)
      CHARACTER(LEN=255) :: TPBC_DIR_EU
      CHARACTER(LEN=255) :: TPBC_DIR_CH
      CHARACTER(LEN=255) :: FLUX_DIR
      CHARACTER(LEN=255) :: OBS_DIR_OCO2
      CHARACTER(LEN=255) :: OBS_DIR_FLASK
      CHARACTER(LEN=255) :: FFILE_PREFIX
      CHARACTER(LEN=255) :: OCOFILE_PREFIX

      ! End of module
      END MODULE DIRECTORY_MOD
