!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: GC_aerosols_module.f90
!
! !DESCRIPTION: This module declares variables needed for the aerosols calculations
!               and contains routines related to its usage.
!\\
!\\
! !INTERFACE: 
!
MODULE GC_aerosols_module
! 
!  !USES:
!
  USE GC_parameters_module, ONLY: GC_maxlayers, aerunit, max_ch_len
  USE GC_variables_module,  ONLY: naer0, taertau0, &
                                  aer_profile, aer_d_profile_dtau,        &
                                  aer_d_profile_dpkh, aer_d_profile_dhfw, &
                                  aer_profile0, taer_profile, &
                                  aer_flags, aer_opdeps, aer_ssalbs,      &
                                  aer_relqext, aer_phfcn, GC_nlayers,     &
                                  heights, nmessages, messages, aer_ctr
  USE GC_read_input_module, ONLY: skip_to_filemark
  USE GC_error_module
!
  IMPLICIT NONE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: read_aer_control_file
  PUBLIC :: aerosol_profiles
!
! !PUBLIC DATA MEMBERS:
!
  ! -------------------
  ! File names for data
  ! -------------------
! 
! !REVISION HISTORY:
!  April 2013 - G. Gonzalez Abad - Initial Version
!
!EOP
!------------------------------------------------------------------------------
  CONTAINS

    SUBROUTINE read_aer_control_file(error)

      IMPLICIT NONE

      LOGICAL, INTENT(INOUT) :: error ! Error variable

      ! ------------------------------------
      ! Strings in the aerosols control file
      ! ------------------------------------
      CHARACTER(LEN=33), PARAMETER :: aer_columnwf_str    = 'Aerosol column weighting function'
      CHARACTER(LEN=13), PARAMETER :: aod_jacobians_str   = 'AOD jacobians'
      CHARACTER(LEN=14), PARAMETER :: assa_jacobians_str  = 'ASSA jacobians'
      CHARACTER(LEN=15), PARAMETER :: aerph_jacobians_str = 'AERPH jacobians'
      CHARACTER(LEN=15), PARAMETER :: aerhw_jacobians_str = 'AERHW jacobians'
      CHARACTER(LEN=14), PARAMETER :: aod_wavelength_str  = 'AOD wavelength'
      CHARACTER(LEN=15), PARAMETER :: aerosols_type_str   = 'Aerosol types'

      ! ---------
      ! File unit
      ! ---------
      INTEGER :: funit

      ! ---------------
      ! Temporal string
      ! ---------------
      CHARACTER(LEN=max_ch_len) :: tmpstr

      ! ---------------
      ! Index variables
      ! ---------------
      INTEGER :: i

      ! ------
      ! Errors
      ! ------
      INTEGER :: ios

      ! ----------------
      ! Code starts here
      ! ----------------
      error = .FALSE.

      ! -----------------
      ! Control file unit
      ! -----------------
      funit = aerunit

      ! ==========================
      ! Read aerosols control file
      ! ==========================

      OPEN(UNIT=funit, FILE = aer_ctr%aerfile, ACTION='READ', status = 'old', IOSTAT = ios)
      IF ( ios /= 0 ) THEN
         CALL write_err_message ( .TRUE., 'ERROR: Unable to open '//aer_ctr%aerfile )
      END IF
      
      ! ----------------------------------
      ! Aerosols column weighting function
      ! ----------------------------------
      REWIND (funit)
      CALL skip_to_filemark (funit, aer_columnwf_str , tmpstr, error )
      IF ( error ) THEN
         CALL write_err_message ( .TRUE., "ERROR: Can't find "//aer_columnwf_str )
      END IF
      READ (UNIT=funit, FMT=*, IOSTAT=ios) aer_ctr%do_aer_columnwf
      
      ! -------------
      ! AOD Jacobians
      ! -------------
      REWIND (funit)
      CALL skip_to_filemark (funit, aod_jacobians_str , tmpstr, error )
      IF ( error ) THEN
         CALL write_err_message ( .TRUE., "ERROR: Can't find "//aod_jacobians_str )
      END IF
      READ (UNIT=funit, FMT=*, IOSTAT=ios) aer_ctr%do_aod_jacobians
      
      ! --------------
      ! ASSA Jacobians
      ! --------------
      REWIND (funit)
      CALL skip_to_filemark (funit, assa_jacobians_str , tmpstr, error )
      IF ( error ) THEN
         CALL write_err_message ( .TRUE., "ERROR: Can't find "//assa_jacobians_str )
      END IF
      READ (UNIT=funit, FMT=*, IOSTAT=ios) aer_ctr%do_assa_jacobians
      
      ! ---------------
      ! AERPH Jacobians
      ! ---------------
      REWIND (funit)
      CALL skip_to_filemark (funit, aerph_jacobians_str , tmpstr, error )
      IF ( error ) THEN
         CALL write_err_message ( .TRUE., "ERROR: Can't find "//aerph_jacobians_str )
      END IF
      READ (UNIT=funit, FMT=*, IOSTAT=ios) aer_ctr%do_aerph_jacobians
      
      ! ---------------
      ! AERHW Jacobians
      ! ---------------
      REWIND (funit)
      CALL skip_to_filemark (funit, aerhw_jacobians_str , tmpstr, error )
      IF ( error ) THEN
         CALL write_err_message ( .TRUE., "ERROR: Can't find "//aerhw_jacobians_str )
      END IF
      READ (UNIT=funit, FMT=*, IOSTAT=ios) aer_ctr%do_aerhw_jacobians
      
      ! --------------
      ! AOD wavelength
      ! --------------
      REWIND (funit)
      CALL skip_to_filemark (funit,  aod_wavelength_str , tmpstr, error )
      IF ( error ) THEN
         CALL write_err_message ( .TRUE., "ERROR: Can't find "//  aod_wavelength_str)
      END IF
      READ (UNIT=funit, FMT=*, IOSTAT=ios) aer_ctr%aer_reflambda
     
      IF (.NOT. aer_ctr%use_aerprof) THEN
         ! ---------------------------------------------------------------
         ! Six aerosol types: dust, sulfate, organic carbon, black carbon, 
         ! sea salt fine, sea salt coarse
         ! ---------------------------------------------------------------
         ! -------------
         ! Aerosol types
         ! -------------
         REWIND (funit)
         CALL skip_to_filemark (funit, aerosols_type_str, tmpstr, error )
         IF ( error ) THEN
            CALL write_err_message ( .TRUE., "ERROR: Can't find "//aerosols_type_str)
         END IF
         READ (UNIT=funit, FMT=*, IOSTAT=ios) aer_ctr%naer
         
         ALLOCATE(aer_ctr%aer_types(1:aer_ctr%naer))
         ALLOCATE(aer_ctr%aer_tau0s(1:aer_ctr%naer))
         ALLOCATE(aer_ctr%aer_z_upperlimit(1:aer_ctr%naer))
         ALLOCATE(aer_ctr%aer_z_lowerlimit(1:aer_ctr%naer))
         ALLOCATE(aer_ctr%aer_z_peakheight(1:aer_ctr%naer))
         ALLOCATE(aer_ctr%aer_half_width(1:aer_ctr%naer))
         ALLOCATE(aer_ctr%aer_relax(1:aer_ctr%naer))

         DO i = 1, aer_ctr%naer
            READ (UNIT=funit, FMT=*, IOSTAT=ios) aer_ctr%aer_types(i)
            READ (UNIT=funit, FMT=*, IOSTAT=ios) aer_ctr%aer_tau0s(i)
            READ (UNIT=funit, FMT=*, IOSTAT=ios) aer_ctr%aer_z_upperlimit(i),   &
                 aer_ctr%aer_z_lowerlimit(i),   &
                 aer_ctr%aer_z_peakheight(i)
            READ (UNIT=funit, FMT=*, IOSTAT=ios) aer_ctr%aer_half_width(i)
         ENDDO
      ENDIF
      
      ! -----------------------------------------
      ! Possible erros when reading from the file
      ! -----------------------------------------
      IF (ios .NE. 0) error = .TRUE.
      
      CLOSE(funit)
      
    END SUBROUTINE read_aer_control_file

    SUBROUTINE aerosol_profiles(error)
   
      IMPLICIT NONE
   
      LOGICAL, INTENT(INOUT) :: error ! Error variable

      ! ---------------
      ! Index variables
      ! ---------------
      INTEGER :: i

      ! ----------------
      ! Code starts here
      ! ----------------
      error = .FALSE.

      IF (aer_ctr%use_aerprof) THEN
         aer_profile = aer_profile0; aer_ctr%naer = naer0
         IF (aer_ctr%naer == 0) THEN
            aer_ctr%do_aerosols        = .FALSE.
            aer_ctr%do_aer_columnwf    = .FALSE.
            aer_ctr%do_aod_Jacobians   = .FALSE.
            aer_ctr%do_assa_Jacobians  = .FALSE.
            aer_ctr%do_aerph_Jacobians = .FALSE.
            aer_ctr%do_aerhw_Jacobians = .FALSE.
         END IF

      ELSE
         DO i = 1, aer_ctr%naer

            IF (aer_ctr%aer_tau0s(i) <= 0.d0) THEN
               CALL write_err_message ( .TRUE., &
                    "ERROR: Aerosol optical depth must be greater than 0")
            ENDIF
            IF (aer_ctr%aer_z_lowerlimit(i) >= aer_ctr%aer_z_upperlimit(i) ) THEN
               CALL write_err_message ( .TRUE., &
                    "ERROR: Aerosol bottom must be lower than aerosol top")
            ENDIF
            IF (aer_ctr%aer_z_peakheight(i) > aer_ctr%aer_z_upperlimit(i)) THEN
               CALL write_err_message ( .TRUE., &
                    "ERROR: Aerosol peak height should be <= upperlimit")
            END IF
            IF (aer_ctr%aer_z_peakheight(i) < aer_ctr%aer_z_lowerlimit(i)) THEN
               CALL write_err_message ( .TRUE., &
                    "ERROR: Aerosol peak height should be >= lowerlimit")
            ENDIF
            
            IF (aer_ctr%aer_z_lowerlimit(i) < heights(GC_nlayers)) &
                 aer_ctr%aer_z_lowerlimit(i) = heights(GC_nlayers)
            IF (aer_ctr%aer_z_upperlimit(i) > heights(0)) &
                 aer_ctr%aer_z_upperlimit(i) = heights(0)
            
            ! ----------------------------------------------
            ! Generate aerosol plume/profiles based on input
            ! ----------------------------------------------
            CALL generate_plume                                            &
                 ( GC_maxlayers, aer_ctr%aer_z_upperlimit(i), &
                 aer_ctr%aer_z_lowerlimit(i), & ! input
                 aer_ctr%aer_z_peakheight(i), &
                 aer_ctr%aer_tau0s(i), aer_ctr%aer_half_width(i),     & ! input
                 GC_nlayers, heights(0:GC_maxlayers),                      & ! input
                 aer_profile(i, 1:GC_maxlayers),                           & ! output
                 aer_d_profile_dtau(i, 1:GC_maxlayers),                    & ! output   
                 aer_d_profile_dpkh(i,1:GC_maxlayers),                     & ! output
                 aer_d_profile_dhfw(i,1:GC_maxlayers),                     & ! output
                 error, messages(nmessages+1) ) 

            ! ------
            ! Errors
            ! ------
            IF (error) CALL write_err_message (.TRUE., messages(nmessages+1))

         END DO
      END IF
      
      DO i = 1, aer_ctr%naer
         WHERE (aer_profile(i, 1:GC_nlayers) > 0.0d0)
            aer_flags(1:GC_nlayers) = .TRUE.
         END WHERE
      END DO
      
      DO i = 1, GC_nlayers
         IF (aer_flags(i)) taer_profile(i) = SUM(aer_profile(1:aer_ctr%naer, i))
      END DO

    END SUBROUTINE aerosol_profiles

END MODULE GC_aerosols_module
