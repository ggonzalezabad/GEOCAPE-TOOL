MODULE GC_clouds_module

  USE GC_parameters_module, ONLY: maxcld, GC_maxlayers, cldunit,         &
                                  max_ch_len, maxgases, maxaer
  USE GC_variables_module,  ONLY: cldfile, do_clouds, do_lambertian_cld, &
                                  use_cldprof, ngases, nmessages,        &
                                  messages, GC_nlayers, pressures,       &
                                  heights, pmaxl, pmaxc, pnc, profids,   &
                                  profile_data, temperatures, aircolumns,&
                                  daircolumns_dT, gas_partialcolumns,    &
                                  naer0, aer_profile0, cld_types,        &
                                  cld_bots, cld_profile, cld_tops,       &
                                  cld_flags, cld_taus, tcld_profile,     &
                                  tcldtau0, cld_uppers, cld_lowers, ncld,&
                                  cfrac, cld_reflambda, use_zcldspec,    &
                                  lambertian_cldalb, do_ctp_Jacobians,   &
                                  do_cfrac_Jacobians, do_cssa_Jacobians, &
                                  do_cod_Jacobians, do_cld_columnwf
  USE GC_read_input_module, ONLY: skip_to_filemark
  USE GC_error_module

  IMPLICIT NONE

  PUBLIC :: read_cld_control_file
  PUBLIC :: cloud_profiles

  CONTAINS

    SUBROUTINE read_cld_control_file(error)

      IMPLICIT NONE

      LOGICAL, INTENT(OUT) :: error

      ! ----------------------------------
      ! Strings in the clouds control file
      ! ----------------------------------
      CHARACTER(LEN=31), PARAMETER :: cld_columnwf_str    = 'Cloud column weighting function'
      CHARACTER(LEN=13), PARAMETER :: cod_jacobians_str   = 'COD jacobians'
      CHARACTER(LEN=14), PARAMETER :: cssa_jacobians_str  = 'CSSA jacobians'
      CHARACTER(LEN=15), PARAMETER :: cfrac_jacobians_str = 'CFRAC jacobians'
      CHARACTER(LEN=13), PARAMETER :: ctp_jacobians_str   = 'CTP jacobians'
      CHARACTER(LEN=23), PARAMETER :: lcld_albedo_str     = 'Lambertian cloud albedo'
      CHARACTER(LEN=12), PARAMETER :: zcldspec_str        = 'Use zcldspec'
      CHARACTER(LEN=14), PARAMETER :: cod_wavelength_str  = 'COD wavelength'
      CHARACTER(LEN=14), PARAMETER :: cld_fraction_str    = 'Cloud fraction'
      CHARACTER(LEN= 9), PARAMETER :: cld_top_str         = 'Cloud top'
      CHARACTER(LEN=16), PARAMETER :: cld_types_str       = 'Number of clouds'
      
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
      funit = cldunit

      ! =======================
      ! Read cloud control file
      ! =======================

      OPEN(UNIT=funit, FILE = cldfile, ACTION='READ', status = 'old', IOSTAT = ios)
      IF ( ios /= 0 ) THEN
         CALL write_err_message ( .TRUE., 'ERROR: Unable to open '//cldfile )
      END IF
      
      ! --------------------------------
      ! Clouds column weighting function
      ! --------------------------------
      REWIND (funit)
      CALL skip_to_filemark (funit, cld_columnwf_str , tmpstr, error )
      IF ( error ) THEN
         CALL write_err_message ( .TRUE., "ERROR: Can't find "//cld_columnwf_str )
      END IF
      READ (UNIT=funit, FMT=*, IOSTAT=ios) do_cld_columnwf

      ! -------------
      ! COD Jacobians
      ! -------------
      REWIND (funit)
      CALL skip_to_filemark (funit, cod_jacobians_str , tmpstr, error )
      IF ( error ) THEN
         CALL write_err_message ( .TRUE., "ERROR: Can't find "//cod_jacobians_str )
      END IF
      READ (UNIT=funit, FMT=*, IOSTAT=ios) do_cod_Jacobians

      ! --------------
      ! CSSA Jacobians
      ! --------------
      REWIND (funit)
      CALL skip_to_filemark (funit, cssa_jacobians_str , tmpstr, error )
      IF ( error ) THEN
         CALL write_err_message ( .TRUE., "ERROR: Can't find "//cssa_jacobians_str )
      END IF
      READ (UNIT=funit, FMT=*, IOSTAT=ios) do_cssa_Jacobians

      ! ---------------
      ! CFRAC Jacobians
      ! ---------------
      REWIND (funit)
      CALL skip_to_filemark (funit, cfrac_jacobians_str , tmpstr, error )
      IF ( error ) THEN
         CALL write_err_message ( .TRUE., "ERROR: Can't find "//cfrac_jacobians_str )
      END IF
      READ (UNIT=funit, FMT=*, IOSTAT=ios) do_cfrac_Jacobians


      ! -------------
      ! CTP Jacobians
      ! -------------
      REWIND (funit)
      CALL skip_to_filemark (funit, ctp_jacobians_str , tmpstr, error )
      IF ( error ) THEN
         CALL write_err_message ( .TRUE., "ERROR: Can't find "//ctp_jacobians_str )
      END IF
      READ (UNIT=funit, FMT=*, IOSTAT=ios) do_ctp_Jacobians

      ! -----------------------
      ! Lambertian cloud albedo
      ! -----------------------
      REWIND (funit)
      CALL skip_to_filemark (funit, lcld_albedo_str , tmpstr, error )
      IF ( error ) THEN
         CALL write_err_message ( .TRUE., "ERROR: Can't find "//lcld_albedo_str )
      END IF
      READ (UNIT=funit, FMT=*, IOSTAT=ios) lambertian_cldalb

      ! ------------
      ! Use ZCLDSPEC
      ! ------------
      REWIND (funit)
      CALL skip_to_filemark (funit, zcldspec_str , tmpstr, error )
      IF ( error ) THEN
         CALL write_err_message ( .TRUE., "ERROR: Can't find "//zcldspec_str )
      END IF
      READ (UNIT=funit, FMT=*, IOSTAT=ios) use_zcldspec

      ! --------------------------
      ! CLOUD reference wavelength
      ! --------------------------
      REWIND (funit)
      CALL skip_to_filemark (funit, cod_wavelength_str , tmpstr, error )
      IF ( error ) THEN
         CALL write_err_message ( .TRUE., "ERROR: Can't find "//cod_wavelength_str )
      END IF
      READ (UNIT=funit, FMT=*, IOSTAT=ios) cld_reflambda

      IF (.NOT. use_cldprof) THEN

         ! --------------
         ! Cloud fraction
         ! --------------
         REWIND (funit)
         CALL skip_to_filemark (funit, cld_fraction_str , tmpstr, error )
         IF ( error ) THEN
            CALL write_err_message ( .TRUE., "ERROR: Can't find "//cld_fraction_str )
         END IF
         READ (UNIT=funit, FMT=*, IOSTAT=ios) cfrac

         ! ---------
         ! Cloud top
         ! ---------
         REWIND (funit)
         CALL skip_to_filemark (funit, cld_top_str , tmpstr, error )
         IF ( error ) THEN
            CALL write_err_message ( .TRUE., "ERROR: Can't find "//cld_top_str )
         END IF
         READ (UNIT=funit, FMT=*, IOSTAT=ios) cld_tops(1)

         IF (.NOT. do_lambertian_cld .AND. cfrac .GT. 0.0d0) THEN
            
            ! -------------
            ! COD Jacobians
            ! -------------
            REWIND (funit)
            CALL skip_to_filemark (funit, cld_types_str , tmpstr, error )
            IF ( error ) THEN
               CALL write_err_message ( .TRUE., "ERROR: Can't find "//cld_types_str )
            END IF
            READ (UNIT=funit, FMT=*, IOSTAT=ios) ncld
            DO i = 1, ncld
               READ (UNIT=funit, FMT=*, IOSTAT=ios) cld_types(i)
               READ (UNIT=funit, FMT=*, IOSTAT=ios) cld_bots(i), cld_tops(i), &
                                                     cld_taus(i)
            END DO
         END IF
      END IF

      ! -----------------------------------------
      ! Possible erros when reading from the file
      ! -----------------------------------------
      IF (ios .NE. 0) error = .TRUE.

    END SUBROUTINE read_cld_control_file

    SUBROUTINE cloud_profiles(error)

      IMPLICIT NONE

      LOGICAL, INTENT(INOUT) :: error ! Error variable

      INTEGER :: i

      ! ----------------
      ! Code starts here
      ! ----------------
      error = .FALSE.

      IF (cfrac .le. 0.0d0) THEN
         ncld               = 0
         cfrac              = 0.0d0
         do_cld_columnwf    = .FALSE.
         do_clouds          = .FALSE.
         do_cod_Jacobians   = .FALSE.
         do_cssa_Jacobians  = .FALSE.
         do_ctp_Jacobians   = .FALSE. 
         do_cfrac_Jacobians = .FALSE.  
      ELSE
         IF (cfrac .ge. 1.0d0) THEN
            cfrac = 1.0d0
            do_cfrac_Jacobians = .FALSE.
         END IF
      
         IF (do_lambertian_cld) THEN
            do_cod_Jacobians  = .FALSE.
            do_cssa_jacobians = .FALSE.
            ncld = 1
         ENDIF
      
         ! ---------------------------------------------------------------------
         ! convert from cloud specification in pressure to altitude if necessary
         ! ---------------------------------------------------------------------
         IF (.NOT. use_zcldspec) THEN
            CALL convert_cldspec_p2z(GC_nlayers, pressures(0:GC_nlayers),      &
                                     heights(0:GC_nlayers), do_lambertian_cld, &
                                     use_cldprof, maxcld, ncld, cld_tops,      &
                                     cld_bots)
         END IF
      
         IF (do_lambertian_cld) THEN
            IF (cld_tops(1) < heights(GC_nlayers)) CALL write_err_message &
                 ( .TRUE., "Lambertian cloud surface below surface!!!")
            
            IF (cld_tops(1) > heights(0))          CALL write_err_message &
                 ( .TRUE., "Lambertian cloud surface above TOA!!!")
         END IF
 
         IF (use_cldprof .and. .not. do_lambertian_cld) THEN
            CALL get_cldprof_from_atmosprof(pmaxl, pmaxc, GC_nlayers, pnc, &
                 profids, profile_data, maxcld, ncld, cld_types,           &
                 cld_profile(1:maxcld, 1:GC_nlayers), cld_lowers, cld_uppers)

         IF (ncld == 0) THEN
            cfrac              = 0.0d0
            do_clouds          = .FALSE.
            do_cod_Jacobians   = .FALSE.
            do_cssa_Jacobians  = .FALSE.
            do_ctp_Jacobians   = .FALSE.  
            do_cfrac_Jacobians = .FALSE.
         END IF
      ELSE
         ! ------------------
         ! Check cloud set up
         ! ------------------
         IF (ncld > maxcld) THEN
            CALL write_err_message ( .TRUE., &
                 "Too many clouds, increase maxcld or reduce ncld!!!")
            error = .TRUE.
         END IF
         IF (.NOT. do_lambertian_cld) THEN
            DO i = 1, ncld
               IF (cld_taus(i) <= 0.d0) THEN 
                  CALL write_err_message ( .TRUE., &
                       "Cloud optical depth must be greater than 0!!!")
                  error = .TRUE.
               END IF
               IF (cld_bots(i) >= cld_tops(i) ) THEN
                  CALL write_err_message ( .TRUE., &
                       "Cloud bottom must be lower than cloud top!!!'")
                  error = .TRUE.
               END IF
               IF (cld_bots(i) < heights(GC_nlayers)) THEN 
                  CALL write_err_message ( .TRUE., &
                       "Cloud bottom below surface!!!")
                  error = .TRUE.
               END IF
               IF (cld_tops(i) > heights(0)) THEN 
                  CALL write_err_message ( .TRUE., &
                       "Cloud top should not be above TOA!!!")
                  error = .TRUE.
               END IF
            END DO
         END IF

         IF (error) CALL error_exit (error)
         
         ! --------------------
         ! Insert clouds levels
         ! --------------------
         call insert_clouds(GC_maxlayers, GC_nlayers, heights, pressures, temperatures, &
              aircolumns, daircolumns_dT, gas_partialcolumns, maxgases, ngases,         &
              do_lambertian_cld, maxcld, ncld, cld_bots, cld_tops, cld_taus,            &
              cld_lowers, cld_uppers, cld_profile, maxaer, naer0, aer_profile0,         &
              error, messages(nmessages+1))

         IF (error) CALL write_err_message ( .TRUE., messages(nmessages+1))

      END IF
      
      DO i = 1, ncld
         WHERE (cld_profile(i, 1:GC_nlayers) > 0.0d0)
            cld_flags(1:GC_nlayers) = .TRUE.
         END WHERE
      END DO

      DO i = 1, GC_nlayers
         IF (cld_flags(i)) tcld_profile(i) = SUM(cld_profile(1:ncld, i))
      END DO
   END IF

   tcldtau0 = SUM(tcld_profile(1:GC_nlayers)) 

 END SUBROUTINE cloud_profiles
       
END MODULE GC_clouds_module
