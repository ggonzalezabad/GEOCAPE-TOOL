!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: GC_read_input_module.f90
!
! !DESCRIPTION: This module declares variables that control the GEOCAPE tool and
!               contains subroutines to read them or derive them from the 
!               previously readed variables. The control file name has to be
!               supplied when calling the main program.
!\\
!\\
! !INTERFACE: 
!
MODULE GC_read_input_module
!
! !USES
  USE GC_parameters_module, ONLY: max_ch_len, maxlambdas, maxflambdas, maxmoms, ctrunit
  USE GC_variables_module,  ONLY: database_dir, results_dir, profile_data_filename,   &
                                  debug_filename, albspectra_fname, aerfile, cldfile, &
                                  do_aerosols, do_clouds, do_lambertian_cld,          &
                                  use_aerprof, use_cldprof,                           &
                                  do_StokesQU_output, idix,                           &
                                  do_AMF_calculation, do_T_Jacobians,                 &
                                  do_sfcprs_Jacobians, do_normalized_WFoutput,        &
                                  do_normalized_radiance,                             &
                                  aercld_nmoments_input,                              &
                                  GC_n_user_altitudes,                                &
                                  GC_user_altitudes, GC_do_user_altitudes,            &
                                  lambda_start, lambda_finish, use_wavelength,        &
                                  do_effcrs, lambda_resolution, lambda_dw,            &
                                  lambda_dfw, nlambdas, lambdas, wavnums, nflambdas,  &
                                  nclambdas, edge_dw, flambdas, fwavnums, clambdas,   &
                                  cwavnums, use_footprint_info, do_debug_geocape_tool,&
                                  use_solar_photons, ngases, which_gases,             &
                                  use_lambertian, do_brdf_surface, use_fixedalbedo,   &
                                  use_albspectra, fixed_albedo, is_ocean, wavlens_um, &
                                  wind_speed, water_rn, water_cn, do_sat_viewcalc,    &
                                  year, month, day, latitude, longitude, satlon,      &
                                  satlat, satalt, utc, geometry_data,                 &
                                  solar_spec_filename
  USE GC_error_module
! 
  IMPLICIT NONE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: read_control_file
  PUBLIC :: check_input
  PUBLIC :: skip_to_filemark
!
! !PUBLIC DATA MEMBERS:
!

!
! !REVISION HISTORY:
!  April 2013 - G. Gonzalez Abad - Initial Version
!
!EOP
!------------------------------------------------------------------------------

  CONTAINS

!BOP
!
! !IROUTINE:  read_control_file
!
! !DESCRIPTION: This routine does something to the input variable and returns
! the result in the output variable.
!\\
!\\
! !INTERFACE: read_control_file
!
    SUBROUTINE read_control_file(error)

      IMPLICIT NONE

!
! !INPUT/OUTPUT PARAMETERS: 
!
      LOGICAL,       INTENT(INOUT) :: error ! Error variable
!
! !REVISION HISTORY: 
!  Apri 2013 - G. Gonzalez Abad - Initial Version
!
!EOP
!BOC

      ! ---------------
      ! Input file name
      ! ---------------
      CHARACTER(LEN=max_ch_len) :: fname
      INTEGER                   :: funit

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

      ! ---------------------------------
      ! Strings in the Control Input file
      ! ---------------------------------
      CHARACTER(LEN=13), PARAMETER :: output_str        = 'Output folder'
      CHARACTER(LEN=18), PARAMETER :: database_str      = 'Database directory'
      CHARACTER(LEN=21), PARAMETER :: profile_str       = 'Profile data filename'
      CHARACTER(LEN=18), PARAMETER :: footprint_str     = 'Use footprint info'
      CHARACTER(LEN=13), PARAMETER :: debug_str         = 'Debug geocape'
      CHARACTER(LEN=13), PARAMETER :: stokes_str        = 'Stokes output'
      CHARACTER(LEN=16), PARAMETER :: airmass_str       = 'Air mass factors'
      CHARACTER(LEN=21), PARAMETER :: tjacobians_str    = 'Temperature jacobians'
      CHARACTER(LEN=26), PARAMETER :: spjacobians_str   = 'Surface pressure jacobians'
      CHARACTER(LEN=20), PARAMETER :: normalizedwf_str  = 'Normalized WF output'
      CHARACTER(LEN=19), PARAMETER :: normalizedrd_str  = 'Normalized radiance'
      CHARACTER(LEN=13), PARAMETER :: solarphotons_str  = 'Solar photons'
      CHARACTER(LEN=24), PARAMETER :: effcross_str      = 'Effective cross sections'
      CHARACTER(LEN= 8), PARAMETER :: spectral_str      = 'Spectral'
      CHARACTER(LEN= 5), PARAMETER :: gases_str         = 'Gases'
      CHARACTER(LEN= 6), PARAMETER :: albedo_str        = 'Albedo'
      CHARACTER(LEN= 8), PARAMETER :: aerosols_str      = 'Aerosols'
      CHARACTER(LEN= 6), PARAMETER :: clouds_str        = 'Clouds'
      CHARACTER(LEN=17), PARAMETER :: moments_str       = 'Number of moments'
      CHARACTER(LEN=13), PARAMETER :: sun_str           = 'Sun positions'
      CHARACTER(LEN=11), PARAMETER :: view_str          = 'View angles'
      CHARACTER(LEN=14), PARAMETER :: azimuth_str       = 'Azimuth angles'
      CHARACTER(LEN=11), PARAMETER :: userlevels_str    = 'User levels'
      CHARACTER(LEN=14), PARAMETER :: useraltitudes_str = 'User altitudes'
      CHARACTER(LEN= 9), PARAMETER :: upwelling_str     = 'Upwelling'

      ! ----------------
      ! Code starts here
      ! ----------------
      error = .FALSE.

      ! -----------------
      ! Control file unit
      ! -----------------
      funit = ctrunit

      ! ---------------------------------------------------
      ! The tool will spect to get the filename as argument
      ! ---------------------------------------------------
      CALL GETARG(1, fname)

      ! -----------------
      ! Open control file
      ! -----------------
      OPEN (UNIT=funit, FILE=TRIM(ADJUSTL(fname)), ACTION='READ', STATUS='OLD', IOSTAT=ios )
      IF ( ios /= 0 ) THEN
         CALL write_err_message ( .TRUE., 'ERROR: Unable to open '//TRIM(ADJUSTL(fname)) )
         CALL error_exit (error)
      END IF

      ! --------------
      ! Results folder
      ! --------------
      REWIND (funit)
      CALL skip_to_filemark ( funit, output_str , tmpstr, error )
      IF ( error ) THEN
         CALL write_err_message ( .TRUE., "ERROR: Can't find "//output_str )
         CALL error_exit (error)
      END IF
      READ (UNIT=funit, FMT='(A)', IOSTAT=ios) results_dir

      ! ---------------
      ! Database folder
      ! ---------------
      REWIND (funit)
      CALL skip_to_filemark ( funit, database_str , tmpstr, error )
      IF ( error ) THEN
         CALL write_err_message ( .TRUE., "ERROR: Can't find "//database_str )
         CALL error_exit (error)
      END IF
      READ (UNIT=funit, FMT='(A)', IOSTAT=ios) database_dir

      ! ---------------------
      ! Profile data filename
      ! ---------------------
      REWIND (funit)
      CALL skip_to_filemark ( funit, profile_str , tmpstr, error )
      IF ( error ) THEN
         CALL write_err_message ( .TRUE., "ERROR: Can't find "//profile_str )
         CALL error_exit (error)
      END IF
      READ (UNIT=funit, FMT='(A)', IOSTAT=ios) profile_data_filename

      ! ------------------
      ! Use footprint info
      ! ------------------
      REWIND (funit)
      CALL skip_to_filemark ( funit, footprint_str , tmpstr, error )
      IF ( error ) THEN
         CALL write_err_message ( .TRUE., "ERROR: Can't find "//footprint_str )
         CALL error_exit (error)
      END IF
      READ (UNIT=funit, FMT=*, IOSTAT=ios) use_footprint_info
      ! ---------
      ! Footprint
      ! ---------
      READ (UNIT=funit, FMT=*, IOSTAT=ios) year, month, day, utc
      READ (UNIT=funit, FMT=*, IOSTAT=ios) longitude, latitude
      READ (UNIT=funit, FMT=*, IOSTAT=ios) satlon, satlat, satalt
      READ (UNIT=funit, FMT=*, IOSTAT=ios) do_sat_viewcalc

      ! -------------
      ! Debug Geocape
      ! -------------
      REWIND (funit)
      CALL skip_to_filemark ( funit, debug_str , tmpstr, error )
      IF ( error ) THEN
         CALL write_err_message ( .TRUE., "ERROR: Can't find "//debug_str )
         CALL error_exit (error)
      END IF
      READ (UNIT=funit, FMT=*, IOSTAT=ios) do_debug_geocape_tool
      READ (UNIT=funit, FMT='(A)', IOSTAT=ios) debug_filename

      ! -------------
      ! Stokes output
      ! -------------
      REWIND (funit)
      CALL skip_to_filemark ( funit, stokes_str , tmpstr, error )
      IF ( error ) THEN
         CALL write_err_message ( .TRUE., "ERROR: Can't find "//stokes_str )
         CALL error_exit (error)
      END IF
      READ (UNIT=funit, FMT=*, IOSTAT=ios) do_StokesQU_output

      ! ---------------
      ! AMF calculation
      ! ---------------
      REWIND (funit)
      CALL skip_to_filemark ( funit, airmass_str , tmpstr, error )
      IF ( error ) THEN
         CALL write_err_message ( .TRUE., "ERROR: Can't find "//airmass_str )
         CALL error_exit (error)
      END IF
      READ (UNIT=funit, FMT=*, IOSTAT=ios) do_AMF_calculation

      ! ---------------------
      ! Temperature jacobians
      ! ---------------------
      REWIND (funit)
      CALL skip_to_filemark ( funit, tjacobians_str , tmpstr, error )
      IF ( error ) THEN
         CALL write_err_message ( .TRUE., "ERROR: Can't find "//tjacobians_str )
         CALL error_exit (error)
      END IF
      READ (UNIT=funit, FMT=*, IOSTAT=ios) do_T_Jacobians

      ! --------------------------
      ! Surface pressure jacobians
      ! --------------------------
      REWIND (funit)
      CALL skip_to_filemark ( funit, spjacobians_str , tmpstr, error )
      IF ( error ) THEN
         CALL write_err_message ( .TRUE., "ERROR: Can't find "//spjacobians_str )
         CALL error_exit (error)
      END IF
      READ (UNIT=funit, FMT=*, IOSTAT=ios) do_sfcprs_Jacobians

      ! --------------------
      ! Normalized WF output
      ! --------------------
      REWIND (funit)
      CALL skip_to_filemark ( funit, normalizedwf_str , tmpstr, error )
      IF ( error ) THEN
         CALL write_err_message ( .TRUE., "ERROR: Can't find "//normalizedwf_str )
         CALL error_exit (error)
      END IF
      READ (UNIT=funit, FMT=*, IOSTAT=ios) do_normalized_WFoutput

      ! -------------------
      ! Normalized radiance
      ! -------------------
      REWIND (funit)
      CALL skip_to_filemark ( funit, normalizedrd_str , tmpstr, error )
      IF ( error ) THEN
         CALL write_err_message ( .TRUE., "ERROR: Can't find "//normalizedrd_str )
         CALL error_exit (error)
      END IF
      READ (UNIT=funit, FMT=*, IOSTAT=ios) do_normalized_radiance
      READ (UNIT=funit, FMT='(A)', IOSTAT=ios) solar_spec_filename

      ! -------------
      ! Solar photons
      ! -------------
      REWIND (funit)
      CALL skip_to_filemark ( funit, solarphotons_str , tmpstr, error )
      IF ( error ) THEN
         CALL write_err_message ( .TRUE., "ERROR: Can't find "//solarphotons_str )
         CALL error_exit (error)
      END IF
      READ (UNIT=funit, FMT=*, IOSTAT=ios) use_solar_photons

      ! ------------------------
      ! Effective cross sections
      ! ------------------------
      REWIND (funit)
      CALL skip_to_filemark ( funit, effcross_str , tmpstr, error )
      IF ( error ) THEN
         CALL write_err_message ( .TRUE., "ERROR: Can't find "//effcross_str )
         CALL error_exit (error)
      END IF
      READ (UNIT=funit, FMT=*, IOSTAT=ios) do_effcrs

      ! --------
      ! Spectral
      ! --------
      REWIND (funit)
      CALL skip_to_filemark ( funit, spectral_str , tmpstr, error )
      IF ( error ) THEN
         CALL write_err_message ( .TRUE., "ERROR: Can't find "//spectral_str )
         CALL error_exit (error)
      END IF
      READ (UNIT=funit, FMT=*, IOSTAT=ios) use_wavelength
      READ (UNIT=funit, FMT=*, IOSTAT=ios) lambda_start, lambda_finish, lambda_resolution
      READ (UNIT=funit, FMT=*, IOSTAT=ios) lambda_dw, lambda_dfw


      ! -----
      ! Gases
      ! -----
      REWIND (funit)
      CALL skip_to_filemark ( funit, gases_str , tmpstr, error )
      IF ( error ) THEN
         CALL write_err_message ( .TRUE., "ERROR: Can't find "//gases_str )
         CALL error_exit (error)
      END IF
      READ (UNIT=funit, FMT=*, IOSTAT=ios) ngases
      READ (UNIT=funit, FMT=*, IOSTAT=ios) (which_gases(i), i=1, ngases)

      ! ------
      ! Albedo
      ! ------
      REWIND (funit)
      CALL skip_to_filemark ( funit, albedo_str , tmpstr, error )
      IF ( error ) THEN
         CALL write_err_message ( .TRUE., "ERROR: Can't find "//albedo_str )
         CALL error_exit (error)
      END IF
      READ (UNIT=funit, FMT=*, IOSTAT=ios) use_lambertian, use_fixedalbedo
      READ (UNIT=funit, FMT=*, IOSTAT=ios) fixed_albedo
      READ (UNIT=funit, FMT=*, IOSTAT=ios) use_albspectra
      READ (UNIT=funit, FMT='(A)', IOSTAT=ios) albspectra_fname
      READ (UNIT=funit, FMT=*, IOSTAT=ios) wind_speed

      ! --------
      ! Aerosols
      ! --------
      REWIND (funit)
      CALL skip_to_filemark ( funit, aerosols_str , tmpstr, error )
      IF ( error ) THEN
         CALL write_err_message ( .TRUE., "ERROR: Can't find "//aerosols_str )
         CALL error_exit (error)
      END IF
      READ (UNIT=funit, FMT=*, IOSTAT=ios) do_aerosols
      READ (UNIT=funit, FMT=*, IOSTAT=ios) use_aerprof
      READ (UNIT=funit, FMT='(A)', IOSTAT=ios) aerfile

      ! ------
      ! Clouds
      ! ------
      REWIND (funit)
      CALL skip_to_filemark ( funit, clouds_str , tmpstr, error )
      IF ( error ) THEN
         CALL write_err_message ( .TRUE., "ERROR: Can't find "//clouds_str )
         CALL error_exit (error)
      END IF
      READ (UNIT=funit, FMT=*, IOSTAT=ios) do_clouds, do_lambertian_cld
      READ (UNIT=funit, FMT=*, IOSTAT=ios) use_cldprof
      READ (UNIT=funit, FMT='(A)', IOSTAT=ios) cldfile

      ! --------------------------
      ! Number of aerosols moments
      ! --------------------------
      REWIND (funit)
      CALL skip_to_filemark ( funit, moments_str , tmpstr, error )
      IF ( error ) THEN
         CALL write_err_message ( .TRUE., "ERROR: Can't find "//moments_str )
         CALL error_exit (error)
      END IF
      READ (UNIT=funit, FMT=*, IOSTAT=ios) aercld_nmoments_input

      ! --------------
      ! User altitudes
      ! --------------
      REWIND (funit)
      CALL skip_to_filemark ( funit, useraltitudes_str , tmpstr, error )
      IF ( error ) THEN
         CALL write_err_message ( .TRUE., "ERROR: Can't find "//useraltitudes_str )
         CALL error_exit (error)
      END IF
      READ (UNIT=funit, FMT=*, IOSTAT=ios) GC_do_user_altitudes
      READ (UNIT=funit, FMT=*, IOSTAT=ios) GC_n_user_altitudes
      READ (UNIT=funit, FMT=*, IOSTAT=ios) (GC_user_altitudes(i), i=1, GC_n_user_altitudes)

      CLOSE(funit)

    END SUBROUTINE read_control_file
!EOC

!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  skip_to_filemark
!
! !DESCRIPTION: This routine does something to the input variable and returns
! the result in the output variable.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE skip_to_filemark ( funit, fmark, lastline, error )

    IMPLICIT NONE
!
! !INPUT PARAMETERS: 
!
    INTEGER,          INTENT(IN)  :: funit ! File unit 
    CHARACTER(LEN=*), INTENT(IN)  :: fmark ! String header
!
! !INPUT/OUTPUT PARAMETERS: 
!
    LOGICAL,          INTENT(INOUT) :: error ! Error variable
!
! !OUTPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(OUT) :: lastline ! Output variable
!
! !REVISION HISTORY: 
!  Apri 2013 - G. Gonzalez Abad - Initial Version
!
!EOP
!------------------------------------------------------------------------------
!BOC
    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER                    :: lmlen, ios
    CHARACTER (LEN=max_ch_len) :: tmpline

    ! -------------------------------------------
    ! Determine the length of the string landmark
    ! -------------------------------------------
    lmlen = LEN(TRIM(ADJUSTL(fmark)))
  
    ! -------------------------------------------------------
    ! Read lines in the file until we either find the string,
    ! reach the end of the file, or reading fails otherwise.
    ! ----------------------------------------------------
    ios = 0
    getlm: DO WHILE ( ios == 0 )
       READ (UNIT=funit, FMT='(A)', IOSTAT=ios) tmpline
       tmpline = TRIM(ADJUSTL(tmpline))
       IF ( ios /= 0 .OR. tmpline(1:lmlen) == fmark ) EXIT getlm
    END DO getlm

    ! ---------------------------------------------------
    ! Return the last line read for the case that we need 
    ! to extract further information from it
    ! ---------------------------------------------------
    lastline = TRIM(ADJUSTL(tmpline))

    ! -----------------------------------------
    ! Set error flag if READ was not successful
    ! -----------------------------------------
    IF ( ios /= 0 ) error = .TRUE.
    IF ( ios ==  0 ) error = .FALSE.

    RETURN
  END SUBROUTINE skip_to_filemark
!EOC

!BOP
!
! !IROUTINE:  check_input
!
! !DESCRIPTION: This routine checks for some incosistencies in the input values and
!               changes values accordingly.
!\\
!\\
! !INTERFACE: check_input(error)
!
  SUBROUTINE check_input(error)

    IMPLICIT NONE
!
! !INPUT/OUTPUT PARAMETERS: 
!
    LOGICAL,       INTENT(INOUT) :: error ! Error variable
!
! !REVISION HISTORY: 
!  Apri 2013 - G. Gonzalez Abad - Initial Version
!
!EOP
!BOC

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER :: i, du, w
    
    error = .FALSE.
      
    DO i = 1, LEN(results_dir)
       IF (results_dir(i:i) == ' ') THEN
          results_dir = results_dir(1:i-1); EXIT
       ENDIF
    ENDDO

    DO i = 1, LEN(database_dir)
       IF (database_dir(i:i) == ' ') THEN
          database_dir = database_dir(1:i-1); EXIT
       ENDIF
    ENDDO

    DO i = 1, LEN(aerfile)
       IF (aerfile(i:i) == ' ') THEN
          aerfile = aerfile(1:i-1); EXIT
       ENDIF
    ENDDO

    DO i = 1, LEN(cldfile)
       IF (cldfile(i:i) == ' ') THEN
          cldfile = cldfile(1:i-1); EXIT
       ENDIF
    ENDDO

    DO i = 1, LEN(profile_data_filename)
       IF (profile_data_filename(i:i) == ' ') THEN
          profile_data_filename = profile_data_filename(1:i-1); EXIT
       ENDIF
    ENDDO

    DO i = 1, LEN(debug_filename)
       IF (debug_filename(i:i) == ' ') THEN
          debug_filename = debug_filename(1:i-1); EXIT
       ENDIF
    ENDDO
    debug_filename = TRIM(ADJUSTL(results_dir)) // debug_filename

    ! ----------------------
    ! Check this input data
    ! ----------------------

    !  No aerosols for a vector calculation
    IF ( do_aerosols .OR. (do_clouds .AND. .NOT. do_lambertian_cld)) THEN

       IF (aercld_nmoments_input .Gt. maxmoms) THEN
          CALL write_err_message ( .FALSE., "Need to increase maxmoms to >= "// &
               "aercld_nmoments_input ")
          error = .TRUE.
       END IF
    END IF

    !  No satellite viewing geometry calculation for downwelling calculation
!!$    IF ( idix(1) .NE. 1 .AND. do_sat_viewcalc) THEN
!!$       CALL write_err_message ( .FALSE., "Cannot calculate viewing geometry for downwelling!!!")
!!$       error = .TRUE.
!!$    END IF
    
    !  Exception handling this section
    IF (error) CALL error_exit(error)
    
    !  set wavelength grid
    IF (lambda_resolution < 0.0 ) THEN
       CALL write_err_message ( .FALSE., "Spectral resolution must be >= 0.0!!!")
       error = .TRUE.
    END IF
    
    IF (lambda_resolution == 0.0 .AND. lambda_dw /= lambda_dfw) THEN
       CALL write_err_message ( .FALSE., "Fine and coarse grids must be the same "// &
                                        "for high-resolution calculation!!!")
       error = .TRUE.
    END IF
    
    IF (lambda_resolution == 0.0) THEN
       do_effcrs = .FALSE.
    END IF
    
    nclambdas = NINT ( (lambda_finish - lambda_start ) / lambda_dw ) + 1
    clambdas(1) = lambda_start
    DO w = 2, nclambdas
       clambdas(w) = clambdas(w-1) + lambda_dw
    END DO
    
    !  For high-resolution grid, need to have extra wavelength space for convolution
    edge_dw = 2.0 * lambda_resolution  
    nflambdas = NINT ( (lambda_finish - lambda_start + 2.0 * edge_dw) / lambda_dfw ) + 1
    flambdas(1) = lambda_start - edge_dw
    DO w = 2, nflambdas
       flambdas(w) = flambdas(w-1) + lambda_dfw
    END DO
    
    IF (.NOT. use_wavelength) THEN
       cwavnums(1:nclambdas) = clambdas(1:nclambdas)
       fwavnums(1:nflambdas) = flambdas(1:nflambdas)
       
       clambdas(1:nclambdas) = 1.0D+07 / cwavnums(1:nclambdas) 
       flambdas(1:nflambdas) = 1.0D+07 / fwavnums(1:nflambdas) 
    ELSE
       cwavnums(1:nclambdas) = 1.0D+07 / clambdas(1:nclambdas) 
       fwavnums(1:nflambdas) = 1.0D+07 / flambdas(1:nflambdas) 
    END IF
    
    ! Fine grid will be used for reading solar irradiance spectrum, and cross sections
    ! no matter whether do_effcrs is true or false
    nlambdas = nflambdas
    lambdas(1:nlambdas) = flambdas(1:nflambdas)
    wavnums(1:nlambdas) = fwavnums(1:nflambdas)
    
    IF (maxlambdas > maxflambdas) THEN
       CALL write_err_message ( .FALSE., "maxlambdas sould be <= maxflambdas!")
       error = .TRUE.
    ENDIF
    
    IF (.NOT. do_effcrs .AND. maxlambdas /= maxflambdas) THEN
       CALL write_err_message ( .FALSE., "Need to set maxlambdas = maxflambdas "// &
            "for fine-grid calculation!")
       error = .TRUE.
    END IF
    
    IF ( nclambdas > maxlambdas ) THEN
       CALL write_err_message ( .FALSE., "Need to increase maxlambdas!")
       error = .TRUE.
    END IF
    
    IF ( nflambdas > maxflambdas ) THEN
       CALL write_err_message ( .FALSE., "Need to increase maxflambdas!")
       error = .TRUE.
   END IF
   
   !  Exception handling this section
   IF (error) CALL error_exit(error)
   
   !  Open debug file if flagged. The unit number is 'du'
   IF ( do_debug_geocape_tool ) THEN
      du=66
      OPEN(du, FILE=debug_filename, STATUS = 'unknown')
   END IF
   
 END SUBROUTINE check_input
 !EOC

END MODULE GC_read_input_module
