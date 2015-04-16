MODULE GC_xsections_module

  USE GC_parameters_module, ONLY: maxgases, maxflambdas, maxgases,   &
                                  GC_maxlayers, ctrunit, max_ch_len, &
                                  maxlambdas
  USE GC_variables_module,  ONLY: nmessages, messages, tmpwaves,       &
                                  GC_nlayers, mid_pressures,           &
                                  mid_temperatures, tmpcwaves,         &
                                  gas_totalcolumns, aircolumns,        &
                                  gas_partialcolumns,                  &
                                  nspec, ngases, database_dir,         &
                                  which_gases, lambdas, nlambdas,      &
                                  use_wavelength, do_effcrs, clambdas, &
                                  nclambdas, wavnums, cwavnums,        &
                                  lambda_resolution, solar_spec_data,  &
                                  xsec_data_filenames, xsec_data_path, &
                                  gas_nxsecs, gas_xsecs_type,          &
                                  gas_xsecs, o3c1_xsecs, o3c2_xsecs,   &
                                  Rayleigh_xsecs, Rayleigh_depols
  USE GC_read_input_module, ONLY: skip_to_filemark
  USE GC_error_module

  IMPLICIT NONE

  CONTAINS

    SUBROUTINE read_xsec_filenames(error)

      IMPLICIT NONE

      LOGICAL, INTENT(INOUT) :: error

      CHARACTER(LEN=14), PARAMETER :: xsec_str = 'Cross sections'

      CHARACTER(LEN=max_ch_len) :: fname
      CHARACTER(LEN=max_ch_len) :: tmpstr, tmpxsecfname, tmpname

      INTEGER :: funit, igas, tmpnumberxsec, tmpxsectype, ios
      LOGICAL :: noxsec

      ! -----------
      ! Code starts
      ! -----------
      error = .FALSE.
      noxsec   = .FALSE.

      ! ---------------------
      ! Trace gas source path
      ! ---------------------
      xsec_data_path = TRIM(ADJUSTL(database_dir)) // '/SAO_crosssections/'

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
      END IF

      ! ----------------------------
      ! Move to cross sections input
      ! ----------------------------
      DO igas = 1, ngases
         REWIND (funit)
         CALL skip_to_filemark ( funit, xsec_str , tmpstr, error )
         IF ( error ) THEN
            CALL write_err_message ( .TRUE., "ERROR: Can't find "//xsec_str )
            CALL error_exit (error)
         END IF
         DO
            READ(UNIT=funit, FMT=*, iostat=ios) tmpname
            IF (ios .GT. 0) THEN
               CALL write_err_message(.TRUE., "Error reading x-sec info from control file")
               noxsec = .TRUE.
               CALL error_exit(error)
            ELSE IF (ios .LT. 0) THEN
               CALL write_err_message(.TRUE., "End of input control file reached without "// &
                                              "finding x-sec for "//which_gases(igas))
               noxsec = .TRUE.
               EXIT
            ELSE
               READ(UNIT=funit, FMT=*) tmpxsecfname
               READ(UNIT=funit, FMT=*) tmpnumberxsec, tmpxsectype
               IF (TRIM(ADJUSTL(which_gases(igas))) .EQ. TRIM(ADJUSTL(tmpname))) THEN
                  xsec_data_filenames(igas) = TRIM(ADJUSTL(tmpxsecfname))
                  gas_nxsecs(igas)          = tmpnumberxsec
                  gas_xsecs_type(igas)      = tmpxsectype
                  EXIT
               END IF
            END IF
         END DO
      END DO

      error = noxsec

    END SUBROUTINE read_xsec_filenames

    SUBROUTINE prepare_xsec(error)

      IMPLICIT NONE

      LOGICAL :: error

      ! ---------------
      ! Local variables
      ! ---------------
      REAL(KIND=8)                                         :: scalex
      REAL(KIND=8),    DIMENSION(maxlambdas, GC_maxlayers) :: temp_xsecs
      INTEGER(KIND=4), DIMENSION(maxflambdas)              :: rvidxs, rcvidxs ! For wavenumbers

      CHARACTER(max_ch_len) :: tmpchar
      INTEGER :: g
      
      ! -------------------------------------------------------------------
      ! Rayleigh calculation is based on the Bodhaine et al., (1999) paper.
      ! -------------------------------------------------------------------
      ! -----------
      ! Code starts
      ! -----------
      error = .FALSE.

      ! -------------------------------------------------------------------------
      ! Call to get the cross-sections
      !   -- 2 additional entries O3 (quadratic parameterization of T-dependence)
      !   -- All Xsec output in [cm^2/mol]
      ! -------------------------------------------------------------------------
      tmpwaves(1:nlambdas) = lambdas(1:nlambdas)
      tmpcwaves(1:nclambdas) = clambdas(1:nclambdas)

      IF (.NOT. use_wavelength) THEN
         CALL reverse(tmpwaves(1:nlambdas), nlambdas)
         CALL reverse(tmpcwaves(1:nclambdas), nclambdas)
      ENDIF

      CALL geocape_xsec_setter_1                                                       &     
           ( maxflambdas, maxgases, GC_maxlayers, xsec_data_path, xsec_data_filenames, & ! Input
           nlambdas, tmpwaves, nclambdas, tmpcwaves, ngases, which_gases,              & ! Input
           gas_partialcolumns, gas_xsecs_type, GC_nlayers, mid_pressures,              & ! Input
           mid_temperatures, do_effcrs, use_wavelength, lambda_resolution,             & ! Input
           solar_spec_data,                                                            & ! Input
           gas_xsecs, o3c1_xsecs, o3c2_xsecs, Rayleigh_xsecs, Rayleigh_depols,         & ! Output
           tmpchar, error )                                                              ! Errors

      IF (.NOT. use_wavelength) THEN

         CALL reverse_idxs(nlambdas, rvidxs(1:nlambdas))
         IF ( ANY(gas_xsecs_type(1:ngases) .eq. 3) .AND. do_effcrs ) &
              CALL reverse_idxs(nclambdas, rcvidxs(1:nclambdas))

         DO g = 1, ngases 
            IF (gas_xsecs_type(g) .eq. 3 .AND. trim(adjustl(xsec_data_filenames(g))) &
                 .eq. 'HITRAN' .AND. do_effcrs) THEN ! Already convolved in get_hitran_crs
               gas_xsecs(1:nclambdas, 1:GC_nlayers, g) = &
                    gas_xsecs(rcvidxs(1:nclambdas), 1:GC_nlayers, g) 
            ELSE
               gas_xsecs(1:nlambdas, 1:GC_nlayers, g) = &
                    gas_xsecs(rvidxs(1:nlambdas), 1:GC_nlayers, g) 
            ENDIF
         ENDDO

         o3c1_xsecs(1:nlambdas)      = o3c1_xsecs(rvidxs(1:nlambdas))
         o3c2_xsecs(1:nlambdas)      = o3c2_xsecs(rvidxs(1:nlambdas))
         Rayleigh_xsecs(1:nlambdas)  = Rayleigh_xsecs(rvidxs(1:nlambdas)) 
         Rayleigh_depols(1:nlambdas) = Rayleigh_depols(rvidxs(1:nlambdas)) 

      END IF
  
      ! --------------------
      ! If there is an error
      ! --------------------
      IF (error) CALL write_err_message (.FALSE., tmpchar)

      ! -------------------------------------------------------------------------
      ! If using effective cross sections, need to convolve high-resolution cross
      ! sections with slit functions
      ! -------------------------------------------------------------------------
      IF (do_effcrs) THEN

         IF (use_wavelength) THEN
            tmpwaves (1:nlambdas)  = lambdas (1:nlambdas)
            tmpcwaves(1:nclambdas) = clambdas(1:nclambdas)
         ELSE
            tmpwaves (1:nlambdas)  = wavnums (1:nlambdas)
            tmpcwaves(1:nclambdas) = cwavnums(1:nclambdas)
         END IF
      
         DO g = 1, ngases

            IF (gas_xsecs_type(g) .ne. 3 .and. trim(adjustl(xsec_data_filenames(g))) .ne. 'HITRAN' ) THEN
               nspec  = gas_nxsecs(g)
               scalex = gas_totalcolumns(g) * 2.68668D16 * 2.0D0
               CALL gauss_f2ci0(tmpwaves(1:nlambdas), gas_xsecs(1:nlambdas, 1:nspec, g), &
                    solar_spec_data(1:nlambdas), nlambdas, nspec, scalex,    &
                    lambda_resolution, tmpcwaves(1:nclambdas),               &
                    temp_xsecs(1:nclambdas, 1:nspec), nclambdas)
               
               gas_xsecs(1:nclambdas, 1:nspec, g) = temp_xsecs(1:nclambdas, 1:nspec)
            ENDIF

            IF (which_gases(g) .eq. 'O3  ' ) THEN
               o3c1_xsecs(1:nclambdas) = gas_xsecs(1:nclambdas, 2, g)
               o3c2_xsecs(1:nclambdas) = gas_xsecs(1:nclambdas, 3, g)
            END IF
         ENDDO

         nspec = 1
         scalex = SUM(aircolumns(1:GC_nlayers)) * 2.0
         CALL gauss_f2ci0(tmpwaves(1:nlambdas), Rayleigh_xsecs(1:nlambdas),     &
                          solar_spec_data(1:nlambdas), nlambdas, nspec, scalex, &
                          lambda_resolution, tmpcwaves(1:nclambdas),            &
                          temp_xsecs(1:nclambdas, 1:1), nclambdas)
         
         rayleigh_xsecs(1:nclambdas) = temp_xsecs(1:nclambdas, 1)
      
         scalex = 1.0
         CALL gauss_f2ci0(tmpwaves(1:nlambdas), Rayleigh_depols(1:nlambdas),    &
                          solar_spec_data(1:nlambdas), nlambdas, nspec, scalex, &
                          lambda_resolution, tmpcwaves(1:nclambdas),            &
                          temp_xsecs(1:nclambdas, 1:1), nclambdas)

         Rayleigh_depols(1:nclambdas) = temp_xsecs(1:nclambdas, 1)
      
         nlambdas = nclambdas
         lambdas(1:nlambdas) = clambdas(1:nclambdas)
         wavnums(1:nlambdas) = cwavnums(1:nclambdas)
      END IF

    END SUBROUTINE prepare_xsec

END MODULE GC_xsections_module
