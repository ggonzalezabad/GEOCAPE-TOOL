MODULE GC_solar_module

  USE GC_parameters_module, ONLY: maxflambdas
  USE GC_variables_module,  ONLY: messages, nmessages, hour, minute, &
                                  second, nspec, use_wavelength,     &
                                  lambdas, nlambdas, database_dir,   &
                                  use_solar_photons, wavnums, month, &
                                  day, year, lambda_resolution,      &
                                  clambdas, nclambdas, cwavnums,     &
                                  jday, jday0, rdis, julday,         &
                                  solar_spec_data, solar_cspec_data, &
                                  solar_spec_filename
  USE GC_error_module

  IMPLICIT NONE

  CONTAINS

    SUBROUTINE prepare_solar_spectrum(error)

      IMPLICIT NONE

      LOGICAL, INTENT(INOUT) :: error

      ! ----------------
      ! Code starts here
      ! ----------------
      error = .FALSE.

      ! ------------------------------------------------------------------------
      ! Flux Factor
      !   -- Set to 1.0 if you want sun-normalized only (solar, NO thermal)
      !   -- Set to PHYSICAL UNIT W/sqm (Same as BLackBody) with Solar + Thermal
      !   -- Not used for thermal only options
      ! ------------------------------------------------------------------------

      !   IF (do_normalized_radiance) THEN
      !      solar_spec_data(1:nlambdas)  = 1.0d0
      !      solar_cspec_data(1:nclambdas) = 1.0d0
      !   ELSE
      

      ! ------------------------
      ! Decide which file to use
      ! The file always should
      ! have units of W/cm2/cm-1
      ! ------------------------
      IF (use_wavelength .AND. lambdas(1) > 201.d0 .AND. &
         lambdas(nlambdas) < 1000.d0) THEN
         solar_spec_filename = TRIM(ADJUSTL(database_dir)) // &
                               '/chance_solarspec_jqsrt2011_bis.dat'
      ELSE
         solar_spec_filename = TRIM(ADJUSTL(database_dir)) // &
                               '/newkur.dat'
      END IF

      ! ---------
      ! Read file
      ! ---------
      CALL solarspec_prep_newkur &
           ( solar_spec_filename, maxflambdas, nlambdas, wavnums,  & ! input
           solar_spec_data, messages(nmessages+1), error )        ! output in W / m ^2 / cm^-1
               
      ! ----------------------------
      ! Convert to photons if needed
      ! ----------------------------
      IF ( use_solar_photons ) THEN  ! photons/cm^2/nm/s
         solar_spec_data(1:nlambdas) = 5.0341174716D+18 * solar_spec_data(1:nlambdas) &
              / lambdas(1:nlambdas)
      END IF
      
      ! ------------------------------
      ! Correct for sun-earth distance
      ! ------------------------------
      jday  = julday ( month, day, year, hour, minute, second)
      jday0 = julday (     1,   1, year,    0,      0,  0.0d0)
      jday  = jday - jday0
      CALL sun_earth_distance (jday, rdis)
      solar_spec_data(1:nlambdas) = solar_spec_data(1:nlambdas) / rdis / rdis
      
      ! -------------------------------
      ! Excpetion handling this section
      ! -------------------------------
      IF (error) CALL write_err_message ( .TRUE., messages(nmessages+1))

      ! -----------------------------------------------------
      ! Get solar irradiance spectrum at coarse spectral grid
      ! -----------------------------------------------------
      nspec = 1
      IF (lambda_resolution > 0.0) THEN
         IF (use_wavelength) THEN
            CALL gauss_f2c (lambdas(1:nlambdas), solar_spec_data(1:nlambdas), &
                            nlambdas, nspec, lambda_resolution, clambdas,     &
                            solar_cspec_data(1:nclambdas), nclambdas)
         ELSE
            CALL gauss_f2c (wavnums(1:nlambdas), solar_spec_data(1:nlambdas), &
                            nlambdas, nspec, lambda_resolution, cwavnums,     &
                            solar_cspec_data(1:nclambdas), nclambdas)
         END IF
      ELSE
         solar_cspec_data(1:nlambdas) = solar_spec_data(1:nlambdas)
      END IF

    END SUBROUTINE prepare_solar_spectrum

END MODULE GC_solar_module
