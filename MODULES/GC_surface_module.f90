MODULE GC_surface_module

  USE GC_parameters_module, ONLY: maxlambdas
  USE GC_variables_module,  ONLY: tmpwaves, surface_data_path, do_clouds, &
                                  do_lambertian_cld, use_lambertian,      &
                                  nlambdas, use_fixedalbedo, fixed_albedo,&
                                  use_wavelength, lambdas, use_albspectra,&
                                  albspectra_fname, database_dir, cfrac,  &
                                  lambertian_cldalb, messages, nmessages, &
                                  latitude, longitude, month, ground_ler, &
                                  wavlens_um, water_rn, water_cn
  USE GC_error_module

  IMPLICIT NONE

  CONTAINS

  SUBROUTINE prepare_surface(error)

    IMPLICIT NONE

    LOGICAL, INTENT(INOUT) :: error

    ! ----------------
    ! Code starts here
    ! ----------------
    ! -----------------------------------------------------------
    ! Only taking into account the Lambertian part of the surface
    ! the BRDF option is done in geocape_tool_v2p6 using the BRDF
    ! supplement.
    ! -----------------------------------------------------------
!!$    IF (do_clouds .AND. do_lambertian_cld .AND. cfrac .GT. 0.0d0 &
!!$         .AND. .NOT. use_lambertian) use_lambertian = .TRUE.
    
    IF (use_lambertian) THEN
       IF (do_clouds .AND. do_lambertian_cld .AND. cfrac .GE. 1.0d0) THEN 
          ground_ler(1:nlambdas) = lambertian_cldalb
       ELSE IF (use_fixedalbedo) THEN  ! Used fixed albedo from input file
          ground_ler(1:nlambdas) = fixed_albedo
       ELSE   
          tmpwaves(1:nlambdas) = lambdas(1:nlambdas)
          IF (.NOT. use_wavelength) CALL reverse(tmpwaves(1:nlambdas), nlambdas)
          IF (use_albspectra) THEN     ! Used albedo reflectance spectra
             
             albspectra_fname = TRIM(ADJUSTL(database_dir)) //'/ReflSpectra/'&
                  //TRIM(ADJUSTL(albspectra_fname))
             CALL geocape_surface_setter_2                             &
                  ( albspectra_fname, maxlambdas, nlambdas, tmpwaves,  & ! inputs
                  ground_ler,  messages(nmessages+1), error ) ! outputs, exception handling
          ELSE                         ! Use TEMIS data
             surface_data_path = TRIM(ADJUSTL(database_dir)) // '/TEMIS_reflectances/TEMIS/' 
             CALL geocape_surface_setter_1 &
                  ( surface_data_path, latitude, longitude, month, & ! inputs
                  maxlambdas, nlambdas, tmpwaves,                  & ! inputs
                  ground_ler,                                      & ! Output
                  messages(nmessages+1), error )                     ! Exception handling
          END IF
          IF (.NOT. use_wavelength) CALL reverse(ground_ler(1:nlambdas), nlambdas)
       END IF
    ELSE
       wavlens_um(1:nlambdas) = lambdas(1:nlambdas) * 1.0D-3
       IF (.NOT. use_wavelength) CALL reverse(wavlens_um(1:nlambdas), nlambdas)
       CALL index_water(wavlens_um(1:nlambdas), water_rn(1:nlambdas), &
            water_cn(1:nlambdas), nlambdas)
       IF (.NOT. use_wavelength) THEN
          CALL reverse(wavlens_um(1:nlambdas), nlambdas)
          CALL reverse(water_rn(1:nlambdas), nlambdas)
          CALL reverse(water_cn(1:nlambdas), nlambdas)
       END IF
    END IF
    
    ! --------------------
    ! If there is an error
    ! --------------------
    IF (error) CALL write_err_message(.FALSE., messages(nmessages+1))
    
  END SUBROUTINE prepare_surface

END MODULE GC_surface_module
