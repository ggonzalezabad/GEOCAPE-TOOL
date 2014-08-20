MODULE GC_netcdf_module

  USE GC_parameters_module, ONLY: max_ch_len
  USE GC_variables_module,  ONLY: latitude, longitude, year, month, day, utc, VLIDORT_ModIn, nlambdas, &
                                  lambdas, GC_nlayers, heights, pressures, temperatures, ngases,       &
                                  which_gases, gas_partialcolumns, aircolumns, cfrac, cld_tops,        &
                                  lambertian_cldalb, taer_profile, tcld_profile, opdeps, ssalbs,       &
                                  aer_opdeps, aer_ssalbs, cld_opdeps, cld_ssalbs, ground_ler,          &
                                  wind_speed, VLIDORT_FixIn, aercld_nmoments_input, lambda_dw,         &
                                  lambda_dfw, use_wavelength, aer_reflambda, cld_reflambda,            &
                                  do_vector_calculation, do_StokesQU_output, do_Jacobians,             &
                                  do_QU_Jacobians, do_AMF_calculation, do_T_Jacobians,                 &
                                  do_sfcprs_Jacobians, do_aod_Jacobians, do_assa_Jacobians,            &
                                  do_cod_Jacobians, do_cssa_Jacobians, do_cfrac_Jacobians,             &
                                  do_aer_columnwf, do_cld_columnwf, do_normalized_WFoutput,            &
                                  do_normalized_radiance, use_lambertian, do_lambertian_cld, do_effcrs,&
                                  use_solar_photons, lambda_resolution, solar_cspec_data, GC_radiances,&
                                  GC_Qvalues, GC_Uvalues, GC_Tracegas_Jacobians, GC_Scattering_Weights,&
                                  GC_AMFs, GC_Tracegas_QJacobians, GC_Tracegas_UJacobians,             &
                                  GC_Temperature_Jacobians, GC_Surfalbedo_Jacobians,                   &
                                  GC_Surfalbedo_QJacobians, GC_Surfalbedo_UJacobians,                  &
                                  GC_Windspeed_Jacobians, GC_WindSpeed_QJacobians,                     &
                                  GC_Windspeed_UJacobians, GC_aod_Jacobians, GC_aod_QJacobians,        &
                                  GC_aod_UJacobians, GC_assa_Jacobians, GC_assa_QJacobians,            &
                                  GC_assa_UJacobians, GC_cod_Jacobians, GC_cod_QJacobians,             &
                                  GC_cod_UJacobians, GC_cssa_Jacobians, GC_cssa_QJacobians,            &
                                  GC_cssa_UJacobians, GC_cfrac_Jacobians, GC_cfrac_QJacobians,         &
                                  GC_cfrac_UJacobians, GC_sfcprs_Jacobians, GC_sfcprs_QJacobians,      &
                                  GC_sfcprs_UJacobians, results_dir, GC_n_sun_positions,               &
                                  GC_n_view_angles, GC_n_azimuths, VLIDORT_Out,                        &
                                  GC_flux, GC_Qflux, GC_Uflux, GC_direct_flux, GC_Qdirect_flux,        &
                                  GC_Udirect_flux, Total_brdf, NSTOKESSQ, do_brdf_surface,             &
                                  OUTPUT_WSABSA, WSA_CALCULATED, BSA_CALCULATED, didx
  USE GC_error_module

  IMPLICIT NONE

CONTAINS

  SUBROUTINE netcdf_output(error)

    IMPLICIT NONE

    ! ------------------
    ! Modified variables
    ! ------------------
    LOGICAL,       INTENT(INOUT) :: error ! Error variable

    ! ---------------
    ! Local variables
    ! ---------------
    CHARACTER(LEN=max_ch_len) :: tmperror
    CHARACTER(LEN=max_ch_len) :: netfname

    ! ----------------
    ! Code starts here
    ! ----------------
    error = .FALSE.

    ! -----------------
    ! NETCDF file write
    ! -----------------
    IF (didx .EQ. 1) THEN
       netfname = TRIM(ADJUSTL(TRIM(ADJUSTL(results_dir)) // 'GC_upwelling_output.nc'))
    ELSEIF (didx .EQ. 2) THEN
       netfname = TRIM(ADJUSTL(TRIM(ADJUSTL(results_dir)) // 'GC_downwelling_output.nc'))
    ELSE
       WRITE(*,*) 'No upwelling or downwelling selected... no output'
       STOP
    END IF
    call netcdf_wrt ( netfname, latitude, longitude, year, month, day, utc,              &
        GC_n_sun_positions, GC_n_view_angles, GC_n_azimuths,                             &
        VLIDORT_Out%Main%TS_N_GEOMETRIES,                                                &
        VLIDORT_ModIn%MSunrays%TS_SZANGLES(1:GC_n_sun_positions),                        &
        VLIDORT_ModIn%MUserVal%TS_USER_VZANGLES_INPUT(1:GC_n_view_angles),               &
        VLIDORT_ModIn%MUserVal%TS_USER_RELAZMS(1:GC_n_azimuths),                         &
        VLIDORT_ModIn%MUserVal%TS_USER_LEVELS(1), nlambdas,                              &
        lambdas(1:nlambdas), GC_nlayers, heights(0:GC_nlayers), pressures(0:GC_nlayers), &
        temperatures(0:GC_nlayers), ngases, which_gases(1:ngases),                       &
        gas_partialcolumns(1:GC_nlayers, 1:ngases), aircolumns(1:GC_nlayers),  cfrac,    &
        cld_tops(1), lambertian_cldalb, taer_profile(1:GC_nlayers), tcld_profile(1:GC_nlayers), &
        opdeps(1:nlambdas, 1:GC_nlayers), ssalbs(1:nlambdas, 1:GC_nlayers),              &
        aer_opdeps(1:nlambdas, 1:GC_nlayers), aer_ssalbs(1:nlambdas, 1:GC_nlayers),      &
        cld_opdeps(1:nlambdas, 1:GC_nlayers), cld_ssalbs(1:nlambdas, 1:GC_nlayers),      &
        ground_ler(1:nlambdas), wind_speed, VLIDORT_FixIn%Cont%TS_NSTREAMS,              &
        aercld_nmoments_input, lambda_dw,                                                &
        lambda_dfw, use_wavelength, aer_reflambda, cld_reflambda, do_vector_calculation, &
        do_StokesQU_output, do_Jacobians, do_QU_Jacobians, do_AMF_calculation,           &
        do_T_Jacobians, do_sfcprs_Jacobians, do_aod_Jacobians, do_assa_Jacobians,        &
        do_cod_Jacobians, do_cssa_Jacobians, do_cfrac_Jacobians, do_aer_columnwf,        &
        do_cld_columnwf, do_normalized_WFoutput, do_normalized_radiance,use_lambertian,  &
        do_lambertian_cld, VLIDORT_FixIn%Bool%TS_DO_UPWELLING, do_effcrs,                &
        use_solar_photons, lambda_resolution,                                            &
        solar_cspec_data(1:nlambdas), GC_radiances(1:nlambdas, 1:VLIDORT_Out%Main%TS_N_GEOMETRIES),     &
        GC_Qvalues(1:nlambdas, 1:VLIDORT_Out%Main%TS_N_GEOMETRIES),                                     &
        GC_Uvalues(1:nlambdas, 1:VLIDORT_Out%Main%TS_N_GEOMETRIES),                                     &
        GC_Tracegas_Jacobians(1:nlambdas, 1:GC_nlayers, 1:VLIDORT_Out%Main%TS_N_GEOMETRIES, 1:ngases),  &
        GC_Scattering_Weights(1:nlambdas, 1:GC_nlayers, 1:VLIDORT_Out%Main%TS_N_GEOMETRIES),            &
        GC_AMFs(1:nlambdas, 1:VLIDORT_Out%Main%TS_N_GEOMETRIES, 1:ngases),                              &
        GC_Tracegas_QJacobians(1:nlambdas, 1:GC_nlayers, 1:VLIDORT_Out%Main%TS_N_GEOMETRIES, 1:ngases), &
        GC_Tracegas_UJacobians(1:nlambdas, 1:GC_nlayers, 1:VLIDORT_Out%Main%TS_N_GEOMETRIES, 1:ngases), &
        GC_Temperature_Jacobians(1:nlambdas, 1:GC_nlayers, 1:VLIDORT_Out%Main%TS_N_GEOMETRIES),         &
        GC_Surfalbedo_Jacobians(1:nlambdas, 1:VLIDORT_Out%Main%TS_N_GEOMETRIES),                        &
        GC_Surfalbedo_QJacobians(1:nlambdas, 1:VLIDORT_Out%Main%TS_N_GEOMETRIES),                       &
        GC_Surfalbedo_UJacobians(1:nlambdas, 1:VLIDORT_Out%Main%TS_N_GEOMETRIES),                       &
        GC_Windspeed_Jacobians(1:nlambdas, 1:VLIDORT_Out%Main%TS_N_GEOMETRIES),                         &
        GC_Windspeed_QJacobians(1:nlambdas, 1:VLIDORT_Out%Main%TS_N_GEOMETRIES),                        &
        GC_Windspeed_UJacobians(1:nlambdas, 1:VLIDORT_Out%Main%TS_N_GEOMETRIES),                        &
        GC_aod_Jacobians(1:nlambdas, 1:GC_nlayers, 1:VLIDORT_Out%Main%TS_N_GEOMETRIES),                 &
        GC_aod_QJacobians(1:nlambdas, 1:GC_nlayers, 1:VLIDORT_Out%Main%TS_N_GEOMETRIES),                &
        GC_aod_UJacobians(1:nlambdas, 1:GC_nlayers, 1:VLIDORT_Out%Main%TS_N_GEOMETRIES),                &
        GC_assa_Jacobians(1:nlambdas, 1:GC_nlayers, 1:VLIDORT_Out%Main%TS_N_GEOMETRIES),                &
        GC_assa_QJacobians(1:nlambdas, 1:GC_nlayers, 1:VLIDORT_Out%Main%TS_N_GEOMETRIES),               &
        GC_assa_UJacobians(1:nlambdas, 1:GC_nlayers, 1:VLIDORT_Out%Main%TS_N_GEOMETRIES),               &
        GC_cod_Jacobians(1:nlambdas, 1:GC_nlayers, 1:VLIDORT_Out%Main%TS_N_GEOMETRIES),                 &
        GC_cod_QJacobians(1:nlambdas, 1:GC_nlayers, 1:VLIDORT_Out%Main%TS_N_GEOMETRIES),                &
        GC_cod_UJacobians(1:nlambdas, 1:GC_nlayers, 1:VLIDORT_Out%Main%TS_N_GEOMETRIES),                &
        GC_cssa_Jacobians(1:nlambdas, 1:GC_nlayers, 1:VLIDORT_Out%Main%TS_N_GEOMETRIES),                &
        GC_cssa_QJacobians(1:nlambdas, 1:GC_nlayers, 1:VLIDORT_Out%Main%TS_N_GEOMETRIES),               &
        GC_cssa_UJacobians(1:nlambdas, 1:GC_nlayers, 1:VLIDORT_Out%Main%TS_N_GEOMETRIES),               &
        GC_cfrac_Jacobians(1:nlambdas, 1:VLIDORT_Out%Main%TS_N_GEOMETRIES),                             &
        GC_cfrac_QJacobians(1:nlambdas, 1:VLIDORT_Out%Main%TS_N_GEOMETRIES),                            &
        GC_cfrac_UJacobians(1:nlambdas, 1:VLIDORT_Out%Main%TS_N_GEOMETRIES),                            &
        GC_sfcprs_Jacobians(1:nlambdas, 1:VLIDORT_Out%Main%TS_N_GEOMETRIES),                            &
        GC_sfcprs_QJacobians(1:nlambdas, 1:VLIDORT_Out%Main%TS_N_GEOMETRIES),                           &
        GC_sfcprs_UJacobians(1:nlambdas, 1:VLIDORT_Out%Main%TS_N_GEOMETRIES),                           &
        GC_flux(1:nlambdas, 1:VLIDORT_ModIn%MSunrays%TS_N_SZANGLES),                                    &
        GC_Qflux(1:nlambdas, 1:VLIDORT_ModIn%MSunrays%TS_N_SZANGLES),                                   &
        GC_Uflux(1:nlambdas, 1:VLIDORT_ModIn%MSunrays%TS_N_SZANGLES),                                   &
        GC_direct_flux(1:nlambdas, 1:VLIDORT_ModIn%MSunrays%TS_N_SZANGLES),                             &
        GC_Qdirect_flux(1:nlambdas, 1:VLIDORT_ModIn%MSunrays%TS_N_SZANGLES),                            &
        GC_Udirect_flux(1:nlambdas, 1:VLIDORT_ModIn%MSunrays%TS_N_SZANGLES),                            &
        Total_brdf(1:NSTOKESSQ,1:GC_n_view_angles,1:GC_n_azimuths,1:GC_n_sun_positions), NSTOKESSQ,     &
        do_brdf_surface, OUTPUT_WSABSA, WSA_CALCULATED, BSA_CALCULATED, tmperror, error)

    IF (error) CALL write_err_message (.TRUE., tmperror)

  END SUBROUTINE netcdf_output
  
END MODULE GC_netcdf_module
