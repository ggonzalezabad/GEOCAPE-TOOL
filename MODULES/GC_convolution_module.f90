MODULE GC_convolution_module

  USE GC_variables_module, ONLY: do_effcrs, lambda_resolution, use_wavelength, &
                                 tmpwaves, flambdas, nflambdas, tmpcwaves,     &
                                 clambdas, nclambdas, fwavnums, cwavnums,      &
                                 nspec, VLIDORT_Out, GC_radiances, temp_rad,   &
                                 do_vector_calculation, do_StokesQU_output,    &
                                 GC_Qvalues, GC_Uvalues, ground_ler,           &
                                 do_Jacobians, ngases, GC_Tracegas_Jacobians,  &
                                 VLIDORT_FixIn, temp_wf, do_AMF_calculation,   &
                                 GC_Scattering_Weights, GC_AMFs,               &
                                 GC_aod_Jacobians, &
                                 GC_assa_jacobians, do_cld_columnwf,           &
                                 do_cod_jacobians, GC_cod_Jacobians,           &
                                 do_cssa_jacobians, GC_cssa_Jacobians,         &
                                 do_cfrac_Jacobians, GC_cfrac_Jacobians,       &
                                 GC_Surfalbedo_Jacobians,                      &
                                 GC_Windspeed_Jacobians, GC_sfcprs_Jacobians,  &
                                 do_QU_Jacobians, GC_Tracegas_QJacobians,      &
                                 GC_Tracegas_UJacobians, GC_aod_QJacobians,    &
                                 GC_aod_UJacobians, GC_assa_QJacobians,        &
                                 GC_assa_UJacobians, GC_cod_QJacobians,        &
                                 GC_cod_UJacobians, GC_cssa_QJacobians,        &
                                 GC_cssa_UJacobians, GC_cfrac_QJacobians,      &
                                 GC_cfrac_UJacobians, GC_Surfalbedo_QJacobians,&
                                 GC_Surfalbedo_UJacobians,                     &
                                 GC_Windspeed_QJacobians,                      &
                                 GC_Windspeed_UJacobians, GC_sfcprs_QJacobians,&
                                 GC_sfcprs_UJacobians, do_T_Jacobians,         &
                                 GC_Temperature_Jacobians, nlambdas, lambdas,  &
                                 clambdas, wavnums, GC_flux, GC_Uflux,         &
                                 GC_Qflux, GC_direct_flux, GC_Qdirect_flux,    &
                                 GC_Udirect_flux, VLIDORT_ModIn, didx, ilev, aer_ctr
                                 
  USE GC_error_module

  IMPLICIT NONE

CONTAINS

  SUBROUTINE convolve_slit(error)

    IMPLICIT NONE

    ! ------------------
    ! Modified variables
    ! ------------------
    LOGICAL,       INTENT(INOUT) :: error ! Error variable

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER :: g,v

    ! ----------------
    ! Code starts here
    ! ----------------
    error = .FALSE.

    IF (use_wavelength) THEN
       tmpwaves (1:nflambdas) = flambdas (1:nflambdas)
       tmpcwaves(1:nclambdas) = clambdas(1:nclambdas)
    ELSE
       tmpwaves (1:nflambdas) = fwavnums (1:nflambdas)
       tmpcwaves(1:nclambdas) = cwavnums(1:nclambdas)
    END IF
    
    nspec = VLIDORT_Out%Main%TS_N_GEOMETRIES
    CALL gauss_f2c (tmpwaves(1:nflambdas), GC_radiances(1:nflambdas, ilev, 1:nspec, didx), nflambdas, &
         nspec, lambda_resolution, tmpcwaves(1:nclambdas), temp_rad(1:nclambdas, 1:nspec), nclambdas)
    GC_radiances(1:nclambdas, ilev, 1:nspec, didx) = temp_rad(1:nclambdas, 1:nspec)

    nspec = VLIDORT_ModIn%MSunrays%TS_N_SZANGLES
    CALL gauss_f2c (tmpwaves(1:nflambdas), GC_flux(1:nflambdas, ilev, 1:nspec, didx), nflambdas, &
         nspec, lambda_resolution, tmpcwaves(1:nclambdas), temp_rad(1:nclambdas, 1:nspec), nclambdas)
    GC_flux(1:nclambdas, ilev, 1:nspec, didx) = temp_rad(1:nclambdas, 1:nspec)

    CALL gauss_f2c (tmpwaves(1:nflambdas), GC_direct_flux(1:nflambdas, ilev, 1:nspec, didx), nflambdas, &
         nspec, lambda_resolution, tmpcwaves(1:nclambdas), temp_rad(1:nclambdas, 1:nspec), nclambdas)
    GC_direct_flux(1:nclambdas, ilev, 1:nspec, didx) = temp_rad(1:nclambdas, 1:nspec)

    IF ( do_vector_calculation .AND. do_StokesQU_output ) THEN

       nspec = VLIDORT_Out%Main%TS_N_GEOMETRIES
       CALL gauss_f2c (tmpwaves(1:nflambdas), GC_Qvalues(1:nflambdas, ilev, 1:nspec, didx), nflambdas, &
            nspec, lambda_resolution, tmpcwaves(1:nclambdas), temp_rad(1:nclambdas, 1:nspec), nclambdas)
       GC_Qvalues(1:nclambdas, ilev, 1:nspec, didx) = temp_rad(1:nclambdas, 1:nspec)
       
       CALL gauss_f2c (tmpwaves(1:nflambdas), GC_Uvalues(1:nflambdas, ilev, 1:nspec, didx), nflambdas, &
            nspec, lambda_resolution, tmpcwaves(1:nclambdas), temp_rad(1:nclambdas, 1:nspec), nclambdas)
       GC_Uvalues(1:nclambdas, ilev, 1:nspec, didx) = temp_rad(1:nclambdas, 1:nspec)

       nspec = VLIDORT_ModIn%MSunrays%TS_N_SZANGLES
       CALL gauss_f2c (tmpwaves(1:nflambdas), GC_Qflux(1:nflambdas, ilev, 1:nspec, didx), nflambdas, &
            nspec, lambda_resolution, tmpcwaves(1:nclambdas), temp_rad(1:nclambdas, 1:nspec), nclambdas)
       GC_Qflux(1:nclambdas, ilev, 1:nspec, didx) = temp_rad(1:nclambdas, 1:nspec)

       CALL gauss_f2c (tmpwaves(1:nflambdas), GC_Uflux(1:nflambdas, ilev, 1:nspec, didx), nflambdas, &
            nspec, lambda_resolution, tmpcwaves(1:nclambdas), temp_rad(1:nclambdas, 1:nspec), nclambdas)
       GC_Uflux(1:nclambdas, ilev, 1:nspec, didx) = temp_rad(1:nclambdas, 1:nspec)
       
       CALL gauss_f2c (tmpwaves(1:nflambdas), GC_Qdirect_flux(1:nflambdas, ilev, 1:nspec, didx), nflambdas, &
            nspec, lambda_resolution, tmpcwaves(1:nclambdas), temp_rad(1:nclambdas, 1:nspec), nclambdas)
       GC_Qdirect_flux(1:nclambdas, ilev, 1:nspec, didx) = temp_rad(1:nclambdas, 1:nspec)

       CALL gauss_f2c (tmpwaves(1:nflambdas), GC_Udirect_flux(1:nflambdas, ilev, 1:nspec, didx), nflambdas, &
            nspec, lambda_resolution, tmpcwaves(1:nclambdas), temp_rad(1:nclambdas, 1:nspec), nclambdas)
       GC_Udirect_flux(1:nclambdas, ilev, 1:nspec, didx) = temp_rad(1:nclambdas, 1:nspec)

    END IF
    
    nspec = 1
    CALL gauss_f2c (tmpwaves(1:nflambdas), ground_ler(1:nflambdas), nflambdas, &
         nspec, lambda_resolution, tmpcwaves(1:nclambdas), temp_rad(1:nclambdas, 1:nspec), nclambdas)
    ground_ler(1:nclambdas) = temp_rad(1:nclambdas, 1)
    
    IF (do_Jacobians) THEN
       DO v = 1, VLIDORT_Out%Main%TS_N_GEOMETRIES
          DO g = 1, ngases
             CALL gauss_f2c (tmpwaves(1:nflambdas), GC_Tracegas_Jacobians(1:nflambdas,                    &
                  1:VLIDORT_FixIn%Cont%TS_NLAYERS, ilev, v, g, didx), nflambdas, VLIDORT_FixIn%Cont%TS_NLAYERS, &
                  lambda_resolution, tmpcwaves(1:nclambdas),                                              &
                  temp_wf(1:nclambdas, 1:VLIDORT_FixIn%Cont%TS_NLAYERS), nclambdas)
             GC_Tracegas_Jacobians(1:nclambdas, 1:VLIDORT_FixIn%Cont%TS_NLAYERS, ilev, v, g, didx) = &
                  temp_wf(1:nclambdas, 1:VLIDORT_FixIn%Cont%TS_NLAYERS)
          END DO
          
          IF (do_AMF_calculation) THEN
             CALL gauss_f2c (flambdas(1:nflambdas), GC_Scattering_Weights(1:nflambdas,                    &
                  1:VLIDORT_FixIn%Cont%TS_NLAYERS, ilev, v, didx), nflambdas, VLIDORT_FixIn%Cont%TS_NLAYERS,    &
                  lambda_resolution, clambdas(1:nclambdas),                                               &
                  temp_wf(1:nclambdas, 1:VLIDORT_FixIn%Cont%TS_NLAYERS), nclambdas)
             GC_Scattering_Weights(1:nclambdas, 1:VLIDORT_FixIn%Cont%TS_NLAYERS, ilev, v, didx) = &
                  temp_wf(1:nclambdas, 1:VLIDORT_FixIn%Cont%TS_NLAYERS)
             
             nspec = 1
             DO g = 1, ngases
                CALL gauss_f2c (flambdas(1:nflambdas), GC_AMFs(1:nflambdas, ilev, v, g, didx), nflambdas, &
                     nspec, lambda_resolution, clambdas(1:nclambdas), temp_wf(1:nclambdas, 1), nclambdas)
                GC_AMFs(1:nclambdas, ilev, v, g, didx) = temp_wf(1:nclambdas, 1)
             END DO
          END IF
          
          IF (aer_ctr%do_aer_columnwf) THEN
             nspec=1
             IF (aer_ctr%do_aod_jacobians) THEN
                CALL gauss_f2c (tmpwaves(1:nflambdas), GC_aod_Jacobians(1:nflambdas, 1, ilev, v, didx), nflambdas, &
                     nspec, lambda_resolution, tmpcwaves(1:nclambdas), temp_wf(1:nclambdas, 1), nclambdas)
                GC_aod_Jacobians(1:nclambdas, 1, ilev, v, didx) = temp_wf(1:nclambdas, 1) 
             END IF
             IF (aer_ctr%do_assa_jacobians) THEN
                CALL gauss_f2c (tmpwaves(1:nflambdas), GC_assa_Jacobians(1:nflambdas, 1, ilev, v, didx), nflambdas, &
                     nspec, lambda_resolution, tmpcwaves(1:nclambdas), temp_wf(1:nclambdas, 1), nclambdas)
                GC_assa_Jacobians(1:nclambdas, 1, ilev, v, didx) = temp_wf(1:nclambdas, 1) 
             END IF
          ELSE
             IF (aer_ctr%do_aod_jacobians) THEN
                CALL gauss_f2c (tmpwaves(1:nflambdas), GC_aod_Jacobians(1:nflambdas,      &
                     1:VLIDORT_FixIn%Cont%TS_NLAYERS, ilev, v, didx), nflambdas,                &
                     VLIDORT_FixIn%Cont%TS_NLAYERS, lambda_resolution,                    &
                     tmpcwaves(1:nclambdas),                                              &
                     temp_wf(1:nclambdas, 1:VLIDORT_FixIn%Cont%TS_NLAYERS), nclambdas)
                GC_aod_Jacobians(1:nclambdas, 1:VLIDORT_FixIn%Cont%TS_NLAYERS, ilev, v, didx) = &
                     temp_wf(1:nclambdas, 1:VLIDORT_FixIn%Cont%TS_NLAYERS) 
             END IF
             IF (aer_ctr%do_assa_jacobians) THEN
                CALL gauss_f2c (tmpwaves(1:nflambdas), GC_assa_Jacobians(1:nflambdas,      &
                     1:VLIDORT_FixIn%Cont%TS_NLAYERS, ilev, v, didx), nflambdas,                 &
                     VLIDORT_FixIn%Cont%TS_NLAYERS, lambda_resolution,                     &
                     tmpcwaves(1:nclambdas), temp_wf(1:nclambdas,                          &
                     1:VLIDORT_FixIn%Cont%TS_NLAYERS), nclambdas)
                GC_assa_Jacobians(1:nclambdas, 1:VLIDORT_FixIn%Cont%TS_NLAYERS, ilev, v, didx) = &
                     temp_wf(1:nclambdas, 1:VLIDORT_FixIn%Cont%TS_NLAYERS) 
             END IF
          ENDIF
          
          IF (do_cld_columnwf) THEN
             nspec=1
             IF (do_cod_jacobians) THEN
                CALL gauss_f2c (tmpwaves(1:nflambdas), GC_cod_Jacobians(1:nflambdas, 1, ilev, v, didx), nflambdas, &
                     nspec, lambda_resolution, tmpcwaves(1:nclambdas), temp_wf(1:nclambdas, 1), nclambdas)
                GC_cod_Jacobians(1:nclambdas, 1, ilev, v, didx) = temp_wf(1:nclambdas, 1) 
             END IF
             IF (do_cssa_jacobians) THEN
                CALL gauss_f2c (tmpwaves(1:nflambdas), GC_cssa_Jacobians(1:nflambdas, 1, ilev, v, didx), nflambdas, &
                     nspec, lambda_resolution, tmpcwaves(1:nclambdas), temp_wf(1:nclambdas, 1), nclambdas)
                GC_cssa_Jacobians(1:nclambdas, 1, ilev, v, didx) = temp_wf(1:nclambdas, 1) 
             END IF
          ELSE
             IF (do_cod_jacobians) THEN
                CALL gauss_f2c (tmpwaves(1:nflambdas), GC_cod_Jacobians(1:nflambdas,      &
                     1:VLIDORT_FixIn%Cont%TS_NLAYERS, ilev, v, didx), nflambdas,                &
                     VLIDORT_FixIn%Cont%TS_NLAYERS, lambda_resolution,                    &
                     tmpcwaves(1:nclambdas),                                              &
                     temp_wf(1:nclambdas, 1:VLIDORT_FixIn%Cont%TS_NLAYERS), nclambdas)
                GC_cod_Jacobians(1:nclambdas, 1:VLIDORT_FixIn%Cont%TS_NLAYERS, ilev, v, didx) = &
                     temp_wf(1:nclambdas, 1:VLIDORT_FixIn%Cont%TS_NLAYERS) 
             END IF
             IF (do_cssa_jacobians) THEN
                CALL gauss_f2c (tmpwaves(1:nflambdas), GC_cssa_Jacobians(1:nflambdas,      &
                     1:VLIDORT_FixIn%Cont%TS_NLAYERS, ilev, v, didx), nflambdas,                 &
                     VLIDORT_FixIn%Cont%TS_NLAYERS, lambda_resolution,                     &
                     tmpcwaves(1:nclambdas), temp_wf(1:nclambdas,                          &
                     1:VLIDORT_FixIn%Cont%TS_NLAYERS), nclambdas)
                GC_cssa_Jacobians(1:nclambdas, 1:VLIDORT_FixIn%Cont%TS_NLAYERS, ilev, v, didx) = &
                     temp_wf(1:nclambdas, 1:VLIDORT_FixIn%Cont%TS_NLAYERS) 
             END IF
          END IF
          
          IF (do_cfrac_Jacobians) THEN
             nspec = 1
             CALL gauss_f2c (tmpwaves(1:nflambdas), GC_cfrac_Jacobians(1:nflambdas, ilev, v, didx), nflambdas, &
                  nspec, lambda_resolution, tmpcwaves(1:nclambdas), temp_wf(1:nclambdas, 1), nclambdas)
             GC_cfrac_Jacobians(1:nclambdas, ilev, v, didx) = temp_wf(1:nclambdas, 1) 
          END IF
          
          IF (VLIDORT_FixIn%Bool%TS_DO_LAMBERTIAN_SURFACE) THEN
             nspec = 1
             CALL gauss_f2c (tmpwaves(1:nflambdas), GC_Surfalbedo_Jacobians(1:nflambdas, ilev, v, didx), nflambdas, &
                  nspec, lambda_resolution, tmpcwaves(1:nclambdas), temp_wf(1:nclambdas, 1), nclambdas)
             GC_Surfalbedo_Jacobians(1:nclambdas, ilev, v, didx) = temp_wf(1:nclambdas, 1)
          ELSE
             nspec = 1
             CALL gauss_f2c (tmpwaves(1:nflambdas), GC_Windspeed_Jacobians(1:nflambdas, ilev, v, didx), nflambdas, &
                  nspec, lambda_resolution, tmpcwaves(1:nclambdas), temp_wf(1:nclambdas, 1), nclambdas)
             GC_Windspeed_Jacobians(1:nclambdas, ilev, v, didx) = temp_wf(1:nclambdas, 1)               
          END IF
          nspec = 1
          CALL gauss_f2c (tmpwaves(1:nflambdas), GC_sfcprs_Jacobians(1:nflambdas, ilev, v, didx), nflambdas, &
               nspec, lambda_resolution, tmpcwaves(1:nclambdas), temp_wf(1:nclambdas, 1), nclambdas)
          GC_sfcprs_Jacobians(1:nclambdas, ilev, v, didx) = temp_wf(1:nclambdas, 1) 
          
       END DO
    END IF
    
    IF ( do_Jacobians .AND. do_QU_Jacobians) THEN
       DO v = 1, VLIDORT_Out%Main%TS_N_GEOMETRIES
          DO g = 1, ngases
             CALL gauss_f2c (tmpwaves(1:nflambdas),                                                            &
                  GC_Tracegas_QJacobians(1:nflambdas, 1:VLIDORT_FixIn%Cont%TS_NLAYERS, ilev, v, g, didx), nflambdas, &
                  VLIDORT_FixIn%Cont%TS_NLAYERS, lambda_resolution, tmpcwaves(1:nclambdas),                    &
                  temp_wf(1:nclambdas, 1:VLIDORT_FixIn%Cont%TS_NLAYERS), nclambdas)
             GC_Tracegas_QJacobians(1:nclambdas, 1:VLIDORT_FixIn%Cont%TS_NLAYERS, ilev, v, g, didx) = &
                  temp_wf(1:nclambdas, 1:VLIDORT_FixIn%Cont%TS_NLAYERS)
             
             CALL gauss_f2c (tmpwaves(1:nflambdas),                                                            &
                  GC_Tracegas_UJacobians(1:nflambdas, 1:VLIDORT_FixIn%Cont%TS_NLAYERS, ilev, v, g, didx), nflambdas, &
                  VLIDORT_FixIn%Cont%TS_NLAYERS, lambda_resolution, tmpcwaves(1:nclambdas),                    &
                  temp_wf(1:nclambdas, 1:VLIDORT_FixIn%Cont%TS_NLAYERS), nclambdas)
             GC_Tracegas_UJacobians(1:nclambdas, 1:VLIDORT_FixIn%Cont%TS_NLAYERS, ilev, v, g, didx) = &
                  temp_wf(1:nclambdas, 1:VLIDORT_FixIn%Cont%TS_NLAYERS)
          END DO
          
          IF (aer_ctr%do_aer_columnwf) THEN
             nspec=1
             IF (aer_ctr%do_aod_Jacobians) THEN 
                CALL gauss_f2c (tmpwaves(1:nflambdas), GC_aod_QJacobians(1:nflambdas, 1, ilev, v, didx), nflambdas, &
                     nspec, lambda_resolution, tmpcwaves(1:nclambdas), temp_wf(1:nclambdas, 1), nclambdas)
                GC_aod_QJacobians(1:nclambdas, 1, ilev, v, didx) = temp_wf(1:nclambdas, 1) 
                CALL gauss_f2c (tmpwaves(1:nflambdas), GC_aod_UJacobians(1:nflambdas, 1, ilev, v, didx), nflambdas, &
                     nspec, lambda_resolution, tmpcwaves(1:nclambdas), temp_wf(1:nclambdas, 1), nclambdas)
                GC_aod_UJacobians(1:nclambdas, 1, ilev, v, didx) = temp_wf(1:nclambdas, 1) 
             END IF
             IF (aer_ctr%do_assa_Jacobians) THEN 
                CALL gauss_f2c (tmpwaves(1:nflambdas), GC_assa_QJacobians(1:nflambdas, 1, ilev, v, didx), nflambdas, &
                     nspec, lambda_resolution, tmpcwaves(1:nclambdas), temp_wf(1:nclambdas, 1), nclambdas)
                GC_assa_QJacobians(1:nclambdas, 1, ilev, v, didx) = temp_wf(1:nclambdas, 1) 
                
                CALL gauss_f2c (tmpwaves(1:nflambdas), GC_assa_UJacobians(1:nflambdas, 1, ilev, v, didx), nflambdas, &
                     nspec, lambda_resolution, tmpcwaves(1:nclambdas), temp_wf(1:nclambdas, 1), nclambdas)
                GC_assa_UJacobians(1:nclambdas, 1, ilev, v, didx) = temp_wf(1:nclambdas, 1) 
             END IF
          ELSE
             IF (aer_ctr%do_aod_Jacobians) THEN 
                CALL gauss_f2c (tmpwaves(1:nflambdas),                                                    &
                     GC_aod_QJacobians(1:nflambdas, 1:VLIDORT_FixIn%Cont%TS_NLAYERS, ilev, v, didx), nflambdas, &
                     VLIDORT_FixIn%Cont%TS_NLAYERS, lambda_resolution, tmpcwaves(1:nclambdas),            &
                     temp_wf(1:nclambdas, 1:VLIDORT_FixIn%Cont%TS_NLAYERS), nclambdas)
                GC_aod_QJacobians(1:nclambdas, 1:VLIDORT_FixIn%Cont%TS_NLAYERS, ilev, v, didx) = &
                     temp_wf(1:nclambdas, 1:VLIDORT_FixIn%Cont%TS_NLAYERS) 
                CALL gauss_f2c (tmpwaves(1:nflambdas),                                                    &
                     GC_aod_UJacobians(1:nflambdas, 1:VLIDORT_FixIn%Cont%TS_NLAYERS, ilev, v, didx), nflambdas, &
                     VLIDORT_FixIn%Cont%TS_NLAYERS, lambda_resolution, tmpcwaves(1:nclambdas),            &
                     temp_wf(1:nclambdas, 1:VLIDORT_FixIn%Cont%TS_NLAYERS), nclambdas)
                GC_aod_UJacobians(1:nclambdas, 1:VLIDORT_FixIn%Cont%TS_NLAYERS, ilev, v, didx) = &
                     temp_wf(1:nclambdas, 1:VLIDORT_FixIn%Cont%TS_NLAYERS) 
             END IF
             IF (aer_ctr%do_assa_Jacobians) THEN 
                CALL gauss_f2c (tmpwaves(1:nflambdas),                                                     &
                     GC_assa_QJacobians(1:nflambdas, 1:VLIDORT_FixIn%Cont%TS_NLAYERS, ilev, v, didx), nflambdas, &
                     VLIDORT_FixIn%Cont%TS_NLAYERS, lambda_resolution, tmpcwaves(1:nclambdas),             &
                     temp_wf(1:nclambdas, 1:VLIDORT_FixIn%Cont%TS_NLAYERS), nclambdas)
                GC_assa_QJacobians(1:nclambdas, 1:VLIDORT_FixIn%Cont%TS_NLAYERS, ilev, v, didx) = &
                     temp_wf(1:nclambdas, 1:VLIDORT_FixIn%Cont%TS_NLAYERS) 
                
                CALL gauss_f2c (tmpwaves(1:nflambdas),                                                     &
                     GC_assa_UJacobians(1:nflambdas, 1:VLIDORT_FixIn%Cont%TS_NLAYERS, ilev, v, didx), nflambdas, &
                     VLIDORT_FixIn%Cont%TS_NLAYERS, lambda_resolution, tmpcwaves(1:nclambdas),             &
                     temp_wf(1:nclambdas, 1:VLIDORT_FixIn%Cont%TS_NLAYERS), nclambdas)
                GC_assa_UJacobians(1:nclambdas, 1:VLIDORT_FixIn%Cont%TS_NLAYERS, ilev, v, didx) = &
                     temp_wf(1:nclambdas, 1:VLIDORT_FixIn%Cont%TS_NLAYERS) 
             END IF
          END IF
          
          IF (do_cld_columnwf) THEN
             nspec=1
             IF (do_cod_Jacobians) THEN 
                CALL gauss_f2c (tmpwaves(1:nflambdas), GC_cod_QJacobians(1:nflambdas, 1, ilev, v, didx), nflambdas, &
                     nspec, lambda_resolution, tmpcwaves(1:nclambdas), temp_wf(1:nclambdas, 1), nclambdas)
                GC_cod_QJacobians(1:nclambdas, 1, ilev, v, didx) = temp_wf(1:nclambdas, 1) 
                CALL gauss_f2c (tmpwaves(1:nflambdas), GC_cod_UJacobians(1:nflambdas, 1, ilev, v, didx), nflambdas, &
                     nspec, lambda_resolution, tmpcwaves(1:nclambdas), temp_wf(1:nclambdas, 1), nclambdas)
                GC_cod_UJacobians(1:nclambdas, 1, ilev, v, didx) = temp_wf(1:nclambdas, 1) 
             END IF
             IF (do_cssa_Jacobians) THEN 
                CALL gauss_f2c (tmpwaves(1:nflambdas), GC_cssa_QJacobians(1:nflambdas, 1, ilev, v, didx), nflambdas, &
                     nspec, lambda_resolution, tmpcwaves(1:nclambdas), temp_wf(1:nclambdas, 1), nclambdas)
                GC_cssa_QJacobians(1:nclambdas, 1, ilev, v, didx) = temp_wf(1:nclambdas, 1) 
                CALL gauss_f2c (tmpwaves(1:nflambdas), GC_cssa_UJacobians(1:nflambdas, 1, ilev, v, didx), nflambdas, &
                     nspec, lambda_resolution, tmpcwaves(1:nclambdas), temp_wf(1:nclambdas, 1), nclambdas)
                GC_cssa_UJacobians(1:nclambdas, 1, ilev, v, didx) = temp_wf(1:nclambdas, 1) 
             END IF
          ELSE
             IF (do_cod_Jacobians) THEN 
                CALL gauss_f2c (tmpwaves(1:nflambdas),                                                    &
                     GC_cod_QJacobians(1:nflambdas, 1:VLIDORT_FixIn%Cont%TS_NLAYERS, ilev, v, didx), nflambdas, &
                     VLIDORT_FixIn%Cont%TS_NLAYERS, lambda_resolution, tmpcwaves(1:nclambdas),            &
                     temp_wf(1:nclambdas, 1:VLIDORT_FixIn%Cont%TS_NLAYERS), nclambdas)
                GC_cod_QJacobians(1:nclambdas, 1:VLIDORT_FixIn%Cont%TS_NLAYERS, ilev, v, didx) = &
                     temp_wf(1:nclambdas, 1:VLIDORT_FixIn%Cont%TS_NLAYERS) 
                CALL gauss_f2c (tmpwaves(1:nflambdas),                                                    &
                     GC_cod_UJacobians(1:nflambdas, 1:VLIDORT_FixIn%Cont%TS_NLAYERS, ilev, v, didx), nflambdas, &
                     VLIDORT_FixIn%Cont%TS_NLAYERS, lambda_resolution, tmpcwaves(1:nclambdas),            &
                     temp_wf(1:nclambdas, 1:VLIDORT_FixIn%Cont%TS_NLAYERS), nclambdas)
                GC_cod_UJacobians(1:nclambdas, 1:VLIDORT_FixIn%Cont%TS_NLAYERS, ilev, v, didx) = &
                     temp_wf(1:nclambdas, 1:VLIDORT_FixIn%Cont%TS_NLAYERS) 
             END IF
             IF (do_cssa_Jacobians) THEN 
                CALL gauss_f2c (tmpwaves(1:nflambdas),                                                     &
                     GC_cssa_QJacobians(1:nflambdas, 1:VLIDORT_FixIn%Cont%TS_NLAYERS, ilev, v, didx), nflambdas, &
                     VLIDORT_FixIn%Cont%TS_NLAYERS, lambda_resolution, tmpcwaves(1:nclambdas),             &
                     temp_wf(1:nclambdas, 1:VLIDORT_FixIn%Cont%TS_NLAYERS), nclambdas)
                GC_cssa_QJacobians(1:nclambdas, 1:VLIDORT_FixIn%Cont%TS_NLAYERS, ilev, v, didx) = &
                     temp_wf(1:nclambdas, 1:VLIDORT_FixIn%Cont%TS_NLAYERS) 
                CALL gauss_f2c (tmpwaves(1:nflambdas),                                                     &
                     GC_cssa_UJacobians(1:nflambdas, 1:VLIDORT_FixIn%Cont%TS_NLAYERS, ilev, v, didx), nflambdas, &
                     VLIDORT_FixIn%Cont%TS_NLAYERS, lambda_resolution, tmpcwaves(1:nclambdas),             &
                     temp_wf(1:nclambdas, 1:VLIDORT_FixIn%Cont%TS_NLAYERS), nclambdas)
                GC_cssa_UJacobians(1:nclambdas, 1:VLIDORT_FixIn%Cont%TS_NLAYERS, ilev, v, didx) = &
                     temp_wf(1:nclambdas, 1:VLIDORT_FixIn%Cont%TS_NLAYERS) 
             END IF
          END IF
          
          IF (do_cfrac_Jacobians) THEN
             nspec = 1
             CALL gauss_f2c (tmpwaves(1:nflambdas), GC_cfrac_QJacobians(1:nflambdas, ilev, v, didx), nflambdas, &
                  nspec, lambda_resolution, tmpcwaves(1:nclambdas), temp_wf(1:nclambdas, 1), nclambdas)
             GC_cfrac_QJacobians(1:nclambdas, ilev, v, didx) = temp_wf(1:nclambdas, 1) 
             CALL gauss_f2c (tmpwaves(1:nflambdas), GC_cfrac_UJacobians(1:nflambdas, ilev, v, didx), nflambdas, &
                  nspec, lambda_resolution, tmpcwaves(1:nclambdas), temp_wf(1:nclambdas, 1), nclambdas)
             GC_cfrac_UJacobians(1:nclambdas, ilev, v, didx) = temp_wf(1:nclambdas, 1) 
          END IF
          
          IF (VLIDORT_FixIn%Bool%TS_DO_LAMBERTIAN_SURFACE) THEN
             nspec = 1
             CALL gauss_f2c (tmpwaves(1:nflambdas), GC_Surfalbedo_QJacobians(1:nflambdas, ilev, v, didx), nflambdas, &
                  nspec, lambda_resolution, tmpcwaves(1:nclambdas), temp_wf(1:nclambdas, 1), nclambdas)
             GC_Surfalbedo_QJacobians(1:nclambdas, ilev, v, didx) = temp_wf(1:nclambdas, 1)
             
             nspec = 1
             CALL gauss_f2c (tmpwaves(1:nflambdas), GC_Surfalbedo_UJacobians(1:nflambdas, ilev, v, didx), nflambdas, &
                  nspec, lambda_resolution, tmpcwaves(1:nclambdas), temp_wf(1:nclambdas, 1), nclambdas)
             GC_Surfalbedo_UJacobians(1:nclambdas, ilev, v, didx) = temp_wf(1:nclambdas, 1)
          ELSE
             nspec = 1
             CALL gauss_f2c (tmpwaves(1:nflambdas), GC_Windspeed_QJacobians(1:nflambdas, ilev, v, didx), nflambdas, &
                  nspec, lambda_resolution, tmpcwaves(1:nclambdas), temp_wf(1:nclambdas, 1), nclambdas)
             GC_Windspeed_QJacobians(1:nclambdas, ilev, v, didx) = temp_wf(1:nclambdas, 1)
             
             nspec = 1
             CALL gauss_f2c (tmpwaves(1:nflambdas), GC_Windspeed_UJacobians(1:nflambdas, ilev, v, didx), nflambdas, &
                  nspec, lambda_resolution, tmpcwaves(1:nclambdas), temp_wf(1:nclambdas, 1), nclambdas)
             GC_Windspeed_UJacobians(1:nclambdas, ilev, v, didx) = temp_wf(1:nclambdas, 1)               
          END IF
          
          nspec = 1
          CALL gauss_f2c (tmpwaves(1:nflambdas), GC_sfcprs_QJacobians(1:nflambdas, ilev, v, didx), nflambdas, &
               nspec, lambda_resolution, tmpcwaves(1:nclambdas), temp_wf(1:nclambdas, 1), nclambdas)
          GC_sfcprs_QJacobians(1:nclambdas, ilev, v, didx) = temp_wf(1:nclambdas, 1) 
          CALL gauss_f2c (tmpwaves(1:nflambdas), GC_sfcprs_UJacobians(1:nflambdas, ilev, v, didx), nflambdas, &
               nspec, lambda_resolution, tmpcwaves(1:nclambdas), temp_wf(1:nclambdas, 1), nclambdas)
          GC_sfcprs_UJacobians(1:nclambdas, ilev, v, didx) = temp_wf(1:nclambdas, 1) 
          
       END DO
    END IF
    
    IF ( do_T_Jacobians ) THEN
       DO v = 1, VLIDORT_Out%Main%TS_N_GEOMETRIES
          CALL gauss_f2c (tmpwaves(1:nflambdas),                                                                 &
               GC_Temperature_Jacobians(1:nflambdas, 1:VLIDORT_FixIn%Cont%TS_NLAYERS, ilev, v, didx), nflambdas, &
               VLIDORT_FixIn%Cont%TS_NLAYERS, lambda_resolution, tmpcwaves(1:nclambdas),                         &
               temp_wf(1:nclambdas, 1:VLIDORT_FixIn%Cont%TS_NLAYERS), nclambdas)
          GC_Temperature_Jacobians(1:nclambdas, 1:VLIDORT_FixIn%Cont%TS_NLAYERS, ilev, v, didx) = &
               temp_wf(1:nclambdas, 1:VLIDORT_FixIn%Cont%TS_NLAYERS)
       END DO
    END IF
    
    nlambdas = nclambdas
    lambdas(1:nlambdas) = clambdas(1:nclambdas)
    wavnums(1:nlambdas) = cwavnums(1:nclambdas)
    
  END SUBROUTINE convolve_slit
  
END MODULE GC_convolution_module
