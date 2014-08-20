MODULE GC_Vlidort_module

  USE VLIDORT_PARS,        ONLY: GISSCOXMUNK_CRI_IDX, VLIDORT_SERIOUS
  USE GC_parameters_module,ONLY: max_ch_len
  USE GC_variables_module, ONLY: VLIDORT_FixIn, ndir, idix, didx, GC_do_user_altitudes,        &
                                 GC_n_user_levels, GC_user_altitudes, GC_user_levels, heights, &
                                 GC_nlayers, VLIDORT_ModIn, nstreams_choice, VLIDORT_LinModIn, &
                                 VLIDORT_LinFixIn, do_vector_calculation, ngksec, nactgkmatc,  &
                                 GC_n_sun_positions, GC_sun_positions, GC_n_view_angles,       &
                                 GC_view_angles, GC_n_azimuths, GC_azimuths, gaswfidx, twfidx, &
                                 aodwfidx, assawfidx, sfcprswfidx, codwfidx, cssawfidx,        &
                                 do_Jacobians, do_T_Jacobians, do_sfcprs_Jacobians,            &
                                 do_aod_Jacobians, do_assa_Jacobians, do_cod_Jacobians,        &
                                 do_cssa_Jacobians,do_lambertian_cld, N_TOTALPROFILE_WFS_wcld, &
                                 do_clouds, N_TOTALPROFILE_WFS_ncld, do_aerosols,              &
                                 aercld_nmoments_input, use_lambertian, wind_speed,            &
                                 VBRDF_Sup_In, VBRDF_LinSup_In, VLIDORT_Out, GC_Radiances,     &
                                 do_StokesQU_output, GC_Qvalues, GC_Uvalues, ngases,           &
                                 gas_partialcolumns, gasabs, GC_Tracegas_Jacobians,            &
                                 do_normalized_WFoutput, ratio, total_gasabs, VLIDORT_LinOut,  &
                                 do_AMF_calculation, GC_Scattering_Weights, w, GC_AMFs,        &
                                 GC_Tracegas_QJacobians, GC_Tracegas_UJacobians,               &
                                 GC_Temperature_Jacobians, mid_temperatures, do_aer_columnwf,  &
                                 total_wf, aer_flags, GC_aod_Jacobians, taertau0, taer_profile,&
                                 do_QU_Jacobians, total_Qwf, total_Uwf, GC_aod_QJacobians,     &
                                 GC_aod_UJacobians, total_aersca, total_aertau, aer_opdeps,    &
                                 aer_ssalbs, total_vaerssa, GC_assa_Jacobians,                 &
                                 GC_assa_QJacobians, GC_assa_UJacobians, do_cld_columnwf,      &
                                 cld_flags, GC_cod_Jacobians, GC_cod_QJacobians,               &
                                 GC_cod_UJacobians, tcldtau0, tcld_profile, total_cldsca,      &
                                 total_cldtau, cld_opdeps, cld_ssalbs, total_cldtau,           &
                                 GC_cssa_Jacobians, GC_cssa_QJacobians, GC_cssa_UJacobians,    &
                                 total_vcldssa, GC_Surfalbedo_Jacobians, ground_ler,           &
                                 GC_Surfalbedo_QJacobians, GC_Surfalbedo_UJacobians,           &
                                 GC_WindSpeed_Jacobians, GC_WindSpeed_QJacobians,              &
                                 GC_WindSpeed_UJacobians, deltp, pressures,                    &
                                 GC_sfcprs_Jacobians, GC_sfcprs_QJacobians,                    &
                                 GC_sfcprs_UJacobians, do_normalized_radiance, do_effcrs,      &
                                 solar_cspec_Data, depol, rayleigh_depols, beta_2, pray_20,    &
                                 pray_20, pray_22, pray_12, pray_32, pray_41, phasmoms_input,  &
                                 database_dir, maxaer, naer, aer_types, solar_spec_data,       &
                                 aer_profile, lambdas, aer_reflambda, aer_relqext, aer_phfcn,  &
                                 maxcld, ncld, cld_types, cld_profile, cld_reflambda,          &
                                 cld_relqext, cld_phfcn, water_rn, water_cn, stokes_clrcld,    &
                                 stokes_flux, stokes_direct_flux,                              &
                                 profilewf_sum, surfacewf_clrcld, ipafrac,                     &
                                 cfrac, cld_uppers, lambertian_cldalb, temp, temp_sq,          &
                                 total_molabs, which_gases, gas_xsecs_type, xsec, gas_xsecs,   &
                                 o3c1_xsecs, o3c2_xsecs, xsec_o3save, total_molsca, aircolumns,&
                                 rayleigh_xsecs, total_moltau, nscatter, scaco_input,          &
                                 total_tau, total_sca, aerscaidx, scaco_input, cldscaidx,      &
                                 omega, opdeps, ssalbs, smallnum, phasmoms_total_input,        &
                                 greekmat_idxs, phasmoms_idxs, abs_o3, dxsec_dt, l_gas, l_air, &
                                 daircolumns_dt, l_phasmoms_total_input, pvar,                 &
                                 do_debug_geocape_tool, du, vlidort_sup, vlidort_linsup,       &
                                 openerrorfileflag, do_cfrac_jacobians, GC_cfrac_Jacobians,    &
                                 GC_cfrac_QJacobians, GC_cfrac_UJacobians, GC_flux,            &
                                 GC_direct_flux, GC_Qflux, GC_Uflux, GC_Qdirect_flux,          &
                                 GC_Udirect_flux, VBRDF_Sup_Out, VBRDF_LinSup_Out,             &
                                 VLIDORT_Sup, VLIDORT_LinSup, VBRDF_Sup_In, Total_brdf, GCM,   &
                                 NSTOKESSQ, do_brdf_surface, OUTPUT_WSABSA, WSA_CALCULATED,    &
                                 BSA_CALCULATED, use_footprint_info, do_sat_viewcalc,          &
                                 GC_n_user_altitudes
  USE GC_error_module

  IMPLICIT NONE

CONTAINS
  SUBROUTINE Vlidort_GC_config(error)
    
    IMPLICIT NONE
    
    ! ------------------
    ! Modified variables
    ! ------------------
    LOGICAL,       INTENT(INOUT) :: error ! Error variable

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER :: i,j,n

    ! ----------------
    ! Code starts here
    ! ----------------
    error = .FALSE.
    
    ! ========================
    ! SECTION 1 (Fixed inputs)
    ! ========================
    
    VLIDORT_FixIn%Bool%TS_DO_QUAD_OUTPUT  = .FALSE. ! No option in Vlidort Control file

    ! =====================================
    ! To clean up the GC control input file
    ! and avoid duplications between it and
    ! Vlidort control file
    ! =====================================
       do_vector_calculation = VLIDORT_FixIn%Bool%TS_DO_FULLRAD_MODE
    
    ! ----------------------------------------
    ! Check for aerosols and scattering clouds
    ! ----------------------------------------
    IF ( do_aerosols .OR. (do_clouds .AND. .NOT. do_lambertian_cld)) THEN
       IF ( aercld_nmoments_input .LT. 2 * VLIDORT_FixIn%Cont%TS_NSTREAMS) THEN
          CALL write_err_message ( .FALSE., "Not enough aerosol/cloud moments"// &
               ", must be at least 2*nstreams")
          error = .TRUE.
       END IF
    END IF
    
    ! Basic control integers
    !  -- Number of layers determined from GC profiles
    !  -- Number of fine-layer subdivisions in single-scattering = 2
    !  -- Number of expansion coefficients = 2     for Rayleigh-only
    !  -- Number of expansion coefficients = Input for Rayleigh+Aerosol
    !     NSTOKES  = Will be set in Section 2
    !     NLAYERS  = Will be set in Section 2
    
    ! Basic  numbers
    !    -- Fourier convergence accuracy (not required)
    !    -- earth radius (fixed here, but could be function of LAT & LONG)
    !    -- Refractive geometry parameter (not required)
    !    -- Geometry specification height = Botton of height grid (set later)
    VLIDORT_ModIn%MChapman%TS_EARTH_RADIUS       = 6378.d0
    VLIDORT_FixIn%Optical%TS_THERMAL_BB_INPUT(:) = 0.0d0
    VLIDORT_FixIn%Optical%TS_SURFACE_BB_INPUT    = 0.0d0
    
    
    ! ========================================================
    ! SECTION 2 (Other VLIDORT inputs depending on GC control)
    ! ========================================================
    IF (VLIDORT_FixIn%Bool%TS_DO_UPWELLING .AND. VLIDORT_FixIn%Bool%TS_DO_DNWELLING) THEN
       ndir = 2; idix(1) = 1; idix(2) = 2
    ELSEIF (VLIDORT_FixIn%Bool%TS_DO_UPWELLING .AND. .NOT. VLIDORT_FixIn%Bool%TS_DO_DNWELLING) THEN
       ndir = 1; idix(1) = 1
    ELSEIF (.NOT. VLIDORT_FixIn%Bool%TS_DO_UPWELLING .AND. VLIDORT_FixIn%Bool%TS_DO_DNWELLING) THEN
       ndir = 1; idix(1) = 2
    END IF
    
    ! --------------------------------------------------------------------------
    ! Output level: Usually it comes from VLIDORT input file. In GC control file
    ! it is possible to use altitude in km to set the output levels.
    ! --------------------------------------------------------------------------
    IF (GC_do_user_altitudes) THEN
       DO i = 1, GC_n_user_altitudes
          IF (GC_user_altitudes(i) >= heights(0)) THEN
             GC_user_levels(i) = 0.0
          ELSE IF (GC_user_altitudes(i) <= heights(GC_nlayers)) THEN
             GC_user_levels(i) = GC_nlayers
             GC_user_altitudes(i) = heights(GC_nlayers)
          ELSE
             DO j = 1, GC_nlayers
                IF (GC_user_altitudes(i) >= heights(j)) THEN
                   GC_user_levels(i) = (heights(j-1) - GC_user_altitudes(i)) &
                        / (heights(j-1) - heights(j)) + j - 1
                   EXIT
                ENDIF
             ENDDO
          ENDIF
       ENDDO
       VLIDORT_FixIn%UserVal%TS_N_USER_LEVELS = GC_n_user_altitudes
       GC_n_user_levels = GC_n_user_altitudes
       DO n = 1, VLIDORT_FixIn%UserVal%TS_N_USER_LEVELS
          VLIDORT_ModIn%MUserVal%TS_USER_LEVELS(n) = GC_user_levels(n)
       END DO
    ELSE
       GC_n_user_levels = VLIDORT_FixIn%UserVal%TS_N_USER_LEVELS
       GC_user_levels(1:GC_n_user_levels) = VLIDORT_ModIn%MUserVal%TS_USER_LEVELS(1:GC_n_user_levels)
       DO i = 1, GC_n_user_levels
          j = INT(GC_user_levels(i))
          GC_user_altitudes(i) = heights(j) - (GC_user_levels(i) - j) &
               * (heights(j+1) - heights(j))
       ENDDO
    ENDIF

    ! Set the Number of Stokes parameters (1 or 3) and layers
    VLIDORT_FixIn%Cont%TS_NSTOKES = 1
    IF ( do_vector_calculation ) VLIDORT_FixIn%Cont%TS_NSTOKES = 3 
    
    ngksec = 1
    IF ( do_vector_calculation ) ngksec = 6 
    
    nactgkmatc = 1
    IF (do_vector_calculation) nactgkmatc = 8
    
    VLIDORT_FixIn%Cont%TS_NLAYERS = GC_nlayers
    
    ! Set the height grid, and Geometry specification height
    VLIDORT_FixIn%Chapman%TS_height_grid(0:VLIDORT_FixIn%Cont%TS_NLAYERS) = &
         heights(0:VLIDORT_FixIn%Cont%TS_NLAYERS)
    VLIDORT_ModIn%MUserVal%TS_GEOMETRY_SPECHEIGHT = heights(VLIDORT_FixIn%Cont%TS_NLAYERS)

    ! Now will be done in the clou/clear loop  
    ! Reduce # of layers for Lambertian cloud surface (now done in clear/cloud loop)
    ! if (do_clouds .and. do_lambertian_cld) NLAYERS = cld_uppers(1) - 1
    
    ! Set the Geometry. Just a straight copy of VLIDORT input into GC_ inputs if we are
    ! not using footprint information. If doing so, the angles have been set before in
    ! GC_profiles_module.f90
    IF (use_footprint_info .OR. do_sat_viewcalc) THEN
       VLIDORT_ModIn%MSunrays%TS_N_SZANGLES  = GC_n_sun_positions
       VLIDORT_ModIn%MSunrays%TS_SZANGLES(1) = GC_sun_positions(1)
       
       VLIDORT_ModIn%MUserVal%TS_N_USER_VZANGLES        = GC_n_view_angles
       VLIDORT_ModIn%MUserVal%TS_USER_VZANGLES_INPUT(1) = GC_view_angles(1) 
       
       VLIDORT_ModIn%MUserVal%TS_N_USER_RELAZMS  = GC_n_azimuths
       VLIDORT_ModIn%MUserVal%TS_USER_RELAZMS(1) = GC_azimuths(1)
    ELSE
       GC_n_sun_positions = VLIDORT_ModIn%MSunrays%TS_N_SZANGLES
       DO n = 1, VLIDORT_ModIn%MSunrays%TS_N_SZANGLES
          GC_sun_positions(n) = VLIDORT_ModIn%MSunrays%TS_SZANGLES(n)
       END DO       
       GC_n_view_angles = VLIDORT_ModIn%MUserVal%TS_N_USER_VZANGLES
       DO n = 1, VLIDORT_ModIn%MUserVal%TS_N_USER_VZANGLES
          GC_view_angles(n) = VLIDORT_ModIn%MUserVal%TS_USER_VZANGLES_INPUT(n)
       END DO
       GC_n_azimuths = VLIDORT_ModIn%MUserVal%TS_N_USER_RELAZMS
       DO n = 1, VLIDORT_ModIn%MUserVal%TS_N_USER_RELAZMS
          GC_azimuths(n) = VLIDORT_ModIn%MUserVal%TS_USER_RELAZMS(n)
       END DO
    ENDIF

    ! -----------------------------------------------------------------------
    ! Set the linearization control for profile Jacobians (cloudy conditions)
    ! Might be different for clear-sky part
    ! -----------------------------------------------------------------------
    ! Linearization inputs
    !    -- No total column Jacobians, no BRDF Jacobians
    !    -- profile and Surface Jacobians, flags set by GC input control
    VLIDORT_LinModIn%MCont%TS_DO_COLUMN_LINEARIZATION = .FALSE.
    VLIDORT_LinFixIn%Cont%TS_N_TOTALCOLUMN_WFS        = 0

    ! Set up do_JACOBIANS variable based in Vlidort input file
    IF ( (VLIDORT_LinFixIn%Cont%TS_do_simulation_only        .EQV. .FALSE.) .AND. &
         (VLIDORT_LinModIn%MCont%TS_do_profile_linearization .EQV. .TRUE. ) ) THEN
       do_JACOBIANS = .TRUE.
       VLIDORT_LinModIn%MCont%TS_do_atmos_linearization = .TRUE.
       VLIDORT_LinModIn%MCont%TS_do_linearization       = .TRUE.
    ELSE
       do_JACOBIANS = .FALSE.
       VLIDORT_LinModIn%MCont%TS_do_atmos_linearization = .FALSE.
       VLIDORT_LinModIn%MCont%TS_do_linearization       = .FALSE.
    ENDIF

    ! Only able to do Stokes output if do_vector_calculation
    IF ( (do_vector_calculation .EQV. .FALSE.) .AND. &
         (do_StokesQU_Output    .EQV. .TRUE. ) ) THEN
       WRITE(*,*) 'Vector calculation set to false and Stokes output to true'
       WRITE(*,*) 'Check input setup in GC_input_file'
       STOP
    ENDIF

    ! If do vector calculation and Stokes output set do_QU_Jacobians to .TRUE.
    IF ( (do_vector_calculation .EQV. .TRUE.) .AND. &
         (do_StokesQU_Output    .EQV. .TRUE. ) ) THEN
       do_QU_Jacobians = .TRUE.
    ELSE
       do_QU_Jacobians = .FALSE.
    ENDIF

    ! If not Jacobians and then it is not possible to do AMF
    IF ( (do_JACOBIANS       .EQV. .FALSE.) .AND. &
          do_AMF_calculation .EQV. .TRUE. ) THEN
       WRITE(*,*) 'You can not do AMFs without doing Jacobians first'
       WRITE(*,*) 'Check Vlidort input and GC_Input'
       STOP
    ENDIF

    IF ( ( do_T_Jacobians   .OR. do_sfcprs_Jacobians .OR.  &
           do_aod_Jacobians .OR. do_assa_Jacobians   .OR.  &
           do_cod_Jacobians .OR. do_cssa_Jacobians)  .AND. &
           (.NOT. do_Jacobians) ) THEN
       WRITE(*,*) do_Jacobians
       WRITE(*,*) 'You have requested some jacobians while the main linearization'
       WRITE(*,*) 'flags in VLIDORT control file are not set for that:'
       WRITE(*,*) 'Do simulation only? F' 
       WRITE(*,*) 'Do atmospheric profile weighting functions? T'
       STOP
    ENDIF

    gaswfidx = 0; twfidx = 0; aodwfidx = 0; assawfidx = 0
    sfcprswfidx = 0; codwfidx = 0; cssawfidx = 0
    IF ( do_JACOBIANS ) THEN

       VLIDORT_LinFixIn%Cont%TS_N_TOTALPROFILE_WFS        = 1
       VLIDORT_LinFixIn%Cont%TS_profilewf_names(1)       = '-Trace Gas Volume Mixing Ratio-'
       gaswfidx = 1
       
       IF ( do_T_Jacobians ) THEN
          VLIDORT_LinFixIn%Cont%TS_N_TOTALPROFILE_WFS    = 2
          VLIDORT_LinFixIn%Cont%TS_profilewf_names(2)    = '------Layer Temperatures-------'
          twfidx = 2
       END IF
       
       IF (do_sfcprs_Jacobians) THEN
          VLIDORT_LinFixIn%Cont%TS_N_TOTALPROFILE_WFS    = VLIDORT_LinFixIn%Cont%TS_N_TOTALPROFILE_WFS + 1
          sfcprswfidx = VLIDORT_LinFixIn%Cont%TS_N_TOTALPROFILE_WFS
          VLIDORT_LinFixIn%Cont%TS_profilewf_names(VLIDORT_LinFixIn%Cont%TS_N_TOTALPROFILE_WFS)  = &
               '-------Surface Pressure--------'
       END IF
       
       IF (do_aod_Jacobians) THEN
          VLIDORT_LinFixIn%Cont%TS_N_TOTALPROFILE_WFS = VLIDORT_LinFixIn%Cont%TS_N_TOTALPROFILE_WFS + 1
          aodwfidx = VLIDORT_LinFixIn%Cont%TS_N_TOTALPROFILE_WFS
          VLIDORT_LinFixIn%Cont%TS_profilewf_names(VLIDORT_LinFixIn%Cont%TS_N_TOTALPROFILE_WFS)  = &
               '-----Aerosol Optical Depth-----'
       END IF
       
       IF (do_assa_Jacobians) THEN
          VLIDORT_LinFixIn%Cont%TS_N_TOTALPROFILE_WFS = VLIDORT_LinFixIn%Cont%TS_N_TOTALPROFILE_WFS + 1
          assawfidx = VLIDORT_LinFixIn%Cont%TS_N_TOTALPROFILE_WFS
          VLIDORT_LinFixIn%Cont%TS_profilewf_names(VLIDORT_LinFixIn%Cont%TS_N_TOTALPROFILE_WFS)  = &
               '---Aerosol Single Sca Albedo---'
       END IF
       
       N_TOTALPROFILE_WFS_ncld = VLIDORT_LinFixIn%Cont%TS_N_TOTALPROFILE_WFS
       
       ! Put cloud related Jacobians at the end
       IF (do_cod_Jacobians) THEN
          VLIDORT_LinFixIn%Cont%TS_N_TOTALPROFILE_WFS = VLIDORT_LinFixIn%Cont%TS_N_TOTALPROFILE_WFS + 1
          codwfidx = VLIDORT_LinFixIn%Cont%TS_N_TOTALPROFILE_WFS
          VLIDORT_LinFixIn%Cont%TS_profilewf_names(VLIDORT_LinFixIn%Cont%TS_N_TOTALPROFILE_WFS)  = &
               '------Cloud Optical Depth------'
       END IF
       
       IF (do_cssa_Jacobians) THEN
          VLIDORT_LinFixIn%Cont%TS_N_TOTALPROFILE_WFS = VLIDORT_LinFixIn%Cont%TS_N_TOTALPROFILE_WFS + 1
          cssawfidx = VLIDORT_LinFixIn%Cont%TS_N_TOTALPROFILE_WFS
          VLIDORT_LinFixIn%Cont%TS_profilewf_names(VLIDORT_LinFixIn%Cont%TS_N_TOTALPROFILE_WFS)  = &
               '----Cloud Single Sca Albedo----'
       END IF
       
       N_TOTALPROFILE_WFS_wcld = VLIDORT_LinFixIn%Cont%TS_N_TOTALPROFILE_WFS
       
       VLIDORT_LinFixIn%Cont%TS_layer_vary_flag(1:VLIDORT_FixIn%Cont%TS_NLAYERS)   = .TRUE.
       VLIDORT_LinFixIn%Cont%TS_layer_vary_number(1:VLIDORT_FixIn%Cont%TS_NLAYERS) = &
            VLIDORT_LinFixIn%Cont%TS_N_TOTALPROFILE_WFS
    ELSE
       VLIDORT_LinModIn%MCont%TS_do_profile_linearization = .FALSE.
       VLIDORT_LinModIn%MCont%TS_do_atmos_linearization   = .FALSE.
       VLIDORT_LinModIn%MCont%TS_do_linearization         = .FALSE.
       VLIDORT_LinFixIn%Cont%TS_N_TOTALPROFILE_WFS        = 0
       n_totalprofile_wfs_ncld  = 0
       n_totalprofile_wfs_wcld  = 0
    END IF
    
    ! Rayleigh only flag, set if no aerosols and scattering clouds
    VLIDORT_ModIn%MBool%TS_DO_RAYLEIGH_ONLY = .NOT. do_aerosols &
         .AND. .NOT. (do_clouds .AND. .NOT. do_lambertian_cld)
    
    ! Rayleigh only, no delta-M scaling, no performance enhancements, nmoms=2
    IF ( VLIDORT_ModIn%MBool%TS_DO_RAYLEIGH_ONLY ) THEN
       VLIDORT_ModIn%MBool%TS_DO_DOUBLE_CONVTEST   = .FALSE.
       VLIDORT_ModIn%MBool%TS_DO_DELTAM_SCALING    = .FALSE.
       VLIDORT_ModIn%MBool%TS_DO_SOLUTION_SAVING   = .FALSE.
       VLIDORT_ModIn%MBool%TS_DO_BVP_TELESCOPING   = .FALSE.
       VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT = 2
    ENDIF
    
    ! Rayleigh+Aerosols, delta-M scaling, performance enhancements, using input moments
    IF ( .not. VLIDORT_ModIn%MBool%TS_DO_RAYLEIGH_ONLY ) THEN
       VLIDORT_ModIn%MBool%TS_DO_DOUBLE_CONVTEST   = .TRUE. ! Change after Xiong settings
       VLIDORT_ModIn%MBool%TS_DO_DELTAM_SCALING    = .TRUE. ! Change after Xiong settings
       VLIDORT_FixIn%Bool%TS_DO_SSCORR_TRUNCATION  = .FALSE.
       VLIDORT_ModIn%MBool%TS_DO_SOLUTION_SAVING   = .FALSE.
       VLIDORT_ModIn%MBool%TS_DO_BVP_TELESCOPING   = .FALSE.
       VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT = aercld_nmoments_input
    ENDIF
    
    ! --------------------------
    ! BRDF or Lambertian surface
    ! --------------------------
    !    -- Kernel settings are automatic, but we'll put them here
    !    -- Actual Albedo comes from the data extraction, see later on.
    ! -----------------------------------------------------------------
    VLIDORT_FixIn%Bool%TS_DO_LAMBERTIAN_SURFACE = use_lambertian 
    !   LAMBERTIAN_ALBEDO     = Will be set in Wavelength loop
    ! Security check, if TS_DO_LAMBERTIAN and BS_DO_BRDF_SURFACE are
    ! not compatible, stop program and print message to screen
    IF ( VLIDORT_FixIn%Bool%TS_DO_LAMBERTIAN_SURFACE .EQV. &
         VBRDF_Sup_In%BS_DO_BRDF_SURFACE ) THEN
       WRITE(*,*) 'Or you have a Lambertian Surface or you have '//&
            'BRDF surface but not both. Check input files.'
       STOP
    ENDIF
    
    ! Set the linearization control for surface (albedo) linearization
    IF (VLIDORT_FixIn%Bool%TS_DO_LAMBERTIAN_SURFACE) THEN
       IF ( do_JACOBIANS ) THEN
          VLIDORT_LinModIn%MCont%TS_DO_SURFACE_LINEARIZATION = .TRUE.
          VLIDORT_LinFixIn%Cont%TS_N_SURFACE_WFS = 1
       ELSE
          VLIDORT_LinModIn%MCont%TS_DO_SURFACE_LINEARIZATION = .FALSE.
       END IF
    ELSE
       IF ( DO_JACOBIANS ) THEN
          VLIDORT_LinModIn%MCont%TS_DO_SURFACE_LINEARIZATION = .TRUE.
          VLIDORT_LinFixIn%Cont%TS_N_SURFACE_WFS = 1
       ELSE
          VLIDORT_LinModIn%MCont%TS_DO_SURFACE_LINEARIZATION = .FALSE.
          VLIDORT_LinFixIn%Cont%TS_N_SURFACE_WFS             = 0
       END IF
    END IF

    ! Simulation only settings
    IF ( .NOT. do_JACOBIANS ) THEN
       VLIDORT_LinFixIn%Cont%TS_N_SURFACE_WFS             = 0
       VLIDORT_LinFixIn%Cont%TS_N_TOTALPROFILE_WFS        = 0
       VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_FLAG           = .FALSE.
       VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_NUMBER         = 0
    END IF
    
  END SUBROUTINE Vlidort_GC_config

  SUBROUTINE save_results(error)

    IMPLICIT NONE

    ! ------------------
    ! Modified variables
    ! ------------------
    LOGICAL,       INTENT(INOUT) :: error ! Error variable

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER :: IB, UM, UA, V, q, g, n, il
    
    ! ----------------
    ! Code starts here
    ! ----------------
    error = .FALSE.

    ! ============
    ! save results
    ! ============

    ! ---------
    ! Save BRDF
    ! ---------
    IF (do_brdf_surface) THEN
       GCM = (VBRDF_Sup_In%BS_which_brdf(1) .eq. 10 .or. &
            VBRDF_Sup_In%BS_which_brdf(2) .eq. 10 .or. &
            VBRDF_Sup_In%BS_which_brdf(3) .eq. 10)
       
       if ( GCM ) then
          NSTOKESSQ = VBRDF_Sup_In%BS_NSTOKES * VBRDF_Sup_In%BS_NSTOKES
       else
          NSTOKESSQ = 1
       endif
       Total_brdf = VLIDORT_Sup%BRDF%TS_EXACTDB_BRDFUNC
       OUTPUT_WSABSA = .FALSE.
       IF (VBRDF_Sup_In%BS_DO_WSABSA_OUTPUT) THEN
         WSA_CALCULATED = VBRDF_Sup_Out%BS_WSA_CALCULATED
         BSA_CALCULATED = VBRDF_Sup_Out%BS_BSA_CALCULATED
         OUTPUT_WSABSA = .TRUE.
      END IF
    ENDIF

    ! ------------------------------------
    ! Save Radiances, flux and direct flux
    ! ------------------------------------
    DO IB = 1, VLIDORT_ModIn%MSunrays%TS_N_SZANGLES
       GC_flux(W,IB)        = VLIDORT_Out%Main%TS_FLUX_STOKES(1,IB,1,didx)
       GC_direct_flux(W,IB) = VLIDORT_Out%Main%TS_FLUX_DIRECT(1,IB,1)
       DO UM = 1, VLIDORT_ModIn%MUserVal%TS_N_USER_VZANGLES
          DO UA = 1, VLIDORT_ModIn%MUserVal%TS_N_USER_RELAZMS
             V = VLIDORT_Out%Main%TS_VZA_OFFSETS(IB,UM) + UA
             GC_Radiances(W,V)   = VLIDORT_Out%Main%TS_STOKES(1,V,1,didx)
          END DO
       END DO
    END DO
    
    ! -------------------------------------------------
    ! Save Q and U values (optional output, if flagged)
    ! -------------------------------------------------
    IF ( do_vector_calculation .AND. do_StokesQU_output ) THEN
       DO IB = 1, VLIDORT_ModIn%MSunrays%TS_N_SZANGLES
                GC_Qflux(W,IB)        = VLIDORT_Out%Main%TS_FLUX_STOKES(1,IB,2,didx)
                GC_Uflux(W,IB)        = VLIDORT_Out%Main%TS_FLUX_STOKES(1,IB,3,didx)
                GC_Qdirect_flux(W,IB) = VLIDORT_Out%Main%TS_FLUX_DIRECT(1,IB,2)
                GC_Udirect_flux(W,IB) = VLIDORT_Out%Main%TS_FLUX_DIRECT(1,IB,3)
          DO UM = 1, VLIDORT_ModIn%MUserVal%TS_N_USER_VZANGLES   
             DO UA = 1, VLIDORT_ModIn%MUserVal%TS_N_USER_RELAZMS
                V = VLIDORT_Out%Main%TS_VZA_OFFSETS(IB,UM) + UA
                GC_Qvalues(W,V)      = VLIDORT_Out%Main%TS_STOKES(1,V,2,didx)
                GC_Uvalues(W,V)      = VLIDORT_Out%Main%TS_STOKES(1,V,3,didx)
             END DO
          END DO
       END DO
    END IF

    ! ------------------------
    ! Save trace gas Jacobians
    ! -----------------------------------------------------
    !    --- VLIDORT output is normalized = GAS.dRad/dGAS
    !    --- Other trace gases, use ratio of cross-sections
    ! -----------------------------------------------------
    IF ( do_Jacobians ) THEN
       q = gaswfidx
       DO g = 1, ngases
          DO n = 1, VLIDORT_FixIn%Cont%TS_NLAYERS
             IF ( gas_partialcolumns(n,g) .EQ. 0.0d0 .OR. gasabs(n, g) .EQ. 0.d0 ) THEN
                GC_Tracegas_Jacobians(w,n,1:VLIDORT_Out%Main%TS_N_GEOMETRIES,g) = 0.0d0
             ELSE
                IF ( do_normalized_WFoutput ) THEN                        !  Normalized output
                   ratio = gasabs(n, g) / total_gasabs(n)
                ELSE                                                     !  Non-normalized output
                   ratio = gasabs(n, g) / gas_partialcolumns(n,g) / total_gasabs(n)                
                END IF
                DO v = 1, VLIDORT_Out%Main%TS_N_GEOMETRIES
                   GC_Tracegas_Jacobians(w,n,v,g) = &
                        VLIDORT_LinOut%Prof%TS_PROFILEWF(q,n,1,v,1,didx) * ratio
                END DO
             END IF
          END DO            
       END DO
       
       ! Compute scattering weights
       IF (do_AMF_calculation) THEN
          DO v = 1, VLIDORT_Out%Main%TS_N_GEOMETRIES
             DO n = 1, VLIDORT_FixIn%Cont%TS_NLAYERS
                GC_Scattering_Weights(w, n, v) = &
                     -VLIDORT_LinOut%Prof%TS_PROFILEWF(q,n,1,v,1,didx) &
                     / total_gasabs(n) / GC_Radiances(w,v)
             END DO
             
             ! Compute Air mass factor (don't need to save for every one)
             DO g = 1, ngases
                GC_AMFs(w, v, g) = &
                     SUM(GC_Scattering_weights(w, 1:VLIDORT_FixIn%Cont%TS_NLAYERS, v) &
                     * gas_partialcolumns(1:VLIDORT_FixIn%Cont%TS_NLAYERS,g)) /       &
                     SUM(gas_partialcolumns(1:VLIDORT_FixIn%Cont%TS_NLAYERS,g))
             END DO
          END DO
       END IF
    END IF
    
    ! -----------------------------------------
    ! Save trace gas Jacobians, Q and U values.
    ! ...OPTIONAL OUTPUT ONLY FOR TRACE GASES..
    ! -----------------------------------------
    !    --- VLIDORT output is normalized = GAS.dRad/dGAS
    !    --- Other trace gases, use ratio of cross-sections
    
    IF ( do_QU_Jacobians) THEN
       q = gaswfidx
       DO g = 1, ngases
          DO n = 1, VLIDORT_FixIn%Cont%TS_NLAYERS
             IF ( gas_partialcolumns(n,g) .EQ. 0.0d0 .OR. gasabs(n, g) == 0.d0) THEN
                GC_Tracegas_QJacobians(w,n,1:VLIDORT_Out%Main%TS_N_GEOMETRIES,g) = 0.0d0
                GC_Tracegas_UJacobians(w,n,1:VLIDORT_Out%Main%TS_N_GEOMETRIES,g) = 0.0d0
             ELSE
                IF ( do_normalized_WFoutput ) THEN                        !  Normalized output
                   ratio = gasabs(n,g) / total_gasabs(n)
                ELSE                                                      !  Non-normalized output
                   ratio = gasabs(n,g) / gas_partialcolumns(n, g) / total_gasabs(n)
                END IF
                DO v = 1, VLIDORT_Out%Main%TS_N_GEOMETRIES
                   GC_Tracegas_QJacobians(w,n,v,g) = ratio * &
                        VLIDORT_LinOut%Prof%TS_PROFILEWF(q,n,1,v,2,didx)
                   GC_Tracegas_UJacobians(w,n,v,g) = ratio * &
                        VLIDORT_LinOut%Prof%TS_PROFILEWF(q,n,1,v,3,didx)
                END DO
             END IF
          END DO            
       END DO
    END IF
    
    ! -------------------------
    ! Save temperature Jacobian
    ! -------------------------
    IF ( do_T_Jacobians ) THEN
       q = twfidx
       DO n = 1, VLIDORT_FixIn%Cont%TS_NLAYERS
          DO v = 1, VLIDORT_Out%Main%TS_N_GEOMETRIES
             IF ( do_normalized_WFoutput ) THEN    !  Normalized output
                GC_Temperature_Jacobians(w,n,v) = &
                     VLIDORT_LinOut%Prof%TS_PROFILEWF(q,n,1,v,1,didx)
             ELSE                                  !  Non-normalized output
                GC_Temperature_Jacobians(w,n,v) = &
                     VLIDORT_LinOut%Prof%TS_PROFILEWF(q,n,1,v,1,didx) / mid_temperatures(n)
             END IF
          END DO
       END DO
    END IF
    
    ! -----------------------------------------
    ! Save aerosols optical thickness Jacobians
    ! -----------------------------------------
    IF ( do_aod_Jacobians ) THEN
       q = aodwfidx
       DO v = 1, VLIDORT_Out%Main%TS_N_GEOMETRIES
          IF (do_aer_columnwf) THEN
             total_wf = 0.0d0
             DO il = 1, VLIDORT_FixIn%Cont%TS_NLAYERS
                IF (aer_flags(il)) THEN
                   total_wf = total_wf + VLIDORT_LinOut%Prof%TS_PROFILEWF(q,il,1,v,1,didx)
                END IF
             END DO
             GC_aod_Jacobians(w, 1, v) = total_wf               ! Normalized output
             IF ( .NOT. do_normalized_WFoutput ) THEN           ! Non-normalized output
                GC_aod_Jacobians(w, 1, v) = total_wf / taertau0 ! Using tau at ref. lambda
             END IF
          ELSE
             GC_aod_Jacobians(w, 1:VLIDORT_FixIn%Cont%TS_NLAYERS, v) = &
                  VLIDORT_LinOut%Prof%TS_PROFILEWF(q,1:VLIDORT_FixIn%Cont%TS_NLAYERS,1,v,1,didx)
             IF ( .NOT. do_normalized_WFoutput ) THEN     
                WHERE( aer_flags(1:VLIDORT_FixIn%Cont%TS_NLAYERS) )
                   GC_aod_Jacobians(w, 1:VLIDORT_FixIn%Cont%TS_NLAYERS, v) =    &
                        GC_aod_Jacobians(w, 1:VLIDORT_FixIn%Cont%TS_NLAYERS, v) &
                        / taer_profile(1:VLIDORT_FixIn%Cont%TS_NLAYERS)
                END WHERE
             END IF
          END IF
       END DO
    END IF
    
    IF (  do_aod_Jacobians .AND. do_QU_Jacobians) THEN
       q = aodwfidx
       DO v = 1, VLIDORT_Out%Main%TS_N_GEOMETRIES
          IF (do_aer_columnwf) THEN
             total_Qwf = 0.0d0; total_Uwf = 0.0d0
             DO il = 1, VLIDORT_FixIn%Cont%TS_NLAYERS
                IF (aer_flags(il)) THEN
                   total_Qwf = total_Qwf + VLIDORT_LinOut%Prof%TS_PROFILEWF(q,il,1,v,2,didx)
                   total_Uwf = total_Uwf + VLIDORT_LinOut%Prof%TS_PROFILEWF(q,il,1,v,3,didx)
                END IF
             END DO
             GC_aod_QJacobians(w, 1, v) = total_Qwf              ! Normalized output
             GC_aod_UJacobians(w, 1, v) = total_Uwf              !
             IF ( .NOT. do_normalized_WFoutput ) THEN            ! Non-normalized output
                GC_aod_QJacobians(w, 1, v) = total_Qwf / taertau0! Using tau at ref. lambda
                GC_aod_UJacobians(w, 1, v) = total_Uwf / taertau0
             END IF
          ELSE
             GC_aod_QJacobians(w, 1:VLIDORT_FixIn%Cont%TS_NLAYERS, v) = &
                  VLIDORT_LinOut%Prof%TS_PROFILEWF(q,1:VLIDORT_FixIn%Cont%TS_NLAYERS,1,v,2,didx)
             GC_aod_UJacobians(w, 1:VLIDORT_FixIn%Cont%TS_NLAYERS, v) = &
                  VLIDORT_LinOut%Prof%TS_PROFILEWF(q,1:VLIDORT_FixIn%Cont%TS_NLAYERS,1,v,3,didx)
             IF ( .NOT. do_normalized_WFoutput ) THEN     
                WHERE( aer_flags(1:VLIDORT_FixIn%Cont%TS_NLAYERS) )
                   GC_aod_QJacobians(w, 1:VLIDORT_FixIn%Cont%TS_NLAYERS, v) =    &
                        GC_aod_QJacobians(w, 1:VLIDORT_FixIn%Cont%TS_NLAYERS, v) &
                        / taer_profile(1:VLIDORT_FixIn%Cont%TS_NLAYERS)
                   GC_aod_UJacobians(w, 1:VLIDORT_FixIn%Cont%TS_NLAYERS, v) =    &
                        GC_aod_UJacobians(w, 1:VLIDORT_FixIn%Cont%TS_NLAYERS, v) &
                        / taer_profile(1:VLIDORT_FixIn%Cont%TS_NLAYERS)
                END WHERE
             END IF
          END IF
       END DO
    END IF
    
    ! ------------------------------------------------
    ! Save aerosols single scattering albedo Jacobians
    ! ------------------------------------------------
    IF ( do_assa_Jacobians ) THEN
       q = assawfidx
       DO v = 1, VLIDORT_Out%Main%TS_N_GEOMETRIES
          IF (do_aer_columnwf) THEN
             total_wf = 0.0d0; total_aersca = 0.0d0; total_aertau = 0.0d0
             DO il = 1, VLIDORT_FixIn%Cont%TS_NLAYERS
                IF (aer_flags(il)) THEN
                   total_wf = total_wf + VLIDORT_LINOUT%PROF%TS_PROFILEWF(q,il,1,v,1,didx) 
                   total_aersca = total_aersca + aer_opdeps(w, il) * aer_ssalbs(w, il) 
                   total_aertau = total_aertau + aer_opdeps(w, il) 
                END IF
             END DO
             total_vaerssa = total_aersca / total_aertau
             GC_assa_Jacobians(w, 1, v) = total_wf       ! Normalized output
             IF ( .NOT. do_normalized_WFoutput ) THEN    ! Non-normalized output
                GC_assa_Jacobians(w, 1, v) = total_wf / total_vaerssa
             END IF
          ELSE
             GC_assa_Jacobians(w, 1:VLIDORT_FixIn%Cont%TS_NLAYERS, v) = &
                  VLIDORT_LINOUT%PROF%TS_PROFILEWF(q,1:VLIDORT_FixIn%Cont%TS_NLAYERS,1,v,1,didx)
             IF ( .NOT. do_normalized_WFoutput ) THEN     
                WHERE( aer_flags(1:VLIDORT_FixIn%Cont%TS_NLAYERS) )
                   GC_assa_Jacobians(w, 1:VLIDORT_FixIn%Cont%TS_NLAYERS, v) = &
                        GC_assa_Jacobians(w, 1:VLIDORT_FixIn%Cont%TS_NLAYERS, v) &
                        / aer_ssalbs(w, 1:VLIDORT_FixIn%Cont%TS_NLAYERS)
                END WHERE
             END IF
          END IF
       END DO
    END IF
    
    IF ( do_assa_Jacobians .AND. do_QU_Jacobians) THEN
       q = assawfidx
       DO v = 1, VLIDORT_Out%Main%TS_N_GEOMETRIES
          IF (do_aer_columnwf) THEN
             total_Qwf = 0.0d0; total_Uwf = 0.0d0; total_aersca = 0.0d0; total_aertau = 0.0d0
             DO il = 1, VLIDORT_FixIn%Cont%TS_NLAYERS
                IF (aer_flags(il)) THEN
                   total_Qwf = total_Qwf + VLIDORT_LINOUT%PROF%TS_PROFILEWF(q,il,1,v,2,didx) 
                   total_Uwf = total_Uwf + VLIDORT_LINOUT%PROF%TS_PROFILEWF(q,il,1,v,3,didx) 
                   total_aersca = total_aersca + aer_opdeps(w, il) * aer_ssalbs(w, il) 
                   total_aertau = total_aertau + aer_opdeps(w, il) 
                END IF
             END DO
             total_vaerssa = total_aersca / total_aertau
             GC_assa_QJacobians(w, 1, v) = total_Qwf     ! Normalized output
             GC_assa_UJacobians(w, 1, v) = total_Uwf     ! Normalized output
             IF ( .NOT. do_normalized_WFoutput ) THEN    ! Non-normalized output
                GC_assa_QJacobians(w, 1, v) = total_Qwf / total_vaerssa
                GC_assa_UJacobians(w, 1, v) = total_Uwf / total_vaerssa
             END IF
          ELSE
             GC_assa_QJacobians(w, 1:VLIDORT_FixIn%Cont%TS_NLAYERS, v) = &
                  VLIDORT_LINOUT%PROF%TS_PROFILEWF(q,1:VLIDORT_FixIn%Cont%TS_NLAYERS,1,v,2,didx)
             GC_assa_UJacobians(w, 1:VLIDORT_FixIn%Cont%TS_NLAYERS, v) = &
                  VLIDORT_LINOUT%PROF%TS_PROFILEWF(q,1:VLIDORT_FixIn%Cont%TS_NLAYERS,1,v,3,didx)
             IF ( .NOT. do_normalized_WFoutput ) THEN     
                WHERE( aer_flags(1:VLIDORT_FixIn%Cont%TS_NLAYERS) )
                   GC_assa_QJacobians(w, 1:VLIDORT_FixIn%Cont%TS_NLAYERS, v) = &
                        GC_assa_QJacobians(w, 1:VLIDORT_FixIn%Cont%TS_NLAYERS, v) &
                        / aer_ssalbs(w, 1:VLIDORT_FixIn%Cont%TS_NLAYERS)
                   GC_assa_UJacobians(w, 1:VLIDORT_FixIn%Cont%TS_NLAYERS, v) = &
                        GC_assa_UJacobians(w, 1:VLIDORT_FixIn%Cont%TS_NLAYERS, v) &
                        / aer_ssalbs(w, 1:VLIDORT_FixIn%Cont%TS_NLAYERS)
                END WHERE
             END IF
          END IF
       END DO
    END IF
    
    ! --------------------------------------
    ! Save cloud optical thickness Jacobians
    ! --------------------------------------
    IF ( do_cod_Jacobians ) THEN
       q = codwfidx
       DO v = 1, VLIDORT_Out%Main%TS_N_GEOMETRIES
          IF (do_cld_columnwf) THEN
             total_wf = 0.0d0
             DO il = 1, VLIDORT_FixIn%Cont%TS_NLAYERS
                IF (cld_flags(il)) THEN
                   total_wf = total_wf + VLIDORT_LINOUT%PROF%TS_PROFILEWF(q,il,1,v,1,didx)
                END IF
             END DO
             GC_cod_Jacobians(w, 1, v) = total_wf        !  Normalized output
             IF ( .NOT. do_normalized_WFoutput ) THEN    !  Non-normalized output
                GC_cod_Jacobians(w, 1, v) = total_wf / tcldtau0
             END IF
          ELSE
             GC_cod_Jacobians(w, 1:VLIDORT_FixIn%Cont%TS_NLAYERS, v) = &
                  VLIDORT_LINOUT%PROF%TS_PROFILEWF(q,1:VLIDORT_FixIn%Cont%TS_NLAYERS,1,v,1,didx)
             IF ( .NOT. do_normalized_WFoutput ) THEN     
                WHERE( cld_flags(1:VLIDORT_FixIn%Cont%TS_NLAYERS) )
                   GC_cod_Jacobians(w, 1:VLIDORT_FixIn%Cont%TS_NLAYERS, v) = &
                        GC_cod_Jacobians(w, 1:VLIDORT_FixIn%Cont%TS_NLAYERS, v) &
                        / tcld_profile(1:VLIDORT_FixIn%Cont%TS_NLAYERS)
                END WHERE
             END IF
          END IF
       END DO
    END IF
    
    IF (  do_cod_Jacobians .AND. do_QU_Jacobians) THEN
       q = codwfidx
       DO v = 1, VLIDORT_Out%Main%TS_N_GEOMETRIES
          IF (do_cld_columnwf) THEN
             total_Qwf = 0.0d0; total_Uwf = 0.0d0
             DO il = 1, VLIDORT_FixIn%Cont%TS_NLAYERS
                IF (cld_flags(il)) THEN
                   total_Qwf = total_Qwf + VLIDORT_LINOUT%PROF%TS_PROFILEWF(q,il,1,v,2,didx)
                   total_Uwf = total_Uwf + VLIDORT_LINOUT%PROF%TS_PROFILEWF(q,il,1,v,3,didx)
                END IF
             END DO
             GC_cod_QJacobians(w, 1, v) = total_Qwf      !  Normalized output
             GC_cod_UJacobians(w, 1, v) = total_Uwf    
             IF ( .NOT. do_normalized_WFoutput ) THEN    !  Non-normalized output
                GC_cod_QJacobians(w, 1, v) = total_Qwf / tcldtau0 
                GC_cod_UJacobians(w, 1, v) = total_Uwf / tcldtau0 
             END IF
          ELSE
             GC_cod_QJacobians(w, 1:VLIDORT_FixIn%Cont%TS_NLAYERS, v) = &
                  VLIDORT_LINOUT%PROF%TS_PROFILEWF(q,1:VLIDORT_FixIn%Cont%TS_NLAYERS,1,v,2,didx)
             GC_cod_UJacobians(w, 1:VLIDORT_FixIn%Cont%TS_NLAYERS, v) = &
                  VLIDORT_LINOUT%PROF%TS_PROFILEWF(q,1:VLIDORT_FixIn%Cont%TS_NLAYERS,1,v,3,didx)
             IF ( .NOT. do_normalized_WFoutput ) THEN     
                WHERE( cld_flags(1:VLIDORT_FixIn%Cont%TS_NLAYERS) )
                   GC_cod_QJacobians(w, 1:VLIDORT_FixIn%Cont%TS_NLAYERS, v) =    &
                        GC_cod_QJacobians(w, 1:VLIDORT_FixIn%Cont%TS_NLAYERS, v) &
                        / tcld_profile(1:VLIDORT_FixIn%Cont%TS_NLAYERS)
                   GC_cod_UJacobians(w, 1:VLIDORT_FixIn%Cont%TS_NLAYERS, v) =    &
                        GC_cod_UJacobians(w, 1:VLIDORT_FixIn%Cont%TS_NLAYERS, v) &
                        / tcld_profile(1:VLIDORT_FixIn%Cont%TS_NLAYERS)
                END WHERE
             END IF
          END IF
       END DO
    END IF
    
    ! ---------------------------------------------
    ! Save cloud single scattering albedo Jacobians
    ! ---------------------------------------------
    IF ( do_cssa_Jacobians ) THEN
       q = cssawfidx
       DO v = 1, VLIDORT_Out%Main%TS_N_GEOMETRIES
          IF (do_cld_columnwf) THEN
             total_wf = 0.0d0; total_cldsca = 0.0d0; total_cldtau = 0.0d0
             DO il = 1, VLIDORT_FixIn%Cont%TS_NLAYERS
                IF (cld_flags(il)) THEN
                   total_wf = total_wf + VLIDORT_LINOUT%PROF%TS_PROFILEWF(q,il,1,v,1,didx) 
                   total_cldsca = total_cldsca + cld_opdeps(w, il)  * cld_ssalbs(w, il) 
                   total_cldtau = total_cldtau + cld_opdeps(w, il) 
                END IF
             END DO
             total_vcldssa = total_cldsca / total_cldtau
             GC_cssa_Jacobians(w, 1, v) = total_wf        ! Normalized output
             IF ( .NOT. do_normalized_WFoutput ) THEN     ! Non-normalized output
                GC_cssa_Jacobians(w, 1, v) = total_wf /  total_vcldssa
             END IF
          ELSE
             GC_cssa_Jacobians(w, 1:VLIDORT_FixIn%Cont%TS_NLAYERS, v) = &
                  VLIDORT_LINOUT%PROF%TS_PROFILEWF(q,1:VLIDORT_FixIn%Cont%TS_NLAYERS,1,v,1,didx)
             IF ( .NOT. do_normalized_WFoutput ) THEN     
                WHERE( cld_flags(1:VLIDORT_FixIn%Cont%TS_NLAYERS) )
                   GC_cssa_Jacobians(w, 1:VLIDORT_FixIn%Cont%TS_NLAYERS, v) =    &
                        GC_cssa_Jacobians(w, 1:VLIDORT_FixIn%Cont%TS_NLAYERS, v) &
                        / cld_ssalbs(w, 1:VLIDORT_FixIn%Cont%TS_NLAYERS)
                END WHERE
             END IF
          END IF
       END DO
    END IF
    
    IF (  do_cssa_Jacobians .AND. do_QU_Jacobians) THEN
       q = cssawfidx
       DO v = 1, VLIDORT_Out%Main%TS_N_GEOMETRIES
          IF (do_cld_columnwf) THEN
             total_Qwf = 0.0d0; total_Uwf = 0.0d0; total_cldsca = 0.0d0; total_cldtau = 0.0d0
             DO il = 1, VLIDORT_FixIn%Cont%TS_NLAYERS
                IF (cld_flags(il)) THEN
                   total_Qwf = total_Qwf + VLIDORT_LINOUT%PROF%TS_PROFILEWF(q,il,1,v,2,didx) 
                   total_Uwf = total_Uwf + VLIDORT_LINOUT%PROF%TS_PROFILEWF(q,il,1,v,3,didx) 
                   total_cldsca = total_cldsca + cld_opdeps(w, il) * cld_ssalbs(w, il) 
                   total_cldtau = total_cldtau + cld_opdeps(w, il) 
                END IF
             END DO
             total_vcldssa = total_cldsca / total_cldtau
             GC_cssa_QJacobians(w, 1, v) = total_Qwf     ! Normalized output
             GC_cssa_UJacobians(w, 1, v) = total_Uwf     ! Normalized output
             IF ( .NOT. do_normalized_WFoutput ) THEN    ! Non-normalized output
                GC_cssa_QJacobians(w, 1, v) = total_Qwf / total_vcldssa
                GC_cssa_UJacobians(w, 1, v) = total_Uwf / total_vcldssa
             END IF
          ELSE
             GC_cssa_QJacobians(w, 1:VLIDORT_FixIn%Cont%TS_NLAYERS, v) = &
                  VLIDORT_LINOUT%PROF%TS_PROFILEWF(q,1:VLIDORT_FixIn%Cont%TS_NLAYERS,1,v,2,didx)
             GC_cssa_UJacobians(w, 1:VLIDORT_FixIn%Cont%TS_NLAYERS, v) = &
                  VLIDORT_LINOUT%PROF%TS_PROFILEWF(q,1:VLIDORT_FixIn%Cont%TS_NLAYERS,1,v,3,didx)
             IF ( .NOT. do_normalized_WFoutput ) THEN     
                WHERE( cld_flags(1:VLIDORT_FixIn%Cont%TS_NLAYERS) )
                   GC_cssa_QJacobians(w, 1:VLIDORT_FixIn%Cont%TS_NLAYERS, v) =    &
                        GC_cssa_QJacobians(w, 1:VLIDORT_FixIn%Cont%TS_NLAYERS, v) &
                        / cld_ssalbs(w, 1:VLIDORT_FixIn%Cont%TS_NLAYERS)
                   GC_cssa_UJacobians(w, 1:VLIDORT_FixIn%Cont%TS_NLAYERS, v) =    &
                        GC_cssa_UJacobians(w, 1:VLIDORT_FixIn%Cont%TS_NLAYERS, v) &
                        / cld_ssalbs(w, 1:VLIDORT_FixIn%Cont%TS_NLAYERS)
                END WHERE
             END IF
          END IF
       END DO
    END IF
    
    IF (VLIDORT_FixIn%Bool%TS_DO_LAMBERTIAN_SURFACE) THEN
       
       ! ----------------------------------------------------------------------
       ! Save Surface albedo Jacobian, *** Note that: always unnormalized  ****
       ! ----------------------------------------------------------------------
       IF ( do_Jacobians ) THEN
          DO v = 1, VLIDORT_Out%Main%TS_N_GEOMETRIES
             IF ( do_normalized_WFoutput ) THEN   !  Normalized output
                GC_Surfalbedo_Jacobians(w,v) = &
                     VLIDORT_LINOUT%SURF%TS_SURFACEWF(1,1,v,1,didx) * ground_ler(w) 
             ELSE                                 !  Non-normalized output
                GC_Surfalbedo_Jacobians(w,v) = &
                     VLIDORT_LINOUT%SURF%TS_SURFACEWF(1,1,v,1,didx) 
             END IF
          END DO
       END IF
       
       ! ----------------------------------------------------------------------------------
       ! Save Surface albedo Jacobian for Q and U, *** Note that: always unnormalized  ****
       ! ----------------------------------------------------------------------------------
       IF ( do_QU_Jacobians ) THEN
          DO v = 1, VLIDORT_Out%Main%TS_N_GEOMETRIES
             IF ( do_normalized_WFoutput ) THEN   !  Normalized output
                GC_Surfalbedo_QJacobians(w,v) = &
                     VLIDORT_LINOUT%SURF%TS_SURFACEWF(1,1,v,2,didx) * ground_ler(w) 
                GC_Surfalbedo_UJacobians(w,v) = &
                     VLIDORT_LINOUT%SURF%TS_SURFACEWF(1,1,v,3,didx) * ground_ler(w) 
             ELSE                                !  Non-normalized output
                GC_Surfalbedo_QJacobians(w,v) = &
                     VLIDORT_LINOUT%SURF%TS_SURFACEWF(1,1,v,2,didx)  
                GC_Surfalbedo_UJacobians(w,v) = &
                     VLIDORT_LINOUT%SURF%TS_SURFACEWF(1,1,v,3,didx)  
             END IF
          END DO
       END IF
       
    ELSE
       
       ! -------------------------------------------------------------
       ! Save wind speed Jacobian - VLIDORT output ALWAYS unnormalized
       ! -------------------------------------------------------------
       IF ( do_Jacobians ) THEN
          DO v = 1, VLIDORT_Out%Main%TS_N_GEOMETRIES
             IF ( do_normalized_WFoutput ) THEN  !  Normalized output
                GC_Windspeed_Jacobians(w,v) = &
                     wind_speed * VLIDORT_LINOUT%SURF%TS_SURFACEWF(1,1,v,1,didx)
             ELSE                                !  Non-normalized output
                GC_Windspeed_Jacobians(w,v) = &
                     VLIDORT_LINOUT%SURF%TS_SURFACEWF(1,1,v,1,didx)
             END IF
          END DO
       END IF
       
       ! -------------------------------------------------------------------------
       ! Save wind speed Jacobian for Q and U - VLIDORT output ALWAYS unnormalized
       ! -------------------------------------------------------------------------
       IF ( do_QU_Jacobians ) THEN
          DO v = 1, VLIDORT_Out%Main%TS_N_GEOMETRIES
             IF ( do_normalized_WFoutput ) THEN   !  Normalized output
                GC_Windspeed_QJacobians(w,v) = &
                     wind_speed * VLIDORT_LINOUT%SURF%TS_SURFACEWF(1,1,v,2,didx)
                GC_Windspeed_UJacobians(w,v) = &
                     wind_speed * VLIDORT_LINOUT%SURF%TS_SURFACEWF(1,1,v,3,didx)
             ELSE                                 !  Non-normalized output
                GC_Windspeed_QJacobians(w,v) = VLIDORT_LINOUT%SURF%TS_SURFACEWF(1,1,v,2,didx)
                GC_Windspeed_UJacobians(w,v) = VLIDORT_LINOUT%SURF%TS_SURFACEWF(1,1,v,3,didx)
             END IF
          END DO
       END IF
       
    END IF
    
    ! ------------------------------
    ! Save surface pressure Jacobian
    ! ------------------------------ 
    IF ( do_sfcprs_Jacobians ) THEN
       q = sfcprswfidx
       DO v = 1, VLIDORT_Out%Main%TS_N_GEOMETRIES
          deltp = pressures(VLIDORT_FixIn%Cont%TS_NLAYERS) - &
                  pressures(VLIDORT_FixIn%Cont%TS_NLAYERS-1)
          IF ( do_normalized_WFoutput ) THEN         !  Normalized output
             GC_sfcprs_Jacobians(w,v) = &
                  VLIDORT_LINOUT%PROF%TS_PROFILEWF(q,VLIDORT_FixIn%Cont%TS_NLAYERS,1,v,1,didx) 
          ELSE                                       !  Non-normalized output
             GC_sfcprs_Jacobians(w,v) = &
                  VLIDORT_LINOUT%PROF%TS_PROFILEWF(q,VLIDORT_FixIn%Cont%TS_NLAYERS,1,v,1,didx) &
                  / deltp
          END IF
       END DO
    END IF
    
    ! ------------------------------------------
    ! Save surface pressure Jacobian for Q and U
    ! ------------------------------------------
    IF ( do_sfcprs_Jacobians .AND. do_QU_Jacobians ) THEN
       q = sfcprswfidx
       DO v = 1, VLIDORT_Out%Main%TS_N_GEOMETRIES
          deltp = pressures(VLIDORT_FixIn%Cont%TS_NLAYERS) - &
                  pressures(VLIDORT_FixIn%Cont%TS_NLAYERS-1)
          IF ( do_normalized_WFoutput ) THEN        !  Normalized output
             GC_sfcprs_QJacobians(w,v) = &
                  VLIDORT_LINOUT%PROF%TS_PROFILEWF(q,VLIDORT_FixIn%Cont%TS_NLAYERS,1,v,2,didx) 
             GC_sfcprs_UJacobians(w,v) = &
                  VLIDORT_LINOUT%PROF%TS_PROFILEWF(q,VLIDORT_FixIn%Cont%TS_NLAYERS,1,v,3,didx) 
          ELSE                                      !  Non-normalized output
             GC_sfcprs_QJacobians(w,v) = &
                  VLIDORT_LINOUT%PROF%TS_PROFILEWF(q,VLIDORT_FixIn%Cont%TS_NLAYERS,1,v,2,didx) &
                  / deltp
             GC_sfcprs_UJacobians(w,v) = &
                  VLIDORT_LINOUT%PROF%TS_PROFILEWF(q,VLIDORT_FixIn%Cont%TS_NLAYERS,1,v,3,didx) &
                  / deltp
          END IF
       END DO
    END IF

    ! -----------------------------
    ! Space for future BRDF results
    ! -----------------------------
    
  END SUBROUTINE save_results
  
  SUBROUTINE Vlidort_set_optical (error)

    IMPLICIT NONE

    ! ------------------
    ! Modified variables
    ! ------------------
    LOGICAL,       INTENT(INOUT) :: error ! Error variable

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER                   :: nw1
    CHARACTER(LEN=max_ch_len) :: tmperrmessage

    ! ----------------
    ! Code starts here
    ! ----------------
    error = .FALSE.

    ! -------------------------------------------------
    ! FLUX_FACTOR = solar_cspec_data(nlambdas-w+1)*1.d4
    ! Always apply a Flux Factor here
    ! -------------------------------------------------
    IF (do_normalized_radiance) then
       VLIDORT_FixIn%Sunrays%TS_FLUX_FACTOR = 1.0d0
    ELSE 
       IF (do_effcrs) then
          VLIDORT_FixIn%Sunrays%TS_FLUX_FACTOR = solar_cspec_data(w)
       ELSE
          VLIDORT_FixIn%Sunrays%TS_FLUX_FACTOR = solar_spec_data(w)
       ENDIF
    ENDIF
    
    ! ---------------------------------------------------
    ! VLIDORT variable | ** Initialize optical properties
    ! ---------------------------------------------------
    VLIDORT_FixIn%Optical%TS_GREEKMAT_TOTAL_INPUT(0:VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT, &
         1:GC_nlayers, :) = 0.0d0
    VLIDORT_FixIn%Optical%TS_DELTAU_VERT_INPUT(1:GC_nlayers)       = 0.0d0
    VLIDORT_ModIn%MOptical%TS_OMEGA_TOTAL_INPUT(1:GC_nlayers)      = 0.0d0
    
    ! --------------------------------------------------------------
    ! VLIDORT variable | ** Initialize linearized optical properties
    ! --------------------------------------------------------------
    IF (VLIDORT_LinModIn%MCont%TS_DO_ATMOS_LINEARIZATION) THEN
       VLIDORT_LinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT(1:n_totalprofile_wfs_wcld, &
            0:VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT, 1:GC_nlayers, :)                     = 0.0d0
       VLIDORT_LinFixIn%Optical%TS_L_DELTAU_VERT_INPUT(1:n_totalprofile_wfs_wcld, 1:GC_nlayers) = 0.0d0
       VLIDORT_LinFixIn%Optical%TS_L_OMEGA_TOTAL_INPUT(1:n_totalprofile_wfs_wcld, 1:GC_nlayers) = 0.0d0
    END IF
    
    ! ---------------------------
    ! Rayleigh phase matrix input
    ! ---------------------------
    depol  = Rayleigh_depols(w)
    beta_2 = (1.0d0 - depol) / (2.0d0 + depol) 
    pRay_20 =1.0d0
    pRay_22 = beta_2
    IF ( VLIDORT_FixIn%Cont%TS_NSTOKES .GT. 1 ) THEN
       pRay_12 = 6.0d0 * beta_2
       pRay_32 = -SQRT(6.0d0) * beta_2
       pRay_41 = 3.0d0 * (1.0d0 - 2.0d0*depol) / (2.0d0 + depol) 
    END IF
    
    phasmoms_input(0, 1, 1) = 1.0d0
    phasmoms_input(2, 1, 1) = pRay_22  
    IF (VLIDORT_FixIn%Cont%TS_NSTOKES > 1) THEN
       phasmoms_input(2, 2, 1) = pRay_12
       phasmoms_input(2, 5, 1) = pRay_32
       phasmoms_input(1, 4, 1) = pRay_41
    END IF
    
    ! ------------------------------------------------------------------
    ! Get aerosol optical properties (use 1 aerosol type for all layers)
    ! ------------------------------------------------------------------
    IF (do_aerosols) THEN
       nw1=1
       CALL prepare_aercld_optical_property(database_dir, maxaer, naer, aer_types, GC_nlayers, &
            aer_profile(1:maxaer, 1:GC_nlayers), aer_flags(1:GC_nlayers), lambdas(w), nw1,     &
            aer_reflambda, VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT, ngksec,                &
            aer_opdeps(w, 1:GC_nlayers), aer_ssalbs(w, 1:GC_nlayers),                          &
            aer_relqext(w, 1:GC_nlayers),                                                      &
            aer_phfcn(1:GC_nlayers, 0:VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT, 1:ngksec),  &
            tmperrmessage, error)

       IF (error) CALL write_err_message (.FALSE., tmperrmessage)
       
    END IF
    
    ! ----------------------------
    ! Get cloud optical properties
    ! ----------------------------
    IF (do_clouds .AND. .NOT. do_lambertian_cld) THEN
       nw1=1
       CALL prepare_aercld_optical_property(database_dir, maxcld, ncld, cld_types, GC_nlayers, &
            cld_profile(1:maxcld, 1:GC_nlayers), cld_flags(1:GC_nlayers), lambdas(w), nw1,     &
            cld_reflambda, VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT, ngksec,                &
            cld_opdeps(w, 1:GC_nlayers), cld_ssalbs(w, 1:GC_nlayers),                          &
            cld_relqext(w, 1:GC_nlayers),                                                      &
            cld_phfcn(1:GC_nlayers, 0:VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT, 1:ngksec),  &
            tmperrmessage, error) 
       
       IF (error) CALL write_err_message (.FALSE., tmperrmessage)

    END IF

!!$    IF (.NOT. VLIDORT_FixIn%Bool%TS_DO_LAMBERTIAN_SURFACE) THEN
!!$       
!!$       ! ------------------------------------------------
!!$       ! VLIDORT variable | ** !  set the BRDF parameters
!!$       ! ------------------------------------------------
!!$       VBRDF_Sup_In%BS_BRDF_PARAMETERS(1,2) = water_rn(w)
!!$       VBRDF_Sup_In%BS_BRDF_PARAMETERS(1,3) = water_cn(w)
!!$       
!!$       !         !  Get BRDF and linearized BRDF functions         
!!$       !         CALL VLIDORT_BRDF_LS_MASTER                                        &
!!$       !              ( DO_USER_VZANGLES, DO_SHADOW_EFFECT, DO_COXMUNK_DBMS,        & ! Inputs
!!$       !              DO_SURFACE_EMISSION, .FALSE., NSTOKES,                        & ! Inputs
!!$       !              2*NSTREAMS-1, N_BRDF_KERNELS, WHICH_BRDF,                     & ! Inputs
!!$       !              LAMBERTIAN_KERNEL_FLAG, NSTREAMS_BRDF, BRDF_FACTORS,          & ! Inputs
!!$       !              N_BRDF_PARAMETERS, BRDF_PARAMETERS, BRDF_NAMES,               & ! Inputs
!!$       !              DO_KERNEL_PARAMS_WFS, DO_KERNEL_FACTOR_WFS, DO_KPARAMS_DERIVS,& ! Inputs
!!$       !              N_KERNEL_FACTOR_WFS, N_KERNEL_PARAMS_WFS, N_SURFACE_WFS,      & ! Inputs
!!$       !              N_SZANGLES, NSTREAMS, N_USER_VZANGLES, N_USER_RELAZMS,        & ! Inputs
!!$       !              SZANGLES, USER_VZANGLES, USER_RELAZMS,                        & ! Inputs
!!$       !              BRDF_F_0, BRDF_F, USER_BRDF_F_0, USER_BRDF_F,                 & ! Outputs
!!$       !              LS_BRDF_F_0, LS_BRDF_F, LS_USER_BRDF_F_0, LS_USER_BRDF_F,     & ! Outputs
!!$       !              EXACTDB_BRDFUNC, LS_EXACTDB_BRDFUNC,                          & ! Outputs
!!$       !              EMISSIVITY, USER_EMISSIVITY, LS_EMISSIVITY, LS_USER_EMISSIVITY )   ! Outputs
!!$    END IF    
    
  END SUBROUTINE Vlidort_set_optical

  SUBROUTINE Vlidort_cloud_and_calculation (error)

    USE VLIDORT_LPS_MASTERS
    USE VLIDORT_AUX

    IMPLICIT NONE
    
    ! ------------------
    ! Modified variables
    ! ------------------
    LOGICAL,       INTENT(INOUT) :: error ! Error variable

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER                   :: ic, n, q, g, j, k

    ! ----------------
    ! Code starts here
    ! ----------------
    error = .FALSE.
    
    ! ==============================================================
    !            Loop for clear and cloudy part
    !  Set up cloud dependent vlidort options and optical properties
    !  Weight results based Indepdent Pixel Approximation (IPA)
    ! ==============================================================
    ! ---------------------------------------------------------------------
    ! Initialize output: stokes, profile & surface jacobians at each lambda
    ! ---------------------------------------------------------------------
    stokes_clrcld(1:VLIDORT_Out%Main%TS_N_GEOMETRIES, 1:VLIDORT_FixIn%Cont%TS_NSTOKES, :) = 0.0d0
    profilewf_sum(1:n_totalprofile_wfs_wcld, 1:GC_nlayers, 1:VLIDORT_Out%Main%TS_N_GEOMETRIES, &
         1:VLIDORT_FixIn%Cont%TS_NSTOKES)                                        = 0.0d0
    surfacewf_clrcld(1:VBRDF_LinSup_in%BS_N_SURFACE_WFS, 1:VLIDORT_Out%Main%TS_N_GEOMETRIES,   &
         1:VLIDORT_FixIn%Cont%TS_NSTOKES, :)                                     = 0.0d0

                 ! --------------------------
    DO ic = 1, 2 ! Beging cloud fraction loop
                 ! --------------------------
       VLIDORT_FixIn%Cont%TS_NLAYERS = GC_nlayers
       
       IF (ic == 1)  THEN     ! Clear
          ipafrac = 1.0d0 - cfrac
       ELSE
          ipafrac = cfrac     ! Cloudy part
          IF (do_lambertian_cld) VLIDORT_FixIn%Cont%TS_NLAYERS = cld_uppers(1) - 1
       END IF
       
       IF (ipafrac == 0.0d0) CYCLE
       
       ! ---------------------------
       ! Set lambertian cloud albedo
       ! ---------------------------
!!$       IF ( VLIDORT_FixIn%Bool%TS_DO_LAMBERTIAN_SURFACE ) THEN
       IF (ic == 2 .AND. do_lambertian_cld) THEN
          VLIDORT_FixIn%Optical%TS_LAMBERTIAN_ALBEDO  = lambertian_cldalb
          VLIDORT_FixIn%Bool%TS_DO_LAMBERTIAN_SURFACE = .TRUE.
!!$       ENDIF
       ELSE
          VLIDORT_FixIn%Optical%TS_LAMBERTIAN_ALBEDO = ground_ler(w)
       END IF
!!$          VBRDF_Sup_In%BS_BRDF_FACTORS(1) = VLIDORT_FixIn%Optical%TS_LAMBERTIAN_ALBEDO
!!$    END IF
       
       ! ------------------------------------------------------------------
       ! Reset Rayleigh only flag, set if no aerosols and scattering clouds
       ! ------------------------------------------------------------------
       IF (.NOT. do_aerosols .AND. .NOT. do_lambertian_cld) THEN
          IF (ic == 1) THEN
             VLIDORT_ModIn%MBool%TS_DO_RAYLEIGH_ONLY     = .TRUE.
             VLIDORT_ModIn%MBool%TS_DO_DOUBLE_CONVTEST   = .FALSE.
             VLIDORT_ModIn%MBool%TS_DO_DELTAM_SCALING    = .FALSE.
             VLIDORT_ModIn%MBool%TS_DO_SOLUTION_SAVING   = .FALSE.
             VLIDORT_ModIn%MBool%TS_DO_BVP_TELESCOPING   = .FALSE.
             VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT = 2
          ELSE
             VLIDORT_ModIn%MBool%TS_DO_RAYLEIGH_ONLY     = .FALSE.
             VLIDORT_ModIn%MBool%TS_DO_DOUBLE_CONVTEST   = .TRUE.
             VLIDORT_ModIn%MBool%TS_DO_DELTAM_SCALING    = .TRUE.
             VLIDORT_FixIn%Bool%TS_DO_SSCORR_TRUNCATION  = .FALSE. ! Change after Xiong settings
             VLIDORT_ModIn%MBool%TS_DO_SOLUTION_SAVING   = .FALSE.
             VLIDORT_ModIn%MBool%TS_DO_BVP_TELESCOPING   = .FALSE.
             VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT = aercld_nmoments_input
          END IF
       END IF
       
       IF (do_Jacobians) THEN
          IF (ic == 1) THEN
             VLIDORT_LinFixIn%Cont%TS_N_TOTALPROFILE_WFS = N_TOTALPROFILE_WFS_ncld
             VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_NUMBER(1:VLIDORT_FixIn%Cont%TS_NLAYERS) = &
                  N_TOTALPROFILE_WFS_ncld
          ELSE
             VLIDORT_LinFixIn%Cont%TS_N_TOTALPROFILE_WFS = N_TOTALPROFILE_WFS_wcld
             VLIDORT_LinFixIn%Cont%TS_LAYER_VARY_NUMBER(1:VLIDORT_FixIn%Cont%TS_NLAYERS) = &
                  N_TOTALPROFILE_WFS_wcld
          END IF
       END IF
       
       ! ----------------------------------------------------------
       ! Here, we just make the Bulk  properties at each wavelength
       ! Need to derive optical properties for ic = 1 or ic = 2 
       ! .and. ipafrac = 1.0) .and. scattering clouds
       ! ----------------------------------------------------------
         
       IF (ic == 1 .OR. (ic == 2 .AND. ipafrac >= 1.0) &
            .OR. (ic == 2 .AND. .NOT. do_lambertian_cld) ) THEN
          
                               ! -----------
          DO n = 1, GC_nlayers ! Layers loop
                               ! -----------

             ! -----------------------------------------
             ! Temperature factors for O3 cross-sections
             ! -----------------------------------------
             temp    = mid_temperatures(n) - 273.15d0
             temp_sq = temp * temp
             
             ! ---------------------
             ! Trace gas absorptions
             ! ---------------------
             total_molabs = 0.0d0
             DO g = 1, ngases
                gasabs(n,g) = 0.0d0
                IF ( which_gases(g) == 'O3  ' .AND. gas_xsecs_type(g) == 2 ) THEN
                   xsec = gas_xsecs(w,1,g) + temp * o3c1_xsecs(w) + temp_sq * o3c2_xsecs(w)
                   xsec_o3save(n) = xsec
                   gasabs(n,g) = gas_partialcolumns(n,g) * xsec
                ELSE IF (gas_xsecs_type(g) == 2) THEN
                   xsec = gas_xsecs(w,1,g) + temp * gas_xsecs(w,2,g) + temp_sq * gas_xsecs(w,3,g)
                   gasabs(n,g) = gas_partialcolumns(n,g) * xsec
                ELSE IF (gas_xsecs_type(g) == 3) THEN
                   gasabs(n,g) = gas_partialcolumns(n,g) * gas_xsecs(w,n,g)
                ELSE IF (gas_xsecs_type(g) == 1) THEN
                   gasabs(n,g) = gas_partialcolumns(n,g) * gas_xsecs(w,1,g)
                END IF
                total_molabs = total_molabs + gasabs(n,g)
             END DO
             total_gasabs(n) = total_molabs
             
             ! ---------------------------------
             ! Rayleigh scattering optical depth
             ! ---------------------------------            
             total_molsca = aircolumns(n) * Rayleigh_xsecs(w)
             
             ! -----------
             ! Total value
             ! -----------            
             total_moltau = total_molabs + total_molsca
             
             nscatter = 1
             scaco_input(1) = total_molsca
             total_tau = total_moltau
             total_sca = total_molsca
             
             IF ( aer_flags(n) ) THEN
                nscatter = nscatter + 1
                aerscaidx = nscatter
                total_aertau = aer_opdeps(w, n)
                total_aersca = total_aertau * aer_ssalbs(w, n)
                total_tau = total_tau + total_aertau
                total_sca = total_sca + total_aersca
                scaco_input(nscatter) = total_aersca
                phasmoms_input(0:VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT, 1:ngksec, nscatter) = &
                     aer_phfcn(n, 0:VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT, 1:ngksec)
             ENDIF
             
             ! --------------------------------------------
             ! Only for cloudy parts with scattering clouds
             ! --------------------------------------------
             IF ( cld_flags(n) .AND. ic == 2) THEN
                nscatter = nscatter + 1
                cldscaidx = nscatter
                total_cldtau = cld_opdeps(w, n)
                total_cldsca = total_cldtau * cld_ssalbs(w, n)
                total_tau = total_tau + total_cldtau
                total_sca = total_sca + total_cldsca
                scaco_input(nscatter) = total_cldsca
                phasmoms_input(0:VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT, 1:ngksec, nscatter) = &
                     cld_phfcn(n, 0:VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT, 1:ngksec)
             END IF
             
             omega = total_sca / total_tau
             opdeps(w, n) = total_tau
             ssalbs(w, n) = omega
             IF (omega .GT. 1.0 - smallnum ) omega = 1.0d0 - smallnum
             IF (omega .LT. smallnum) omega = smallnum
             VLIDORT_ModIn%MOptical%TS_OMEGA_TOTAL_INPUT(n) = omega
             VLIDORT_FixIn%Optical%TS_DELTAU_VERT_INPUT(n)  = total_tau
             
             DO j = 0, VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT
                DO k = 1, ngksec
                   phasmoms_total_input(j, k) = SUM(phasmoms_input(j, k, 1:nscatter) &
                        * scaco_input(1:nscatter)) / total_sca
                END DO
             END DO
             
             VLIDORT_FixIn%Optical%TS_GREEKMAT_TOTAL_INPUT                                 &
                  (0:VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT, n,                       &
                  greekmat_idxs(1:nactgkmatc)) =                                           &
                  phasmoms_total_input(0:VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT,      &
                  phasmoms_idxs(1:nactgkmatc))

!xliu, 6/6/2013, based on discussion with Rob and Vjay/ gga 1/24/2014
!!$             IF ( nactgkmatc > 1 ) VLIDORT_FixIn%Optical%TS_GREEKMAT_TOTAL_INPUT           &
!!$                  (0:VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT, n, 15) =                 &
!!$                  -VLIDORT_FixIn%Optical%TS_GREEKMAT_TOTAL_INPUT(                          &
!!$                  0:VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT, n, 15)
         
             IF ( nactgkmatc > 1 ) THEN
                VLIDORT_FixIn%Optical%TS_GREEKMAT_TOTAL_INPUT           &
                  (0:VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT, n, 2) =                 &
                  -VLIDORT_FixIn%Optical%TS_GREEKMAT_TOTAL_INPUT(                         &
                  0:VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT, n, 2)

                VLIDORT_FixIn%Optical%TS_GREEKMAT_TOTAL_INPUT           &
                  (0:VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT, n, 5) =                 &
                  -VLIDORT_FixIn%Optical%TS_GREEKMAT_TOTAL_INPUT(                          &
                  0:VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT, n, 5)

                VLIDORT_FixIn%Optical%TS_GREEKMAT_TOTAL_INPUT                              &
                  (0:VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT, n, 12) =                 &
                  -VLIDORT_FixIn%Optical%TS_GREEKMAT_TOTAL_INPUT(                          &
                  0:VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT, n, 12)
             ENDIF
    
             ! ----------------------------------------------------------------------------------
             ! This should always be 1, but may be slightly different due to numerical truncation
             ! ----------------------------------------------------------------------------------
             VLIDORT_FixIn%Optical%TS_GREEKMAT_TOTAL_INPUT(0, n, 1) = 1.d0  
             
             ! -----------------------------------------------
             ! Set linearized VLIDORT inputs
             ! First Jacobian  ( q = 1 ), use total absorption
             ! -----------------------------------------------
             IF ( do_jacobians ) THEN
                q = gaswfidx
                DO g = 1, ngases
                   IF ( which_gases(g) == 'O3  ' ) abs_o3(n) = gasabs(n,g)
                END DO
                  !ratio = abs_o3(n) / total_tau
                ratio = total_gasabs(n) / total_tau
                VLIDORT_LinFixIn%Optical%TS_L_DELTAU_VERT_INPUT(q,n) =   ratio
                VLIDORT_LinFixIn%Optical%TS_L_OMEGA_TOTAL_INPUT(q,n) = - ratio
             END IF
             
             ! ---------------------------------------------------------------
             ! Second Jacobian ( q = 2 ). Temperature. Only if flagged.
             ! L_gas is the contribution from the T-dependent O3 cross-section
             ! L_air is the contribution from the Air columns
             ! ---------------------------------------------------------------
             
             IF ( do_T_jacobians ) THEN
                q = twfidx
                DO g = 1, ngases
                   ! ------------------------------------------------------
                   ! So far only temperature dependence in O3 is considered
                   ! ------------------------------------------------------
                   IF ( which_gases(g) == 'O3  ' ) THEN
                      ratio = gasabs(n,g) / total_tau
                      xsec  = gas_xsecs(w,1,g) + temp * o3c1_xsecs(w) + temp_sq*o3c2_xsecs(w)
                      dxsec_dT = dxsec_dT + o3c1_xsecs(w) + 2.0d0 * temp * o3c2_xsecs(w)
                      L_gas    = L_gas + ratio * dxsec_dT / xsec 
                      EXIT
                   END IF
                END DO
                
                ratio = total_moltau / total_tau
                L_air = ratio * daircolumns_dT(n) / aircolumns(n)
                VLIDORT_LinFixIn%Optical%TS_L_DELTAU_VERT_INPUT(q,n) = &
                     mid_temperatures(n) * ( L_gas + L_air )
                
                ratio = total_molsca / total_sca
                L_air = mid_temperatures(n) * ratio * daircolumns_dT(n) / aircolumns(n)
                VLIDORT_LinFixIn%Optical%TS_L_OMEGA_TOTAL_INPUT(q,n) = &
                     L_air - VLIDORT_LinFixIn%Optical%TS_L_DELTAU_VERT_INPUT(q,n) 
                
                VLIDORT_LinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT(q,             &
                     0:VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT, n, :) = 0.0d0
                DO j = 0, VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT
                   DO k = 1, ngksec
                      IF (phasmoms_total_input(j, k) /= 0.d0) THEN
                         l_phasmoms_total_input(q, j, k) = ( phasmoms_input(j, k, 1) &
                              - phasmoms_total_input(j, k) ) / phasmoms_total_input(j, k) * L_air
                      ELSE
                         l_phasmoms_total_input(q, j, k) = 0.0d0
                      END IF
                   END DO
                END DO
                VLIDORT_LinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT(q,         &
                     0:VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT, n,        &
                     greekmat_idxs(1:nactgkmatc)) = l_phasmoms_total_input(q, &
                     0:VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT,           &
                     phasmoms_idxs(1:nactgkmatc))

!xliu, 6/6/2013, based on discussion with Rob and Vjay/ gga 1/24/2014
!!$                IF ( nactgkmatc > 1 )  VLIDORT_LinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT(q, &
!!$                     0:VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT, n, 15)                   &
!!$                     = -VLIDORT_LinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT(q,                &
!!$                     0:VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT, n, 15)                             

                IF ( nactgkmatc > 1 )  THEN

                   VLIDORT_LinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT(q,                     &
                     0:VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT, n, 2)                   &
                     = -VLIDORT_LinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT(q,                &
                     0:VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT, n, 2)

                   VLIDORT_LinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT(q,                     &
                     0:VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT, n, 5)                   &
                     = -VLIDORT_LinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT(q,                &
                     0:VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT, n, 5)      

                   VLIDORT_LinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT(q,                     &
                     0:VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT, n, 12)                   &
                     = -VLIDORT_LinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT(q,                &
                     0:VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT, n, 12)   
                ENDIF

             END IF
             
             ! ------------------------------------------------------------
             ! Ground surface pressure Jacobians, only if flagged
             ! Rayleigh only, decouple Rayleigh from gases, aerosol, clouds
             ! ------------------------------------------------------------
             IF (do_sfcprs_Jacobians .AND. n == GC_nlayers) THEN
                q = sfcprswfidx
                pvar = total_molsca / total_tau
                VLIDORT_LinFixIn%Optical%TS_L_OMEGA_TOTAL_INPUT(q, n) =     &
                     (1.0 - VLIDORT_ModIn%MOptical%TS_OMEGA_TOTAL_INPUT(n)) &
                     * total_molsca / total_tau / VLIDORT_ModIn%MOptical%TS_OMEGA_TOTAL_INPUT(n)
                VLIDORT_LinFixIn%Optical%TS_L_DELTAU_VERT_INPUT(q, n) = pvar
                VLIDORT_LinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT(q,   &
                     0:VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT, n, :) = 0.0d0
                DO j = 0, VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT
                   DO k = 1, ngksec
                      IF (phasmoms_total_input(j, k) /= 0.0) THEN
                         l_phasmoms_total_input(q, j, k) = (phasmoms_input(j, k, 1) - &
                              phasmoms_total_input(j, k)) / phasmoms_total_input(j, k) &
                              * total_molsca / total_sca
                      ELSE
                         l_phasmoms_total_input(q, j, k) = 0.0
                      END IF
                   END DO
                END DO
             END IF
             
             ! ------------------------------------------
             ! Aerosol optical thickness, only if flagged
             ! ------------------------------------------
             IF (do_aod_Jacobians .AND. aer_flags(n) ) THEN
                q = aodwfidx
                
                VLIDORT_LinFixIn%Optical%TS_L_DELTAU_VERT_INPUT(q, n) = &
                     + total_aertau / total_tau
                VLIDORT_LinFixIn%Optical%TS_L_OMEGA_TOTAL_INPUT(q, n) = &
                     (total_aersca / total_aertau - VLIDORT_ModIn%MOptical%TS_OMEGA_TOTAL_INPUT(n)) &
                     / VLIDORT_ModIn%MOptical%TS_OMEGA_TOTAL_INPUT(n)  * total_aertau / total_tau
                VLIDORT_LinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT(q,   &
                     0:VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT, n, :) = 0.0d0
                
                DO j = 0, VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT
                   DO k = 1, ngksec
                      IF (phasmoms_total_input(j, k) /= 0.d0) THEN
                         l_phasmoms_total_input(q, j, k) = ( phasmoms_input(j, k, aerscaidx) &
                              - phasmoms_total_input(j, k) ) / phasmoms_total_input(j, k) &
                              * total_aersca / total_sca
                      ELSE
                         l_phasmoms_total_input(q, j, k) = 0.0d0
                      END IF
                   END DO
                END DO
                
                VLIDORT_LinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT(q,                         &
                     0:VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT, n,                        &
                     greekmat_idxs(1:nactgkmatc)) =                                           &
                     l_phasmoms_total_input(q, 0:VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT, &
                     phasmoms_idxs(1:nactgkmatc))

!xliu, 6/6/2013, changes based on discussion with Rob and Vijay / gga 1/24/2014
!!$                IF ( nactgkmatc > 1 )  VLIDORT_LinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT(q, &
!!$                     0:VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT, n, 15)                   &
!!$                     = -VLIDORT_LinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT(q,                &
!!$                     0:VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT, n, 15)
                
               IF ( nactgkmatc > 1 )  THEN
                   VLIDORT_LinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT(q, &
                     0:VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT, n, 2)                   &
                     = -VLIDORT_LinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT(q,                &
                     0:VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT, n, 2)

                   VLIDORT_LinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT(q, &
                     0:VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT, n, 5)                   &
                     = -VLIDORT_LinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT(q,                &
                     0:VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT, n, 5)

                   VLIDORT_LinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT(q, &
                     0:VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT, n, 12)                   &
                     = -VLIDORT_LinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT(q,                &
                     0:VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT, n, 12)
                ENDIF
 
             END IF
             
             ! ------------------------------------------------------
             ! Aerosol SSA thickness, only if flagged (AOD unchanged)
             ! ------------------------------------------------------
             IF (do_assa_Jacobians .AND. aer_flags(n)) THEN
                q = assawfidx
                VLIDORT_LinFixIn%Optical%TS_L_DELTAU_VERT_INPUT(q,n) = 0.0d0
                VLIDORT_LinFixIn%Optical%TS_L_OMEGA_TOTAL_INPUT(q,n) = total_aersca / &
                     VLIDORT_ModIn%MOptical%TS_OMEGA_TOTAL_INPUT(n) / total_tau
                VLIDORT_LinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT(q, &
                     0:VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT, n, :) = 0.0d0
                DO j = 0, VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT
                   DO k = 1, ngksec
                      IF (phasmoms_total_input(j, k) /= 0.d0) THEN
                         l_phasmoms_total_input(q, j, k) = ( phasmoms_input(j, k, aerscaidx) &
                              - phasmoms_total_input(j, k) ) / phasmoms_total_input(j, k) &
                              * total_aersca / total_sca
                      ELSE
                         l_phasmoms_total_input(q, j, k) = 0.d0
                      END IF
                   END DO
                END DO

                VLIDORT_LinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT(q,                                 &
                     0:VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT, n, greekmat_idxs(1:nactgkmatc)) = &
                     l_phasmoms_total_input(q, 0:VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT,         &
                     phasmoms_idxs(1:nactgkmatc))

!xliu, 6/6/2013, changes based on discussion with Rob and Vijay / gga 1/24/2014
!!$                IF ( nactgkmatc > 1 )  VLIDORT_LinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT(q, &
!!$                     0:VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT, n, 15)                   &
!!$                     = - VLIDORT_LinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT(q,               &
!!$                     0:VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT, n, 15)

                IF ( nactgkmatc > 1 )  THEN
                   VLIDORT_LinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT(q, &
                     0:VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT, n, 2)                   &
                     = - VLIDORT_LinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT(q,               &
                     0:VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT, n, 2)

                   VLIDORT_LinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT(q, &
                     0:VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT, n, 5)                   &
                     = - VLIDORT_LinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT(q,               &
                     0:VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT, n, 5)

                   VLIDORT_LinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT(q, &
                     0:VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT, n, 12)                   &
                     = - VLIDORT_LinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT(q,               &
                     0:VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT, n, 12)
                ENDIF

             END IF
             
             ! ----------------------------------------
             ! Cloud optical thickness, only if flagged
             ! ----------------------------------------
             IF (do_cod_Jacobians .AND. cld_flags(n) .AND. ic == 2 ) THEN
                q = codwfidx
                
                VLIDORT_LinFixIn%Optical%TS_L_DELTAU_VERT_INPUT(q, n) = + total_cldtau / total_tau
                VLIDORT_LinFixIn%Optical%TS_L_OMEGA_TOTAL_INPUT(q, n) = (total_cldsca / total_cldtau &
                     - VLIDORT_ModIn%MOptical%TS_OMEGA_TOTAL_INPUT(n))                               &
                     / VLIDORT_ModIn%MOptical%TS_OMEGA_TOTAL_INPUT(n) * total_cldtau / total_tau
                VLIDORT_LinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT(q, &
                     0:VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT, n, :) = 0.0d0
                DO j = 0, VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT
                   DO k = 1, ngksec
                      IF (phasmoms_total_input(j, k) /= 0.d0) THEN
                         l_phasmoms_total_input(q, j, k) = ( phasmoms_input(j, k, cldscaidx) &
                              - phasmoms_total_input(j, k) ) / phasmoms_total_input(j, k)    &
                              * total_cldsca / total_sca
                      ELSE
                         l_phasmoms_total_input(q, j, k) = 0.d0
                      END IF
                   END DO
                END DO
                VLIDORT_LinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT(q,                         &
                     0:VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT, n,                        &
                     greekmat_idxs(1:nactgkmatc)) =                                           &
                     l_phasmoms_total_input(q, 0:VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT, &
                     phasmoms_idxs(1:nactgkmatc))

!xliu, 6/6/2013, changes based on discussion with Rob and Vijay / gga 1/24/2014
!!$                IF ( nactgkmatc > 1 )  VLIDORT_LinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT(q,  &
!!$                     0:VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT, n, 15)                    &
!!$                     = -VLIDORT_LinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT(q,                 &
!!$                     0:VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT, n, 15)

                IF ( nactgkmatc > 1 ) THEN
                   VLIDORT_LinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT(q,  &
                     0:VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT, n, 2)                    &
                     = -VLIDORT_LinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT(q,                 &
                     0:VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT, n, 2)

                   VLIDORT_LinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT(q,  &
                     0:VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT, n, 5)                    &
                     = -VLIDORT_LinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT(q,                 &
                     0:VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT, n, 5)

                   VLIDORT_LinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT(q,  &
                     0:VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT, n, 12)                    &
                     = -VLIDORT_LinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT(q,                 &
                     0:VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT, n, 12)
                ENDIF
                
             END IF
             
             ! ----------------------------------------------------
             ! Cloud SSA thickness, only if flagged (AOD unchanged)
             ! ----------------------------------------------------
             IF (do_cssa_Jacobians .AND. cld_flags(n) .AND. ic == 2) THEN
                q = cssawfidx
                VLIDORT_LinFixIn%Optical%TS_L_DELTAU_VERT_INPUT(q,n) = 0.0d0
                VLIDORT_LinFixIn%Optical%TS_L_OMEGA_TOTAL_INPUT(q,n) = total_cldsca / &
                     VLIDORT_ModIn%MOptical%TS_OMEGA_TOTAL_INPUT(n) / total_tau
                VLIDORT_LinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT(q, &
                     0:VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT, n, :) = 0.0d0
                DO j = 0, VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT
                   DO k = 1, ngksec
                      IF (phasmoms_total_input(j, k) /= 0.d0) THEN
                         l_phasmoms_total_input(q, j, k) = ( phasmoms_input(j, k, cldscaidx) &
                              - phasmoms_total_input(j, k) ) / phasmoms_total_input(j, k) &
                              * total_cldsca / total_sca
                      ELSE
                         l_phasmoms_total_input(q, j, k) = 0.d0
                      END IF
                   END DO
                END DO
                VLIDORT_LinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT(q, &
                     0:VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT, n,&
                     greekmat_idxs(1:nactgkmatc)) =                   &
                     l_phasmoms_total_input(q,                        &
                     0:VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT,   &
                     phasmoms_idxs(1:nactgkmatc))

!xliu, 6/6/2013, changes based on discussion with Rob and Vijay / gga 1/24/2014
!!$                IF ( nactgkmatc > 1 )  VLIDORT_LinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT(q, &
!!$                     0:VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT, n, 15)                   &
!!$                     = - VLIDORT_LinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT(q,               &
!!$                     0:VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT, n, 15)

                IF ( nactgkmatc > 1 ) THEN
                   VLIDORT_LinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT(q, &
                     0:VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT, n, 2)                   &
                     = - VLIDORT_LinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT(q,               &
                     0:VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT, n, 2)

                   VLIDORT_LinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT(q, &
                     0:VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT, n, 5)                   &
                     = - VLIDORT_LinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT(q,               &
                     0:VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT, n, 5)

                   VLIDORT_LinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT(q, &
                     0:VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT, n, 12)                   &
                     = - VLIDORT_LinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT(q,               &
                     0:VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT, n, 12)
                ENDIF

             END IF
             
             VLIDORT_LinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT(:, 0, n, 1) = 0.d0
             
                 ! --------------
          END DO ! End layer loop
                 ! --------------
       END IF
       
       !  Debug optical properties
       
       if ( do_debug_geocape_tool ) then
          write(du,'(/a,i4,a, i4)')'Debug optical properties for wavelength # ',w, ' ic = ', ic
          if ( do_Jacobians ) then
             total_tau = 0.0d0
             do n = 1, VLIDORT_FixIn%Cont%TS_NLAYERS
                total_tau = total_tau + VLIDORT_FixIn%Optical%TS_DELTAU_VERT_INPUT(n)
                write(du,'(i3,1p12e20.10)') n, VLIDORT_FixIn%Optical%TS_DELTAU_VERT_INPUT(n), &
                     VLIDORT_ModIn%MOptical%TS_OMEGA_TOTAL_INPUT(n),                          &
                     VLIDORT_FixIn%Optical%TS_GREEKMAT_TOTAL_INPUT(1,n,1),                    &
                     VLIDORT_FixIn%Optical%TS_GREEKMAT_TOTAL_INPUT(2,n,1),                    &
                     VLIDORT_LinFixIn%Optical%TS_L_DELTAU_VERT_INPUT(1,n),                    &
                     VLIDORT_LinFixIn%Optical%TS_L_OMEGA_TOTAL_INPUT(1,n),                    & 
                     VLIDORT_LinFixIn%Optical%TS_L_DELTAU_VERT_INPUT(2,n),                    &
                     VLIDORT_LinFixIn%Optical%TS_L_OMEGA_TOTAL_INPUT(2,n),                    &
                     VLIDORT_LinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT(1,1,n,1),             &
                     VLIDORT_LinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT(2,1,n,1),             &
                     VLIDORT_LinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT(1,2,n,1),             &
                     VLIDORT_LinFixIn%Optical%TS_L_GREEKMAT_TOTAL_INPUT(2,2,n,1)
             enddo
             write(du,'(a,f10.5)')'--------Total optical depth = ',total_tau
          else
             total_tau = 0.0d0
             do n = 1, VLIDORT_FixIn%Cont%TS_NLAYERS
                total_tau = total_tau + VLIDORT_FixIn%Optical%TS_DELTAU_VERT_INPUT(n)
                write(du,'(i3,1p4e20.10)')n,VLIDORT_FixIn%Optical%TS_DELTAU_VERT_INPUT(n), &
                     VLIDORT_ModIn%MOptical%TS_OMEGA_TOTAL_INPUT(n),                       &
                     VLIDORT_FixIn%Optical%TS_GREEKMAT_TOTAL_INPUT(1,n,1),                 &
                     VLIDORT_FixIn%Optical%TS_GREEKMAT_TOTAL_INPUT(2,n,1)
             enddo
             write(du,'(a,f10.5)')'--------Total optical depth = ',total_tau
          endif
       endif
       
       ! ------------
       ! VLIDORT call
       ! ------------
       ! ----------------
       ! Now call VLIDORT
       ! -----------------
       WRITE(*,*)'doing VLIDORT calculation # ', w, didx
       CALL VLIDORT_LPS_MASTER ( &
            VLIDORT_FixIn,     &
            VLIDORT_ModIn,     &
            VLIDORT_Sup,       &
            VLIDORT_Out,       &
            VLIDORT_LinFixIn,  &
            VLIDORT_LinModIn,  &
            VLIDORT_LinSup,    &
            VLIDORT_LinOut )

       ! ------------------
       ! Exception handling
       ! ------------------
       OPENERRORFILEFLAG = .FALSE.
       CALL VLIDORT_WRITE_STATUS (            &
            'V2p6_VLIDORT_LPS_Execution.log', &
            35, OPENERRORFILEFLAG,            &
            VLIDORT_Out%Status )
       
       ! --------------
       ! Stop if failed
       ! --------------
       
       IF ( VLIDORT_Out%Status%TS_STATUS_INPUTCHECK .EQ. VLIDORT_SERIOUS ) THEN
          WRITE(*,*)'VLIDORT input abort, thread number = ', w
          STOP
       ELSE IF ( VLIDORT_Out%Status%TS_STATUS_INPUTCHECK .NE. VLIDORT_SERIOUS .AND. &
            VLIDORT_Out%Status%TS_STATUS_CALCULATION .EQ. VLIDORT_SERIOUS ) THEN
          WRITE(*,*)'VLIDORT calculation abort, thread number = ', w
          STOP
       END IF

       ! --------------
       ! Copy radiances
       ! --------------
       stokes_clrcld(1:VLIDORT_Out%Main%TS_N_GEOMETRIES, 1:VLIDORT_FixIn%Cont%TS_NSTOKES, ic) = &
            VLIDORT_Out%Main%TS_STOKES(1, 1:VLIDORT_Out%Main%TS_N_GEOMETRIES, &
            1:VLIDORT_FixIn%Cont%TS_NSTOKES, didx)

       ! ---------
       ! Copy flux
       ! ---------
       stokes_flux(1:VLIDORT_ModIn%MSunrays%TS_N_SZANGLES, 1:VLIDORT_FixIn%Cont%TS_NSTOKES, ic) = &
            VLIDORT_Out%Main%TS_FLUX_STOKES(1, 1:VLIDORT_ModIn%MSunrays%TS_N_SZANGLES, &
            1:VLIDORT_FixIn%Cont%TS_NSTOKES, didx)
       
       ! -----------------------------------
       ! Copy direct flux (only downwelling)
       ! -----------------------------------
       stokes_direct_flux(1:VLIDORT_ModIn%MSunrays%TS_N_SZANGLES, 1:VLIDORT_FixIn%Cont%TS_NSTOKES, ic) = &
            VLIDORT_Out%Main%TS_FLUX_DIRECT(1, 1:VLIDORT_ModIn%MSunrays%TS_N_SZANGLES, &
            1:VLIDORT_FixIn%Cont%TS_NSTOKES)
       
       IF (do_jacobians) THEN
          
          profilewf_sum(1:VLIDORT_LinFixIn%Cont%TS_N_TOTALPROFILE_WFS,               &
               1:VLIDORT_FixIn%Cont%TS_NLAYERS,                             &
               1:VLIDORT_Out%Main%TS_N_GEOMETRIES,                          &
               1:VLIDORT_FixIn%Cont%TS_NSTOKES) =                           &
               profilewf_sum(1:VLIDORT_LinFixIn%Cont%TS_N_TOTALPROFILE_WFS, &
               1:VLIDORT_FixIn%Cont%TS_NLAYERS,                             &
               1:VLIDORT_Out%Main%TS_N_GEOMETRIES,                          &
               1:VLIDORT_FixIn%Cont%TS_NSTOKES) +                           &
               VLIDORT_LinOut%Prof%TS_PROFILEWF                             &
               (1:VLIDORT_LinFixIn%Cont%TS_N_TOTALPROFILE_WFS,              &
               1:VLIDORT_FixIn%Cont%TS_NLAYERS, 1,                          &
               1:VLIDORT_Out%Main%TS_N_GEOMETRIES,                          &
               1:VLIDORT_FixIn%Cont%TS_NSTOKES, didx)*ipafrac
          surfacewf_clrcld(1:VLIDORT_LinFixIn%Cont%TS_N_SURFACE_WFS,                 &
               1:VLIDORT_Out%Main%TS_N_GEOMETRIES,                          &
               1:VLIDORT_FixIn%Cont%TS_NSTOKES, ic) =                       &
               VLIDORT_LinOut%Surf%TS_SURFACEWF                             &
               (1:VLIDORT_LinFixIn%Cont%TS_N_SURFACE_WFS, 1,                &
               1:VLIDORT_Out%Main%TS_N_GEOMETRIES,                          &
               1:VLIDORT_FixIn%Cont%TS_NSTOKES, didx)
       END IF
       
           ! ----------------------- 
    END DO ! end cloud fraction loop
           ! -----------------------
    
    VLIDORT_FixIn%Cont%TS_NLAYERS = GC_nlayers  ! use the most # of layers
    
    ! ------------------------------------------------------------------
    ! Independent pixel approximation for radiance, flux and direct flux
    ! ------------------------------------------------------------------
    VLIDORT_Out%Main%TS_STOKES(1, 1:VLIDORT_Out%Main%TS_N_GEOMETRIES, &
         1:VLIDORT_FixIn%Cont%TS_NSTOKES, didx) =                     &
         stokes_clrcld(1:VLIDORT_Out%Main%TS_N_GEOMETRIES,            &
         1:VLIDORT_FixIn%Cont%TS_NSTOKES, 1) * (1.0 - cfrac )         &
         + stokes_clrcld(1:VLIDORT_Out%Main%TS_N_GEOMETRIES,          &
         1:VLIDORT_FixIn%Cont%TS_NSTOKES, 2) * cfrac

    VLIDORT_Out%Main%TS_FLUX_STOKES(1, 1:VLIDORT_ModIn%MSunrays%TS_N_SZANGLES,   &
         1:VLIDORT_FixIn%Cont%TS_NSTOKES, didx) =                          &
         stokes_flux(1:VLIDORT_ModIn%MSunrays%TS_N_SZANGLES,                     &
         1:VLIDORT_FixIn%Cont%TS_NSTOKES, 1) * (1.0 - cfrac )              &
         + stokes_flux(1:VLIDORT_ModIn%MSunrays%TS_N_SZANGLES,                   &
         1:VLIDORT_FixIn%Cont%TS_NSTOKES, 2) * cfrac

    VLIDORT_Out%Main%TS_FLUX_DIRECT(1, 1:VLIDORT_ModIn%MSunrays%TS_N_SZANGLES,   &
         1:VLIDORT_FixIn%Cont%TS_NSTOKES) =                                &
         stokes_direct_flux(1:VLIDORT_ModIn%MSunrays%TS_N_SZANGLES,              &
         1:VLIDORT_FixIn%Cont%TS_NSTOKES, 1) * (1.0 - cfrac )              &
         + stokes_direct_flux(1:VLIDORT_ModIn%MSunrays%TS_N_SZANGLES,            &
         1:VLIDORT_FixIn%Cont%TS_NSTOKES, 2) * cfrac
    
    IF (do_jacobians) THEN
       VLIDORT_LinOut%Prof%TS_PROFILEWF(1:VLIDORT_LinFixIn%Cont%TS_N_TOTALPROFILE_WFS, &
            1:GC_nlayers, 1, 1:VLIDORT_Out%Main%TS_N_GEOMETRIES,                       &
            1:VLIDORT_FixIn%Cont%TS_NSTOKES, didx) =                                   &
            profilewf_sum(1:VLIDORT_LinFixIn%Cont%TS_N_TOTALPROFILE_WFS, 1:GC_nlayers, &
            1:VLIDORT_Out%Main%TS_N_GEOMETRIES,                                        &
            1:VLIDORT_FixIn%Cont%TS_NSTOKES)  
       
       IF (cfrac <= 0.0d0) THEN
          VLIDORT_LinOut%Surf%TS_SURFACEWF(1:VLIDORT_LinFixIn%Cont%TS_N_SURFACE_WFS, 1,     &
               1:VLIDORT_Out%Main%TS_N_GEOMETRIES, 1:VLIDORT_FixIn%Cont%TS_NSTOKES, didx) = &
               surfacewf_clrcld(1:VLIDORT_LinFixIn%Cont%TS_N_SURFACE_WFS,                   &
               1:VLIDORT_Out%Main%TS_N_GEOMETRIES, 1:VLIDORT_FixIn%Cont%TS_NSTOKES, 1)  
       ELSE IF (cfrac >= 1.0d0) THEN
          VLIDORT_LinOut%Surf%TS_SURFACEWF(1:VLIDORT_LinFixIn%Cont%TS_N_SURFACE_WFS, 1, &
               1:VLIDORT_Out%Main%TS_N_GEOMETRIES,                                      &
               1:VLIDORT_FixIn%Cont%TS_NSTOKES, didx) = &
               SURFACEWF_CLRCLD(1:VLIDORT_LinFixIn%Cont%TS_N_SURFACE_WFS,               &
               1:VLIDORT_Out%Main%TS_N_GEOMETRIES, 1:VLIDORT_FixIn%Cont%TS_NSTOKES, 2) 
       ELSE IF (do_lambertian_cld .AND. VLIDORT_FixIn%Bool%TS_DO_LAMBERTIAN_SURFACE) THEN ! from clear-sky only
          VLIDORT_LinOut%Surf%TS_SURFACEWF(1:VLIDORT_LinFixIn%Cont%TS_N_SURFACE_WFS, 1, &
               1:VLIDORT_Out%Main%TS_N_GEOMETRIES, 1:VLIDORT_FixIn%Cont%TS_NSTOKES,     &
               didx) = surfacewf_clrcld(1:VLIDORT_LinFixIn%Cont%TS_N_SURFACE_WFS,       &
               1:VLIDORT_Out%Main%TS_N_GEOMETRIES, 1:VLIDORT_FixIn%Cont%TS_NSTOKES, 1) *&
               ( 1.0-cfrac)
       ELSE IF (do_lambertian_cld .AND. .NOT. VLIDORT_FixIn%Bool%TS_DO_LAMBERTIAN_SURFACE) THEN ! Should not happen
          PRINT *, 'Should not happen!!!'; STOP
       ELSE
          VLIDORT_LinOut%Surf%TS_SURFACEWF(1:VLIDORT_LinFixIn%Cont%TS_N_SURFACE_WFS, 1,     &
               1:VLIDORT_Out%Main%TS_N_GEOMETRIES, 1:VLIDORT_FixIn%Cont%TS_NSTOKES, didx) = &
               surfacewf_clrcld(1:VLIDORT_LinFixIn%Cont%TS_N_SURFACE_WFS,                   &
               1:VLIDORT_Out%Main%TS_N_GEOMETRIES, 1:VLIDORT_FixIn%Cont%TS_NSTOKES, 1) *    &
               ( 1.0-cfrac) + surfacewf_clrcld(1:VLIDORT_LinFixIn%Cont%TS_N_SURFACE_WFS,    &
               1:VLIDORT_Out%Main%TS_N_GEOMETRIES, 1:VLIDORT_FixIn%Cont%TS_NSTOKES, 2) * cfrac
       END IF
       
       ! Compute cloud fraction weighting function if selected
       IF (do_cfrac_jacobians) THEN
          IF (cfrac .GT. 0.0d0 .AND. cfrac .LT. 1.0d0 ) THEN
             GC_cfrac_jacobians(w, 1:VLIDORT_Out%Main%TS_N_GEOMETRIES) =   &
                  stokes_clrcld(1:VLIDORT_Out%Main%TS_N_GEOMETRIES, 1, 2)  &
                  - stokes_clrcld(1:VLIDORT_Out%Main%TS_N_GEOMETRIES, 1, 1)
          ELSE
             GC_cfrac_Jacobians(w, 1:VLIDORT_Out%Main%TS_N_GEOMETRIES) = 0.0d0
          END IF
       END IF
       
       IF (do_cfrac_jacobians .AND. do_QU_Jacobians) THEN
          IF (cfrac .GT. 0.0d0 .AND. cfrac .LT. 1.0d0 ) THEN
             GC_cfrac_QJacobians(w, 1:VLIDORT_Out%Main%TS_N_GEOMETRIES) = &
                  stokes_clrcld(1:VLIDORT_Out%Main%TS_N_GEOMETRIES, 2, 2) &
                  - stokes_clrcld(1:VLIDORT_Out%Main%TS_N_GEOMETRIES, 2, 1)
             GC_cfrac_UJacobians(w, 1:VLIDORT_Out%Main%TS_N_GEOMETRIES) = &
                  stokes_clrcld(1:VLIDORT_Out%Main%TS_N_GEOMETRIES, 3, 2) &
                  - stokes_clrcld(1:VLIDORT_Out%Main%TS_N_GEOMETRIES, 3, 1)
          ELSE
             GC_cfrac_QJacobians(w, 1:VLIDORT_Out%Main%TS_N_GEOMETRIES) = 0.0d0
             GC_cfrac_UJacobians(w, 1:VLIDORT_Out%Main%TS_N_GEOMETRIES) = 0.0d0
          END IF
       END IF
    END IF
    
  END SUBROUTINE Vlidort_cloud_and_calculation

  SUBROUTINE VBRDF_TO_VLIDORT(error)

    IMPLICIT NONE

    ! ------------------
    ! Modified variables
    ! ------------------
    LOGICAL,       INTENT(INOUT) :: error ! Error variable

    ! ----------------
    ! Code starts here
    ! ----------------
    error = .FALSE.
    
    VLIDORT_Sup%BRDF%TS_BRDF_F_0        = VBRDF_Sup_Out%BS_BRDF_F_0
    VLIDORT_Sup%BRDF%TS_BRDF_F          = VBRDF_Sup_Out%BS_BRDF_F
    VLIDORT_Sup%BRDF%TS_USER_BRDF_F_0   = VBRDF_Sup_Out%BS_USER_BRDF_F_0
    VLIDORT_Sup%BRDF%TS_USER_BRDF_F     = VBRDF_Sup_Out%BS_USER_BRDF_F
    VLIDORT_Sup%BRDF%TS_EXACTDB_BRDFUNC = VBRDF_Sup_Out%BS_EXACTDB_BRDFUNC
    VLIDORT_Sup%BRDF%TS_EMISSIVITY      = VBRDF_Sup_Out%BS_EMISSIVITY
    VLIDORT_Sup%BRDF%TS_USER_EMISSIVITY = VBRDF_Sup_Out%BS_USER_EMISSIVITY
    
    VLIDORT_LinSup%BRDF%TS_LS_BRDF_F_0        = VBRDF_LinSup_Out%BS_LS_BRDF_F_0
    VLIDORT_LinSup%BRDF%TS_LS_BRDF_F          = VBRDF_LinSup_Out%BS_LS_BRDF_F
    VLIDORT_LinSup%BRDF%TS_LS_USER_BRDF_F_0   = VBRDF_LinSup_Out%BS_LS_USER_BRDF_F_0
    VLIDORT_LinSup%BRDF%TS_LS_USER_BRDF_F     = VBRDF_LinSup_Out%BS_LS_USER_BRDF_F
    VLIDORT_LinSup%BRDF%TS_LS_EXACTDB_BRDFUNC = VBRDF_LinSup_Out%BS_LS_EXACTDB_BRDFUNC
    VLIDORT_LinSup%BRDF%TS_LS_USER_EMISSIVITY = VBRDF_LinSup_Out%BS_LS_USER_EMISSIVITY
    VLIDORT_LinSup%BRDF%TS_LS_EMISSIVITY      = VBRDF_LinSup_Out%BS_LS_EMISSIVITY

  END SUBROUTINE VBRDF_TO_VLIDORT
  
END MODULE GC_Vlidort_module
