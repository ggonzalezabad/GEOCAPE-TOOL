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
                                  OUTPUT_WSABSA, WSA_CALCULATED, BSA_CALCULATED, didx, W
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
    call netcdf_wrt ( netfname, GC_n_sun_positions, VLIDORT_Out%Main%TS_N_GEOMETRIES, &
         GC_nlayers, ngases, opdeps(W, 1:GC_nlayers), ssalbs(W, 1:GC_nlayers),        &
         aer_opdeps(W, 1:GC_nlayers), aer_ssalbs(W, 1:GC_nlayers),                    &
         cld_opdeps(W, 1:GC_nlayers), cld_ssalbs(W, 1:GC_nlayers),                    &
         ground_ler(W), do_vector_calculation, do_StokesQU_output, do_Jacobians,      &
         do_QU_Jacobians, do_AMF_calculation, do_T_Jacobians, do_sfcprs_Jacobians,    &
         do_aod_Jacobians, do_assa_Jacobians, do_cod_Jacobians, do_cssa_Jacobians,    &
         do_cfrac_Jacobians, do_aer_columnwf, do_cld_columnwf, use_lambertian,        &
         solar_cspec_data(W), GC_radiances(W, 1:VLIDORT_Out%Main%TS_N_GEOMETRIES),    &
         GC_Qvalues(W, 1:VLIDORT_Out%Main%TS_N_GEOMETRIES),                           &
         GC_Uvalues(W, 1:VLIDORT_Out%Main%TS_N_GEOMETRIES),                           &
         GC_Tracegas_Jacobians(W, 1:GC_nlayers, 1:VLIDORT_Out%Main%TS_N_GEOMETRIES, 1:ngases),  &
         GC_Scattering_Weights(W, 1:GC_nlayers, 1:VLIDORT_Out%Main%TS_N_GEOMETRIES),            &
         GC_AMFs(W, 1:VLIDORT_Out%Main%TS_N_GEOMETRIES, 1:ngases),                              &
         GC_Tracegas_QJacobians(W, 1:GC_nlayers, 1:VLIDORT_Out%Main%TS_N_GEOMETRIES, 1:ngases), &
         GC_Tracegas_UJacobians(W, 1:GC_nlayers, 1:VLIDORT_Out%Main%TS_N_GEOMETRIES, 1:ngases), &
         GC_Temperature_Jacobians(W, 1:GC_nlayers, 1:VLIDORT_Out%Main%TS_N_GEOMETRIES),         &
         GC_Surfalbedo_Jacobians(W, 1:VLIDORT_Out%Main%TS_N_GEOMETRIES),                        &
         GC_Surfalbedo_QJacobians(W, 1:VLIDORT_Out%Main%TS_N_GEOMETRIES),                       &
         GC_Surfalbedo_UJacobians(W, 1:VLIDORT_Out%Main%TS_N_GEOMETRIES),                       &
         GC_Windspeed_Jacobians(W, 1:VLIDORT_Out%Main%TS_N_GEOMETRIES),                         &
         GC_Windspeed_QJacobians(W, 1:VLIDORT_Out%Main%TS_N_GEOMETRIES),                        &
         GC_Windspeed_UJacobians(W, 1:VLIDORT_Out%Main%TS_N_GEOMETRIES),                        &
         GC_aod_Jacobians(W, 1:GC_nlayers, 1:VLIDORT_Out%Main%TS_N_GEOMETRIES),                 &
         GC_aod_QJacobians(W, 1:GC_nlayers, 1:VLIDORT_Out%Main%TS_N_GEOMETRIES),                &
         GC_aod_UJacobians(W, 1:GC_nlayers, 1:VLIDORT_Out%Main%TS_N_GEOMETRIES),                &
         GC_assa_Jacobians(W, 1:GC_nlayers, 1:VLIDORT_Out%Main%TS_N_GEOMETRIES),                &
         GC_assa_QJacobians(W, 1:GC_nlayers, 1:VLIDORT_Out%Main%TS_N_GEOMETRIES),               &
         GC_assa_UJacobians(W, 1:GC_nlayers, 1:VLIDORT_Out%Main%TS_N_GEOMETRIES),               &
         GC_cod_Jacobians(W, 1:GC_nlayers, 1:VLIDORT_Out%Main%TS_N_GEOMETRIES),                 &
         GC_cod_QJacobians(W, 1:GC_nlayers, 1:VLIDORT_Out%Main%TS_N_GEOMETRIES),                &
         GC_cod_UJacobians(W, 1:GC_nlayers, 1:VLIDORT_Out%Main%TS_N_GEOMETRIES),                &
         GC_cssa_Jacobians(W, 1:GC_nlayers, 1:VLIDORT_Out%Main%TS_N_GEOMETRIES),                &
         GC_cssa_QJacobians(W, 1:GC_nlayers, 1:VLIDORT_Out%Main%TS_N_GEOMETRIES),               &
         GC_cssa_UJacobians(W, 1:GC_nlayers, 1:VLIDORT_Out%Main%TS_N_GEOMETRIES),               &
         GC_cfrac_Jacobians(W, 1:VLIDORT_Out%Main%TS_N_GEOMETRIES),                             &
         GC_cfrac_QJacobians(W, 1:VLIDORT_Out%Main%TS_N_GEOMETRIES),                            &
         GC_cfrac_UJacobians(W, 1:VLIDORT_Out%Main%TS_N_GEOMETRIES),                            &
         GC_sfcprs_Jacobians(W, 1:VLIDORT_Out%Main%TS_N_GEOMETRIES),                            &
         GC_sfcprs_QJacobians(W, 1:VLIDORT_Out%Main%TS_N_GEOMETRIES),                           &
         GC_sfcprs_UJacobians(W, 1:VLIDORT_Out%Main%TS_N_GEOMETRIES),                           &
         GC_flux(W, 1:VLIDORT_ModIn%MSunrays%TS_N_SZANGLES),                                    &
         GC_Qflux(W, 1:VLIDORT_ModIn%MSunrays%TS_N_SZANGLES),                                   &
         GC_Uflux(W, 1:VLIDORT_ModIn%MSunrays%TS_N_SZANGLES),                                   &
         GC_direct_flux(W, 1:VLIDORT_ModIn%MSunrays%TS_N_SZANGLES),                             &
         GC_Qdirect_flux(W, 1:VLIDORT_ModIn%MSunrays%TS_N_SZANGLES),                            &
         GC_Udirect_flux(W, 1:VLIDORT_ModIn%MSunrays%TS_N_SZANGLES),                            &
         W, tmperror, error)

    IF (error) CALL write_err_message (.TRUE., tmperror)

  END SUBROUTINE netcdf_output

  SUBROUTINE Create_netcdf_output_file (error)

    IMPLICIT NONE
    INCLUDE '../netcdf/netcdf.inc'
    ! ------------------
    ! Modified variables
    ! ------------------
    LOGICAL,       INTENT(INOUT) :: error ! Error variable

    ! ---------------
    ! Local variables
    ! ---------------
    CHARACTER(LEN=max_ch_len) :: tmperror
    CHARACTER(LEN=max_ch_len) :: netfname

    integer :: ncid, rcode, nlen, i
    integer :: szadim, vzadim, azadim, gasdim, lvldim, laydim, wavdim, geodim, nsqdim, onedim
    integer :: szaid, vzaid, azaid, gasid, lvlid, wavid, psid, tsid, airid, &
         aer0id, cld0id, radid, qid, uid, irradid, sfcid, sfcwfid, wswfid, cfracwfid, &
         sfcprswfid, aodwfid, assawfid, codwfid, cssawfid, sfcqwfid, wsqwfid, cfracqwfid, &
         sfcprsqwfid, aodqwfid, assaqwfid, codqwfid, cssaqwfid, sfcuwfid, wsuwfid, cfracuwfid,&
         sfcprsuwfid, aoduwfid, assauwfid, coduwfid, cssauwfid, gascolid, aodsid, assasid, &
         codsid, cssasid, scatwtid, odsid, ssasid, amfid, gaswfid, gasqwfid, &
         gasuwfid, tempwfid, fluxid, dfluxid, qfluxid, ufluxid, qdfluxid, udfluxid, brdfid, &
         wsaid, bsaid
    integer, dimension(2)      :: gascol_dims, wavalt_dims, wavgas_dims, wavgeo_dims, wavsza_dims
    integer, dimension(3)      :: wavaltgeo_dims, wavgeogas_dims
    integer, dimension(4)      :: wavaltgeogas_dims, brdfdim
    integer, dimension(ngases) :: gas_indices
    integer, dimension(1)      :: ndimstart1, ndimcount1
    integer, dimension(2)      :: ndimstart2, ndimcount2
    integer, dimension(4)      :: ndimstart4, ndimcount4
    character(len=256)         :: gasstr
    character(len=30)          :: irrad_unitc, rad_unitc
    integer(kind=4)            :: ngeom

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

    error = .false.; tmperror = ' '
    
    ncid = nccre(trim(netfname), ncclob, rcode)
    if (rcode .eq. -1) then
       tmperror =  ' error in netcdf_out: nccre failed'
       error = .true.; return
    endif

    ngeom = GC_n_sun_positions * GC_n_view_angles * GC_n_azimuths
    ! Create the dimensions of the dataset:
    onedim  = ncddef (ncid, 'one', 1,  rcode)
    szadim  = ncddef (ncid, 'nsza', GC_n_sun_positions, rcode)
    vzadim  = ncddef (ncid, 'nvza', GC_n_view_angles, rcode)
    azadim  = ncddef (ncid, 'naza', GC_n_azimuths, rcode)
    gasdim  = ncddef (ncid, 'ngas', ngases, rcode)
    lvldim  = ncddef (ncid, 'nlevel', GC_nlayers+1, rcode)
    laydim  = ncddef (ncid, 'nlayer', GC_nlayers, rcode)
    wavdim  = ncddef (ncid, 'nw', nlambdas, rcode)
    geodim  = ncddef (ncid, 'ngeom',  ngeom, rcode)
    nsqdim  = ncddef (ncid, 'nstokessq', NSTOKESSQ, rcode)
    
    !=============================================================================
    ! Create the coordinate (aka independent) variables:
    !============================================================================= 
    szaid  = ncvdef (ncid, 'Solarzenithangle',     ncfloat, 1, szadim, rcode)
    vzaid  = ncvdef (ncid, 'Viewingzenithangle',   ncfloat, 1, vzadim, rcode)
    azaid  = ncvdef (ncid, 'Relativeazimuthangle', ncfloat, 1, azadim, rcode)
    gasid  = ncvdef (ncid, 'gas',                  nclong,  1, gasdim, rcode)
    lvlid  = ncvdef (ncid, 'zs',                   ncfloat, 1, lvldim, rcode)
    wavid  = ncvdef (ncid, 'Wavelength',           ncfloat, 1, wavdim, rcode)
    
    !=============================================================================
    ! Create a vector containing the often referred to spacetime coordinates:
    !=============================================================================
    gascol_dims(1) = laydim; gascol_dims(2) = gasdim
    wavalt_dims(1) = wavdim; wavalt_dims(2) = laydim
    wavgas_dims(1) = wavdim; wavgas_dims(2) = gasdim
    wavgeo_dims(1) = wavdim; wavgeo_dims(2) = geodim
    wavsza_dims(1) = wavdim; wavsza_dims(2) = szadim
    
    wavaltgeo_dims(1) = wavdim; wavaltgeo_dims(2) = laydim; wavaltgeo_dims(3) = geodim
    wavgeogas_dims(1) = wavdim; wavgeogas_dims(2) = geodim; wavgeogas_dims(3) = gasdim
    
    wavaltgeogas_dims(1) = wavdim; wavaltgeogas_dims(2) = laydim
    wavaltgeogas_dims(3) = geodim; wavaltgeogas_dims(4) = gasdim
    
    brdfdim(1) = nsqdim; brdfdim(2) = vzadim; brdfdim(3)=azadim; brdfdim(4) = szadim

    !=============================================================================
    ! Create the dependent variables:
    !=============================================================================
    ! ps, ts
    psid = ncvdef(ncid, 'ps', ncfloat, 1, lvldim, rcode)
    tsid = ncvdef(ncid, 'ts', ncfloat, 1, lvldim, rcode)
    
    ! aircol, aods0, cods0
    airid  = ncvdef(ncid, 'aircol', ncfloat, 1, laydim, rcode)
    aer0id = ncvdef(ncid, 'aods0',  ncfloat, 1, laydim, rcode)
    cld0id = ncvdef(ncid, 'cods0',  ncfloat, 1, laydim, rcode)
        
    ! varaibles with 1D, wavdim
    irradid = ncvdef(ncid, 'irradiance', ncfloat, 1, wavdim, rcode)
    sfcid   = ncvdef(ncid, 'surfalb',    ncfloat, 1, wavdim, rcode)
    
    ! variables with 2D, laydim, gasdim
    gascolid = ncvdef(ncid,  'gascol', ncdouble, 2, gascol_dims, rcode)
    
    ! variables with 2D, wavdim, laydim
    odsid   = ncvdef(ncid, 'ods',   ncfloat, 2, wavalt_dims, rcode)
    ssasid  = ncvdef(ncid, 'ssas',  ncfloat, 2, wavalt_dims, rcode)
    aodsid  = ncvdef(ncid, 'aods',  ncfloat, 2, wavalt_dims, rcode)
    assasid = ncvdef(ncid, 'assas', ncfloat, 2, wavalt_dims, rcode)
    codsid  = ncvdef(ncid, 'cods',  ncfloat, 2, wavalt_dims, rcode)
    cssasid = ncvdef(ncid, 'cssas', ncfloat, 2, wavalt_dims, rcode)
    
    ! variables with 2D, wavdim, geodim
    radid   = ncvdef(ncid, 'radiance',    ncfloat, 2, wavgeo_dims, rcode)
    fluxid  = ncvdef(ncid, 'flux',        ncfloat, 2, wavsza_dims, rcode)
    dfluxid = ncvdef(ncid, 'direct_flux', ncfloat, 2, wavsza_dims, rcode)
    if (do_vector_calculation .and. do_StokesQU_output) then
       qid      = ncvdef(ncid, 'q',            ncfloat, 2, wavgeo_dims, rcode)
       uid      = ncvdef(ncid, 'u',            ncfloat, 2, wavgeo_dims, rcode)
       qfluxid  = ncvdef(ncid, 'qflux',        ncfloat, 2, wavsza_dims, rcode)
       ufluxid  = ncvdef(ncid, 'uflux',        ncfloat, 2, wavsza_dims, rcode)
       qdfluxid = ncvdef(ncid, 'qdirect_flux', ncfloat, 2, wavsza_dims, rcode)
       udfluxid = ncvdef(ncid, 'udirect_flux', ncfloat, 2, wavsza_dims, rcode)
    endif
    
    if (do_Jacobians) then
       
       if (use_lambertian) then
          sfcwfid = ncvdef(ncid, 'surfalb_jac',   ncfloat, 2, wavgeo_dims, rcode)
       else 
          wswfid  = ncvdef(ncid, 'windspeed_jac', ncfloat, 2, wavgeo_dims, rcode)
       endif
       
       if (do_cfrac_Jacobians) &
            cfracwfid  = ncvdef(ncid, 'cfrac_jac',     ncfloat, 2, wavgeo_dims, rcode)
       if (do_sfcprs_Jacobians) &
            sfcprswfid = ncvdef(ncid, 'sfcprs_jac',    ncfloat, 2, wavgeo_dims, rcode)
       if (do_aer_columnwf .and. do_aod_Jacobians) &  !No column jacobians yet
            aodwfid    = ncvdef(ncid, 'aodcolwf_jac',  ncfloat, 2, wavgeo_dims, rcode)
       if (do_aer_columnwf .and. do_assa_Jacobians) & !No column jacobians yet
            assawfid   = ncvdef(ncid, 'assacolwf_jac', ncfloat, 2, wavgeo_dims, rcode)
       if (do_cld_columnwf .and. do_cod_Jacobians) &  !No column jacobians yet
            codwfid    = ncvdef(ncid, 'codcolwf_jac',  ncfloat, 2, wavgeo_dims, rcode)
       if (do_cld_columnwf .and. do_cssa_Jacobians) & !No column jacobians yet
            cssawfid   = ncvdef(ncid, 'cssacolwf_jac', ncfloat, 2, wavgeo_dims, rcode)
       
    endif
    
    if (do_QU_Jacobians) then
       
       if (use_lambertian) then
          sfcqwfid = ncvdef(ncid, 'surfalb_qjac',  ncfloat, 2, wavgeo_dims, rcode)
          sfcuwfid = ncvdef(ncid, 'surfalb_ujac',  ncfloat, 2, wavgeo_dims, rcode)
       else 
          wsqwfid = ncvdef(ncid, 'windspeed_qjac', ncfloat, 2, wavgeo_dims, rcode)
          wsuwfid = ncvdef(ncid, 'windspeed_ujac', ncfloat, 2, wavgeo_dims, rcode)
       endif
       
       if (do_cfrac_Jacobians) then
          cfracqwfid = ncvdef(ncid, 'cfrac_qjac',  ncfloat, 2, wavgeo_dims, rcode)
          cfracuwfid = ncvdef(ncid, 'cfrac_ujac',  ncfloat, 2, wavgeo_dims, rcode)
       endif
       if (do_sfcprs_Jacobians) then
          sfcprsqwfid = ncvdef(ncid, 'sfcprs_qjac',  ncfloat, 2, wavgeo_dims, rcode)
          sfcprsuwfid = ncvdef(ncid, 'sfcprs_ujac',  ncfloat, 2, wavgeo_dims, rcode)
       endif
       if (do_aer_columnwf .and. do_aod_Jacobians) then  !No column jacobians yet
          aodqwfid = ncvdef(ncid, 'aodcol_qjac',  ncfloat, 2, wavgeo_dims, rcode)
          aoduwfid = ncvdef(ncid, 'aodcol_ujac',  ncfloat, 2, wavgeo_dims, rcode)
       endif
       if (do_aer_columnwf .and. do_assa_Jacobians) then !No column jacobians yet
          assaqwfid = ncvdef(ncid, 'assacol_qjac', ncfloat, 2, wavgeo_dims, rcode)
          assauwfid = ncvdef(ncid, 'assacol_ujac', ncfloat, 2, wavgeo_dims, rcode)
       endif
       if (do_cld_columnwf .and. do_cod_Jacobians) then  !No column jacobians yet
          codqwfid = ncvdef(ncid, 'codcol_qjac',  ncfloat, 2, wavgeo_dims, rcode)
          coduwfid = ncvdef(ncid, 'codcol_ujac',  ncfloat, 2, wavgeo_dims, rcode)
       endif
       if (do_cld_columnwf .and. do_cssa_Jacobians) then !No column jacobians yet
          cssaqwfid = ncvdef(ncid, 'cssacol_qjac', ncfloat, 2, wavgeo_dims, rcode)
          cssauwfid = ncvdef(ncid, 'cssacol_ujac', ncfloat, 2, wavgeo_dims, rcode)
       endif
       
    endif
    
    ! variables with 3D, wavdim, laydim, geodim
    if (do_AMF_calculation) then
       scatwtid = ncvdef(ncid, 'scatweights', ncfloat, 3, wavaltgeo_dims, rcode)
    endif
    
    if (do_Jacobians) then
       if (do_T_Jacobians) &
            tempwfid = ncvdef(ncid, 't_jac',    ncfloat, 3, wavaltgeo_dims, rcode)
       if (.not. do_aer_columnwf .and. do_aod_Jacobians) &
            aodwfid  = ncvdef(ncid, 'aod_jac',  ncfloat, 3, wavaltgeo_dims, rcode)  
       if (.not. do_aer_columnwf .and. do_assa_Jacobians) &
            assawfid = ncvdef(ncid, 'assa_jac', ncfloat, 3, wavaltgeo_dims, rcode)   
       if (.not. do_cld_columnwf .and. do_cod_Jacobians) &
            codwfid  = ncvdef(ncid, 'cod_jac',  ncfloat, 3, wavaltgeo_dims, rcode)  
       if (.not. do_cld_columnwf .and. do_cssa_Jacobians) &
            cssawfid = ncvdef(ncid, 'cssa_jac', ncfloat, 3, wavaltgeo_dims, rcode)  
    endif
    
    if (do_QU_Jacobians) then
       if (.not. do_aer_columnwf .and. do_aod_Jacobians) then 
          aodqwfid = ncvdef(ncid, 'aod_qjac', ncfloat, 3, wavaltgeo_dims, rcode)  
          aoduwfid = ncvdef(ncid, 'aod_ujac', ncfloat, 3, wavaltgeo_dims, rcode)  
       endif
       if (.not. do_aer_columnwf .and. do_assa_Jacobians) then
          assaqwfid = ncvdef(ncid, 'assa_qjac', ncfloat, 3, wavaltgeo_dims, rcode)   
          assauwfid = ncvdef(ncid, 'assa_ujac', ncfloat, 3, wavaltgeo_dims, rcode)   
       endif
       if (.not. do_cld_columnwf .and. do_cod_Jacobians) then
          codqwfid = ncvdef(ncid, 'cod_qjac', ncfloat, 3, wavaltgeo_dims, rcode)  
          coduwfid = ncvdef(ncid, 'cod_ujac', ncfloat, 3, wavaltgeo_dims, rcode)  
       endif
       if (.not. do_cld_columnwf .and. do_cssa_Jacobians) then
          cssaqwfid = ncvdef(ncid, 'cssa_qjac', ncfloat, 3, wavaltgeo_dims, rcode)  
          cssauwfid = ncvdef(ncid, 'cssa_ujac', ncfloat, 3, wavaltgeo_dims, rcode) 
       endif
    endif
    
    ! variables with 3D, wavdim, geodim, gasdim
    if (do_AMF_calculation) then
       amfid = ncvdef(ncid, 'amf', ncfloat, 3, wavgeogas_dims, rcode)
    endif
    
    ! variables with 4D, wavdim, laydim, geodim, gasdim
    if (do_Jacobians) then
       gaswfid = ncvdef(ncid,  'gas_jac', ncfloat, 4, wavaltgeogas_dims, rcode)
       if (do_QU_Jacobians) then
          gasqwfid = ncvdef(ncid,  'gas_qjac', ncfloat, 4, wavaltgeogas_dims, rcode)
          gasuwfid = ncvdef(ncid,  'gas_ujac', ncfloat, 4, wavaltgeogas_dims, rcode)
       endif
    endif
    
    ! BRDF variables onedim, nsqdim, nvza, naza, nsza
    if (do_brdf_surface) then
    ! WSA, BSA amd BRDF
       wsaid = ncvdef(ncid, 'WSA', ncfloat, 1, onedim, rcode)
       bsaid = ncvdef(ncid, 'BSA', ncfloat, 1, onedim, rcode)
       brdfid = ncvdef(ncid, 'BRDF', ncfloat, 4, brdfdim, rcode)
    endif
    
    !=============================================================================
    ! Assign attributes (meta-data) to the various variables:
    !============================================================================= 
    call ncapt (ncid, ncglobal, 'lon',          ncfloat, 1, real(longitude, kind=4), rcode)
    call ncapt (ncid, ncglobal, 'lat',          ncfloat, 1, real(latitude, kind=4), rcode)
    call ncapt (ncid, ncglobal, 'month',        nclong,  1, month, rcode)
    call ncapt (ncid, ncglobal, 'year',         nclong,  1, year, rcode)
    call ncapt (ncid, ncglobal, 'day',          nclong,  1, day, rcode)
    call ncapt (ncid, ncglobal, 'utc',          ncfloat, 1, real(utc, kind=4), rcode)
    call ncapt (ncid, ncglobal, 'cfrac',        ncfloat, 1, real(cfrac, kind=4), rcode)
    call ncapt (ncid, ncglobal, 'zcldtop',      ncfloat, 1, real(cld_tops, kind=4), rcode)
    call ncapt (ncid, ncglobal, 'cldalb',       ncfloat, 1, real(lambertian_cldalb, kind=4), rcode)
    call ncapt (ncid, ncglobal, 'userlvl',      ncfloat, 1, real(VLIDORT_ModIn%MUserVal%TS_USER_LEVELS, kind=4), rcode)
    call ncapt (ncid, ncglobal, 'windspeed',    ncfloat, 1, real(wind_speed, kind=4), rcode)
    call ncapt (ncid, ncglobal, 'lam_reso',     ncfloat, 1, real(lambda_resolution, kind=4), rcode)
    call ncapt (ncid, ncglobal, 'aer_reflam',   ncfloat, 1, real(aer_reflambda, kind=4), rcode)
    call ncapt (ncid, ncglobal, 'cld_reflam',   ncfloat, 1, real(cld_reflambda, kind=4), rcode)
    call ncapt (ncid, ncglobal, 'lam_dw',       ncfloat, 1, real(lambda_dw, kind=4), rcode)
    call ncapt (ncid, ncglobal, 'lam_dfw',      ncfloat, 1, real(lambda_dfw, kind=4), rcode)
    call ncapt (ncid, ncglobal, 'nstreams',     nclong,  1, VLIDORT_FixIn%Cont%TS_NSTREAMS, rcode)
    call ncapt (ncid, ncglobal, 'aercld_nmoms', nclong,  1, aercld_nmoments_input, rcode)
    call ncapt (ncid, ncglobal, 'nsza',         nclong,  1, GC_n_sun_positions, rcode)
    call ncapt (ncid, ncglobal, 'nvza',         nclong,  1, GC_n_view_angles, rcode)
    call ncapt (ncid, ncglobal, 'naza',         nclong,  1, GC_n_azimuths, rcode)
    call ncapt (ncid, ncglobal, 'ngeometries',  nclong,  1, ngeom, rcode)
    call ncapt (ncid, ncglobal, 'nwavelengths', nclong,  1, nlambdas, rcode)
    call ncapt (ncid, ncglobal, 'nlayers',      nclong,  1, GC_nlayers, rcode)
    call ncapt (ncid, ncglobal, 'ngas',         nclong,  1, ngases, rcode)
    
    if (use_lambertian) then
       call ncapt (ncid, ncglobal, 'user_lambertian', ncbyte, 1, 1, rcode)
    else
       call ncapt (ncid, ncglobal, 'user_lambertian', ncbyte, 1, 0, rcode)
    endif
    if (do_lambertian_cld) then
       call ncapt (ncid, ncglobal, 'do_lambertian_cld', ncbyte, 1, 1, rcode)
    else
       call ncapt (ncid, ncglobal, 'do_lambertian_cld', ncbyte, 1, 0, rcode)
    endif
    if (do_effcrs) then
       call ncapt (ncid, ncglobal, 'do_effcrs', ncbyte, 1, 1, rcode)
    else
       call ncapt (ncid, ncglobal, 'do_effcrs', ncbyte, 1, 0, rcode)
    endif
    if (use_wavelength) then
       call ncapt (ncid, ncglobal, 'use_wavelength', ncbyte, 1, 1, rcode)
    else
       call ncapt (ncid, ncglobal, 'use_wavelength', ncbyte, 1, 0, rcode)
    endif
    if (VLIDORT_FixIn%Bool%TS_DO_UPWELLING) then
       call ncapt (ncid, ncglobal, 'do_upwelling', ncbyte, 1, 1, rcode)
    else
       call ncapt (ncid, ncglobal, 'do_upwelling', ncbyte, 1, 0, rcode)
    endif
    if (do_normalized_WFoutput) then
       call ncapt (ncid, ncglobal, 'do_norm_WFoutput', ncbyte, 1, 1, rcode)
    else
       call ncapt (ncid, ncglobal, 'do_norm_WFoutput', ncbyte, 1, 0, rcode)
    endif
    if (do_normalized_radiance) then
       call ncapt (ncid, ncglobal, 'do_norm_radiance', ncbyte, 1, 1, rcode)
    else
       call ncapt (ncid, ncglobal, 'do_norm_radiance', ncbyte, 1, 0, rcode)
    endif
    if (use_solar_photons) then
       call ncapt (ncid, ncglobal, 'use_solar_photons', ncbyte, 1, 1, rcode)
    else
       call ncapt (ncid, ncglobal, 'use_solar_photons', ncbyte, 1, 0, rcode)
    endif
    if (do_vector_calculation) then
       call ncapt (ncid, ncglobal, 'do_vector_calc', ncbyte, 1, 1, rcode)
    else
       call ncapt (ncid, ncglobal, 'do_vector_calc', ncbyte, 1, 0, rcode)
    endif
    if (do_StokesQU_output) then
       call ncapt (ncid, ncglobal, 'do_QU_output', ncbyte, 1, 1, rcode)
    else
       call ncapt (ncid, ncglobal, 'do_QU_output', ncbyte, 1, 0, rcode)
    endif
    if (do_Jacobians) then
       call ncapt (ncid, ncglobal, 'do_Jacobian', ncbyte, 1, 1, rcode)
    else
       call ncapt (ncid, ncglobal, 'do_Jacobian', ncbyte, 1, 0, rcode)
    endif
    if (do_QU_Jacobians) then
       call ncapt (ncid, ncglobal, 'do_QU_Jacobian', ncbyte, 1, 1, rcode)
    else
       call ncapt (ncid, ncglobal, 'do_QU_Jacobian', ncbyte, 1, 0, rcode)
    endif
    if (do_AMF_calculation) then
       call ncapt (ncid, ncglobal, 'do_AMF_calc', ncbyte, 1, 1, rcode)
    else
       call ncapt (ncid, ncglobal, 'do_AMF_calc', ncbyte, 1, 0, rcode)
    endif
    if (do_T_Jacobians) then
       call ncapt (ncid, ncglobal, 'do_T_Jacobian', ncbyte, 1, 1, rcode)
    else
       call ncapt (ncid, ncglobal, 'do_T_Jacobian', ncbyte, 1, 0, rcode)
    endif
    if (do_sfcprs_Jacobians) then
       call ncapt (ncid, ncglobal, 'do_sfcprs_Jacobian', ncbyte, 1, 1, rcode)
    else
       call ncapt (ncid, ncglobal, 'do_sfcprs_Jacobian', ncbyte, 1, 0, rcode)
    endif
    if (do_aod_Jacobians) then
       call ncapt (ncid, ncglobal, 'do_aod_Jacobian', ncbyte, 1, 1, rcode)
    else
       call ncapt (ncid, ncglobal, 'do_aod_Jacobian', ncbyte, 1, 0, rcode)
    endif
    if (do_assa_Jacobians) then
       call ncapt (ncid, ncglobal, 'do_assa_Jacobian', ncbyte, 1, 1, rcode)
    else
       call ncapt (ncid, ncglobal, 'do_assa_Jacobian', ncbyte, 1, 0, rcode)
    endif
    if (do_cod_Jacobians) then
       call ncapt (ncid, ncglobal, 'do_cod_Jacobian', ncbyte, 1, 1, rcode)
    else
       call ncapt (ncid, ncglobal, 'do_cod_Jacobian', ncbyte, 1, 0, rcode)
    endif
    if (do_cssa_Jacobians) then
       call ncapt (ncid, ncglobal, 'do_cssa_Jacobian', ncbyte, 1, 1, rcode)
    else
       call ncapt (ncid, ncglobal, 'do_cssa_Jacobian', ncbyte, 1, 0, rcode)
    endif
    if (do_cfrac_Jacobians) then
       call ncapt (ncid, ncglobal, 'do_cfrac_Jacobian', ncbyte, 1, 1, rcode)
    else
       call ncapt (ncid, ncglobal, 'do_cfrac_Jacobian', ncbyte, 1, 0, rcode)
    endif
    if (do_aer_columnwf) then
       call ncapt (ncid, ncglobal, 'do_aer_columnwf', ncbyte, 1, 1, rcode)
    else
       call ncapt (ncid, ncglobal, 'do_aer_columnwf', ncbyte, 1, 0, rcode)
    endif
    if (do_cld_columnwf) then
       call ncapt (ncid, ncglobal, 'do_cld_columnwf', ncbyte, 1, 1, rcode)
    else
       call ncapt (ncid, ncglobal, 'do_cld_columnwf', ncbyte, 1, 0, rcode)
    endif
    if (do_brdf_surface) then
       call ncapt (ncid, ncglobal, 'do_brdf', ncbyte, 1, 1, rcode)
    else
       call ncapt (ncid, ncglobal, 'do_brdf', ncbyte, 1, 0, rcode)
    endif
    if (OUTPUT_WSABSA) then
       call ncapt (ncid, ncglobal, 'do_output_wsabsa', ncbyte, 1, 1, rcode)
    else
       call ncapt (ncid, ncglobal, 'do_output_wsabsa', ncbyte, 1, 0, rcode)
    endif
    
    ! write the list of gases in one string
    write(gasstr, '(10(I2,A1,A4,A1))') (i,':',which_gases(i),',', i=1, ngases)
    nlen=len_trim(gasstr); gasstr=gasstr(1:nlen-1); nlen=nlen-1
    call ncaptc (ncid, ncglobal, 'gases', ncchar,  nlen, trim(gasstr), rcode)
    
    ! Units for variables
    call ncaptc (ncid, ncglobal, 'windspeed_units', ncchar, 3, 'm/s', rcode)
    call ncaptc (ncid, ncglobal, 'lonlat_units',    ncchar, 7, 'degrees', rcode)
    call ncaptc (ncid, ncglobal, 'zcldtop_units',   ncchar, 2, 'km', rcode)
    call ncaptc (ncid, szaid,    'units', ncchar, 7, 'degrees', rcode)
    call ncaptc (ncid, vzaid,    'units', ncchar, 7, 'degrees', rcode)
    call ncaptc (ncid, azaid,    'units', ncchar, 7, 'degrees', rcode)
    call ncaptc (ncid, lvlid,    'units', ncchar, 2, 'km', rcode)
    if (use_wavelength) then
       call ncaptc(ncid,wavid,'units', ncchar,2, 'nm',rcode)
    else 
       call ncaptc(ncid,wavid,'units', ncchar,4, 'cm-1',rcode)
    endif
    if (use_solar_photons) then
       irrad_unitc='photons/cm2/nm/s'
    else
       irrad_unitc='W/m2/cm-1'
    endif
    call ncaptc(ncid, irradid, 'units', ncchar, len_trim(irrad_unitc), irrad_unitc, rcode)
    if (do_normalized_radiance) then
       rad_unitc='unitless'
    else
       rad_unitc=TRIM(irrad_unitc)//'/sr'
    endif
    call ncaptc(ncid, radid,    'units', ncchar, len_trim(rad_unitc), rad_unitc, rcode) 
    call ncaptc(ncid, gascolid, 'units', ncchar, 13, 'molecules/cm2',            rcode)
    call ncaptc(ncid, airid,    'units', ncchar, 13, 'molecules/cm2',            rcode)
    call ncaptc(ncid, psid,     'units', ncchar, 3,  'hPa',                      rcode)
    call ncaptc(ncid, tsid,     'units', ncchar, 1,  'K',                        rcode)    

    !=============================================================================
    ! Get out of 'define' mode, and into 'data' mode
    !=============================================================================
    call ncendf (ncid, rcode)

    !=============================================================================
    ! Fill in the values of the non time varying variables.
    ! Remember, we're still in "initialization" here - this is only done the
    ! first time through.
    !=============================================================================
    call ncvpt (ncid, szaid, 1, GC_n_sun_positions, &
         real(VLIDORT_ModIn%MSunrays%TS_SZANGLES(1:GC_n_sun_positions), kind=4), rcode)
    call ncvpt (ncid, vzaid, 1, GC_n_view_angles, &
         real(VLIDORT_ModIn%MUserVal%TS_USER_VZANGLES_INPUT(1:GC_n_view_angles), kind=4), rcode)
    call ncvpt (ncid, azaid, 1, GC_n_azimuths, &
         real(VLIDORT_ModIn%MUserVal%TS_USER_RELAZMS(1:GC_n_azimuths), kind=4), rcode)
    call ncvpt (ncid, lvlid, 1, GC_nlayers+1, &
         real(heights(0:GC_nlayers), kind=4), rcode)
    call ncvpt (ncid, wavid, 1, nlambdas, &
         real(lambdas(1:nlambdas), kind=4), rcode)
    
    do i = 1, ngases
       gas_indices(i) = i
    enddo
    call ncvpt (ncid, gasid, 1, ngases, gas_indices, rcode)

    !=============================================================================
    ! Define the START and COUNT arrays for each of the array variables.
    ! Fill in the values for other variables
    !=============================================================================
    ! write ps and ts
    ndimstart1 = (/ 1 /)
    ndimcount1 = (/ GC_nlayers+1 /)  
    call ncvpt (ncid, psid, ndimstart1, ndimcount1, real(pressures(0:GC_nlayers), kind=4), rcode)
    call ncvpt (ncid, tsid, ndimstart1, ndimcount1, real(temperatures(0:GC_nlayers), kind=4), rcode)
    
    ! write 1D with laydim aircol
    ndimstart1 = (/ 1 /)
    ndimcount1 = (/ GC_nlayers /)  
    call ncvpt (ncid, airid,  ndimstart1, ndimcount1, real(aircolumns(1:GC_nlayers), kind=4), rcode)
    call ncvpt (ncid, aer0id, ndimstart1, ndimcount1, real(taer_profile(1:GC_nlayers), kind=4),  rcode)
    call ncvpt (ncid, cld0id, ndimstart1, ndimcount1, real(tcld_profile(1:GC_nlayers), kind=4),  rcode)

    ! 2D, variables, laydim, gasdim
    ndimstart2 = (/ 1, 1 /)
    ndimcount2 = (/ GC_nlayers, ngases /)  
    call ncvpt (ncid, gascolid, ndimstart2, ndimcount2, gas_partialcolumns(1:GC_nlayers, 1:ngases), rcode)

    ! 4d variables, nsqdim, vzadim, azadim, szadim
    ndimstart4 = (/ 1, 1, 1, 1 /)
    ndimcount4 = (/ NSTOKESSQ, GC_n_view_angles, GC_n_azimuths, GC_n_sun_positions /) 
    if (do_brdf_surface) then
       call ncvpt (ncid, brdfid, ndimstart4, ndimcount4, &
            real(Total_brdf(1:NSTOKESSQ,1:GC_n_view_angles,1:GC_n_azimuths,1:GC_n_sun_positions), kind=4), rcode)  
       if (OUTPUT_WSABSA) then
          ndimstart1 = (/ 1 /)
          ndimcount1 = (/ 1 /)  
          call ncvpt (ncid, wsaid, ndimstart1, ndimcount1, real(WSA_CALCULATED, kind=4), rcode)
          call ncvpt (ncid, bsaid, ndimstart1, ndimcount1, real(BSA_CALCULATED, kind=4), rcode)
       endif
    endif

    !==============================================================================
    ! CLOSE the NetCDF file
    !==============================================================================
    
    call ncclos (ncid, rcode)
    
  END SUBROUTINE Create_netcdf_output_file
  
END MODULE GC_netcdf_module
