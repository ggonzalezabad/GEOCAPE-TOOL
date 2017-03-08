MODULE GC_netcdf_module

  USE GC_parameters_module, ONLY: max_ch_len
  USE GC_variables_module,  ONLY: latitude, longitude, year, month, day, utc, VLIDORT_ModIn, nlambdas, &
                                  lambdas, GC_nlayers, heights, pressures, temperatures, ngases,       &
                                  which_gases, gas_partialcolumns, aircolumns, cfrac, cld_tops,        &
                                  lambertian_cldalb, taer_profile, tcld_profile, opdeps, ssalbs,       &
                                  aer_opdeps, aer_ssalbs, cld_opdeps, cld_ssalbs, ground_ler,          &
                                  wind_speed, VLIDORT_FixIn, lambda_dw, clambdas, nclambdas,           &
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
                                  OUTPUT_WSABSA, WSA_CALCULATED, BSA_CALCULATED, didx, W,              &
                                  GC_user_altitudes, GC_n_user_levels, ilev, GC_user_levels
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
    call netcdf_wrt( netfname, tmperror, error)

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
    integer :: szadim, vzadim, azadim, gasdim, lvldim, laydim, wavdim, geodim, nsqdim, onedim, olvdim
    integer :: szaid, vzaid, azaid, gasid, lvlid, wavid, psid, tsid, airid, &
         aer0id, cld0id, radid, qid, uid, irradid, sfcid, sfcwfid, wswfid, cfracwfid, &
         sfcprswfid, aodwfid, assawfid, codwfid, cssawfid, sfcqwfid, wsqwfid, cfracqwfid, &
         sfcprsqwfid, aodqwfid, assaqwfid, codqwfid, cssaqwfid, sfcuwfid, wsuwfid, cfracuwfid,&
         sfcprsuwfid, aoduwfid, assauwfid, coduwfid, cssauwfid, gascolid, aodsid, assasid, &
         codsid, cssasid, scatwtid, odsid, ssasid, amfid, gaswfid, gasqwfid, &
         gasuwfid, tempwfid, fluxid, dfluxid, qfluxid, ufluxid, qdfluxid, udfluxid, brdfid, &
         wsaid, bsaid, outlevid
    integer, dimension(2)      :: gascol_dims, wavalt_dims
    integer, dimension(3)      :: wavgeolev_dims, wavszalev_dims
    integer, dimension(4)      :: wavaltgeolev_dims, wavgeogaslev_dims, brdfdim
    integer, dimension(5)      :: wavaltgeogaslev_dims
    integer, dimension(ngases) :: gas_indices
    integer, dimension(1)      :: ndimstart1, ndimcount1
    integer, dimension(2)      :: ndimstart2, ndimcount2
    integer, dimension(4)      :: ndimstart4, ndimcount4
    character(len=256)         :: gasstr
    character(len=30)          :: irrad_unitc, rad_unitc
    integer(kind=4)            :: ngeom

    ! -----------
    ! Fill values
    ! ----------
    INTEGER, parameter :: FILL_MODE   = 0 !Set to 1 to disable fill mode
 
   ! ---------------
    ! Chunking arrays
    ! ---------------
    INTEGER, dimension(1) :: chunk_1d
    INTEGER, dimension(2) :: chunk_2d
    INTEGER, dimension(3) :: chunk_3d
    INTEGER, dimension(4) :: chunk_4d

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
    
    rcode = NF_CREATE(netfname, OR(NF_CLOBBER,NF_NETCDF4), ncid)
    IF (rcode .NE. NF_NOERR) then
       tmperror =  ' error in netcdf_out: nccre failed'
       error = .true.; return
    endif

    ngeom = GC_n_sun_positions * GC_n_view_angles * GC_n_azimuths
    ! Create the dimensions of the dataset:
    rcode = NF_DEF_DIM (ncid, 'one', 1, onedim)
    rcode = NF_DEF_DIM (ncid, 'nsza', GC_n_sun_positions, szadim)
    rcode = NF_DEF_DIM (ncid, 'nvza', GC_n_view_angles, vzadim)
    rcode = NF_DEF_DIM (ncid, 'naza', GC_n_azimuths, azadim)
    rcode = NF_DEF_DIM (ncid, 'ngas', ngases, gasdim)
    rcode = NF_DEF_DIM (ncid, 'nlevel', GC_nlayers+1, lvldim)
    rcode = NF_DEF_DIM (ncid, 'noutputlevel', GC_n_user_levels, olvdim)
    rcode = NF_DEF_DIM (ncid, 'nlayer', GC_nlayers, laydim)
    IF ( .NOT. do_effcrs .AND. lambda_resolution /= 0.0d0 ) THEN
       rcode = NF_DEF_DIM (ncid, 'nw', nclambdas, wavdim)
    ELSE
       rcode = NF_DEF_DIM (ncid, 'nw', nlambdas, wavdim)
    END IF
    rcode = NF_DEF_DIM (ncid, 'ngeom',  ngeom, geodim)
    IF (do_brdf_surface) rcode = NF_DEF_DIM (ncid, 'nstokessq', NSTOKESSQ, nsqdim)
    
    !=============================================================================
    ! Create the coordinate (aka independent) variables:
    !============================================================================= 
    rcode = NF_DEF_VAR (ncid, 'Solarzenithangle',     NF_FLOAT, 1, szadim, szaid)
    rcode = NF_DEF_VAR_FILL(ncid, szaid, FILL_MODE, NF_FILL_REAL)
    chunk_1d(1) = szadim

    rcode = NF_DEF_VAR (ncid, 'Viewingzenithangle',   NF_FLOAT, 1, vzadim, vzaid)
    rcode = NF_DEF_VAR_FILL(ncid, vzaid, FILL_MODE, NF_FILL_REAL)

    rcode = NF_DEF_VAR (ncid, 'Relativeazimuthangle', NF_FLOAT, 1, azadim, azaid)
    rcode = NF_DEF_VAR_FILL(ncid, azaid, FILL_MODE, NF_FILL_REAL)

    rcode = NF_DEF_VAR (ncid, 'gas',                  NF_INT,   1, gasdim, gasid)
    rcode = NF_DEF_VAR_FILL(ncid, gasid, FILL_MODE, NF_FILL_INT)

    rcode = NF_DEF_VAR (ncid, 'zs',                   NF_FLOAT, 1, lvldim, lvlid)
    rcode = NF_DEF_VAR_FILL(ncid, lvlid, FILL_MODE, NF_FILL_REAL)

    rcode = NF_DEF_VAR (ncid, 'Wavelength',           NF_FLOAT, 1, wavdim, wavid)
    rcode = NF_DEF_VAR_FILL(ncid, wavid, FILL_MODE, NF_FILL_REAL)

    rcode = NF_DEF_VAR (ncid, 'outputlevels',         NF_FLOAT, 1, olvdim, outlevid)
    rcode = NF_DEF_VAR_FILL(ncid, outlevid, FILL_MODE, NF_FILL_REAL)
    
    !=============================================================================
    ! Create a vector containing the often referred to spacetime coordinates:
    !=============================================================================
    gascol_dims(1) = laydim; gascol_dims(2) = gasdim
    wavalt_dims(1) = wavdim; wavalt_dims(2) = laydim
    wavgeolev_dims(1) = wavdim; wavgeolev_dims(2) = geodim; wavgeolev_dims(3) = olvdim
    wavszalev_dims(1) = wavdim; wavszalev_dims(2) = szadim; wavszalev_dims(3) = olvdim
    
    wavaltgeolev_dims(1) = wavdim; wavaltgeolev_dims(2) = laydim; wavaltgeolev_dims(3) = geodim; wavaltgeolev_dims(4) = olvdim
    wavgeogaslev_dims(1) = wavdim; wavgeogaslev_dims(2) = geodim; wavgeogaslev_dims(3) = gasdim; wavgeogaslev_dims(4) = olvdim
    
    wavaltgeogaslev_dims(1) = wavdim; wavaltgeogaslev_dims(2) = laydim
    wavaltgeogaslev_dims(3) = geodim; wavaltgeogaslev_dims(4) = gasdim; wavaltgeogaslev_dims(5) = olvdim
    
    brdfdim(1) = nsqdim; brdfdim(2) = vzadim; brdfdim(3)=azadim; brdfdim(4) = szadim

    !=============================================================================
    ! Create the dependent variables:
    !=============================================================================
    ! ps, ts
    rcode = NF_DEF_VAR(ncid, 'ps', NF_FLOAT, 1, lvldim, psid)
    rcode = NF_DEF_VAR_FILL(ncid, psid, FILL_MODE, NF_FILL_REAL)
    rcode = NF_DEF_VAR(ncid, 'ts', NF_FLOAT, 1, lvldim, tsid)
    rcode = NF_DEF_VAR_FILL(ncid, tsid, FILL_MODE, NF_FILL_REAL)
    
    ! aircol, aods0, cods0
    rcode = NF_DEF_VAR(ncid, 'aircol', NF_FLOAT, 1, laydim, airid)
    rcode = NF_DEF_VAR_FILL(ncid, airid, FILL_MODE, NF_FILL_REAL)
    rcode = NF_DEF_VAR(ncid, 'aods0',  NF_FLOAT, 1, laydim, aer0id)
    rcode = NF_DEF_VAR_FILL(ncid, aer0id, FILL_MODE, NF_FILL_REAL)
    rcode = NF_DEF_VAR(ncid, 'cods0',  NF_FLOAT, 1, laydim, cld0id)
    rcode = NF_DEF_VAR_FILL(ncid, cld0id, FILL_MODE, NF_FILL_REAL)
        
    ! varaibles with 1D, wavdim
    rcode = NF_DEF_VAR(ncid, 'irradiance', NF_FLOAT, 1, wavdim, irradid)
    rcode = NF_DEF_VAR_FILL(ncid, irradid, FILL_MODE, NF_FILL_REAL)
    rcode = NF_DEF_VAR(ncid, 'surfalb',    NF_FLOAT, 1, wavdim, sfcid)
    rcode = NF_DEF_VAR_FILL(ncid, sfcid, FILL_MODE, NF_FILL_REAL)

    ! variables with 2D, laydim, gasdim
    rcode = NF_DEF_VAR(ncid,  'gascol', NF_DOUBLE, 2, gascol_dims, gascolid)
    rcode = NF_DEF_VAR_FILL(ncid, gascolid, FILL_MODE, NF_FILL_DOUBLE)
    
    ! variables with 2D, wavdim, laydim
    rcode = NF_DEF_VAR(ncid, 'ods',   NF_FLOAT, 2, wavalt_dims, odsid)
    rcode = NF_DEF_VAR_FILL(ncid, odsid, FILL_MODE, NF_FILL_REAL)
    rcode = NF_DEF_VAR(ncid, 'ssas',  NF_FLOAT, 2, wavalt_dims, ssasid)
    rcode = NF_DEF_VAR_FILL(ncid, ssasid, FILL_MODE, NF_FILL_REAL)
    rcode = NF_DEF_VAR(ncid, 'aods',  NF_FLOAT, 2, wavalt_dims, aodsid)
    rcode = NF_DEF_VAR_FILL(ncid, aodsid, FILL_MODE, NF_FILL_REAL)
    rcode = NF_DEF_VAR(ncid, 'assas', NF_FLOAT, 2, wavalt_dims, assasid)
    rcode = NF_DEF_VAR_FILL(ncid, assasid, FILL_MODE, NF_FILL_REAL)
    rcode = NF_DEF_VAR(ncid, 'cods',  NF_FLOAT, 2, wavalt_dims, codsid)
    rcode = NF_DEF_VAR_FILL(ncid, codsid, FILL_MODE, NF_FILL_REAL)
    rcode = NF_DEF_VAR(ncid, 'cssas', NF_FLOAT, 2, wavalt_dims, cssasid)
    rcode = NF_DEF_VAR_FILL(ncid, cssasid, FILL_MODE, NF_FILL_REAL)
    
    ! variables with 3D, wavdim, geodim, olvdim
    rcode = NF_DEF_VAR(ncid, 'radiance',    NF_FLOAT, 3, wavgeolev_dims, radid)
    rcode = NF_DEF_VAR_FILL(ncid, radid, FILL_MODE, NF_FILL_REAL)
    rcode = NF_DEF_VAR(ncid, 'flux',        NF_FLOAT, 3, wavszalev_dims, fluxid)
    rcode = NF_DEF_VAR_FILL(ncid, fluxid, FILL_MODE, NF_FILL_REAL)
    rcode = NF_DEF_VAR(ncid, 'direct_flux', NF_FLOAT, 3, wavszalev_dims, dfluxid)
    rcode = NF_DEF_VAR_FILL(ncid, dfluxid, FILL_MODE, NF_FILL_REAL)
    if (do_vector_calculation .and. do_StokesQU_output) then
       rcode = NF_DEF_VAR(ncid, 'q',            NF_FLOAT, 3, wavgeolev_dims, qid)
    rcode = NF_DEF_VAR_FILL(ncid, qid, FILL_MODE, NF_FILL_REAL)
       rcode = NF_DEF_VAR(ncid, 'u',            NF_FLOAT, 3, wavgeolev_dims, uid)
    rcode = NF_DEF_VAR_FILL(ncid, uid, FILL_MODE, NF_FILL_REAL)
       rcode = NF_DEF_VAR(ncid, 'qflux',        NF_FLOAT, 3, wavszalev_dims, qfluxid)
    rcode = NF_DEF_VAR_FILL(ncid, qfluxid, FILL_MODE, NF_FILL_REAL)
       rcode = NF_DEF_VAR(ncid, 'uflux',        NF_FLOAT, 3, wavszalev_dims, ufluxid)
    rcode = NF_DEF_VAR_FILL(ncid, ufluxid, FILL_MODE, NF_FILL_REAL)
       rcode = NF_DEF_VAR(ncid, 'qdirect_flux', NF_FLOAT, 3, wavszalev_dims, qdfluxid)
    rcode = NF_DEF_VAR_FILL(ncid, qdfluxid, FILL_MODE, NF_FILL_REAL)
       rcode = NF_DEF_VAR(ncid, 'udirect_flux', NF_FLOAT, 3, wavszalev_dims, udfluxid)
    rcode = NF_DEF_VAR_FILL(ncid, udfluxid, FILL_MODE, NF_FILL_REAL)
    endif
    
    if (do_Jacobians) then
       
       if (use_lambertian) then
          rcode = NF_DEF_VAR(ncid, 'surfalb_jac',   NF_FLOAT, 3, wavgeolev_dims, sfcwfid)
          rcode = NF_DEF_VAR_FILL(ncid, sfcwfid, FILL_MODE, NF_FILL_REAL)
       else 
          rcode = NF_DEF_VAR(ncid, 'windspeed_jac', NF_FLOAT, 3, wavgeolev_dims, wswfid)
          rcode = NF_DEF_VAR_FILL(ncid, wswfid, FILL_MODE, NF_FILL_REAL)
       endif
       
       if (do_cfrac_Jacobians) then
            rcode = NF_DEF_VAR(ncid, 'cfrac_jac',     NF_FLOAT, 3, wavgeolev_dims, cfracwfid)
            rcode = NF_DEF_VAR_FILL(ncid, cfracwfid, FILL_MODE, NF_FILL_REAL)
         endif
       if (do_sfcprs_Jacobians) then
            rcode = NF_DEF_VAR(ncid, 'sfcprs_jac',    NF_FLOAT, 3, wavgeolev_dims, sfcprswfid)
            rcode = NF_DEF_VAR_FILL(ncid, sfcprswfid, FILL_MODE, NF_FILL_REAL)
         endif
       if (do_aer_columnwf .and. do_aod_Jacobians) then  !No column jacobians yet
          rcode = NF_DEF_VAR(ncid, 'aodcolwf_jac',  NF_FLOAT, 3, wavgeolev_dims, aodwfid)
          rcode = NF_DEF_VAR_FILL(ncid, aodwfid, FILL_MODE, NF_FILL_REAL)
       endif
       if (do_aer_columnwf .and. do_assa_Jacobians) then !No column jacobians yet
          rcode = NF_DEF_VAR(ncid, 'assacolwf_jac', NF_FLOAT, 3, wavgeolev_dims, assawfid)
          rcode = NF_DEF_VAR_FILL(ncid, assawfid, FILL_MODE, NF_FILL_REAL)
       endif
       if (do_cld_columnwf .and. do_cod_Jacobians) then  !No column jacobians yet
          rcode = NF_DEF_VAR(ncid, 'codcolwf_jac',  NF_FLOAT, 3, wavgeolev_dims, codwfid)
          rcode = NF_DEF_VAR_FILL(ncid, codwfid, FILL_MODE, NF_FILL_REAL)
       endif
       if (do_cld_columnwf .and. do_cssa_Jacobians) then !No column jacobians yet
          rcode = NF_DEF_VAR(ncid, 'cssacolwf_jac', NF_FLOAT, 3, wavgeolev_dims, cssawfid)
          rcode = NF_DEF_VAR_FILL(ncid, cssawfid, FILL_MODE, NF_FILL_REAL)
       endif
       
    endif
    
    if (do_QU_Jacobians) then
       
       if (use_lambertian) then
          rcode = NF_DEF_VAR(ncid, 'surfalb_qjac',  NF_FLOAT, 3, wavgeolev_dims, sfcqwfid)
          rcode = NF_DEF_VAR_FILL(ncid, sfcqwfid, FILL_MODE, NF_FILL_REAL)
          rcode = NF_DEF_VAR(ncid, 'surfalb_ujac',  NF_FLOAT, 3, wavgeolev_dims, sfcuwfid)
          rcode = NF_DEF_VAR_FILL(ncid, sfcuwfid, FILL_MODE, NF_FILL_REAL)
       else 
          rcode = NF_DEF_VAR(ncid, 'windspeed_qjac', NF_FLOAT, 3, wavgeolev_dims, wsqwfid)
          rcode = NF_DEF_VAR(ncid, 'windspeed_ujac', NF_FLOAT, 3, wavgeolev_dims, wsuwfid)
       endif
       
       if (do_cfrac_Jacobians) then
          rcode = NF_DEF_VAR(ncid, 'cfrac_qjac',  NF_FLOAT, 3, wavgeolev_dims, cfracqwfid)
          rcode = NF_DEF_VAR_FILL(ncid, cfracqwfid, FILL_MODE, NF_FILL_REAL)
          rcode = NF_DEF_VAR(ncid, 'cfrac_ujac',  NF_FLOAT, 3, wavgeolev_dims, cfracuwfid)
          rcode = NF_DEF_VAR_FILL(ncid, cfracuwfid, FILL_MODE, NF_FILL_REAL)
       endif
       if (do_sfcprs_Jacobians) then
          rcode = NF_DEF_VAR(ncid, 'sfcprs_qjac',  NF_FLOAT, 3, wavgeolev_dims, sfcprsqwfid)
          rcode = NF_DEF_VAR_FILL(ncid, sfcprsqwfid, FILL_MODE, NF_FILL_REAL)
          rcode = NF_DEF_VAR(ncid, 'sfcprs_ujac',  NF_FLOAT, 3, wavgeolev_dims, sfcprsuwfid)
          rcode = NF_DEF_VAR_FILL(ncid, sfcprsuwfid, FILL_MODE, NF_FILL_REAL)
       endif
       if (do_aer_columnwf .and. do_aod_Jacobians) then  !No column jacobians yet
          rcode = NF_DEF_VAR(ncid, 'aodcol_qjac',  NF_FLOAT, 3, wavgeolev_dims, aodqwfid)
          rcode = NF_DEF_VAR_FILL(ncid, aodqwfid, FILL_MODE, NF_FILL_REAL)
          rcode = NF_DEF_VAR(ncid, 'aodcol_ujac',  NF_FLOAT, 3, wavgeolev_dims, aoduwfid)
          rcode = NF_DEF_VAR_FILL(ncid, aoduwfid, FILL_MODE, NF_FILL_REAL)
       endif
       if (do_aer_columnwf .and. do_assa_Jacobians) then !No column jacobians yet
          rcode = NF_DEF_VAR(ncid, 'assacol_qjac', NF_FLOAT, 3, wavgeolev_dims, assaqwfid)
          rcode = NF_DEF_VAR_FILL(ncid, assaqwfid, FILL_MODE, NF_FILL_REAL)
          rcode = NF_DEF_VAR(ncid, 'assacol_ujac', NF_FLOAT, 3, wavgeolev_dims, assauwfid)
          rcode = NF_DEF_VAR_FILL(ncid, assauwfid, FILL_MODE, NF_FILL_REAL)
       endif
       if (do_cld_columnwf .and. do_cod_Jacobians) then  !No column jacobians yet
          rcode = NF_DEF_VAR(ncid, 'codcol_qjac',  NF_FLOAT, 3, wavgeolev_dims, codqwfid)
          rcode = NF_DEF_VAR_FILL(ncid, codqwfid, FILL_MODE, NF_FILL_REAL)
          rcode = NF_DEF_VAR(ncid, 'codcol_ujac',  NF_FLOAT, 3, wavgeolev_dims, coduwfid)
          rcode = NF_DEF_VAR_FILL(ncid, coduwfid, FILL_MODE, NF_FILL_REAL)
       endif
       if (do_cld_columnwf .and. do_cssa_Jacobians) then !No column jacobians yet
          rcode = NF_DEF_VAR(ncid, 'cssacol_qjac', NF_FLOAT, 3, wavgeolev_dims, cssaqwfid)
          rcode = NF_DEF_VAR_FILL(ncid, cssaqwfid, FILL_MODE, NF_FILL_REAL)
          rcode = NF_DEF_VAR(ncid, 'cssacol_ujac', NF_FLOAT, 3, wavgeolev_dims, cssauwfid)
          rcode = NF_DEF_VAR_FILL(ncid, cssauwfid, FILL_MODE, NF_FILL_REAL)
       endif
       
    endif
    
    ! variables with 4D, wavdim, laydim, geodim, olvdim
    if (do_AMF_calculation) then
       rcode = NF_DEF_VAR(ncid, 'scatweights', NF_FLOAT, 4, wavaltgeolev_dims, scatwtid)
       rcode = NF_DEF_VAR_FILL(ncid, scatwtid, FILL_MODE, NF_FILL_REAL)
    endif
    
    if (do_Jacobians) then
       if (do_T_Jacobians) then
            rcode = NF_DEF_VAR(ncid, 't_jac',    NF_FLOAT, 4, wavaltgeolev_dims, tempwfid)
            rcode = NF_DEF_VAR_FILL(ncid, tempwfid, FILL_MODE, NF_FILL_REAL)
         endif
       if (.not. do_aer_columnwf .and. do_aod_Jacobians) then
            rcode = NF_DEF_VAR(ncid, 'aod_jac',  NF_FLOAT, 4, wavaltgeolev_dims, aodwfid)  
            rcode = NF_DEF_VAR_FILL(ncid, aodwfid, FILL_MODE, NF_FILL_REAL)
         endif
       if (.not. do_aer_columnwf .and. do_assa_Jacobians) then
            rcode = NF_DEF_VAR(ncid, 'assa_jac', NF_FLOAT, 4, wavaltgeolev_dims, assawfid)   
            rcode = NF_DEF_VAR_FILL(ncid, assawfid, FILL_MODE, NF_FILL_REAL)
         endif
       if (.not. do_cld_columnwf .and. do_cod_Jacobians) then
            rcode = NF_DEF_VAR(ncid, 'cod_jac',  NF_FLOAT, 4, wavaltgeolev_dims, codwfid)  
            rcode = NF_DEF_VAR_FILL(ncid, codwfid, FILL_MODE, NF_FILL_REAL)
         endif
       if (.not. do_cld_columnwf .and. do_cssa_Jacobians) then
            rcode = NF_DEF_VAR(ncid, 'cssa_jac', NF_FLOAT, 4, wavaltgeolev_dims, cssawfid)  
            rcode = NF_DEF_VAR_FILL(ncid, cssawfid, FILL_MODE, NF_FILL_REAL)
         endif
    endif
    
    if (do_QU_Jacobians) then
       if (.not. do_aer_columnwf .and. do_aod_Jacobians) then 
          rcode = NF_DEF_VAR(ncid, 'aod_qjac', NF_FLOAT, 4, wavaltgeolev_dims, aodqwfid)  
          rcode = NF_DEF_VAR_FILL(ncid, aodqwfid, FILL_MODE, NF_FILL_REAL)
          rcode = NF_DEF_VAR(ncid, 'aod_ujac', NF_FLOAT, 4, wavaltgeolev_dims, aoduwfid)  
          rcode = NF_DEF_VAR_FILL(ncid, aoduwfid, FILL_MODE, NF_FILL_REAL)
       endif
       if (.not. do_aer_columnwf .and. do_assa_Jacobians) then
          rcode = NF_DEF_VAR(ncid, 'assa_qjac', NF_FLOAT, 4, wavaltgeolev_dims, assaqwfid)   
          rcode = NF_DEF_VAR_FILL(ncid, assaqwfid, FILL_MODE, NF_FILL_REAL)
          rcode = NF_DEF_VAR(ncid, 'assa_ujac', NF_FLOAT, 4, wavaltgeolev_dims, assauwfid)   
          rcode = NF_DEF_VAR_FILL(ncid, assauwfid, FILL_MODE, NF_FILL_REAL)
       endif
       if (.not. do_cld_columnwf .and. do_cod_Jacobians) then
          rcode = NF_DEF_VAR(ncid, 'cod_qjac', NF_FLOAT, 4, wavaltgeolev_dims, codqwfid)  
          rcode = NF_DEF_VAR_FILL(ncid, codqwfid, FILL_MODE, NF_FILL_REAL)
          rcode = NF_DEF_VAR(ncid, 'cod_ujac', NF_FLOAT, 4, wavaltgeolev_dims, coduwfid)  
          rcode = NF_DEF_VAR_FILL(ncid, coduwfid, FILL_MODE, NF_FILL_REAL)
       endif
       if (.not. do_cld_columnwf .and. do_cssa_Jacobians) then
          rcode = NF_DEF_VAR(ncid, 'cssa_qjac', NF_FLOAT, 4, wavaltgeolev_dims, cssaqwfid)  
          rcode = NF_DEF_VAR_FILL(ncid, cssaqwfid, FILL_MODE, NF_FILL_REAL)
          rcode = NF_DEF_VAR(ncid, 'cssa_ujac', NF_FLOAT, 4, wavaltgeolev_dims, cssauwfid) 
          rcode = NF_DEF_VAR_FILL(ncid, cssauwfid, FILL_MODE, NF_FILL_REAL)
       endif
    endif
    
    ! variables with 4D, wavdim, geodim, gasdim, olvdim
    if (do_AMF_calculation) then
       rcode = NF_DEF_VAR(ncid, 'amf', NF_FLOAT, 4, wavgeogaslev_dims, amfid)
       rcode = NF_DEF_VAR_FILL(ncid, amfid, FILL_MODE, NF_FILL_REAL)
    endif
    
    ! variables with 4D, wavdim, laydim, geodim, gasdim, olvdim
    if (do_Jacobians) then
       rcode = NF_DEF_VAR(ncid,  'gas_jac', NF_FLOAT, 5, wavaltgeogaslev_dims, gaswfid)
       rcode = NF_DEF_VAR_FILL(ncid, gaswfid, FILL_MODE, NF_FILL_REAL)
       if (do_QU_Jacobians) then
          rcode = NF_DEF_VAR(ncid,  'gas_qjac', NF_FLOAT, 5, wavaltgeogaslev_dims, gasqwfid)
          rcode = NF_DEF_VAR_FILL(ncid, gasqwfid, FILL_MODE, NF_FILL_REAL)
          rcode = NF_DEF_VAR(ncid,  'gas_ujac', NF_FLOAT, 5, wavaltgeogaslev_dims, gasuwfid)
          rcode = NF_DEF_VAR_FILL(ncid, gasuwfid, FILL_MODE, NF_FILL_REAL)
       endif
    endif
    
    ! BRDF variables onedim, nsqdim, nvza, naza, nsza
    if (do_brdf_surface) then
    ! WSA, BSA amd BRDF
       rcode = NF_DEF_VAR(ncid, 'WSA', NF_FLOAT, 1, onedim, wsaid)
       rcode = NF_DEF_VAR_FILL(ncid, wsaid, FILL_MODE, NF_FILL_REAL)
       rcode = NF_DEF_VAR(ncid, 'BSA', NF_FLOAT, 1, onedim, bsaid)
       rcode = NF_DEF_VAR_FILL(ncid, bsaid, FILL_MODE, NF_FILL_REAL)
       rcode = NF_DEF_VAR(ncid, 'BRDF', NF_FLOAT, 4, brdfdim, brdfid)
       rcode = NF_DEF_VAR_FILL(ncid, brdfid, FILL_MODE, NF_FILL_REAL)
    endif

    !=============================================================================
    ! Assign attributes (meta-data) to the various variables:
    !============================================================================= 
    rcode = NF_PUT_ATT_REAL (ncid, NF_GLOBAL, 'lon',          NF_FLOAT, 1, real(longitude, kind=4))
    rcode = NF_PUT_ATT_REAL (ncid, NF_GLOBAL, 'lat',          NF_FLOAT, 1, real(latitude, kind=4))
    rcode = NF_PUT_ATT_INT (ncid, NF_GLOBAL, 'month',        NF_INT,  1, month)
    rcode = NF_PUT_ATT_INT (ncid, NF_GLOBAL, 'year',         NF_INT,  1, year)
    rcode = NF_PUT_ATT_INT (ncid, NF_GLOBAL, 'day',          NF_INT,  1, day)
    rcode = NF_PUT_ATT_REAL (ncid, NF_GLOBAL, 'utc',          NF_FLOAT, 1, real(utc, kind=4))
    rcode = NF_PUT_ATT_REAL (ncid, NF_GLOBAL, 'cfrac',        NF_FLOAT, 1, real(cfrac, kind=4))
    rcode = NF_PUT_ATT_REAL (ncid, NF_GLOBAL, 'zcldtop',      NF_FLOAT, 1, real(cld_tops, kind=4))
    rcode = NF_PUT_ATT_REAL (ncid, NF_GLOBAL, 'cldalb',       NF_FLOAT, 1, real(lambertian_cldalb, kind=4))
    rcode = NF_PUT_ATT_REAL (ncid, NF_GLOBAL, 'userlvl',      NF_FLOAT, 1, real(VLIDORT_ModIn%MUserVal%TS_USER_LEVELS, kind=4))
    rcode = NF_PUT_ATT_REAL (ncid, NF_GLOBAL, 'windspeed',    NF_FLOAT, 1, real(wind_speed, kind=4))
    rcode = NF_PUT_ATT_REAL (ncid, NF_GLOBAL, 'lam_reso',     NF_FLOAT, 1, real(lambda_resolution, kind=4))
    rcode = NF_PUT_ATT_REAL (ncid, NF_GLOBAL, 'aer_reflam',   NF_FLOAT, 1, real(aer_reflambda, kind=4))
    rcode = NF_PUT_ATT_REAL (ncid, NF_GLOBAL, 'cld_reflam',   NF_FLOAT, 1, real(cld_reflambda, kind=4))
    rcode = NF_PUT_ATT_REAL (ncid, NF_GLOBAL, 'lam_dw',       NF_FLOAT, 1, real(lambda_dw, kind=4))
    rcode = NF_PUT_ATT_REAL (ncid, NF_GLOBAL, 'lam_dfw',      NF_FLOAT, 1, real(lambda_dfw, kind=4))
    rcode = NF_PUT_ATT_INT (ncid, NF_GLOBAL, 'nstreams',     NF_INT,  1, VLIDORT_FixIn%Cont%TS_NSTREAMS)
    rcode = NF_PUT_ATT_INT (ncid, NF_GLOBAL, 'aercld_nmoms', NF_INT,  1, VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT)
    rcode = NF_PUT_ATT_INT (ncid, NF_GLOBAL, 'nsza',         NF_INT,  1, GC_n_sun_positions)
    rcode = NF_PUT_ATT_INT (ncid, NF_GLOBAL, 'nvza',         NF_INT,  1, GC_n_view_angles)
    rcode = NF_PUT_ATT_INT (ncid, NF_GLOBAL, 'naza',         NF_INT,  1, GC_n_azimuths)
    rcode = NF_PUT_ATT_INT (ncid, NF_GLOBAL, 'ngeometries',  NF_INT,  1, ngeom)
    rcode = NF_PUT_ATT_INT (ncid, NF_GLOBAL, 'noutputlevels',NF_INT,  1, VLIDORT_FixIn%UserVal%TS_N_USER_LEVELS)

    IF ( .NOT. do_effcrs .AND. lambda_resolution /= 0.0d0 ) THEN
       rcode = NF_PUT_ATT_INT (ncid, NF_GLOBAL, 'nwavelengths', NF_INT,  1, nclambdas)
    ELSE
       rcode = NF_PUT_ATT_INT (ncid, NF_GLOBAL, 'nwavelengths', NF_INT,  1, nlambdas)
    END IF
    rcode = NF_PUT_ATT_INT (ncid, NF_GLOBAL, 'nlayers',      NF_INT,  1, GC_nlayers)
    rcode = NF_PUT_ATT_INT (ncid, NF_GLOBAL, 'ngas',         NF_INT,  1, ngases)
    
    if (use_lambertian) then
       rcode = NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'user_lambertian', NF_INT1, 1, 1)
    else
       rcode = NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'user_lambertian', NF_INT1, 1, 0)
    endif
    if (do_lambertian_cld) then
       rcode = NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'do_lambertian_cld', NF_INT1, 1, 1)
    else
       rcode = NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'do_lambertian_cld', NF_INT1, 1, 0)
    endif
    if (do_effcrs) then
       rcode = NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'do_effcrs', NF_INT1, 1, 1)
    else
       rcode = NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'do_effcrs', NF_INT1, 1, 0)
    endif
    if (use_wavelength) then
       rcode = NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'use_wavelength', NF_INT1, 1, 1)
    else
       rcode = NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'use_wavelength', NF_INT1, 1, 0)
    endif
    if (VLIDORT_FixIn%Bool%TS_DO_UPWELLING) then
       rcode = NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'do_upwelling', NF_INT1, 1, 1)
    else
       rcode = NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'do_upwelling', NF_INT1, 1, 0)
    endif
    if (do_normalized_WFoutput) then
       rcode = NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'do_norm_WFoutput', NF_INT1, 1, 1)
    else
       rcode = NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'do_norm_WFoutput', NF_INT1, 1, 0)
    endif
    if (do_normalized_radiance) then
       rcode = NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'do_norm_radiance', NF_INT1, 1, 1)
    else
       rcode = NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'do_norm_radiance', NF_INT1, 1, 0)
    endif
    if (use_solar_photons) then
       rcode = NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'use_solar_photons', NF_INT1, 1, 1)
    else
       rcode = NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'use_solar_photons', NF_INT1, 1, 0)
    endif
    if (do_vector_calculation) then
       rcode = NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'do_vector_calc', NF_INT1, 1, 1)
    else
       rcode = NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'do_vector_calc', NF_INT1, 1, 0)
    endif
    if (do_StokesQU_output) then
       rcode = NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'do_QU_output', NF_INT1, 1, 1)
    else
       rcode = NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'do_QU_output', NF_INT1, 1, 0)
    endif
    if (do_Jacobians) then
       rcode = NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'do_Jacobian', NF_INT1, 1, 1)
    else
       rcode = NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'do_Jacobian', NF_INT1, 1, 0)
    endif
    if (do_QU_Jacobians) then
       rcode = NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'do_QU_Jacobian', NF_INT1, 1, 1)
    else
       rcode = NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'do_QU_Jacobian', NF_INT1, 1, 0)
    endif
    if (do_AMF_calculation) then
       rcode = NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'do_AMF_calc', NF_INT1, 1, 1)
    else
       rcode = NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'do_AMF_calc', NF_INT1, 1, 0)
    endif
    if (do_T_Jacobians) then
       rcode = NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'do_T_Jacobian', NF_INT1, 1, 1)
    else
       rcode = NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'do_T_Jacobian', NF_INT1, 1, 0)
    endif
    if (do_sfcprs_Jacobians) then
       rcode = NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'do_sfcprs_Jacobian', NF_INT1, 1, 1)
    else
       rcode = NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'do_sfcprs_Jacobian', NF_INT1, 1, 0)
    endif
    if (do_aod_Jacobians) then
       rcode = NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'do_aod_Jacobian', NF_INT1, 1, 1)
    else
       rcode = NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'do_aod_Jacobian', NF_INT1, 1, 0)
    endif
    if (do_assa_Jacobians) then
       rcode = NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'do_assa_Jacobian', NF_INT1, 1, 1)
    else
       rcode = NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'do_assa_Jacobian', NF_INT1, 1, 0)
    endif
    if (do_cod_Jacobians) then
       rcode = NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'do_cod_Jacobian', NF_INT1, 1, 1)
    else
       rcode = NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'do_cod_Jacobian', NF_INT1, 1, 0)
    endif
    if (do_cssa_Jacobians) then
       rcode = NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'do_cssa_Jacobian', NF_INT1, 1, 1)
    else
       rcode = NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'do_cssa_Jacobian', NF_INT1, 1, 0)
    endif
    if (do_cfrac_Jacobians) then
       rcode = NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'do_cfrac_Jacobian', NF_INT1, 1, 1)
    else
       rcode = NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'do_cfrac_Jacobian', NF_INT1, 1, 0)
    endif
    if (do_aer_columnwf) then
       rcode = NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'do_aer_columnwf', NF_INT1, 1, 1)
    else
       rcode = NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'do_aer_columnwf', NF_INT1, 1, 0)
    endif
    if (do_cld_columnwf) then
       rcode = NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'do_cld_columnwf', NF_INT1, 1, 1)
    else
       rcode = NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'do_cld_columnwf', NF_INT1, 1, 0)
    endif
    if (do_brdf_surface) then
       rcode = NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'do_brdf', NF_INT1, 1, 1)
    else
       rcode = NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'do_brdf', NF_INT1, 1, 0)
    endif
    if (OUTPUT_WSABSA) then
       rcode = NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'do_output_wsabsa', NF_INT1, 1, 1)
    else
       rcode = NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'do_output_wsabsa', NF_INT1, 1, 0)
    endif
    
    ! write the list of gases in one string
    write(gasstr, '(10(I2,A1,A4,A1))') (i,':',which_gases(i),',', i=1, ngases)
    nlen=LEN(trim(gasstr)) ; gasstr=gasstr(1:nlen-1)
    rcode = NF_PUT_ATT_TEXT (ncid, NF_GLOBAL, 'gases', ncchar,  len(trim(gasstr)), trim(gasstr))

    ! Units for variables
    rcode = NF_PUT_ATT_TEXT (ncid, NF_GLOBAL, 'windspeed_units', ncchar, 3, 'm/s')
    rcode = NF_PUT_ATT_TEXT (ncid, NF_GLOBAL, 'lonlat_units',    ncchar, 7, 'degrees')
    rcode = NF_PUT_ATT_TEXT (ncid, NF_GLOBAL, 'zcldtop_units',   ncchar, 2, 'km')
    rcode = NF_PUT_ATT_TEXT (ncid, szaid,    'units', ncchar, 7, 'degrees')
    rcode = NF_PUT_ATT_TEXT (ncid, vzaid,    'units', ncchar, 7, 'degrees')
    rcode = NF_PUT_ATT_TEXT (ncid, azaid,    'units', ncchar, 7, 'degrees')
    rcode = NF_PUT_ATT_TEXT (ncid, lvlid,    'units', ncchar, 2, 'km')
    rcode = NF_PUT_ATT_TEXT (ncid, outlevid, 'units', ncchar, 8, 'unitless')
    if (use_wavelength) then
       rcode = NF_PUT_ATT_TEXT(ncid,wavid,'units', ncchar,2, 'nm',rcode)
    else 
       rcode = NF_PUT_ATT_TEXT(ncid,wavid,'units', ncchar,4, 'cm-1',rcode)
    endif
    if (use_solar_photons) then
       irrad_unitc='photons/cm2/nm/s'
    else
       irrad_unitc='W/m2/cm-1'
    endif
    rcode = NF_PUT_ATT_TEXT(ncid, irradid, 'units', ncchar, len_trim(irrad_unitc), irrad_unitc)
    if (do_normalized_radiance) then
       rad_unitc='unitless'
    else
       rad_unitc=TRIM(irrad_unitc)//'/sr'
    endif
    rcode = NF_PUT_ATT_TEXT(ncid, radid,    'units', ncchar, len_trim(rad_unitc), rad_unitc) 
    rcode = NF_PUT_ATT_TEXT(ncid, gascolid, 'units', ncchar, 13, 'molecules/cm2')
    rcode = NF_PUT_ATT_TEXT(ncid, airid,    'units', ncchar, 13, 'molecules/cm2')
    rcode = NF_PUT_ATT_TEXT(ncid, psid,     'units', ncchar, 3,  'hPa')
    rcode = NF_PUT_ATT_TEXT(ncid, tsid,     'units', ncchar, 1,  'K')    

    !=============================================================================
    ! Get out of 'define' mode, and into 'data' mode
    !=============================================================================
    rcode = NF_ENDDEF (ncid)

    !=============================================================================
    ! Fill in the values of the non time varying variables.
    ! Remember, we're still in "initialization" here - this is only done the
    ! first time through.
    !=============================================================================
    rcode = NF_PUT_VARA_REAL (ncid, szaid, 1, GC_n_sun_positions, &
         real(VLIDORT_ModIn%MSunrays%TS_SZANGLES(1:GC_n_sun_positions), kind=4))
    rcode = NF_PUT_VARA_REAL (ncid, vzaid, 1, GC_n_view_angles, &
         real(VLIDORT_ModIn%MUserVal%TS_USER_VZANGLES_INPUT(1:GC_n_view_angles), kind=4))
    rcode = NF_PUT_VARA_REAL (ncid, azaid, 1, GC_n_azimuths, &
         real(VLIDORT_ModIn%MUserVal%TS_USER_RELAZMS(1:GC_n_azimuths), kind=4))
    rcode = NF_PUT_VARA_REAL (ncid, lvlid, 1, GC_nlayers+1, &
         real(heights(0:GC_nlayers), kind=4))
    IF ( .NOT. do_effcrs .AND. lambda_resolution /= 0.0d0) THEN
       rcode = NF_PUT_VARA_REAL (ncid, wavid, 1, nclambdas, &
            real(clambdas(1:nclambdas), kind=4))
    ELSE
       rcode = NF_PUT_VARA_REAL (ncid, wavid, 1, nlambdas, &
            real(lambdas(1:nlambdas), kind=4))
    END IF
    
    do i = 1, ngases
       gas_indices(i) = i
    enddo
    rcode = NF_PUT_VARA_INT (ncid, gasid, 1, ngases, gas_indices)
    rcode = NF_PUT_VARA_REAL (ncid, outlevid, 1, VLIDORT_FixIn%UserVal%TS_N_USER_LEVELS, &
         real(VLIDORT_ModIn%MUserVal%TS_USER_LEVELS(1:VLIDORT_FixIn%UserVal%TS_N_USER_LEVELS), kind=4))

    !=============================================================================
    ! Define the START and COUNT arrays for each of the array variables.
    ! Fill in the values for other variables
    !=============================================================================
    ! write ps and ts
    ndimstart1 = (/ 1 /)
    ndimcount1 = (/ GC_nlayers+1 /)  
    rcode = NF_PUT_VARA_REAL (ncid, psid, ndimstart1, ndimcount1, real(pressures(0:GC_nlayers), kind=4))
    rcode = NF_PUT_VARA_REAL (ncid, tsid, ndimstart1, ndimcount1, real(temperatures(0:GC_nlayers), kind=4))
    
    ! write 1D with laydim aircol
    ndimstart1 = (/ 1 /)
    ndimcount1 = (/ GC_nlayers /)  
    rcode = NF_PUT_VARA_REAL (ncid, airid,  ndimstart1, ndimcount1, real(aircolumns(1:GC_nlayers), kind=4))
    rcode = NF_PUT_VARA_REAL (ncid, aer0id, ndimstart1, ndimcount1, real(taer_profile(1:GC_nlayers), kind=4),  rcode)
    rcode = NF_PUT_VARA_REAL (ncid, cld0id, ndimstart1, ndimcount1, real(tcld_profile(1:GC_nlayers), kind=4),  rcode)

    ! 2D, variables, laydim, gasdim
    ndimstart2 = (/ 1, 1 /)
    ndimcount2 = (/ GC_nlayers, ngases /)  
    rcode = NF_PUT_VARA_DOUBLE (ncid, gascolid, ndimstart2, ndimcount2, gas_partialcolumns(1:GC_nlayers, 1:ngases))

    ! 4d variables, nsqdim, vzadim, azadim, szadim
    if (do_brdf_surface) then
       ndimstart4 = (/ 1, 1, 1, 1 /)
       ndimcount4 = (/ NSTOKESSQ, GC_n_view_angles, GC_n_azimuths, GC_n_sun_positions /) 
       rcode = NF_PUT_VARA_REAL (ncid, brdfid, ndimstart4, ndimcount4, &
            real(Total_brdf(1:NSTOKESSQ,1:GC_n_view_angles,1:GC_n_azimuths,1:GC_n_sun_positions), kind=4))  
       if (OUTPUT_WSABSA) then
          ndimstart1 = (/ 1 /)
          ndimcount1 = (/ 1 /)  
          rcode = NF_PUT_VARA_REAL (ncid, wsaid, ndimstart1, ndimcount1, real(WSA_CALCULATED, kind=4))
          rcode = NF_PUT_VARA_REAL (ncid, bsaid, ndimstart1, ndimcount1, real(BSA_CALCULATED, kind=4))
       endif
    endif

    !==============================================================================
    ! CLOSE the NetCDF file
    !==============================================================================
    
    rcode = NF_CLOSE(ncid)
    
  END SUBROUTINE Create_netcdf_output_file
  
END MODULE GC_netcdf_module
