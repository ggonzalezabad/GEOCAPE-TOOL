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
    call netcdf_wrt( netfname)

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
    CHARACTER(LEN=max_ch_len)    :: netfname
    CHARACTER(LEN=25), parameter :: location='Create_netcdf_output_file'

    integer :: ncid, nlen, i
    integer :: szadim, vzadim, azadim, gasdim, lvldim, laydim, wavdim, geodim, nsqdim, onedim, olvdim
    integer :: szadimid, vzadimid, azadimid, gasdimid, lvldimid, laydimid, wavdimid, &
         geodimid, nsqdimid, onedimid, olvdimid
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
    INTEGER, parameter, dimension(1) :: members_chunk_1 = 4096
    INTEGER, parameter, dimension(1) :: members_chunk_2 = 2048
    INTEGER, parameter, dimension(1) :: members_chunk_4 = 1024
    INTEGER, parameter, dimension(1) :: members_chunk_8 =  512
    INTEGER, dimension(1) :: chunk_1d
    INTEGER, dimension(2) :: chunk_2d
    INTEGER, dimension(3) :: chunk_3d
    INTEGER, dimension(4) :: chunk_4d
    INTEGER, dimension(5) :: chunk_5d
    INTEGER :: max_chunk

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
    
    CALL netcdf_handle_error(location,NF_CREATE(netfname, OR(NF_CLOBBER,NF_NETCDF4), ncid))

    ngeom = GC_n_sun_positions * GC_n_view_angles * GC_n_azimuths
    ! Create the dimensions of the dataset:
    CALL netcdf_handle_error(location,NF_DEF_DIM (ncid, 'one', 1, onedimid)); onedim = 1
    CALL netcdf_handle_error(location,NF_DEF_DIM (ncid, 'nsza', GC_n_sun_positions, szadimid)); szadim = GC_n_sun_positions
    CALL netcdf_handle_error(location,NF_DEF_DIM (ncid, 'nvza', GC_n_view_angles, vzadimid)); vzadim = GC_n_view_angles
    CALL netcdf_handle_error(location,NF_DEF_DIM (ncid, 'naza', GC_n_azimuths, azadimid)); azadim = GC_n_azimuths
    CALL netcdf_handle_error(location,NF_DEF_DIM (ncid, 'ngas', ngases, gasdimid)); gasdim = ngases
    CALL netcdf_handle_error(location,NF_DEF_DIM (ncid, 'nlevel', GC_nlayers+1, lvldimid)); lvldim = GC_nlayers+1
    CALL netcdf_handle_error(location,NF_DEF_DIM (ncid, 'noutputlevel', GC_n_user_levels, olvdimid)); olvdim = GC_n_user_levels
    CALL netcdf_handle_error(location,NF_DEF_DIM (ncid, 'nlayer', GC_nlayers, laydimid)); laydim = GC_nlayers
    IF ( .NOT. do_effcrs .AND. lambda_resolution /= 0.0d0 ) THEN
       CALL netcdf_handle_error(location,NF_DEF_DIM (ncid, 'nw', nclambdas, wavdimid)); wavdim = nclambdas
    ELSE
       CALL netcdf_handle_error(location,NF_DEF_DIM (ncid, 'nw', nlambdas, wavdimid)); wavdim = nlambdas
    END IF
    CALL netcdf_handle_error(location,NF_DEF_DIM (ncid, 'ngeom',  ngeom, geodimid)); geodim = ngeom
    IF (do_brdf_surface) THEN
       CALL netcdf_handle_error(location,NF_DEF_DIM (ncid, 'nstokessq', NSTOKESSQ, nsqdimid)); nsqdim = NSTOKESSQ
    END IF
    
    !=============================================================================
    ! Create the coordinate (aka independent) variables:
    !============================================================================= 
    CALL netcdf_handle_error(location,NF_DEF_VAR (ncid, 'Solarzenithangle',     NF_FLOAT, 1, szadimid, szaid))
    CALL netcdf_handle_error(location,NF_DEF_VAR_FILL(ncid, szaid, FILL_MODE, NF_FILL_REAL))

    CALL netcdf_handle_error(location,NF_DEF_VAR (ncid, 'Viewingzenithangle',   NF_FLOAT, 1, vzadimid, vzaid))
    CALL netcdf_handle_error(location,NF_DEF_VAR_FILL(ncid, vzaid, FILL_MODE, NF_FILL_REAL))

    CALL netcdf_handle_error(location,NF_DEF_VAR (ncid, 'Relativeazimuthangle', NF_FLOAT, 1, azadimid, azaid))
    CALL netcdf_handle_error(location,NF_DEF_VAR_FILL(ncid, azaid, FILL_MODE, NF_FILL_REAL))

    CALL netcdf_handle_error(location,NF_DEF_VAR (ncid, 'gas',                  NF_INT,   1, gasdimid, gasid))
    CALL netcdf_handle_error(location,NF_DEF_VAR_FILL(ncid, gasid, FILL_MODE, NF_FILL_INT))

    CALL netcdf_handle_error(location,NF_DEF_VAR (ncid, 'zs',                   NF_FLOAT, 1, lvldimid, lvlid))
    CALL netcdf_handle_error(location,NF_DEF_VAR_FILL(ncid, lvlid, FILL_MODE, NF_FILL_REAL))

    CALL netcdf_handle_error(location,NF_DEF_VAR (ncid, 'Wavelength',           NF_FLOAT, 1, wavdimid, wavid))
    CALL netcdf_handle_error(location,NF_DEF_VAR_FILL(ncid, wavid, FILL_MODE, NF_FILL_REAL))
    IF (wavdim .LT. members_chunk_4(1)) THEN
       chunk_1d(1) = wavdim
    ELSE
       chunk_1d(1) = members_chunk_4(1)
    END IF
    CALL netcdf_handle_error(location,NF_DEF_VAR_CHUNKING(ncid, wavid,NF_CHUNKED,chunk_1d))

    CALL netcdf_handle_error(location,NF_DEF_VAR (ncid, 'outputlevels',         NF_FLOAT, 1, olvdimid, outlevid))
    CALL netcdf_handle_error(location,NF_DEF_VAR_FILL(ncid, outlevid, FILL_MODE, NF_FILL_REAL))
    
    !=============================================================================
    ! Create a vector containing the often referred to spacetime coordinates:
    !=============================================================================
    gascol_dims(1) = laydimid; gascol_dims(2) = gasdimid
    wavalt_dims(1) = wavdimid; wavalt_dims(2) = laydimid
    wavgeolev_dims(1) = wavdimid; wavgeolev_dims(2) = geodimid; wavgeolev_dims(3) = olvdimid
    wavszalev_dims(1) = wavdimid; wavszalev_dims(2) = szadimid; wavszalev_dims(3) = olvdimid
    
    wavaltgeolev_dims(1) = wavdimid; wavaltgeolev_dims(2) = laydimid
    wavaltgeolev_dims(3) = geodimid; wavaltgeolev_dims(4) = olvdimid

    wavgeogaslev_dims(1) = wavdimid; wavgeogaslev_dims(2) = geodimid 
    wavgeogaslev_dims(3) = gasdimid; wavgeogaslev_dims(4) = olvdimid
    
    wavaltgeogaslev_dims(1) = wavdimid; wavaltgeogaslev_dims(2) = laydimid
    wavaltgeogaslev_dims(3) = geodimid; wavaltgeogaslev_dims(4) = gasdimid; wavaltgeogaslev_dims(5) = olvdimid
    
    brdfdim(1) = nsqdimid; brdfdim(2) = vzadimid; brdfdim(3)=azadimid; brdfdim(4) = szadimid

    !=============================================================================
    ! Create the dependent variables:
    !=============================================================================
    ! ps, ts
    CALL netcdf_handle_error(location,NF_DEF_VAR(ncid, 'ps', NF_FLOAT, 1, lvldimid, psid))
    CALL netcdf_handle_error(location,NF_DEF_VAR_FILL(ncid, psid, FILL_MODE, NF_FILL_REAL))
    CALL netcdf_handle_error(location,NF_DEF_VAR(ncid, 'ts', NF_FLOAT, 1, lvldimid, tsid))
    CALL netcdf_handle_error(location,NF_DEF_VAR_FILL(ncid, tsid, FILL_MODE, NF_FILL_REAL))
    
    ! aircol, aods0, cods0
    CALL netcdf_handle_error(location,NF_DEF_VAR(ncid, 'aircol', NF_FLOAT, 1, laydimid, airid))
    CALL netcdf_handle_error(location,NF_DEF_VAR_FILL(ncid, airid, FILL_MODE, NF_FILL_REAL))
    CALL netcdf_handle_error(location,NF_DEF_VAR(ncid, 'aods0',  NF_FLOAT, 1, laydimid, aer0id))
    CALL netcdf_handle_error(location,NF_DEF_VAR_FILL(ncid, aer0id, FILL_MODE, NF_FILL_REAL))
    CALL netcdf_handle_error(location,NF_DEF_VAR(ncid, 'cods0',  NF_FLOAT, 1, laydimid, cld0id))
    CALL netcdf_handle_error(location,NF_DEF_VAR_FILL(ncid, cld0id, FILL_MODE, NF_FILL_REAL))
        
    ! varaibles with 1D, wavdim
    CALL netcdf_handle_error(location,NF_DEF_VAR(ncid, 'irradiance', NF_FLOAT, 1, wavdimid, irradid))
    CALL netcdf_handle_error(location,NF_DEF_VAR_FILL(ncid, irradid, FILL_MODE, NF_FILL_REAL))
    CALL netcdf_handle_error(location,NF_DEF_VAR_CHUNKING(ncid, irradid,NF_CHUNKED,chunk_1d))
    ! Work out chunking
    CALL netcdf_handle_error(location,NF_DEF_VAR(ncid, 'surfalb',    NF_FLOAT, 1, wavdimid, sfcid))
    CALL netcdf_handle_error(location,NF_DEF_VAR_FILL(ncid, sfcid, FILL_MODE, NF_FILL_REAL))
    CALL netcdf_handle_error(location,NF_DEF_VAR_CHUNKING(ncid, sfcid,NF_CHUNKED,chunk_1d))

    ! variables with 2D, laydim, gasdim
    CALL netcdf_handle_error(location,NF_DEF_VAR(ncid,  'gascol', NF_DOUBLE, 2, gascol_dims, gascolid))
    CALL netcdf_handle_error(location,NF_DEF_VAR_FILL(ncid, gascolid, FILL_MODE, NF_FILL_DOUBLE))
    
    ! variables with 2D, wavdim, laydim
    ! Work out chunking
    max_chunk = FLOOR(REAL(members_chunk_4(1),KIND=4) / REAL(laydim,KIND=4))
    IF (wavdim .LT. max_chunk) THEN
       chunk_2d(1) = wavdim
    ELSE
       chunk_2d(1) = max_chunk
    END IF
    chunk_2d(2) = laydim

    CALL netcdf_handle_error(location,NF_DEF_VAR(ncid, 'ods',   NF_FLOAT, 2, wavalt_dims, odsid))
    CALL netcdf_handle_error(location,NF_DEF_VAR_FILL(ncid, odsid, FILL_MODE, NF_FILL_REAL))
    CALL netcdf_handle_error(location,NF_DEF_VAR_CHUNKING(ncid, odsid,NF_CHUNKED,chunk_2d))
    CALL netcdf_handle_error(location,NF_DEF_VAR(ncid, 'ssas',  NF_FLOAT, 2, wavalt_dims, ssasid))
    CALL netcdf_handle_error(location,NF_DEF_VAR_FILL(ncid, ssasid, FILL_MODE, NF_FILL_REAL))
    CALL netcdf_handle_error(location,NF_DEF_VAR_CHUNKING(ncid, ssasid,NF_CHUNKED,chunk_2d))
    CALL netcdf_handle_error(location,NF_DEF_VAR(ncid, 'aods',  NF_FLOAT, 2, wavalt_dims, aodsid))
    CALL netcdf_handle_error(location,NF_DEF_VAR_FILL(ncid, aodsid, FILL_MODE, NF_FILL_REAL))
    CALL netcdf_handle_error(location,NF_DEF_VAR_CHUNKING(ncid, aodsid,NF_CHUNKED,chunk_2d))
    CALL netcdf_handle_error(location,NF_DEF_VAR(ncid, 'assas', NF_FLOAT, 2, wavalt_dims, assasid))
    CALL netcdf_handle_error(location,NF_DEF_VAR_FILL(ncid, assasid, FILL_MODE, NF_FILL_REAL))
    CALL netcdf_handle_error(location,NF_DEF_VAR_CHUNKING(ncid, assasid,NF_CHUNKED,chunk_2d))
    CALL netcdf_handle_error(location,NF_DEF_VAR(ncid, 'cods',  NF_FLOAT, 2, wavalt_dims, codsid))
    CALL netcdf_handle_error(location,NF_DEF_VAR_FILL(ncid, codsid, FILL_MODE, NF_FILL_REAL))
    CALL netcdf_handle_error(location,NF_DEF_VAR_CHUNKING(ncid, codsid,NF_CHUNKED,chunk_2d))
    CALL netcdf_handle_error(location,NF_DEF_VAR(ncid, 'cssas', NF_FLOAT, 2, wavalt_dims, cssasid))
    CALL netcdf_handle_error(location,NF_DEF_VAR_FILL(ncid, cssasid, FILL_MODE, NF_FILL_REAL))
    CALL netcdf_handle_error(location,NF_DEF_VAR_CHUNKING(ncid, cssasid,NF_CHUNKED,chunk_2d))
    
    ! variables with 3D, wavdim, geodim, olvdim
    ! Work out chunking
    max_chunk = FLOOR(REAL(members_chunk_4(1),KIND=4) / REAL(geodim,KIND=4))
    IF (wavdim .LT. max_chunk) THEN
       chunk_3d(1) = wavdim
    ELSE
       chunk_3d(1) = max_chunk
    END IF
    chunk_3d(2) = geodim; chunk_3d(3) = 1

    CALL netcdf_handle_error(location,NF_DEF_VAR(ncid, 'radiance',    NF_FLOAT, 3, wavgeolev_dims, radid))
    CALL netcdf_handle_error(location,NF_DEF_VAR_FILL(ncid, radid, FILL_MODE, NF_FILL_REAL))
    CALL netcdf_handle_error(location,NF_DEF_VAR_CHUNKING(ncid, radid,NF_CHUNKED,chunk_3d))
    CALL netcdf_handle_error(location,NF_DEF_VAR(ncid, 'flux',        NF_FLOAT, 3, wavszalev_dims, fluxid))
    CALL netcdf_handle_error(location,NF_DEF_VAR_FILL(ncid, fluxid, FILL_MODE, NF_FILL_REAL))
    CALL netcdf_handle_error(location,NF_DEF_VAR_CHUNKING(ncid, fluxid,NF_CHUNKED,chunk_3d))
    CALL netcdf_handle_error(location,NF_DEF_VAR(ncid, 'direct_flux', NF_FLOAT, 3, wavszalev_dims, dfluxid))
    CALL netcdf_handle_error(location,NF_DEF_VAR_FILL(ncid, dfluxid, FILL_MODE, NF_FILL_REAL))
    CALL netcdf_handle_error(location,NF_DEF_VAR_CHUNKING(ncid, dfluxid,NF_CHUNKED,chunk_3d))
    if (do_vector_calculation .and. do_StokesQU_output) then
       CALL netcdf_handle_error(location,NF_DEF_VAR(ncid, 'q',            NF_FLOAT, 3, wavgeolev_dims, qid))
       CALL netcdf_handle_error(location,NF_DEF_VAR_FILL(ncid, qid, FILL_MODE, NF_FILL_REAL))
       CALL netcdf_handle_error(location,NF_DEF_VAR(ncid, 'u',            NF_FLOAT, 3, wavgeolev_dims, uid))
       CALL netcdf_handle_error(location,NF_DEF_VAR_FILL(ncid, uid, FILL_MODE, NF_FILL_REAL))
       CALL netcdf_handle_error(location,NF_DEF_VAR(ncid, 'qflux',        NF_FLOAT, 3, wavszalev_dims, qfluxid))
       CALL netcdf_handle_error(location,NF_DEF_VAR_FILL(ncid, qfluxid, FILL_MODE, NF_FILL_REAL))
       CALL netcdf_handle_error(location,NF_DEF_VAR(ncid, 'uflux',        NF_FLOAT, 3, wavszalev_dims, ufluxid))
       CALL netcdf_handle_error(location,NF_DEF_VAR_FILL(ncid, ufluxid, FILL_MODE, NF_FILL_REAL))
       CALL netcdf_handle_error(location,NF_DEF_VAR(ncid, 'qdirect_flux', NF_FLOAT, 3, wavszalev_dims, qdfluxid))
       CALL netcdf_handle_error(location,NF_DEF_VAR_FILL(ncid, qdfluxid, FILL_MODE, NF_FILL_REAL))
       CALL netcdf_handle_error(location,NF_DEF_VAR(ncid, 'udirect_flux', NF_FLOAT, 3, wavszalev_dims, udfluxid))
       CALL netcdf_handle_error(location,NF_DEF_VAR_FILL(ncid, udfluxid, FILL_MODE, NF_FILL_REAL))
    endif
    
    if (do_Jacobians) then
       
       if (use_lambertian) then
          CALL netcdf_handle_error(location,NF_DEF_VAR(ncid, 'surfalb_jac',   NF_FLOAT, 3, wavgeolev_dims, sfcwfid))
          CALL netcdf_handle_error(location,NF_DEF_VAR_FILL(ncid, sfcwfid, FILL_MODE, NF_FILL_REAL))
          CALL netcdf_handle_error(location,NF_DEF_VAR_CHUNKING(ncid, sfcwfid,NF_CHUNKED,chunk_3d))
       else 
          CALL netcdf_handle_error(location,NF_DEF_VAR(ncid, 'windspeed_jac', NF_FLOAT, 3, wavgeolev_dims, wswfid))
          CALL netcdf_handle_error(location,NF_DEF_VAR_FILL(ncid, wswfid, FILL_MODE, NF_FILL_REAL))
          CALL netcdf_handle_error(location,NF_DEF_VAR_CHUNKING(ncid, wswfid,NF_CHUNKED,chunk_3d))
       endif
       
       if (do_cfrac_Jacobians) then
            CALL netcdf_handle_error(location,NF_DEF_VAR(ncid, 'cfrac_jac',     NF_FLOAT, 3, wavgeolev_dims, cfracwfid))
            CALL netcdf_handle_error(location,NF_DEF_VAR_FILL(ncid, cfracwfid, FILL_MODE, NF_FILL_REAL))
            CALL netcdf_handle_error(location,NF_DEF_VAR_CHUNKING(ncid, cfracwfid,NF_CHUNKED,chunk_3d))
         endif
       if (do_sfcprs_Jacobians) then
            CALL netcdf_handle_error(location,NF_DEF_VAR(ncid, 'sfcprs_jac',    NF_FLOAT, 3, wavgeolev_dims, sfcprswfid))
            CALL netcdf_handle_error(location,NF_DEF_VAR_FILL(ncid, sfcprswfid, FILL_MODE, NF_FILL_REAL))
            CALL netcdf_handle_error(location,NF_DEF_VAR_CHUNKING(ncid, sfcprswfid,NF_CHUNKED,chunk_3d))
         endif
       if (do_aer_columnwf .and. do_aod_Jacobians) then  !No column jacobians yet
          CALL netcdf_handle_error(location,NF_DEF_VAR(ncid, 'aodcolwf_jac',  NF_FLOAT, 3, wavgeolev_dims, aodwfid))
          CALL netcdf_handle_error(location,NF_DEF_VAR_FILL(ncid, aodwfid, FILL_MODE, NF_FILL_REAL))
          CALL netcdf_handle_error(location,NF_DEF_VAR_CHUNKING(ncid, aodwfid,NF_CHUNKED,chunk_3d))
       endif
       if (do_aer_columnwf .and. do_assa_Jacobians) then !No column jacobians yet
          CALL netcdf_handle_error(location,NF_DEF_VAR(ncid, 'assacolwf_jac', NF_FLOAT, 3, wavgeolev_dims, assawfid))
          CALL netcdf_handle_error(location,NF_DEF_VAR_FILL(ncid, assawfid, FILL_MODE, NF_FILL_REAL))
          CALL netcdf_handle_error(location,NF_DEF_VAR_CHUNKING(ncid, assawfid,NF_CHUNKED,chunk_3d))
       endif
       if (do_cld_columnwf .and. do_cod_Jacobians) then  !No column jacobians yet
          CALL netcdf_handle_error(location,NF_DEF_VAR(ncid, 'codcolwf_jac',  NF_FLOAT, 3, wavgeolev_dims, codwfid))
          CALL netcdf_handle_error(location,NF_DEF_VAR_FILL(ncid, codwfid, FILL_MODE, NF_FILL_REAL))
          CALL netcdf_handle_error(location,NF_DEF_VAR_CHUNKING(ncid, codwfid,NF_CHUNKED,chunk_3d))
       endif
       if (do_cld_columnwf .and. do_cssa_Jacobians) then !No column jacobians yet
          CALL netcdf_handle_error(location,NF_DEF_VAR(ncid, 'cssacolwf_jac', NF_FLOAT, 3, wavgeolev_dims, cssawfid))
          CALL netcdf_handle_error(location,NF_DEF_VAR_FILL(ncid, cssawfid, FILL_MODE, NF_FILL_REAL))
          CALL netcdf_handle_error(location,NF_DEF_VAR_CHUNKING(ncid, cssawfid,NF_CHUNKED,chunk_3d))
       endif
       
    endif
    
    if (do_QU_Jacobians) then
       
       if (use_lambertian) then
          CALL netcdf_handle_error(location,NF_DEF_VAR(ncid, 'surfalb_qjac',  NF_FLOAT, 3, wavgeolev_dims, sfcqwfid))
          CALL netcdf_handle_error(location,NF_DEF_VAR_FILL(ncid, sfcqwfid, FILL_MODE, NF_FILL_REAL))
          CALL netcdf_handle_error(location,NF_DEF_VAR_CHUNKING(ncid, sfcqwfid,NF_CHUNKED,chunk_3d))
          CALL netcdf_handle_error(location,NF_DEF_VAR(ncid, 'surfalb_ujac',  NF_FLOAT, 3, wavgeolev_dims, sfcuwfid))
          CALL netcdf_handle_error(location,NF_DEF_VAR_FILL(ncid, sfcuwfid, FILL_MODE, NF_FILL_REAL))
          CALL netcdf_handle_error(location,NF_DEF_VAR_CHUNKING(ncid, sfcuwfid,NF_CHUNKED,chunk_3d))
       else 
          CALL netcdf_handle_error(location,NF_DEF_VAR(ncid, 'windspeed_qjac', NF_FLOAT, 3, wavgeolev_dims, wsqwfid))
          CALL netcdf_handle_error(location,NF_DEF_VAR_FILL(ncid, wsqwfid, FILL_MODE, NF_FILL_REAL))
          CALL netcdf_handle_error(location,NF_DEF_VAR_CHUNKING(ncid, wsqwfid,NF_CHUNKED,chunk_3d))
          CALL netcdf_handle_error(location,NF_DEF_VAR(ncid, 'windspeed_ujac', NF_FLOAT, 3, wavgeolev_dims, wsuwfid))
          CALL netcdf_handle_error(location,NF_DEF_VAR_FILL(ncid, wsuwfid, FILL_MODE, NF_FILL_REAL))
          CALL netcdf_handle_error(location,NF_DEF_VAR_CHUNKING(ncid, wsuwfid,NF_CHUNKED,chunk_3d))
       endif
       
       if (do_cfrac_Jacobians) then
          CALL netcdf_handle_error(location,NF_DEF_VAR(ncid, 'cfrac_qjac',  NF_FLOAT, 3, wavgeolev_dims, cfracqwfid))
          CALL netcdf_handle_error(location,NF_DEF_VAR_FILL(ncid, cfracqwfid, FILL_MODE, NF_FILL_REAL))
          CALL netcdf_handle_error(location,NF_DEF_VAR_CHUNKING(ncid, cfracqwfid,NF_CHUNKED,chunk_3d))
          CALL netcdf_handle_error(location,NF_DEF_VAR(ncid, 'cfrac_ujac',  NF_FLOAT, 3, wavgeolev_dims, cfracuwfid))
          CALL netcdf_handle_error(location,NF_DEF_VAR_FILL(ncid, cfracuwfid, FILL_MODE, NF_FILL_REAL))
          CALL netcdf_handle_error(location,NF_DEF_VAR_CHUNKING(ncid, cfracuwfid,NF_CHUNKED,chunk_3d))
       endif
       if (do_sfcprs_Jacobians) then
          CALL netcdf_handle_error(location,NF_DEF_VAR(ncid, 'sfcprs_qjac',  NF_FLOAT, 3, wavgeolev_dims, sfcprsqwfid))
          CALL netcdf_handle_error(location,NF_DEF_VAR_FILL(ncid, sfcprsqwfid, FILL_MODE, NF_FILL_REAL))
          CALL netcdf_handle_error(location,NF_DEF_VAR_CHUNKING(ncid, sfcprsqwfid,NF_CHUNKED,chunk_3d))
          CALL netcdf_handle_error(location,NF_DEF_VAR(ncid, 'sfcprs_ujac',  NF_FLOAT, 3, wavgeolev_dims, sfcprsuwfid))
          CALL netcdf_handle_error(location,NF_DEF_VAR_FILL(ncid, sfcprsuwfid, FILL_MODE, NF_FILL_REAL))
          CALL netcdf_handle_error(location,NF_DEF_VAR_CHUNKING(ncid, sfcprsuwfid,NF_CHUNKED,chunk_3d))
       endif
       if (do_aer_columnwf .and. do_aod_Jacobians) then  !No column jacobians yet
          CALL netcdf_handle_error(location,NF_DEF_VAR(ncid, 'aodcol_qjac',  NF_FLOAT, 3, wavgeolev_dims, aodqwfid))
          CALL netcdf_handle_error(location,NF_DEF_VAR_FILL(ncid, aodqwfid, FILL_MODE, NF_FILL_REAL))
          CALL netcdf_handle_error(location,NF_DEF_VAR_CHUNKING(ncid, aodqwfid,NF_CHUNKED,chunk_3d))
          CALL netcdf_handle_error(location,NF_DEF_VAR(ncid, 'aodcol_ujac',  NF_FLOAT, 3, wavgeolev_dims, aoduwfid))
          CALL netcdf_handle_error(location,NF_DEF_VAR_FILL(ncid, aoduwfid, FILL_MODE, NF_FILL_REAL))
          CALL netcdf_handle_error(location,NF_DEF_VAR_CHUNKING(ncid, aoduwfid,NF_CHUNKED,chunk_3d))
       endif
       if (do_aer_columnwf .and. do_assa_Jacobians) then !No column jacobians yet
          CALL netcdf_handle_error(location,NF_DEF_VAR(ncid, 'assacol_qjac', NF_FLOAT, 3, wavgeolev_dims, assaqwfid))
          CALL netcdf_handle_error(location,NF_DEF_VAR_FILL(ncid, assaqwfid, FILL_MODE, NF_FILL_REAL))
          CALL netcdf_handle_error(location,NF_DEF_VAR_CHUNKING(ncid, assaqwfid,NF_CHUNKED,chunk_3d))
          CALL netcdf_handle_error(location,NF_DEF_VAR(ncid, 'assacol_ujac', NF_FLOAT, 3, wavgeolev_dims, assauwfid))
          CALL netcdf_handle_error(location,NF_DEF_VAR_FILL(ncid, assauwfid, FILL_MODE, NF_FILL_REAL))
          CALL netcdf_handle_error(location,NF_DEF_VAR_CHUNKING(ncid, assauwfid,NF_CHUNKED,chunk_3d))
       endif
       if (do_cld_columnwf .and. do_cod_Jacobians) then  !No column jacobians yet
          CALL netcdf_handle_error(location,NF_DEF_VAR(ncid, 'codcol_qjac',  NF_FLOAT, 3, wavgeolev_dims, codqwfid))
          CALL netcdf_handle_error(location,NF_DEF_VAR_FILL(ncid, codqwfid, FILL_MODE, NF_FILL_REAL))
          CALL netcdf_handle_error(location,NF_DEF_VAR_CHUNKING(ncid, codqwfid,NF_CHUNKED,chunk_3d))
          CALL netcdf_handle_error(location,NF_DEF_VAR(ncid, 'codcol_ujac',  NF_FLOAT, 3, wavgeolev_dims, coduwfid))
          CALL netcdf_handle_error(location,NF_DEF_VAR_FILL(ncid, coduwfid, FILL_MODE, NF_FILL_REAL))
          CALL netcdf_handle_error(location,NF_DEF_VAR_CHUNKING(ncid, coduwfid,NF_CHUNKED,chunk_3d))
       endif
       if (do_cld_columnwf .and. do_cssa_Jacobians) then !No column jacobians yet
          CALL netcdf_handle_error(location,NF_DEF_VAR(ncid, 'cssacol_qjac', NF_FLOAT, 3, wavgeolev_dims, cssaqwfid))
          CALL netcdf_handle_error(location,NF_DEF_VAR_FILL(ncid, cssaqwfid, FILL_MODE, NF_FILL_REAL))
          CALL netcdf_handle_error(location,NF_DEF_VAR_CHUNKING(ncid, cssaqwfid,NF_CHUNKED,chunk_3d))
          CALL netcdf_handle_error(location,NF_DEF_VAR(ncid, 'cssacol_ujac', NF_FLOAT, 3, wavgeolev_dims, cssauwfid))
          CALL netcdf_handle_error(location,NF_DEF_VAR_FILL(ncid, cssauwfid, FILL_MODE, NF_FILL_REAL))
          CALL netcdf_handle_error(location,NF_DEF_VAR_CHUNKING(ncid, cssauwfid,NF_CHUNKED,chunk_3d))
       endif       
    endif
    
    ! variables with 4D, wavdim, laydim, geodim, olvdim
    ! Work out chunking
    max_chunk = FLOOR(REAL(members_chunk_4(1),KIND=4) / REAL(geodim*laydim,KIND=4))
    IF (wavdim .LT. max_chunk) THEN
       chunk_4d(1) = wavdim
    ELSE
       chunk_4d(1) = max_chunk
    END IF
    chunk_4d(2) = laydim; chunk_4d(3) = geodim; chunk_4d(4) = 1

    if (do_AMF_calculation) then
       CALL netcdf_handle_error(location,NF_DEF_VAR(ncid, 'scatweights', NF_FLOAT, 4, wavaltgeolev_dims, scatwtid))
       CALL netcdf_handle_error(location,NF_DEF_VAR_FILL(ncid, scatwtid, FILL_MODE, NF_FILL_REAL))
       CALL netcdf_handle_error(location,NF_DEF_VAR_CHUNKING(ncid, scatwtid,NF_CHUNKED,chunk_4d))
    endif
    
    if (do_Jacobians) then
       if (do_T_Jacobians) then
            CALL netcdf_handle_error(location,NF_DEF_VAR(ncid, 't_jac',    NF_FLOAT, 4, wavaltgeolev_dims, tempwfid))
            CALL netcdf_handle_error(location,NF_DEF_VAR_FILL(ncid, tempwfid, FILL_MODE, NF_FILL_REAL))
            CALL netcdf_handle_error(location,NF_DEF_VAR_CHUNKING(ncid, tempwfid,NF_CHUNKED,chunk_4d))
         endif
       if (.not. do_aer_columnwf .and. do_aod_Jacobians) then
            CALL netcdf_handle_error(location,NF_DEF_VAR(ncid, 'aod_jac',  NF_FLOAT, 4, wavaltgeolev_dims, aodwfid))
            CALL netcdf_handle_error(location,NF_DEF_VAR_FILL(ncid, aodwfid, FILL_MODE, NF_FILL_REAL))
            CALL netcdf_handle_error(location,NF_DEF_VAR_CHUNKING(ncid, aodwfid,NF_CHUNKED,chunk_4d))
         endif
       if (.not. do_aer_columnwf .and. do_assa_Jacobians) then
            CALL netcdf_handle_error(location,NF_DEF_VAR(ncid, 'assa_jac', NF_FLOAT, 4, wavaltgeolev_dims, assawfid))
            CALL netcdf_handle_error(location,NF_DEF_VAR_FILL(ncid, assawfid, FILL_MODE, NF_FILL_REAL))
            CALL netcdf_handle_error(location,NF_DEF_VAR_CHUNKING(ncid, assawfid,NF_CHUNKED,chunk_4d))
         endif
       if (.not. do_cld_columnwf .and. do_cod_Jacobians) then
            CALL netcdf_handle_error(location,NF_DEF_VAR(ncid, 'cod_jac',  NF_FLOAT, 4, wavaltgeolev_dims, codwfid))
            CALL netcdf_handle_error(location,NF_DEF_VAR_FILL(ncid, codwfid, FILL_MODE, NF_FILL_REAL))
            CALL netcdf_handle_error(location,NF_DEF_VAR_CHUNKING(ncid, codwfid,NF_CHUNKED,chunk_4d))
         endif
       if (.not. do_cld_columnwf .and. do_cssa_Jacobians) then
            CALL netcdf_handle_error(location,NF_DEF_VAR(ncid, 'cssa_jac', NF_FLOAT, 4, wavaltgeolev_dims, cssawfid))
            CALL netcdf_handle_error(location,NF_DEF_VAR_FILL(ncid, cssawfid, FILL_MODE, NF_FILL_REAL))
            CALL netcdf_handle_error(location,NF_DEF_VAR_CHUNKING(ncid, cssawfid,NF_CHUNKED,chunk_4d))
         endif
    endif
    
    if (do_QU_Jacobians) then
       if (.not. do_aer_columnwf .and. do_aod_Jacobians) then 
          CALL netcdf_handle_error(location,NF_DEF_VAR(ncid, 'aod_qjac', NF_FLOAT, 4, wavaltgeolev_dims, aodqwfid))
          CALL netcdf_handle_error(location,NF_DEF_VAR_FILL(ncid, aodqwfid, FILL_MODE, NF_FILL_REAL))
          CALL netcdf_handle_error(location,NF_DEF_VAR_CHUNKING(ncid, aodqwfid,NF_CHUNKED,chunk_4d))
          CALL netcdf_handle_error(location,NF_DEF_VAR(ncid, 'aod_ujac', NF_FLOAT, 4, wavaltgeolev_dims, aoduwfid))
          CALL netcdf_handle_error(location,NF_DEF_VAR_FILL(ncid, aoduwfid, FILL_MODE, NF_FILL_REAL))
          CALL netcdf_handle_error(location,NF_DEF_VAR_CHUNKING(ncid, aoduwfid,NF_CHUNKED,chunk_4d))
       endif
       if (.not. do_aer_columnwf .and. do_assa_Jacobians) then
          CALL netcdf_handle_error(location,NF_DEF_VAR(ncid, 'assa_qjac', NF_FLOAT, 4, wavaltgeolev_dims, assaqwfid))
          CALL netcdf_handle_error(location,NF_DEF_VAR_FILL(ncid, assaqwfid, FILL_MODE, NF_FILL_REAL))
          CALL netcdf_handle_error(location,NF_DEF_VAR_CHUNKING(ncid, assaqwfid,NF_CHUNKED,chunk_4d))
          CALL netcdf_handle_error(location,NF_DEF_VAR(ncid, 'assa_ujac', NF_FLOAT, 4, wavaltgeolev_dims, assauwfid))
          CALL netcdf_handle_error(location,NF_DEF_VAR_FILL(ncid, assauwfid, FILL_MODE, NF_FILL_REAL))
          CALL netcdf_handle_error(location,NF_DEF_VAR_CHUNKING(ncid, assauwfid,NF_CHUNKED,chunk_4d))
       endif
       if (.not. do_cld_columnwf .and. do_cod_Jacobians) then
          CALL netcdf_handle_error(location,NF_DEF_VAR(ncid, 'cod_qjac', NF_FLOAT, 4, wavaltgeolev_dims, codqwfid))
          CALL netcdf_handle_error(location,NF_DEF_VAR_FILL(ncid, codqwfid, FILL_MODE, NF_FILL_REAL))
          CALL netcdf_handle_error(location,NF_DEF_VAR_CHUNKING(ncid, codqwfid,NF_CHUNKED,chunk_4d))
          CALL netcdf_handle_error(location,NF_DEF_VAR(ncid, 'cod_ujac', NF_FLOAT, 4, wavaltgeolev_dims, coduwfid))
          CALL netcdf_handle_error(location,NF_DEF_VAR_FILL(ncid, coduwfid, FILL_MODE, NF_FILL_REAL))
          CALL netcdf_handle_error(location,NF_DEF_VAR_CHUNKING(ncid, coduwfid,NF_CHUNKED,chunk_4d))
       endif
       if (.not. do_cld_columnwf .and. do_cssa_Jacobians) then
          CALL netcdf_handle_error(location,NF_DEF_VAR(ncid, 'cssa_qjac', NF_FLOAT, 4, wavaltgeolev_dims, cssaqwfid))
          CALL netcdf_handle_error(location,NF_DEF_VAR_FILL(ncid, cssaqwfid, FILL_MODE, NF_FILL_REAL))
          CALL netcdf_handle_error(location,NF_DEF_VAR_CHUNKING(ncid, cssaqwfid,NF_CHUNKED,chunk_4d))
          CALL netcdf_handle_error(location,NF_DEF_VAR(ncid, 'cssa_ujac', NF_FLOAT, 4, wavaltgeolev_dims, cssauwfid))
          CALL netcdf_handle_error(location,NF_DEF_VAR_FILL(ncid, cssauwfid, FILL_MODE, NF_FILL_REAL))
          CALL netcdf_handle_error(location,NF_DEF_VAR_CHUNKING(ncid, cssauwfid,NF_CHUNKED,chunk_4d))
       endif
    endif
    
    ! variables with 4D, wavdim, geodim, gasdim, olvdim
    max_chunk = FLOOR(REAL(members_chunk_4(1),KIND=4) / REAL(geodim,KIND=4))
    IF (wavdim .LT. max_chunk) THEN
       chunk_4d(1) = wavdim
    ELSE
       chunk_4d(1) = max_chunk
    END IF
    chunk_4d(2) = geodim; chunk_4d(3) = 1; chunk_4d(4) = 1

    if (do_AMF_calculation) then
       CALL netcdf_handle_error(location,NF_DEF_VAR(ncid, 'amf', NF_FLOAT, 4, wavgeogaslev_dims, amfid))
       CALL netcdf_handle_error(location,NF_DEF_VAR_FILL(ncid, amfid, FILL_MODE, NF_FILL_REAL))
       CALL netcdf_handle_error(location,NF_DEF_VAR_CHUNKING(ncid, amfid,NF_CHUNKED,chunk_4d))
    endif
    
    ! variables with 5D, wavdim, laydim, geodim, gasdim, olvdim
    max_chunk = FLOOR(REAL(members_chunk_4(1),KIND=4) / (REAL(laydim,KIND=4)*REAL(geodim,KIND=4)) )
    IF (wavdim .LT. max_chunk) THEN
       chunk_5d(1) = wavdim
    ELSE
       chunk_5d(1) = max_chunk
    END IF
    chunk_5d(2) = laydim; chunk_5d(3) = geodim; chunk_5d(4) = 1; chunk_5d(5) = 1

    if (do_Jacobians) then
       CALL netcdf_handle_error(location,NF_DEF_VAR(ncid,  'gas_jac', NF_FLOAT, 5, wavaltgeogaslev_dims, gaswfid))
       CALL netcdf_handle_error(location,NF_DEF_VAR_FILL(ncid, gaswfid, FILL_MODE, NF_FILL_REAL))
       CALL netcdf_handle_error(location,NF_DEF_VAR_CHUNKING(ncid, gaswfid,NF_CHUNKED,chunk_5d))
       if (do_QU_Jacobians) then
          CALL netcdf_handle_error(location,NF_DEF_VAR(ncid,  'gas_qjac', NF_FLOAT, 5, wavaltgeogaslev_dims, gasqwfid))
          CALL netcdf_handle_error(location,NF_DEF_VAR_FILL(ncid, gasqwfid, FILL_MODE, NF_FILL_REAL))
          CALL netcdf_handle_error(location,NF_DEF_VAR_CHUNKING(ncid, gasqwfid,NF_CHUNKED,chunk_5d))
          CALL netcdf_handle_error(location,NF_DEF_VAR(ncid,  'gas_ujac', NF_FLOAT, 5, wavaltgeogaslev_dims, gasuwfid))
          CALL netcdf_handle_error(location,NF_DEF_VAR_FILL(ncid, gasuwfid, FILL_MODE, NF_FILL_REAL))
          CALL netcdf_handle_error(location,NF_DEF_VAR_CHUNKING(ncid, gasuwfid,NF_CHUNKED,chunk_5d))
       endif
    endif
    
    ! BRDF variables onedim, nsqdim, nvza, naza, nsza
    if (do_brdf_surface) then
    ! WSA, BSA amd BRDF
       CALL netcdf_handle_error(location,NF_DEF_VAR(ncid, 'WSA', NF_FLOAT, 1, onedim, wsaid))
       CALL netcdf_handle_error(location,NF_DEF_VAR_FILL(ncid, wsaid, FILL_MODE, NF_FILL_REAL))
       CALL netcdf_handle_error(location,NF_DEF_VAR(ncid, 'BSA', NF_FLOAT, 1, onedim, bsaid))
       CALL netcdf_handle_error(location,NF_DEF_VAR_FILL(ncid, bsaid, FILL_MODE, NF_FILL_REAL))
       CALL netcdf_handle_error(location,NF_DEF_VAR(ncid, 'BRDF', NF_FLOAT, 4, brdfdim, brdfid))
       CALL netcdf_handle_error(location,NF_DEF_VAR_FILL(ncid, brdfid, FILL_MODE, NF_FILL_REAL))
    endif

    !=============================================================================
    ! Assign attributes (meta-data) to the various variables:
    !============================================================================= 
    CALL netcdf_handle_error(location,NF_PUT_ATT_REAL (ncid, NF_GLOBAL, 'lon', NF_FLOAT, 1, real(longitude, kind=4)))
    CALL netcdf_handle_error(location,NF_PUT_ATT_REAL (ncid, NF_GLOBAL, 'lat', NF_FLOAT, 1, real(latitude, kind=4)))
    CALL netcdf_handle_error(location,NF_PUT_ATT_INT (ncid, NF_GLOBAL, 'month', NF_INT,  1, month))
    CALL netcdf_handle_error(location,NF_PUT_ATT_INT (ncid, NF_GLOBAL, 'year', NF_INT,  1, year))
    CALL netcdf_handle_error(location,NF_PUT_ATT_INT (ncid, NF_GLOBAL, 'day', NF_INT,  1, day))
    CALL netcdf_handle_error(location,NF_PUT_ATT_REAL (ncid, NF_GLOBAL, 'utc', NF_FLOAT, 1, real(utc, kind=4)))
    CALL netcdf_handle_error(location,NF_PUT_ATT_REAL (ncid, NF_GLOBAL, 'cfrac', NF_FLOAT, 1, real(cfrac, kind=4)))
    CALL netcdf_handle_error(location,NF_PUT_ATT_REAL (ncid, NF_GLOBAL, 'zcldtop', NF_FLOAT, 1, real(cld_tops, kind=4)))
    CALL netcdf_handle_error(location,NF_PUT_ATT_REAL (ncid, NF_GLOBAL, 'cldalb', NF_FLOAT, 1, real(lambertian_cldalb, kind=4)))
    CALL netcdf_handle_error(location,NF_PUT_ATT_REAL (ncid, NF_GLOBAL, 'userlvl', NF_FLOAT, 1, &
         real(VLIDORT_ModIn%MUserVal%TS_USER_LEVELS, kind=4)))
    CALL netcdf_handle_error(location,NF_PUT_ATT_REAL (ncid, NF_GLOBAL, 'windspeed', NF_FLOAT, 1, real(wind_speed, kind=4)))
    CALL netcdf_handle_error(location,NF_PUT_ATT_REAL (ncid, NF_GLOBAL, 'lam_reso', NF_FLOAT, 1, real(lambda_resolution, kind=4)))
    CALL netcdf_handle_error(location,NF_PUT_ATT_REAL (ncid, NF_GLOBAL, 'aer_reflam', NF_FLOAT, 1, real(aer_reflambda, kind=4)))
    CALL netcdf_handle_error(location,NF_PUT_ATT_REAL (ncid, NF_GLOBAL, 'cld_reflam', NF_FLOAT, 1, real(cld_reflambda, kind=4)))
    CALL netcdf_handle_error(location,NF_PUT_ATT_REAL (ncid, NF_GLOBAL, 'lam_dw', NF_FLOAT, 1, real(lambda_dw, kind=4)))
    CALL netcdf_handle_error(location,NF_PUT_ATT_REAL (ncid, NF_GLOBAL, 'lam_dfw', NF_FLOAT, 1, real(lambda_dfw, kind=4)))
    CALL netcdf_handle_error(location,NF_PUT_ATT_INT (ncid, NF_GLOBAL, 'nstreams', NF_INT,  1, VLIDORT_FixIn%Cont%TS_NSTREAMS))
    CALL netcdf_handle_error(location,NF_PUT_ATT_INT (ncid, NF_GLOBAL, 'aercld_nmoms', NF_INT,  1, &
         VLIDORT_ModIn%MCont%TS_NGREEK_MOMENTS_INPUT))
    CALL netcdf_handle_error(location,NF_PUT_ATT_INT (ncid, NF_GLOBAL, 'nsza', NF_INT,  1, GC_n_sun_positions))
    CALL netcdf_handle_error(location,NF_PUT_ATT_INT (ncid, NF_GLOBAL, 'nvza', NF_INT,  1, GC_n_view_angles))
    CALL netcdf_handle_error(location,NF_PUT_ATT_INT (ncid, NF_GLOBAL, 'naza', NF_INT,  1, GC_n_azimuths))
    CALL netcdf_handle_error(location,NF_PUT_ATT_INT (ncid, NF_GLOBAL, 'ngeometries', NF_INT,  1, ngeom))
    CALL netcdf_handle_error(location,NF_PUT_ATT_INT (ncid, NF_GLOBAL, 'noutputlevels',NF_INT,  1, &
         VLIDORT_FixIn%UserVal%TS_N_USER_LEVELS))

    IF ( .NOT. do_effcrs .AND. lambda_resolution /= 0.0d0 ) THEN
       CALL netcdf_handle_error(location,NF_PUT_ATT_INT (ncid, NF_GLOBAL, 'nwavelengths', NF_INT,  1, nclambdas))
    ELSE
       CALL netcdf_handle_error(location,NF_PUT_ATT_INT (ncid, NF_GLOBAL, 'nwavelengths', NF_INT,  1, nlambdas))
    END IF
    CALL netcdf_handle_error(location,NF_PUT_ATT_INT (ncid, NF_GLOBAL, 'nlayers', NF_INT,  1, GC_nlayers))
    CALL netcdf_handle_error(location,NF_PUT_ATT_INT (ncid, NF_GLOBAL, 'ngas', NF_INT,  1, ngases))
    
    if (use_lambertian) then
       CALL netcdf_handle_error(location,NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'user_lambertian', NF_INT1, 1, 1))
    else
       CALL netcdf_handle_error(location,NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'user_lambertian', NF_INT1, 1, 0))
    endif
    if (do_lambertian_cld) then
       CALL netcdf_handle_error(location,NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'do_lambertian_cld', NF_INT1, 1, 1))
    else
       CALL netcdf_handle_error(location,NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'do_lambertian_cld', NF_INT1, 1, 0))
    endif
    if (do_effcrs) then
       CALL netcdf_handle_error(location,NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'do_effcrs', NF_INT1, 1, 1))
    else
       CALL netcdf_handle_error(location,NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'do_effcrs', NF_INT1, 1, 0))
    endif
    if (use_wavelength) then
       CALL netcdf_handle_error(location,NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'use_wavelength', NF_INT1, 1, 1))
    else
       CALL netcdf_handle_error(location,NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'use_wavelength', NF_INT1, 1, 0))
    endif
    if (VLIDORT_FixIn%Bool%TS_DO_UPWELLING) then
       CALL netcdf_handle_error(location,NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'do_upwelling', NF_INT1, 1, 1))
    else
       CALL netcdf_handle_error(location,NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'do_upwelling', NF_INT1, 1, 0))
    endif
    if (do_normalized_WFoutput) then
       CALL netcdf_handle_error(location,NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'do_norm_WFoutput', NF_INT1, 1, 1))
    else
       CALL netcdf_handle_error(location,NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'do_norm_WFoutput', NF_INT1, 1, 0))
    endif
    if (do_normalized_radiance) then
       CALL netcdf_handle_error(location,NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'do_norm_radiance', NF_INT1, 1, 1))
    else
       CALL netcdf_handle_error(location,NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'do_norm_radiance', NF_INT1, 1, 0))
    endif
    if (use_solar_photons) then
       CALL netcdf_handle_error(location,NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'use_solar_photons', NF_INT1, 1, 1))
    else
       CALL netcdf_handle_error(location,NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'use_solar_photons', NF_INT1, 1, 0))
    endif
    if (do_vector_calculation) then
       CALL netcdf_handle_error(location,NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'do_vector_calc', NF_INT1, 1, 1))
    else
       CALL netcdf_handle_error(location,NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'do_vector_calc', NF_INT1, 1, 0))
    endif
    if (do_StokesQU_output) then
       CALL netcdf_handle_error(location,NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'do_QU_output', NF_INT1, 1, 1))
    else
       CALL netcdf_handle_error(location,NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'do_QU_output', NF_INT1, 1, 0))
    endif
    if (do_Jacobians) then
       CALL netcdf_handle_error(location,NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'do_Jacobian', NF_INT1, 1, 1))
    else
       CALL netcdf_handle_error(location,NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'do_Jacobian', NF_INT1, 1, 0))
    endif
    if (do_QU_Jacobians) then
       CALL netcdf_handle_error(location,NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'do_QU_Jacobian', NF_INT1, 1, 1))
    else
       CALL netcdf_handle_error(location,NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'do_QU_Jacobian', NF_INT1, 1, 0))
    endif
    if (do_AMF_calculation) then
       CALL netcdf_handle_error(location,NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'do_AMF_calc', NF_INT1, 1, 1))
    else
       CALL netcdf_handle_error(location,NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'do_AMF_calc', NF_INT1, 1, 0))
    endif
    if (do_T_Jacobians) then
       CALL netcdf_handle_error(location,NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'do_T_Jacobian', NF_INT1, 1, 1))
    else
       CALL netcdf_handle_error(location,NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'do_T_Jacobian', NF_INT1, 1, 0))
    endif
    if (do_sfcprs_Jacobians) then
       CALL netcdf_handle_error(location,NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'do_sfcprs_Jacobian', NF_INT1, 1, 1))
    else
       CALL netcdf_handle_error(location,NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'do_sfcprs_Jacobian', NF_INT1, 1, 0))
    endif
    if (do_aod_Jacobians) then
       CALL netcdf_handle_error(location,NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'do_aod_Jacobian', NF_INT1, 1, 1))
    else
       CALL netcdf_handle_error(location,NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'do_aod_Jacobian', NF_INT1, 1, 0))
    endif
    if (do_assa_Jacobians) then
       CALL netcdf_handle_error(location,NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'do_assa_Jacobian', NF_INT1, 1, 1))
    else
       CALL netcdf_handle_error(location,NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'do_assa_Jacobian', NF_INT1, 1, 0))
    endif
    if (do_cod_Jacobians) then
       CALL netcdf_handle_error(location,NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'do_cod_Jacobian', NF_INT1, 1, 1))
    else
       CALL netcdf_handle_error(location,NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'do_cod_Jacobian', NF_INT1, 1, 0))
    endif
    if (do_cssa_Jacobians) then
       CALL netcdf_handle_error(location,NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'do_cssa_Jacobian', NF_INT1, 1, 1))
    else
       CALL netcdf_handle_error(location,NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'do_cssa_Jacobian', NF_INT1, 1, 0))
    endif
    if (do_cfrac_Jacobians) then
       CALL netcdf_handle_error(location,NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'do_cfrac_Jacobian', NF_INT1, 1, 1))
    else
       CALL netcdf_handle_error(location,NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'do_cfrac_Jacobian', NF_INT1, 1, 0))
    endif
    if (do_aer_columnwf) then
       CALL netcdf_handle_error(location,NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'do_aer_columnwf', NF_INT1, 1, 1))
    else
       CALL netcdf_handle_error(location,NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'do_aer_columnwf', NF_INT1, 1, 0))
    endif
    if (do_cld_columnwf) then
       CALL netcdf_handle_error(location,NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'do_cld_columnwf', NF_INT1, 1, 1))
    else
       CALL netcdf_handle_error(location,NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'do_cld_columnwf', NF_INT1, 1, 0))
    endif
    if (do_brdf_surface) then
       CALL netcdf_handle_error(location,NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'do_brdf', NF_INT1, 1, 1))
    else
       CALL netcdf_handle_error(location,NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'do_brdf', NF_INT1, 1, 0))
    endif
    if (OUTPUT_WSABSA) then
       CALL netcdf_handle_error(location,NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'do_output_wsabsa', NF_INT1, 1, 1))
    else
       CALL netcdf_handle_error(location,NF_PUT_ATT_INT1 (ncid, NF_GLOBAL, 'do_output_wsabsa', NF_INT1, 1, 0))
    endif
    
    ! write the list of gases in one string
    write(gasstr, '(10(I2,A1,A4,A1))') (i,':',which_gases(i),',', i=1, ngases)
    nlen=LEN(trim(gasstr)) ; gasstr=gasstr(1:nlen-1)
    CALL netcdf_handle_error(location,NF_PUT_ATT_TEXT (ncid, NF_GLOBAL, 'gases', nlen, gasstr))

    ! Units for variables
    CALL netcdf_handle_error(location,NF_PUT_ATT_TEXT (ncid, NF_GLOBAL, 'windspeed_units', 3, 'm/s'))
    CALL netcdf_handle_error(location,NF_PUT_ATT_TEXT (ncid, NF_GLOBAL, 'lonlat_units',    7, 'degrees'))
    CALL netcdf_handle_error(location,NF_PUT_ATT_TEXT (ncid, NF_GLOBAL, 'zcldtop_units',   2, 'km'))
    CALL netcdf_handle_error(location,NF_PUT_ATT_TEXT (ncid, szaid,    'units', 7, 'degrees'))
    CALL netcdf_handle_error(location,NF_PUT_ATT_TEXT (ncid, vzaid,    'units', 7, 'degrees'))
    CALL netcdf_handle_error(location,NF_PUT_ATT_TEXT (ncid, azaid,    'units', 7, 'degrees'))
    CALL netcdf_handle_error(location,NF_PUT_ATT_TEXT (ncid, lvlid,    'units', 2, 'km'))
    CALL netcdf_handle_error(location,NF_PUT_ATT_TEXT (ncid, outlevid, 'units', 8, 'unitless'))
    if (use_wavelength) then
       CALL netcdf_handle_error(location,NF_PUT_ATT_TEXT(ncid,wavid,'units', 2, 'nm'))
    else 
       CALL netcdf_handle_error(location,NF_PUT_ATT_TEXT(ncid,wavid,'units', 4, 'cm-1'))
    endif
    if (use_solar_photons) then
       irrad_unitc='photons/cm2/nm/s'
    else
       irrad_unitc='W/m2/cm-1'
    endif
    CALL netcdf_handle_error(location,NF_PUT_ATT_TEXT(ncid, irradid, 'units', len_trim(irrad_unitc), irrad_unitc))
    if (do_normalized_radiance) then
       rad_unitc='unitless'
    else
       rad_unitc=TRIM(irrad_unitc)//'/sr'
    endif
    CALL netcdf_handle_error(location,NF_PUT_ATT_TEXT(ncid, radid,    'units', len_trim(rad_unitc), rad_unitc))
    CALL netcdf_handle_error(location,NF_PUT_ATT_TEXT(ncid, gascolid, 'units', 13, 'molecules/cm2'))
    CALL netcdf_handle_error(location,NF_PUT_ATT_TEXT(ncid, airid,    'units', 13, 'molecules/cm2'))
    CALL netcdf_handle_error(location,NF_PUT_ATT_TEXT(ncid, psid,     'units', 3,  'hPa'))
    CALL netcdf_handle_error(location,NF_PUT_ATT_TEXT(ncid, tsid,     'units', 1,  'K')) 

    !=============================================================================
    ! Get out of 'define' mode, and into 'data' mode
    !=============================================================================
    CALL netcdf_handle_error(location,NF_ENDDEF (ncid))

    !=============================================================================
    ! Fill in the values of the non time varying variables.
    ! Remember, we're still in "initialization" here - this is only done the
    ! first time through.
    !=============================================================================
    CALL netcdf_handle_error(location,NF_PUT_VARA_REAL (ncid, szaid, 1, GC_n_sun_positions, &
         real(VLIDORT_ModIn%MSunrays%TS_SZANGLES(1:GC_n_sun_positions), kind=4)))
    CALL netcdf_handle_error(location,NF_PUT_VARA_REAL (ncid, vzaid, 1, GC_n_view_angles, &
         real(VLIDORT_ModIn%MUserVal%TS_USER_VZANGLES_INPUT(1:GC_n_view_angles), kind=4)))
    CALL netcdf_handle_error(location,NF_PUT_VARA_REAL (ncid, azaid, 1, GC_n_azimuths, &
         real(VLIDORT_ModIn%MUserVal%TS_USER_RELAZMS(1:GC_n_azimuths), kind=4)))
    CALL netcdf_handle_error(location,NF_PUT_VARA_REAL (ncid, lvlid, 1, GC_nlayers+1, &
         real(heights(0:GC_nlayers), kind=4)))
    IF ( .NOT. do_effcrs .AND. lambda_resolution /= 0.0d0) THEN
       CALL netcdf_handle_error(location,NF_PUT_VARA_REAL (ncid, wavid, 1, nclambdas, &
            real(clambdas(1:nclambdas), kind=4)))
    ELSE
       CALL netcdf_handle_error(location,NF_PUT_VARA_REAL (ncid, wavid, 1, nlambdas, &
            real(lambdas(1:nlambdas), kind=4)))
    END IF
    
    do i = 1, ngases
       gas_indices(i) = i
    enddo
    CALL netcdf_handle_error(location,NF_PUT_VARA_INT (ncid, gasid, 1, ngases, gas_indices))
    CALL netcdf_handle_error(location,NF_PUT_VARA_REAL (ncid, outlevid, 1, VLIDORT_FixIn%UserVal%TS_N_USER_LEVELS, &
         real(VLIDORT_ModIn%MUserVal%TS_USER_LEVELS(1:VLIDORT_FixIn%UserVal%TS_N_USER_LEVELS), kind=4)))

    !=============================================================================
    ! Define the START and COUNT arrays for each of the array variables.
    ! Fill in the values for other variables
    !=============================================================================
    ! write ps and ts
    ndimstart1 = (/ 1 /)
    ndimcount1 = (/ GC_nlayers+1 /)  
    CALL netcdf_handle_error(location,NF_PUT_VARA_REAL (ncid, psid, ndimstart1, ndimcount1, &
         real(pressures(0:GC_nlayers), kind=4)))
    CALL netcdf_handle_error(location,NF_PUT_VARA_REAL (ncid, tsid, ndimstart1, ndimcount1, &
         real(temperatures(0:GC_nlayers), kind=4)))
    
    ! write 1D with laydim aircol
    ndimstart1 = (/ 1 /)
    ndimcount1 = (/ GC_nlayers /)  
    CALL netcdf_handle_error(location,NF_PUT_VARA_REAL (ncid, airid,  ndimstart1, ndimcount1, &
         real(aircolumns(1:GC_nlayers), kind=4)))
    CALL netcdf_handle_error(location,NF_PUT_VARA_REAL (ncid, aer0id, ndimstart1, ndimcount1, &
         real(taer_profile(1:GC_nlayers), kind=4)))
    CALL netcdf_handle_error(location,NF_PUT_VARA_REAL (ncid, cld0id, ndimstart1, ndimcount1, &
         real(tcld_profile(1:GC_nlayers), kind=4)))

    ! 2D, variables, laydim, gasdim
    ndimstart2 = (/ 1, 1 /)
    ndimcount2 = (/ GC_nlayers, ngases /)  
    CALL netcdf_handle_error(location,NF_PUT_VARA_DOUBLE (ncid, gascolid, ndimstart2, &
         ndimcount2, gas_partialcolumns(1:GC_nlayers, 1:ngases)))

    ! 4d variables, nsqdim, vzadim, azadim, szadim
    if (do_brdf_surface) then
       ndimstart4 = (/ 1, 1, 1, 1 /)
       ndimcount4 = (/ NSTOKESSQ, GC_n_view_angles, GC_n_azimuths, GC_n_sun_positions /) 
       CALL netcdf_handle_error(location,NF_PUT_VARA_REAL (ncid, brdfid, ndimstart4, ndimcount4, &
            real(Total_brdf(1:NSTOKESSQ,1:GC_n_view_angles,1:GC_n_azimuths,1:GC_n_sun_positions), kind=4)))
       if (OUTPUT_WSABSA) then
          ndimstart1 = (/ 1 /)
          ndimcount1 = (/ 1 /)  
          CALL netcdf_handle_error(location,NF_PUT_VARA_REAL (ncid, wsaid, ndimstart1, ndimcount1, real(WSA_CALCULATED, kind=4)))
          CALL netcdf_handle_error(location,NF_PUT_VARA_REAL (ncid, bsaid, ndimstart1, ndimcount1, real(BSA_CALCULATED, kind=4)))
       endif
    endif

    !==============================================================================
    ! CLOSE the NetCDF file
    !==============================================================================
    
    CALL netcdf_handle_error(location,NF_CLOSE(ncid))
    
  END SUBROUTINE Create_netcdf_output_file

  SUBROUTINE netcdf_handle_error(location,status)

    IMPLICIT NONE
    INCLUDE '../netcdf/netcdf.inc'

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER, INTENT(IN) :: status
    CHARACTER(*)        :: location

    ! ---------------
    ! Local variables
    !----------------
    CHARACTER(NF_MAX_NAME) :: message

    ! Code starts here
    IF (status .NE. NF_NOERR) THEN
       message = 'Error in '//TRIM(location)//': '//NF_STRERROR(status)
       CALL write_err_message(.TRUE., message)
       CALL error_exit(.TRUE.)
    ENDIF

  END SUBROUTINE netcdf_handle_error
  
END MODULE GC_netcdf_module
